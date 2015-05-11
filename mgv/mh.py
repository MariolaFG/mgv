import sys
import pysam
from collections import defaultdict

from pymongo import MongoClient
import readview

mongo = MongoClient()
snp_db =  mongo.metagenvar_filtered.snp
sites_db =  mongo.metagenvar_filtered.sites
samples_db = mongo.metagenvar_filtered.samples
gene_db = mongo.ref9.genes
contig_db = mongo.ref9.contigs

print "Loading reference sequence"
try:
    DEBUG = int(sys.argv[2])
except:
    DEBUG = False
else:
    print "DEBUG MODE, process reads:", DEBUG
    
if sys.argv[1] == "falk":
    FALK = True
    CONTIG="S4__C0_L=313228;"
    TAXID=""
    CONTIG_NAME = CONTIG
    from Bio import SeqIO
    for record in SeqIO.parse("/g/bork1/hildebra/SNP/GNMassSimu//simulated_metaG2SI_3//assemblies/metag/scaffolds.fasta",  "fasta"):
        if record.id == CONTIG:
            REFSEQ = str(record.seq)
            break
    BAMFILE =  "/g/bork1/hildebra/SNP/GNMassSimu//simulated_metaG2SI_3/mapping/Align_ment-smd.bam"
else:
    FALK = False
    CONTIG = "58253.NC_009614"
    TAXID = 435590    
    CONTIG_NAME = "%s.%s" %(TAXID, CONTIG)
    REFSEQ = contig_db.find_one({"sp":TAXID, "c":CONTIG}, {"nt"})["nt"]
    BAMFILE = sys.argv[1]

    
print "Loading reads"
samfile = pysam.AlignmentFile(BAMFILE, "rb" )
paired_reads = {}
total_reads = 0

indel_reads = 0 
for read in samfile.fetch(CONTIG_NAME):
    if len(read.blocks) > 1:
        indel_reads += 1
        continue
    if FALK:
        paired_reads.setdefault(read.qname, []).append(read)
    else:
        paired_reads.setdefault(read.qname[:-3], []).append(read)
    total_reads += 1

print "Total reads", total_reads 
print "Paired_reads", len(paired_reads)
print "read with Indels", indel_reads
print len([p for p in paired_reads.itervalues() if len(p) == 1]), "== 1"
print len([p for p in paired_reads.itervalues() if len(p) == 2]), "== 2"
print len([p for p in paired_reads.itervalues() if len(p) > 2]), " >2 "

# calling paired alleles
alleles = defaultdict(list)

def sort_alleles(alist, a, reverse=True):
    def _sort(x, y):
        r = cmp(len(x), len(y))
        if r == 0:
            r = cmp(len(a[x]), len(a[y]))
        return r
    return sorted(alist, _sort, reverse=reverse)

snps = defaultdict(int)
multi_alg = 0
MIN_QUALITY = 13
for nread, pair in enumerate(paired_reads.itervalues()):
    if nread % 1000 == 0:
        if nread % 100000 == 0:
            #temp = sort_alleles(alleles, alleles)
            #if temp:
            #    for x in temp[:10]:
            #        print "       ", x, alleles[x]
            print "\r% 8d Alleles % 5d, max_size: % 2d, snps: % 7d" %(nread, len(alleles), max([0] + [len(s) for s in alleles]), len(snps)),
        else:
            print "\r% 8d" %(nread),
        sys.stdout.flush()
   
    hap = []
    covered_regions = []
    for read in pair:
        seq = REFSEQ[read.reference_start:read.reference_end]
        var = []
        gvar = []
        for i, x in enumerate(read.query_alignment_sequence):
            if x != seq[i] and read.query_qualities[i] >= MIN_QUALITY:
                if 0:
                    var.append((i, x))
                gcoord = i + read.reference_start
                gvar.append((gcoord, x))
                snps[(gcoord, x)] += 1
        if gvar:
            if 0: 
                seq = [s for s in seq]
                for i, x in var:
                    seq[i] = "*"
                print ''.join(seq)
            hap.extend(gvar)
        covered_regions.append((read.reference_start, read.reference_end, read.qname))
    if len(hap)>1:
        hap_key = frozenset(hap)
        alleles[hap_key].append(covered_regions)
       
    if DEBUG and nread == DEBUG:
        break
print "\r% 8d Alleles % 5d, max_size: % 2d, snps: % 7d" %(nread, len(alleles), max([0] + [len(s) for s in alleles]), len(snps))

MIN_COMMON_SNPS_TO_MERGE = 1
MIN_SNP_COV = 1

def are_compatible_hap(a, regs_a, b, regs_b):
    def overalp_reg(ra, rb):
        return min(ra[1], rb[1]) - max(ra[0], rb[0]) >= 0
    
    common_snps = a & b
    common_snps = [pos for pos in common_snps if snps[pos] >= MIN_SNP_COV]
    
    if len(common_snps) < MIN_COMMON_SNPS_TO_MERGE:
        return False
    else:
        flat_regs_b = []
        [flat_regs_b.extend(breg) for breg in regs_b]
        flat_regs_a = []
        [flat_regs_a.extend(areg) for areg in regs_a]
        
        # check for incompatibilities
        for unique_a in a - b:
            for rb in flat_regs_b:
                #print (unique_a[0], unique_a[0]), rb
                if overalp_reg((unique_a[0], unique_a[0]), rb):
                    #print rb, "oooooooooooooooooooops"
                    return False
        for unique_b in b - a:
            for ra in flat_regs_a:
                #print (unique_b[0], unique_b[0]), ra
                if overalp_reg((unique_b[0], unique_b[0]), ra):
                    #print ra, "oooooooooooooooooooops"
                    return False
    return True

print "\nMerging alleles"
def print_alleles(alleles_group, alleles):
    all_pos = set()
    for a in alleles_group:
        all_pos.update(a)
        
    for pos in all_pos:
        if pos in alleles_group[0]:
            sys.stdout.write(pos[1])
        else:
            sys.stdout.write("-") 
    #print "  % 2d" %len(alleles[alleles_group[0]])
    print 
    # this are expected to be the merged alleles
    for b in alleles_group[1:]:
        for pos in all_pos:
            if pos in b:
                sys.stdout.write(pos[1])
            else:
                sys.stdout.write("-")
        #print "  % 2d" %len(alleles[b])
        print 
    for pos in all_pos:
        sys.stdout.write(pos[1])
    print "  MERGED" 
    print all_pos

   
merged_alleles = []

s_alleles = sort_alleles(alleles, alleles, reverse=False)
while s_alleles:
    a = s_alleles.pop()
    a_regs = alleles[a]
    temp_a = a
    temp_a_regs = a_regs
    
    compatible_b = []
    nocomp = False
    while nocomp == False:
        nocomp = True
        size = len(s_alleles)
        for index, b in enumerate(reversed(s_alleles)):
            b_regs = alleles[b]
            if are_compatible_hap(temp_a, temp_a_regs, b, b_regs):
                temp_a = temp_a | b
                temp_a_regs = temp_a_regs + b_regs
                compatible_b.append(b)
                del s_alleles[(size - index)-1]
                nocomp = False
                break
                
    if compatible_b:
        merged_alleles.append([a] + compatible_b)
        #print_alleles(merged_alleles[-1], alleles)
        print "\r % 6d -% 6d" %(len(s_alleles), len(merged_alleles[-1])),
    else:
        merged_alleles.append([a])
        print "\r % 6d -% 6d" %(len(s_alleles), 1),

print     
print "Total merged", len(merged_alleles)
sizes = [len(m) for m in merged_alleles]
lengths = [len(set([val for sublist in list_of_lists for val in sublist])) for list_of_lists in merged_alleles]

def global_region(target, alleles):
    min_pos, max_pos = None, 0
    for a in target:
        for hap in alleles[a]:
            for reg in hap:
                if min_pos is None or reg[0] < min_pos:
                    min_pos = reg[0]
                if reg[1] > max_pos:
                    max_pos = reg[1]
    return min_pos, max_pos            

def view_reads(target, alleles, ref, groups=None):
    break_points = {}
    if groups:
        target = []
        for gr in groups:
            target.extend(gr)
            break_points[len(target)-1] = str(gr)
    print
    
    print len(break_points), break_points.keys(), "break_points"
    print len(groups)
    raw_input()

    min_pos, max_pos = global_region(target, alleles)
    refseq = ref[min_pos:max_pos]
    
    lines = []
    group_lines = []
    group_hap = set()
    for aindex, a in enumerate(target):
        hap_variants = dict(a)
        for hap in alleles[a]:
            read_seq = [" "] * len(refseq)
            min_read_start = None
            for reg in hap:
                if min_read_start is None or reg[0] < min_read_start:
                    min_read_start = reg[0]
                for i in xrange(min_pos, max_pos):
                    if i>= reg[0] and i <= reg[1]:
                        if i in hap_variants:
                            read_seq[i - min_pos] = hap_variants[i]
                        else:
                            read_seq[i - min_pos] = "-"
            group_lines.append([min_read_start, ''.join(read_seq)])
        group_hap.add(a)
        if aindex in break_points:
            #group_lines.append([len(refseq), '____%s'%sorted(group_hap)+'_'*len(refseq)])
            group_lines.append([len(refseq), '_'*len(refseq)])
            group_lines.sort()
            lines.append(group_lines)
            for x in group_lines:
                print x
            raw_input()
            group_lines = []
            group_hap = set()

    if group_lines:
        group_lines.sort()
        lines.append(group_lines)
        
    lines.sort()
    flat_lines = [ln[1] for sublist in lines for ln in sublist]                     
    hap_length = len(set([val for sublist in m for val in sublist]))
    title = "%d:%d (%dbp), read-pairs:%d, hap_length:%d" %(min_pos, max_pos, max_pos-min_pos, len(m), hap_length)
    
    readview.view([refseq] + flat_lines, title)

    
import numpy
print "min", numpy.min(sizes)
print "max", numpy.max(sizes)
print "mean", numpy.mean(sizes)
print "median", numpy.median(sizes)
print ">1 ", len([s for s in sizes if s>1])

print "min", numpy.min(lengths)
print "max", numpy.max(lengths)
print "mean", numpy.mean(lengths)
print "median", numpy.median(lengths)
print ">1 ", len([s for s in lengths if s>1])


# Merging merged

merged_alleles = [_m for  _l,_s,_m in sorted(zip(lengths, sizes, merged_alleles), reverse=True)]
compatibles = []
used = set()
for i in xrange(len(merged_alleles)):
    print '\r 6d', i, 
    ma = set()
    [ma.update(a) for a in merged_alleles[i]]
    comp = []
    for j in xrange(i, len(merged_alleles)):
        if i != j:
            mb = set()
            [mb.update(a) for a in merged_alleles[j]]
            if len(ma & mb) > 0:
                comp.append(merged_alleles[j])
                
    if len(comp) >= 7:
        #print_alleles(merged_alleles[i]+comp, alleles)
        groups =  [b for b in comp]
        groups.append(merged_alleles[i])
        print 
        print len(groups), "groups"
        raw_input()
        view_reads(target=None, groups = groups, alleles=alleles, ref=REFSEQ)
        print
        
    compatibles.append(len(comp))
print

print numpy.max(compatibles)
print numpy.min(compatibles)
print numpy.mean(compatibles)
print len([_ for _ in compatibles if _==0])
print len([_ for _ in compatibles if _==1])
print len([_ for _ in compatibles if _>1])
raw_input()

# Visualize
for m in merged_alleles:
    minpos, maxpos = global_region(m, alleles)
    print minpos, maxpos, maxpos-minpos
    print "alleles", len(m)
    print "length", len(set([val for sublist in m for val in sublist]))
    print_alleles(m, alleles)
    print 
    view_reads(m, alleles, REFSEQ)
    print
    
                
                

    
    
