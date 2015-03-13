import sys
import os
from glob import glob
from collections import defaultdict
from argparse import ArgumentParser
import time
import cPickle
from pprint import pprint
#import numpy
from pymongo import MongoClient
import pysam

# Database pointers
mongo = MongoClient('pc-bork62')
snp_db =  mongo.metagenvar_filtered.snp
sites_db =  mongo.metagenvar_filtered.sites
samples_db = mongo.metagenvar_filtered.samples
gene_db = mongo.ref9.genes
contig_db = mongo.ref9.contigs

def iter_reads_by_contig(bamfile, target_taxa=None, target_region=None):
    print "Loading reads..."
    bam = pysam.AlignmentFile(bamfile, "rb" )
    if target_region:
        contigs = [target_region[0]]
        start = target_region[1]
        end = target_region[2]
    else:
        contigs = bam.references
        start = None
        end = None
        

    for i, ctg in enumerate(contigs):
        taxid, cid = ctg.split(".", 1)
        # Skip contigs from unwanted species
        if target_taxa and taxid not in target_taxa:
            continue
                
        print "Loading % 40s reads (task: %d/%d)" %(ctg, i, len(contigs))
        paired_reads = {}
        indel_reads = 0
        total_reads = 0
        # Collect and parse reads
        print "reads:"
        print bam.count(ctg, start, end)
        sys.stdout.flush()
        for read in bam.fetch(ctg, start, end):
            # if the read matches several blocks, skip it for now: probably
            # pointing to an indel
            if len(read.blocks) > 1:
                indel_reads += 1
                continue
            # store read with its pair
            paired_reads.setdefault(read.qname[:-3], []).append(read)
            total_reads += 1
            
        if paired_reads:
            yield taxid, cid, paired_reads, total_reads 
    
def get_paired_haplotypes(paired_reads, refseq, min_site_quality, min_snps_in_haplotype):
    haplotypes = defaultdict(list)
    snp_coverage = defaultdict(int)
    for nread, read_pair in enumerate(paired_reads.itervalues()):
        # print progress
        if nread % 1000 == 0:
            if nread % 100000 == 0:
                print "\r   reads:% 8d, haplotypes:% 5d, max-haplotype-size:% 2d, unique-snps:% 7d" %\
                    (nread, len(haplotypes), max([0] + [len(s) for s in haplotypes]), len(snp_coverage)),
            else:
                print "\r   reads:% 8d" %(nread),
            sys.stdout.flush()

        # build paired haplotypes    
        read_snps = []
        covered_regions = []
        for read in read_pair:
            seq = refseq[read.reference_start:read.reference_end]
            for i, nt in enumerate(read.query_alignment_sequence):
                if nt != seq[i] and read.query_qualities[i] >= min_site_quality:
                    genomic_coord = i + read.reference_start
                    read_snps.append( (genomic_coord, nt) )
                    snp_coverage[ (genomic_coord, nt) ] += 1
            covered_regions.append((read.reference_start, read.reference_end, read.qname))
            
        if len(read_snps) > min_snps_in_haplotype:
            hap_key = frozenset(read_snps)
            # haplotype coverage could be guessed by len(haplotypes[hap_key])
            # I need to keep convered regions for subsequent analyses
            haplotypes[hap_key].append(covered_regions)
        #else:
        #    haplotypes[frozenset(read_snps)].append(covered_regions)

    print "\r   reads:% 8d, haplotypes:% 5d, max-haplotype-size:% 2d, unique-snps:% 7d" %\
        (nread, len(haplotypes), max([0] + [len(s) for s in haplotypes]), len(snp_coverage)),
    print
    return haplotypes, snp_coverage
    
def merge_haplotypes(haplotype2regions, snp2coverage, min_snp_coverage, min_anchoring_snps):
    # haplotype2regions = { set(snp1, snp2, ....): regions covered by original reads } 
    
    def sort_haplotypes_by_size(hap1, hap2):
        # by number of snps in haplotype
        r = cmp(len(hap1), len(hap2))
        if r == 0:
            s1 = min([s[0] for s in hap1])
            s2 = min([s[0] for s in hap2])
            r = cmp(s1, s2) * -1
            # otherwise by number of read_pairs covering 
            #r = cmp(len(haplotype2regions[hap1]), len(haplotype2regions[hap1]))
        return r
        
    def sort_haplotypes_by_pos(hap1, hap2):
        # by number of snps in haplotype
        #s1, e1 = global_region(hap1, haplotype2regions)
        #s2, e2 = global_region(hap2, haplotype2regions)
        s1 = min([s[0] for s in hap1])
        s2 = min([s[0] for s in hap2])
        
        return cmp(s1, s2) * -1
        
    haplotypes_list = haplotype2regions.keys()
    haplotypes_list.sort(sort_haplotypes_by_size, reverse=True)

    print 'Flattening genomic regions'
    flattened_regions = []
    for hap in haplotypes_list:
        flat_regs = []
        [flat_regs.extend(reg) for reg in haplotype2regions[hap]]
        flattened_regions.append(flat_regs)

    visited = set()
    merged_haplotypes = 0
    hap_sizes = [0]
    t1 = time.time()
    scanned = 0
    log_scanned = 0
   
    for i, current_hap in enumerate(haplotypes_list):
        if i in visited: continue
        
        visited.add(i)
        
        compatible_haps = [i]
        query_hap = set(current_hap)
        query_hap_regs = list(flattened_regions[i])
        # for j, target_hap in enumerate(haplotypes_list[i:]):
        #     orig_pos = j + i
        #     if orig_pos in visited: continue
        target_visited = set()
        last_pos = 0
        while last_pos+1 < len(haplotypes_list):
            for orig_pos, target_hap in enumerate(haplotypes_list):
                last_pos = orig_pos
                if orig_pos == i: continue # do not merged with itself
                if orig_pos in target_visited: continue 
                if orig_pos in visited: continue ###################

                target_hap_regs = flattened_regions[orig_pos]
                if mergeable_haplotypes(query_hap, query_hap_regs,
                                        target_hap, target_hap_regs,
                                        snp2coverage, 
                                        min_snp_coverage,
                                        min_anchoring_snps):

                    # if compatible, build a merged haplotype and continue
                    # searching. Add also this hap to the visited list, so it will
                    # not be used as seed hap.
                    visited.add(orig_pos)
                    query_hap.update(target_hap)
                    query_hap_regs.extend(target_hap_regs)
                    compatible_haps.append(orig_pos)

                    target_visited.add(orig_pos)
                    break

                
        merged_haplotypes += 1
        hap_sizes.append(len(query_hap))
        
        yield [haplotypes_list[x] for x in compatible_haps], query_hap, query_hap_regs
        
        if scanned == log_scanned:
            if scanned == 0:
                etime = time.time()-t1 # time used to process one hap
                log_scanned = int(5 / etime)
                print "max scan time:", etime, "logging time adjusted to show progress every ~5 secs. ", log_scanned, "haps" 
            else: 
                etime = (scanned/(time.time()-t1)) 
            print "     scanned:%d merged:%d, sizes:%d-%d median-size:%d time:%d/s" %(len(visited),
                                                                                       merged_haplotypes,
                                                                                       max(hap_sizes),
                                                                                       min(hap_sizes),
                                                                                       0,
                                                                                       etime)
            sys.stdout.flush()
            t1 = time.time()
            scanned = 0
            #print_haplotypes(compatible_haps, haplotype2regions)
        scanned += 1

    print "     scanned:%d merged:%d, sizes:%d-%d median-size:%d time:%d/s" %(len(visited),
                                                                              merged_haplotypes,
                                                                              max(hap_sizes),
                                                                              min(hap_sizes),
                                                                              0,
                                                                              scanned / (time.time()-t1))

    #crappyhist(snp_count, max(snp_count))
    #return merged_haplotypes


def flatten_snps(items, haplotypes):
    unique_snps = set()
    map(lambda i: unique_snps.update(haplotypes[i]), items)
    return unique_snps
   
    
def crappyhist(a, bins):
    '''Draws a crappy text-mode histogram of an array'''
    import numpy as np
    import string
    from math import log10

    h,b = np.histogram(a, bins)

    for i in range (0, bins-1):
        print string.rjust(`b[i]`, 7)[:int(log10(np.amax(b)))+5], '| ', '#'*int(70*h[i-1]/np.amax(h))
        #print string.rjust(`b[bins]`, 7)[:int(log10(np.amax(b)))+5] 
        
        
    
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

def snp_contined_in_region(snp_pos, reg):
    if snp_pos < reg[0]: return False
    if snp_pos > reg[1]: return False
    return True

def overlap_reg(ra, rb):
    '''Returns the overalp of two genomic regions'''    
    return (min(ra[1], rb[1]) - max(ra[0], rb[0])) >= 0
  
    
def mergeable_haplotypes(a, flat_regs_a, b, flat_regs_b, snp2coverage, min_snp_coverage, min_anchoring_snps):
    '''Returns True if two haplotypes (composed of a series of SNP coords) could be
    merged, meaning that they are probably from the same DNA fragment.
    '''
    
    # builds a list of snps present both in a and b haplotypes
    common_snps = [snp for snp in a & b if snp2coverage[snp] >= min_snp_coverage]
    if len(common_snps) and (len(b) == len(common_snps) or len(a) == len(common_snps)):
        #if one haplotype is contained in the other, the can safely be merged
        return True        
    elif len(common_snps) < min_anchoring_snps:
        # if they contribute different snps, but share less and the anchoring thr, cannot be merged 
        return False
    else:
        #flat_regs_b = []
        #[flat_regs_b.extend(breg) for breg in regs_b]
        #flat_regs_a = []
        #[flat_regs_a.extend(areg) for areg in regs_a]
        
        # does any SNP in one profile overlap with the sampled region of the
        # other haplotype?
        for unique_a in a - b:
            for rb in flat_regs_b:
                if snp_contined_in_region(unique_a[0], rb):
                    return False
                #if overlap_reg((unique_a[0], unique_a[0]), rb):
                #    return False
                    
        for unique_b in b - a:
            for ra in flat_regs_a:
                if snp_contined_in_region(unique_b[0], ra):
                    return False
                #if overlap_reg((unique_b[0], unique_b[0]), ra):
                #    return False
    return True

def print_haplotypes(alleles_group, alleles):
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
    print all_pos

def view_reads(target, alleles, ref, groups=None):
    import readview
    
    break_points = {}
    if groups:
        target = []
        for gr in groups:
            target.extend(gr)
            break_points[len(target)-1] = str(gr)
    print
    

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
            group_lines = []
            group_hap = set()

    if group_lines:
        group_lines.sort()
        lines.append(group_lines)
        
    lines.sort()
    flat_lines = [ln[1] for sublist in lines for ln in sublist]                     
    hap_length = len(set([val for sublist in target for val in sublist]))
    title = "%d:%d (%dbp), read-pairs:%d, hap_length:%d" %(min_pos, max_pos, max_pos-min_pos, len(target), hap_length)
    readview.view([refseq] + flat_lines, title)

def print2(hap_groups):
    all_snps = set()
    for hap in hap_groups:
        all_snps.update(hap)
    for hap in hap_groups:
        for snp in sorted(all_snps):
            if snp in hap:
                print "% 5d-%s" %(snp[0], snp[1]),
            else:
                print "% 5d- " %(snp[0]),
        print

def merged_snps(hap_groups):
    all_snps = set()
    for hap in hap_groups:
        all_snps.update(hap)
    return sorted(all_snps)

def clean_lines(lines):
    new_lines = [""]*len(lines)
    skiplable = set(['-', ' '])
    
    for pos in xrange(len(lines[0])):
        column_chars = set([ln[pos] for ln in lines[1:] if not ln.startswith(">")])
        if column_chars - skiplable:
            for i, newl in enumerate(new_lines):
                new_lines[i] += lines[i][pos]
    return new_lines
    
def view_reads2(hap_blocks, contigseq, show_grouped=True):
    import readview

    global_start = len(contigseq)
    global_end = 0

    for hap_group, hap2regions, title in hap_blocks:
        min_pos, max_pos = global_region(hap_group, hap2regions)
        global_start = min(min_pos, global_start)
        global_end = max(global_end, max_pos)
    print global_start, global_end
    
    refseq = contigseq[global_start:global_end + 1]
    all_lines = [refseq]
    for hap_group, hap2regions, title in hap_blocks:
        title = '> %s (%d read pairs) - (%d snps)' %(title, len(hap_group), len(merged_snps(hap_group)))
        hap_group_lines = []
        if show_grouped:
            line = [' '] * len(refseq)
        for hap in hap_group:
            for reg in hap2regions[hap]:
                if not show_grouped:
                    line = [' '] * len(refseq)
                for read_pair in reg:
                    read_start = read_pair[0] - global_start
                    read_end = read_pair[1] - global_start
                    line[read_start:read_end] = ['-'] * (read_end-read_start)
                for pos, letter in hap:
                    line[pos - global_start] = letter
                
                if not show_grouped:
                    hap_group_lines.append(''.join(line))
        if show_grouped:
            hap_group_lines.append(''.join(line))
        #print2(hap_group)
        hap_group_lines.sort()
        title = [title+str('_'*(len(refseq)-len(title)))]
        all_lines.extend(hap_group_lines)
    all_lines = clean_lines(all_lines)
    readview.view(sorted(all_lines, reverse=True), "blablbas")
    return 
   

def main(args):
    if len(args.bamfiles) == 1:
        target_bamfiles = glob(args.bamfiles[0])
    else:
        target_bamfiles = args.bamfiles

    target_region_haps = []
    
    for bamfile in target_bamfiles:
        for taxid, contigid, paired_reads, total_reads in iter_reads_by_contig(bamfile, target_taxa=args.taxa, target_region=args.target_region):
            contigseq = contig_db.find_one({"sp":int(taxid), "c":contigid}, {"nt"})["nt"]
            print "Contig lentgh:", len(contigseq)

            gene_coords = None
            if args.gene_snps_only:
                print 'loading gene coords...'
                gene_coords = []
                for g in gene_db.find({"sp":int(taxid), "c":contigid}, {"n":1, "s":1, "e":1}):
                    gene_coords.append([g["s"], g["e"], g["n"]])
                print len(gene_coords)
            haplotype2regions, snp_cov = get_paired_haplotypes(paired_reads, contigseq, 
                                                              min_site_quality=args.min_site_quality,
                                                              min_snps_in_haplotype=args.min_snps_in_haplotype)

            if gene_coords:
                print 'filtering snps outside gene regions'
                new_hap2regions = {}
                for i, hap in enumerate(haplotype2regions):
                    selected_snps = []
                    print '\r% 6d/% 6d' %(i, len(haplotype2regions)),
                    sys.stdout.flush()
                    for s in hap:
                        for g in gene_coords:
                            if s[0] < g[0]: 
                                break # gene_coord should be sorted in ascending orther
                            if snp_contined_in_region(s[0], g[0:2]):
                                selected_snps.append(s)
                                break
                        if len(selected_snps) == len(hap):
                            new_hap2regions[hap] = haplotype2regions[hap]

                print '\n%d haplotypes selected out of %d' %(len(new_hap2regions), len(haplotype2regions))
                haplotype2regions = new_hap2regions
                                
            merged = []
            for hap_group, merged_snps, flatten_regs in merge_haplotypes(haplotype2regions, snp_cov,
                                                                         min_anchoring_snps=args.min_anchoring_snps,
                                                                         min_snp_coverage=args.min_snp_coverage):
                merged.append([hap_group, merged_snps, flatten_regs])

                if args.view and len(merged_snps) >= args.view:
                    #print_haplotypes(hap_group, refseq)
                    #raw_input()
                    if args.target_region:
                        target_region_haps.append([hap_group, haplotype2regions, bamfile])
                    else:
                        view_reads(hap_group, haplotype2regions, contigseq)
                
            if args.output:
                cPickle.dump(merged, open(args.output, "w"))
                
    if args.target_region and args.view:
        view_reads2(target_region_haps, contigseq)
    

if __name__ == '__main__':
    # GENERAL CONFIG
    parser = ArgumentParser()
    parser.add_argument('--debug', dest="debug")
    parser.add_argument('--taxa', dest="taxa", default=[None], nargs="*")
    parser.add_argument('--bam', dest="bamfiles", nargs='+')
    parser.add_argument('--target_region', dest='target_region', type=str)
    parser.add_argument('--contig', dest="contig")
    parser.add_argument('--min_snp_coverage', dest='min_snp_coverage', default=5, type=int)
    parser.add_argument('--min_site_quality', dest='min_site_quality', default=13, type=int)
    parser.add_argument('--min_anchoring_snps', dest='min_anchoring_snps', default=2, type=int)
    parser.add_argument('--min_snps_in_haplotype', dest='min_snps_in_haplotype', default=2, type=int)
    parser.add_argument('--output', dest='output')
    parser.add_argument('--view', dest='view', type=int)
    parser.add_argument('--gene_snps_only', dest='gene_snps_only', action='store_true')
    args = parser.parse_args()

    if args.target_region:
        try:
            contig, raw_pos = args.target_region.split(':')
            args.target_region = [contig.strip()] + map(int, raw_pos.split('-'))
        except Exception:
            print >>sys.stderr, 'ERROR: Invalid contig region. Use contigid:start-end syntax\n'
            raise
        print args.target_region
        
    if args.output and os.path.exists(args.output):
        print 'Output file exits. Skipping'
    else:
        main(args)


    

