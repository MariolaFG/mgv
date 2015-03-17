import ipdb

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
import logging
logging.basicConfig()
log = logging.getLogger()

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
        if target_region[0] is None:
            contigs = bam.references
            start = None
            end = None
        else:
            contigs = [target_region[0]]
            start = target_region[1]
            end = target_region[2]
    else:
        contigs = bam.references
        start = None
        end = None
    print contigs
    
    for i, ctg in enumerate(contigs):
        if target_taxa:
            taxid, cid = ctg.split(".", 1)
        else:
            taxid, cid = None, ctg
            
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

        log.info("Loading site coverage...")
        cov = defaultdict(int)
        for col in bam.pileup(ctg, start, end):
            for pileupread in col.pileups:
                nt = pileupread.alignment.query_sequence[pileupread.query_position]  
                cov[(col.pos, nt)] += 1
        log.info("Loading read content...")
       
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
            yield taxid, cid, paired_reads, total_reads, cov 
    
def get_paired_haplotypes(paired_reads, refseq, coverage, min_snp_coverage, min_site_quality, min_snps_in_haplotype):
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

        # build paired haplotypes. Pysam loaded read 0-based positions, so we
        # can directly compare to refseq string
        read_snps = []
        covered_regions = []
        for read in read_pair:
            seq = refseq[read.reference_start:read.reference_end]
            for i, nt in enumerate(read.query_alignment_sequence):
                if nt != seq[i] and read.query_qualities[i] >= min_site_quality:
                    genomic_coord = i + read.reference_start
                    snp_coverage[ (genomic_coord, nt) ] += 1
                    if coverage[(genomic_coord, nt)] >= min_snp_coverage:
                        read_snps.append( (genomic_coord, nt) )
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
        log.debug('using %s as seed' %(haplotypes_list[i]))
        compatible_haps = [i]
        query_hap = set(current_hap)
        query_hap_regs = list(flattened_regions[i])
        target_visited = set()
        retry = True
        while retry:
            retry = False
            for j, target_hap in enumerate(haplotypes_list):
                if j == i: continue # do not merged with itself
                if j in target_visited: continue 
                if j in visited: continue ###################

                target_hap_regs = flattened_regions[j]
                if mergeable_haplotypes(query_hap, query_hap_regs,
                                        target_hap, target_hap_regs,
                                        snp2coverage, 
                                        min_snp_coverage,
                                        min_anchoring_snps):
                    log.debug('  Merged with %s' %target_hap)
                    # if compatible, build a merged haplotype and continue
                    # searching. Add also this hap to the visited list, so it will
                    # not be used as seed hap.
                    visited.add(j)
                    query_hap.update(target_hap)
                    query_hap_regs.extend(target_hap_regs)
                    compatible_haps.append(j)

                    target_visited.add(j)
                    log.debug('  found mergeable %s ' %(haplotypes_list[j]))
                    retry = True
                    #break
                    
                        
        merged_haplotypes += 1
        hap_sizes.append(len(query_hap))
        hap_group = [haplotypes_list[x] for x in compatible_haps]
        #if len(hap_group) > 1:
        #    print_hap_group(hap_group)
        yield hap_group, query_hap, query_hap_regs
        
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
            t1 = time.time() - 0.00001
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
    log.debug(" Comparing %s vs %s" %(a, b))
    # builds a list of snps present both in a and b haplotypes
    common_snps = [snp for snp in a & b if snp2coverage[snp]] # >= min_snp_coverage]
    
    #if len(common_snps) and (len(b) == len(common_snps) or len(a) == len(common_snps)):
        #if one haplotype is contained in the other, they can safely be merged
    #    log.debug(" ACCEPT: compatible NO new SNPs")
    #    return True        
    if len(common_snps) < min_anchoring_snps and a^b:
        log.debug(" DISCARD: Not enough anchoring points to join %s" %len(common_snps))
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
                    log.debug(" DISCARD: contradicting (B) regions: %s in %s" %(unique_a[0], rb))
                    return False
                #if overlap_reg((unique_a[0], unique_a[0]), rb):
                #    return False
                    
        for unique_b in b - a:
            for ra in flat_regs_a:
                if snp_contined_in_region(unique_b[0], ra):
                    log.debug(" DISCARD: contradicting (A) regions: %s in %s" %(unique_b[0], ra))
                    return False
                #if overlap_reg((unique_b[0], unique_b[0]), ra):
                #    return False

    log.debug(" ACCEPT: new snps added")
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

def print_hap_group(hap_groups):
    all_snps = set()
    for hap in hap_groups:
        all_snps.update(hap)
    for hap in hap_groups:
        for snp in sorted(all_snps):
            if snp in hap:
                #print "% 5d-%s" %(snp[0], snp[1]),
                print snp[1],
            else:
                print " ",
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

    master_line = []
    for pos in xrange(len(new_lines[0])):
        column_chars = set([ln[pos] for ln in new_lines[1:] if not ln.startswith(">") and ln[pos] != ' '])
        if len(column_chars) > 1:
            master_line.append(str(len(column_chars)))
        else:
            master_line.append('.')
    #new_lines.append(''.join(master_line))
    return new_lines, ''.join(master_line)
    
def view_reads2(hap_blocks, contigseq, show_grouped=True, snp_only=False, cid=""):
    import readview

    global_start = len(contigseq)
    global_end = 0

    for hap_group, hap2regions, title in hap_blocks:
        min_pos, max_pos = global_region(hap_group, hap2regions)
        global_start = min(min_pos, global_start)
        global_end = max(global_end, max_pos)

    print global_start, global_end
    
    refseq = contigseq[global_start:global_end+1]
    all_lines = [refseq]

    if show_grouped:
        for hap_group, hap2regions, title in hap_blocks:
            title = '> %s (%d read pairs) - (%d snps)' %(title, len(hap_group), len(merged_snps(hap_group)))
            line = [' '] * len(refseq)
            for hap in hap_group:
                for reg in hap2regions[hap]:
                    for read_pair in reg:
                        for i in xrange(read_pair[0], read_pair[1] + 1):
                            if i >= global_start and i<= global_end:
                                if line[i - global_start] == ' ':
                                    line[i - global_start] = '-'
                for pos, nt in hap:
                    line[pos - global_start] = nt
            all_lines.append(''.join(line))
    else:
        for hap_group, hap2regions, title in hap_blocks:
            title = '> %s (%d read pairs) - (%d snps)' %(title, len(hap_group), len(merged_snps(hap_group)))
            hap_group_lines = []
            for hap in hap_group:
                for reg in hap2regions[hap]:
                    line = [' '] * len(refseq)
                    for read_pair in reg:
                        read_start = read_pair[0] - global_start
                        read_end = read_pair[1] - global_start
                        line[read_start:read_end] = ['-'] * (read_end - read_start + 1)

                    for pos, letter in hap:
                        line[pos - global_start] = letter
                    hap_group_lines.append(''.join(line))

            hap_group_lines.sort()
            hap_group_lines.append('>'+'_'*len(line))
            all_lines.extend(hap_group_lines)
        
    if snp_only: 
        new_lines, master_line = clean_lines(all_lines)
        header = [new_lines[0], master_line]
        target_lines = new_lines[1:]
    else:
        ruler = ''.join(["|" if not x%10 else " " for x in xrange(1, len(all_lines[0]))])
        header = [ruler, all_lines[0]]
        target_lines = all_lines[1:]
       
    if show_grouped:
        target_lines.sort(reverse=True)
    readview.view(header + target_lines, "%s:%d-%d: %d merged read pairs, %d total SNPs." %(cid,global_start, global_end, len(hap_group), len(merged_snps(hap_group))))

    return 
   
def main(args):
    if len(args.bamfiles) == 1:
        target_bamfiles = glob(args.bamfiles[0])
    else:
        target_bamfiles = args.bamfiles

    for reg in args.target_regions:
        target_region_haps = defaultdict(list)
        for bamfile in target_bamfiles:
            # Iterate over groups of paired reads
            
            for taxid, contigid, paired_reads, total_reads, cov in iter_reads_by_contig(bamfile, target_taxa=args.taxa, target_region=reg):
                print 'contig:', contigid, taxid
                if args.refseqs:
                    contigseq = args.refseqs.get_seq(contigid)
                   
                else:
                    contigseq = contig_db.find_one({"sp":int(taxid), "c":contigid}, {"nt"})["nt"]
                print "Contig lentgh:", len(contigseq)
                
                
                # generate basic haplotypes based on paired reads containing a minimum number of snps. 
                haplotype2regions, snp_cov = get_paired_haplotypes(paired_reads, contigseq, cov,
                                                                   min_snp_coverage=args.min_snp_coverage,
                                                                   min_site_quality=args.min_site_quality,
                                                                   min_snps_in_haplotype=args.min_snps_in_haplotype)

                gene_coords = None
                if args.gene_snps_only:
                    print 'Loading gene coords...'
                    gene_coords = []
                    for g in gene_db.find({"sp":int(taxid), "c":contigid}, {"n":1, "s":1, "e":1}):
                        gene_coords.append([g["s"], g["e"], g["n"]])
                    print len(gene_coords)
                
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

                # merges paried-haplotypes based non-conflicting overlaping snps
                # from different reads.
                merged = []
                for hap_group, merged_snps, flatten_regs in merge_haplotypes(haplotype2regions, snp_cov,
                                                                             min_anchoring_snps=args.min_anchoring_snps,
                                                                             min_snp_coverage=args.min_snp_coverage):

                    if len(hap_group):
                        merged.append([hap_group, merged_snps, flatten_regs])
                        if args.view and len(merged_snps) >= args.view:
                            view_reads2([[hap_group, haplotype2regions, bamfile]], contigseq, show_grouped=False, snp_only=True, cid=contigid)
                        if args.groupedview and len(merged_snps) >= args.groupedview:
                            target_region_haps[(len(contigseq), contigid, contigseq)].append([hap_group, haplotype2regions, bamfile])
            
                if args.output:
                    cPickle.dump(merged, open(args.output, "w"))
                
            if args.groupedview:
                for (clen, cid, cseq), region_haps in sorted(target_region_haps.items(), reverse=True, cmp=lambda x,y: cmp(x[0][0], y[0][0])):
                    #get_master_haplotype(region_haps)
                    view_reads2(region_haps, cseq, show_grouped=True, cid=cid, snp_only=True)

def get_master_haplotype(hap_groups_in_region):
    #    ipdb.set_trace()
    counter = defaultdict(set)
    flat = []
    for hap_group, regions, _ in hap_groups_in_region:
        flat.extend(hap_group)
        for hap in hap_group:
            snp_pos = set([snp[0] for snp in hap])
            for regs in regions[hap]:
                for r in regs:
                    for i in xrange(r[0], r[1]+1):
                        if i not in snp_pos:
                            counter[i].add('-')
            for snp in hap:
                counter[snp[0]].add(snp[1])
    print_hap_group(flat)
    
    master_line = ''.join(["." if pos[1]==1 else str(pos[1]) for pos in sorted([(pos, len(v)) for pos,v in counter.iteritems() if len(v)>1 or '-' not in v ])])
    print master_line
    return counter
      

                    

def test():
    test_contig = "GAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG"    
    test_hap2regions = {
        frozenset([(0, "C"), (9, "T"), (19, "G")]): [ [(0, 13), (17, 24)] , [(1, 32)] ],
        frozenset([(0, "C"), (9, "T"), (39, "C")]): [ [(0, 19), (22, 40)] ],
        frozenset([(30, "C"), (35, "G")]): [ [(20, 60), (21, 54)] , [(20, 60)] ]
    }
    test_snp2coverage = defaultdict(lambda:100)
    target_region_haps = defaultdict(list)
    for hap_group, merged_snps, flat_regs in  merge_haplotypes(test_hap2regions, test_snp2coverage, 10, 2):
        #view_reads2([[hap_group, test_hap2regions, "test example"]],
        #            test_contig, show_grouped=False, snp_only=False, cid="test")
        target_region_haps[(len(test_contig), "test", test_contig)].append([hap_group, test_hap2regions, "test contig"])

    raw_input()
        
    for (clen, cid, cseq), region_haps in sorted(target_region_haps.items(), reverse=True, cmp=lambda x,y: cmp(x[0][0], y[0][0])):
        view_reads2(region_haps, cseq, show_grouped=False, cid=cid, snp_only=False)
        #view_reads2(region_haps, cseq, show_grouped=True, cid=cid, snp_only=True)
    
                    
if __name__ == '__main__':
    
    if len(sys.argv) == 2 and sys.argv[1] == "test":
        log.setLevel(level=logging.DEBUG)
        test()
        sys.exit(0)
    
    # GENERAL CONFIG
    parser = ArgumentParser()
    parser.add_argument('--debug', dest="debug", action="store_true")
    parser.add_argument('--taxa', dest="taxa", default=[], nargs="*")
    parser.add_argument('--bam', dest="bamfiles", nargs='+')
    parser.add_argument('--target_regions', dest='target_regions', type=str, nargs="+")
    parser.add_argument('--target_genes', dest='target_genes', type=str, nargs="+")
    parser.add_argument('--contig', dest="contig")
    parser.add_argument('--min_snp_coverage', dest='min_snp_coverage', default=5, type=int)
    parser.add_argument('--min_site_quality', dest='min_site_quality', default=13, type=int)
    parser.add_argument('--min_anchoring_snps', dest='min_anchoring_snps', default=2, type=int)
    parser.add_argument('--min_snps_in_haplotype', dest='min_snps_in_haplotype', default=2, type=int)
    parser.add_argument('--output', dest='output')
    parser.add_argument('--view', dest='view', type=int)
    parser.add_argument('--groupedview', dest='groupedview', type=int)
    parser.add_argument('--gene_snps_only', dest='gene_snps_only', action='store_true')
    parser.add_argument('--refseqs', dest='refseqs', type=str)
    args = parser.parse_args()

    if args.debug:
        log.setLevel(level=logging.DEBUG)
    
    if args.target_regions:
        regions = []
        for reg in args.target_regions:
            try:
                contig, raw_pos = reg.split(':')
                if not raw_pos:
                    start, end = None, None
                else:
                    start, end = map(int, raw_pos.split('-'))
                regions.append([contig.strip(), start, end])
            except Exception:
                print >>sys.stderr, 'ERROR: Invalid contig region. Use contigid:start-end syntax\n'
                raise
        args.target_regions = regions
        
    if args.target_genes:
        # extract region of genes and add it to target_regions
        pass

    if not args.target_regions:
        # If not regions requested, scan all contigs completely 
        args.target_regions = [ [None, None, None] ]

    if args.refseqs:
        from ete2 import SeqGroup
        args.refseqs = SeqGroup(args.refseqs)
                
    if args.output and os.path.exists(args.output):
        print 'Output file exits. Skipping'
    else:
        main(args)



