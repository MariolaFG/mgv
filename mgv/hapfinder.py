import sys
from glob import glob
from collections import defaultdict
from argparse import ArgumentParser

import numpy
from pymongo import MongoClient
import pysam

import readview

# Database pointers
mongo = MongoClient()
snp_db =  mongo.metagenvar_filtered.snp
sites_db =  mongo.metagenvar_filtered.sites
samples_db = mongo.metagenvar_filtered.samples
gene_db = mongo.ref9.genes
contig_db = mongo.ref9.contigs

def iter_reads_by_contig(bamfile, target_taxa=None):
    print "Loading reads..."
    bam = pysam.AlignmentFile(bamfile, "rb" )
    contigs = bam.references
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
        for read in bam.fetch(ctg):
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

    print "\r   reads:% 8d, haplotypes:% 5d, max-haplotype-size:% 2d, unique-snps:% 7d" %\
        (nread, len(haplotypes), max([0] + [len(s) for s in haplotypes]), len(snp_coverage)),
    print
    return haplotypes, snp_coverage
    
def merge_haplotypes(haplotype2regions, snp2coverage, min_snp_coverage, min_anchoring_snps):
    # haplotype2regions = { set(snp1, snp2, ....): regions covered by original reads } 
    
    def sort_haplotypes(hap1, hap2):
        # by number of snps in haplotype
        r = cmp(len(hap1), len(hap2))
        if r == 0:
            # otherwise by number of read_pairs covering 
            r = cmp(len(haplotype2regions[hap1]), len(haplotype2regions[hap1]))
        return r
  
    merged_haplotypes = []
    haplotypes_list = haplotype2regions.keys()
    haplotypes_list.sort(sort_haplotypes, reverse=True)
    cycles = 0
    visited = set()
    while haplotypes_list:
        cycles += 1
        
        # extract first (smallest haplotype)
        #current_hap = haplotypes_list.pop()
        current_hap = haplotypes_list[0]
        haplotypes_list = haplotypes_list[1:]
        
        current_hap_regs = haplotype2regions[current_hap]
        #print current_hap, len(current_hap_regs)
        
        compatible_haps = []
        search_compatibles = True
        query_hap = current_hap
        query_hap_regs = current_hap_regs
        while search_compatibles:
            cycles += 1
            search_compatibles = False
            size = len(haplotypes_list)
            for index, target_hap in enumerate(haplotypes_list): # reverse?
                target_hap_regs = haplotype2regions[target_hap]
                if mergeable_haplotypes(query_hap, query_hap_regs,
                                        target_hap, target_hap_regs,
                                        snp2coverage, 
                                        min_snp_coverage,
                                        min_anchoring_snps):
                    # if compatible, build a merged haplotype and continue
                    # searching
                    query_hap = query_hap | target_hap
                    query_hap_regs = query_hap_regs + target_hap_regs
                    compatible_haps.append(target_hap)
                    #del haplotypes_list[(size - index)-1]
                    #del haplotypes_list[index]
                    haplotypes_list = [e for j, e in enumerate(haplotypes_list) if j != index]
                    
                    search_compatibles = True
                    break
                    
        compatible_haps.append(current_hap)
        merged_haplotypes.append(compatible_haps)
        if cycles % 1000 == 0:
            snp_count = [len(get_merged(x)) for x in merged_haplotypes]
            print "\r   Haps left:%d merged:%d, sizes:%d-%d median-size:%d" %(len(haplotypes_list),
                                                                              len(merged_haplotypes),
                                                                              max(snp_count),
                                                                              min(snp_count),
                                                                              numpy.median(snp_count))
            sys.stdout.flush()                                                                    
            #print_haplotypes(compatible_haps, haplotype2regions)

            
    snp_count = [len(get_merged(x)) for x in merged_haplotypes]
    print "\r   Haps left:%d merged:%d, sizes:%d-%d median-size:%d" %(len(haplotypes_list),
                                                                      len(merged_haplotypes),
                                                                      max(snp_count),
                                                                      min(snp_count),
                                                                      numpy.median(snp_count))

    crappyhist(snp_count, max(snp_count))
    return merged_haplotypes


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

    
def view_reads(target, alleles, ref, groups=None):
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

    
def get_merged(haplotypes):
    unique_snps = set()
    for hap in haplotypes:
        unique_snps.update(hap)
    return unique_snps

    
def mergeable_haplotypes(a, regs_a, b, regs_b, snp2coverage, min_snp_coverage, min_anchoring_snps):
    '''Returns True if two haplotypes (composed of a series of SNP coords) could be
    merged, meaning that they are probably from the same DNA fragment.
    '''
    
    def overlap_reg(ra, rb):
        '''Returns the overalp of two genomic regions''' 
        return min(ra[1], rb[1]) - max(ra[0], rb[0]) >= 0

    # builds a list of snps present both in a and b haplotypes
    common_snps = [snp for snp in (a & b) if snp2coverage[snp] >= min_snp_coverage]
        
    if len(common_snps) < min_anchoring_snps:
        return False
    else:
        flat_regs_b = []
        [flat_regs_b.extend(breg) for breg in regs_b]
        flat_regs_a = []
        [flat_regs_a.extend(areg) for areg in regs_a]
        
        # check for incompatibilities between the 
        for unique_a in a - b:
            for rb in flat_regs_b:
                if overlap_reg((unique_a[0], unique_a[0]), rb):
                    return False
                    
        for unique_b in b - a:
            for ra in flat_regs_a:
                if overlap_reg((unique_b[0], unique_b[0]), ra):
                    return False
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



def main(args):
    for bamfile in glob(args.bamfiles):
        for taxid, contigid, paired_reads, total_reads in iter_reads_by_contig(bamfile, target_taxa=args.taxa):
            refseq = contig_db.find_one({"sp":int(taxid), "c":contigid}, {"nt"})["nt"]
            haplotype2regions, snp_cov = get_paired_haplotypes(paired_reads, refseq, 
                                                              min_site_quality=args.min_site_quality,
                                                              min_snps_in_haplotype=args.min_snps_in_haplotype)
            merged = merge_haplotypes(haplotype2regions, snp_cov,
                             min_anchoring_snps=args.min_anchoring_snps,
                             min_snp_coverage=args.min_snp_coverage)


            if args.output:
                import cPickle
                cPickle.dump(merged, open(args.output, "w"))

            # Visualize
            # for m in merged:
            #     print_haplotypes(m, refseq)
            #     raw_input()
            #     view_reads(m, haplotype2regions, refseq)
            #     print
    
                

            

if __name__ == '__main__':
    # GENERAL CONFIG
    parser = ArgumentParser()
    parser.add_argument('--debug', dest="debug")
    parser.add_argument('--taxa', dest="taxa", default=[None], nargs="*")
    parser.add_argument('--bam', dest="bamfiles")
    parser.add_argument('--contig', dest="contig")
    parser.add_argument('--min_snp_coverage', dest='min_snp_coverage', default=4, type=int)
    parser.add_argument('--min_site_quality', dest='min_site_quality', default=13, type=int)
    parser.add_argument('--min_anchoring_snps', dest='min_anchoring_snps', default=2, type=int)
    parser.add_argument('--min_snps_in_haplotype', dest='min_snps_in_haplotype', default=2, type=int)
    parser.add_argument('--output', dest='output')
    args = parser.parse_args()
    main(args)


    

