from pprint import pprint
import pymongo
import numpy
import cPickle
import sys
import argparse
from scipy import stats
from collections import defaultdict

import zlib
import base64

from pymongo import MongoClient

def get_global_coverage(cov, header_names, sample_indexes):
    full_cov = [-1] * len(sample_indexes)
    for i, sname in enumerate(header_names):
        full_cov[sample_indexes[sname]] = cov[i]
    return full_cov

def dump_snp_matrix(sorted_sp):
    # iter over the list of species sorted by number of pos with snps
    for percent, (nsnps, sp) in enumerate(reversed(sorted_sp)):
        if TARGET_SPECIES and sp not in TARGET_SPECIES:
            continue

        #nsnps = snp.find({'sp':sp}, {'c':1}).count()
        print 'processing', sp, 'with', nsnps, 'snps. ', percent+1, "/", len(sorted_sp)

        positions = [ ]
        bygene_positions = {}
        genes = set()
        gene2pos2cov = {}
        
        for i, site in enumerate(sites_db.find({'sp':sp}, {'c':1, 'pos':1, 'cov':1, 'ref':1, 'g':1}).sort([("c", 1), ("pos", 1)])):
            if args.testcase and i == args.testcase:
                    break
                    
            if i % 1000 == 0:

                print '% 6d % 6d genes:%d\r' %(i, len(positions), len(genes)), 
                sys.stdout.flush()

            if FILTERED_DB:
                species_samples = samples_db.find_one({'sp':sp}, {'samples':1})["samples"]
                full_site_cov = get_global_coverage(site['cov'], species_samples, SAMPLE_INDEX)
            else:
                species_samples = SAMPLE_NAMES
                full_site_cov = site['cov']

            dpos = [None] * SAMPLE_SIZE
            dcov = [None] * SAMPLE_SIZE
            data_array = [0] * SAMPLE_SIZE
            
            if args.genes_only and not site["g"]:
                continue

            genes.add((site["c"], site["g"]))
            
            
            for j, snp in enumerate(snp_db.find({'sp':sp, 'c':site['c'], 'pos':site['pos']}, {'snp':1, 'cov':1, 'syn':1})):
                if FILTERED_DB:
                    full_snp_cov = get_global_coverage(snp['cov'], species_samples, SAMPLE_INDEX)
                else:
                    full_snp_cov = snp['cov']

                if not snp["syn"].startswith("N"):
                    continue
                    
                for k in xrange(SAMPLE_SIZE):
                    _cov_site = full_site_cov[k]
                    _cov_snp = full_snp_cov[k]
                    
                    cov_counter = gene2pos2cov.setdefault((k, site["c"], site["g"]), {}).setdefault((site['pos'], _cov_site), [])
                    
                    if args.atype == "nt":
                        ref_symbol = site["ref"]
                        snp_symbol = snp["snp"]
                    elif args.atype == "pa":
                        ref_symbol = "0"
                        snp_symbol = "1"
                        
                    if _cov_site == -1:
                        dpos[k] = "-"
                        dcov[k] = -1
                    elif _cov_site < MIN_SITE_COV:
                        dpos[k] = "-" # if there is not enough coverage, we don't know if snp or ref site
                        dcov[k] = _cov_snp
                    elif _cov_snp < MIN_SNP_COV:
                        dpos[k] = ref_symbol
                        dcov[k] = _cov_snp
                    elif dcov[k] is None:
                        dpos[k] = snp_symbol
                        dcov[k] = _cov_snp
                        data_array[k] = 1
                        cov_counter.append((_cov_snp, snp["snp"], 'snp'))
                    elif dpos[k] == snp_symbol:
                        print
                        print (k, site["c"], site["g"], site['pos'], snp['snp']), dpos[k]
                        raw_input('')
                    elif _cov_snp > dcov[k]:
                        dpos[k] = snp_symbol
                        dcov[k] = _cov_snp
                        cov_counter.append((_cov_snp, snp["snp"], 'snp'))
                    elif _cov_snp <= dcov[k]:
                        cov_counter.append((_cov_snp, snp["snp"], 'snp'))
                        pass # ignore less frequent snp
                    else:
                        print
                        print 'ref', site['ref'], _cov_site
                        print 'snp', snp['snp'], _cov_snp
                        print 'current_cov', dcov[k]
                        print 'current_sit', dpos[k]
                        print MIN_SITE_COV
                        print MIN_SNP_COV
                        raise ValueError('This should never occur!')

                for k in xrange(SAMPLE_SIZE):
                    _cov_site = full_site_cov[k]
                    if _cov_site == -1:
                        continue
                   
                    cov_counter = gene2pos2cov.setdefault((k, site["c"], site["g"]), {}).setdefault((site['pos'], _cov_site), [])
                    total_snp_cov = sum(cov[0] for cov in cov_counter)
                    #print total_snp_cov, _cov_site
                    if _cov_site != total_snp_cov:
                        cov_counter.append((_cov_site - total_snp_cov, site["ref"], 'ref'))
                
            ndata = sum(data_array)
            if ndata >= MIN_COL_SNP:
                positions.append([''.join(dpos), ndata])
                bygene_positions.setdefault((site["c"], site["g"]), []).append(''.join(dpos))

        if positions:
            print 'dumping', len(positions), 'informative columns'
            OUT = open('%s.%s.fa' %(OUTFILE, sp), 'w')
            INFO = open('%s.%s.fa.info' %(OUTFILE, sp), 'w')
            column_data = [ndata for p, ndata in positions]
            print >>INFO, "# MIN_SITE_COV: ", MIN_SITE_COV
            print >>INFO, "# MIN_SNP_COV: ", MIN_SNP_COV
            print >>INFO, "# MIN_COL_SNP: ", MIN_COL_SNP
            print >>INFO, "# nsites: ", len(column_data)
            print >>INFO, "# avg snps per column: ", numpy.mean(column_data)
            print >>INFO, "# std deviation: ", numpy.std(column_data)
            print >>INFO, "# Num samples: ", len(species_samples)
            print >>INFO, "# samples: ", ' '.join(species_samples)

            # for i, sname in enumerate(SAMPLE_NAMES):
            #     sample_seq = ''.join([p[i] for p, ndata in positions])
            #     print >>OUT, '>%s\n%s' %(sname, sample_seq)
            for i, sname in enumerate(SAMPLE_NAMES):
                sample_seq = ""
                for gene, snppos in bygene_positions.iteritems():
                    sample_seq += ''.join([p[i] for p in snppos]) + " "
                if args.remove_all_gaps and not (set(sample_seq) - set("- ")):
                    continue
                print >>OUT, '>%s\n%s' %(sname, sample_seq)
            
            OUT.close()
            INFO.close()
            biga, bigb = [], []
            binom_genes = defaultdict(int)
            binom_samples = defaultdict(int)
            
            for gene in gene2pos2cov:
                #print SAMPLE_NAMES[gene[0]], gene
                show = False
                a, b = [], []
                for (pos, poscov), values in gene2pos2cov[gene].iteritems():
                    if poscov > 0:
                        #print pos, poscov 
                        values.sort()
                        #pprint(values)
                        if len(values) == 2:
                            a.append(values[0][0])
                            b.append(values[1][0])
                            
                if a and b:
                    biga.extend(a)
                    bigb.extend(b)
                    d, pvalue = stats.ks_2samp(a, b) 
                    if d==1 and pvalue <= 0.00001:
                        print "GENE", gene, SAMPLE_NAMES[gene[0]]
                        gene_info = gene_db.find({'sp':sp, 'n':gene[2]}, {"nt":1, "c":1, "s":1, "e":1}).next()
                        contig_seq = contig_db.find({'sp':sp, 'c':gene_info["c"]}, {"nt":1}).next()
                        print gene
                        print gene_info["nt"]

                        gstart, gend = gene_info["s"], gene_info["e"]
                        ref_seq = contig_seq["nt"][gstart-1:gend]
                        var_seq = ref_seq
                       
                        for (pos, poscov), values in sorted(gene2pos2cov[gene].items()):
                            pass
                            print "   POS: %06d %06d    REF:%s   %s " %(pos, pos-gstart, ref_seq[pos-gstart], values)
                        print "KS test: ", d, pvalue, a, b

                        # extract reference gene sfgequence
                        binom_genes[gene[2]] += 1
                        binom_samples[SAMPLE_NAMES[gene[0]]] +=1 
                        raw_input('enter to continue')

            print len(binom_genes), 'bimodel genes',  sorted(binom_genes.items(), key=lambda x: x[1], reverse=True)[0:10]
            print len(binom_samples), 'bimodel samples',  sorted(binom_samples.items(), key=lambda x: x[1], reverse=True)[0:10]
            print stats.ks_2samp(biga, bigb)
            print "--"
            
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", choices=["full", "filtered"], dest="db", required=True)
    parser.add_argument("-o", dest="outfile", required=True)
    parser.add_argument("--min_site_cov", dest="min_site_cov", default=5, type=int)
    parser.add_argument("--min_snp_cov", dest="min_snp_cov", default=2, type=int)
    parser.add_argument("--min_snp_samples", dest="min_snp_samples", default=1, type=int)
    parser.add_argument("--target_species", dest="target_species", nargs="*", type=int)
    parser.add_argument("--atype", dest='atype', choices=["pa", "nt"])
    parser.add_argument("--remove_all_gaps", dest='remove_all_gaps', action="store_true")
    parser.add_argument("--genes_only", dest='genes_only', action="store_true")
    parser.add_argument("--test", dest='testcase', type=int)

    args = parser.parse_args()

    MIN_SITE_COV = args.min_site_cov
    MIN_SNP_COV = args.min_snp_cov
    MIN_COL_SNP = args.min_snp_samples
    FILTERED_DB = args.db == "filtered"
    OUTFILE = args.outfile

    TARGET_SPECIES = set(args.target_species) if args.target_species else None
    print TARGET_SPECIES

    _NN = {
        'A': 'A',
        'C': 'C',
        'T': 'T',
        'G': 'G',
        'N': 'N',
        'AG' : 'R',
        'CT' : 'Y',
        'GC' : 'S',
        'AT' : 'W',
        'GT' : 'K',
        'AC' : 'M',
        'CGT' : 'B',
        'AGT' : 'D',
        'ACT' : 'H',
        'ACG' : 'V',
        'ACTG' : 'N'}

    NN = dict([(tuple(sorted(list(k))), v) for k,v in _NN.iteritems()])
    NN2 = dict([(v, k) for k,v in NN.iteritems()])

    SAMPLE_NAMES = [line.strip() for line in open('1086_sample_names')]
    SAMPLE_INDEX = dict([(name, i) for i, name in enumerate(SAMPLE_NAMES)])

    mongo = MongoClient()
    if FILTERED_DB:
        snp_db =  mongo.metagenvar_filtered.snp
        sites_db =  mongo.metagenvar_filtered.sites
        samples_db = mongo.metagenvar_filtered.samples
    else:
        snp_db =  mongo.metagenvar.snp
        sites_db =  mongo.metagenvar.sites
        samples_db = mongo.metagenvar.samples

    gene_db = mongo.ref9.genes
    contig_db = mongo.ref9.contigs
            
    SP_INFO_CACHE = 'species_%s.pkl' %args.db 

    try:
        sorted_sp = cPickle.load(open(SP_INFO_CACHE))
        print 'loaded from file'
    except Exception:
        species =  snp_db.distinct('sp')
        sorted_sp = []
        for sp in species:
            nsnps = snp_db.find({'sp':sp}, {'c':1}).count()
            sorted_sp.append([nsnps, sp])
        sorted_sp.sort()
        cPickle.dump(sorted_sp, open(SP_INFO_CACHE, 'w'))
    print sorted_sp
    SAMPLE_SIZE = len(SAMPLE_NAMES)

    dump_snp_matrix(sorted_sp)

