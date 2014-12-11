from pprint import pprint
import pymongo
import numpy
import cPickle
import sys
import argparse

from pymongo import MongoClient

def get_global_coverage(cov, header_names, sample_indexes):
    full_cov = [-1] * len(sample_indexes)
    for i, sname in enumerate(header_names):
        full_cov[sample_indexes[sname]] = cov[i]
    return full_cov

parser = argparse.ArgumentParser()
parser.add_argument("--db", choices=["full", "filtered"], dest="db", required=True)
parser.add_argument("-o", dest="outfile", required=True)
parser.add_argument("--min_site_cov", dest="min_site_cov", default=5, type=int)
parser.add_argument("--min_snp_cov", dest="min_snp_cov", default=2, type=int)
parser.add_argument("--min_snp_samples", dest="min_snp_samples", default=1, type=int)
parser.add_argument("--target_species", dest="target_species", nargs="*", type=int)

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
    
sample_size = len(SAMPLE_NAMES)
pos2info = {}
# iter over the list of species sorted by number of pos with snps
for percent, (nsnps, sp) in enumerate(reversed(sorted_sp)):
    if TARGET_SPECIES and sp not in TARGET_SPECIES:
        continue
        
    #nsnps = snp.find({'sp':sp}, {'c':1}).count()
    print 'processing', sp, 'with', nsnps, 'snps. ', percent+1, "/", len(sorted_sp)

    positions = [ ]
    for i, site in enumerate(sites_db.find({'sp':sp}, {'c':1, 'pos':1, 'cov':1, 'ref':1})):
        if i % 1000 == 0:
            #if i == 5000: break
            print '% 6d %d\r' %(i, len(positions)),  
            sys.stdout.flush()
            
        if FILTERED_DB:
            species_samples = samples_db.find_one({'sp':sp}, {'samples':1})["samples"]
            full_site_cov = get_global_coverage(site['cov'], species_samples, SAMPLE_INDEX)
        else:
            species_samples = SAMPLE_NAMES
            full_site_cov = site['cov']
            
        dpos = [None] * sample_size
        dcov = [None] * sample_size
        data_array = [0] * sample_size
        for j, snp in enumerate(snp_db.find({'sp':sp, 'c':site['c'], 'pos':site['pos']}, {'snp':1, 'cov':1})):
            if FILTERED_DB:
                full_snp_cov = get_global_coverage(snp['cov'], species_samples, SAMPLE_INDEX)
            else:
                full_snp_cov = snp['cov']
                
            for k in xrange(sample_size):
                _cov_site = full_site_cov[k]
                _cov_snp = full_snp_cov[k]
                if _cov_site == -1:
                    dpos[k] = "-"
                    dcov[k] = -1
                elif _cov_site < MIN_SITE_COV:
                    dpos[k] = str(site['ref'])
                    dcov[k] = _cov_snp
                elif _cov_snp < MIN_SNP_COV:
                    dpos[k] = str(site['ref'])
                    dcov[k] = _cov_snp
                elif dcov[k] is None:
                    dpos[k] = snp['snp']
                    dcov[k] = _cov_snp
                    data_array[k] = 1
                elif _cov_snp > dcov[k]:
                    dpos[k] = snp['snp']
                    dcov[k] = _cov_snp
                elif _cov_snp <= dcov[k]:
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
                    
        ndata = sum(data_array)
        if ndata >= MIN_COL_SNP:
            positions.append([''.join(dpos), ndata])

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
        
        for i, sname in enumerate(SAMPLE_NAMES):
            sample_seq = ''.join([p[i] for p, ndata in positions])
            print >>OUT, '>%s\n%s' %(sname, sample_seq)
        OUT.close()
        INFO.close()
        


