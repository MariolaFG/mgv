from pprint import pprint
import pymongo
import numpy
import cPickle
import sys

MIN_SITE_COV = 10
MIN_SNP_COV = 5
MIN_COL_SNP = 2

OUTFILE = sys.argv[1]

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

from pymongo import MongoClient
mongo = MongoClient()

snp =  mongo.metagenvar.snp
sites =  mongo.metagenvar.sites

try:
    sorted_sp = cPickle.load(open('species.pkl'))
    print 'loaded from file'
except Exception:
    species =  snp.distinct('sp')
    sorted_sp = []
    for sp in species:
        nsnps = snp.find({'sp':sp}, {'c':1}).count()
        sorted_sp.append([nsnps, sp])
    sorted_sp.sort()
    cPickle.dump(sorted_sp, open('species.pkl', 'w'))
    
sample_size = 1086
pos2info = {}
# iter over the list of species sorted by number of pos with snps
for percent, (nsnps, sp) in enumerate(reversed(sorted_sp)):
    #nsnps = snp.find({'sp':sp}, {'c':1}).count()
    print 'processing', sp, 'with', nsnps, 'snps. ', percent+1, "/", len(sorted_sp)

    positions = []
    for i, site in enumerate(sites.find({'sp':sp}, {'c':1, 'pos':1, 'cov':1, 'ref':1})):
        if i % 1000 == 0:
            #if i == 5000: break
            print '% 6d %d\r' %(i, len(positions)),  
            sys.stdout.flush()
        dpos = [str(site['ref'])] * sample_size
        dcov = [0] * sample_size
        ndata = 0
        for j, snp_site in enumerate(snp.find({'sp':sp, 'c':site['c'], 'pos':site['pos']}, {'snp':1, 'cov':1})):
            for k in xrange(sample_size):
                if site['cov'][k] >= MIN_SITE_COV and snp_site['cov'][k] >= MIN_SNP_COV and snp_site['cov'][k] >= dcov[k]:
                    dcov[k] = snp_site['cov'][k]
                    dpos[k] = snp_site['snp']
                    ndata += 1
        if ndata:
            positions.append([''.join(dpos), ndata])

    if positions:
        print 'dumping', len(positions), 'informative columns'
        OUT = open('%s.%s.fa' %(OUTFILE, sp), 'w')
        column_data = [ndata for p, ndata in positions if ndata >= MIN_COL_SNP]
        print >>OUT, "# MIN_COL_SNP: ", MIN_COL_SNP
        print >>OUT, "# MIN_SITE_COV: ", MIN_SITE_COV
        print >>OUT, "# MIN_SNP_COV: ", MIN_SNP_COV
        print >>OUT, "# nsites: ", len(column_data)
        print >>OUT, "# avg snps per column: ", numpy.mean(column_data)
        print >>OUT, "# std deviation: ", numpy.std(column_data)
        for i, sname in enumerate(SAMPLE_NAMES):
            sample_seq = ''.join([p[i] for p, ndata in positions if ndata >= MIN_COL_SNP])
            print >>OUT, '>%s\n%s' %(sname, sample_seq)
        OUT.close()        
        


