import sys
from string import strip
import json

def as_json(d):
    sep = (',', ':')
    return json.dumps(d, separators=sep)

SNP = open('snp.json', 'w')
SITE = open('sites.json', 'w')
for fname in sys.argv[1:]:
    print >>sys.stderr, fname
    for line in open(fname):
        contig, gene, pos, base, cov_per_sample, snp_data = map(strip, line.split('\t'))
        pos = int(pos)
        taxa, contig_name = contig.split('.', 1)
        if gene != '-':
            taxa2, gene_name = gene.split('.', 1)
        else:
            taxa2 = taxa
            gene_name = ''
        taxa, taxa2 = int(taxa), int(taxa2)
        if taxa != taxa2:
            raise ValueError('taxa names are not the same')
        coverage = map(int, cov_per_sample.split('|'))

        for nsnp, snp_track in enumerate(snp_data.split(',')):
            snp_track = snp_track.split('|')
            nreads, snp, syn, snp_coverage =  int(snp_track[0]), snp_track[1], snp_track[2], map(int, snp_track[3:])

            if len(snp_coverage) != len(coverage):
                raise ValueError("coverages are not the same size")
            snp_entry = {"nr":nreads, "sp":taxa, "snp":snp, "pos":pos, "g":gene_name, "c":contig_name, "cov":snp_coverage, "syn":syn}
            print >>SNP, as_json(snp_entry)
        print >>SITE, as_json({"pos":pos, "sp":taxa, "g":gene_name, "c":contig_name, "cov":coverage, "ref":base, 'nsnp':nsnp+1})
SNP.close()
SITE.close()
