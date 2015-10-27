import sys
from Bio import SeqIO

with open(sys.argv[1]) as f:
    for rec in SeqIO.parse(f, 'fasta'):
        with open("chr%s.fa" % rec.id, 'w') as o:
            print "Writing chr%s" % rec.id
            o.write(">chr%s\n%s\n" % (rec.id, rec.seq))

