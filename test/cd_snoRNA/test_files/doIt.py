import os
import sys
import pandas as pd
from Bio import SeqIO

unm = pd.read_table(sys.argv[1], header=None)[0].tolist()
print unm[:10]
seqs = []
with open(os.path.expanduser("~/snoRNAs/hybridSearch/Disa/data/Disa_426-4_fib_hek/unmapped_sequences.fa")) as f:
    for rec in SeqIO.parse(f, 'fasta'):
        seqs.append([str(rec.id), str(rec.seq)])
seqs[:5]
df_seqs = pd.DataFrame(seqs)
df_seqs.head()
df_seqs.columns = ["id", "seq"]
df_seqs.head()
df_seqs.index = df_seqs.id
df_seqs.head()
df_seqs.ix[unm].head()
ndf = df_seqs.ix[unm]
ndf.head()
with open(sys.argv[1] + ".fa", 'w') as out:
    for idx, row in ndf.iterrows():
        out.write(">%s\n%s\n" % (row.id, row.seq))

