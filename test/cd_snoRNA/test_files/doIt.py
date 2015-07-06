unm = pd.read_table("rna18s_unmapped.tab", header=None)
unm.head()
unm = pd.read_table("rna18s_unmapped.tab", header=None)[0].tolist()
unm[:10]
seqs = []
with open("/import/bc2/home/zavolan/gumiennr/snoRNAs/hybridSearch/Disa/data/Disa_426-4_fib_hek/unmapped_sequences.fa") as f:
    for rec in SeqIO.parse(f, 'fasta'):
        seqs.append([str(rec.id), str(rec.seq)])
seqs[:5]
df_seqs = pd.DataFrame(seqs)
df_seqs.head()
df_seqs.columns = ["id", "seq"]
df_seqs.head()
df.index = df.id
df_seqs.index = df_seqs.id
df_seqs.head()
df_seqs.ix[unm].head()
ndf = df_seqs.ix[unm]
ndf.head()
with open("rna18s_unmapped.fa", 'w') as out:
    for idx, row in ndf.iterrows():
        out.write(">%s\n%s\n" % (row.id, row.seq))

