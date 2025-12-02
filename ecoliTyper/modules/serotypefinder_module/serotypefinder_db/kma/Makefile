CFLAGS ?= -Wall -O3
CFLAGS += -std=c99
LIBS = align.o alnfrags.o ankers.o assembly.o chain.o cmp.o compdna.o compkmers.o compress.o conclave.o db.o decon.o dist.o ef.o filebuff.o frags.o hashmap.o hashmapcci.o hashmapkma.o hashmapkmers.o hashtable.o index.o kma.o kmapipe.o kmeranker.o kmers.o kmmap.o loadupdate.o makeindex.o matrix.o merge.o middlelayer.o mt1.o nw.o pherror.o printconsensus.o qc.o qseqs.o qualcheck.o reassign.o runinput.o runkma.o sam.o savekmers.o seq2fasta.o seqmenttree.o seqparse.o shm.o sparse.o spltdb.o stdnuc.o stdstat.o trim.o tmp.o tsv.o update.o updateindex.o updatescores.o valueshash.o vcf.o xml.o
PROGS = kma kma_index kma_shm kma_update

.c .o:
	$(CC) $(CFLAGS) -c -o $@ $<

all: $(PROGS)

kma: main.c libkma.a
	$(CC) $(CFLAGS) -o $@ main.c libkma.a -lm -lpthread -lz $(LDFLAGS)

kma_index: kma_index.c libkma.a
	$(CC) $(CFLAGS) -o $@ kma_index.c libkma.a -lm -lz $(LDFLAGS)

kma_shm: kma_shm.c libkma.a
	$(CC) $(CFLAGS) -o $@ kma_shm.c libkma.a $(LDFLAGS)

kma_update: kma_update.c libkma.a
	$(CC) $(CFLAGS) -o $@ kma_update.c libkma.a $(LDFLAGS)

libkma.a: $(LIBS)
	$(AR) -csr $@ $(LIBS)

clean:
	$(RM) $(LIBS) $(PROGS) libkma.a

align.o: align.h chain.h compdna.h hashmapcci.h nw.h stdnuc.h stdstat.h
alnfrags.o: alnfrags.h align.h ankers.h compdna.h hashmapcci.h qseqs.h threader.h updatescores.h
ankers.o: ankers.h compdna.h pherror.h qseqs.h
assembly.o: assembly.h align.h filebuff.h hashmapcci.h kmapipe.h pherror.h stdnuc.h stdstat.h threader.h
chain.o: chain.h penalties.h pherror.h stdstat.h
cmp.o: cmp.h hashmapkma.h kmmap.h pherror.h tmp.h version.h
compdna.o: compdna.h pherror.h stdnuc.h
compkmers.o: compkmers.h pherror.h
compress.o: compress.h hashmap.h hashmapkma.h pherror.h valueshash.h
conclave.o: conclave.h frags.h pherror.h qseqs.h stdnuc.h
db.o: db.h hashmapkma.h pherror.h stdstat.h
decon.o: decon.h compdna.h filebuff.h hashmapkma.h seqparse.h stdnuc.h qseqs.h updateindex.h
dist.o: dist.h hashmapkma.h matrix.h pherror.h
ef.o: ef.h assembly.h stdnuc.h vcf.h version.h
filebuff.o: filebuff.h pherror.h qseqs.h
frags.o: frags.h filebuff.h pherror.h qseqs.h
hashmap.o: hashmap.h hashtable.h pherror.h stdstat.h
hashmapcci.o: hashmapcci.h pherror.h stdnuc.h stdstat.h
hashmapkma.o: hashmapkma.h pherror.h stdstat.h
hashmapkmers.o: hashmapkmers.h pherror.h
hashtable.o: hashtable.h hashmapkma.h hashmapkmers.h pherror.h
index.o: index.h compress.h decon.h hashmap.h hashmapkma.h loadupdate.h makeindex.h pherror.h stdstat.h version.h
kma.o: kma.h ankers.h assembly.h chain.h hashmapkma.h kmers.h mt1.h 
penalties.h pherror.h qc.h qseqs.h runinput.h runkma.h savekmers.h sparse.h spltdb.h tmp.h version.h kmapipe.o: kmapipe.h pherror.h
kmeranker.o: kmeranker.h penalties.h
kmers.o: kmers.h ankers.h compdna.h hashmapkma.h kmapipe.h pherror.h qseqs.h savekmers.h spltdb.h
kmmap.o: kmmap.h hashmapkma.h
loadupdate.o: loadupdate.h pherror.h hashmap.h hashmapkma.h updateindex.h
makeindex.o: makeindex.h compdna.h filebuff.h hashmap.h pherror.h qseqs.h seqparse.h updateindex.h
matrix.o: matrix.h pherror.h
merge.o: merge.h hashmapkma.h kmmap.h middlelayer.h pherror.h stdstat.h tmp.h
middlelayer.o: middlelayer.h hashmapkma.h pherror.h
mt1.o: mt1.h assembly.h chain.h filebuff.h hashmapcci.h kmapipe.h nw.h penalties.h pherror.h printconsensus.h qseqs.h runkma.h stdstat.h tsv.h vcf.h
nw.o: nw.h pherror.h stdnuc.h penalties.h
pherror.o: pherror.h
printconsensus.o: printconsensus.h assembly.h
qc.o: qc.h pherror.h
qseqs.o: qseqs.h pherror.h
qualcheck.o: qualcheck.h compdna.h hashmap.h pherror.h stdnuc.h stdstat.h
reassign.o: reassign.h assembly.h compdna.h hashmapkma.h kmmap.h pherror.h qseqs.h runkma.h stdnuc.h
runinput.o: runinput.h compdna.h filebuff.h pherror.h qseqs.h seqparse.h
runkma.o: runkma.h align.h alnfrags.h assembly.h chain.h compdna.h ef.h filebuff.h frags.h hashmapcci.h kmapipe.h nw.h pherror.h printconsensus.h qseqs.h reassign.h stdnuc.h stdstat.h tmp.h tsv.h vcf.h
sam.o: sam.h nw.h pherror.h qseqs.h runkma.h
savekmers.o: savekmers.h ankers.h compdna.h hashmapkma.h kmeranker.h penalties.h pherror.h qseqs.h stdnuc.h stdstat.h threader.h
seq2fasta.o: seq2fasta.h pherror.h qseqs.h runkma.h stdnuc.h
seqmenttree.o: seqmenttree.h pherror.h
seqparse.o: seqparse.h filebuff.h qseqs.h
shm.o: shm.h pherror.h hashmapkma.h version.h
sparse.o: sparse.h compkmers.h hashtable.h kmapipe.h pherror.h qseqs.h qc.h runinput.h savekmers.h stdnuc.h stdstat.h
spltdb.o: spltdb.h align.h alnfrags.h assembly.h chain.h compdna.h ef.h filebuff.h frags.h hashmapcci.h kmapipe.h nw.h pherror.h printconsensus.h qseqs.h runkma.h stdnuc.h stdstat.h tmp.h tsv.h vcf.h
stdnuc.o: stdnuc.h
stdstat.o: stdstat.h
trim.o: trim.h compdna.h pherror.h runinput.h qc.h qseqs.h
tmp.o: tmp.h pherror.h threader.h
tsv.o: tsv.h assembly.h
update.o: update.h hashmapkma.h pherror.h stdnuc.h
updateindex.o: updateindex.h compdna.h hashmap.h hashmapcci.h pherror.h qualcheck.h stdnuc.h pherror.h
updatescores.o: updatescores.h qseqs.h
valueshash.o: valueshash.h pherror.h
vcf.o: vcf.h assembly.h filebuff.h stdnuc.h stdstat.h version.h
xml.o: xml.h pherror.h version.h
