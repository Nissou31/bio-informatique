from Bio import SeqIO
from Bio import Entrez


"""
recordGB = SeqIO.read("data/sequence-2.gb", "genbank")
recordFASTA = SeqIO.read("data/sequence.fasta", "fasta")

#print(recordFASTA)

string_seq_gb = str(recordGB.seq)
string_seq_fasta = str(recordFASTA.seq)



#print("Sequence alphabet %s" % string_seq_fasta)
#print("Sequence alphabet %s" % string_seq_gb)

if string_seq_gb == string_seq_fasta :
	print("The sequences in both format are identical.")

# recordFASTA_rc = recordFASTA.reverse_complement(id=recordFASTA.id + "_rc")


# recordGB_rc = recordGB.seq.reverse_complement()


# print("Reverse sequence alphabet %s" % str(recordGB_rc))


#for f in recordGB.features:
#	print(f)

# first_feature = recordGB.features[3]

# print(first_feature)





handle = Entrez.elink(dbfrom="nuccore",db="gene", id="NM_007389")
record = Entrez.read(handle)
handle.close()
print(record)
linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
print(linked)



handle = Entrez.esummary(db="gene", id="11435")
record = Entrez.read(handle)
handle.close()

dictElem = dict(record)['DocumentSummarySet']
premierElem = dict(dictElem)['DocumentSummary']
info = []

for e in premierElem:
	for stuff in dict(e).items():
		print(stuff)
		if stuff[0] == 'GenomicInfo':
			info = stuff[1]

print(info)
print(info[0]['ChrStart'])




handle = Entrez.efetch(db='nucleotide', id='NC_000068.8', rettype='fasta', strand=-1,retmode='text', seq_start=73393625, seq_stop=73410682)
record = handle.readline().strip()
genomic_record = SeqIO.read(handle, "fasta")

handle.close()

print(genomic_record.seq)


from Bio import motifs

with open("data/motifs/MA0114.4.jaspar") as handle:
    m = motifs.read(handle, "jaspar")
 	
print(m.counts.normalize(0.1))
pssm = motifs.matrix.PositionWeightMatrix.log_odds(m.counts.normalize(0.1))

print(pssm)
print(pssm[0][0])
"""

Entrez.email = "aness306@gmail.com"

"""
handle = Entrez.elink(dbfrom="nuccore",db="gene", id="NM_007389")
record = Entrez.read(handle)
handle.close()
print(record)
linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
print(linked)


"""
handle = Entrez.esummary(db="gene", id="11435")
record = Entrez.read(handle)
handle.close()
dictElem = dict(record)['DocumentSummarySet']
premierElem = dict(dictElem)['DocumentSummary']
info = []

for e in premierElem:
	for stuff in dict(e).items():
		if stuff[0] == 'GenomicInfo':
			info = stuff[1]



seq_start=info[0]['ChrStart']
seq_stop= info[0]['ChrStop']

if int(seq_start) < int(seq_stop):
	handle = Entrez.efetch(db='nucleotide', id='NC_000068.8', rettype='fasta', strand=1,retmode='text', seq_start=seq_start, seq_stop=seq_stop)
else:
	handle = Entrez.efetch(db='nucleotide', id='NC_000068.8', rettype='fasta', strand=2,retmode='text', seq_start=seq_start, seq_stop=seq_stop)

genomic_record = SeqIO.read(handle, "fasta")
SeqIO.write(genomic_record, "test_1000.gb", "fasta")
handle.close()


def find_cds(seqrecord):
	
	result = ""
	

	for f in seqrecord.features :
		if f.type == "gene":
			resulat= f.gene

	return result



print(find_cds(genomic_record))

