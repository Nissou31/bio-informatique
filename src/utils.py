from Bio import SeqIO
from Bio import Entrez
#import time
#from tqdm import tqdm
import os


dbfrom="nucleotide"
db="gene"


"""
This was an idea to make the algorithm work faster by not calculating all the points each time, and use the previous window data.
NOT USED.

class Window:
	#A class used to store data about previous calculated windows.

	#The goal is to reduce time of calculating score for a window, by storing 
	#data of the current window , when calculating a score for the next one,
	#jumpt into position that are greater then the end of the previous window,
	#and then add the previously calculated score.

	#Exemple :
		
	#	We have a window starting at position 0 and ending at 40,
	#	first we calculate a score for this one ( By the score_window function),
	#	then we store the start, end and score for this window into a the global class
	#	window. Next, when scanning the next window at position , by instence, 15 and ending at 55,
	#	all the score between 15 and 40 already being calculated, so we only work on 41 to 55,
	#	and add to them the previous data.

	def __init__(self,start,end):
		self.start = start
		self.end = end
		self._data_window_ = dict()

	def set_start(self,start):
		self.start = start
	
	def set_end(self,end):
		self.end = end

	def get_start(self):
		return self.start

	def get_end(self):
		return self.end

	def store_data_window(self,cp,np,pscore):
		pos_couple = (cp,np)
		self._data_window_ = dict()
		self._data_window_[pos_couple] = pscore 

	def get_data(self):
		return self._data_window_
"""


def find_cds(seqrecord):
	"""Fonction used to find the coding area on a gene.
	
	Returns the start and end location of the coding area on a gene.
	"""
	result = []
	couplePosition = tuple()

	for f in seqrecord.features :
		if f.type == "CDS":
			couplePosition = (f.location._start.position,f.location._end.position)
			result.append(couplePosition)

	return result


def mrna_to_gene(num_accession):
	"""Search for the id of gene given the number of accession of a mrna.
	
	Returns the id of the gene from the gene data base.
	"""
	
	handle = Entrez.elink(dbfrom=dbfrom,db=db, id=num_accession)
	record = Entrez.read(handle)
	handle.close()
	
	try:
		linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
		return linked[0]
	except ValueError:
		print("No id found...")

def get_genom_id_accession(id_gene):
	"""Search for the id accession fo the chromosom given a gene id.

	Returns the id accession of chromosomes.
	"""
	
	handle = Entrez.esummary(db=db, id=id_gene)
	record = Entrez.read(handle)
	handle.close()

	dictElem = dict(record)['DocumentSummarySet']
	premierElem = dict(dictElem)['DocumentSummary']
	info = []

	for e in premierElem:
		for stuff in dict(e).items():
			if stuff[0] == 'GenomicInfo':
				info = stuff[1]

	return info[0]['ChrAccVer']


def get_seq_start_stop(id_gene):
	"""Search for the start and end sequence from a gene.

	Given an id of a gene, search for the start and end position of it sequence on the gene.
	"""
	seq_start = 0
	seq_stop = 0

	handle = Entrez.esummary(db=db, id=id_gene)
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

	return seq_start, seq_stop


def upstream_gene_seq(id_gene, seq_length):
	"""Find the upstream sequence of a gene.
	Return the promotor sequence of a gene.
    """
	
	seq_start, seq_stop = get_seq_start_stop(id_gene)
	genom_accession = get_genom_id_accession(id_gene)
	
	if seq_start < seq_stop:
		handle = Entrez.efetch(db=dbfrom, id=genom_accession, rettype='fasta', strand=1,retmode='text', seq_start=str(int(seq_start)-seq_length), seq_stop=seq_start)
		genomic_record = SeqIO.read(handle, "fasta")
	else :
		handle = Entrez.efetch(db=dbfrom, id=genom_accession, rettype='fasta', strand=-1,retmode='text', seq_start=str(int(seq_start)+2), seq_stop=str(int(seq_start)+seq_length))
		genomic_record = SeqIO.read(handle, "fasta")

	
	handle.close()


	return genomic_record


"""
def download_promotors_tqdm(ids_mrna_list,seq_length,out_put_dir):

	print("Downloading files, please wait :")
	for id_mrna in ids_mrna_list:
		print("Fetching data for gene : "+id_mrna)
		for i in tqdm(range(100)):
			time.sleep(0.009)
		id_gene = mrna_to_gene(id_mrna)
		seq = upstream_gene_seq(id_gene, seq_length)
		filename = os.path.join(out_put_dir, id_mrna+"_"+str(seq_length)+".fa")
		SeqIO.write(seq, filename, 'fasta')
	print("Download done.")
"""


def download_promotors(ids_mrna_list,seq_length,out_put_dir):
	"""Downlaod promotors sequences for a list of MRNA as a fasta files.
	"""

	print("Downloading files, please wait :")
	for id_mrna in ids_mrna_list:
		print("Fetching data for gene : "+id_mrna)
		id_gene = mrna_to_gene(id_mrna)
		seq = upstream_gene_seq(id_gene, seq_length)
		filename = os.path.join(out_put_dir, id_mrna+"_"+str(seq_length)+".fa")
		SeqIO.write(seq, filename, 'fasta')
	print("Download done.")


def Seq_obj_from_files(files_list):
	"""Take a list of promotor sequences from fasta files.

	Return a list of SeqRecord objects.
	"""
	
	seq_obj_list = list()
	
	for file in files_list:
		seq_object = SeqIO.read(file,"fasta")
		seq_obj_list.append(seq_object)

	return seq_obj_list
	
"""
if __name__ == "__main__":


	

	#Entrez.email = "Your.Name.Here@example.org"
	
	
	recordGB = SeqIO.read("../data/sequence-2.gb", "genbank")

	
	handle = Entrez.efetch(db="nucleotide", id="NM_007389", rettype="gb", retmode="text")

	recordFetch = SeqIO.read(handle, "genbank")

	string_seq_gb = str(recordGB.seq)
	string_seq_fetch = str(recordFetch.seq)
	if string_seq_gb == string_seq_fetch :
		print("The sequences in both records are identical.")


	localres = find_cds(recordGB)
	fetchres = find_cds(recordFetch)
	if localres == fetchres :
		print("The CDS location in both records are identical.")

	
	id_gene = mrna_to_gene("NM_007389")
	print("Id gene for NM_007389 is %s"  %id_gene)
	
	seq_start , seq_stop = get_seq_start_stop(id_gene)
	gen_accession = get_genom_id_accession(id_gene)

	print("Id accession is %s"  %gen_accession)

	seq = upstream_gene_seq(id_gene, 100)
	print(seq)
	
	handle.close()
	
	
	id_gene = mrna_to_gene("NM_001267550")
	print("Id gene for NM_001267550 is %s"  %id_gene)
	
	seq_start , seq_stop = get_seq_start_stop(id_gene)
	gen_accession = get_genom_id_accession(id_gene)
	print((seq_start,seq_stop), gen_accession)
	print(upstream_gene_seq(id_gene, 1000).seq)

	#handle.close()

	
	#list_mRNA = ["NM_007389", "NM_079420", "NM_001267550", "NM_002470", "NM_003279", "NM_005159", "NM_003281", "NM_002469", "NM_004997", "NM_004320", "NM_001100", "NM_006757"]
	
	id_gene = mrna_to_gene("NM_007389")
	print("Id gene for NM_007389 is %s"  %id_gene)
	
	seq_start , seq_stop = get_seq_start_stop(id_gene)
	gen_accession = get_genom_id_accession(id_gene)
	print((seq_start,seq_stop), gen_accession)
	print(upstream_gene_seq(id_gene, 1000).seq)

	handle.close()



	list_mRNA = ["NM_007389", "NM_079420", "NM_001267550"]
	download_promotors(list_mRNA,1000,"../data/")
"""

	
