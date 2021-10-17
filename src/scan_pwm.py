# ./scan_pwm.py <matrix.jaspar> <genbank_mRNA_ID> <pseudo-count> <scorethreshold>

from utils import *
from pwm import *

seq_length = 100

if __name__ == "__main__":

	accession = sys.argv[1]
	path = sys.argv[2]
	psw = float(sys.argv[3])
	threshold = float(sys.argv[4])
	
	pssm = pwm2pssm(path,psw)

	id_gene = mrna_to_gene(accession)
	print("Id gene is %s ",id_gene)
	
	seq_start , seq_stop = get_seq_start_stop(id_gene)
	gen_accession = get_genom_id_accession(id_gene)

	seq = upstream_gene_seq(gen_accession, seq_length , seq_start, seq_stop)

	print(pssm)
	print(seq)
	print(scan_sequence(pssm,seq,threshold))