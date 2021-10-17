from Bio import motifs
from utils import *


Entrez.email = "Your.Name.Here@example.org"


import sys
import glob


def pwm2pssm(FPM,psw):
	"""Return TF name and his PositionSpecificScoringMatrix object from a .jaspar file and 
	a given PositionWeigthMAtrix object.
	"""
	
	FPM+=".jaspar" # Add .jaspar the jaspar filename.
	path=os.path.join("../data/motifs",FPM) # Create a valid path to the target file.

	with open(path) as handle:
		m = motifs.read(handle, "jaspar") 
	
	# Get the PSSM from PSW.
	pssm = motifs.matrix.PositionWeightMatrix.log_odds(m.counts.normalize(psw))

	return pssm, m.name

def scan_sequence(pssm, seq, scorethreshold):
	"""Return a liste of couple position/score of hits in the promotor sequences with the psw that are
	greater or equal to the given scorethreshold.

	Note: both=False refere to only to positions on the positive srand.
	"""
	return list(pssm.search(seq,scorethreshold,both=True)) #-> [(pos,score),..]

def scan_all_sequences(pssm, list_seq, scorethreshold):
	"""Create a dictionnary of couple position/score for multiple sequences.
	
	It uses scan_sequence(pssm, seq, scorethreshold) compute for each seq the position and score of the TBFS hits,
	store them into a dictionnay with the key being the sequence id and values the list of (position/score). 
	"""

	data_struct = dict()
	seq_obj_list = Seq_obj_from_files(list_seq) # Returns a list of SeqIO object from a list of.

	for seq_object,mrna in zip(seq_obj_list,mRNAs):
		data_struct[mrna] = scan_sequence(pssm, seq_object.seq, scorethreshold)
	return data_struct # -> {id1:[(pos,score,),...],id2:[(pos2,score2),..],....}

def score_window(sas,start,end):
	"""Scan a window between a start and end positions ( points ).

	This function uses the sliding window algorithm to analyze a certain window between a starting and ending 
	points. The main algorithm ( approach ) used to calculate the score, is to check if hits (occurences) of a set
	of genes are near to each other; For instance, we can say that to points are in a proximity if their distance is 
	close to zero. So we calculate the diffrence like this : diff = abs(pointA - pointB), then we campute a score on the difference:
			 A -> 1 ( 100 %)
			 B -> x ( ?% )
			 diff = abs ( A - B )
			 x = diff/B 
	Here x means the percentage of the distance between A and B, if x is less or equal then the third of the window size / 100 ( to get 
	the percentage ), then A is near to B, else A is not in proximty of B, those occurences are far from each other and the score of this window 
	is not the preferead score ( see best_window ).
	It returns start, end , sum(scores), and widow info. The window info contains the hits of the sequences. 
	"""


	samples = list() # A list to hold samples from each sequences.
	percentage = list()  # A list to hold scores of distances between elements of consecutive rows.
	occurences = list()# A list to hold all the scores that are less or equal than a percentage threshold.
	i = 0 # Index.
	r = 0 # Rows.
	window_info = dict()

	# Iterate througth our dict of seq object.
	for item in sas:
		# For each sequence, create a row.
		row = []
		# Iterate througth pos,score in each sequence.
		for position, score in sas[item]:
			# If the position is neg, it means we are in the - srand.
			window_info[item] = position,score
			if position < 0:
				# Calculate the correcte position.
				# position = prom_len + position + 1
				# If the algorithm finds a occurence in the window, then 
				# append it to the row.
				if (prom_len + position + 1) >= start and (prom_len + position + 1) < end:
					row.append(prom_len + position + 1)
			else:
				# Same as above but in the + srand.
				neg_srand = False
				if position >= start and position < end:
					row.append(position)
		# If we have non-empty row,
		# add it to the samples.
		if row != []:
			samples.append(row)
			r+=1 # Increment the rows.
		else:
			pass
	
	# If non-empty samples.
	if samples != []:
		# For each row, compare each element of it with
		# the next row, and so on.
		for i in range(r-1):
			curr_row = samples[i-1]
			next_row = samples[i]
			for e in curr_row:
				for n in next_row:
					# Calculate the distence between two positions.
					diff = abs(e-n)
					# Give a percentage like score for the difference,
					# Then add it to the percentage list.
					percentage.append(diff/n)
			# For each elements of percentage, take only those with
			# a value less or equal then threshold, i.e that all this values
			# are in proximity.
			occurences = [ x for x in percentage if x < seuil ]	
	# Sum up the score for this window.		
	
	

	return start, end, sum(occurences), window_info


def best_window(sas,window_size):
	"""Calculate the best window score by scanning a set of sequences.

	Start from the head of sequences, and then slide througth the set, between each starting point in a sequence 
	and the size of the window, calls scan_window to compute a score on this locations, then check if the score is 
	lower or equal to the window threshold given as argument.

	Returns windows infos, a dictionnary with the score, start and end position and the window info on this specific location.
	"""

	global seuil # Define a global threshlod for the score_window condition.
	seuil = (window_size/3)/100 # Take the percentage of the third of the window size.

	start_bw = 0 # Starting point for the best window.
	end_bw = 0 # End point for the best window.
	start = 0 # Start the scan at position 0 for all the sequences.
	bws= 0 # Define a bws for the scans.
	tmp = 0
	end = window_size # Set the end for the window scan.
	window_index = 1
	
	slide_step = 7 # Arbitrary value to indicate the slide step of the window

	windows_info = dict()

	# Scann all the sequences entil reaching the end.
	while start < (prom_len - window_size ):
		# Get the score of the current window between start and end positions.
		start_curr, end_curr, curr_bws, window_info = score_window(sas,start,end)
		# Check if the current score is valid.
		if curr_bws > 0 and w_threshold > curr_bws:
			# To get ride of same value.
			if tmp != curr_bws:
				# Set it start and end window as the best start window.
				start_bw = start_curr 
				end_bw = end_curr
				windows_info[window_index] = (start_bw,end_bw,curr_bws,window_info)
				window_index+=1
			tmp = curr_bws
		# Take a step further.
		start+= slide_step 
		end += slide_step 
	return windows_info

	

def search_luncher(list_mRNA,jaspar_matrix_name,psw,threshold,len_prom,window_size,window_threshlod):
	""" Act as a luncher for the putative_TFBS.py script.

	This function is like a main, takes all the arguments from the input from users, process 
	the search, and returns TF name and a dictionnary of the winodws informations.
	"""
	global prom_len
	prom_len = len_prom

	global w_threshold
	w_threshold = window_threshlod

	global mRNAs
	mRNAs = list_mRNA

	pssm, tf_name = pwm2pssm(jaspar_matrix_name,psw)	

	file_path_sequences="../data/sequences" # directory to store downlaoded promotors.
	
	# Check if theirs no .fa file, so the program start dowloading data, else just use the data on the dir.
	if os.path.exists(file_path_sequences) and os.path.isdir(file_path_sequences):
		if not any(fname.endswith('.fa') for fname in os.listdir(file_path_sequences)):
			print("Empty directory.")
			download_promotors(mRNAs,prom_len,file_path_sequences) # dowload promotors.
		else:
			print("Using files from directory.")
	else:
		print("Given directory doesn't exist")
		print("Creating the '"+file_path_sequences+"' directory")
		os.makedirs(file_path_sequences)
		download_promotors(mRNAs,prom_len,file_path_sequences) # dowload promotors.

	
	file_list = glob.glob(os.path.join(os.getcwd(), file_path_sequences, "*.fa")) # create a list with all the files in the dowloand directory.
	dict_seq = scan_all_sequences(pssm,file_list,threshold)


	return tf_name, best_window(dict_seq,window_size)


"""
if __name__ == "__main__":

	accession = sys.argv[1]
	path = sys.argv[2]
	psw = float(sys.argv[3])
	threshold = float(sys.argv[4])
	seq_length = int(sys.argv[5])
	
	pssm, motif_name = pwm2pssm(path,psw)

	global prom_len
	prom_len = seq_length

	global w_threshold
	w_threshold = 0.2



	file_path_sequences='../data/sequences' # directory to store downlaoded promotors.
	
	global list_mRNA
	list_mRNA = ["NM_007389", "NM_079420", "NM_001267550","NM_002470", "NM_003279", "NM_005159", "NM_003281", "NM_002469", "NM_004997", "NM_004320", "NM_001100", "NM_006757"]
	
	#list_mRNA = ["NM_007389", "NM_079420", "NM_001267550"]
	
	if os.path.exists('../data/sequences') and os.path.isdir('../data/sequences'):
		if not any(fname.endswith('.fa') for fname in os.listdir('../data/sequences')):
			download_promotors(list_mRNA,seq_length,file_path_sequences) # dowload promotors.
		else:
			print("Directory is not empty")
	else:
		print("Given directory doesn't exist")

	file_list = glob.glob(os.path.join(os.getcwd(), file_path_sequences, "*.fa")) # create a list with all the files in the dowloand directory.
	dict_seq = scan_all_sequences(pssm,file_list,threshold)

	
	for item in dict_seq:
		for pos, score in dict_seq[item]:
			print("%s %d %f" %(motif_name, pos, score))
	
	


	#score_window(dict_seq,0,40)


	wsi = best_window(dict_seq,40)
	
	for i in wsi :
		print("Window id : "+str(i)+" | TF : "+motif_name+" | Window pos : [%d,%d] | Window score : %f" %(wsi[i][0],wsi[i][1],wsi[i][2]))
		wi = wsi[i][3]
		for j in wi :
			print(""+str(i)+" "+motif_name+"  %s  %d  %f" %(j,wi[j][0],wi[j][1]))
		
"""



