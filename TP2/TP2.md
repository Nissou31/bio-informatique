## TP2 

# Binôme :
- AMRI Anes SaadEddine 
- NASRI Ahmed-Faez

## Réponses :

# Lécture de fichiers :

1. Récuperation de séquences des entrées :
		
		from Bio import SeqIO

		recordFASTA = SeqIO.read("data/sequence.fasta", "fasta")
		recordGB = SeqIO.read("data/sequence.gb", "genbank")
		
		print("Sequence alphabet %s" % recordFASTA.seq)
		print("Sequence alphabet %s" % recordGB.seq)

2. Récuperation de séquences sous forme str :
		
		string_seq_gb = str(recordGB.seq)
		string_seq_fasta = str(recordFASTA.seq)
		
		# Check if both string are equals.
		if string_seq_gb == string_seq_fasta :
			print("The sequences in both format are identical.")

3. Pour obtenir la séquence complémentaire inverse , on utilise la fonction `reverse_complement`:
		
		recordFASTA_rc = recordFASTA.reverse_complement(id=recordFASTA.id + "_rc")
		print("Reverse sequence alphabet %s" % str(recordFASTA_rc.seq))
	
	Mais dans le format GenBank, on obtient une `ValueError: DNA found, RNA expected`, on peut contourner cette erreur on passant a la fonction
	comme parametre la sequence :
		
		recordGB_rc = recordGB.seq.reverse_complement()
		print("Reverse sequence alphabet %s" % str(recordGB_rc))

4. Les attributs pour les entrées sont :
		
		print(recordFASTA.features)

		for f in recordGB.features :
			print(f)
		
	- type
	- location
	- qualifiers{ chromosome ,db_xref ,ap, mol_type, organism, strain }
    		
    Ces attributs on les retrouvent seulement dans le format GenBank, car dans le format FASTA on obtient un vecteur vide.

5. On peut accéder a ses positions de début et de fin :
	
		first_feature = recordGB.features[0]

		print(first_feature.location)

# Fonction :

* Fonction `find_cds`:

		def find_cds(seqrecord):
	
			result = []
			couplePosition = tuple()

			for f in seqrecord.features :
				if f.type == "CDS":
					couplePosition = (list(f.location)[0],list(f.location)[-1])
					result.append(couplePosition)

			return result

# API du NCBI avec Biopython

### efetch :

1. Utilisez la méthode efetch pour récupérer l'entrée NM_007389 au format Genbank puis au format FASTA :
		
		Entrez.email = "email@example.org"
		
		handleGB = Entrez.efetch(db="nucleotide", id="NM_007389", rettype="gb", retmode="text")
		handleFASTA = Entrez.efetch(db="nucleotide", id="NM_007389", rettype="fasta", retmode="text")

		recordFetchGB = SeqIO.read(handleGB, "genbank")
		recordFetchFASTA = SeqIO.read(handleFASTA, "fasta")

2. Vérifiez que les séquences des entrées obtenues sont bien identiques
	
		string_seq_gb = str(recordGB.seq)
		string_seq_fetch = str(recordFetch.seq)
		if string_seq_gb == string_seq_fetch :
			print("The sequences in both records are identical.")

3. En utilisant votre méthode find_cds vérifiez que les CDS sont bien identiques.

		localres = find_cds(recordGB)
		fetchres = find_cds(recordFetch)
		if localres == fetchres :
			print("The CDS location in both records are identical.")

### elink :


1. Utilisation de la méthode elink pour récupérer le gène correspondant à l'entrée NM_007389 :
	
		handle = Entrez.elink(dbfrom="nuccore",db="gene", id="NM_007389")
		record = Entrez.read(handle)

2. Identifiant du gène :
	
		linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]

		for obj in linked :
			print(obj)

3. Méthode mrna_to_gene qui prenne un numéro d'accession d'un ARNm et qui renvoie l'identifiant du gène correspondant (ou qui lève une exception ValueError en cas de problème) :
		
		def mrna_to_gene(num_accession):
	
			handle = Entrez.elink(dbfrom=dbfrom,db=db, id=num_accession)
			record = Entrez.read(handle)
			handle.close()
	
			try:
				linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
				return linked[0]
			except ValueError:
				print("No id found...")



### Récupération de la portion amont d’un gène :


1. Méthode esummary :

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

		print(info)

2. À partir de l’identifiant de ce chromosome,on obtient la séquence en amont du gène comme ce ci :
	
		handle = Entrez.efetch(db='nucleotide', id='NC_000068.8', rettype='fasta', strand=-1,retmode='text', seq_start=73393625, seq_stop=73410682)
		record = handle.readline().strip()
		genomic_record = SeqIO.read(handle, "fasta")

		handle.close()

		print(genomic_record.seq)


3. Fonction upstream_gene_seq :

		def upstream_gene_seq(id_gene, seq_length, seq_start, seq_stop):

			handle = Entrez.efetch(db=dbfrom, id=id_gene, rettype='fasta', strand=-1,retmode='text', seq_start=seq_start, seq_stop=seq_stop)
			genomic_record = SeqIO.read(handle, "fasta")

			handle.close()


			return genomic_record.seq[:seq_length]


### Calcul de score à partir de matrices de fréquences :

1. Nombre de matrices lues : 1

2. Pour obtenir la matrice *pwm* :

		from Bio import motifs

		with open("data/motifs/MA0114.4.jaspar") as handle:
			m = motifs.read(handle, "jaspar")

		pmw = m.counts.normalize(0.1)
		
		print(pmw)

	On met la valeur pour le pseudo-count a **O.1**.

3. On obtient une *pssm* a partir d'une *pwm* de la façon suivante :

		pssm = motifs.matrix.PositionWeightMatrix.log_odds(pmw)

4. Avec la fonction `search()`.























