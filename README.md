# Bioinformatique : cis-regulating TFBS searcher.

# Description

The project consists in solving the problem of finding cis-regulating TFBS of a set of genes.
Using the `Biopython` module and the interaction with the **NCBI** and **JASPAR** databases, we realized a program in `python` which allows to look for TFBS on a set of genes supposed co-regulated by a TF.
The method used is the sliding window, and our approach is the of calculation of the distance between the occurrences (hits) of TFBS on promoter sequences (see slide.).

# Project tree :

    ├── README.md          
    │
    ├── data
    │    ├── motifs         
    │    └── sequences     
    │
    ├── src
    │    ├── utiles.py
    │    ├── pwm.py
    │    ├── putative_TFBS.py
    │    └── scan_pwm.py
    │
    ├── TP2  
    |	├── TP2.md
    |	└── TP2_bin.py
    │       
    └── test.sh


# Getting started :

Before cloning the project into your local disk, please check the requierments.

## Requierments :

- `python 3.8`
- `valide SSL CERTIFICATE for python`
- `biopython`

## Usage :

    git clone https://github.com/Nissou31/bio-informatique.git
    cd bio-informatique
    ./test.sh    

The bash script will run automatically the program with differents settings.

To test only the **putative_TFBS.py**

	cd src
	python3 putative_TFBS.py --help

Results : 

	usage: putative_TFBS.py [-h] -t T [-l L] [-w W] -m M -s S [-p P] [mrna [mrna ...]]

	Find the index of the best window for a TFBS.

	positional arguments:
		mrna                  list of mRNA id.

	optional arguments:
		-h, --help            show this help message and exit
		-t T, --threshold T   score threshold for pssm
		-l L, --promotor-length L
                        	  promotor sequence length
        -w W, --window-size W
                        	  sliding window size
        -m M, --pfm M         jaspar matrix name
        -s S, --window-threshold S
                       		  window score threshold
        -p P, --pseudocount P
                         	  pseudo-count weight

# Autors :

- AMRI Anes Saad Eddine.

