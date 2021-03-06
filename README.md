[![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com) &nbsp;
<img src="./data/img/using-biopython.svg"> &nbsp;
[![forthebadge](https://forthebadge.com/images/badges/open-source.svg)](https://forthebadge.com) &nbsp;


# Bioinformatique : cis-regulating TFBS searcher ๐งฌ๐ฌ

# Description ๐ :

The project consists in solving the problem of finding cis-regulating TFBS of a set of genes.
Using the `Biopython` module and the interaction with the **NCBI** and **JASPAR** databases, we realized a program in `python` which allows to look for TFBS on a set of genes supposed co-regulated by a TF.
The method used is the sliding window, and our approach is the calculation of the distance between the occurrences (hits) of TFBS on promoter sequences (see slide.).

# Project tree ๐ฒ :

    โโโ README.md          
    โ
    โโโ data
    โ    โโโ motifs         
    โ    โโโ sequences     
    โ
    โโโ src
    โ    โโโ utiles.py
    โ    โโโ pwm.py
    โ    โโโ putative_TFBS.py
    โ    โโโ scan_pwm.py
    โ
    โ       
    โโโ test.sh


# Getting started ๐ :

Before cloning the project into your local disk, please check the requierments.

## Requierments ๐งพ :

- `python 3.8`
- `valide SSL CERTIFICATE for python`
- `biopython`

## Usage ๐ป :

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

# Autors โ๏ธ :

- AMRI Anes Saad Eddine.

