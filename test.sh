#!/bin/sh

cd src

echo "Run search for MA0006" 
python3 putative_TFBS.py NM_007389 NM_079420 NM_001267550 NM_002470 NM_003279 NM_005159 NM_003281 NM_002469 NM_004997 NM_004320 NM_001100 NM_006757 -m MA0006.1 -p 0.1 -t 3.0 -l 1000 -w 40 -s 0.3

echo "Run search for MA0114"
python3 putative_TFBS.py NM_007389 NM_079420 NM_001267550 NM_002470 NM_003279 NM_005159 NM_003281 NM_002469 NM_004997 NM_004320 NM_001100 NM_006757 -m MA0114.4 -p 0.1 -t '-2.0' -l 1000 -w 40 -s 0.3

echo "Run search for MA0083"
python3 putative_TFBS.py NM_007389 NM_079420 NM_001267550 NM_002470 NM_003279 NM_005159 NM_003281 NM_002469 NM_004997 NM_004320 NM_001100 NM_006757 -m MA0083.3 -p 0.1 -t 1.0 -l 1000 -w 40 -s 0.3

echo "Run search for MA0057" 
python3 putative_TFBS.py NM_007389 NM_079420 NM_001267550 NM_002470 NM_003279 NM_005159 NM_003281 NM_002469 NM_004997 NM_004320 NM_001100 NM_006757 -m MA0057.1 -p 0.1 -t '-1.0' -l 1000 -w 40 -s 0.3

echo "Run search for MA0056"
python3 putative_TFBS.py NM_007389 NM_079420 NM_001267550 NM_002470 NM_003279 NM_005159 NM_003281 NM_002469 NM_004997 NM_004320 NM_001100 NM_006757 -m MA0056.2 -p 0.1 -t 0.8 -l 1000 -w 40 -s 0.3