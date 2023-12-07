#!/bin/bash
#BSUB -J R 
#BSUB -o R.%J.out
#BSUB -e R.%J.error
#BSUB -n 15
#BSUB -M 10000
#BSUB -R "span[hosts=1] rusage [mem=10000]" 

Python fig6_pycaret_var100_AI_swapping.py