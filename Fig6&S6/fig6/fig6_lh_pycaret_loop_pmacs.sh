#!/bin/bash
#BSUB -J R 
#BSUB -o R.%J.out
#BSUB -e R.%J.error
#BSUB -n 15
#BSUB -M 131072
#BSUB -R "span[hosts=1] rusage [mem=131072]" 

jupyter nbconvert --to notebook --execute pycaret-BLA-loop-sampling.ipynb