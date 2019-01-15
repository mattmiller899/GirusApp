#!/bin/bash
#PBS -l select=1:ncpus=28:mem=168gb:ngpus=1
#PBS -l walltime=12:00:00

cd $PBS_O_WORKDIR
source activate pyenv
module load singularity
singularity run ./keras_gpu3.img ./save_models.py -k $KMER -i $IN -o $OUT -w $WORKER 
