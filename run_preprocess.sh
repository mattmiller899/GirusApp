#!/bin/bash
#PBS -l select=1:ncpus=7:mem=42gb
#PBS -l walltime=24:00:00

KMER=$PBS_ARRAY_INDEX
cd $PBS_O_WORKDIR
#module load singularity
source activate pyenv
INPUT_DIR="/rsgrps/bhurwitz/mattmiller899/girus_work/FASTA_files/${ORG}_unsplit/${CONTIG}_contigs"
OUTPUT_DIR="/rsgrps/bhurwitz/mattmiller899/girus_work/COMPILED_FEATURES_test/${ORG}/${CONTIG}_contigs/${KMER}"
init_dir "$OUTPUT_DIR"
#singularity run ./keras_gpu3.img pipeline.py -i $INPUT_DIR -o $OUTPUT_DIR -e $EXT -k $KMER
echo "python pipeline.py -i $INPUT_DIR -o $OUTPUT_DIR -e $EXT -k $KMER -m 1"
python pipeline.py -i $INPUT_DIR -o $OUTPUT_DIR -e $EXT -k $KMER -m 1
