#!/bin/bash

set -u

source /rsgrps/bhurwitz/mattmiller899/girus_work/config.sh
CONTIGS=(3000)
KMERS=(2)
export STDERR_DIR="/rsgrps/bhurwitz/mattmiller899/GirusApp/err"
export STDOUT_DIR="/rsgrps/bhurwitz/mattmiller899/GirusApp/out"
init_dir "$STDERR_DIR" "$STDOUT_DIR"
#WORKERS=(2 4 6 8 10 12 14 16 18 20 22 24 26 28)
for x in ${CONTIGS[@]}; do
    for y in ${KMERS[@]}; do
            export CONTIG=$x
            export KMER=$y
            export WORKER=22
            export IN="/rsgrps/bhurwitz/mattmiller899/GirusApp/hdf_test/hdf_${CONTIG}_${KMER}_0_train.h5"
            export OUT="/rsgrps/bhurwitz/mattmiller899/GirusApp/models/model_${CONTIG}_${KMER}.h5"
            #export HDF="/rsgrps/bhurwitz/mattmiller899/girus_work/deep_learning/HDF_files/unsplit/hdf_mpi_${ORGS}_PBSIM_${KMER}.h5"
            #export OUT="/rsgrps/bhurwitz/mattmiller899/girus_work/deep_learning/${ORGS}_PBSIM_${KMER}_hdf5_output_mpi${WORKER}.txt"
            ARGS="-q standard -W group_list=bhurwitz -M mattmiller899@email.arizona.edu -m a"
            export NAME="model_${CONTIG}_${KMER}"
            JOB_ID=`qsub $ARGS -v CONTIG,KMER,WORKER,IN,OUT -N $NAME -e $STDERR_DIR -o $STDOUT_DIR ./run_save_models.sh`
            #JOB_ID=`qsub $ARGS -v ORGS,KMER,WORKER,HDF,OUT -N test_mpi ./run_create_hdf5_files_mpi.sh`
            if [ "${JOB_ID}x" != "x" ]; then
                echo Job: \"$JOB_ID\"
            else
                echo Problem submitting job. Job terminated.
                exit 1
            fi
            echo "job successfully submitted"
    done
done
