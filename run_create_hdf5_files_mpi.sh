#!/bin/bash
#PBS -l select=1:ncpus=28:mem=168gb
#PBS -l walltime=12:00:00

cd $PBS_O_WORKDIR
source activate pyenv
module load singularity/2/2.6.1
cat $PBS_NODEFILE > $HOME/hosts/hosts_${CONTIG}_${KMER}
export PBS_NODEFILE=$HOME/hosts/hosts_${CONTIG}_${KMER}
#echo "PBS_NODEFILE = ${PBS_NODEFILE}"
#if [[ "$SPLIT" == "True" ]]; then
#    echo "hello"
#    singularity run ./h5py_mpi.img $WORKER create_hdf5_files_mpi.py -n $ORGS -k $KMER -c $CONTIG -d $HDF -o $OUT -f $SPLIT_FILE -z
#else
#    echo "howdy"
#    singularity run ./h5py_mpi.img $WORKER create_hdf5_files_mpi.py -n $ORGS -k $KMER -c $CONTIG -d $HDF -o $OUT -f $SPLIT_FILE
#fi
singularity run ./h5py_mpi.img $WORKER create_hdf5_files_mpi.py -k $KMER -c $CONTIG -d $HDF -p $POS_DIR -n $NEG_DIR -f $FOLDS


#singularity run ./h5py_mpi.img $WORKER create_hdf5_files_pbsim_mpi.py -n $ORGS -k $KMER -d $HDF -o $OUT
