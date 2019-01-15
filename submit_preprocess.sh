#!/bin/bash

set -u
ORGS=(girus virus)
CONTIGS=(1000)
#export STDERR_DIR="./err"
#export STDOUT_DIR="./out"
#init_dir "$STDERR_DIR" "$STDOUT_DIR"
export EXT="fasta"
for x in ${ORGS[@]}; do
    for y in ${CONTIGS[@]}; do
        export ORG=$x
        export CONTIG=$y
        ARGS="-q standard -W group_list=bhurwitz -M mattmiller899@email.arizona.edu -m a"
        echo "$ORG $CONTIG"
        #JOB_ID=`qsub $ARGS -v ORG -N gc_content_${x} -e "$STDERR_DIR" -o "$STDOUT_DIR" -J 1-5 $GIRUS_DIR/gc_content.sh`
        #JOB_ID=`qsub $ARGS -v ORG,CONTIG,EXT -N test_preprocess -e "$STDERR_DIR" -o "$STDOUT_DIR" -J 2-6 ./run_preprocess.sh`
        JOB_ID=`qsub $ARGS -v ORG,CONTIG,EXT -N test_preprocess -J 2-6 ./run_preprocess.sh`

        if [ "${JOB_ID}x" != "x" ]; then
            echo Job: \"$JOB_ID\"
        else
            echo Problem submitting job. Job terminated.
            exit 1
        fi
        echo "job successfully submitted"
    done
done
