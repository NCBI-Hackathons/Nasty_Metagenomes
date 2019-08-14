#!/bin/bash
set -e
set -u
set -o pipefail

while getopts ':q::h:' OPTION; do
       case "$OPTION" in
	       q)    queryfile="$OPTARG";;
	       h)    hmmdir="$OPTARG";;
               ?)    echo "Script usage: [-h] [dir containing hmm profiles] [-q] [gzipped fastq file]"
       esac
done       

count=0
python read_translator.py $queryfile -testmode
for hmmfile in $hmmdir/*.HMM;
do
	time hmmsearch --noali -T 50 --acc --incT 50 -o "hmmresult${count}" $hmmfile translated_reads.fa
	count=`expr $count + 1`
done

mkdir hmm_results_temp
mv hmmresult* hmm_results_temp/

python do_things_to_hmm_output #add this next
#after that add local assembly.... then delete all the temporary hmmresult files.
