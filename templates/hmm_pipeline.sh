#!/bin/bash
set -e
set -u
set -o pipefail

while getopts ':q::h::a:' OPTION; do
       case "$OPTION" in
	       q)    queryfile="$OPTARG";;
	       a)    threshold="$OPTARG";;
	       h)    hmmdir="$OPTARG";;
               ?)    echo "Script usage: [-h] [dir containing hmm profiles] [-q] [fastq file] [-a] [threshold bitscore for hmmer]"
       esac
done      

if (( $OPTIND == 1)); then
	echo "Incorrect arguments supplied."
	echo "To run hmm_pipeline.sh, supply the following arguments:"
	echo "[-h] [directory containing hmm profiles] [-q] [fastq file]"
	echo "[-a] [bitscore threshold for read to be considered a hit]"
	echo "in any order."
	exit 1
fi

#Start by translating each read in all 6 possible reading frames using
#the read_translator.py script.
count=0
python read_translator.py $queryfile -testmode
echo "reads translated; beginning HMMsearch; time for each profile will be printed"
#Use hmmsearch to bounce all 562 HMMs against the fasta file containing the translated reads
for hmmfile in $hmmdir/*.HMM;
do
	hmmsearch --noali -T $threshold --acc --incT $threshold -o "hmmresult${count}" $hmmfile translated_reads.fa
	count=`expr $count + 1`
done

echo "hmm searches complete"
#Move all the results files to a temporary directory we will delete later.
mkdir hmm_results_temp
mv hmmresult* hmm_results_temp/

#Use the hits_compiler script to retrieve a representative nucleotide sequence
#corresponding to each hit from the nf_to_seq fasta file and output to
#hit_nucleotide_seqs.fa, which will be the target for the guided assembler.
python hits_compiler.py $queryfile

#Remove all the temporary files
rm hit_reads.fa translated_reads.fa
rm -r hmm_results_temp
echo "HMM analysis complete; moving to guided assembly."
