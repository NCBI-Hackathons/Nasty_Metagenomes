#!/usr/bin/env python3

import Bio, os
from Bio import SeqIO


current_dir = os.getcwd()
file_list = "${hmmresult}".split()
hit_list = []

#Loop through all files in the temporary hmm_results_temp directory.
#Parse each file, looking to see whether it contains any results listed
#under the '------' line. If so, add the accession number stored early
#in the file to a list of hits.
for inp_file in file_list:
    file_handle = open(inp_file)
    collect_outputs = False
    for line in file_handle:
        if 'Accession:' in line:
            accession_id = line.strip().split()[1]
        if collect_outputs:
            if len(line) > 1:
                if '[No hits detected that satisfy' in line or 'Domain annotation' in line:
                    break
                else:
                    if '_forward_' in line or '_reverse_' in line:
                        hit_list.append(accession_id)
                        break
        if '------' in line:
            collect_outputs = True
    file_handle.close()

#Now go to the nf_to_seq.fa file, which contains a reprentative nucleotide
#sequence for each NF accession number. For each accession number that had a
#hit, retrieve the appropriate nucleotide sequence and write it to
#hit_nucleotide_seqs.fa, which will be our output to the next stage of
#the pipeline.
hit_set = set(hit_list)
os.chdir(current_dir)
hit_nucleotide_seqs = []
with open('hit_nucleotide_seqs.fa', 'w+') as hit_file:
    with open("${nf_to_seq}") as nf_map_file:
        for record in SeqIO.parse(nf_map_file, 'fasta'):
            if record.id.split('|')[0] in hit_set:
                hit_nucleotide_seqs.append(record)
        SeqIO.write(hit_nucleotide_seqs, hit_file, 'fasta')
