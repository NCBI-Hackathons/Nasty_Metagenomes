#!/usr/bin/env python3


import Bio, sys, multiprocessing, subprocess
from multiprocessing import Process
from Bio import SeqIO
from datetime import datetime

startTime = datetime.now()

fastq = "$fastq"

#This function writes a record to the output file.
def write_record(record_id, orfseq, output_file, translation_num, orientation):
    output_line = ''.join(['>', record_id, '_%s_'%orientation, str(translation_num), '\n', orfseq, '\n', '\n'])
    output_file.write(output_line)

#Loop through the read file supplied as an argument.
#For each read, translate in all 6 possible reading frames and break up
#into ORFs. For all ORFs with lengths greater than 25 aa,
#write to the query file that we will search the HMMs against in the
#next step.
def process_subfile(subinput_name, suboutput_name):
    with open(subinput_name, "r") as in_file:
        with open(suboutput_name, 'w+') as out_file:
            for num_reads, record in enumerate(SeqIO.parse(in_file, 'fasta')):
                counter = 0
                rev_comp = record.seq.reverse_complement()
                for i in range(0,3):
                    translation_length = len(record) - (len(record[i:])%3)
                    forward_pass = record.seq[i:translation_length].translate()
                    reverse_pass = rev_comp[i:translation_length].translate()
                    if '*' in forward_pass:
                        for orf in str(forward_pass).split('*'):
                            if len(orf) > 30:
                                write_record(record.id, orf, out_file, counter, 'forward')
                                counter += 1
                    else:
                        write_record(record.id, str(forward_pass), out_file, counter, 'forward')
                        counter += 1
                    if '*' in reverse_pass:
                        for orf in str(reverse_pass).split('*'):
                            if len(orf) > 30:
                                write_record(record.id, orf, out_file, counter, 'reverse')
                                counter += 1
                    else:
                        write_record(record.id, str(reverse_pass), out_file, counter, 'reverse')
                        counter += 1
                if num_reads % 250000 == 0 and num_reads > 0:
                    print('%s reads complete'%num_reads)

if __name__ == "__main__":
    line_count = 0
    with open(fastq) as inpt:
        for line in inpt:
            line_count += 1
    
    print('counting records complete')

    twentypercent_split = (line_count / 4) // 5
    num_to_terminate = twentypercent_split
    subset_files = [open('temp_sub_file%s.fa'%i, 'w+') for i in range(0,5)]
    current_file = 0
    for num_reads, record in enumerate(SeqIO.parse(open(fastq), 'fastq')):
        subset_files[current_file].write(''.join(['>',record.id, '\n', str(record.seq), '\n']))
        if num_reads > num_to_terminate and current_file < 4:
            current_file += 1
            num_to_terminate += twentypercent_split
    print('%s reads processed'%num_reads)
    for subset_file in subset_files:
        subset_file.close()
    jobs = [Process(target=process_subfile, args=('temp_sub_file%s.fa'%i, 'temp_prot%s.fa'%i)) for i in range(0,5)]
    for job in jobs:
        job.start()
    for job in jobs:
        job.join()

    subprocess.call('cat temp_prot0.fa temp_prot1.fa temp_prot2.fa temp_prot3.fa temp_prot4.fa > translated_reads.fa', shell=True)
