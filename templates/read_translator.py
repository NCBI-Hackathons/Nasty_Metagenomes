import Bio, sys
from Bio import SeqIO
from datetime import datetime

startTime = datetime.now()

if len(sys.argv) not in [2,3]:
    print('Incorrect # of arguments! Please supply path to FASTQ as argument to script.')
    exit()
#If in testmode, only look at the first 100,000 reads.
testmode=False

#This function writes a record to the output file.
def write_record(record_id, orfseq, output_file, translation_num, orientation):
    output_line = ''.join(['>', record_id, '_%s_'%orientation,
                        str(translation_num), '\n', orfseq, '\n', '\n'])
    output_file.write(output_line)

#Loop through the read file supplied as an argument.
#For each read, translate in all 6 possible reading frames and break up
#into ORFs. For all ORFs with lengths greater than 25 aa,
#write to the query file that we will search the HMMs against in the
#next step.
with open(sys.argv[1]) as in_file:
    with open('translated_reads.fa', 'w+') as out_file:
        for num_reads, record in enumerate(SeqIO.parse(in_file, 'fastq')):
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
                if testmode == True:
                    print('running in test mode; 25,000 reads processed.')
                    break

print('This script processed %s reads'%num_reads)
print('Time elapsed: %s'%(datetime.now() - startTime))
