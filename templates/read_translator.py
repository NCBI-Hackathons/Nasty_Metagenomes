import gzip, Bio, sys
from Bio import SeqIO
from datetime import datetime

startTime = datetime.now()

if len(sys.argv) is not in [2,3]:
    print('Incorrect # of arguments! Please supply path to gzipped FASTQ as argument to script.')
    exit()

testmode=False
if sys.argv[2] == '-testmode':
    testmode = True

def write_record(record_id, orfseq, output_file, translation_num, orientation):
    output_line = ''.join(['>', record_id, '_%s_'%orientation,
                        str(translation_num), '\n', orfseq, '\n', '\n'])
    output_file.write(output_line)


with gzip.open(sys.argv[1], 'rt') as in_file:
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
        if num_reads % 50000 == 0:
            print('%s reads processed'%num_reads)
            if testmode == True:
                print('Running in test mode; first 50,000 reads processed.')
                exit()

print('This script processed %s reads'%num_reads)
print('Time elapsed: %s'%(datetime.now() - startTime))
