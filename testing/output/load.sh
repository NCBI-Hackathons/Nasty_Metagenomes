cat ~/Nasty_Metagenomes/testing/output/ERR1600439.ga.fa_vs_Bacteria_RefSeq_nr_blastn_parsed.out | while read line; do echo -e "ERR1600439\t$line"; done > /tmp/blast.out
bq load --source_format=CSV --field_delimiter=tab nasty.contig_blast /tmp/blast.out run:STRING,contig:STRING,ref:STRING,species:STRING,origin:STRING,alignment_score:INTEGER,identity:FLOAT

cat contigs.bed | while read line; do echo -e "ERR1600439\t$line"; done > contigs.acc_bed
bq load --source_format=CSV --field_delimiter=tab nasty.beds contigs.acc_bed run:STRING,contig:STRING,start:INTEGER,stop:INTEGER,AMR:STRING,score:INTEGER,strand:STRING


cat calls_amr.vcf | while read line ; do echo -e "ERR1600439\t$line"; done | cut -f1-9 > calls_amr_acc.vcf 
zaluninvv@hackathon2:/data/ERR1600439/magicblast_output$ bq load --source_format=CSV --field_delimiter=tab nasty.snp calls_amr_acc.vcf run:STRING,contig:STRING,pos:INTEGER,id:STRING,ref:STRING,alt:STRING,qual:STRING,filter:STRING,info:STRING
