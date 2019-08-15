# Nasty_Metagenomes
Antimicrobial Resistance Characterization in Metagenomes

### Previous Hackathon Work
* [NastyBugs](https://github.com/NCBI-Hackathons/MetagenomicAntibioticResistance)
* [Bugs_And_Drugs](https://github.com/NCBI-Hackathons/Bugs_And_Drugs) (Based on MagicBlast)

### Goals:
* Compare alignment of FASTQ files 1) MagicBlast 2) HMM-er and 3) MASH using reference AMR Finder database
* Generate a "hits" file with aligned reads to each of the publically available AMR databases
* Determine species and plasmid contextual information about AMR genes by creating a merged chromosome and plasmid database and BLAST-ting "hits" 

### Dependencies & Tools
* [MagicBlast](https://ncbi.github.io/magicblast/)
* [HHM-er](https://github.com/EddyRivasLab/hmmer)
* [MASH](https://github.com/marbl/Mash)
* [SamTools](https://github.com/samtools)
* [Skesa](https://github.com/ncbi/SKESA)
* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/)
* [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
* [Biopython](https://anaconda.org/anaconda/biopython)

### AMR Database
* [AMR Finder](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047), which includes:
  * [CARD](https://card.mcmaster.ca/)
  * [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/)
  * [Lahey](https://externalwebapps.lahey.org/studies/)
  * [Pasteur Institute Beta Lactamases](https://bigsdb.pasteur.fr/klebsiella/klebsiella.html)

## Workflow Diagram
![workflow](https://github.com/NCBI-Hackathons/Nasty_Metagenomes/blob/master/images/Workflow.JPG)

## Workflow Steps
* Input: SRA ID
* Output files: AMR hits file, AMR by species, AMR's on plasmids

1. [Create Blast Databases](#Step-1)
2. [Use SamTools](#Step-2)
3. [Input SRA/FASTQ to MagicBlast or HHM-er](#Step-3)
4. [SKESA guided assembly](#Step-4)
5. [Species identification, plasmid identification](#Step-5)

### Step 1. 
# ------------------ 
    # Download [RefSeq Plasmid Database](https://www.ncbi.nlm.nih.gov/refseq/)
    ## Use [FTP](ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/) to download plasmid database, and concatinate into one file
    
    # Download [AMR Finder Database](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047)
    ## Use webserver to download database [AMR_CDS](ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinder/data/2019-04-29.1/)
    
    # Download [bacterial chromosome Databases](https://www.ncbi.nlm.nih.gov/assembly)
    ## Search assemblies all[sb]
    ## Download Assembly: Bacteria, Latest RefSeq, Assembly from Type
    ## Download Assembly: Bacteria, Latest RefSeq, Reference
    
    # Merge Plasmid, Assembly from Type, Reference Databses
    ## Create non-redudant list of accession numbers
    
    # Create Blast Databases
    ## AMR blast db
    ## Merged Plasmid, Assembly from Type, Reference
      
    #**command line**
    /opt/ncbi-blast-2.9.0+/bin/makeblastdb -in /data/DBs/Bacteria_type_rep_plasmid_cat_nr.fa -parse_seqids -input_type fasta -dbtype nucl -out /data/DBs/Bacteria_type_rep_plasmid_refseq_nr.blastdb -max_file_sz 4GB 
    # Create Mash sketches:
    ## mash sketch -i AMR.fa
    ## mash sketch -i -p 12 Bacteria_rep.fna
    ## mash sketch -i -p 12 Bacteria_typ.fna 
#### (Back to [Workflow](#Workflow-Steps))
    
    
### Step 2.
# ------------------
    # Use Sam tools

    #**command line**
    
#### (Back to [Workflow](#Workflow-Steps))


### Step 3.
# ------------------
    # MagicBlast
    # Use AMR finder database as reference
    # Use SamTools to sort and create depth summary file
    # Run Cov_dep_cal.pl for coverage depth and average sequence coverage
      
    #**command line**
    /opt/magicblast/ncbi-magicblast-1.4.0/bin/magicblast -sra ERR1600439 -db /data/AMR_CDS.blastdb -outfmt sam -out ERR1600439_v_AMR_CDS_magicblast_sam.out -num_threads 8 -paired -no_unaligned
      
    # HHM-er
    # Use AMR finder database as a reference; use a user-specified bitscore as a threshold
    # to filter hmm hits.
    # Translate each read into protein in all six possible reading frames
    # Break translated reads into ORFs; discard any length 25 aa or less
    # Use hmmsearch to run 562 HMM profiles against the translated reads
    # For each hmm profile that scored a hit, extract a representative nucleotide sequence
    # Output the representative nucleotide sequences for guided assembly
      
    #**command line**
    ./hmm_pipeline.sh -q [fastq file path] -h [path to hmm_databases directory] -a 
    [threshold bitscore]
    
    # MASH
    # Given a fasta file with AMR genes we build a MASH sketch and screen the reads against the sketch. 
    # This produces kmer distances between the read set and each AMR, which then is used to extract only 
    # the AMRs that are close the readset.
    
    #**command line**
    # Screen the reads against the AMR sequences with the minimum score 0.85:
    ## mash screen -p 12 -w AMR.fa.msh ERR1600439*.fastq | awk '$1>0.85' > ERR1600439.amr.screen
    # produce a list of candidate AMRs:
    ## cut -f 5 ERR1600439.amr.screen
#### (Back to [Workflow](#Workflow-Steps))


### Step 4.
# ------------------
    # SKESA Guided Assembly
    # Guided assembly allows to assemble contigs based on some known sequences used as baits. The assembler stacks kmers and extends the ends of each guide sequence optionally output a list of variants assembled. 
    # In our analysis we use AMR fasta sequences selected in previous steps as guides for assembly. The resutls are presented as contigs fasta file. 

    #**command line**
    # to assemble contigs with AMR_CDS_by_ERR1600439_ref.fasta as guides:
    # guidedassembler --cores 8 --sra_run ERR1600439 --targets /data/ERR1600439/magicblast_output/AMR_CDS_by_ERR1600439_ref.fasta --contigs_out ERR1600439.ga.fa --fraction 0.1
    # to assemble contigs and print out all variants of contigs:
    # guidedassembler_graph --targets ../../AMR_CDS_norm.fasta --consensus ERR1600439.amr.contigs.fa --all_variants ERR1600439.amr.all-contigs.fa --gfa /dev/null --sra_run ERR1600439
#### (Back to [Workflow](#Workflow-Steps))


### Step 5.
# ------------------
    # Species and Plasmid Identification
    # Blast AMR hits lists against combined database
    # Parse for Species level and Plasmid identification
    
    #**command line**
    sudo /opt/ncbi-blast-2.9.0+/bin/blastn -query /data/ERR1600439/magicblast_output/ERR1600439.ga.fa -task blastn -db /data/DBs/Bacteria_type_rep_plasmid_refseq_nr.blastdb -outfmt 6 -evalue 1e-6 -out /data/ERR1600439/magicblast_output/ERR1600439.ga.fa_vs_Bacteria_RefSeq_nr_blastn.out -max_target_seqs 100
#### (Back to [Workflow](#Workflow-Steps))


### Step 6.
# ------------------
    # Alignment and stats
    
    # build blast db for contigs
    # makeblastdb -parse_seqids -in ERR1600439.ga.fa -input_type fasta -dbtype nucl -out ERR1600439.ga.blastdb
    # align reads onto the contigs
    # magicblast -db ERR1600439.ga.blastdb -query ERR1600439*.fastq | samtools view -Sb - | samtools sort - > ERR1600439_amr_contigs.bam

## Authors
* Xin Huang
  * National Human Genome Research Institute, National Institutes of Health, Bethesda, MD 20851
* Inês Mendes
  * Instituto de Microbiologia, Instituto de Medicina Molecular, Faculdade de Medicina, Universidade de Lisboa, Lisboa, Portugal; University of Groningen, University Medical Center Groningen, Department of Medical Microbiology and Infection Prevention, Groningen, The Netherlands
* Jonathan Parkinson
  * Qpex Biopharma, Inc., San Diego, CA 92121
* Samantha Sevilla
  * Cancer Genomics Research Laboratory, Division of Cancer Epidemiology and Genetics, National Institutes of Health, Leidos Biomedical, Inc., Gaithersburg, MD 20877
* Vadim Zalunin
  * National Center for Biotechnology Information, National Library of Medicine, National Institutes of Health, Bethesda, MD 20894, USA
