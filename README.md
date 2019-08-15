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
    ## 
    ## Create non-redudant list of accession numbers
    
    
    # Create Blast Databases
    ## AMR blast db
    ## Merged Plasmid, Assembly from Type, Reference
      
    #**command line**
    /opt/ncbi-blast-2.9.0+/bin/makeblastdb -in /data/DBs/Bacteria_type_rep_plasmid_cat.fa -parse_seqids -input_type fasta -dbtype nucl -out Bacteria_type_rep_plasmid_refseq.blastdb 

      
### Step 2.
# ------------------  
    # Use Sam tools

    #**command line**

### Step 3.
# ------------------
    # MagicBlast
    # Use AMR finder database as reference
    # Use SamTools to sort and create depth summary file
    # Run Cov_dep_cal.pl for coverage depth and average sequence coverage
      
    #**command line**
    /opt/magicblast/ncbi-magicblast-1.4.0/bin/magicblast -sra ERR1600439 -db /data/AMR_CDS.blastdb -outfmt sam -out ERR1600439_v_AMR_CDS_magicblast_sam.out -num_threads 8 -paired -no_unaligned
      
    # HHM-er
    # Use AMR finder database as a reference
    # Translate all 6 reading frames
    # Create hit lits of representative AMR genes
      
    #**command line**
    
    # MASH
    # 
    
    #**command line**


### Step 4.
# ------------------
    # SKESA Guided Assembly

    #**command line**

### Step 5.
# ------------------
    # Species and Plasmid Identification
    # Blast AMR hits lists against combined database
    # Parse for Species level and Plasmid identification
    
    #**command line**
    /opt/ncbi-blast-2.9.0+/bin/blastn -query /data/ERR1600439/magicblast_output/ERR1600439.ga.fa -task blastn -db /data/DBs/Bacteria_type_rep_plasmid_refseq.blastdb -outfmt 6 -evalue 1e-6 -out /data/ERR1600439/magicblast_output/ERR1600439.ga.fa_vs_Bacteria_RefSeq_blastn.out

## Authors
Xin Huang
National Human Genome Research Institute, National Institutes of Health, Bethesda, MD 20851

Inês Mendes
Instituto de Microbiologia, Instituto de Medicina Molecular, Faculdade de Medicina, Universidade de Lisboa, Lisboa, Portugal; University of Groningen, University Medical Center Groningen, Department of Medical Microbiology and Infection Prevention, Groningen, The Netherlands

Jonathan Parkinson
Qpex Biopharma, Inc., San Diego, CA 92121

Samantha Sevilla
Cancer Genomics Research Laboratory, Division of Cancer Epidemiology and Genetics, National Institutes of Health, Leidos Biomedical, Inc., Gaithersburg, MD 20877

Vadim Zalunin
National Center for Biotechnology Information, National Library of Medicine, National Institutes of Health, Bethesda, MD 20894, USA
