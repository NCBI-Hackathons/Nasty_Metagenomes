# Nasty_Metagenomes
Antimicrobial Resistance Characterization in Metagenomes

### Previous Hackathon Work
* [NastyBugs](https://github.com/NCBI-Hackathons/MetagenomicAntibioticResistance)
* [Bugs_And_Drugs](https://github.com/NCBI-Hackathons/Bugs_And_Drugs) (Based on MagicBlast)

### Goals:
* Compare alignment of FASTQ files from SRA Id's to MagicBlast and HMM-er using reference AMR databased (CARD, Resfinder, and ARG-ANNOT)
* Generate a "hits" file with aligned reads to each of the publically available AMR databases
* Determine contextual information about AMR genes by 1) BLAST against RefSeq database to determine species level information 2) determine whether the AMR's are on a plasmid

### Tools
* [MagicBlast](https://ncbi.github.io/magicblast/)
* [HHM-er](https://github.com/EddyRivasLab/hmmer)
* [SamTools](https://github.com/samtools)
* [Skesa](https://github.com/ncbi/SKESA)
* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/)

### AMR Databases
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

1. Create Blast Databases
2. Use SamTools
3. Input SRA/FASTQ to MagicBlast or HHM-er
4. SKESA guided assembly
5. Species identification, plasmid identification

### Step 1.
# ------------------
    # Download Plasmid Database
    
    # Download AMR Finder Database
    
    # Download [bacterial chromosome Databases] (https://www.ncbi.nlm.nih.gov/assembly)
    ## Search assemblies all[sb]
    ## Download Assembly: Bacteria, Latest RefSeq, Assembly from Type
    ## Download Assembly: Bacteria, Latest RefSeq, Reference
    
    # Merge Plasmid, Assembly from Type, Reference Databses
    
    
    # Create Blast Databases
    ## AMR blast db
    ## Merged Plasmid, Assembly from Type, Reference
      
    #**command line**
    makeblastdb -in .fna -dbtype nucl -parse_seqids -out 
    makeblastdb -in all_plasmid_genomic.fna -dbtype nucl -parse_seqids -out plasmid_db

      
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

      
    # HHM-er
    # Use AMR finder database as a reference
      
    #**command line**

### Step 4.
# ------------------
    # SKESA Guided Assembly

    #**command line**

### Step 5.
# ------------------
    # Species Identification
    # Input file to blast against RefSeq database to create species ID list
    # Cross reference list to AMR genes
    
    #**command line**
    
    # Plasmid Identification
    # input file to determine which AMR's are on plasmid's
    
    #**command line**
