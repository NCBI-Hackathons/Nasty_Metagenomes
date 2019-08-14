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
* [CARD](https://card.mcmaster.ca/)
* [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/)
* [ARG-ANNOT](https://omictools.com/arg-annot-tool)

## Workflow Diagram
![workflow](https://github.com/NCBI-Hackathons/Nasty_Metagenomes/blob/master/images/Workflow.JPG)

## Workflow Steps
* Input: SRA ID
* Output files: AMR hits file, AMR by species, AMR's on plasmids

* STEP 1: Input SRA ID to:
  * MagicBlast
   * use AMR database as reference
   * Sort .sam file
   * Create depth of coverage file
   * Run Cov_dep_cal.pl for coverage depth and average sequence coverage
  * HHM-er
* STEP 2: Take output "hits" file, and 
  * Species Identification
   * Input file to blast against RefSeq database to create species ID list
   * Cross reference list to AMR genes
  * Plasmid Identification
   * input file to determine which AMR's are on plasmid's
