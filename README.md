# Nasty_Metagenomes
Antimicrobial Resistance Characterization in Metagenomes

### Previous Pipelines
* [NastyBugs](https://github.com/NCBI-Hackathons/MetagenomicAntibioticResistance)
* [Bugs_And_Drugs](https://github.com/NCBI-Hackathons/Bugs_And_Drugs) (Based on MagicBlast)

### Goals:
* Compare alignment of FASTQ files to: MagicBlast and HMM-er
* Generate a "hits" file with aligned reads to each of the publically available AMR databases: [CARD](https://card.mcmaster.ca/), [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/), [ARG-ANNOT](https://omictools.com/arg-annot-tool)
* Determine contextual information about AMR genes by 1) BLAST against RefSeq database to determine species level information 2) determine whether the AMR's are on a plasmid

### Tools
* Nextflow
* Docker
* 

### AMR Databases
* [CARD](https://card.mcmaster.ca/)
* [Resfinder](https://cge.cbs.dtu.dk/services/ResFinder/)
* [ARG-ANNOT](https://omictools.com/arg-annot-tool)

### Alignment Tools
* [MegaBlast](https://ncbi.github.io/magicblast/)
* [HHM-er](https://github.com/EddyRivasLab/hmmer)

## Workflow Diagram
![workflow](https://github.com/NCBI-Hackathons/Nasty_Metagenomes/blob/master/images/Workflow.JPG)

## Workflow Steps
* Input: SRA ID
* Output files: AMR hits file, AMR by species, AMR's on plasmids

1) Input SRA ID to MagicBlast OR Input SRA ID to HHM-er
2) Take output "hits" file, and 
--input file to blast against RefSeq database to determine which ARM's are assigned to which species
--input file to determine which AMR's are on plasmid's
