

# OM-HiC-scaffolding 

Here the scripts are used in the paper "Improving and correcting the contiguity of long-read genome assemblies of three plant species using optical mapping and chromatin conformation capture data". http://genome.cshlp.org/content/27/5/778


## Description
These scripts are used to do genome assembly scaffolding in three ways:

1)  scaffolding using optical consensus maps and PacBio read assembly contigs from Falcon and PBcR

2)  scaffolding using Dovetail Hi-C reads and PacBio read assembly contigs from Falcon and PBcR

3)  scaffolding using optical consensus maps, Dovetail Hi-C reads and PacBio read assembly contigs from Falcon and PBcR


## System requirements

Perl(version >5.10); Bioperl 

IrysSolve scripts http://www.bnxinstall.com/scripts/v3692/

RefAligner http://www.bnxinstall.com/RefalignerAssembler/Linux/SSE/

BLAST (makeblastdb, blastn)

HiRise (https://github.com/DovetailGenomics/HiRise_July2015_GR)

### requirements for HiRise        
    Python 3 
  
    BWA (version: 0.7.15-r1140)
  
    samtools (version:1.2)
  
    samblaster (https://github.com/GregoryFaust/samblaster)
  
    Meraculous (http://downloads.sourceforge.net/project/meraculous20/release-2.0.5.tgz)
  

## Usage
Please read the "README.md" and "run.sh" in each scaffolding worflow
