# OM-HiC-scaffolding
A workflow for correcting and scaffolding long-read(such as PacBio, nanopore) assemblies using optical mapping, Dovetail Hi-C and illumina mate-pair data. Scripts were used in our manuscript.


## filter raw PacBio reads using SMRTAnalysis


## run PacBio assembly using Falcon and PBcR, respectively
Falcon:fc_run.py fc_run.cfg



## do assembly polishing using filtered PacBio reads by QUIVER

## map illumina whole genome shotgun sequencing reads to the assemblies, call SNPs and InDels, and further correct the assemblies


## intergrated Optical map scaffolding
step 1:
step 2:
step 3:

## intergrated Dovetail Hi-C scaffolding


## intergrated OM and Hi-C scaffolding


## assembly validation using illumina mate-pairs

## assembly validation using GBS genetic maps

