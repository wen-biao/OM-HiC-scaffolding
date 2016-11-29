# Improved Dovetail Hi-C Scaffolding workflow

Working directory: $your_working_directory/OM-HiC-scaffolding/HiCint

The final scaffolds are in the file "HiCint/2-falconBroken/scafSeq/hybrid.map.scaffolds.fasta"

This workflow mainly contains steps as below:

1) mapping Hi-C reads to input sequences of Falcon assembly contigs and PBcR assembly contigs, respectively.

2) run HiRise scaffolding for Falcon contigs and PBcR contigs

3) do in-silico Optical Mapping scaffolding using scaffolds from Falcon contigs' scaffolding and in-silico cmap of scaffold sequences from HiRise scaffolding of PBcR contigs

