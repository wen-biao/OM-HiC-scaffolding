Requirement: Please make sure that the system Perl is running version 5.10 or version 5.14.

Please make sure that the following directory structure is kept, and that you have full permission. 
top directory:
rerun_coassembly_pipeline_021615.pl 
hybridScaffold_v1.pl
hybridScaffold_v2.pl
RefAligner	- note: an executable

scripts/
AssignAlignType.pl  
calc_cmap_stats.pl  
calc_xmap_stats.pl  
fa2cmap.pl  
fa_key_convert.pl  
findBNG.pl  
findNGS.pl  
MergeNGS_BN.pl
perl5/

subdirectories:
input_data/
mar3_NA12878.scf.fasta	- note: first PacBio sequence assembly (Celera)
primary_tigs_c_copy.fa	- note: second PacBio sequence assembly (Falcon)
EXP_REFINEFINAL1.cmap	- note: BioNano assembly

xml/
hybridScaffold_config_v1.xml	- note: put the absolute path of the RefAligner executable file in the refaligner field
hybridScaffold_config_v2.xml	- note: put the absolute path of the RefAligner executable file in the refaligner field

Before running the pipeline, be sure to change refaligner field to the absolute path of RefAligner executable. 
To run the pipelines, please run the rerun_coassembly_pipeline_020615.pl script using the command
perl rerun_coassembly_pipeline_021615.pl

There are two output directories v1 and v2, corresponding to the results of the two rounds of hybrid scaffolding: first PacBio sequence assembly
with BioNano assembly, second v1 hybrid scaffolds with second PacBio assembly.

output directory:
v1/output_results/
v2/output_results/

After an hybrid scaffold run is finished, it should have the following files
<input_file1>_<input_file2>_BNGcontigs_HYBRID_SCAFFOLD.xmap  - note: alignment between the input file 1 fragments and the hybrid scaffolds
<input_file1>_<input_file2>_BNGcontigs_HYBRID_SCAFFOLD_q.cmap
<input_file1>_<input_file2>_BNGcontigs_HYBRID_SCAFFOLD_r.cmap
<input_file1>_<input_file2>_NGScontigs_HYBRID_SCAFFOLD.xmap	- note: alignment between the input file 2 fragments and the hybrid scaffolds
<input_file1>_<input_file2>_NGScontigs_HYBRID_SCAFFOLD_q.cmap
<input_file1>_<input_file2>_NGScontigs_HYBRID_SCAFFOLD_r.cmap
hybridScaffold_config_v1.xml
output_results_log.txt	- note: contains a summary of the hybrid scaffold run
align0/
align1/
align_final/
assignAlignType/
fa2cmap/
mergeNGS_BN/


The file output_results_log.txt should contain the running messages, as well as some summary statistics. 

Note, if you choose to try the beta pipeline for rendering fastas from hybrid scaffolds (instead of the Irsyview GUI), note that 
we have used an updated version of RefAligner (called RefAligner2 in this repository as called by "replace_contigs_for_final_scaffold.sh").
For best results in NEW studies we recommend using RefAligner2 for ALL scaffolding steps.  

