####Integrated Optical Mapping Scaffolding

## step 1:  align CMAPs to sequence assemblies with low and high initial alignment p-value, respectively
mkdir OMint
mkdir OMint/1-falcon OMint/1-pbcr

perl hybridScaffold.pl -n data/Cp.falcon.contig.fasta -b data/Cp.OM.cmap -c xml/hybridScaffold_config_out5.xml -o ./OMint/1-falcon/out5
perl hybridScaffold.pl -n data/Cp.pbcr.contig.fasta -b data/Cp.OM.cmap -c xml/hybridScaffold_config_out5.xml -o ./OMint/1-pbcr/out5

perl hybridScaffold.pl -n data/Cp.falcon.contig.fasta -b data/Cp.OM.cmap -c xml/hybridScaffold_config_out7.xml -o ./OMint/1-falcon/out7
perl hybridScaffold.pl -n data/Cp.pbcr.contig.fasta -b data/Cp.OM.cmap -c xml/hybridScaffold_config_out7.xml -o ./OMint/1-pbcr/out7

## step 2: find conflict alignment regions between sequence assembly and optical consensus maps
cd OMint
 perl ./OM.find.conflict.breaks.pl -k ../data/Cp.falcon.contig_BspQI_0Kb_0labels_key.txt -d ./1-falcon/out5/ -o ./1-falcon/out5/breaks 
 perl ./OM.find.conflict.breaks.pl -k ../data/Cp.falcon.contig_BspQI_0Kb_0labels_key.txt -d ./1-falcon/out7/ -o ./1-falcon/out7/breaks

 perl ./OM.find.conflict.breaks.pl -k ../data/Cp.pbcr.contig_BspQI_0Kb_0labels_key.txt -d ./1-pbcr/out7/ -o ./1-pbcr/out7/breaks 
 perl ./OM.find.conflict.breaks.pl -k ../data/Cp.pbcr.contig_BspQI_0Kb_0labels_key.txt -d ./1-pbcr/out5/ -o ./1-pbcr/out5/breaks

## step 3: Merge conflicting alignment regions found from previous aligments under low and high initial alignment p-value 
 perl ./OM.merge.low-high.breaks.pl -a ./1-falcon/out5/breaks/seq.assembly.conflicts.breaks -b ./1-falcon/out7/breaks/seq.assembly.conflicts.breaks -c ./1-falcon/out5/breaks/OM.assembly.conflicts.breaks -d ./1-falcon/out7/breaks/OM.assembly.conflicts.breaks -o ./2-falconBroken/mergebreaks 

 perl ./OM.merge.low-high.breaks.pl -a ./1-pbcr/out5/breaks/seq.assembly.conflicts.breaks -b ./1-pbcr/out7/breaks/seq.assembly.conflicts.breaks -c ./1-pbcr/out5/breaks/OM.assembly.conflicts.breaks -d ./1-pbcr/out7/breaks/OM.assembly.conflicts.breaks -o ./2-pbcrBroken/mergebreaks 


## step 4: Assembly cross check and define the misassembled regions
mkdir 1_pbcr-falcon_breaks
perl  ./OM.falcon.pbcr.misass.cross-check.pl -p ./2-pbcrBroken/mergebreaks/seq.breaks4 -x ./2-pbcrBroken/mergedbreaks/OM.breaks4 -f ./2-falconBroken/mergebreaks/seq.breaks4 -y ./2-falconBroken/mergebreaks/OM.breaks4 -pX ./1-pbcr/out5/align1/align1.xmap -fX 1-falcon/out5/align1/align1.xmap -pK ../data/Cp.pbcr.contig_BspQI_0Kb_0labels_key.txt -fK ../data/Cp.falcon.contig_BspQI_0Kb_0labels_key.txt -o ./1_pbcr-falcon_breaks/ >1_pbcr-falcon_breaks/pbcr-falcon.brk.log 

## step 5: split the contig sequences based on the result of misassembled regions from previous step
perl ./OM.split-misassembly.pl -k ../data/Cp.falcon.contig_BspQI_0Kb_0labels_key.txt -f ../data/Cp.falcon.contig.fasta -s ./1_pbcr-falcon_breaks/falcon.seq.breaks5 -m ./1_pbcr-falcon_breaks/falcon.BNG.breaks5 -o ./1-falcon/out5 -d ./2-falconBroken/ > 2-falconBroken/split.log 

## step 6: Realign and scaffolding
mkdir -p 2-falconBroken 2-pbcrBroken
cd ../

perl hybridScaffold.pl -n OMint/2-falconBroken/seq.assembly.conflicts.broken.fasta -b OMint/2-falconBroken/BNG.assembly.conflicts.broken.cmap -c ./xml/hybridScaffold_config_out54.xml -o OMint/2-falconBroken/output1 -f >OMint/2-falconBroken/out1.log 

perl hybridScaffold.pl -n OMint/2-pbcrBroken/seq.assembly.conflicts.broken.fasta -b OMint/2-pbcrBroken/BNG.assembly.conflicts.broken.cmap -c ./xml/hybridScaffold_config_out54.xml -o OMint/2-pbcrBroken/output1 -f > OMint/2-pbcrBroken/out1.log 

cd OMint/
RefAligner -f -ref ./2-falconBroken/output1/align_final/step2.hybrid.cmap -i ./2-falconBroken/seq.assembly.conflicts.broken_BspQI_0Kb_0labels.cmap -o 2-falconBroken/bCtgSeq_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite 2-falconBroken/bCtgSeq_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 
## scaffolding
perl ./OM.scaffolding.pl -x 2-falconBroken/bCtgSeq_HybridMap.xmap -r 2-falconBroken/bCtgSeq_HybridMap_r.cmap -q 2-falconBroken/bCtgSeq_HybridMap_q.cmap -key 2-falconBroken//seq.broken.new.key.txt -outdir 2-falconBroken/scafSeq -seq 2-falconBroken/seq.assembly.conflicts.broken.fasta -sID scaf >2-falconBroken/scaffolding.log

#pbcr 
RefAligner -f -ref ./2-pbcrBroken/output1/align_final/step2.hybrid.cmap -i ./2-pbcrBroken/seq.assembly.conflicts.broken_BspQI_0Kb_0labels.cmap -o ./2-pbcrBroken/bCtgSeq_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite ./2-pbcrBroken/bCtgSeq_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 
## scaffolding
perl ./OM.scaffolding.pl -x 2-pbcrBroken/bCtgSeq_HybridMap.xmap -r 2-pbcrBroken/bCtgSeq_HybridMap_r.cmap -q 2-pbcrBroken/bCtgSeq_HybridMap_q.cmap -key 2-pbcrBroken/seq.broken.new.key.txt -outdir 2-pbcrBroken/scafSeq -seq 2-pbcrBroken/seq.assembly.conflicts.broken.fasta -sID scaf >2-pbcrBroken/scaffolding.log

## step 7: falcon-scaffolds align with pbcr-hybridMap 
mkdir -p 3-pbcrHmap-falconScaf
cp 2-pbcrBroken/output1/mergeNGS_BN/step2.hybrid.cmap 3-pbcrHmap-falconScaf/pbcr.hybrid.BN.naive.cmap
grep -v '#' 2-pbcrBroken/output1/mergeNGS_BN/step1.BN.naive.cmap >>3-pbcrHmap-falconScaf/pbcr.hybrid.BN.naive.cmap
cp 2-falconBroken/scafSeq/hybrid.map.scaffolds.fasta ./3-pbcrHmap-falconScaf/

cd ../
perl hybridScaffold.pl -n OMint/3-pbcrHmap-falconScaf/hybrid.map.scaffolds.fasta -b OMint/3-pbcrHmap-falconScaf/pbcr.hybrid.BN.naive.cmap -c xml/hybridScaffold_config_out5.xml -o OMint/3-pbcrHmap-falconScaf/output -f >OMint/3-pbcrHmap-falconScaf/output.log

## step 8: find the conflicting alignment and break the scaffold seq and v1.hybrid-cmap
cd OMint
perl ./OM.find.conflict-regions.pl -k 3-pbcrHmap-falconScaf/hybrid.map.scaffolds_BspQI_0Kb_0labels_key.txt -d 3-pbcrHmap-falconScaf/output/ -o 3-pbcrHmap-falconScaf/output/breaks
perl ./OM.split-misassembly.pl -k ./3-pbcrHmap-falconScaf/hybrid.map.scaffolds_BspQI_0Kb_0labels_key.txt -f ./3-pbcrHmap-falconScaf/hybrid.map.scaffolds.fasta -s ./3-pbcrHmap-falconScaf/output/breaks/seq.assembly.conflicts.mrg.breaks -m ./3-pbcrHmap-falconScaf/output/breaks/OM.assembly.conflicts.mrg.breaks -o ./3-pbcrHmap-falconScaf/output -d ./3-pbcrHmap-falconScaf/output/breaks > 3-pbcrHmap-falconScaf/output/breaks/split.log

## step 9: 2nd round hybrid scaffolding 
cd ../
 perl hybridScaffold.pl -n OMint/3-pbcrHmap-falconScaf/output/breaks/seq.assembly.conflicts.broken.fasta -b OMint/3-pbcrHmap-falconScaf/output/breaks/BNG.assembly.conflicts.broken.cmap -c xml/hybridScaffold_config_out54.xml -o OMint/3-pbcrHmap-falconScaf/outSplitRealign -f > OMint/3-pbcrHmap-falconScaf/outSplitRealign.log

cd OMint/3-pbcrHmap-falconScaf/outSplitRealign
RefAligner -f -ref align_final/step2.hybrid.cmap -i ../output/breaks/seq.assembly.conflicts.broken_BspQI_0Kb_0labels.cmap -o bScafSeq_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite bScafSeq_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 

perl ../../OM.scaffolding.pl -x bScafSeq_HybridMap.xmap -r bScafSeq_HybridMap_r.cmap -q bScafSeq_HybridMap_q.cmap -key ../output/breaks/seq.broken.new.key.txt -outdir scafSeq -seq ../output/breaks/seq.assembly.conflicts.broken.fasta -sID xcaf > scaffolding.log

