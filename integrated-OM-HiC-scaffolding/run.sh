####Integrated Optical-Mapping and Dovetail-HiC Scaffolding

## step 0: find misassembled regions using Hi-C data
mkdir 0-falconHiC 0-pbcrHiC
cd 0-falconHiC 
bwa index falcon.contigs.fasta  ##copy your falcon assembly contigs here
mkdir HiCMap shortMap

cd HiCMap
#ln -s xx.fq ./
ctg1=../falcon.contigs.fasta 
  bwa aln $ctg1 ./reads_1.fq -f reads1_1.sai -t 20
  bwa aln $ctg1 ./reads_2.fq -f reads1_2.sai -t 20
  bwa sampe -a 300000 $ctg1 reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 

cd ../shortMap
#ln -s xx.fq ./
  bwa aln $ctg1 ./reads_1.fq -f reads1_1.sai -t 20
  bwa aln $ctg1 ./reads_2.fq -f reads1_2.sai -t 20
  bwa sampe -a 300000 $ctg1 reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 

#  run HiRise scaffolding
cd ../scaffolding
run hirise_commands.pa.sh ../HiCMap/reads.srt.rmdup.bam ../falcon.contigs.fasta ../shortMap/reads.srt.rmdup.bam 10 30

cd ../../0-pbcrHiC
bwa index pbcr.contigs.fasta  ##copy your pbcr assembly contigs here

mkdir HiCMap shortMap
cd HiCMap
#ln -s xx.fq ./
ctg2=../pbcr.contigs.fasta 
  bwa aln $ctg2 ./reads_1.fq -f reads1_1.sai -t 20
  bwa aln $ctg2 ./reads_2.fq -f reads1_2.sai -t 20
  bwa sampe -a 300000 $ctg2 reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 

cd ../shortMap
#ln -s xx.fq ./
  bwa aln $ctg2 ./reads_1.fq -f reads1_1.sai -t 20
  bwa aln $ctg2 ./reads_2.fq -f reads1_2.sai -t 20
  bwa sampe -a 300000 $ctg2 reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 

#  run HiRise scaffolding
cd ../scaffolding
run hirise_commands.pa.sh ../HiCMap/reads.srt.rmdup.bam ../pbcr.contigs.fasta ../shortMap/reads.srt.rmdup.bam 10 30

cd ../../
mkdir 1-falcon 1-pbcr 2-falconBroken 2-pbcrBroken
perl ./hirise.breaks.table.to.bed.pl 0-falconHiC/scaffolding/hirise_iter_broken_3.gapclosed.table 1-falcon/falcon.hirise.breaks.bed
perl ./hirise.breaks.table.to.bed.pl 0-pbcrHiC/scaffolding/hirise_iter_broken_3.gapclosed.table 1-pbcr/pbcr.hirise.breaks.bed

## step 1:  align CMAPs to sequence assemblies with low and high initial alignment p-value, respectively

#low p-value Falcon
perl hybridScaffold.pl -n Aalpina.falcon.assembly.fasta -b Aalpina.denovo.cmap -c xml/hybridScaffold_config_out5.xml -o OM-HiC_int/1-falcon/output5 -f > OM-HiC_int/1-falcon/output5.log 
#high p-value Falcon
perl hybridScaffold.pl -n Aalpina.falcon.assembly.fasta -b Aalpina.denovo.cmap -c xml/hybridScaffold_config_out7.xml -o OM-HiC_int/1-falcon/output7 -f > OM-HiC_int/1-falcon/output7.log 

#low p-value PBcR
perl hybridScaffold.pl -n Aalpina.pbcr.assembly.fasta -b Aalpina.denovo.cmap -c xml/hybridScaffold_config_out5.xml -o OM-HiC_int/1-pbcr/output5 -f > OM-HiC_int/1-pbcr/output5.log 
#high p-value PBcR
perl hybridScaffold.pl -n Aalpina.pbcr.assembly.fasta -b Aalpina.denovo.cmap -c xml/hybridScaffold_config_out7.xml -o OM-HiC_int/1-pbcr/output7 -f > OM-HiC_int/1-pbcr/output7.log 

## step 2: find breaks

cd OM-HiC_int
perl OM.find.conflict-regions.pl -k ../Aalpina.falcon.assembly_BspQI_0Kb_0labels_key.txt -d ./1-falcon/output5/ -o ./1-falcon/output5/breaks >1-falcon/output5/breaks.log 
perl OM.find.conflict-regions.pl -k ../Aalpina.falcon.assembly_BspQI_0Kb_0labels_key.txt -d ./1-falcon/output7/ -o ./1-falcon/output7/breaks >1-falcon/output7/breaks.log 

perl OM.find.conflict-regions.pl -k ../Aalpina.pbcr.assembly_BspQI_0Kb_0labels_key.txt -d ./1-pbcr/output7/ -o ./1-pbcr/output7/breaks >1-pbcr/output7/breaks.log 
perl OM.find.conflict-regions.pl -k ../Aalpina.pbcr.assembly_BspQI_0Kb_0labels_key.txt -d ./1-pbcr/output5/ -o ./1-pbcr/output5/breaks >1-pbcr/output5/breaks.log 


## step 3: merge breaks from low and high pvalue runnings
perl OM.merge.low-high.breaks.pl -a ./1-falcon/output5/breaks/seq.assembly.conflicts.breaks -b ./1-falcon/output7/breaks/seq.assembly.conflicts.breaks -c ./1-falcon/output5/breaks/OM.assembly.conflicts.breaks -d ./1-falcon/output7/breaks/OM.assembly.conflicts.breaks -o ./2-falconBroken/mergeBreaks >2-falconBroken/mergeBreaks.log 
perl OM.merge.low-high.breaks.pl -a ./1-pbcr/output5/breaks/seq.assembly.conflicts.breaks -b ./1-pbcr/output7/breaks/seq.assembly.conflicts.breaks -c ./1-pbcr/output5/breaks/OM.assembly.conflicts.breaks -d ./1-pbcr/output7/breaks/OM.assembly.conflicts.breaks -o ./2-pbcrBroken/mergeBreaks >2-pbcrBroken/mergeBreaks.log 


## step 4: windowBed btx OM-breaks, HR-breaks
perl ./OM.conflicts.add.HiRise.pl -s 2-falconBroken/mergeBreaks/seq.breaks4 -m 1-falcon/output5/align1/align1_r.cmap -r ./falcon.hirise.breaks.bed  -o 2-falconBroken/mergeBreaks >>2-falconBroken/mergeBreaks.log 
perl ./OM.conflicts.add.HiRise.pl -s 2-pbcrBroken/mergeBreaks/seq.breaks4 -m 1-pbcr/output5/align1/align1_r.cmap -r ./pbcr.hirise.breaks.bed  -o 2-pbcrBroken/mergeBreaks >>2-pbcrBroken/mergeBreaks.log 


## step 5: check whether conflict CMAP region is fullly aligned by another seq CMAP
perl ./OM.conflicts.falcon-pbcr.cross-check.pl -c 2-falconBroken/mergeBreaks/OM.breaks4 -s 2-falconBroken/mergeBreaks/seq.breaks4.HR.bed -x 1-pbcr/output7/align1/align1.xmap -k ../Aalpina.pbcr.assembly_BspQI_0Kb_0labels_key.txt -o 2-falconBroken/mergeBreaks >> 2-falconBroken/mergeBreaks.log

perl ./OM.conflicts.falcon-pbcr.cross-check.pl -c 2-pbcrBroken/mergeBreaks/OM.breaks4 -s 2-pbcrBroken/mergeBreaks/seq.breaks4.HR.bed -x 1-falcon/output7/align1/align1.xmap -k ../Aalpina.falcon.assembly_BspQI_0Kb_0labels_key.txt -o 2-pbcrBroken/mergeBreaks >> 2-pbcrBroken/mergeBreaks.log

## step 6: filter Seq and OM breaks based on 2-CMAP, HR-breaks and falcon-pbcr cross-check
perl ./OM.breaks.filter.pl -c 2-falconBroken/mergeBreaks/OM.breaks4.OMGoodAln.bed -s 2-falconBroken/mergeBreaks/seq.breaks4.HR.OMGoodAln.bed -o 2-falconBroken/mergeBreaks >2-falconBroken/mergeBreaks.log
perl ./OM.breaks.filter.pl -c 2-pbcrBroken/mergeBreaks/OM.breaks4.OMGoodAln.bed -s 2-pbcrBroken/mergeBreaks/seq.breaks4.HR.OMGoodAln.bed -o 2-pbcrBroken/mergeBreaks >2-pbcrBroken/mergeBreaks.log

## step 7: split the Seq and OM 
perl ./OM.breaks.ctg.pl -s  2-falconBroken/mergeBreaks/seq.breaks4.flt.manual.check.bed -f ../Aalpina.falcon.assembly.fasta -m 2-falconBroken/mergeBreaks/OM.breaks4.flt.manual.check.bed -c 1-falcon/output5/align1/align1_q.cmap -o 2-falconBroken/split > 2-falconBroken/split.log

perl ./OM.breaks.ctg.pl -s  2-pbcrBroken/mergeBreaks/seq.breaks4.flt.manual.check.bed -f ../Aalpina.pbcr.assembly.fasta -m 2-pbcrBroken/mergeBreaks/OM.breaks4.flt.manual.check.bed -c 1-pbcr/output5/align1/align1_q.cmap -o 2-pbcrBroken/split > 2-pbcrBroken/split.log

## step 8: realign and ## step 9: do hybrid scaffolding 
#lower pvalue for assigning conflict alignment, to avoid wrong conflict due to partial/multiple alignment
cd ../
perl hybridScaffold.pl -n OM-HiC_int/2-falconBroken/split/seq.broken.fasta -b  OM-HiC_int/2-falconBroken/split/OM.broken.cmap -c xml/hybridScaffold_config_out54.xml -o  OM-HiC_int/2-falconBroken/output54 -f >  OM-HiC_int/2-falconBroken/output54.log 
perl hybridScaffold.pl -n  OM-HiC_int/2-pbcrBroken/split/seq.broken.fasta -b  OM-HiC_int/2-pbcrBroken/split/OM.broken.cmap -c xml/hybridScaffold_config_out54.xml -o  OM-HiC_int/2-pbcrBroken/output54 -f >  OM-HiC_int/2-pbcrBroken/output54.log 

cd OM-HiC_int/2-falconBroken/output54
RefAligner -f -ref ./align_final/step2.hybrid.cmap -i ../split/seq.broken_BspQI_0Kb_0labels.cmap -o brkCtg_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite brkCtg_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50
 perl ../../OM.scaffolding.pl -x brkCtg_HybridMap.xmap -r brkCtg_HybridMap_r.cmap -q brkCtg_HybridMap_q.cmap -key ../split/seq.broken_BspQI_0Kb_0labels_key.txt -outdir scafSeq -seq ../split/seq.broken.fasta -sID scaf > scaffolding.log&

cd ../../../OM-HiC_int/2-pbcrBroken/output54
RefAligner -f -ref ./align_final/step2.hybrid.cmap -i ../split/seq.broken_BspQI_0Kb_0labels.cmap -o brkCtg_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite brkCtg_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50
 perl ../../OM.scaffolding.pl -x brkCtg_HybridMap.xmap -r brkCtg_HybridMap_r.cmap -q brkCtg_HybridMap_q.cmap -key ../split/seq.broken_BspQI_0Kb_0labels_key.txt -outdir scafSeq -seq ../split/seq.broken.fasta -sID scaf > scaffolding.log&

echo "start 2nd-round hybrid scaffolding"
## step 10: align the pbcr-hybrid cmap to falcon-hybrid scaffolds
cd ../../
mkdir 3-falconOMint
cp 2-pbcrBroken/output54/mergeNGS_BN/step2.hybrid.cmap ./3-falconOMint/pbcr.hybrid.BN.naive.cmap
grep -v '#' 2-pbcrBroken/output54/mergeNGS_BN/step1.BN.naive.cmap >>./3-falconOMint/pbcr.hybrid.BN.naive.cmap

cp 2-falconBroken/output54/scafSeq/hybrid.map.scaffolds.fasta 3-falconOMint/hybrid.map.scaffolds.fasta
cd ../

perl hybridScaffold.pl -n OM-HiC_int/3-falconOMint/hybrid.map.scaffolds.fasta -b OM-HiC_int/3-falconOMint/pbcr.hybrid.BN.naive.cmap -c xml/hybridScaffold_config_out5.xml -o OM-HiC_int/3-falconOMint/output5 -f >OM-HiC_int/3-falconOMint/output5.log 


## step 11: find the conflicting alignment and break the scaffold sequences and hybrid cmap
cd OM-HiC_int
perl ./OM.find.conflict-regions.pl -k 3-falconOMint/hybrid.map.scaffolds_BspQI_0Kb_0labels_key.txt -d 3-falconOMint/output/ -o 3-falconOMint/output/breaks
###check the conflicts 
perl ./OM.split-misassembly.pl  -k ./3-falconOMint/hybrid.map.scaffolds_BspQI_0Kb_0labels_key.txt -f ./3-falconOMint/hybrid.map.scaffolds.fasta  -s ./3-falconOMint/output/breaks/seq.assembly.conflicts.mrg.breaks -m ./3-falconOMint/output/breaks/OM.assembly.conflicts.mrg.breaks -o ./3-falconOMint/output -d ./3-falconOMint/output/breaks > 3-falconOMint/output/breaks/split.log

## step 12: further hybrid scaffolding 
perl hybridScaffold.pl -n OM-HiC_int/3-falconOMint/output/breaks/seq.assembly.conflicts.broken.fasta -bOM-HiC_int/3-falconOMint/output/breaks/BNG.assembly.conflicts.broken.cmap -c xml/hybridScaffold_config_out54.xml -o OM-HiC_int/3-falconOMint/outSplitRealign -f > OM-HiC_int/3-falconOMint/outSplitRealign.log
cd OM-HiC_int/3-falconOMint/outSplitRealign

RefAligner -f -ref align_final/step2.hybrid.cmap -i ../output/breaks/seq.assembly.conflicts.broken_BspQI_0Kb_0labels.cmap -o bScafSeq_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite bScafSeq_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 
perl ../../OM.scaffolding.pl -x bScafSeq_HybridMap.xmap -r bScafSeq_HybridMap_r.cmap -q bScafSeq_HybridMap_q.cmap -key ../output/breaks/seq.broken.new.key.txt -outdir scafSeq -seq ../output/breaks/seq.assembly.conflicts.broken.fasta -sID xcaf > scaffolding.log

###intergrate Hi-C data
# step 13: run HiRise to do assembly scaffolding using Dovetail Hi-C data
cd ../../../
mkdir 4-HiC
cd 4-HiC; mkdir HiCMap shortMap scaffolding
cp ../../3-falconOMint/hybrid.map.scaffolds.fasta ./

# step 14: map Dovetail Hi-C reads and illumina reads from short fragments library to the above assembled scaffolds
bwa index hybrid.map.scaffolds.fasta
cd HiCMap
#ln -s xx.fq ./
ref=../hybrid.map.scaffolds.fasta
  bwa aln $ref ./reads_1.fq -f reads1_1.sai -t 20
  bwa aln $ref ./reads_2.fq -f reads1_2.sai -t 20
  bwa sampe -a 300000 $ref reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 

cd ../shortMap
#ln -s xx.fq ./
  bwa aln $ref ./reads_1.fq -f reads1_1.sai -t 20
  bwa aln $ref ./reads_2.fq -f reads1_2.sai -t 20
  bwa sampe -a 300000 $ref reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 

# step 15: run HiRise scaffolding
cd ../scaffolding
run hirise_commands.pa.sh ../HiCMap/reads.srt.rmdup.bam ../hybrid.map.scaffolds.fasta ../shortMap/reads.srt.rmdup.bam 10 30

# step 16: interative OM scaffolding
cd ../../
mkdir 5-OMiter
cp 4-HiC/scaffolding/hirise_iter_broken_3.gapclosed.fasta ./5-OMiter/scaffolds.fasta
cp 3-OMint/outSplitRealign/step2.hybrid.cmap ./5-OMinter/OMint.hybrid.BN.naive.cmap
grep -v '#' 3-OMint/outSplitRealign/mergeNGS_BN/step1.BN.naive.cmap >> ./5-OMinter/OMint.hybrid.BN.naive.cmap
cd ../
perl hybridScaffold.pl -n OM-HiC_int/5-OMiter/hybrid.map.scaffolds.fasta -b OM-HiC_int/5-OMiter/OMint.hybrid.BN.naive.cmap -c xml/hybridScaffold_config_out54.xml -o OM-HiC_int/5-OMiter/output54 -f >OM-HiC_int/5-OMiter/output54.log 

cd OM-HiC_int/
perl ./OM.find.conflict-regions.pl -k 5-OMiter/scaffolds_BspQI_0Kb_0labels_key.txt -d 5-OMiter/output54/ -o 5-OMiter/output54/breaks

perl ./OM.split-misassembly.pl  -k 5-OMiter/scaffolds_BspQI_0Kb_0labels_key.txt -f 5-OMiter/scaffolds.fasta  -s 5-OMiter/output54/breaks/seq.assembly.conflicts.mrg.breaks -m 5-OMiter/output54/breaks/OM.assembly.conflicts.mrg.breaks -o 5-OMiter/output54 -d 5-OMiter/output54/breaks > 5-OMiter/output54/breaks/split.log

cd ../

perl hybridScaffold.pl -n OM-HiC_int/5-OMiter/output54/breaks/seq.assembly.conflicts.broken.fasta -b OM-HiC_int/5-OMiter/OMint.hybrid.BN.naive.cmap -c xml/hybridScaffold_config_out55.xml -o OM-HiC_int/5-OMiter/output55 -f >OM-HiC_int/5-OMiter/output55.log

cd OM-HiC_int/5-OMiter/output55
RefAligner -f -ref align_final/step2.hybrid.cmap -i ../output54/breaks/seq.assembly.conflicts.broken_BspQI_0Kb_0labels.cmap -o bScafSeq_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite bScafSeq_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50
perl ../../OM.scaffolding.pl -x bScafSeq_HybridMap.xmap -r bScafSeq_HybridMap_r.cmap -q bScafSeq_HybridMap_q.cmap -key ../output54/breaks/seq.broken.new.key.txt -outdir scafSeq -seq ../output54/breaks/seq.assembly.conflicts.broken.fasta -sID xcaf > scaffolding.log


