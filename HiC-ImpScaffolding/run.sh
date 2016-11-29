###  Improved Dovetail Hi-C Scaffolding workflow
# Working directory: $your_working_directory/OM-HiC-scaffolding/HiCint
# The final scaffolds are in the file "HiCint/2-falconBroken/scafSeq/hybrid.map.scaffolds.fasta"

## mapping Hi-C reads to input sequences of Falcon assembly contigs and PBcR assembly contigs, respectively.
fCtg=$1
pCtg=$2
read1=$3
read2=$4
sread1=$5
sread2=$6

mkdir -p 0-falconHiC
mkdir -p 0-pbcrHiC

cd 0-falconHiC 
bwa index $fCtg
bwa index $pCtg

mkdir HiCMap shortMap
cd HiCMap
  bwa aln $fCtg $read1 -f reads1_1.sai -t 20
  bwa aln $fCtg $read2 -f reads1_2.sai -t 20
  bwa sampe -a 300000 $fCtg reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 
cd ../shortMap
  bwa aln $fCtg $sread1 -f reads1_1.sai -t 20
  bwa aln $fCtg $sread2 -f reads1_2.sai -t 20
  bwa sampe -a 300000 $fCtg reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 

# run HiRise scaffolding
cd ../scaffolding
run hirise_commands.pa.sh ../HiCMap/reads.srt.rmdup.bam $fCtg ../shortMap/reads.srt.rmdup.bam 10 30

cd ../../0-pbcrHiC
mkdir HiCMap shortMap
cd HiCMap
  bwa aln $pCtg $read1 -f reads1_1.sai -t 20
  bwa aln $pCtg $read2 -f reads1_2.sai -t 20
  bwa sampe -a 300000 $pCtg reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 
cd ../shortMap
  bwa aln $pCtg $sread1 -f reads1_1.sai -t 20
  bwa aln $pCtg $sread2 -f reads1_2.sai -t 20
  bwa sampe -a 300000 $pCtg reads_1.sai reads_2.sai reads_1.fq reads_2.fq |samblaster --addMateTags |samtools view -bS - -o reads.bam
  samtools sort reads.bam reads.srt
  samtools rmdup reads.srt.bam reads.srt.rmdup.bam
  samtools index reads.srt.rmdup.bam 

# run HiRise scaffolding
cd ../scaffolding
run hirise_commands.pa.sh ../HiCMap/reads.srt.rmdup.bam $pCtg ../shortMap/reads.srt.rmdup.bam 10 30

cd ../../

## in-silico Optical Mapping scaffolding
mkdir 1-falconScaf 1-pbcrScaf 2-falconBroken 2-pbcrBroken
cp 0-falconHiC/scaffolding/hirise_iter_broken_3.gapclosed.fasta falcon.hirise.scaf.bwa.fasta
cp 0-pbcrHiC/scaffolding/hirise_iter_broken_3.gapclosed.fasta pbcr.hirise.scaf.bwa.fasta

#generate in-silico optical maps using Falcon or PBcR assembly contigs, use these maps as reference maps, do the initial alignment between falcon sequences and pbcr in-silico optical maps
perl ../scripts/fa2cmap.pl -i pbcr.hirise.scaf.bwa.fasta -n BspQI -m 5 -M 50 -o ./ 

cd ../
perl hybridScaffold.pl -n ./HiCint/falcon.hirise.scaf.bwa.fasta -b ./HiCint/pbcr.hirise.scaf.bwa_BspQI_50Kb_5labels.cmap -c xml/hybridScaffold_config_out5.xml -o ./HiCint/1-falconScaf/output -f 

# find conflicts and split misassemblies
perl ./HiCint/OM.find.conflict-regions.for.HiC.improved.pl -k HiCint/falcon.hirise.scaf.bwa_BspQI_0Kb_0labels_key.txt -d HiCint/1-falconScaf/output/ -o HiCint/2-falconBroken/split > HiCint/2-falconBroken/breaks.log

perl ./HiCint/OM.split-misassembly.pl -key ./HiCint/falcon.hirise.scaf.bwa_BspQI_0Kb_0labels_key.txt -seq ./HiCint/falcon.hirise.scaf.bwa.fasta -seqMis ./HiCint/2-falconBroken/split/seq.assembly.conflicts.breaks3.flt -aln ./HiCint/1-falconScaf/output -OMmis ./HiCint/2-falconBroken/split/OM.assembly.conflicts.breaks3.flt -outDir ./HiCint/2-falconBroken/split/ > HiCint/2-falconBroken/split.log 

# realign to get hybrid consensus maps and do scaffolding
perl hybridScaffold.pl -n HiCint/2-falconBroken/split/seq.assembly.conflicts.broken.fasta -b HiCint/2-falconBroken/split/BNG.assembly.conflicts.broken.cmap -c xml/hybridScaffold_config_out54.xml -o HiCint/2-falconBroken/output -f > HiCint/2-falconBroken/output.log 

cd HiCint
RefAligner -f -ref ./2-falconBroken/output/align_final/step2.hybrid.cmap -i ./2-falconBroken/split/seq.assembly.conflicts.broken_BspQI_0Kb_0labels.cmap -o ./2-falconBroken/bScafSeq_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite ./2-falconBroken/bScafSeq_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 

cd 2-falconBroken
perl ../OM.scaffolding.pl -x bScafSeq_HybridMap.xmap -r bScafSeq_HybridMap_r.cmap -q bScafSeq_HybridMap_q.cmap -key split/seq.broken.new.key.txt -outdir scafSeq -seq split/seq.assembly.conflicts.broken.fasta -sID xcaf >scaf.log

#the final scaffolds are in the file "HiCint/2-falconBroken/scafSeq/hybrid.map.scaffolds.fasta"

