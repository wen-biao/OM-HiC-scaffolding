
####Improved Dovetail Hi-C Scaffolding workflow
##step1: mapping Hi-C reads to input sequences of Falcon assembly contigs and PBcR assembly contigs, respectively.
  bwa aln $ref ./reads1_1.fq -f reads1_1.sai -t $t
  bwa aln $ref ./reads1_2.fq -f reads1_2.sai -t $t

  bwa aln $ref ./reads2_1.fq -f reads2_1.sai -t $t
  bwa aln $ref ./reads2_2.fq -f reads2_2.sai -t $t

  bwa sampe -a 300000 $ref reads1_1.sai reads1_2.sai reads1_1.fq reads1_2.fq |samblaster --addMateTags |samtools view -bS - -o reads1.bam
  bwa sampe -a 300000 $ref reads2_1.sai reads2_2.sai reads2_1.fq reads2_2.fq |samblaster --addMateTags |samtools view -bS - -o reads2.bam
  samtools merge merged.bam reads1.bam reads2.bam
  samtools sort merged.bam merged.srt
  samtools rmdup merged.srt.bam merged.srt.rmdup.bam
  samtools index merged.srt.rmdup.bam 

##step2: run HiRise scaffolding 
  bash hirise_commands.pa.sh 10 30
 #Note: We used BWA to do reads mapping. Because BWA's mapping quality scale is different from SNAP, we modified the default cufoff of the minimal mapping quality in the HiRise scripts. 

##step3: in-silico Optical Mapping scaffolding
#step3.0: generate in-silico optical maps using Falcon or PBcR assembly contigs, use these maps as reference maps, do the initial alignment between falcon sequences and pbcr in-silico optical maps
  export PERLLIB=/opt/share/software/packages/perl/5.10.1d/lib/5.10.1d:/opt/share/software/packages/perl/5.10.1d/lib/5.10.1d/x86_64-linux-thread-multi
/opt/share/software/packages/perl/5.10.1d/bin/perl $bindir/IrysScaffoldingKSU/scripts/HybridScaffold/scripts/fa2cmap.pl -i ./pbcr.hirise.scaf.bwa.fasta -n BspQI -m 5 -M 50 -o ./ &

  /opt/share/software/packages/perl/5.10.1d/bin/perl hybridScaffold.pl -n ./HiCint/falcon.hirise.scaf.bwa.fasta -b ./HiCint/fa2cmap/pbcr.hirise.scaf.bwa_BspQI_50Kb_5labels.cmap -c xml/hybridScaffold_config_out5.xml -o ./HiCint/1-falconScaf/output -f &

#step3.1: find conflicts and split misassemblies
  perl ./HiCint/OM.find.conflict-regions.for.HiC.improved.pl -k HiCint/falcon.hirise.scaf.bwa_BspQI_0Kb_0labels_key.txt -d HiCint/1-falconScaf/output/ -o HiCint/2-falconBroken/split > HiCint/2-falconBroken/breaks.log

  perl HiCint/OM.split-misassembly.pl -key ./HiCint/falcon.hirise.scaf.bwa_BspQI_0Kb_0labels_key.txt -seq ./HiCint/falcon.hirise.scaf.bwa.fasta -seqMis ./HiCint/2-falconBroken/split/seq.assembly.conflicts.breaks3.flt -aln ./HiCint/1-falconScaf/output -OMmis ./HiCint/2-falconBroken/split/OM.assembly.conflicts.breaks3.flt -outDir ./HiCint/2-falconBroken/split/ > HiCint/2-falconBroken/split.log 


#step3.2: realign to get hybrid consensus maps and do scaffolding
export PERLLIB=/opt/share/software/packages/perl/5.10.1d/lib/5.10.1d:/opt/share/software/packages/perl/5.10.1d/lib/5.10.1d/x86_64-linux-thread-multi
nohup /opt/share/software/packages/perl/5.10.1d/bin/perl hybridScaffold.pl -n HiCint/2-falconBroken/split/seq.assembly.conflicts.broken.fasta -b HiCint/2-falconBroken/split/BNG.assembly.conflicts.broken.cmap -c xml/hybridScaffold_config_out54.xml -o HiCint/2-falconBroken/output -f > HiCint/2-falconBroken/output.log &

source $HOME/.bashrc
cd HiCint
/projects/dep_coupland/grp_nordstrom/bin/IrysScaffoldingKSU/tools/RefAligner -f -ref ./2-falconBroken/output/align_final/step2.hybrid.cmap -i ./2-falconBroken/split/seq.assembly.conflicts.broken_BspQI_0Kb_0labels.cmap -o ./2-falconBroken/bScafSeq_HybridMap -endoutlier 1e-3 -outlier 1e-4 -extend 0 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 0 -resSD 0.75 -mres 0 -A 5 -biaswt 0 -M 1 -Mfast 0 -maxmem 128 -maxthreads 20 -deltaX 9 -deltaY 9 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-9 -BestRef 1 -nosplit 2 -XmapStatWrite ./2-falconBroken/bScafSeq_HybridMap -stdout -stderr -hashgen 5 3 2.4 1.5 0.05 5.0 1 1 1 -hash -hashdelta 50 &

cd 2-falconBroken
nohup perl ../OM.scaffolding.pl -x bScafSeq_HybridMap.xmap -r bScafSeq_HybridMap_r.cmap -q bScafSeq_HybridMap_q.cmap -key split/seq.broken.new.key.txt -outdir scafSeq -seq split/seq.assembly.conflicts.broken.fasta -sID xcaf >scaf.log&
