merSize=16
mhap=-k 16 --num-hashes 512 --num-min-matches 3 --threshold 0.04 --weighted

useGrid=0
scriptOnGrid=0
gridEngine=LSF

assembleCoverage=25

ovlMemory=32
ovlStoreMemory=32000
threads=32
ovlConcurrency=1
cnsConcurrency=8
merylThreads=32
merylMemory=32000
ovlRefBlockSize=20000
frgCorrThreads = 16
frgCorrBatchSize = 100000
ovlCorrBatchSize = 100000


sgeScript = -pe threads 1
sgeConsensus = -pe threads 8
sgeOverlap = -pe threads 15 –l mem=2GB
sgeCorrection = -pe threads 15 –l mem=2GB
sgeFragmentCorrection = -pe threads 16 –l mem=2GB
sgeOverlapCorrection = -pe threads 1 –l mem=16GB
