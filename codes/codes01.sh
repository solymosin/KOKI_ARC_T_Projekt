
cd /data/sn/HErik/testosteron_220403/FASTQ_Generation_2022-04-03_20_27_02Z-548161623

mv *ds*/*.gz .
rm -rf *ds*

export PATH=$PATH:/data/tools/FastQC

mkdir qc01

for f in *fastq.gz 
do 
  fastqc -t 38 $f -o qc01
done 

multiqc qc01 . 


for f in *_L001_R1_001.fastq.gz
do
  zcat ${f/'_L001_R1_001.fastq.gz'/'*'} > ../${f/'_L001_R1_001.fastq.gz'/'.fastq'}
done

TRIM='/data/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'

for f in *.fastq
do 
  o=${f/'.fastq'/'_trimmed.fastq'}
  java -jar $TRIM SE -threads 38 \
    $f $o \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:30 \
    MINLEN:50
done

export PATH=/data/tools/STAR-2.7.9a/bin/Linux_x86_64
export PATH=/data/tools/subread-2.0.2-Linux-x86_64/bin:$PATH


idx='/data/dbs/STAR/GRCm39_104'

for f in *_trimmed.fastq
do
  root=${f/'.fastq'/''}
  STAR -- genomeDir $idx \
       -- readFilesIn $f \
       -- outFileNamePrefix 'GRCm39_104_'$root \
       -- outFilterMultimapNmax 1 \
       -- outReadsUnmapped Fastx \
       -- outSAMtype BAM SortedByCoordinate \
       -- twopassMode Basic \
       -- runThreadN 14 
done
 
featureCounts -O \
   -a $idx/Mus_musculus.GRCm39.104.gtf \
   -o featureCounts_GRCm39.104_O.txt \
   GRCm39_104*bam 


