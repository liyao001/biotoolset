#PBS -N {TISSUE}-{ERRCODE}
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -q high
READS1={ERRCODE}_1.fastq.gz
READS2={ERRCODE}_2.fastq.gz
path=/public/home/yl/normal/ebi/{TISSUE}
cd $path
gunzip $READS1
gunzip $READS2
READS1=${READS1/.gz/}
READS2=${READS2/.gz/}
toolPath=/public/home/yl/bin
cd $toolPath
./hisat2 -q -x /public/home/yl/genomes/grch38/genome -1 $path/$READS1 -2 $path/$READS2 -S $path/{EXP}.sam -p 12 --dta-cufflinks
./samtools view -@ 12 -bS $path/{EXP}.sam -o $path/$EXP.bam
./samtools sort -@ 12 $path/{EXP}.bam -o $path/$EXP.sorted.bam
./cufflinks -p 12 -G /public/home/yl/genomes/Homo_sapiens.GRCh381.79.gtf -o $path/{EXP} $path/{EXP}.sorted.bam
