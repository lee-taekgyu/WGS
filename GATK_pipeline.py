# Install SRA data from ncbi
/BiO/Install/sratoolkit.2.9.6-ubuntu64/bin/fastq-dump --gzip --split-3 SRR1002940

# Get first n lines of the file if wanted (if file too big)
head -n 40000 filename.fq > filename.temp.fq

# Run FastQC with downloaded file
/BiO/Install/FastQC_0.10.1/fastqc -t 4 --nogroup SRR1002940.r1.temp.fq

# Trimming using Trimmomatic
java -jar /BiO/Install/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 SRR1002940.r1.temp.fq SRR1002940.r2.temp.fq SRR1002940.r1.trim.fq SRR1002940.r1.unpair.fq SRR1002940.r2.trim.fq SRR1002940.r2.unpair.fq ILLUMINACLIP:/BiO/Install/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:151:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Run FastQC again with trimmed data (just to compare with previous FastQC)
# /BiO/Install/FastQC_0.10.1/fastqc -t 4 --nogroup SRR1002940.r1.trim.fq

# Map reads to reference sequence using BWA mem
bwa mem -t 4 -R ‘@RG\tPL:Illumina\tID:YUHL\tSM:SRR1002940\tLB:HiSeq’ /BiO/Education/WGS/REF/hg19.fa SRR1002940.r1.trim.fq SRR1002940.r2.trim.fq > SRR1002940.sam

# Mark Duplicates (PICARD)
## 1. Duplication tagging
mkdir TEMP_PICARD
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar AddOrReplaceReadGroups TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I=SRR1002940.sam O=SRR1002940_sorted.bam RGID=SRR1002940 RGLB=HiSeq RGPL=Illumina RGPU=unit1 RGSM=SRR1002940 CREATE_INDEX=true

## 2. Remove duplicates
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar MarkDuplicates TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT I=SRR1002940_sorted.bam O=SRR1002940_dedup.sam M=SRR1002940.duplicate_metrics REMOVE_DUPLICATES=true AS=true

## 3. Sorting
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar SortSam TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I=SRR1002940_dedup.sam O=SRR1002940_dedup.bam CREATE_INDEX=true

### Converting from bam to sam (and vice versa)
### samtools view SRR1002940_sorted.bam > test.sam
### samtools view -Sb SRR1002940.sam > test.bam

# Base Quality Score Recalibration
## First pass
### 1. Analyze covariation pattern in sequencing data (find the region for recalibration and store in table)
java –Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_dedup.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf -O SRR1002940_recal_pass1.table
### 2. Recalibration
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ApplyBQSR -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_dedup.bam --bqsr-recal-file SRR1002940_recal_pass1.table -O SRR1002940_recal_pass1.bam

## Second Pass
### Repeat the same exact steps as in the first pass except now use SRR1002940_recal_pass1.bam instead of SRR1002940_dedup.bam as the input. Similarly, the outputs will be SRR1002940_recal_pass2.table and SRR1002940_recal_pass2.bam

# Calling variants for all samples (extract variants from sequencing data with HaplotypeCaller)
java –Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar HaplotypeCaller -R /BiO/Education/WGS/REF/hg19.fa -I SRR1002940_recal_pass2.bam -O SRR1002940.rawVariants.g.vcf -ERC GVCF --standard-min-confidence-threshold-for-calling 20

