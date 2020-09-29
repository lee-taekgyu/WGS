sample=$1

/BiO/Install/FastQC_0.10.1/fastqc -t 4 ${sample}.r1.temp.fq

/BiO/Install/FastQC_0.10.1/fastqc -t 4 ${sample}.r2.temp.fq

java -jar /BiO/Install/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 ${sample}.r1.temp.fq ${sample}.r2.temp.fq ${sample}.r1.trim.fq ${sample}.r1.unpair.fq ${sample}.r2.trim.fq ${sample}.r2.unpair.fq ILLUMINACLIP:/BiO/Install/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:151:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

bwa mem -t 4 -R '@RG\tPL:Illumina\tID:YUHL\tSM:${sample}\tLB:HiSeq' /BiO/Education/WGS/REF/hg19.fa ${sample}.r1.trim.fq ${sample}.r2.trim.fq > ${sample}.sam

java -jar /BiO/Install/picard-tools-2.22.3/picard.jar AddOrReplaceReadGroups TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I= ${sample}.sam O= ${sample}_sorted.bam RGID=${sample} RGLB=HiSeq RGPL=Illumona RGPU=unit1 RGSM=${sample} CREATE_INDEX=true

##process two
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar MarkDuplicates TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT I= ${sample}_sorted.bam O= ${sample}_dedup.sam M= ${sample}.duplicate_metrics REMOVE_DUPLICATES=true AS=true

##process three
java -jar /BiO/Install/picard-tools-2.22.3/picard.jar SortSam TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I= ${sample}_dedup.sam O=${sample}_dedup.bam CREATE_INDEX=true

#Base Quality Score Recalibration (GATK)
##Analyze pattern
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I ${sample}_dedup.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf -O ${sample}_recal_pass1.table

##Recalibration
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ApplyBQSR -R /BiO/Education/WGS/REF/hg19.fa -I ${sample}_dedup.bam --bqsr-recal-file ${sample}_recal_pass1.table -O ${sample}_recal_pass1.bam

##Do it again
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I ${sample}_recal_pass1.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf -O ${sample}_recal_pass2.table

##Do it again
java -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ApplyBQSR -R /BiO/Education/WGS/REF/hg19.fa -I ${sample}_recal_pass1.bam -bqsr ${sample}_recal_pass2.table -O ${sample}_recal_pass2.bam

#Calling Variants - HaplotypeCaller
java -Xmx16g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar HaplotypeCaller -R /BiO/Education/WGS/REF/hg19.fa -I ${sample}_recal_pass2.bam -O ${sample}.rawVariants.g.vcf -ERC GVCF --standard-min-confidence-threshold-for-calling 20

java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar GenotypeGVCFs -R /BiO/Education/WGS/REF/hg19.fa -V ${sample}.rawVariants.g.vcf -O ${sample}_genotype.vcf

#Extracting Variants (SNP.INDEL) - SelectVariants
##SNP
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants -R /BiO/Education/WGS/REF/hg19.fa -V ${sample}_genotype.vcf --select-type-to-include SNP -O ${sample}.rawSNPs.vcf

##INDEL
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants -R /BiO/Education/WGS/REF/hg19.fa -V ${sample}_genotype.vcf --select-type-to-include INDEL -O ${sample}.rawINDELs.vcf

#Filtration - VariantFiltration
##SNP
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration -R /BiO/Education/WGS/REF/hg19.fa -V ${sample}.rawSNPs.vcf -O ${sample}.rawSNPs.filtered.vcf --filter-name "." --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || haplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0"

##INDEL
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration -R /BiO/Education/WGS/REF/hg19.fa -V ${sample}.rawINDELs.vcf -O ${sample}.rawINDELs.filtered.vcf --filter-name "." --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

#Merge Variants Information 
##SortVcf - gatk
java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SortVcf -I ${sample}.rawSNPs.filtered.vcf -I ${sample}.rawINDELs.filtered.vcf -O ${sample}.Filtered_gatk.variant.vcf


#Annotation - Annovar
egrep "^#|PASS" ${sample}.Filtered_gatk.variant.vcf > ${sample}.Filtered.variant.PASS.vcf


perl /BiO/Install/annovar/table_annovar.pl ${sample}.Filtered.variant.PASS.vcf /BiO/Education/WGS/humandb/ -buildver hg19 -out ${sample} -remove -protocol refGene,cytoBand,avsnp138,clinvar_20190305 -operation g,r,f,f -nastring . -vcfinput
