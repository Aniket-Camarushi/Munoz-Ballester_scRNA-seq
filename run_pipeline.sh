# Complete RNA-Seq pipeline: FASTQ -> Kallisto -> STAR -> GATK -> VCF
# Terminate on error
set -e

echo "[$(date)] Starting pipeline" >> pipeline.log

# 1. Quality Control with FastQC
echo "Executing FastQC..."
mkdir -p FastQC
fastqc ./Data/*.gz -t 20 -o FastQC/

# 2. Kallisto Quantification
echo "Constructing kallisto index..."
mkdir -p Kallisto
kallisto index -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index ./Data/Homo_sapiens.GRCh38.cdna.all.fa

echo "Running kallisto quantifications..."

# Begin with the healthy subject (H5)
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS01 -t 26 --single -l 250 -s 30 ./Data/SRR8668755.fastq.gz &> ./Kallisto/HS01.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS02 -t 26 --single -l 250 -s 30 ./Data/SRR8668756.fastq.gz &> ./Kallisto/HS02.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS03 -t 26 --single -l 250 -s 30 ./Data/SRR8668757.fastq.gz &> ./Kallisto/HS03.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS04 -t 26 --single -l 250 -s 30 ./Data/SRR8668758.fastq.gz &> ./Kallisto/HS04.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/HS05 -t 26 --single -l 250 -s 30 ./Data/SRR8668759.fastq.gz &> ./Kallisto/HS05.log

# Proceed with the cutaneous leishmaniasis (CL) patients
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL08 -t 26 --single -l 250 -s 30 ./Data/SRR8668769.fastq.gz &> ./Kallisto/CL08.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL10 -t 26 --single -l 250 -s 30 ./Data/SRR8668771.fastq.gz &> ./Kallisto/CL10.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL11 -t 26 --single -l 250 -s 30 ./Data/SRR8668772.fastq.gz &> ./Kallisto/CL11.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL12 -t 26 --single -l 250 -s 30 ./Data/SRR8668773.fastq.gz &> ./Kallisto/CL12.log
kallisto quant -i ./Kallisto/Homo_sapiens.GRCh38.cdna.all_my.index  -o ./Kallisto/CL13 -t 26 --single -l 250 -s 30 ./Data/SRR8668774.fastq.gz &> ./Kallisto/CL13.log

# 3. STAR Alignment Setup
echo "Preparing STAR genome index..."
mkdir -p STAR
STAR --runThreadN 26 \
     --runMode genomeGenerate \
     --genomeDir ./STAR/GRCh38_index \
     --genomeFastaFiles ./Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile ./Data/Homo_sapiens.GRCh38.108.gtf \
     --sjdbOverhang 100	


# 4. STAR Alignment
echo "Executing STAR alignments..."

# Set directories and parameters
FASTQ_DIR="./Data"
INDEX="./STAR/GRCh38_index"
THREADS=26
BAM_RAM_LIMIT=30000000000  # 30GB for BAM sorting

# Process each FASTQ file
for fastq in $FASTQ_DIR/*.fastq.gz; do    
	SAMPLE_NAME=$(basename ${fastq%.fastq.gz})        
	
	echo "Processing sample: $SAMPLE_NAME"        

	STAR --genomeDir $INDEX \
	     --readFilesIn $fastq \
	     --readFilesCommand zcat \
	     --outFileNamePrefix ./STAR/${SAMPLE_NAME}_ \
	     --runThreadN $THREADS \
	     --limitBAMsortRAM $BAM_RAM_LIMIT \
	     --outSAMtype BAM SortedByCoordinate \
	     --quantMode TranscriptomeSAM GeneCounts \
	     --outSAMattrRGline ID:${SAMPLE_NAME} SM:${SAMPLE_NAME} LB:lib1 PL:ILLUMINA PU:unit1 \
	     --outFilterType BySJout \ 
	     --outFilterMultimapNmax 20 \ 
	     --alignSJoverhangMin 8 \
	     --alignSJDBoverhangMin 1 \
	     --outFilterMismatchNmax 999 \
	     --outFilterMismatchNoverLmax 0.04 \
	     --alignIntronMin 20 \
	     --alignIntronMax 1000000 \
	     --alignMatesGapMax 1000000 \
	     --sjdbScore 1 \
	     --twopassMode Basic \
	     --outSAMattributes NH HI AS nM NM MD XS
done

# 5. Post-alignment Processing
echo "Processing BAM files..."

for bam in ./STAR/*Aligned.sortedByCoord.out.bam; do    
	base=${bam%Aligned.sortedByCoord.out.bam}    
	picard MarkDuplicates \
       I=$bam \
       O=${base}_dedup.bam \
       M=${base}_dedup_metrics.txt    
	samtools index ${base}_dedup.bam
done


# 6. GATK Variant Calling
echo "Executing GATK variant calling..."
mkdir -p GATK

# Begin by creating a sequence dictionary for the reference
picard CreateSequenceDictionary \
    R=./Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    O=./Data/Homo_sapiens.GRCh38.dna.primary_assembly.dict
      
# Index the reference
samtools faidx ./Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Perform HaplotypeCaller for each sample
for bam in ./STAR/*_dedup.bam; do    
	base=$(basename $bam _dedup.bam)    
	gatk HaplotypeCaller \
	    -R ./Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	    -I $bam \
	    -O ./GATK/${base}.g.vcf.gz \
	    -ERC GVCF

done

# Merge GVCFs (UPDATED: include all samples using a loop)
GVCFs=$(ls ./GATK/*.g.vcf.gz | tr '\n' ' ')  # Gather all GVCFs
gatk CombineGVCFs \
    -R ./Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    $(printf -- "--variant %s " $GVCFs) \
    -O ./GATK/combined.g.vcf.gz

# Joint genotyping
gatk GenotypeGVCFs \
    -R ./Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -V ./GATK/combined.g.vcf.gz \
    -O ./GATK/joint_variants.vcf.gz


# 7. HARD FILTERING (substitute for VQSR)
echo "Applying hard filters..."
gatk VariantFiltration \
    -R ./Data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -V ./GATK/joint_variants.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "GATK_RNAseq_filter" \
    -O ./GATK/filtered_variants.vcf.gz


# 8. MultiQC Summary
echo "Generating MultiQC summary..."
mkdir -p MultiQC
multiqc . -o MultiQC/

echo "Pipeline completed successfully!"
