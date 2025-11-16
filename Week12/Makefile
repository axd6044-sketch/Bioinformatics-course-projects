# Aim of makefile: calling variants for normal and tumor samples 
# things this makefile does;
# 1. Download reference genome and extract chromosome 17
# 2. Download BAM files for normal and tumor samples and extract reads for region chr17	:7670000-7680000
# 3. Call variants using bcftools for both normal and tumor samples
# 4. Merge and intersect VCF files
# 5. Compare tumor variants to gold standard DeepVariant calls

SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# Define variables
REF_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam
TUMOR_BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam
CHR=chr17
REGION=chr17:7670000-7680000
REF=refs/${CHR}.fasta
BAM=bam/TP53-N-D.bam
TUMOR_BAM=bam/TP53-T.bam
VCF=vcf/
VCF_N=vcf/TP53-N-D-variants.vcf.gz
VCF_T=vcf/TP53-T-D-variants.vcf.gz
VCF_GS=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz
INTERSECT=vcf/intersected_variants/

# Help target to explain the makefile usage
help:
	@echo "Makefile for variant calling in normal and tumor samples"
	@echo ""
	@echo "Available targets:"
	@echo "  getref               - Download and prepare the reference genome"
	@echo "  getbam               - Download and extract alignments for normal sample"
	@echo "  callvariants        - Call variants for normal sample"
	@echo "  gettumorbam         - Download and extract alignments for tumor sample"
	@echo "  calltumorvariants   - Call variants for tumor sample"
	@echo "  bcftools_merge      - Merge VCF files from normal and tumor samples"
	@echo "  bcftools_intersect  - Intersect VCF files from normal and tumor samples"
	@echo "  compare_to_gs       - Compare tumor variants to gold standard DeepVariant calls"
	@echo "  alignments		 - Run alignment extraction for both normal and tumor samples"
	@echo "  variants		 - Run variant calling and merging for both normal and tumor samples"
	@echo "  all                 - Run all steps in order"


## Getting the reference genome 
#URL of the reference genome
REF_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz

#Download the reference genome
REF=refs/GRCh38.fasta.gz

#The chromosome of interest.
CHR=chr17

#usage: make getref 
getref:
	#Create the reference directory
	mkdir -p refs/
	#Download the reference genome
	curl -L ${REF_URL} > ${REF}
	#Index the reference genome for use with samtools.
	samtools faidx ${REF}
	ls -lh ${REF}
	#Extract the chromosome of interest
	samtools faidx ${REF} ${CHR} > refs/${CHR}.fasta

#alignments for the region of interest
#usage: make getbam
getbam:
	#Create the bam directory
	mkdir -p bam/
	# Extract the reads for the region of interest.
	samtools view -b ${BAM_URL} ${REGION} > ${BAM}
	# Index the BAM file.
	samtools index ${BAM}
	ls -lh ${BAM}
	# Get statistics on the BAM file.
	samtools flagstat ${BAM}

#call variants for the region of interest using bcftools
#usage: make callvariants
callvariants:
	#Create the output directory
	mkdir -p ${VCF}
	#Call variants using bcftools mpileup and bcftools call
	bcftools mpileup -f ${REF} ${BAM} | bcftools call -mv -Oz -o ${VCF_N}
	#Index the VCF file
	bcftools index ${VCF_N}
	#View the VCF file
	bcftools view ${VCF_N}
	#Get statistics on the VCF file
	bcftools stats ${VCF_N} > ${VCF_N}.stats.txt
	cat ${VCF_N}.stats.txt
	bcftools view ${VCF_N} | grep -v '^#' | wc -l
	bcftools index --tbi ${VCF_N}


#repeat the process for the tumor sample
#usage: make gettumorbam
gettumorbam:
	#Create the bam directory
	mkdir -p bam/
	# Extract the reads for the region of interest.
	samtools view -b ${TUMOR_BAM_URL} ${REGION} > ${TUMOR_BAM}
	# Index the tumor BAM file.
	samtools index ${TUMOR_BAM}
	ls -lh ${TUMOR_BAM}
	# Get statistics on the tumor BAM file.
	samtools flagstat ${TUMOR_BAM}

#call variants for the tumor sample
#usage: make calltumorvariants
calltumorvariants:
	#Create the output directory
	mkdir -p ${VCF}
	#Call variants using bcftools mpileup and bcftools call
	bcftools mpileup -f ${REF} ${TUMOR_BAM} | bcftools call -mv -Oz -o ${VCF_T}
	#Index the VCF file
	bcftools index ${VCF_T}
	#View the VCF file
	bcftools view ${VCF_T}
	#Get statistics on the VCF file
	bcftools stats ${VCF_T} > ${VCF_T}.stats.txt
	cat ${VCF_T}.stats.txt
	bcftools view ${VCF_T} | grep -v '^#' | wc -l
	bcftools index --tbi ${VCF_T}

#merge all vcf files 
bcftools_merge: 
	#Merge the normal and tumor VCF files
	bcftools merge -Oz -o vcf/merged_variants.vcf.gz ${VCF_N} ${VCF_T}
	#Index the merged VCF file
	bcftools index vcf/merged_variants.vcf.gz
	#Get statistics on the merged VCF file
	bcftools stats vcf/merged_variants.vcf.gz > vcf/merged_variants.stats.txt
	cat vcf/merged_variants.stats.txt
	bcftools view vcf/merged_variants.vcf.gz | grep -v '^#' | wc -l


#intersect the normal and tumor variants
bcftools_intersect: 
	mkdir -p ${INTERSECT}
	#Intersect the normal and tumor variants using bcftools isec
	bcftools isec -p ${INTERSECT} ${VCF_N} ${VCF_T}
	#The number of variants only in normal sample 
	grep -v '^#' vcf/intersected_variants/0000.vcf | wc -l
	#The number of variants only in tumor sample
	grep -v '^#' vcf/intersected_variants/0001.vcf | wc -l
	#The number of variants common in both normal and tumor samples
	grep -v '^#' vcf/intersected_variants/0002.vcf | wc -l
	grep -v '^#' vcf/intersected_variants/0003.vcf | wc -l

#Compare to gold standard DeepVariant calls
VCF_GS=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/Ultima_DeepVariant-somatic-SNV-INDEL_20240822/HG008_edv_AF_recalibration.result.PASS.vcf.gz
compare_to_gs:
	#Download the gold standard VCF file
	curl -L ${VCF_GS} > DeepVariant-HG008.vcf.gz
	#Index the gold standard VCF file
	tabix -p vcf DeepVariant-HG008.vcf.gz
	bcftools view ${VCF_GS} ${REGION} | bcftools stats
	bcftools view ${VCF_T} ${REGION} | bcftools stats
	#Compare the tumor variants to the gold standard using bcftools isec
	bcftools isec -p vcf/gs_comparison/ DeepVariant-HG008.vcf.gz ${VCF_T}
	#The number of variants only in gold standard
	grep -v '^#' vcf/gs_comparison/0000.vcf | wc -l
	#The number of variants only in tumor sample
	grep -v '^#' vcf/gs_comparison/0001.vcf | wc -l
	#The number of variants common in both gold standard and tumor samples
	grep -v '^#' vcf/gs_comparison/0002.vcf | wc -l

alignments: getbam gettumorbam	
variants: callvariants calltumorvariants bcftools_merge
.PHONY: clean help getref getbam callvariants gettumorbam calltumorvariants bcftools_merge bcftools_intersect compare_to_gs all
all: getref getbam callvariants gettumorbam calltumorvariants bcftools_merge bcftools_intersect compare_to_gs alignments variants
