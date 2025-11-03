#This is a makefile for RNA-seq analysis
#Process both single and paired-end reads using pattern-based rules
#In this assignment, you will extend your previous alignment pipeline to include variant calling. 
#Extend your existing Makefile to call variants on a single sample of your choice.

#Use your existing BAM file as input
#Generate a VCF file for this sample
#Follow best practices for variant calling
#Visualize the resulting VCF file data alongside the BAM file

SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

#variables
REF = ref
species ?= Zika
genome ?= GCF_000882815.3 / GCA_000882815.1
srr ?= SRR3191545
READS=reads
fastqcreports ?= $(reads)/fastqc_reports
genome_fa ?= $(ref)/$(species)genome.fa
bam ?= bam
PRJNA ?= PRJNA313294


# ------ No changes in code beyond this point -----

all: design get_genome get_fastq genome_index alignreads bigwig call_variants
	echo ">> All steps completed successfully."

# Display help information
help:
	@echo "Available targets:"
	@echo "  design        - Generate design.csv from BioProject"
	@echo "  get_genome    - Download reference genome"
	@echo "  get_fastq     - Download and process FASTQ files (uses -X 100000 to limit data size)"
	@echo "  genome_index  - Build genome index (supports human or Zika)"
	@echo "  alignreads    - Align reads to genome using BWA (supports single- or paired-end)"
	@echo "  bigwig        - Convert BAM to BedGraph + BigWig for IGV/Genome Browser"
	@echo "  call_variants - Call variants using bcftools"
	@echo ""
	@echo "Usage:"
	@echo "  make design PRJNA=PRJNA313294"
	@echo "  make get_genome species=human GCF=GCF_009914755.1 REF=ref ACC=NC_060948.1"
	@echo "  make get_fastq srr=SRR3191545 READS=reads fastqcreports=reads/fastqc_reports"
	@echo "  make genome_index species=human REF=ref/human_genome.fa"
	@echo "  make alignreads species=human REF=ref/human_genome.fa READS=reads srr=SRR3191545 bam=bam"
	@echo "  make bigwig species=human srr=SRR3191545 bam=bam REF=ref/human_genome.fa"
	@echo "  make call_variants species=human REF=ref/human_genome.fa bam=bam srr=SRR3191545 vcf=vcf"



# Generate design.csv from BioProject
# Usage: make design PRJNA=PRJNA313294
make design:
	bio search -H --csv ${PRJNA} > design.csv
	cat design.csv

# Download reference genome
# Usage:
#   make get_genome species=human  GCF=GCF_000882815.3 REF=ref ACC=NC_060948.1
#   make get_genome species=zika GCF=GCF_000882815.3 REF=ref ACC=NC_012532.1

get_genome:
	mkdir -p $(dir ${REF})
	bio fetch ${ACC} --format fasta > ${REF}/${species}_genome.fa
	@echo ">> Downloaded reference genome for $(species) ($(ACC))."
	#download genome and annotation from NCBI Datasets
	datasets download genome accession ${GCF} --include genome,gtf
	unzip -n ncbi_dataset.zip -d ${REF}
	mv ${REF}/ncbi_dataset/data/$(GCF)/* ${REF}/
	rm -r ${REF}/ncbi_dataset
	rm ncbi_dataset.zip
	@echo ">> Downloaded and extracted genome and annotation for $(species) ($(GCF))."

# Download and process FASTQ files (handles both single and paired-end)
# Usage: make get_fastq srr=SRR3191545 READS=reads fastqcreports=reads/fastqc_reports

get_fastq:
	@echo ">> Downloading reads for ${srr}..."
	mkdir -p ${READS}
	mkdir -p ${fastqcreports}

	@echo ">> Fetching first 100,000 reads for test..."
	fastq-dump --split-files --gzip -X 100000 ${srr} -O ${READS}/

	@echo ">> Running FastQC..."
	fastqc -o ${fastqcreports} ${READS}/*.fastq.gz

	@echo ">> Generating read statistics..."
	seqkit stats ${READS}/*.fastq.gz > ${READS}/${srr}_all_stats.txt
	cat ${READS}/${srr}_all_stats.txt

	@echo ">> Done! FASTQ and QC files saved under ${READS}/"


# Build genome index (supports human or Zika)
# Usage:
#   make genome_index species=human REF=ref/human_genome.fa
#   make genome_index species=zika  REF=ref/zika_genome.fa

genome_index:
	@echo ">> Building genome index for $(species) using reference $(REF)..."
ifeq ($(species),human)
	bwa index ${REF}
	samtools faidx ${REF}
	@echo ">> BWA index built successfully for human genome."
else ifeq ($(species),zika)
	bwa index ${REF}
	samtools faidx ${REF}
	@echo ">> BWA index built successfully for Zika virus genome."
else
	@echo "Error: Unsupported species '$(species)'. Please use species=human or species=zika."
	exit 1
endif


# Align reads to genome using BWA (supports single- or paired-end)
# Usage:
#   make alignreads species=human REF=ref/human_genome.fa READS=reads srr=SRR3191545 bam=bam
#   make alignreads species=zika  REF=ref/zika_genome.fa  READS=reads srr=SRR3191545 bam=bam

alignreads:
	@echo ">> Aligning $(srr) using BWA index $(REF)"
	mkdir -p $(bam)
	FASTQS=($(wildcard $(READS)/$(srr)*.fastq*))
	@if [ $${#FASTQS[@]} -eq 2 ]; then \
		echo ">> Detected paired-end reads: $${FASTQS[@]}"; \
		bwa mem -t 8 "$(REF)" "$${FASTQS[0]}" "$${FASTQS[1]}" | samtools sort -@ 4 -o "$(bam)/$(srr)_$(species).bam"; \
	elif [ $${#FASTQS[@]} -eq 1 ]; then \
		echo ">> Detected single-end read: $${FASTQS[0]}"; \
		bwa mem -t 8 "$(REF)" "$${FASTQS[0]}" | samtools sort -@ 4 -o "$(bam)/$(srr)_$(species).bam"; \
	else \
		echo "!! ERROR: No FASTQ files found for $(srr) in $(READS)"; \
		exit 1; \
	fi
	samtools index "$(bam)/$(srr)_$(species).bam"
	samtools flagstat "$(bam)/$(srr)_$(species).bam" > "$(bam)/$(srr)_$(species)_align_stats.txt"
	@echo ">> Alignment complete: $(bam)/$(srr)_$(species).bam"



# Convert BAM to BedGraph + BigWig for IGV/Genome Browser
# Works for RNA-seq (spliced alignments) or genomic reads
# Usage:
#   make bigwig species=human srr=SRR3191545 bam=bam REF=ref/human_genome.fa
#   make bigwig species=zika  srr=SRR3191545 bam=bam REF=ref/zika_genome.fa

bigwig:
	@echo ">> Generating BedGraph and BigWig for ${srr} (${species})..."
	mkdir -p ${bam}

	@echo ">> Step 1: Creating raw BedGraph..."
	bedtools genomecov -ibam ${bam}/${srr}_${species}.bam -bg -split > ${bam}/${srr}_${species}_raw.bedgraph

	@echo ">> Step 2: Sorting BedGraph..."
	sort -k1,1 -k2,2n ${bam}/${srr}_${species}_raw.bedgraph > ${bam}/${srr}_${species}_sorted.bedgraph

	@echo ">> Step 3: Preparing chromosome size file..."
	cut -f1,2 ${REF}.fai > ${bam}/${species}_chrom.sizes

	@echo ">> Step 4: Filtering valid chromosomes..."
	awk 'NR==FNR {a[$$1]; next} $$1 in a' ${bam}/${species}_chrom.sizes ${bam}/${srr}_${species}_sorted.bedgraph > ${bam}/${srr}_${species}_filtered.bedgraph

	@echo ">> Step 5: Creating BigWig file..."
	bedGraphToBigWig ${bam}/${srr}_${species}_filtered.bedgraph ${bam}/${species}_chrom.sizes ${bam}/${srr}_${species}.bw

	@echo ">> Step 6: Compressing and indexing BedGraph..."
	bgzip -f ${bam}/${srr}_${species}_filtered.bedgraph
	tabix -p bed ${bam}/${srr}_${species}_filtered.bedgraph.gz

	@echo ">> Done! BigWig created: ${bam}/${srr}_${species}.bw"


# Call variants using bcftools
# Usage:
#   make call_variants species=human REF=ref/human_genome.fa bam=bam srr=SRR3191545 vcf=vcf
#   make call_variants species=zika  REF=ref/zika_genome.fa bam=bam srr=SRR3191545 vcf=vcf

call_variants:
	@echo ">> Calling variants for ${srr} (${species})..."
	mkdir -p ${vcf}

	@echo ">> Step 1: Generating genotype likelihoods and calling variants..."
	bcftools mpileup -d 100 --annotate INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP -O u -f ${REF} ${bam}/${srr}_${species}.bam | \
	bcftools call --ploidy 2 --annotate FORMAT/GQ -mv -O u | \
	bcftools norm -f ${REF} -d all -O u | \
	bcftools sort -O z -o ${vcf}/${srr}_${species}.vcf.gz

	@echo ">> Step 2: Indexing VCF..."
	bcftools index -t ${vcf}/${srr}_${species}.vcf.gz

	@echo ">> Done! Variant calls saved to: ${vcf}/${srr}_${species}.vcf.gz"


.PHONY: help design get_genome get_fastq genome_index alignreads bigwig call_variants all