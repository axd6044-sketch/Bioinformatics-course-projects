#This is a makefile for RNA-seq analysis
#Process both single and paired-end reads using pattern-based rules

SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

#variables
ref ?= ref/genome
species ?= Zika
genome ?= GCF_000882815.3 / GCA_000882815.1
srr ?= SRR3194430
reads ?= reads
fastqcreports ?= $(reads)/fastqc_reports
genome_fa ?= $(ref)/$(species)genome.fa
bam ?= bam
PRJNA ?= PRJNA313294


# ------ No changes in code beyond this point -----

all: design get_genome get_fastq genome_index alignreads bigwig 
	echo ">> All steps completed successfully."

# Display help information
help:
	@echo "Available targets:"
	@echo "  design          - Generate design.csv from BioProject"
	@echo "  get_genome        - Download reference genome"
	@echo "  get_fastq        - Download and process FASTQ files (uses -X 100000 to limit data size)"
	@echo "Usage:"
	@echo "  make design PRJNA=PRJNA313294"
	@echo "  make get_genome ref=ref genome=chr22 species=human"
	@echo "  make get_fastq srr=SRR3194430 reads=reads fastqcreports=reads/fastqc_reports"
	@echo "  make genome_index genome_fa=ref/hmanchr22.fa"
	@echo "  make alignreads genome_fa=ref/hmanchr22 reads=reads srr=SRR3194430 bam=bam"
	@echo "  make bigwig srr=SRR3194430 bam=bam genome_fa=ref/hmanchr22"


# Generate design.csv from BioProject
# Usage: make design PRJNA=PRJNA313294
make design:
	bio search -H --csv ${PRJNA} > design.csv
	cat design.csv

# Download reference genome for human analysis (chromosome 22 subset)
# Usage: make get_genome ref=ref genome=chr22 species=human
get_genome:
	@echo ">> Downloading human reference genome (chr22, hg38 build)..."
	mkdir -p $(ref)
	curl -L -o $(ref)/$(genome).fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$(genome).fa.gz
	gunzip -f $(ref)/$(genome).fa.gz
	mv $(ref)/$(genome).fa $(ref)/hmanchr22.fa
	samtools faidx $(ref)/hmanchr22.fa

# Download and process FASTQ files (handles both single and paired-end)
# Usage: make get_fastq srr=SRR3194430 reads=reads fastqcreports=reads/fastqc_reports
get_fastq:
	mkdir -p $(reads)
	mkdir -p $(fastqcreports)
	fastq-dump --split-3 --gzip -X 100000 $(srr) -O $(reads)/
	gunzip $(reads)/*.fastq.gz || true
	# fastqc and seqkit will operate on all FASTQ files in the reads directory
	fastqc -o $(fastqcreports) $(reads)/*.fastq
	seqkit stats $(reads)/*.fastq > $(reads)/${srr}_all_stats.txt
	ls -l $(reads)

# Index genome for BWA and samtools
# Usage: make genome_index genome_fa=ref/hmanchr22.fa
#bwa index -a is -p $(basename $(genome_fa)) $(genome_fa)
#samtools faidx $(genome_fa)
genome_index:
	# build HISAT2 index (spliced aligner)
	hisat2-build $(genome_fa) $(basename $(genome_fa))

# Align reads to genome (works for single- or paired-end FASTQ)
# Usage: make alignreads genome_fa=ref/hmanchr22 reads=reads srr=SRR3194430 bam=bam
alignreads:
	@echo ">> Aligning $(srr) using HISAT2 index $(genome_fa)"
	mkdir -p $(bam)
	FASTQS=($(wildcard $(reads)/$(srr)*.fastq*))
	@if [ $${#FASTQS[@]} -eq 2 ]; then \
		echo ">> Detected paired-end reads: $${FASTQS[@]}"; \
		hisat2 -p 4 -x "$(genome_fa)" -1 "$${FASTQS[0]}" -2 "$${FASTQS[1]}" | samtools sort -o "$(bam)/$(srr).bam"; \
	elif [ $${#FASTQS[@]} -eq 1 ]; then \
		echo ">> Detected single-end reads: $${FASTQS[0]}"; \
		hisat2 -p 4 -x "$(genome_fa)" -U "$${FASTQS[0]}" | samtools sort -o "$(bam)/$(srr).bam"; \
	else \
		echo "!! ERROR: No FASTQ files found for $(srr) in $(reads)"; \
		exit 1; \
	fi
	samtools index "$(bam)/$(srr).bam"
	samtools flagstat "$(bam)/$(srr).bam" > "$(bam)/$(srr)_align_stats.txt"
	@echo ">> Alignment complete: $(bam)/$(srr).bam"


# Convert BAM to BedGraph + BigWig for IGV/Genome Browser
# Works for RNA-seq (spliced alignments) or genomic reads
# Usage: make bigwig srr=SRR3194430 bam=bam genome_fa=ref/hmanchr22.fa
bigwig:
	mkdir -p $(bam)
	bedtools genomecov -ibam "$(bam)/$(srr).bam" -bg -split | awk '$$3>$$2' > "$(bam)/$(srr)_raw.bedgraph"
	sort -k1,1 -k2,2n "$(bam)/$(srr)_raw.bedgraph" -o "$(bam)/$(srr)_sorted.bedgraph"
	cut -f1,2 "$(genome_fa).fai" > "$(bam)/genome.chrom.sizes"
	awk 'NR==FNR {a[$$1]; next} $$1 in a' "$(bam)/genome.chrom.sizes" "$(bam)/$(srr)_sorted.bedgraph" > "$(bam)/$(srr)_filtered.bedgraph"
	bedGraphToBigWig "$(bam)/$(srr)_filtered.bedgraph" "$(bam)/genome.chrom.sizes" "$(bam)/$(srr).bw"
	bgzip -f "$(bam)/$(srr)_filtered.bedgraph"
	tabix -p bed "$(bam)/$(srr)_filtered.bedgraph.gz"



.PHONY: help design get_genome get_fastq genome_index alignreads bigwig all