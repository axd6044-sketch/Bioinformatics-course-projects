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
srr ?= SRR3194430
READS=reads
fastqcreports ?= $(reads)/fastqc_reports
genome_fa ?= $(ref)/$(species)genome.fa
bam ?= bam
PRJNA ?= PRJNA313294
HG38= GCF_000001405.40


# ------ No changes in code beyond this point -----

all: design get_genome get_fastq genome_index alignreads bigwig 
	echo ">> All steps completed successfully."

# Display help information
help:
	@echo "Available targets:"
	@echo "  design      - Generate design.csv from BioProject"
	@echo "  get_genome  - Download reference genome"
	@echo "  get_fastq   - Download and process FASTQ files (uses -X 100000 to limit data size)"

	@echo "Usage:"
	@echo "  make design PRJNA=PRJNA313294"
	@echo "  make get_genome REF=ref genome=chr22 species=human"
	@echo "  make get_fastq srr=SRR3194430 reads=reads fastqcreports=reads/fastqc_reports"
	@echo "  make genome_index genome_fa=ref/hmanchr22.fa"
	@echo "  make alignreads genome_fa=ref/hmanchr22 reads=reads srr=SRR3194430 bam=bam"
	@echo "  make bigwig srr=SRR3194430 bam=bam genome_fa=ref/hmanchr22"


# Generate design.csv from BioProject
# Usage: make design PRJNA=PRJNA313294
make design:
	bio search -H --csv ${PRJNA} > design.csv
	cat design.csv

# Download reference genome
# Usage examples:
# Download reference genome for human or Zika
# Usage examples:
#   make get_genome species=human genome=chr22 REF=ref
#   make get_genome species=human genome=hg38_ucsc REF=ref
#   make get_genome species=zika genome=GCF_000882815.3 REF=ref

get_genome:
ifeq ($(species),human)
	@echo ">> Downloading human reference genome ($(genome), hg38 build)..."
	mkdir -p $(REF)
# Option 1: Download specific chromosome (e.g., chr22)
ifneq (,$(findstring chr,$(genome)))
	curl -L -o $(REF)/$(genome).fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$(genome).fa.gz
	gunzip -f $(REF)/$(genome).fa.gz
	mv $(REF)/$(genome).fa $(REF)/human_$(genome).fa
	samtools faidx $(REF)/human_$(genome).fa
# Option 2: Download entire UCSC hg38 assembly (compact file)
else ifeq ($(genome),hg38_ucsc)
	curl -L -o $(REF)/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
    curl -L -o ref/human_hg38.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz
	gunzip -f ref/human_hg38.gtf.gz
	tabix -p gtf ref/human_hg38.gtf.gz
	mv $(REF)/hg38.fa $(REF)/human_hg38.fa
	samtools faidx $(REF)/human_hg38.fa
else
	@echo "!! Unknown genome option for human. Use 'chrN' (e.g., chr22) or 'hg38_ucsc'."
	exit 1
endif

else ifeq ($(species),zika)
	@echo ">> Downloading Zika virus reference genome ($(genome))..."
	mkdir -p $(REF)
	datasets download genome accession $(genome) --include genome,gff3 --filename $(REF)/${genome}_genome.zip
	unzip -n -o $(REF)/${genome}_genome.zip -d $(REF)
	rm $(REF)/${genome}_genome.zip
	mv $(REF)/ncbi_dataset/data/$(genome)/*.fna $(REF)/ZikaVirus_genome.fa
	mv $(REF)/ncbi_dataset/data/$(genome)/*.gff $(REF)/ZikaVirus_genome.gff
	rm -rf $(REF)/ncbi_dataset
	samtools faidx $(REF)/ZikaVirus_genome.fa

else
	@echo "Error: Unsupported species '$(species)'. Please use species=human or species=zika."
	exit 1
endif


# Download and process FASTQ files (handles both single and paired-end)
# Usage: make get_fastq srr=SRR3194430 READS=reads fastqcreports=reads/fastqc_reports
get_fastq:
	mkdir -p $(READS)
	mkdir -p $(fastqcreports)
	fastq-dump --split-3 --gzip -X 100000 $(srr) -O $(READS)/
	gunzip $(READS)/*.fastq.gz || true
	# fastqc and seqkit will operate on all FASTQ files in the reads directory
	fastqc -o $(fastqcreports) $(READS)/*.fastq
	seqkit stats $(READS)/*.fastq > $(reads)/${srr}_all_stats.txt
	ls -l $(READS)

# Index genome for BWA and samtools
# Usage: make genome_index REF=ref/hmanchr22.fa
#bwa index -a is -p $(basename $(REF)) $(REF)
#samtools faidx $(REF)
# Build genome index depending on species
# Usage examples:
#   make genome_index species=human REF=ref/human_chr22.fa
#   make genome_index species=human REF=ref/human_hg38.fa  ## for full human genome
#   make genome_index species=zika REF=ref/ZikaVirus_genome.fa

genome_index:
ifeq ($(species),human)
	@echo ">> Building HISAT2 index for human reference genome..."
	hisat2-build $(REF) $(basename $(REF))
	@echo ">> HISAT2 index for human built successfully."
else ifeq ($(species),zika)
	@echo ">> Building BWA index for Zika virus reference genome..."
	bwa index $(REF)
	@echo ">> BWA index for Zika virus built successfully."
else
	@echo "Error: Unsupported species '$(species)'. Please use species=human or species=zika."
	exit 1
endif

# Align reads to genome (supports human or Zika, single- or paired-end)
# Usage examples:
#   make alignreads species=human REF=ref/human_hg38 READS=reads srr=SRR3194430 bam=bam
#   make alignreads species=human REF=ref/human_chr22 READS=reads srr=SRR3194430 bam=bam
#   make alignreads species=zika  REF=ref/ZikaVirus_genome READS=reads srr=SRR3194430 bam=bam

alignreads:
	@echo ">> Aligning $(srr) for $(species) using reference $(REF)"
	mkdir -p $(bam)
	FASTQS=($(wildcard $(READS)/$(srr)*.fastq*))

ifeq ($(species),human)
	@if [ $${#FASTQS[@]} -eq 2 ]; then \
		echo ">> Detected paired-end reads: $${FASTQS[@]}"; \
		hisat2 -p 8 --dta -x "$(REF)" -1 "$${FASTQS[0]}" -2 "$${FASTQS[1]}" \
		| samtools sort -@ 4 -o "$(bam)/$(srr)_human.bam"; \
	elif [ $${#FASTQS[@]} -eq 1 ]; then \
		echo ">> Detected single-end reads: $${FASTQS[0]}"; \
		hisat2 -p 8 --dta -x "$(REF)" -U "$${FASTQS[0]}" \
		| samtools sort -@ 4 -o "$(bam)/$(srr)_human.bam"; \
	else \
		echo "!! ERROR: No FASTQ files found for $(srr) in $(READS)"; \
		exit 1; \
	fi

else ifeq ($(species),zika)
	@if [ $${#FASTQS[@]} -eq 2 ]; then \
		echo ">> Detected paired-end reads: $${FASTQS[@]}"; \
		bwa mem -t 8 "$(REF).fa" "$${FASTQS[0]}" "$${FASTQS[1]}" \
		| samtools sort -@ 4 -o "$(bam)/$(srr)_zika.bam"; \
	elif [ $${#FASTQS[@]} -eq 1 ]; then \
		echo ">> Detected single-end reads: $${FASTQS[0]}"; \
		bwa mem -t 8 "$(REF).fa" "$${FASTQS[0]}" \
		| samtools sort -@ 4 -o "$(bam)/$(srr)_zika.bam"; \
	else \
		echo "!! ERROR: No FASTQ files found for $(srr) in $(reads)"; \
		exit 1; \
	fi

else
	@echo "Error: Unsupported species '$(species)'. Please use species=human or species=zika."
	exit 1
endif

	@echo ">> Indexing and collecting alignment statistics..."
	samtools index "$(bam)/$(srr)_$(species).bam"
	samtools flagstat "$(bam)/$(srr)_$(species).bam" > "$(bam)/$(srr)_$(species)_align_stats.txt"
	@echo ">> Alignment complete: $(bam)/$(srr)_$(species).bam"



# Convert BAM to BedGraph + BigWig for IGV/Genome Browser
# Works for RNA-seq (spliced alignments) or genomic reads
# Usage: make bigwig srr=SRR3194430 bam=bam REF=ref/hmanchr22.fa
# Usage: make bigwig srr=SRR3194430 bam=bam REF=ref/human_hg38.fa species=human
# Usage: make bigwig srr=SRR17310612 bam=bam REF=ref/ZikaVirus_genome.fa
bigwig:
	mkdir -p $(bam)
	bedtools genomecov -ibam "$(bam)/$(srr)_$(species).bam" -bg -split | awk '$$3>$$2' > "$(bam)/$(srr)_$(species)_raw.bedgraph"
	sort -k1,1 -k2,2n "$(bam)/$(srr)_$(species)_raw.bedgraph" -o "$(bam)/$(srr)_$(species)_sorted.bedgraph"
	cut -f1,2 "$(REF).fai" > "$(bam)/genome.chrom.sizes"
	awk 'NR==FNR {a[$$1]; next} $$1 in a' "$(bam)/genome.chrom.sizes" "$(bam)/$(srr)_$(species)_sorted.bedgraph" > "$(bam)/$(srr)_$(species)_filtered.bedgraph"
	bedGraphToBigWig "$(bam)/$(srr)_$(species)_filtered.bedgraph" "$(bam)/genome.chrom.sizes" "$(bam)/$(srr)_$(species).bw"
	bgzip -f "$(bam)/$(srr)_$(species)_filtered.bedgraph"
	tabix -p bed "$(bam)/$(srr)_$(species)_filtered.bedgraph.gz"

#calling variants using bcftools
#Usage:
#   make call_variants species=human REF=ref/human_hg38.fa bam=bam srr=SRR3194430 vcf=vcf 
#   make call_variants species=zika  REF=ref/ZikaVirus_genome.fa bam=bam srr=SRR17310612 vcf=vcf
#   make call_variants species=human REF=ref/human_chr22.fa bam=bam srr=SRR3194430 vcf=vcf

call_variants:
	@echo ">> Calling variants for $(srr) ($(species))..."
	mkdir -p $(vcf)
	F1="-d 100 --annotate INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP"
	F2="--ploidy 2 --annotate FORMAT/GQ"
	BAM="$(bam)/$(srr)_$(species).bam"
	VCF="$(vcf)/$(srr)_$(species).vcf.gz"
	bcftools mpileup $${F1} -O u -f "$(REF)" "$${BAM}" | \
	bcftools call $${F2} -mv -O u | \
	bcftools norm -f "$(REF)" -d all -O u | \
	bcftools sort -O z -o "$${VCF}"
	@echo ">> Step 5: Indexing VCF..."
	bcftools index -t "$${VCF}"
	@echo ">> Done! Variant calls saved to: $${VCF}"

.PHONY: help design get_genome get_fastq genome_index alignreads bigwig call_variants all