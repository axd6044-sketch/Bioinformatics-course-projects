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
genome ?= GCF_000882815.3
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
	@echo "  make get_genome genome=GCF_000882815.3 ref=ref/genome species=Zika"
	@echo "  make get_fastq srr=SRR3194430" "reads=reads fastqcreports=reads/fastqc_reports"
	@echo "  make genome_index ref=ref/genome species=Zika"
	@echo "  make alignreads genome_fa=ref/Zikagenome.fa reads=reads srr=SRR3191544 bam=bam"
	@echo "  make bigwig srr=SRR3191544 bam=bam genome_fa=ref/Zikagenome.fa"


# Generate design.csv from BioProject
# Usage: make design PRJNA=PRJNA313294
make design:
	bio search -H --csv ${PRJNA} > design.csv
	cat design.csv

# Download reference genome - only keeps fa and gff files
# Usage: make get_genome genome=GCF_000882815.3 ref=ref/genome species=Zika
get_genome:
	mkdir -p ${dir ${ref}}
	datasets download genome accession $(genome) --include genome,gff3 --filename ${genome}_genome.zip
	unzip -n -o ${genome}_genome.zip -d ${ref}
	mv ${ref}/ncbi_dataset/data/$(genome)/*.fna ${ref}/${species}genome.fa
	mv ${ref}/ncbi_dataset/data/$(genome)/*.gff ${ref}/${species}genome.gff
	rm -rf ${genome}_genome.zip ${ref}/ncbi_dataset ${ref}/README.md ${ref}/md5sum.txt
	echo "Genome saved to ${ref}/${species}genome.fa | GFF saved to ${ref}/${species}genome.gff"

# Download and process FASTQ files (handles both single and paired-end)
# Usage: make get_fastq srr=SRR3194430 reads=reads fastqcreports=reads/fastqc_reports
get_fastq:
	mkdir -p $(reads)
	mkdir -p $(fastqcreports)
	fastq-dump --split-3 --gzip -X 100000 $(srr) -O $(reads)/
	gunzip $(reads)/*.fastq.gz || true
	# fastqc and seqkit will operate on all FASTQ files in the reads directory
	fastqc -o $(fastqcreports) $(reads)/*.fastq
	seqkit stats $(reads)/*.fastq > $(reads)/all_stats.txt
	ls -l $(reads)

# Index genome for BWA and samtools
# Usage: make genome_index genome_fa=ref/Zikagenome.fa
genome_index:
	bwa index -a is -p $(basename $(genome_fa)) $(genome_fa)
	samtools faidx $(genome_fa)


# Align reads to genome (works for single- or paired-end FASTQ)
# rely on make's wildcard expansion to collect either one (single-end)
# or two (paired-end) FASTQ files. Ordering is lexicographic so _1 before _2.
# Usage: make alignreads genome_fa=ref/Zikagenome.fa reads=reads srr=SRR3191544 bam=bam
alignreads:
	mkdir -p $(bam)
	FASTQS=$(wildcard $(reads)/$(srr)*.fastq*)
	bwa mem -t 4 $(basename $(genome_fa)) $${FASTQS} | samtools sort -o $(bam)/$(srr).bam
	samtools index $(bam)/$(srr).bam
	samtools flagstat $(bam)/$(srr).bam > $(bam)/$(srr)_align_stats.txt


# Convert BAM to BedGraph + BigWig for IGV/Genome Browser
# Usage: make bigwig srr=SRR3191544 bam=bam genome_fa=ref/Zikagenome.fa
bigwig:
	bedtools genomecov -ibam $(bam)/$(srr).bam -bg > $(bam)/$(srr).bedgraph
	sort -k1,1 -k2,2n $(bam)/$(srr).bedgraph -o $(bam)/$(srr)_sorted.bedgraph
	awk '{print $$1, $$2, $$3, $$4}' OFS='\t' $(bam)/$(srr)_sorted.bedgraph > $(bam)/$(srr)_final.bedgraph
	cut -f1,2 $(genome_fa).fai > $(bam)/genome.chrom.sizes
	grep -Fwf <(cut -f1 $(bam)/genome.chrom.sizes) $(bam)/$(srr)_final.bedgraph > $(bam)/$(srr)_filtered.bedgraph
	bedGraphToBigWig $(bam)/$(srr)_filtered.bedgraph $(bam)/genome.chrom.sizes $(bam)/$(srr).bw
	bgzip -f $(bam)/$(srr)_filtered.bedgraph
	tabix -p bed $(bam)/$(srr)_filtered.bedgraph.gz


.PHONY: help design get_genome get_fastq genome_index alignreads bigwig all