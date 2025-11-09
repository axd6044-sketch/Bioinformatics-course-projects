# ============================
# Makefile for RNA-seq / mapping
# ============================

SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# ----------- user-ish variables -------------
SPECIES ?= zika           # override: make get_genome SPECIES=human
PRJNA   ?= PRJNA313294
SRR     ?= SRR3191542

REF     ?= ref
READS   ?= reads
FASTQC  ?= $(READS)/fastqc_reports
BAM     ?= bam
VCF     ?= vcf
GFF	 ?= $(REF)/$(SPECIES)_genome.gff
GCF	 ?= GCF_009914755.1
GCFILE ?=ref/ncbi_dataset/data/$(GCF)/*.gff
N       ?= 100000        # number of reads to download (for testing)

# per-species NC_ accessions (the ones efetch can give us)
ACC_human := NC_060948.1
ACC_zika  := NC_012532.1
GCF_human := GCF_009914755.1
GCF_zika  := GCF_000882815.3

# derived filenames (don’t depend on species logic at parse-time)
REF_FA  := $(REF)/$(SPECIES)_genome.fa
REF_GFF := $(REF)/$(SPECIES)_genome.gff
BAMFILE := $(BAM)/$(SRR)_$(SPECIES).bam
VCFILE  := $(VCF)/$(SRR)_$(SPECIES).vcf.gz

# SnpEff variables
NAME := $(strip $(SPECIES))
SNPEFF_IDX   := $(CURDIR)/idx/snpEff
SNPEFF_NAME  := $(SPECIES)
SNPEFF_DIR   := snpEFF_ann
SNPEFF_DB    := $(SNPEFF_IDX)/$(SNPEFF_NAME)/snpEffectPredictor.bin
SNPEFF_ANN  ?= snpEFF_ann/$(SRR)_$(SPECIES).snpeff.vcf.gz
SNPEFF_HTML ?= snpEFF_ann/$(SRR)_$(SPECIES).snpeff.html
SNPEFF_CSV  ?= snpEFF_ann/$(SRR)_$(SPECIES).snpeff.csv
SNPEFF_GFF=ref/human_genome.NC_060948.1.gff
# default
all: design get_genome get_fastq genome_index alignreads bigwig call_variants snpeff
	@echo ">> All steps completed successfully."

help:
	@echo "Targets:"
	@echo "  make design PRJNA=PRJNA313294"
	@echo "  make get_genome SPECIES=human|zika"
	@echo "  make genome_index SPECIES=human|zika"
	@echo "  make get_fastq SRR=SRR3191545"
	@echo "  make alignreads SPECIES=human|zika SRR=SRR..."
	@echo "  make bigwig SPECIES=human|zika SRR=SRR..."
	@echo "  make call_variants SPECIES=human|zika SRR=SRR..."
	@echo "Variables:"
	@echo "  REF=$(REF)  READS=$(READS)  BAM=$(BAM)  VCF=$(VCF)"


# A. DESIGN: 

design:
	bio search -H --csv ${PRJNA} > design.csv
	cat design.csv


# B. GET GENOME: needs species 
#usage: make get_genome SPECIES=human ACC_human=NC_060948.1 GCF_human=GCF_009914755.1
#usage: make get_genome SPECIES=zika  ACC_zika=NC_012532.1  GCF_zika=GCF_000882815.3

get_genome:
	@if [ "$(SPECIES)" = "human" ]; then \
		ACC="$(ACC_human)"; \
		GCF="$(GCF_human)"; \
	elif [ "$(SPECIES)" = "zika" ]; then \
		ACC="$(ACC_zika)"; \
		GCF="$(GCF_zika)"; \
	else \
		echo "Unsupported SPECIES '$(SPECIES)'. Use SPECIES=human or SPECIES=zika."; \
		exit 1; \
	fi; \
	mkdir -p "$(REF)"; \
	echo ">> Downloading FASTA for $$ACC ..."; \
	bio fetch $$ACC --format fasta > "$(REF_FA)"; \
	echo ">> Saved FASTA to $(REF_FA)"; \
	echo ">> Downloading GFF for $$GCF ..."; \
	datasets download genome accession $$GCF --include gff3 >/dev/null; \
	unzip -n ncbi_dataset.zip -d "$(REF)" >/dev/null; \
	cp "$(REF)/ncbi_dataset/data/$$GCF/"*.gff "$(REF_GFF)"; \
	rm -rf "$(REF)/ncbi_dataset" ncbi_dataset.zip; \
	echo ">> Saved GFF to $(REF_GFF)"


# C. FASTQ download + QC
#usage: make get_fastq SRR=SRR3191542 READS=reads FASTQC=reads/fastqc_reports N=100000


STAT := $(READS)/$(SRR)_all_stats.txt
get_fastq:
	@echo ">> Downloading reads for ${SRR}..."
	mkdir -p ${READS} ${FASTQC}
	fastq-dump --split-files --gzip -X ${N} ${SRR} -O ${READS}/
	@echo ">> Running FastQC..."
	fastqc -o ${FASTQC} ${READS}/*.fastq.gz
	@echo ">> Read stats..."
	seqkit stats ${READS}/*.fastq.gz > ${STAT}
	cat ${STAT}
	@echo ">> Done."


# D. GENOME INDEX
# usage: make genome_index SPECIES=human REF_FA=ref/human_genome.fa
# usage: make genome_index SPECIES=zika  REF_FA=ref/zika_genome.fa

genome_index: get_genome
	@echo ">> Building genome index for $(SPECIES) using $(REF_FA)"
	bwa index $(REF_FA)
	samtools faidx $(REF_FA)
	@echo ">> Index built."


# E. ALIGN READS
# usage: make alignreads SPECIES=human SRR=SRR3191542 REF_FA=ref/human_genome.fa BAM=bam
# usage: make alignreads SPECIES=zika  SRR=SRR3191542 REF_FA=ref/zika_genome.fa  BAM=bam

alignreads:
	@echo ">> Aligning ${SRR} to $(REF_FA) ..."
	mkdir -p $(BAM)
	FASTQS=($(wildcard $(READS)/$(SRR)*.fastq*))
	@if [ $${#FASTQS[@]} -eq 2 ]; then \
		echo ">> Paired-end: $${FASTQS[@]}"; \
		bwa mem -t 8 "$(REF_FA)" "$${FASTQS[0]}" "$${FASTQS[1]}" | samtools sort -@ 4 -o "$(BAMFILE)"; \
	elif [ $${#FASTQS[@]} -eq 1 ]; then \
		echo ">> Single-end: $${FASTQS[0]}"; \
		bwa mem -t 8 "$(REF_FA)" "$${FASTQS[0]}" | samtools sort -@ 4 -o "$(BAMFILE)"; \
	else \
		echo "!! ERROR: no FASTQ for $(SRR) in $(READS)"; \
		exit 1; \
	fi
	samtools index "$(BAMFILE)"
	samtools flagstat "$(BAMFILE)" > "$(BAM)/$(SRR)_$(SPECIES)_align_stats.txt"
	@echo ">> Alignment complete: $(BAMFILE)"


# F. BIGWIG
# usage: make bigwig SPECIES=human SRR=SRR3191542 REF_FA=ref/human_genome.fa BAM=bam
# usage: make bigwig SPECIES=zika  SRR=SRR3191542 REF_FA=ref/zika_genome.fa  BAM=bam

bigwig:
	@echo ">> Generating BigWig for $(SRR) ($(SPECIES))"
	mkdir -p $(BAM)
	# 1) bedgraph
	bedtools genomecov -ibam $(BAMFILE) -bg -split > $(BAM)/$(SRR)_$(SPECIES)_raw.bedgraph
	# 2) sort
	sort -k1,1 -k2,2n $(BAM)/$(SRR)_$(SPECIES)_raw.bedgraph > $(BAM)/$(SRR)_$(SPECIES)_sorted.bedgraph
	# 3) chrom sizes
	cut -f1,2 $(REF_FA).fai > $(BAM)/$(SPECIES)_chrom.sizes
	# 4) filter
	awk 'NR==FNR {a[$$1]; next} $$1 in a' $(BAM)/$(SPECIES)_chrom.sizes $(BAM)/$(SRR)_$(SPECIES)_sorted.bedgraph > $(BAM)/$(SRR)_$(SPECIES)_filtered.bedgraph
	# 5) bigwig
	bedGraphToBigWig $(BAM)/$(SRR)_$(SPECIES)_filtered.bedgraph $(BAM)/$(SPECIES)_chrom.sizes $(BAM)/$(SRR)_$(SPECIES).bw
	# 6) compress + index bedgraph
	bgzip -f $(BAM)/$(SRR)_$(SPECIES)_filtered.bedgraph
	tabix -p bed $(BAM)/$(SRR)_$(SPECIES)_filtered.bedgraph.gz
	@echo ">> BigWig written to $(BAM)/$(SRR)_$(SPECIES).bw"


# G. CALL VARIANTS
# usage: make call_variants SPECIES=human SRR=SRR3191542 REF_FA=ref/human_genome.fa BAM=bam VCF=vcf
# usage: make call_variants SPECIES=zika  SRR=SRR3191542 REF_FA=ref/zika_genome.fa  BAM=bam VCF=vcf
call_variants:
	@echo ">> Calling variants for $(SRR) ($(SPECIES))"
	mkdir -p $(VCF)
	bcftools mpileup -d 100 -f $(REF_FA) $(BAMFILE) \
	  | bcftools call --ploidy 2 -mv -O u \
	  | bcftools norm -f $(REF_FA) -d all -O u \
	  | bcftools sort -O z -o $(VCFILE)
	bcftools index -t $(VCFILE)
	@echo ">> VCF: $(VCFILE)"

# H. SnpEff ANNOTATE VARIANTS

# snpEFF

$(SNPEFF_DB): $(SNPEFF_GFF) $(REF_FA)
	# make genome folder
	mkdir -p "$(SNPEFF_IDX)/$(NAME)"

	# Copy the files to the snpEff folder.
	cp -f "$(REF_FA)" "$(SNPEFF_IDX)/$(NAME)/sequences.fa"
	cp -f "$(SNPEFF_GFF)" "$(SNPEFF_IDX)/$(NAME)/genes.gff"

	# Make the configuration file.
	echo "$(NAME).genome : $(NAME)" > $(SNPEFF_DIR)/snpeff.config

	# Build the database.
	snpEff build -dataDir "$(SNPEFF_IDX)" -v "$(NAME)"

# Create the annotated VCF file.
$(SNPEFF_ANN) $(SNPEFF_HTML) $(SNPEFF_CSV): $(VCFILE) $(SNPEFF_DB)
	mkdir -p "$(dir $(SNPEFF_ANN))" "$(dir $(SNPEFF_HTML))" "$(dir $(SNPEFF_CSV))"
	snpEff ann \
	    -csvStats "$(SNPEFF_CSV)" \
		-s "$(SNPEFF_HTML)" \
		-dataDir "$(SNPEFF_IDX)" \
		-v "$(NAME)" \
		"$(VCFILE)" | bcftools view -O z -o "$(SNPEFF_ANN)"
	bcftools index "$(SNPEFF_ANN)"

#usage: make snpeff_build SPECIES=human SRR=SRR3191542 
#usage: make snpeff_run SPECIES=zika  SRR=SRR3191542 
snpeff_build: $(SNPEFF_DB)
	@ls -lh "$(SNPEFF_DB)"

snpeff_run: $(SNPEFF_ANN) $(SNPEFF_HTML) $(SNPEFF_CSV)
	@ls -lh "$(SNPEFF_ANN)"
	@ls -lh "$(SNPEFF_HTML)"
	@ls -lh "$(SNPEFF_CSV)"


snpeff: snpeff_build snpeff_run
getvcf: get_fastq alignreads bigwig call_variants

.PHONY: all help design get_genome get_fastq genome_index alignreads bigwig call_variants snpeff_build snpeff_run snpeff snpeff_all getvcf
