
SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules


# Zika Virus Genome Pipeline
# Purpose: Automate genome download, indexing, alignment, feature analysis, and ~10x read subsampling
# Usage: make all

# ---- VARIABLES ----
ACC_GCF = GCF_000882815.3
ACC_GCA = GCA_000882815.1
TAXON = "Zika virus"
SRR_ID = SRR3194430
NREADS = 1500
OUTDIR = ncbi_dataset
READDIR = reads
GENOME_FA = Zikagenome.fa
GFF = Zikagenome.gff
INDEX_PREFIX = Zikagenome
ALIGN_PREFIX = alignments

# ---- RULES ----
all: genome features accessions reads subsample index align
	@echo "All pipeline steps completed successfully."

# ---- 1. Download genome and annotation ----
$(GENOME_FA) $(GFF):
	@echo ">> Downloading genome and annotation..."
	datasets summary genome accession $(ACC_GCA) | jq
	datasets download genome accession $(ACC_GCF) --include genome,gff3,gtf
	unzip -o ncbi_dataset.zip -d $(OUTDIR) -x README.md
	mv $(OUTDIR)/ncbi_dataset/data/$(ACC_GCF)/*.fna $(GENOME_FA)
	mv $(OUTDIR)/ncbi_dataset/data/$(ACC_GCF)/*.gff $(GFF)
	@echo " Genome and annotation downloaded and renamed."

genome: $(GENOME_FA) $(GFF)
	@echo ">> Computing genome size..."
	@GENOME_SIZE=$$(grep -v ">" $(GENOME_FA) | tr -d '\n' | wc -c); \
	echo "Genome size (bp): $$GENOME_SIZE"

# ---- 2. Count features in GFF ----
features: $(GFF)
	@echo ">> Counting genomic features..."
	@echo "Genomic details:"
	@cut -f1 $(GFF) | sort | uniq
	@echo ""
	@echo "Genome features are:"
	@grep -v '#' $(GFF) | cut -f3 | sort | uniq -c | sort -nr
	@echo ""
	@echo "The genes being coded are:"
	@awk '$$3=="gene" {print $$9}' $(GFF)
	@echo " Feature summary complete."

# ---- 3. List other accessions ----
accessions:
	@echo ">> Listing other genome accessions for $(TAXON)..."
	@datasets summary genome taxon $(TAXON) | grep -oE 'GC[FA]_[0-9]+\.[0-9]' | sort -u | head -10
	@echo " Accessions listed."

# ---- 4. Download sequencing reads ----
$(READDIR)/$(SRR_ID).fastq:
	@echo ">> Prefetching $(SRR_ID)..."
	prefetch $(SRR_ID)
	mkdir -p $(READDIR)
	fasterq-dump $(SRR_ID) -O $(READDIR)/
	@echo " Reads downloaded."

reads: $(READDIR)/$(SRR_ID).fastq
	@seqkit stats $(READDIR)/$(SRR_ID).fastq

# ---- 5. Subsample reads (~10x coverage) ----
reads/$(SRR_ID)_subsample.fastq: reads
	@echo ">> Subsampling $(NREADS) reads (~10x coverage)..."
	seqkit sample -n $(NREADS) $(READDIR)/$(SRR_ID).fastq > $(READDIR)/$(SRR_ID)_subsample.fastq
	seqkit stats $(READDIR)/$(SRR_ID)_subsample.fastq
	@echo " Subsampling complete."

subsample: reads/$(SRR_ID)_subsample.fastq

# ---- 6. Index genome ----
index: $(GENOME_FA)
	@echo ">> Indexing genome with BWA..."
	bwa index -p $(INDEX_PREFIX) $(GENOME_FA)
	@echo " Genome indexing complete."

# ---- 7. Align reads to genome ----
align: index reads/$(SRR_ID)_subsample.fastq
	@echo ">> Aligning reads to genome..."
	mkdir -p $(ALIGN_PREFIX)
	bwa mem $(INDEX_PREFIX) $(READDIR)/$(SRR_ID)_subsample.fastq > $(ALIGN_PREFIX)/$(SRR_ID).sam
	@echo ">> Converting SAM to sorted BAM..."
	samtools view -bS $(ALIGN_PREFIX)/$(SRR_ID).sam | samtools sort -o $(ALIGN_PREFIX)/$(SRR_ID).bam
	samtools index $(ALIGN_PREFIX)/$(SRR_ID).bam
	rm $(ALIGN_PREFIX)/$(SRR_ID).sam
	@echo " Alignment complete: $(ALIGN_PREFIX)/$(SRR_ID).bam"


.PHONY: all genome features accessions reads subsample index align clean
