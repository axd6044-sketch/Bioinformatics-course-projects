# CONFIGURATION
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# File paths
REF        := refs/chr22.genome.fa
GTF        := refs/chr22.gtf
READS_DIR  := reads
BAM_DIR    := bam
DESIGN     := design.csv

COUNTS_TXT := counts.txt
COUNTS_FMT := countsformatted.csv
COUNTS_CSV := counts.csv

# Samples from design file: skip header, take first column
SAMPLES := $(shell tail -n +2 $(DESIGN) 2>/dev/null | cut -d, -f1)
BAMS    := $(addprefix $(BAM_DIR)/,$(addsuffix .bam,$(SAMPLES)))

.PHONY: all getref stats design index alignall counts strandcounts formatcounts tx2gene counts_csv clean

# Run the whole pipeline
all: getref stats design index alignall counts counts_csv

# no changes beyond this point

# 1. DOWNLOAD DATA & REFS
# usage: make getref

getref:
	wget -nc http://data.biostarhandbook.com/data/uhr-hbr.tar.gz
	tar xzvf uhr-hbr.tar.gz

# 2. STATS (REF + FASTQ)
# usage: make stats

stats:
	seqkit stats refs/*.fa
	seqkit stats $(READS_DIR)/*.fq


# 3. DESIGN FILE
# usage: make design

design: $(DESIGN)

$(DESIGN): $(READS_DIR)/*_R1.fq
	ls $(READS_DIR)/*_R1.fq \
	  | sed 's/$(READS_DIR)\/\(.*\)_R1\.fq/\1/' \
	  | sort \
	  | awk 'BEGIN{OFS=","; print "sample,group"} {split($$1,a,"_"); print $$1,a[1]}' \
	  > $@


# 4. HISAT2 INDEX
# usage: make index

index:
	make -f src/run/hisat2.mk index REF=$(REF)


# 5. ALIGNMENTS
# usage: make alignall

alignall: $(DESIGN)
	mkdir -p $(BAM_DIR)
	cat $(DESIGN) | \
	  parallel --colsep , --header : --eta --lb --jobs 2 \
	    make -f src/run/hisat2.mk \
	      REF=$(REF) \
	      R1=$(READS_DIR)/{sample}_R1.fq \
	      BAM=$(BAM_DIR)/{sample}.bam \
	      run


# 6. FEATURECOUNTS (COUNT MATRIX)
# usage: make counts

counts: $(COUNTS_TXT)

$(COUNTS_TXT): $(BAMS) $(GTF)
	featureCounts -a $(GTF) -o $@ $(BAMS)

# Strand-specific 
# usage: make strandcounts
strandcounts: $(BAMS) $(GTF)
	featureCounts -a $(GTF) -o counts_strand_s2.txt -s 2 bam/HBR_1.bam


# 7. REFORMAT COUNTS (R SCRIPTS)
# usage: make formatcounts

formatcounts: $(COUNTS_TXT)
	Rscript src/r/format_featurecounts.r -c $(COUNTS_TXT) -o $(COUNTS_FMT)

# 8. TX2GENE MAPPING & COUNTS CSV
# usage: make tx2gene
tx2gene:
	Rscript src/r/create_tx2gene.r -s > names.txt
	Rscript src/r/create_tx2gene.r -d hsapiens_gene_ensembl

# Final counts CSV
# usage: make counts_csv
counts_csv: $(COUNTS_TXT) tx2gene
	Rscript src/r/format_featurecounts.r -c $(COUNTS_TXT) -t tx2gene.csv -o $(COUNTS_CSV)
