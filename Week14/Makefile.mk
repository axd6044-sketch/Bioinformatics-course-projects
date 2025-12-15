############################################################
# Week14 Makefile
# Differential expression visualization and enrichment
# Biostar workflows
############################################################

# --------------------
# Input files
# --------------------
COUNTS = edger.csv
DESIGN = design.csv

# --------------------
# Output files
# --------------------
PCA = pca.pdf
HEATMAP = heatmap.pdf
GPROF = gprofiler.csv
ENRICH = enrichr.csv

# --------------------
# Scripts
# --------------------
PCA_SCRIPT = src/r/plot_pca.r
HEATMAP_SCRIPT = src/r/plot_heatmap.r

# --------------------
# Default target
# --------------------
all: $(PCA) $(HEATMAP) $(GPROF) $(ENRICH)

############################################################
# PCA plot
############################################################
$(PCA): $(COUNTS) $(DESIGN) $(PCA_SCRIPT)
	@echo "Generating PCA plot"
	$(PCA_SCRIPT) -c $(COUNTS) -d $(DESIGN) -o $(PCA)

############################################################
# Heatmap
############################################################
$(HEATMAP): $(COUNTS) $(DESIGN) $(HEATMAP_SCRIPT)
	@echo "Generating heatmap"
	$(HEATMAP_SCRIPT) -c $(COUNTS) -d $(DESIGN) -o $(HEATMAP)

############################################################
# g:Profiler enrichment
############################################################
$(GPROF): $(COUNTS)
	@echo "Running g:Profiler enrichment"
	bio gprofiler -c $(COUNTS) -d hsapiens

############################################################
# Enrichr enrichment
############################################################
$(ENRICH): $(COUNTS)
	@echo "Running Enrichr enrichment"
	bio enrichr -c $(COUNTS)

############################################################
# Housekeeping
############################################################
clean:
	rm -f $(PCA) $(HEATMAP) $(GPROF) $(ENRICH)

