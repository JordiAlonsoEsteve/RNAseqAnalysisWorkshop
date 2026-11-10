library(DESeq2)
library(stringr)
setwd("/home/neggor/PhD/RNA_seq_analysis")
source("code/Statistical_Analysis/utils.R")
dataset = "PDAMP"
## WARNING (take into account what has to be remove and upper/lower cases) TODO
file <- str_glue("data/{dataset}/counts_genes_x_samples.tsv")
my_colData <- build_coldata(file, ref_treatment ="mock", drop_treatments = c("untreated", "Untreated", "Room1", "Room2", "flg22", "Pep1"))
counts  <- load_counts(file, samples = rownames(my_colData))
print(my_colData)
# NOTE! We get the reference of the factors the baseline treatment, mock, temperature 21 and time 0.5.
# The time is not very important, but is easier to understand if we set it at the earliest!
print(my_colData$temperature)
print(my_colData$treatment)
print(my_colData$time)

print(counts)

# ---- design ----
design_formula <- ~ treatment * time * temperature

# full rank check: fails if any treatment x time x temperature cell is empty
mm <- model.matrix(design_formula, my_colData)
stopifnot(qr(mm)$rank == ncol(mm))
message(ncol(mm), " coefficients for ", nrow(mm), " samples")

### Fit
dds <- DESeqDataSetFromMatrix(countData = counts, colData = my_colData, design = design_formula)
smallestGroupSize <- min(table(interaction(my_colData$treatment, my_colData$time, my_colData$temperature)))

# remove genes with little counts
# keep genes that reach at lplot_gene_trajectoryeast 10 counts in as many samples as the smallest group has, so a gene expressed in only one condition isn’t thrown out.
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]
nrow(dds)
# This is the object that you might want to store (Many tables of pvalues and lfc)
dds <- DESeq(dds)   # size factors, dispersions, Wald tests
# save dds fitted object
saveRDS(dds, file = str_glue("data/{dataset}/dds_object.rds"))
# Now this for gene-by-gene analysis
run_gene_viewer(dds)
plot_gene_diff_from_mock_heatmap(dds, gene = gene)
# Group analysis


