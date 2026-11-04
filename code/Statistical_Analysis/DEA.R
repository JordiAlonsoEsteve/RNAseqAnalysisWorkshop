setwd("/home/neggor/PhD/RNA_seq_analysis")
source("code/Statistical_Analysis/utils.R")
library(DESeq2)
file    <- "data/PDAMP/counts_genes_x_samples.tsv"
colData <- build_coldata(file)
counts  <- load_counts(file, samples = rownames(colData))
print(colData)
# NOTE! We get the reference of the factors the baseline treatment, mock, temperature 21 and time 0.5.
# The time is not very important, but is easier to understand if we set it at the earliest!
print(colData$temperature)
print(colData$treatment)
print(colData$time)

print(counts)

### design
# ---- design ----
design_formula <- ~ treatment * time * temperature

# full rank check: fails if any treatment x time x temperature cell is empty
mm <- model.matrix(design_formula, colData)
stopifnot(qr(mm)$rank == ncol(mm))
message(ncol(mm), " coefficients for ", nrow(mm), " samples")

### Fit
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = design_formula)
smallestGroupSize <- min(table(interaction(colData$treatment, colData$time, colData$temperature)))

# remove genes with little counts
# keep genes that reach at least 10 counts in as many samples as the smallest group has, so a gene expressed in only one condition isn’t thrown out.
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep, ]
nrow(dds)

dds <- DESeq(dds)   # size factors, dispersions, Wald tests