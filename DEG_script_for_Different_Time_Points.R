# Load full count matrix
counts_raw <- read.csv("counts.csv")

# Combine Chr and Start as unique row names (or use gene ID if available)
rownames(counts_raw) <- counts_raw$Geneid

# Remove non-sample columns (assuming 1 to 5 are metadata)
counts <- counts_raw[, -(1:5)] 

head(counts)

# Optional: Clean sample names to match sample info (remove ".sorted.bam")
colnames(counts) <- gsub("\\.sorted\\.bam", "", colnames(counts))


sample_info <- read.csv("sample_info.csv", stringsAsFactors = FALSE)



colnames(counts)
sample_info$sample

rownames(sample_info) <- sample_info$sample
sample_info$sample <- NULL

# Convert time_point to factor with correct levels
sample_info$time_point <- factor(sample_info$time_point, levels = c("0hr", "6hr", "12hr"))

all(colnames(counts) == rownames(sample_info))  # should return TRUE


length(colnames(counts))  # Number of columns (samples) in counts
nrow(sample_info)         # Number of rows (samples) in sample_info

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ time_point)

dds <- DESeq(dds)

res <- results(dds, contrast = c("time_point", "0hr", "12hr"))
head(res)
write.csv(as.data.frame(res), "DESeq2_results_0hr_vs_12hr.csv")
