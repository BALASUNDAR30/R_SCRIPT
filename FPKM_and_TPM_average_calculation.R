# Read TSV replicate files
df1 <- read.delim("SRR14717729.sorted_gene_abundance.tsv", header = TRUE)
df2 <- read.delim("SRR14717728.sorted_gene_abundance.tsv", header = TRUE)
df3 <- read.delim("SRR14717719.sorted_gene_abundance.tsv", header = TRUE)

# Calculate average FPKM and TPM
df1$FPKM <- (df1$FPKM + df2$FPKM + df3$FPKM) / 3
df1$TPM  <- (df1$TPM  + df2$TPM  + df3$TPM)  / 3


# Keep relevant columns
output_df <- df1[, c("Gene.ID", "Gene.Name", "Reference", "Strand", "Start", "End", "Coverage", "FPKM", "TPM")]

# Rename columns if needed (if Gene.ID should be "Gene ID" etc.)
colnames(output_df) <- c("Gene ID", "Gene Name", "Reference", "Strand", "Start", "End", "Coverage", "FPKM", "TPM")

# Write to TSV file
write.table(output_df, file = "0hr_fpkm_tpm.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
