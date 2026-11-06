# Bulk RNA-seq differential expression with DESeq2 (airway dataset)
# Run from the repo folder. Output goes to the working directory.
library(DESeq2)
library(airway)
library(ggplot2)

data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds$dex <- relevel(dds$dex, ref = "untrt")

keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]

dds <- DESeq(dds)
res <- results(dds, contrast = c("dex", "trt", "untrt"))
res$symbol <- rowData(dds)$symbol
summary(res)

vsd <- vst(dds, blind = FALSE)
p_pca <- plotPCA(vsd, intgroup = c("dex", "cell"))
ggsave("pca.png", p_pca, width = 7, height = 5, dpi = 150)

df <- as.data.frame(res)
df <- df[!is.na(df$padj), ]
df$sig <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1,
                 "significant", "not significant")
p_volcano <- ggplot(df, aes(log2FoldChange, -log10(padj), colour = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c("grey60", "red")) +
  labs(title = "Dex-treated vs untreated",
       x = "log2 fold change", y = "-log10 adjusted p-value") +
  theme_minimal()
ggsave("volcano.png", p_volcano, width = 7, height = 5, dpi = 150)

out <- as.data.frame(res)
out$gene_id <- rownames(out)
write.csv(out[order(out$padj), ], "airway_deseq2_results.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
