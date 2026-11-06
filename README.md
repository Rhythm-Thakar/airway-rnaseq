# Bulk RNA-seq differential expression (DESeq2)

Differential expression analysis of the Bioconductor `airway` dataset: human airway smooth muscle cells, untreated vs dexamethasone-treated, 4 cell lines (8 samples, paired).

## Method
- Design `~ cell + dex`, so the treatment effect is estimated within each cell line; `untrt` set as the reference level
- Filtered to genes with at least 10 reads in at least 4 samples (63,677 -> 16,139 genes)
- DESeq2 1.52.0 on R 4.6.2; variance-stabilised counts for PCA
- Code: `scripts/analysis.R`

## Results
- At adjusted p < 0.1: 2,607 genes up (16%) and 2,298 down (14%) in treated samples
- 938 genes pass adjusted p < 0.05 and |log2 fold change| > 1
- Top hits include SPARCL1, CACNB2, DUSP1, SAMHD1, MAOA
- Known glucocorticoid-response genes are strongly up: FKBP5 (log2FC 4.04) and PER1 (log2FC 3.18)

Files: `results/airway_deseq2_results.csv`, `results/sessionInfo.txt`, `figures/pca.png`, `figures/volcano.png`.

## What I learned
- Setting the reference level matters, since factor levels default to alphabetical order
- A paired design needs the pairing variable in the model
- A sanity check against known biology (FKBP5, PER1) builds trust in the result
