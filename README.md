Proteodata:
Creates volcano plot, showing significantly up- and downregulated proteins/genes.
The csv input should contain:
- accession of any kind
- fold change data (will be log2 transformed)
- p value (will be -log10 transformed)

Correlation:
Creates a scatterplot showing correlation between e.g. proteomic and RNAseq data.
The csv input should contain:
- accession of any kind
- fold change data for x axis (e.g. RNAseq)
- fold change data for y axis (e.g. proteomic)
- the fold change data should bt log2 transformed

sPLS-DA-PCA:
Performs sPLS-DA and PCA analyses of quatitative proteomic data.
The csv input should contain:
- log2 transformed proteomic data normalized by Z-score for each individual sample
- optionally, the data should be imputed, but NaN values should be ok
- the data I've used were transformed, imputed, ANOVA-filtered and Z-score normalized
