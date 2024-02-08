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
