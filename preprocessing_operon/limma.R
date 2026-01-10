if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs <- c("GEOquery","limma","Biobase","sva","org.Hs.eg.db","AnnotationDbi","pheatmap","ggplot2")
BiocManager::install(pkgs, ask=FALSE, update=FALSE)

library(GEOquery)
library(limma)
library(Biobase)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)

gse <- getGEO("GSE28551", GSEMatrix=TRUE)
gset <- gse[[1]]

exprs_raw <- exprs(gset)
qx <- as.numeric(quantile(exprs_raw, c(0.,0.25,0.5,0.75,0.99,1), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) #appears to be already log2 transformed

par(mfrow=c(1,2))
boxplot(exprs_raw, las=2) #seems okay

gset@phenoData@data[["data_processing"]] #were background corrected + quantile normalized + variability less than 0.1 were removed from dataset, so okay I guess

#"The pre-normalized intensity values were extracted per spot from the data files. 
#Probes  with signal to noise ratio <3.0, or signal means under 9.5, or variability less than 0.1 were removed from dataset. 
#Control spots were filtered out. Arrays were quantile normalized (R preprocessCore library). 
#Missing expression levels were inferred by imputed by KNN (R impute library). Rscript available from authour."

pca <- prcomp(t(exprs(gset)), scale=TRUE)
pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], 
                     trimester=pheno$`trimester:ch1`)
ggplot(pca_df, aes(PC1, PC2, color=trimester)) + geom_point(size=3) #seems like no batch effect detected

pheno <- pData(gset)
group <- factor(pheno$`trimester:ch1`, levels=c("first","third"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(gset, design)
contrast.matrix <- makeContrasts(first_vs_third = first - third, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTab <- topTable(fit2, coef="first_vs_third", number=Inf, adjust.method="BH")
sum(is.na(topTab$Gene.Symbol) | topTab$Gene.Symbol == "") #5897 deleted unannotated probes
topTab <- topTab[!is.na(topTab$Gene.Symbol) & topTab$Gene.Symbol != "", ]

exprs_mat <- exprs(gset)
first_idx <- which(group == "first")
third_idx <- which(group == "third")

probe_mean_first <- rowMeans(exprs_mat[, first_idx, drop=FALSE], na.rm=TRUE)
probe_mean_third <- rowMeans(exprs_mat[, third_idx, drop=FALSE], na.rm=TRUE)

probe_means_df <- data.frame(ID = rownames(exprs_mat),
                             First.Trimester.Average = probe_mean_first,
                             Third.Trimester.Average = probe_mean_third,
                             stringsAsFactors = FALSE)

topTab <- subset(topTab, select=c("ID", "Gene.Symbol", "GENE", "Gene.Name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
topTab <- merge(topTab, probe_means_df, by = "ID", all.x = TRUE)

topTab <- rename(topTab,
    X = ID,
    SYMBOL = Gene.Symbol,
    ENTREZID = GENE,
    GENENAME = Gene.Name
  )

write.csv(topTab, "GSE28551_limma_third_vs_first.csv", row.names=TRUE)
