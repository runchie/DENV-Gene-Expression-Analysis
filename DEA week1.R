library(DESeq2)
library(tidyverse)

library(readxl)
D0_1 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D0_1.xlsx", sheet = "input")
D0_2 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D0_2.xlsx", sheet = "input")
D1_1 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D1_1.xlsx", sheet = "input")
D1_2 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D1_2.xlsx", sheet = "input")
D2_1 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D2_1.xlsx", sheet = "input")
D2_2 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D2_2.xlsx", sheet = "input")
D3_1 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D3_1.xlsx", sheet = "input")
D3_2 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D3_2.xlsx", sheet = "input")
D4_1 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D4_1.xlsx", sheet = "input")
D4_2 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/D4_2.xlsx", sheet = "input")

View(D3_2)

mock1 <- as.matrix(D0_1) # turn the xlsx into matrix
mock_1 <- mock1[,-3] # omit the length column
View(mock_1)

mock2 <- as.matrix(D0_2)
mock_2 <- mock2[,-3] # omit the length column

denv1.1 <- as.matrix(D1_1)
denv1_1 <- denv1.1[,-3] # omit the length column

denv1.2 <- as.matrix(D1_2)
denv1_2 <- denv1.2[,-3] # omit the length column

denv2.1 <- as.matrix(D2_1)
denv2_1 <- denv2.1[,-3] # omit the length column

denv2.2 <- as.matrix(D2_2)
denv2_2 <- denv2.2[,-3] # omit the length column

denv3.1 <- as.matrix(D3_1)
denv3_1 <- denv3.1[,-3] # omit the length column

denv3.2 <- as.matrix(D3_2)
denv3_2 <- denv3.2[,-3] # omit the length column

denv4.1 <- as.matrix(D4_1)
denv4_1 <- denv4.1[,-3] # omit the length column

denv4.2 <- as.matrix(D4_2)
denv4_2 <- denv4.2[,-3] # omit the length column

View(denv4_1)

colnames(mock_1)[colnames(mock_1) == "expected_count"] <- "Mock1"
View(mock_1)

colnames(mock_2)[colnames(mock_2) == "expected_count"] <- "Mock2"
View(mock_2)

colnames(denv1_1)[colnames(denv1_1) == "expected_count"] <- "DENV1_1"
View(denv1_1)

colnames(denv1_2)[colnames(denv1_2) == "expected_count"] <- "DENV1_2"
View(denv1_2)

colnames(denv2_1)[colnames(denv2_1) == "expected_count"] <- "DENV2_1"
View(denv2_1)

#colnames(denv2_2)[colnames(denv2_2) == "expected_count_denv2_2"] <- "DENV2_2"
#denv2_2

#colnames(mock_1)[colnames(mock_1) == "expected_count"] <- "mock1"
#View(mock_1)

#colnames(mock_2)[colnames(mock_2) == "expected_count"] <- "mock2"
#View(mock_2)

#colnames(denv2_1)[colnames(denv2_1) == "expected_count"] <- "d2.1"
#View(denv2_1)

#colnames(denv2_2)[colnames(denv2_2) == "expected_count"] <- "d2.2"
#View(denv2_2)

colnames(denv2_2)[colnames(denv2_2) == "expected_count"] <- "DENV2_2"
View(denv2_2)

colnames(denv3_1)[colnames(denv3_1) == "expected_count"] <- "DENV3_1"
View(denv3_1)

colnames(denv3_2)[colnames(denv3_2) == "expected_count"] <- "DENV3_2"
View(denv3_2)

colnames(denv4_1)[colnames(denv4_1) == "expected_count"] <- "DENV4_1"
View(denv4_1)

colnames(denv4_2)[colnames(denv4_2) == "expected_count"] <- "DENV4_2"
View(denv4_2)

merged_denv1_1 <- merge(mock_1, denv1_1, by = "SymbolID")
#colnames(coldv1)[colnames(coldv1) == "expected_count.x.x"] <- "Mock 1"
View(merged_denv1_1)

#############

merge_m1d2 <- merge(mock_1, denv2_1, by = "SymbolID")
View(merge_m1d2)


merge_m2d2 <- merge(mock_2, denv2_2, by = "SymbolID")
View(merge_m2d2)

merge_GA <- merge(merge_m1d2, merge_m2d2, by = "SymbolID")



#############

merged_denv1_2 <- merge(mock_2, denv1_2, by = "SymbolID")
#colnames(coldv2)[colnames(coldv2) == "expected_count.x.x"] <- "Mock 1"
View(merged_denv1_2)

merged_denv2_1 <- merge(mock_1, denv2_1, by = "SymbolID")
#colnames(merged_denv2_1)[colnames(merged_denv2_2) == "expected_count.x.x"] <- "Mock 2"
View(merged_denv2_1)

merged_denv2_2 <- merge(mock_2, denv2_2, by = "SymbolID")
#colnames(merged_denv2_2)[colnames(merged_denv2_2) == "expected_count.x.x"] <- "Mock 2"
View(merged_denv2_2)

merged_denv3_1 <- merge(mock_1, denv3_1, by = "SymbolID")
#colnames(merged_denv3_1)[colnames(merged_denv3_1) == "expected_count.x.x"] <- "Mock 1"
View(merged_denv3_1)

merged_denv3_2 <- merge(mock_2, denv3_2, by = "SymbolID")
#colnames(merged_denv3_2)[colnames(merged_denv3_2) == "expected_count.x.x"] <- "Mock 2"
View(merged_denv3_2)

merged_denv4_1 <- merge(mock_1, denv4_1, by = "SymbolID")
#colnames(merged_denv4_1)[colnames(merged_denv4_1) == "expected_count.x.x"] <- "Mock 1"
View(merged_denv4_1)

merged_denv4_2 <- merge(mock_2, denv4_2, by = "SymbolID")
#colnames(merged_denv4_2)[colnames(merged_denv4_2) == "expected_count.x.x"] <- "Mock 2"
View(merged_denv4_2)

#merged_denv1_1 <- merge(mock_1, denv1_1, by = "SymbolID", all = TRUE)
#colnames(merged_denv2_1)[colnames(merged_denv2_1) == "expected_count.x.x"] <- "Mock 1"
#View(merged_denv2_1)

merged_denv1 <- merge(merged_denv1_1, merged_denv1_2, by = "SymbolID")
View(merged_denv1) # merged mock and denv1 = countData

merged_denv2 <- merge(merged_denv2_1, merged_denv2_2, by = "SymbolID")
View(merged_denv2) # merged mock and denv2 = countData

merged_denv3 <- merge(merged_denv3_1, merged_denv3_2, by = "SymbolID")
View(merged_denv3) # merged mock and denv3 = countData

merged_denv4<- merge(merged_denv4_1, merged_denv4_2, by = "SymbolID")
View(merged_denv4) # merged mock and denv3 = countData

mgdv1 <- merged_denv1[complete.cases(merged_denv1), ] # countdata with omitted NA
View(mgdv1)

mgdv2 <- merged_denv2[complete.cases(merged_denv2), ] # countdata with omitted NA
View(mgdv2)

mgdv3 <- merged_denv3[complete.cases(merged_denv3), ] # countdata with omitted NA
View(mgdv3)

mgdv4 <- merged_denv4[complete.cases(merged_denv4), ] # countdata with omitted NA
View(mgdv4)

merge_GA <- merge_GA[complete.cases(merge_GA), ]

#remove duplicated Akr1b10 column

write_xlsx(mgdv1, "/Users/runchie/Desktop/Ph.D./DEA_runchana/mgdv1.xlsx") #export file from R to selected path

write_xlsx(mgdv2, "/Users/runchie/Desktop/Ph.D./DEA_runchana/mgdv2.xlsx") #export file from R to selected path

write_xlsx(mgdv3, "/Users/runchie/Desktop/Ph.D./DEA_runchana/mgdv3.xlsx") #export file from R to selected path

write_xlsx(mgdv4, "/Users/runchie/Desktop/Ph.D./DEA_runchana/mgdv4.xlsx") #export file from R to selected path

mgdv1 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/mgdv1.xlsx")
mgdv2 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/mgdv2.xlsx")
mgdv3 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/mgdv3.xlsx")
mgdv4 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/mgdv4.xlsx")

View(mgdv1)
View(mgdv2)
View(mgdv3)
View(mgdv4)



library(tidyverse)
coldv1 <- mgdv1 %>% remove_rownames %>% column_to_rownames(var="SymbolID")
View(coldv1) # colData for DENV1 set

coldv2 <- mgdv2 %>% remove_rownames %>% column_to_rownames(var="SymbolID")
View(coldv2) # colData for DENV2 set

coldv3 <- mgdv3 %>% remove_rownames %>% column_to_rownames(var="SymbolID")
View(coldv3) # colData for DENV3 set

coldv4 <- mgdv4 %>% remove_rownames %>% column_to_rownames(var="SymbolID")
View(coldv4) # colData for DENV4 set

merge_GA <- merge_GA %>% remove_rownames %>% column_to_rownames(var="SymbolID")
View(merge_GA)

#colnames(mgdv2)[colnames(mgdv2) == "mock2"] <- "Mock2"
#View(mgdv2)

#colnames(mgdv3)[colnames(mgdv3) == "mock2"] <- "Mock2"
#View(mgdv3)

#coldv3 <- mgdv3 %>% remove_rownames %>% column_to_rownames(var="SymbolID")
#coldv3 # colData for DENV3 set

#coldv4 <- mgdv4 %>% remove_rownames %>% column_to_rownames(var="SymbolID")
#coldv4 # colData for DENV4 set

# expected_count.x.x = mock 1
# expected_count.y.x = DENV1_1
# expected_count.x.y = mock 2
# expected_count.y.y = DENV1_2


#merged_denv1_1 # merged mock1 and denv1_1 = countData

#merged_denv1_2 <- merge(mock_2, denv1_2, by = "SymbolID", all = TRUE)
#merged_denv1_2 # merged mock2 and denv1_2 = countData

#merged_denv1 <- merge(merged_denv1_1, merged_denv1_2, by = "SymbolID", all = TRUE)
#merged_denv1 # merged mock and denv1 = countData


mgdv1 <- merged_denv1[complete.cases(merged_denv1), ] # countdata with omitted NA

coldv2 <- mgdv2[complete.cases(mgdv2), ] # countdata with omitted NA





#row.names(mgd_v1) <- 1:nrow(mg_dv1)
#mg_dv1


#MGDV1 <- colnames(mgdv1, do.NULL = TRUE, prefix ="SymbolID")
#MGDV1 <- rownames(mgdv1, var = "SymbolID")


Meta_DENV2 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/Meta_DENV2.xlsx", sheet = "Sheet1")

View(Metadata_1)
View(mgdv1)
mode(coldv1) <- "integer"

mode(merge_GA) <- "integer"

MGA <- as.integer(unlist(merge_GA))

#library(WriteXLS)
#write_xlsx(mg_dv1, "/Users/runchie/Desktop/Ph.D./DEA_runchana/mg_dv1.xlsx") #export file from R to selected path

View(Metadata_2)

library(WriteXLS)
library(writexl)
Write_xlsx(resSig, "/Users/runchie/Desktop/Ph.D./DEA_runchana/resSig_dv1.xlsx")

write.csv( as.data.frame(resSig), file = "resultsig_denv1.csv")


resSig <- subset(res, padj < 0.01)
resSig # input for ggplot2
summary(resSig)

library(ggrepel)
vol_plot_denv1 <- resSig
ggplot(aes(x = log2FoldChange,
           y = padj)) +
  geom_point()

df1=as.data.frame(resSig)


ggplot(data=resSig, aes(x=log2FoldChange,y=padj)+
         geom_point()+
         theme_minimal()+ 
         geom_text_repel()+
         scale_color_manual(values=c('blue','black','red'))+
         theme(text=element_text(size=20)))

Metadata_1 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/Metadata.xlsx")
View(Metadata_1)

Metadata_2 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/Metadata_DENV2.xlsx")
View(Metadata_2)

Metadata_3 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/Metadata_DENV3.xlsx")
View(Metadata_3)

Metadata_4 <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/Metadata_DENV4.xlsx")
View(Metadata_4)

------
####DENV2####





DENV2_count <- coldv2 %>% remove_rownames %>% column_to_rownames(var="SymbolID")
View(DENV2_count) # colData for DENV2 set

View(merge_GA_NoRowName)

MGAdf <- as.data.frame(merge_GA)

MD2 <- as.data.frame(Metadata_2)

dds_2 <- DESeqDataSetFromMatrix(countData = round(coldvD2), colData = MD2, design= ~ batch + condition)

write.csv( as.data.frame(DENV2_count), file = "DENV2_count.csv")
countD2GA <- read_excel("/Users/runchie/Desktop/Ph.D./DEA_runchana/DENV2_count.xlsx", sheet = "DENV2_count")
coldvD2 <- countD2GA %>% remove_rownames %>% column_to_rownames(var="...1")
View(coldvD2)

dds_2$condition <- relevel(dds_2$condition, ref = "untreated")

dds_2 <- DESeq(dds_2)
resultsNames(dds_2)

resD2 <- results(dds_2, name="condition_treated_vs_untreated")
view(resD2)

resSigD2 <- subset(resD2, padj < 0.01)
resSigD2 # input for ggplot2
summary(resSigD2)


##############

view(Metadata_1)

dds_1 <- counts(dds, normalized = T)
  
  
------
  
dds <- DESeqDataSetFromMatrix(countData = round(coldv1),
                              colData = Metadata_1,
                              design= ~ batch + condition)


dds$condition <- relevel(dds$condition, ref = "untreated")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name="condition_treated_vs_untreated")
view(res)

resSig <- subset(res, padj < 0.01)
resSig # input for ggplot2
summary(resSig)

view(Metadata_1)

dds_1 <- counts(dds, normalized = T)


write.csv( as.data.frame(resSig), file = "resultsig_denv1.1.csv")

df1=as.data.frame(resSig)

vsd1 <- varianceStabilizingTransformation(dds)
#plotPCA(vsd1, intgroup=c("condition", "sample"))

sampleDists1 <- dist(t(assay(vsd1)))
library("RColorBrewer")
sampleDistMatrix1 <- as.matrix(sampleDists1)
rownames(sampleDistMatrix1) <- paste(vsd1$condition, vsd1$type, sep="-")
colnames(sampleDistMatrix1) <- NULL
colors1 <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix1,
         clustering_distance_rows=sampleDists1,
         clustering_distance_cols=sampleDists1,
         col=colors1)

df1=as.data.frame(resSig)
df1$diffexpressed <- "NO"
head(df1)
df1$delabel <- NA

df1$diffexpressed[df1$log2FoldChange>1 & df1$pvalue<0.01] <- "Upregulated"
df1$diffexpressed[df1$log2FoldChange<1 & df1$pvalue<0.01] <- "Downregulated"
df1$diffexpressed[df1$log2FoldChange>-1 & df1$log2FoldChange<1] <- "Non-significant"

plot1 <- ggplot(data=df1,aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + geom_point(shape = 21) + theme_minimal() + geom_text() + geom_text_repel() + scale_color_manual(values =c('blue','grey','red')) + theme(text=element_text(size=20)) + geom_vline(xintercept=c(-1, 1),linetype = "dashed", col='grey') + geom_hline(yintercept=-log(0.001),linetype = "dashed", col='grey') + ggtitle("DENV-1")
------
  
dds2 <- DESeqDataSetFromMatrix(countData = round(coldv2),
                              colData = Metadata_2,
                              design= ~ batch + condition)

View(coldv2)

dds2$condition <- relevel(dds2$condition, ref = "untreated")

dds2 <- DESeq(dds2)
resultsNames(dds2)

res2 <- results(dds2, name="condition_treated_vs_untreated")
view(res2)

resSig2 <- subset(res2, padj < 0.01)
resSig2 # input for ggplot2
summary(resSig2)


write.csv( as.data.frame(resSig2), file = "resultsig_denv2.1.csv")

ggplot(data=df2,aes(x=log2FoldChange,y=-log10(pvalue))) + 
  geom_point()

df2=as.data.frame(resSig2)
df2$diffexpressed <- "NO"
head(df2)
df2$delabel <- NA

df2$diffexpressed[df2$log2FoldChange>1 & df2$pvalue<0.01] <- "Upregulated"
df2$diffexpressed[df2$log2FoldChange<1 & df2$pvalue<0.01] <- "Downregulated"
df2$diffexpressed[df2$log2FoldChange>-1 & df2$log2FoldChange<1] <- "Non-significant"

plot2 <- ggplot(data=df2,aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + geom_point(shape = 21) + theme_minimal() + geom_text() + geom_text_repel() + scale_color_manual(values =c('blue','grey','red')) + theme(text=element_text(size=20)) + geom_vline(xintercept=c(-1, 1),linetype = "dashed", col='grey') + geom_hline(yintercept=-log(0.001),linetype = "dashed", col='grey') + ggtitle("DENV-2")



-----------------
  
dds3 <- DESeqDataSetFromMatrix(countData = round(coldv3),
                                 colData = Metadata_3,
                                 design= ~ batch + condition)


dds3$condition <- relevel(dds3$condition, ref = "untreated")

dds3 <- DESeq(dds3)
resultsNames(dds3)

res3 <- results(dds3, name="condition_treated_vs_untreated")
view(res3)

resSig3 <- subset(res3, padj < 0.01)
resSig3 # input for ggplot2
summary(resSig3)
View(resSig3)

write.csv( as.data.frame(resSig3), file = "resultsig_denv3.1.csv")

df3=as.data.frame(resSig3)
df3$diffexpressed <- "NO"
head(df3)
df3$delabel <- NA

df3$diffexpressed[df3$log2FoldChange>1 & df3$pvalue<0.01] <- "Upregulated"
df3$diffexpressed[df3$log2FoldChange<1 & df3$pvalue<0.01] <- "Downregulated"
df3$diffexpressed[df3$log2FoldChange>-1 & df3$log2FoldChange<1] <- "Non-significant"


plot3 <- ggplot(data=df3,aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + geom_point(shape = 21) + theme_minimal() + geom_text() + geom_text_repel() + scale_color_manual(values =c('blue','grey','red')) + theme(text=element_text(size=20)) + geom_vline(xintercept=c(-1, 1),linetype = "dashed", col='grey') + geom_hline(yintercept=-log(0.001),linetype = "dashed", col='grey') + ggtitle("DENV-3")

---------
  
dds4 <- DESeqDataSetFromMatrix(countData = round(coldv4),
                                 colData = Metadata_4,
                                 design= ~ batch + condition)


dds4$condition <- relevel(dds4$condition, ref = "untreated")

dds4 <- DESeq(dds4)
resultsNames(dds4)

res4 <- results(dds4, name="condition_treated_vs_untreated")
view(res4)

resSig4 <- subset(res4, padj < 0.01)
resSig4 # input for ggplot2
summary(resSig4)


write.csv( as.data.frame(resSig4), file = "resultsig_denv4.1.csv")


df4=as.data.frame(resSig4)
df4$diffexpressed <- "NO"
head(df4)
df4$delabel <- NA

df4$diffexpressed[df4$log2FoldChange>1 & df4$pvalue<0.01] <- "Upregulated"
df4$diffexpressed[df4$log2FoldChange<1 & df4$pvalue<0.01] <- "Downregulated"
df4$diffexpressed[df4$log2FoldChange>-1 & df4$log2FoldChange<1] <- "Non-significant"


plot4 <- ggplot(data=df4,aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + geom_point(shape = 21) + theme_minimal() + geom_text() + geom_text_repel() + scale_color_manual(values =c('blue','grey','red')) + theme(text=element_text(size=20)) + geom_vline(xintercept=c(-1, 1),linetype = "dashed", col='grey') + geom_hline(yintercept=-log(0.001),linetype = "dashed", col='grey') + ggtitle("DENV-4")

