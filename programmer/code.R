SRR1177981 = read.table("SRR1177981_Aligned", sep = "")
SRR1177981 = SRR1177981[-1,]
SRR1177981 = SRR1177981[, 7]

SRR1177982 = read.table("SRR1177982_Aligned", sep = "")
SRR1177982 = SRR1177982[-1,]
SRR1177982 = SRR1177982[, 7]

SRR1177983 = read.table("SRR1177983_Aligned", sep = "")
SRR1177983 = SRR1177983[-1,]
SRR1177983 = SRR1177983[, 7]

SRR1178008 = read.table("SRR1178008_Aligned", sep = "")
SRR1178008 = SRR1178008[-1,]
SRR1178008 = SRR1178008[, 7]

SRR1178009 = read.table("SRR1178009_Aligned", sep = "")
SRR1178009 = SRR1178009[-1,]
SRR1178009 = SRR1178009[, 7]

SRR1178010 = read.table("SRR1178010_Aligned", sep = "")
SRR1178010 = SRR1178010[-1,]
SRR1178010 = SRR1178010[, 7]

SRR1178014 = read.table("SRR1178014_Aligned", sep = "")
SRR1178014 = SRR1178014[-1,]
SRR1178014 = SRR1178014[, 7]

SRR1178021 = read.table("SRR1178021_Aligned", sep = "")
SRR1178021 = SRR1178021[-1,]
SRR1178021 = SRR1178021[, 7]

SRR1178047 = read.table("SRR1178047_Aligned", sep = "")
SRR1178047 = SRR1178047[-1,]
SRR1178047 = SRR1178047[, 7]

geneid = read.table("SRR1177981_Aligned", sep = "")
geneid = geneid[, 1]
geneid = geneid[-1]

count_matrix = data.frame(geneid, SRR1177981, SRR1177982, SRR1177983, SRR1178008, SRR1178009, SRR1178010, SRR1178014, SRR1178021, SRR1178047)
write.csv(count_matrix, "count_matrix.csv")

#SRR1177981 = as.integer(paste(SRR1177981))
#SRR1177981 = as.integer(levels(SRR1177981))[SRR1177981]

graphmatrix = count_matrix
graphmatrix[graphmatrix==0] = 1
boxplot(as.integer(paste(graphmatrix$SRR1177981)), 
        as.integer(paste(graphmatrix$SRR1177982)), 
        as.integer(paste(graphmatrix$SRR1177983)), 
        as.integer(paste(graphmatrix$SRR1178008)),
        as.integer(paste(graphmatrix$SRR1178009)),
        as.integer(paste(graphmatrix$SRR1178010)),
        as.integer(paste(graphmatrix$SRR1178014)),
        as.integer(paste(graphmatrix$SRR1178021)),
        as.integer(paste(graphmatrix$SRR1178047)),
        main = "feature counts for SRR1177981", names = names(count_matrix)[-1], log = "y", xlab = "Samples", ylab = "Feature Counts")

count_matrixnew = read.csv("count_matrix.csv", header = TRUE, row.names=1)
row.names(count_matrixnew) = count_matrixnew[,1]
count_matrixnew = count_matrixnew[, -1]
control = read.csv("control_counts.csv", header = TRUE)
control = control[,-1]
control = control[,c("SRR1178050", "SRR1178061", "SRR1178063", "SRR1178004", "SRR1178006", "SRR1178013")]
AhRsamples = cbind(count_matrixnew[,c("SRR1178008", "SRR1178009", "SRR1178010")], control[,c("SRR1178050", "SRR1178061", "SRR1178063")])
CARPXRsamples = cbind(count_matrixnew[,c("SRR1178014", "SRR1178021", "SRR1178047")], control[,c("SRR1178050", "SRR1178061", "SRR1178063")])
DNADMGsamples = cbind(count_matrixnew[,c("SRR1177981", "SRR1177982", "SRR1177983")], control[,c("SRR1178004", "SRR1178006", "SRR1178013")])

library(DESeq2)


# load counts for AhR
AhRcnts = AhRsamples

# filter out rows that have any zeros for funzies
AhRcnts = subset(AhRcnts,rowSums(AhRcnts==0)==0)

# sample information
info = read.csv('toxgroup_3_rna_info.csv', header = TRUE, row.names=1)
AhRcntsinfo = info[c("SRR1178008", "SRR1178009", "SRR1178010", "SRR1178050", "SRR1178061", "SRR1178063"), ]
# create the DESeq object
dds1 = DESeqDataSetFromMatrix(
  countData = AhRcnts,
  colData = AhRcntsinfo,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds1$mode_of_action = relevel(dds1$mode_of_action, ref='Control')

# run DESeq
dds1 = DESeq(dds1)
res1 = results(dds1, contrast=c('mode_of_action','AhR','Control'))
res1 = lfcShrink(dds1, coef=2)
# write out DE results
write.csv(res1,'AhR_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds1,normalized=TRUE),'AhR_deseq_norm_counts.csv')

#For CARPXR
CARPXRcnts = CARPXRsamples
CARPXRcnts = subset(CARPXRcnts,rowSums(CARPXRcnts==0)==0)
CARPXRinfo = info[c("SRR1178014", "SRR1178021", "SRR1178047", "SRR1178050", "SRR1178061", "SRR1178063"), ]
dds2 = DESeqDataSetFromMatrix(
  countData = CARPXRcnts,
  colData = CARPXRinfo,
  design= ~ mode_of_action
)
dds2 = DESeq(dds2)
res2 = results(dds2, contrast=c('mode_of_action','CAR/PXR','Control'))
res2 = lfcShrink(dds2, coef=2)
write.csv(res2,'CARPXR_deseq_results.csv')
write.csv(counts(dds2,normalized=TRUE),'CARPXR_deseq_norm_counts.csv')


#For DNADMG
DNADMGcnts = DNADMGsamples
DNADMGcnts = subset(DNADMGcnts,rowSums(DNADMGcnts==0)==0)
DNADMGinfo = info[c("SRR1177981", "SRR1177982", "SRR1177983", "SRR1178004", "SRR1178006", "SRR1178013"), ]
dds3 = DESeqDataSetFromMatrix(
  countData = DNADMGcnts,
  colData = DNADMGinfo,
  design= ~ mode_of_action
)
dds3 = DESeq(dds3)
res3 = results(dds3, contrast=c('mode_of_action','DNA_Damage','Control'))
res3 = lfcShrink(dds3, coef=2)
write.csv(res3,'DNADMG_deseq_results.csv')
write.csv(counts(dds3,normalized=TRUE),'DNADMG_deseq_norm_counts.csv')

res1_sorted = res1[order(res1$padj),]
res1_sorted = res1_sorted[complete.cases(res1_sorted), ]
res2_sorted = res2[order(res2$padj),]
res2_sorted = res2_sorted[complete.cases(res2_sorted), ]
res3_sorted = res3[order(res3$padj),]
res3_sorted = res3_sorted[complete.cases(res3_sorted), ]

AhR_sig_genecounts = dim(res1_sorted[res1_sorted$padj<0.05,])[1]
AhR_top10 = row.names(res1_sorted)[1:10]
CARPXR_sig_genecounts = dim(res2_sorted[res2_sorted$padj<0.05,])[1]
CARPXR_top10 = row.names(res2_sorted)[1:10]
DNADMG_sig_genecounts = dim(res3_sorted[res3_sorted$padj<0.05,])[1]
DNADMG_top10 = row.names(res3_sorted)[1:10]

Result = data.frame(col1=c(AhR_sig_genecounts, AhR_top10), col2=c(CARPXR_sig_genecounts, CARPXR_top10), col3=c(DNADMG_sig_genecounts, DNADMG_top10))
names(Result) = c("AhR", "CARPXR", "DNADMG")
rownames(Result) = c("Significant genes count", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th", "9th", "10th")
write.csv(Result, "Toxgroup significant gene count and top 10 most significant genes.csv")

hist(res1_sorted[c(1:AhR_sig_genecounts),"log2FoldChange"], main = "Histogram of Log2 fold change in expression for AhR samples", xlab = "Fold change", ylab = "", breaks=20)
hist(res2_sorted[c(1:CARPXR_sig_genecounts),"log2FoldChange"], main = "Histogram of Log2 fold change in expression for CARPXR samples", xlab = "Fold change", ylab = "", breaks=20)
hist(res3_sorted[c(1:DNADMG_sig_genecounts),"log2FoldChange"], main = "Histogram of Log2 fold change in expression for DNADMG samples", xlab = "Fold change", ylab = "", breaks=20)

library(ggplot2)
library(ggrepel)
#volcano plot for AhR
res1_sorted_foldch = res1[order(res1$log2FoldChange),]
res1_sorted_foldch = data.frame(res1_sorted_foldch)
AhRlabels = rbind(head(res1_sorted_foldch), tail(res1_sorted_foldch))
AhRlabels$gene = rownames(AhRlabels)
ggplot(res1_sorted_foldch, aes(x=log2FoldChange, y = -log(padj)))+
  geom_point()+
  geom_hline(yintercept = -log(0.05))+
  geom_point(data=AhRlabels, col = "red")+
  geom_text_repel(data=AhRlabels, aes(label=gene))

#volcano plot for CARPXR
res2_sorted_foldch = res2[order(res2$log2FoldChange),]
res2_sorted_foldch = data.frame(res2_sorted_foldch)
CARPXRlabels = rbind(head(res2_sorted_foldch), tail(res2_sorted_foldch))
CARPXRlabels$gene = rownames(CARPXRlabels)
ggplot(res2_sorted_foldch, aes(x=log2FoldChange, y = -log(padj)))+
  geom_point()+
  geom_hline(yintercept = -log(0.05))+
  geom_point(data=CARPXRlabels, col = "red")+
  geom_text_repel(data=CARPXRlabels, aes(label=gene))

#volcano plot for DNADMG
res3_sorted_foldch = res3[order(res3$log2FoldChange),]
res3_sorted_foldch = data.frame(res3_sorted_foldch)
DNADMGlabels = rbind(head(res3_sorted_foldch), tail(res3_sorted_foldch))
DNADMGlabels$gene = rownames(DNADMGlabels)
ggplot(res3_sorted_foldch, aes(x=log2FoldChange, y = -log(padj)))+
  geom_point()+
  geom_hline(yintercept = -log(0.05))+
  geom_point(data=DNADMGlabels, col = "red")+
  geom_text_repel(data=DNADMGlabels, aes(label=gene))
