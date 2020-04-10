library(limma)
setwd("/projectnb2/bf528/users/group2/project3/analyst")
# sample info dataframe with array_id and chemical columns
samples <- read.csv('group_2_mic_info.csv',as.is=TRUE)
BNAPH_subsamples <- samples[which(samples$chemical == 'BETA-NAPHTHOFLAVONE' | (samples$chemical == 'Control' & samples$vehicle == 'CMC_.5_%')),]
ECONA_subsamples <- samples[which(samples$chemical == 'ECONAZOLE' | (samples$chemical == 'Control'& samples$vehicle == "CORN_OIL_100_%")),]
THIOACE_subsamples <- samples[which(samples$chemical == 'THIOACETAMIDE' | (samples$chemical == 'Control'& samples$vehicle == "SALINE_100_%")),]

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)

# subset the full expression matrix to just those in this comparison
rma.BNAPHsubset <- rma[paste0('X',BNAPH_subsamples$array_id)]
rma.ECONAsubset <- rma[paste0('X',ECONA_subsamples$array_id)]
rma.THIOACEsubset <- rma[paste0('X',THIOACE_subsamples$array_id)]

# construct a design matrix modeling treatment vs control for use by limma
BNAPHdesign <- model.matrix(
  ~factor(
    BNAPH_subsamples$chemical,
    levels=c('Control','BETA-NAPHTHOFLAVONE')
  )
)

ECONAdesign <- model.matrix(
  ~factor(
    ECONA_subsamples$chemical,
    levels=c('Control','ECONAZOLE')
  )
)
THIOACEdesign <- model.matrix(
  ~factor(
    THIOACE_subsamples$chemical,
    levels=c('Control','THIOACETAMIDE')
  )
)
colnames(BNAPHdesign) <- c('Intercept','BETA-NAPHTHOFLAVONE')
colnames(ECONAdesign) <- c('Intercept','ECONAZOLE')
colnames(THIOACEdesign) <- c('Intercept','THIOACETAMIDE')

# run limma
BNAPHfit <- lmFit(rma.BNAPHsubset, BNAPHdesign)
BNAPHfit <- eBayes(BNAPHfit)
BNAPHt <- topTable(BNAPHfit, coef=2, n=nrow(rma.BNAPHsubset), adjust='BH')
orderBNAPHt <- BNAPHt[order(BNAPHt$adj.P.Val),]

ECONAfit <- lmFit(rma.ECONAsubset, ECONAdesign)
ECONAfit <- eBayes(ECONAfit)
ECONAt <- topTable(ECONAfit, coef=2, n=nrow(rma.ECONAsubset), adjust='BH')
orderECONAt <- ECONAt[order(ECONAt$adj.P.Val),]

THIOACEfit <- lmFit(rma.THIOACEsubset, THIOACEdesign)
THIOACEfit <- eBayes(THIOACEfit)
THIOACEt <- topTable(THIOACEfit, coef=2, n=nrow(rma.THIOACEsubset), adjust='BH')
orderTHIOACEt<- THIOACEt[order(THIOACEt$adj.P.Val),]

# write out the results to file
write.csv(orderBNAPHt,'BNAPH_limma_results.csv')
write.csv(orderECONAt,'ECONA_limma_results.csv')
write.csv(orderTHIOACEt,'THIOACE_limma_results.csv')

# get significant results
sigorderBNAPHt<-orderBNAPHt[which(orderBNAPHt$adj.P.Val<0.05),]
sigorderECONAt<-orderECONAt[which(orderECONAt$adj.P.Val<0.05),]
sigorderTHIOACEt<-orderTHIOACEt[which(orderTHIOACEt$adj.P.Val<0.05),]

library(ggplot2)

pdf("B_NAPH.pdf")

print(qplot(10^sigorderBNAPHt$logFC,
      geom="histogram",  
      main = "Histogram for FoldCount", 
      xlab = "FoldCount",  
      fill=I("blue"), 
      col=I("black"),
      xlim = c(0,8)))


dev.off()

pdf("ECONA.pdf")

print(qplot(10^sigorderECONAt$logFC,
      geom="histogram",
      main = "Histogram for FoldCount",
      xlab = "FoldCount",
      fill=I("blue"),
      col=I("black"),
      xlim = c(0,20),
      ))


dev.off()

pdf("THIOACE.pdf")

print(qplot(10^sigorderTHIOACEt$logFC,
      geom="histogram",
      main = "Histogram for FoldCount",
      xlab = "FoldCount",
      fill=I("blue"),
      col=I("black"),
      xlim = c(0,30)
))


dev.off()

pdf("B_NAPH_scatter.pdf")

print(qplot(x = 10^sigorderBNAPHt$logFC,
      y = sigorderBNAPHt$P.Value,
      main = "Scatterplot of FoldCount vs P-Value",
      xlab = "FoldCount",
      ylab = "P-Value",
      xlim = c(0,8)))


dev.off()

pdf("ECONA_scatter.pdf")

print(qplot(x = 10^sigorderECONAt$logFC,  
      y = sigorderECONAt$P.Value,
      main = "Scatterplot of FoldCount vs P-Value",
      xlab = "FoldCount",
      ylab = "P-Value",
      xlim = c(0,20)))


dev.off()


pdf("THIOACE_scatter.pdf")

print(qplot(x = 10^sigorderTHIOACEt$logFC,
      y = sigorderTHIOACEt$P.Value,
      main = "Scatterplot of FoldCount vs P-Value",
      xlab = "FoldCount",
      ylab = "P-Value",
      xlim = c(0,40)))


dev.off()

