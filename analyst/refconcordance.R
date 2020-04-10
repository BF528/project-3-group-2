refconcordance<- function(refset,probeset,mapset){
	lenref <- nrow(refset)
	lenprobe <- nrow(probeset)
	x<-c("numref","numtot","numtotshared","numconc")
	concframe<-data.frame(matrix(ncol = 4, nrow = 0))
	colnames(concframe)<-x
	numtot <- 0
	numtotshared<-0
	numconc <-0
	for (i in 1:nrow(refset)){
		imapped = mapset[which(mapset$REFSEQ==refset$X[i][1]),]
		if(nrow(imapped)>0){
			for(j in 1:nrow(probeset)){ 
				mapped=imapped[which(!is.na(imapped$PROBEID) & imapped$PROBEID==row.names(probeset[j,])[1]),]
				if(nrow(mapped)>0){
					print(i)
					print(j)
					if((refset$log2FoldChange[i][1]>0 & probeset$logFC[j][1]>0)|(refset$log2FoldChange[i][1]<0 & probeset$logFC[j][1]<0)){
						numconc<-numconc+nrow(mapped)
					}
					numtotshared<-numtotshared + nrow(mapped)
				}
				numtot<-numtot+1
				concframe<-rbind(concframe,list(i,numtot,numtotshared,numconc))
			}
		}else{
			numtot<-numtot+lenprobe
			concframe<-rbind(concframe,list(i,numtot,numtotshared,numconc))
		}
	}
	concframe$concordance<- ((numconc*numtot)-(i*(numtot%i))/(numconc+numtot-i-(numtot%i)))
	return(concframe)
}
