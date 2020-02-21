tpmAnalysis <- function(tpm,sampleGroup=NULL,geneGroup=NULL,deGenes=NULL,hmRun=FALSE,hmUserOrder=NULL,dendroRow=TRUE,scaleHM=NULL,corrRun=FALSE,corrUserOrder=NULL,scaleMin=NULL,pcaRun=FALSE,pcs=2,pcgseRun=FALSE){
"tpmAnalysis.R - The purpose of this function is to take a single matrix of l2tpm values and easily analyze by:
    1. creating a heatmap and clustering
    2. creating a correlation matrix and clustering
    3. performing PCA
Within these, the user can easily specify which samples and genes to consider in each analysis."

#User inputs - note that only required argument is tpm, but practically also one of hmRun, corrRun, or pcaRun or else code won't do anything
#tpm is data matrix with column order GQ, SYM, trimmed sample names
#subset_string is string by which to subset data; if none provided, will use all samples
#subset_genes is list of genes by which to subset data; if none provided, will use all genes
#deGenes is data frame of DESeq counts output from GoSeq for heatmap of DE genes only within a given condition  

#hmRun is whether to plot heatmap for selected data; default is FALSE
#hmUserOrder is user-supplied list of columns to be used for order of heatmap; Must be possible through flipping nodes (which means practically
  #you must run with automatic clustering first to determine desired order for rearrangement); Default is automatic clustering
#dendroRow defaults to clustering rows of heatmap but can be set to FALSE  
  
#corrRun is whether to plot correlation matrix for selected data; default is FALSE
#corrUserOrder is user-supplied list of columns to be used for order of heatmap; Must be possible through flipping nodes; Default is 
  #automatic clustering
#scaleMin is minimum for scale on corr mat heatmap. Default is NULL, so scaling is done automatically.
  
#pcaRun is whether to run PCA for selected data; default is FALSE
#pcs is number of principal components to save loadings and PCGSE results for; default is 2
#pcgseRun is whether to run PCGSE on PCA results; default is FALSE and will not run unless pcaRUN set to TRUE
  

# Part 1. Library load & data processing -------------------------------------------------
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(ggfortify)
  library(gplots)
  require(gdata)
  require(ggdendro)
  library(viridis)
  library(dendextend)

#Subset data by genes and/or column names specified by user; if no genes specified, set rownames as GQ.ID and remove GQ.ID/SYM columns
tpm<-na.omit(tpm)
#If a geneGroup is specified, this section will find the gene symbols of that group and subset data
if (!isTRUE(is.null(geneGroup))){  
  geneOnto<-read.delim("C:\\Users\\User1\\Dropbox (MIT)\\00 Tools\\GOslim_lookup.txt")
  geneSub<-subset(geneOnto,grepl(geneGroup,geneOnto$termName)); subset_genes<-as.character(geneSub[,1])
  idx.genes=match(subset_genes,tpm[,2]); tpmNum<-tpm[idx.genes,c(-1,-2)]; rownames(tpmNum)<-make.names(tpm[idx.genes,2],unique=TRUE)
#Else the full gene list will be processed and used
} else { tpmNum<-tpm[,c(-1,-2)]; rownames(tpmNum)<-tpm[,1]}  
#Subset the columns (samples) based on specified string
if (!isTRUE(is.null(sampleGroup))){
  tpmNum<-subset(tpmNum,select=grepl(sampleGroup,colnames(tpmNum)))
  } 
#Limit to only DEGenes based on provided GOSeq output
if (!is.null(deGenes)){
  deGenes<-subset(deGenes,select=grepl(sampleGroup,colnames(deGenes)))
  deGenesonly<-deGenes[rowSums(deGenes) != 0, ]
  idx.DE=match(rownames(deGenesonly),rownames(tpmNum))
  tpmNum<-tpmNum[idx.DE,]
}

# Part 2a. Heatmap of l2TPMs -------------------------------------------------------
#Create heatmap if hmRun=TRUE and cluster accordingly based on other inputs


if(isTRUE(hmRun)){
  print("hello heatmap")

  
  if (!isTRUE(is.null(scaleHM))){scale.hm=scaleHM; breaks=seq(-1.5,1.5,by=(3/50))} else {scale.hm="none"; breaks=NULL}
  
  if(isTRUE(is.null(hmUserOrder))){   #If hmUserOrder specified, will order columns by this and create separate dendrogram
    if(isTRUE(dendroRow)){dend="both"} else {dend="column"}
    pdf("Heatmap1.pdf",width=18,height=12)
    heatmap.2(as.matrix(tpmNum),trace="none",col=viridis(50),Colv=TRUE,Rowv=dendroRow,dendrogram=dend,scale=scale.hm,breaks=breaks,margins=c(12,8))
    dev.off()
  } else {
    tpmNum.hm<-tpmNum[,hmUserOrder]
    if(isTRUE(dendroRow)){dend="row"} else {dend="none"}
    pdf("Heatmap1.pdf",width=9,height=6)
    heatmap.2(as.matrix(tpmNum.hm),trace="none",col=viridis(50),Colv=FALSE,Rowv=dendroRow,dendrogram=dend,scale=scale.hm,breaks=breaks,margins=c(12,8))
    dev.off()
    
    distance.col=dist(t(tpmNum.hm),method="euclidean")
    cluster.col=hclust(distance.col,method="ward.D2")
    lookup<-unlist(cluster.col["order"])
    idx.dend=match(colnames(tpmNum.hm),colnames(tpmNum.hm)[lookup])
    pdf("DendroHeat1.pdf",width=6,height=3)
    plot(rotate(as.dendrogram(cluster.col),idx.dend))
    dev.off()
  }
}

# Part 2b. Correlation matrix and dendrogram ------------------------------------------------------
#Create correlation matrix if corrRun=TRUE and cluster accordingly based on other inputs

if(isTRUE(corrRun)){
  print("hello corrmat")
  
  breaks=NULL     #Allows user to specify scale for corr plot, but defaults to automatic if no scaleMin entered
  if (!isTRUE(is.null(scaleMin))){
    breaks=seq(scaleMin,1,by=(1-scaleMin)/20)
  }
  if(isTRUE(is.null(corrUserOrder))){
      pdf("Corrplot1.pdf",width=6,height=6)
      heatmap.2(cor(tpmNum),trace="none",col=viridis(20),breaks=breaks)
      dev.off()
    } else {
      tpmNum.cp<-tpmNum[,corrUserOrder]
      pdf("Corrplot1.pdf",width=6,height=6)
      heatmap.2(cor(tpmNum.cp),trace="none",col=viridis(20),breaks=breaks,Colv=FALSE,Rowv=FALSE,dendrogram="none")
      dev.off()
      
      distance.col=dist(t(tpmNum.cp),method="euclidean")
      cluster.col=hclust(distance.col,method="ward.D2")
      lookup<-unlist(cluster.col["order"])
      idx.dend=match(colnames(tpmNum.cp),colnames(tpmNum.cp)[lookup])
      pdf("DenddroCorr1.pdf",width=6,height=3)
      plot(rotate(as.dendrogram(cluster.col),idx.dend))
      dev.off()
    }
}

# Part 2c. Principal component analysis -----------------------------------------------------------
#Conduct PCA and PCGSE depending on user specification
if (isTRUE(pcaRun)){
  library(MASS)  #Required for PCGSE
  library(PCGSE) #PCGSE package
  library(WriteXLS) #To save results to xlsx
  print("hello pca")
  
  #Transpose data and load gene.sets, but eliminate any genes not present in dataset (likely due to NAs)
  tpmNum_trans<-setNames(data.frame(t(tpmNum)),rownames(tpmNum))
  pcaRes<-prcomp(tpmNum_trans)
  pdf("PCAplot1.pdf",width=9,height=6,useDingbats=FALSE)
  print(autoplot(pcaRes,label=TRUE,label.size=4))   
  dev.off()
  if (isTRUE(pcgseRun)){
    gene.sets_raw<-as.matrix(read.delim("C:\\Users\\User1\\Dropbox (MIT)\\00 Tools\\GOslim_matrix.txt"))
    idx.pca<-na.omit(match(colnames(tpmNum_trans),colnames(gene.sets_raw)))
    gene.sets<-gene.sets_raw[,idx.pca]
    rownames(gene.sets)<-gene.sets_raw[,1]
    
    ## Execute PCGSE using Fisher-transformed correlation coefficients as the gene-level statistics,
    ## the standardized mean difference as the gene set statistic and a correlation-adjusted t-test
    ## for the determination of statistical significance.
    
    pcgse.results = pcgse(data=tpmNum_trans,
                          prcomp.output=pcaRes,
                          pc.indexes=1:pcs,
                          gene.sets=gene.sets,
                          gene.statistic="loading",
                          transformation="none",
                          gene.set.statistic="mean.diff",
                          gene.set.test="cor.adj.parametric")
    #for (i in 1:pcs) {
    #  pcgse.results$p.values[,i] = p.adjust(pcgse.results$p.values[,i], method="bonferroni")
    #}
    pcaLoad<-as.data.frame(pcaRes$rotation[,1:pcs])
    pcaLoad$sym<-tpm[match(rownames(pcaLoad),tpm[,1]),2]
    pcaLoad$GQid<-rownames(pcaLoad)
    pcaLoad<-pcaLoad[,c(pcs+2,pcs+1,1:pcs)]
    
    p.adj<-data.frame(pcgse.results$p.values) 
    p.adj$GOterms<-rownames(p.adj)
    p.adj<-p.adj[,c(pcs+1,1:pcs)]
    statistics<-data.frame(pcgse.results$statistics)
    statistics$GOterms<-rownames(statistics)
    statistics<-statistics[,c(pcs+1,1:pcs)]
    
    WriteXLS(c("p.adj","statistics","pcaLoad"), ExcelFileName = "PCGSE1.xls", SheetNames = c("adj.p","stats","loadings"))
    }
}
}