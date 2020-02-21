deseqAnalysis<-function(pvals,lfcs,fname=NULL,upsetRun=FALSE,cond1=NULL,goseqRun=FALSE,fgseaRun=FALSE,pcutoff=NULL,plotGeneSet=NULL,volcanoRun=FALSE){
"deseqAnalysis.R - The purpose of this function is to analyze deseqMultifactor.R outputs (SignPADJ and LFC files) by:
  1. Perform intersection analysis and create Upset plots (higher dimension Venn diagram information)
  2. Conduct GOSeq analysis
  3. Conduct FGSEA analysis"
  
#User inputs: padj and lfc data is required, along with setting one of the three analyses to run. All else are optional
#resSignPADJ is the SignPADJ output from deseqMultifactor.R
#resLFC is the LFC output from deseqMultifactor.R

#upsetRun is whether to create upset plots and defaults to FALSE
  
#goseqRun is whether to run GOSeq and defaults to FALSE
  
#fgseaRun is whether to run FGSEA and defaults to FALSE

if(is.null(fname)){fname="DE"}

# Upset plot generation ---------------------------------------------------
if(isTRUE(upsetRun)){
  library(UpSetR)
  degenesLFC<-apply(lfcs,2, function(x) rownames(subset(lfcs,abs(x)>1)))
  degenesPVAL<-apply(pvals,2,function(x) rownames(subset(pvals,abs(x)<0.05)))
  
  if(length(degenesLFC)==length(degenesPVAL)){lenUpset=length(degenesLFC)} else {print("Warning: lfcs and pvals do not have equal number of columns!")}
  
  temp<-as.list(matrix('',lenUpset,1))
  names(temp)<-names(degenesPVAL)
  
  for (ii in c(1:lenUpset)){
    temp[[ii]]<-as.list(intersect(degenesLFC[[ii]],degenesPVAL[[ii]]))
    if (length(temp[[ii]])==0){temp[[ii]]="NULL";print(paste(names(temp)[ii]," contains no DE genes.",sep=""))}
  }
  
  if(is.null(cond1)){
    pdf(paste("UpSet.",fname,"-1.pdf",sep=""),width=14,height=7,useDingbats=FALSE)
    upset(fromList(temp), point.size = 3.5, line.size = 2,
        mainbar.y.label = "DE Gene Intersections", sets.x.label = "Total DE Genes",
        text.scale=c(4,4,4,4,4,4),show.numbers = FALSE,
        keep.order=FALSE, cutoff=min(c(length(degenesLFC),5)), matrix.color="black",
        group.by = "sets",order.by="freq")
    dev.off()
  } else {
  for (con1 in cond1){
    pdf(paste("UpSet.",fname,".",con1,".pdf",sep=""),width=14,height=7,onefile=FALSE,useDingbats=FALSE)
    upset(fromList(temp[grepl(con1,names(temp))]), point.size = 3.5, line.size = 2,
        mainbar.y.label = "DE Gene Intersections", sets.x.label = "Total DE Genes",
        text.scale=c(4,4,4,4,4,4),show.numbers = FALSE,
        keep.order=FALSE, cutoff=5, matrix.color="black",
        group.by = "sets",order.by="freq")
    dev.off()
  }
  }
}

# GOSeq gene ontology analysis --------------------------------------------
if(isTRUE(goseqRun)){
  library("rafalib")
  library("gplots")
  library("stringi")
  library("goseq")
  
  ontology<-as.data.frame(read.delim("C:\\Users\\User1\\Dropbox (MIT)\\00 Tools\\GOslim_goseq.txt"))
  gene_len_raw<-read.csv("C:\\Users\\User1\\Dropbox (MIT)\\00 Tools\\gene_lengths_salmon.csv",header=FALSE)
  
  pmethod="bonferroni"
  if(is.null(pcutoff)){pcutoff=0.05}
  
  lfcs[lfcs>=1 | lfcs<= -1]<- 1
  lfcs[lfcs<1 & lfcs> -1]<- 0
  lfcs[is.na(lfcs)]<-0
  pvals[is.na(pvals)]<-0
  
  goseq_res=setNames(data.frame(matrix(ncol = ncol(pvals), nrow = length(unique(ontology$GO)))),colnames(pvals))
  goseq_res[is.na(goseq_res)]<-""
  deseq_cnts<-setNames(data.frame(matrix(ncol=ncol(pvals),nrow=nrow(pvals))),colnames(pvals))
  rownames(deseq_cnts)<-rownames(pvals)
  ###Loop through each DESeq column to check for enrichment
  for (ii in 1:ncol(pvals)){
    genes.p=as.integer(p.adjust(abs(pvals[,ii]),method=pmethod)<pcutoff)
    genes=genes.p*lfcs[,ii]
    names(genes)=row.names(pvals)
    deseq_cnts[,ii]<-genes
    
    gene_len<-as.vector(gene_len_raw[match(names(genes),gene_len_raw[,1]),2])
    
    png(filename=paste(colnames(pvals)[ii],"goseqQC.png",sep=""),width=1500,height=900,units="px",bg="transparent")
    layout(matrix(c(1,2), 1, 2, byrow = TRUE))
    pwf=nullp(genes,bias.data=gene_len) #Gene length bias correction determination
    
    go.Wall=goseq(pwf,gene2cat=ontology) #Wallenius method gene set enrichment
    dev.off()
    #Note: the below lines are commented out because 1. after many runs, I notice Wallenius is just as good and 2. it's much faster
    #go.Samp=goseq(pwf,gene2cat=ontology,method="Sampling",repcnt=1000) #Sampling method gene set enrichment
    
    
    #plot(log10(go.Wall[,2]),log10(go.Samp[match(go.Wall[,1],go.Samp[,1]),2]), #QC plot to ensure two methods give consistent results
    #     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
    #     xlim=c(-3,0))
    #abline(0,1,col=3,lty=2)
    #dev.off()
    
    idx=match(go.Wall$category[p.adjust(go.Wall$over_represented_pvalue,method=pmethod)<pcutoff],go.Wall$category)
    goseq_res[idx,ii]<-paste(go.Wall$category[idx],format(go.Wall$over_represented_pvalue[idx],digits=3,scientific=TRUE),sep="_")
  }
 
  write.csv(as.data.frame(goseq_res), file=paste(fname,"_GoSeq_",pmethod,pcutoff,".csv",sep=""))
  write.csv(deseq_cnts, file=paste(fname,"_DESeqCnts.csv",sep=""))
}
# FGSEA gene ontology analysis --------------------------------------------
if(isTRUE(fgseaRun)){
  library(fgsea)
  library(stringr)
  library(data.table)
  library(ggplot2)
  
  pathways <- gmtPathways("C:\\Users\\User1\\Dropbox (MIT)\\00 Tools\\GOSlim.gmt")[-3]
  if(is.null(pcutoff)){pcutoff=0.05}
  
  fgsea_filter=setNames(data.frame(matrix(ncol = ncol(pvals), nrow = length(pathways))),colnames(pvals))
  fgsea_filter[is.na(fgsea_filter)]<-""
  rownames(fgsea_filter)<-names(pathways)
  fgsea_pvals=fgsea_filter
  ###Loop through comparisons
  for (jj in 1:ncol(pvals)){
    comp_label=colnames(lfcs)[jj]
    
    lfcs.ord<-lfcs[order(-lfcs[,jj]), ]
    ranks<-na.omit(lfcs.ord[,jj])
    names(ranks)<-rownames(na.omit(lfcs.ord))
    
    fgseaRes<-fgsea(pathways,ranks,minSize=8,maxSize=500,nperm=500)
    
    fwrite(fgseaRes, file=paste(comp_label,"_fgseaRes.txt",sep=""), sep="\t", sep2=c("", " ", ""))
    
    ####Plot GSEA results in table containing sig gene sets and rank plots
    #Filter out significant pathways and obtain their index in results
    topPathways<-fgseaRes[padj<pcutoff][(order(NES)),pathway]
    idxFGSEA1=match(topPathways,fgseaRes$pathway)
    #Add results to output files/plot if there are any significant pathways
    if (length(idxFGSEA1)>0){
      idxFGSEA2=match(topPathways,rownames(fgsea_filter))
      fgsea_filter[idxFGSEA2,jj]=fgseaRes$NES[idxFGSEA2]
      
      png(filename=paste(comp_label,"_fgseaTable.png",sep=""),width=1500,height=900,units="px",bg="transparent")
      plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5)
    }
    idxFGSEA3=match(rownames(fgsea_pvals),fgseaRes$pathway)
    fgsea_pvals[idxFGSEA3,jj]<-fgseaRes$padj*fgseaRes$NES/abs(fgseaRes$NES)
    if(!is.null(plotGeneSet)){
      for (kk in 1:length(plotGeneSet)){
        enrichName=paste(substr(plotGeneSet[kk],1,10),"-",colnames(fgsea_filter)[jj],"-1.pdf",sep="")
        plotEnrichment(pathways[[plotGeneSet[kk]]], ranks)
        ggsave(enrichName,last_plot(),device="pdf",width=7,height=7.5,units="in",dpi=300)
      }
    }
  }
  
  graphics.off()
  write.csv(as.data.frame(fgsea_filter), file=paste(fname,"_FGSEA_",pcutoff,"_SigNES.csv",sep=""))
  write.csv(as.data.frame(fgsea_pvals), file=paste(fname,"_FGSEA_",pcutoff,"_SignPADJ.csv",sep=""))
}
  
# Upset plot generation ---------------------------------------------------
if(isTRUE(volcanoRun)){
  library(ggplot2)
  library(ggfortify)
  source("C:\\Users\\User1\\Dropbox (MIT)\\00 Tools\\R\\ggtheme.R")
  
  sym <- read.delim("C:\\Users\\User1\\Dropbox (MIT)\\00 Tools\\R\\sym_lookup.txt", header=TRUE)
  #sym <- as.data.frame(rownames(pvals))
  #sym$two <- as.vector(rownames(pvals))
  sym <- sym[match(rownames(pvals),sym[,1]),]
  
  comps_num=substr(colnames(lfcs),1,regexpr(".v.",colnames(lfcs))-1)
  comps_denom=substr(colnames(lfcs),regexpr(".v.",colnames(lfcs))+3,regexpr("*$",colnames(lfcs)))
  ii=1
  temp=as.data.frame(sym[,2])
  for (ii in 1:ncol(lfcs)){
    temp[,2]=lfcs[,ii]
    temp[,3]=abs(pvals[,ii])
    names(temp)=c('sym','log2FoldChange','padj')
    
    pdf(paste("Volcano.",colnames(lfcs)[ii],".pdf",sep=""),width=15,height=9,useDingbats = F)
    par(mar=c(7,9,4,4)+0.1,mgp=c(6,2.5,0))
    with(temp, plot(log2FoldChange, -log10(padj), pch=19, xlim=c(-4,4),ylim=c(0,5),cex.lab=4,cex.axis=4,cex.main=4))
    title(bquote(.(comps_num[ii]) ~ "                                             "), col.main = "magenta",cex.main=4)
    title("vs", col.main = "black",cex.main=4)
    title(bquote("                                             " ~ .(comps_denom[ii])), col.main = "blue",cex.main=4)
    
    # Add colored points: colored if padj<0.01 & >1 log2foldchange)
    up.count=nrow(subset(temp, padj<.01 & log2FoldChange>1))
    down.count=nrow(subset(temp, padj<.01 & log2FoldChange< -1))
    
    if (up.count>0){
      with(subset(temp, padj<.01 & log2FoldChange>1), points(log2FoldChange, -log10(padj), pch=19, col="magenta"))
      with(subset(temp, padj<.01 & log2FoldChange>1), text(x=log2FoldChange, y=-log10(padj), labels = sym, col = "magenta", cex = 1.5, pos = 1 ))
    }
    
    if (down.count>0){
      with(subset(temp, padj<.01 & log2FoldChange< -1), points(log2FoldChange, -log10(padj), pch=19, col="blue"))
      with(subset(temp, padj<.01 & log2FoldChange< -1), text(x=log2FoldChange, y=-log10(padj), labels = sym, col = "blue", cex = 1.5, pos = 1 ))
    }
    dev.off()
  }
}
  }