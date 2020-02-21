setwd("C:\\Users\\User1\\Dropbox (MIT)\\00 Tools\\R")
sym <- read.delim("sym_lookup.txt", header=TRUE)

setwd("C:\\Users\\User1\\Dropbox (MIT)\\Promoter study\\0Manuscript\\Analysis\\")
res <- read.csv("X14.v.15_fgseaRes.txt", header=TRUE)

head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, xlim=c(-4,4)))
title(expression("Y11430" * phantom(" v. X33)")), col.main = "magenta")
title(expression(phantom("Y11430") * "v. X33"), col.main = "blue")

# Add colored points: colored if padj<0.01 & >1 log2foldchange)
with(subset(res, padj<.01 & log2FoldChange>1), points(log2FoldChange, -log10(pvalue), pch=20, col="magenta"))
with(subset(res, padj<.01 & log2FoldChange>1), text(x=log2FoldChange, y=-log10(pvalue), labels = sym, col = "magenta", cex = 0.75, pos = 1 ))

with(subset(res, padj<.01 & log2FoldChange< -1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & log2FoldChange< -1), text(x=log2FoldChange, y=-log10(pvalue), labels = sym, col = "blue", cex = 0.75, pos = 1 ))


###From multiple DESeq output files
setwd("C:\\Users\\User1\\Dropbox (MIT)\\Promoter study\\0Manuscript\\Analysis\\")
sym <- read.delim("sym_lookup.txt", header=TRUE)
res <- read.delim("ATG30.Foreign.v.Native_LFC.txt", header=TRUE)
pval <- read.delim("ATG30.Foreign.v.Native_SignPADJ.txt", header=TRUE)

sym <- sym[match(pval[,1],sym[,1]),]

comps_num=substr(colnames(res[,-1]),1,regexpr("vs",colnames(res[,-1]))-1)
comps_denom=substr(colnames(res[,-1]),regexpr("vs",colnames(res[,-1]))+2,regexpr("*$",colnames(res[,-1])))

temp=as.data.frame(sym[,2])
for (ii in 1:(ncol(res)-1)){
  temp[,2]=res[,ii+1]
  temp[,3]=abs(pval[,ii+1])
  names(temp)=c('sym','log2FoldChange','padj')
  
  with(temp, plot(log2FoldChange, -log10(padj), pch=20, xlim=c(-4,4)))
  title(bquote(.(comps_num[ii]) ~ "                                             "), col.main = "magenta")
  title("vs", col.main = "black")
  title(bquote("                                             " ~ .(comps_denom[ii])), col.main = "blue")
  
  # Add colored points: colored if padj<0.01 & >1 log2foldchange)
  up.count=nrow(subset(temp, padj<.01 & log2FoldChange>1))
  down.count=nrow(subset(temp, padj<.01 & log2FoldChange< -1))
  
  if (up.count>0){
  with(subset(temp, padj<.01 & log2FoldChange>1), points(log2FoldChange, -log10(padj), pch=20, col="magenta"))
  with(subset(temp, padj<.01 & log2FoldChange>1), text(x=log2FoldChange, y=-log10(padj), labels = sym, col = "magenta", cex = 0.75, pos = 1 ))
  }
  
  if (down.count>0){
  with(subset(temp, padj<.01 & log2FoldChange< -1), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
  with(subset(temp, padj<.01 & log2FoldChange< -1), text(x=log2FoldChange, y=-log10(padj), labels = sym, col = "blue", cex = 0.75, pos = 1 ))
  }
}
