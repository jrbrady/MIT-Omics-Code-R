deseqMultifactor <- function(intct,designTable,outputName=NULL,designName,cond1=NULL,num,denom){
"deseqMultifactor.R - The purpose of this function is to perform multi-factor DESeq2 analysis on Integer Count RNASeq data
This is achieved by providing a designTable that categorizes samples for dispersion calcs, then looping through desired contrasts"


#User inputs - note that everything initialized with a default value is an optional argument, and all others are required
#intct is a matrix with rownames GQ.ID and columns sample data
#designTable is a data frame where each rowname is a different sample matching the column names of intct and columns are categorical
    #IMPORTANT: one of the columns should be the designName (see below) and should be a concatenation of two or more other columns, 
    #separated by a period (ex. strain: Y11430, timept: 24M, strain.tpt: Y11430.24M)

#outputName is string denoting name for output result files - defaults to designName
#designName is the name of the column in designTable by which to group samples for dispersion calcs preceded by ~, ex. ~strain 

#cond1 will determine the outermost 'for' loop and is the name of a column in designTable 
    #Contrasts will only be made within a single condition (e.g. only strains at the same timept will be compared to each other)
    #In the event that only one condition was tested, this argument is optional and will not be used.

#num specifies a list of identifiers (such as strain name) to be in the numerator of the contrast
#denom specifies a list of identifiers (such as strain name, same identifiers as num) to be in the denominator of the contrast

#Example, designName = "strain.tpt", cond1="timept", num=c("Y11430","Y7556"), denom=c("X33") will 
    #compare Y11430.24M_vs_X33.24M and Y7556.24M_vs_X33.24M and so on for other timepts (0G, 24G in this example)
    


# DDS ---------------------------------------------------------------------
if (is.null(outputName)){outputName=as.character(designName)[-1]}
  
options(warn=-1)
library("DESeq2")
library("ashr")
library("rafalib")

dds<-DESeqDataSetFromMatrix(intct, colData=designTable,designName)
dds<-DESeq(dds)
#resultsNames(dds)


# LFCShrink ---------------------------------------------------------------
#Initialize data frames for results
resMF<-as.data.frame(results(dds))
res_SignPADJ<-resMF[0]
res_LFC<-resMF[0]
counter=0

#If condition 1 entered, loop over unique conditions as listed in designTable; else cond1 is blank string
if (is.null(cond1)){con1Vec=""} else {con1Vec=paste(".",cond1,sep="")}

#Loop through all pairwise combinations specified by cond1, num, and denom
for (con1 in con1Vec){      #Loop over condition 1
  for (den1 in denom){      #Loop over denominators
    for (num1 in num){      #Loop over numerators
      if (den1==num1){
        next
      }
      counter=counter+1
      #Perform lfcShrink & contrast on all pairwise combos
      temp<-as.data.frame(lfcShrink(dds,contrast=c(as.character(designName)[-1],paste(num1,con1, sep=""),paste(den1,con1, sep="")),type="normal"))
    
      #Store LFC and padj*(sign of LFC) in separate data frames for heatmap visualization and subsequent analysis
      res_LFC[,counter]<-temp$log2FoldChange
      colnames(res_LFC)[counter]<-paste(num1,".v.",den1,con1, sep="")
    
      res_SignPADJ[,counter]<-sign(temp$log2FoldChange)*temp$padj
      colnames(res_SignPADJ)[counter]<-paste(num1,".v.",den1,con1, sep="")
    }
  }  
}

write.table(res_SignPADJ, file=paste(outputName,"_SignPADJ.txt",sep=""),sep="\t",col.names=NA)
write.table(res_LFC, file=paste(outputName,"_LFC.txt",sep=""),sep="\t",col.names=NA)
}