## differential gene expression
## Please use gene symbol
## Biological replicates for each sample; no techinical ones, if there are techinical replicates, please combine them. 


args <- commandArgs(trailingOnly=TRUE)
expression.file <- args[1]
class.file <- args[2]
exp.format <- args[3]
out.path <- args[4]
err.file <- args[5]

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
expression.data <- read.table(expression.file,header=T,row.names=1,sep="\t",check.names=FALSE)
myclass = read.table(class.file,row.names=1,check.names=FALSE)   

if (dim(expression.data)[2] != dim(myclass)[1]){
	write("Error: Differential gene expression analysis - number of samples in expression does not qual to number of class labels...",file=err.file)
	stop("Number of samples in expression does not qual to number of class labels...")
}
if (length(unique(myclass[,1])) !=2 || length(grep("benign",myclass[,1],ignore.case=TRUE)) == 0){
	write("Error: Differential gene expression analysis - sample classes are not correct...",file=err.file)
	stop("Sample classes are not correct...")
}

myfactors = data.frame(Tissue=myclass[colnames(expression.data),1])
if(exp.format == "rpkm"){

	if(is.installed("NOISeq")){	
	}else{
		source("http://bioconductor.org/biocLite.R")
		biocLite("NOISeq")
	}

	suppressMessages(library(NOISeq))
	mydata <- readData(data=expression.data,factors=myfactors)
	
	if (length(myfactors[,1]) == 2){
		myresults <- noiseq(mydata, factor = "Tissue", k = NULL, norm = "n", replicates = "no")
		mynoiseq.deup = rownames(degenes(myresults,q=0.9,M="up"))
		mynoiseq.dedown = rownames(degenes(myresults,q=0.9,M="down"))
		pdf(file=paste(out.path,"/DE.pdf",sep=""))
		DE.plot(myresults,q=0.9,graphic="expr",log.scale=TRUE)
		dev.off()
	}else{
		myresults = noiseqbio(mydata,factor="Tissue",norm="n",r=20, k=0, a0per = NULL)     
		mynoiseq.deup = rownames(degenes(myresults,q=0.95,M="up"))
		mynoiseq.dedown = rownames(degenes(myresults,q=0.95,M="down"))
		pdf(file=paste(out.path,"/DE.pdf",sep=""))
		DE.plot(myresults,q=0.95,graphic="expr",log.scale=TRUE)
		dev.off()
	}

	condition.1=gsub("_mean","",colnames(as.matrix(myresults@results[[1]]))[1])
	condition.2=gsub("_mean","",colnames(as.matrix(myresults@results[[1]]))[2])
	
	output = rbind(cbind(mynoiseq.deup,rep(condition.1, length(mynoiseq.deup))), cbind(mynoiseq.dedown, rep(condition.2,length(mynoiseq.dedown))))
	write.table(output, file=paste(out.path,"/DE.gene.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
	
}else{
	if(is.installed("DESeq")){
	}else{
	    source("http://bioconductor.org/biocLite.R")
    	biocLite("DESeq")
	}
	suppressMessages(library(DESeq))
	condition = factor(myfactors[,"Tissue"])
	mydata <- newCountDataSet(expression.data,condition)
	mydata = estimateSizeFactors( mydata )

	if (length(myfactors[,1]) == 2){
		mydata = estimateDispersions( mydata, method="blind", sharingMode="fit-only" )
		myresults = nbinomTest(mydata,unique(condition)[1],unique(condition)[2])
		pdf(file=paste(out.path,"/DE.pdf",sep=""))
			plotMA(myresults, col = ifelse(myresults$padj>=0.1, "gray32", "red3"))
		dev.off()
		resSig = myresults[ (myresults$padj <= 0.1 & is.na(myresults$padj)==FALSE), ]
		my.deup = resSig[resSig$foldChange<1,]
		my.dedown = resSig[resSig$foldChange>1,]
	}else{
		mydata = estimateDispersions( mydata )
		myresults = nbinomTest(mydata,unique(condition)[1],unique(condition)[2])
		pdf(file=paste(out.path,"/DE.pdf",sep=""))
			plotMA(myresults, col = ifelse(myresults$padj>=0.05, "gray32", "red3"))
		dev.off()
		resSig = myresults[ (myresults$padj <= 0.05 & is.na(myresults$padj)==FALSE), ]
		my.deup = resSig[resSig$foldChange<1,]
		my.dedown = resSig[resSig$foldChange>1,]
	}
	
	output = rbind(cbind(my.deup$id,as.character(rep(unique(condition)[1], dim(my.deup)[1])),my.deup$foldChange,my.deup$padj), cbind(my.dedown$id, as.character(rep(unique(condition)[2],dim(my.dedown)[1])),my.dedown$foldChange,my.dedown$padj))
	write.table(output, file=paste(out.path,"/DE.gene.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
}


