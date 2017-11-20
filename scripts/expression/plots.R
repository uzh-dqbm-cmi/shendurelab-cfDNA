#!/usr/bin/env Rscript

library("optparse")

opt = parse_args(OptionParser(option_list=list(
    make_option(
        c("-p", "--protein_atlas"),
        type="character",
        default=NULL,
        help="Protein Atlas RNA table (tsv)",
        metavar="P"
    ),
    make_option(
        c("-l", "--labels"),
        type="character",
        default=NULL,
        help="Expression labels",
        metavar="L"
    ),
    make_option(
        c("-m", "--mask"),
        type="character",
        default=NULL,
        help="WPS path and file mask (with %s)",
        metavar="M"
    ),
    make_option(
        c("-r", "--reference"),
        type="character",
        default=NULL,
        help="ID of the reference individual",
        metavar="R"
    ),
    make_option(
        c("-i", "--individuals"),
        type="character",
        default=NULL,
        help="IDs of other individuals",
        metavar="I"
    ),
    make_option(
        c("-o", "--output_directory"),
        type="character",
        default=NULL,
        help="Output directory",
        metavar="O"
    )
)))

reference = opt[["reference"]]
individuals = strsplit(opt[["individuals"]], ",")[[1]]
samples = c(reference, individuals)

##################################

proteinAtlas <- read.table(opt[["protein_atlas"]],header=T,as.is=T,sep="\t")
rownames(proteinAtlas) <- proteinAtlas$GeneID

ndata <- proteinAtlas[,-1]
ndata[ndata == 0] <- NA
# at least 3 tissues with non-zero values
ndata <- ndata[apply(ndata,1,FUN=function(x) { length(which(is.na(x))) < length(x)-2 } ),] 
ndata[is.na(ndata)] <- 0.04
logndata <- log2(ndata)
dim(logndata)

##################################

tLabels <- read.table(opt[["labels"]],header=T,as.is=T,sep="\t",quote="\"")

fftColumns <- 29:52 # 160-222
selFreq <- c("193","196","199")

##################################

allfreq = paste(opt[["output_directory"]], "allfreq_correlation.pdf", sep="/")
pdf(allfreq, width=8, height=6)

for (sample in samples)
{
    fdata <- read.table(sprintf(opt[["mask"]],sample),as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    logndata2 <- logndata[fdata[,1],]
 
    res <- cor(fdata[,fftColumns],logndata2,use="pairwise.complete.obs")
 
    matplot(as.numeric(sub("X","",names(fdata[,fftColumns]))),res,type="l",xlab="Frequency",ylab="Correlation",col="darkgrey",lwd=1,lty=1,main=sprintf("%s: Correlation of intensities across tissues",sample))
    lines(as.numeric(sub("X","",names(fdata[,fftColumns]))),cor(fdata[,fftColumns],logndata2[,"NB.4"],method="pearson",use="pairwise.complete.obs"),col="black",lwd=2,type="b",pch=19,cex=0.6)
    legend("topright","NB-4",fill="black")
}
dev.off()

##################################

out_mask = paste(opt[["output_directory"]], "%s_ave193-199bp_correlation.txt", sep="/")

for (sample in samples)
{
    fdata <- read.table(sprintf(opt[["mask"]],sample),as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    logndata2 <- logndata[fdata[,1],]
 
    res <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")
    res <- data.frame(category=tLabels$Category,description=tLabels$Type,tissue=colnames(res),correlation=as.numeric(res))
    sink(sprintf(out_mask, sample))
    print(sample)
    print(res[order(res$correlation),])
    sink()
}

##################################

fdata <- read.table(sprintf(opt[["mask"]],opt[["reference"]]),as.is=T,sep="\t",header=T,comment.char="~")
colnames(fdata) <- sub("X","",colnames(fdata))
rownames(fdata) <- fdata[,1]
fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
logndata2 <- logndata[fdata[,1],]
refCorrelation <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")

out_mask = paste(opt[["output_directory"]], "%s_ave193-199bp_correlation_rank.txt", sep="/")

for (sample in samples)
{
    fdata <- read.table(sprintf(opt[["mask"]],sample),as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    logndata2 <- logndata[fdata[,1],]
 
    res <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")
    res <- data.frame(category=tLabels$Category,description=tLabels$Type,tissue=colnames(res),correlation=as.numeric(res),rankDiff=rank(refCorrelation)-rank(res))
    par(mfrow=c(2,1))
 
    sink(sprintf(out_mask, sample))
    sprintf("By correlation rank difference: %s (vs. normal)",sample)
    print(head(res[order(res$rankDiff,decreasing=T),],15))
    sprintf("By correlation: %s",sample)
    print(head(res[order(res$correlation),],15))
    sink()
}
