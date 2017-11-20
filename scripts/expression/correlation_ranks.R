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
        c("-i", "--individual"),
        type="character",
        default=NULL,
        help="ID of the target individual",
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
sample = opt[["individual"]]

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

fdata <- read.table(sprintf(opt[["mask"]],opt[["reference"]]),as.is=T,sep="\t",header=T,comment.char="~")
colnames(fdata) <- sub("X","",colnames(fdata))
rownames(fdata) <- fdata[,1]
fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
logndata2 <- logndata[fdata[,1],]
refCorrelation <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")

fdata <- read.table(sprintf(opt[["mask"]],sample),as.is=T,sep="\t",header=T,comment.char="~")
colnames(fdata) <- sub("X","",colnames(fdata))
rownames(fdata) <- fdata[,1]
fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
logndata2 <- logndata[fdata[,1],]

res <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")
res <- data.frame(category=tLabels$Category,description=tLabels$Type,tissue=colnames(res),correlation=as.numeric(res),rankDiff=rank(refCorrelation)-rank(res))
par(mfrow=c(2,1))

file_mask = paste("%s_by_", opt[["reference"]], "_ave193-199bp_correlation_rank.txt", sep="")
out_mask = paste(opt[["output_directory"]], file_mask, sep="/")

sink(sprintf(out_mask, sample))
sprintf("By correlation rank difference: %s (vs. normal)",sample)
print(head(res[order(res$rankDiff,decreasing=T),],15))
sprintf("By correlation: %s",sample)
print(head(res[order(res$correlation),],15))
sink()
