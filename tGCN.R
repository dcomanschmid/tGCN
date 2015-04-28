###############################################################################################################
# t_GCN.r script for computing targeted gene co-expression networks #
# 1. "guide-query GCN" --- computes correlations between guide and query genes #
# We are grateful to Marc W. Schmid, UZH, schmid.m@access.uzh.ch for improving the speed of this step #
# #
# 2. "guide-query groups GCN" --- computes correlations between guide genes and groups of query genes #
###############################################################################################################

# Required packages and libraries

source("http://bioconductor.org/biocLite.R")
biocLite("multtest")
library(Biobase)
library(multtest)
library(xtable)

################################################################################
# 1. "guide-query GCN" --- computes correlations between guide and query genes #
################################################################################

# set the threshold for P-values

alpha <- 0.05

# read in the file containing the expression data for the guide genes

workDir <- "your_working_directory"

guide_expression <-read.csv(file.path(workDir, "guide_expression.csv"),header=TRUE, row.names=1)
guide_expression_names <- colnames(guide_expression)

# read in the file containing the expression data for the query genes

query_expression <- read.csv(file.path(workDir,"query_expression.csv"),header=TRUE, row.names=1)
m <- ncol(query_expression)
n <- nrow(query_expression)
cornames <- colnames(query_expression)
data <- cbind(guide_expression, query_expression)

# calculate the Pearson correlation coefficients between guide genes and each group of query genes
# to use other similarity measures change for example "pearson" to "spearman"

corData <- cor(data, method ="pearson")[1:ncol(guide_expression),(ncol(guide_expression)+1):(m+ncol(guide_expression))]

# asses the statistical significance of positive Pearson correlation coefficients

all.pValue <- apply(corData, 1, function(x) 1-pnorm(sqrt(n-3)*0.5*(log(1+x)-log(1-x))))

# to assess the statistical significance of negative Pearson correlation coefficients run the entire script with the
# following modificantion: all.pValue <-apply(corData, 1, function(x) pnorm(sqrt(n-3)*0.5*(log(1+x)-log(1-x))))

all.pValue <- as.vector(all.pValue)

# Apply multiple testing correction (Bonferroni, Holm, FDR)

bonf <- mt.rawp2adjp(all.pValue, proc="Bonferroni")
pValue.bonf <- bonf$adjp[,2][order(bonf$index)]
pValue.bonf <- matrix(pValue.bonf,nrow=m)
pValue.bonf <-rbind(pValue.bonf,colSums((pValue.bonf<=alpha)))

holm <- mt.rawp2adjp(all.pValue, proc="Holm")
pValue.holm <- holm$adjp[,2][order(holm$index)]
pValue.holm <- matrix(pValue.holm,nrow=m)
pValue.holm <-rbind(pValue.holm,colSums((pValue.holm<=alpha)))

fdr <- mt.rawp2adjp(all.pValue, proc="BH")
pValue.fdr <- fdr$adjp[,2][order(fdr$index)]
pValue.fdr <- matrix(pValue.fdr,nrow=m)
pValue.fdr <-rbind(pValue.fdr,colSums((pValue.fdr<=alpha)))

colnames(pValue.bonf) <- colnames(pValue.holm) <- colnames(pValue.fdr) <- guide_expression_names
rownames(pValue.bonf) <- rownames(pValue.holm) <- rownames(pValue.fdr) <- c(cornames,"nr. signif.")

# save the results as HTML tables in your working directory)

pValue.bonf <- xtable(pValue.bonf,digits=4)
print(pValue.bonf, type="html", file=paste(workDir,'guide_query_GCN',"bonf",".html",sep=""))

pValue.holm <- xtable(pValue.holm,digits=4)
print(pValue.holm, type="html", file=paste(workDir,'guide_query_GCN',"holm",".html",sep=""))

pValue.fdr <- xtable(pValue.fdr,digits=4)
print(pValue.fdr, type="html", file=paste(workDir,'guide_query_GCN',"fdr",".html",sep=""))

#######################################################################################################
# 2. "guide-query groups GCN" --- computes correlations between guide genes and groups of query genes #
#######################################################################################################

# set the threshold for P-values

alpha <- 0.05

# read in the file containg the names of the query genes groups

query_groups_names <- read.csv(file.path(workDir,"query_groups_names.csv"),header=TRUE)
query_groups_names <- as.matrix(query_groups_names)
query_groups_names_L <- length(query_groups_names)

# read in the file containing the expression data for the guide genes

guide_expression <-read.csv(file.path(workDir,"guide_expression.csv"),header=TRUE, row.names=1)
guide_expression_names <- colnames(guide_expression)

HC.vec.all <- sig.num.bonf <- sig.num.holm <- sig.num.fdr <- numeric(0)

for(j in 1:query_groups_names_L){
  queryname <- query_groups_names[j]
  query_group_expression <- read.csv(file=paste(workDir,queryname,".csv", sep=""),header=TRUE, row.names=1)
  m <- ncol(query_group_expression)
  n <- nrow(query_group_expression)
  cornames <- colnames(query_group_expression)
  all.pValue <- numeric(0)
  pV.HC <- numeric(ncol(guide_expression))
  for(i in 1:ncol(guide_expression)){
    data <- cbind(guide_expression, query_group_expression)
    # calculate the Pearson correlation coefficients between guide genes and each group of query genes
    cordata <- cor(data, method = "pearson")[i,(ncol(guide_expression)+1):(m+ncol(guide_expression))] # to use other similarity measures change for example "pearson" to "spearman"
    # asses the statistical significance of positive Pearson correlation coefficients
    pValue <- (1-pnorm(sqrt(n-3)*0.5*(log(1+cordata)-log(1-cordata))))
    ##############################################################################################################################################
    # to assess the statistical significance of negative Pearson correlation coefficients run the entire script with the following modificantion #
    # pValue <- pnorm(sqrt(n-3)*0.5*(log(1+cordata)-log(1-cordata))) #
    ##############################################################################################################################################
    N <- length(pValue)
    I <- sum(pValue<=alpha)
    pV.HC[i] <- sqrt(N)*(I/N-0.05)/(sqrt(0.05*(1-0.05)))
    all.pValue <- c(all.pValue,pValue)
  }
  # Apply multiple testing correction (Bonferroni, Holm, FDR)
  bonf <- mt.rawp2adjp(all.pValue, proc="Bonferroni")
  pValue.bonf <- bonf$adjp[,2][order(bonf$index)]
  pValue.bonf <- matrix(pValue.bonf,nrow=m)
  sig.num.bonf <- c(sig.num.bonf,colSums((pValue.bonf<=alpha)))
  pValue.bonf <-rbind(pValue.bonf,colSums((pValue.bonf<=alpha)),pV.HC)
  
  holm <- mt.rawp2adjp(all.pValue, proc="Holm")
  pValue.holm <- holm$adjp[,2][order(holm$index)]
  pValue.holm <- matrix(pValue.holm,nrow=m)
  sig.num.holm <- c(sig.num.holm,colSums((pValue.holm<=alpha)))
  pValue.holm <-rbind(pValue.holm,colSums((pValue.holm<=alpha)),pV.HC)

  fdr <- mt.rawp2adjp(all.pValue, proc="BH")
  pValue.fdr <- fdr$adjp[,2][order(fdr$index)]
  pValue.fdr <- matrix(pValue.fdr,nrow=m)
  sig.num.fdr <- c(sig.num.fdr,colSums((pValue.fdr<=alpha)))
  pValue.fdr <-rbind(pValue.fdr,colSums((pValue.fdr<=alpha)),pV.HC)
  colnames(pValue.bonf) <- colnames(pValue.holm) <- colnames(pValue.fdr) <- guide_expression_names
  rownames(pValue.bonf) <- rownames(pValue.holm) <- rownames(pValue.fdr) <- c(cornames,"nr. signif.","Tukey test")
  
  # save the results as HTML tables in your working directory
  pValue.bonf.x <- xtable(pValue.bonf,digits=4)
  print(pValue.bonf.x, type="html", file=paste(workDir,queryname,"bonf",".html",sep=""))
  pValue.holm.x <- xtable(pValue.holm,digits=4)
  print(pValue.holm.x, type="html", file=paste(workDir,queryname,"holm",".html",sep=""))
  pValue.fdr.x <- xtable(pValue.fdr,digits=4)
  print(pValue.fdr.x, type="html", file=paste(workDir,queryname,"fdr",".html",sep=""))
}
