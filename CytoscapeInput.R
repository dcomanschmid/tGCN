##########################################################################
# prepare the output files from t_GCN.R for Cytoscape input              #
# use any t_GCN.R output html file                                       #
# convert to Cytoscape format:                                           #
#     column1 = target gene identifier                                   #
#     column2 = interaction (i.e. adjusted p-value calculated by t.GCN.R #
#     column3 = query gene identifier                                    #
#                                                                        #
# by Diana Coman Schmid, Eawag, 2015                                     #
##########################################################################

library(XML)

# set your working directory and read in any output file from t_GCN.R (e.g. guide_query_GCNfdr.html)

workDir <- "your_working_directory"
gcn <- readHTMLTable(file.path(workDir,"guide_query_GCNfdr.html"), header = TRUE)

gcn.df <- as.data.frame(gcn[[1]][1:(nrow(gcn[[1]])-1),])#ignore the last row containing the total number of significant correlations
row.names(gcn.df) <- gcn.df[,1]
gcn.df[,1] <- NULL

# convert to Cytoscape format
gcn.cyt.l <- list()
for (i in 1:ncol(gcn.df)){
  gcn.cyt <- matrix(0,nrow= nrow(gcn.df), ncol=3,byrow=T)
  gcn.cyt[,1] <- rownames(gcn.df)
  gcn.cyt[,3] <- rep(colnames(gcn.df)[i],each=nrow(gcn.df))
  gcn.cyt[,2] <- as.character(gcn.df[,i])
  gcn.cyt <- gsub(" ","",gcn.cyt)
  gcn.cyt.df <- as.data.frame(gcn.cyt)
  gcn.cyt.l[[i]] <- gcn.cyt.df
}

gcn.cyt.in <- do.call("rbind",gcn.cyt.l)
colnames(gcn.cyt.in)<- c("Target","CorrPadj","Query")

# Optional: select only significant correlations (i.e. adjusted p-values <= 0.05)
# flexibly replace "0.05" with your desired threshold

gcn.cytSig.in<-gcn.cyt.in[apply(gcn.cyt.in,1, function(x) (x[2]<=0.05)),]

# save the the Cytoscape input file

write.table(gcn.cyt.in, file.path(paste(workDir,"cytoscapeIN.txt")), row.names=FALSE,quote=FALSE)
write.table(gcn.cytSig.in, file.path(paste(workDir,"cytoscapeSigIN.txt")), row.names=FALSE,quote=FALSE)


