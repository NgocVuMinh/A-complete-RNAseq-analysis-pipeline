# R script for combine featureCounts outputs into one read count matrix in R
# Adapted from Ahmed Alhendi: https://github.com/AAlhendi1707/htseq-merge/blob/master/htseq-merge_all.R 
# Usage example:  Rscript merge_counts.R


setwd("../featureCounts")
path <- "../featureCounts"
files <- list.files(path=path, pattern="*clean2.txt")
print(sprintf("Files to be merged: "))
print(files)

# using perl to manpulate file names by trimming file extension
labs <- paste("", gsub("\\.txt", "", files, perl=TRUE), sep="")
labs

# load all files to list object, use paste to return the trimpping parts to file name
print(sprintf("File read START %s", format(Sys.time(),"%b %d %Y-%H-%M")))

cov <- list()
filepath <- file.path(path,paste(labs[1],"clean2.txt",sep=""))

for (i in labs) {
  filepath <- file.path(path,paste(i,".txt",sep=""))
  cov[[i]] <- read.table(filepath,sep = "\t", header=T, stringsAsFactors=FALSE)
  colnames(cov[[i]]) <- c("Geneid", substr(i, 1, 3))
}

print(sprintf("File read END %s", format(Sys.time(),"%b %d %Y-%H-%M")))

## construct one data frame from list of data.frames using reduce function
print(sprintf("Merging START %s", format(Sys.time(),"%b %d %Y-%H-%M")))
df <- Reduce(function(x,y) merge(x = x, y = y, by ="Geneid"), cov)
print(sprintf("Merging END %s", format(Sys.time(),"%b %d %Y-%H-%M")))

write.table(df, paste(path, "featureCounts_merged",".txt",sep=""), sep="\t", quote= F, row.names = F)

print(sprintf("Merged file saved %s", format(Sys.time(),"%b %d %Y-%H-%M")))
