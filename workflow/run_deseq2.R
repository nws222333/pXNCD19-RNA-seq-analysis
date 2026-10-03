suppressPackageStartupMessages(library(DESeq2))
args_all <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args_all[grepl('^--file=',args_all)][1])
root <- normalizePath(file.path(dirname(script),'..'),winslash='/')
outdir <- file.path(root,'outputs','DESeq2')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
for(cmp in c('C600_vs_CD19C','TH2_vs_TH2C')) {
  x <- read.delim(file.path(root,'inputs','counts',paste0(cmp,'_host_counts.tsv')),check.names=FALSE,stringsAsFactors=FALSE)
  mat <- as.matrix(x[,-1,drop=FALSE]);storage.mode(mat)<-'integer';rownames(mat)<-x[[1]]
  group <- factor(c(rep('recipient',3),rep('transconjugant',3)),levels=c('recipient','transconjugant'))
  coldata <- data.frame(row.names=colnames(mat),group=group)
  dds <- DESeqDataSetFromMatrix(countData=mat,colData=coldata,design=~group)
  dds <- dds[rowSums(counts(dds))>0,]
  # These are the defaults used by the historical DESeq(dds,quiet=TRUE) call.
  dds <- DESeq(dds,test='Wald',fitType='parametric',sfType='ratio',betaPrior=FALSE,quiet=TRUE)
  res <- results(dds,contrast=c('group','transconjugant','recipient'),alpha=.01,
                 independentFiltering=TRUE,pAdjustMethod='BH',lfcThreshold=0)
  write.table(data.frame(gene_id=rownames(res),as.data.frame(res)),file.path(outdir,paste0(cmp,'_DESeq2_results.tsv')),
              sep='\t',quote=FALSE,row.names=FALSE,na='NA')
  write.table(data.frame(sample=colnames(mat),group=as.character(group),size_factor=sizeFactors(dds),mapped_host_counts=colSums(mat)),
              file.path(outdir,paste0(cmp,'_DESeq2_size_factors.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(outdir,'R_sessionInfo.txt'))
