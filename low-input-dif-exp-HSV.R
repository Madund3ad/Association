library("org.Mm.eg.db")
library(plyr)
library(dplyr)
library(DESeq2)
library(DEGreport)
library(tximport)

#krab proteins
KRABS=read.csv("E:/DATA/mm10/KRABS.csv",header = F)
KRABS
getgenename<-function(ids){
  tab=read.csv("E:/DATA/mm10/ENSEMBL-gene-transcripts-info.tsv",sep="\t")
  genenames=tab$Gene.name[match(ids,tab$Gene.stable.ID)]
  return(genenames)
}
getgenenameV<-function(ids){
  tab=read.csv("E:/DATA/mm10/ENSEMBL-gene-transcripts-info.tsv",sep="\t")
  genenames=tab$Gene.name[match(ids,tab$Gene.stable.ID.version)]
  return(genenames)
}
getgenenameT<-function(ids){
  tab=read.csv("E:/DATA/mm10/ENSEMBL-gene-transcripts-info.tsv",sep="\t")
  genenames=tab$Gene.name[match(ids,tab$Transcript.stable.ID.version)]
  return(genenames)
}
tab=read.csv("E:/DATA/mm10/ENSEMBL-gene-transcripts-info.tsv",sep="\t")
t1=tab[grep(tab$Gene.name, pattern = "mt-"),]
t=t1[c(4,11),]
t1=t1[c(4,11,which(t1$Gene.type=="protein_coding")),]

tx2gene<-read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)]

mit_genes=t1[!duplicated(t1$Gene.stable.ID),]
#mit_genes=rbind(mit_genes,t)

makedata<-function(filenames){
  mm_time=vector()
  mm_virus=vector()
  mm_virus[grep("HSV",filenames)]="HSV"
  mm_virus[grep("IAV",filenames)]="IAV"
  mm_virus[grep("untreat",filenames)]="untreated"
  ab=vector()
  ab[grep("FLAG",filenames)]="FLAG"
  ab[grep("IgG",filenames)]="IgG"
  ab[grep("J2",filenames)]="J2"
  ab[grep("Z22",filenames)]="Z22"

  ab[grep("_input_",filenames)]="input"
  
  mm_time[grep("12h",filenames)]="12h"
  mm_time[grep("8h",filenames)]="8h"
  mm_time[grep("untreat",filenames)]="0h"
  exp=vector()
  exp[grep("low",filenames)]="low"
  exp[grep("total",filenames)]="total_rna"
  #nested=as.factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15))
  length(exp)
  length(mm_time)
  length(ab)
  length(mm_virus)
  colData=data.frame(as.factor(ab),as.factor(mm_time),as.factor(mm_virus),as.factor(exp))
  row.names(colData)=filenames
  colnames(colData)=c("protein","time","virus","exp")
  rownames(colData)=gsub(pattern = "(?:.*salmon/|/quant.genes.sf)", replacement = "", rownames(colData))
  return(colData)
}
#IgG vs Z22 for every group##################################################################################################
htimes=c("8h")
conditions=c("HSV")
condition=conditions[1]
htime=htimes[1]

library("tximport")
lf=list.files("E:/projects/IAV/host/", recursive = T, full.names = T)
lf=list.files("E:/projects/IAV/host/", "*.sf", recursive = T, full.names = T)

fdir= "E:/projects/IAV/host/results/IgGvsZ22/"
#WriteIgGvsZ22ForEveryGroup########################################################
WriteIgGvsZ22ForEveryGroup<-function(fdir, lf, usebl){
  htimes=c("untreat","8h")
  conditions=c("IAV")
  condition=conditions[1]
  htime=htimes[1]
  for (condition in conditions){
    for (htime in htimes){
      print(htime)
      print(condition)
      #print(lf)
      filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|Z22).*quant.sf"),lf)]
      #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|Z22).*quant.genes"),lf)]
      filenames=lf[grep(pattern = paste0(".*",htime,".*(?:IgG|Z22).*quant.genes"),lf)]
      lf
      filenames
      colData=makedata(filenames)
      colData
      ss=which(colData$virus==condition & colData$time==htime& colData$protein!="Z22")
      sscolData=colData
      txi = tximport(filenames, txIn = FALSE, type="salmon", geneIdCol = "Name")
      #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
      ss=grep(paste0(condition,".",htime,".*(?:IgG|Z22)"),rownames(colData))
      
      
      dds <- DESeqDataSetFromTximport(txi,
                                      colData = sscolData, #design = ~virus+protein)
                                      design = ~protein)
      #KRABS1######################################################
      
      KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
      KRABS=unique(KRABS[,1])
      
      #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
      #KRABS[,2]
      #KRABS=KRABS[,2]
      KRABS=gsub(pattern = "\\s",replacement="",KRABS)
      
      counts <- counts(dds, normalized = FALSE)
      acounts=rownames(counts)
      acounts=gsub(pattern = "\\..*",replacement="",acounts)
      mcounts=match(KRABS,acounts)
      #counts[na.omit(m),]
      
      
      #################################
      if (usebl){
        
        blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
        blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
        blackliist=unique(blackliist)
        
        bl2=match(blackliist,rownames(counts(dds)))
        bl2=na.omit(bl2)
        
        dds <- dds[-na.omit(bl2),]
      }
      #################################################
      nm <- assays(dds)[["avgTxLength"]]
      #sf <- estimateSizeFactorsForMatrix(counts(dds))
      dim(counts(dds))
      #blackliist
      rownames(counts(dds))
      
      #dds
      keep <- rowSums(counts(dds)) > 10
      #keep <- rowMeans(counts(dds)) > 10
      dds <- dds[keep,]
      dds$protein <- relevel(dds$protein, ref = "IgG")
      dds <- DESeq(dds)
      #plotDispEsts(dds)
      #title(paste0(condition,"  ",htime))
      
      res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
      ####################################
      a=rownames(res)
      a=gsub(pattern = "\\..*",replacement="",a)
      m=match(KRABS,a)
      #m
      #res[na.omit(m),]
      KRABres=res[na.omit(m),c(2,6)]
      
      
      a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
      dim(counts[na.omit(mcounts),])
      colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
      ##################################  
      summary(res)
      #plotMA(res, ylim=c(-10,10))
      res=res[order(res$padj),]
      ensemblsIDS=row.names(res)
      #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
      res$symbols=getgenename(ensemblsIDS)
      res 
      # counts=counts[order(counts[,1]),]
      #df=sort(counts[,1],decreasing = T)
      #getgenename(names(df)[1:20])
      #names(df)[1:20]
      # colSums(counts)
      # tail(counts)
      counts <- counts(dds, normalized = FALSE)
      counts=data.frame(counts)
      #df[with(df, order(1)), ]
      #
      #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
      #mitcounts$names=getgenename(rownames(mitcounts))
      #mitcounts
      counts <- counts(dds, normalized = TRUE)
      counts=data.frame(counts)
      #KRABS2##########################
      counts <- counts(dds, normalized = T)
      
      acounts=rownames(counts)
      acounts=gsub(pattern = "\\..*",replacement="",acounts)
      mcounts=match(KRABS,acounts)
      #counts[na.omit(mcounts),]
      row.names(a)=a$Row.names
      a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
      a$Row.names=NULL
      dim(counts[na.omit(mcounts),])
      row.names(a)=a$Row.names
      a$Row.names=getgenename(a$Row.names)
      #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
      #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
      ################################
      mit_genes
      mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
      mitcounts=data.frame(counts[mc,])
      mitcounts$names=getgenename(rownames(mitcounts))
      colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
      print(mitcounts)
      
      
      mitres=res[mit_genes$Gene.stable.ID,]
      mitres=na.omit(mitres)
      mitres$log2FoldChange=round(mitres$log2FoldChange,2)
      
      #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
      mitres
      #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.Z22_R1","IAV.Z22_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
      paste0(fdir,"/",htime,"mittable005.csv")
      write.csv(mitres,paste0(fdir,"/",htime,"mittable005.csv"),row.names = F)
      #############
      library(gprofiler2)
      lg=rownames(res[which(res$padj<0.05),])
      lg=gsub(lg, pattern = "\\..*", replacement="")
      if (length(lg)>0){  
        gostres <- gost(query = lg,
                        organism = "mmusculus",user_threshold = 0.05,evcodes = T,
                        correction_method ="fdr")
        if (!is.null(gostres)){
          #png(paste0(fdir,condition,"/",htime,"KEGG.png"))
          if (dim(gostres$result[gostres$result$source=="KEGG",])[1]){
            if (dim(gostres$result[gostres$result$source=="KEGG",])[1]<20){
              publish_gosttable(gostres, 
                                highlight_terms = gostres$result[gostres$result$source=="KEGG",],
                                use_colors = TRUE, 
                                show_columns = c("source", "term_name", "term_size", "intersection_size"),
                                filename =paste0(fdir,"/",htime,"KEGG.png") )}
            else{
              publish_gosttable(gostres, 
                                highlight_terms = gostres$result[gostres$result$source=="KEGG",][1:20,],
                                use_colors = TRUE, 
                                show_columns = c("source", "term_name", "term_size", "intersection_size"),
                                filename =paste0(fdir,"/",htime,"KEGG.png") )
            }
          }
          if (dim(gostres$result[gostres$result$source=="REAC",])[1]){
            if (dim(gostres$result[gostres$result$source=="REAC",])[1]<20){
              publish_gosttable(gostres, 
                                highlight_terms = gostres$result[gostres$result$source=="REAC",],
                                use_colors = TRUE, 
                                show_columns = c("source", "term_name", "term_size", "intersection_size"),
                                filename = paste0(fdir,"/",htime,"REAC.png"))
            }
            else{
              t=gostres$result[gostres$result$source=="REAC",]
              publish_gosttable(gostres, 
                                highlight_terms = gostres$result[gostres$result$source=="REAC",][1:20,],
                                use_colors = TRUE, 
                                show_columns = c("source", "term_name", "term_size", "intersection_size"),
                                filename = paste0(fdir,"/",htime,"REAC.png"))
            }
          }
          dfgo=data.frame(gostres$result$term_name,gostres$result$p_value,gostres$result$source,gostres$result$query_size,gostres$result$intersection_size,gostres$result$intersection)
          #colnames(dfgo)=c("Term name","FDR","source","DE genes in term")
          write.csv(dfgo ,paste0(fdir,"/",htime,"GO.csv"),row.names = F) 
        }
      }
      
      # counts
      # summary(counts)
      # counts[c(which(counts[,3]==max(counts[,3])),which(counts[,3]==max(counts[,3]))+1),]
      # row.names(counts[which(counts[,3]==max(counts[,3])),])
      #geneCounts<-plotCounts(dds_lrt_time, gene=rownames(res)[1], intgroup = c("protein","time"),returnData = T)
      
      #dim(res)
      #res
      # counts=round(counts[rownames(res),],2)
      # colnames(counts)=gsub(pattern = "HSV.|IAV.|1_",replacement = "",colnames(counts))
      # colnames(counts)=gsub(pattern = "mock",replacement = "",colnames(counts))
      #res=cbind(res[,c(2,5,6,7)],counts[rownames(res),])
      #write.csv(res,paste0(fdir,condition,"/",htime,"table_all.csv"))
      
###################################################################
      write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,"/",htime,"table005.csv"))
      write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
    }
  }
}
lf=list.files("E:/projects/IAV/host/", recursive = T, full.names = T)
fdir= "E:/projects/IAV/results/IgGvsZ22/"

#lf=list.files("E:/DATA/mm10/reseq/HSV/", recursive = T, full.names = T)
#fdir="E:/DATA/mm10/reseq/results/IgGvsZ22_for_every_group/"
WriteIgGvsZ22ForEveryGroup(fdir=fdir, lf=lf, usebl = F)
#lf=list.files("E:/DATA/mm10/reseq/untreated/", recursive = T, full.names = T)
#lf
#WriteIgGvsZ22ForEveryGroup(fdir=fdir, lf=lf, usebl = F)





########################
#############################################

cluster_rlog <-assay(rlog(dds, blind=TRUE))[rownames(res[which(res$padj<0.05),]),]


#cluster_rlog =cluster_rlog [,c(5,6,3,4,1,2,11,12,9,10,7,8)]
#a[rownames(res), ]
sscolData=colData
sscolData$time=ordered(sscolData$time,levels=c("0h",  "8h",  "12h"))

#colData[colData=="0h"]="00h"
#colData$time=factor(colData$time, levels = c("0h", "8h", "12h"))
#clusters <- degPatterns(cluster_rlog, metadata = sscolData, time="time", col="protein")
clusters <- degPatterns(cluster_rlog, metadata =sscolData[c(5,6,3,4,1,2,11,12,9,10,7,8),], time="time", col="protein")

cluster_groups <- clusters$df
for (j in unique(cluster_groups$cluster)){
  group <- clusters$df %>%
    filter(cluster == j)
  dir.create(paste0(fdir,"/",tcontrasts[i],"/group",j))
  ensemblsIDS=row.names(group)
  symbols <-getgenename(ensemblsIDS)
  #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
  group$symbols=symbols
  group
  
  write.csv(res[rownames(group),],paste0(fdir,"/",tcontrasts[i],"/group",j,".csv"))
  
  lg=ensemblsIDS
  
  if (length(lg)>0){  
    gostres <- gost(query = lg,
                    organism = "mmusculus",user_threshold = 0.05,evcodes = T,
                    correction_method ="fdr")
    if (!is.null(gostres)){
      #png(paste0(fdir,"/",tcontrasts[i],"/group",i,"KEGG.png"))
      if (dim(gostres$result[gostres$result$source=="KEGG",])[1]){
        if (dim(gostres$result[gostres$result$source=="KEGG",])[1]<20){
          publish_gosttable(gostres, 
                            highlight_terms = gostres$result[gostres$result$source=="KEGG",],
                            use_colors = TRUE, 
                            show_columns = c("source", "term_name", "term_size", "intersection_size"),
                            filename =paste0(fdir,"/",tcontrasts[i],"/group",j,"/KEGG.png") )}
        else{
          publish_gosttable(gostres, 
                            highlight_terms = gostres$result[gostres$result$source=="KEGG",][1:20,],
                            use_colors = TRUE, 
                            show_columns = c("source", "term_name", "term_size", "intersection_size"),
                            filename =paste0(fdir,"/",tcontrasts[i],"/group",j,"/KEGG.png") )
        }
      }
      if (dim(gostres$result[gostres$result$source=="REAC",])[1]){
        if (dim(gostres$result[gostres$result$source=="REAC",])[1]<20){
          publish_gosttable(gostres, 
                            highlight_terms = gostres$result[gostres$result$source=="REAC",],
                            use_colors = TRUE, 
                            show_columns = c("source", "term_name", "term_size", "intersection_size"),
                            filename = paste0(fdir,"/",tcontrasts[i],"/group",j,"/REAC.png"))
        }
        else{
          t=gostres$result[gostres$result$source=="REAC",]
          publish_gosttable(gostres, 
                            highlight_terms = gostres$result[gostres$result$source=="REAC",][1:20,],
                            use_colors = TRUE, 
                            show_columns = c("source", "term_name", "term_size", "intersection_size"),
                            filename = paste0(fdir,"/",tcontrasts[i],"/group",j,"/REAC.png"))
        }
      }
      dfgo=data.frame(gostres$result$term_name,gostres$result$p_value,gostres$result$source,gostres$result$query_size,gostres$result$intersection_size,gostres$result$intersection)
      #colnames(dfgo)=c("Term name","FDR","source","DE genes in term")
      write.csv(dfgo ,paste0(fdir,"/",tcontrasts[i],"/group",j,"/GO.csv"),row.names = F) 
    }
  }
  
}

##########################
lf=list.files("E:/DATA/mm10/input/mm39/newdata/", recursive = T, full.names = T)
fdir="E:/DATA/mm10/input/mm39/newdata/IAV/IgGvsZ22/"
#############################


#input#################################################
lf=list.files("E:/DATA/mm10/reseq/HSV", recursive = T, full.names = T,".sf")
lf
WriteInputvsZ22ForEveryGroup<-function(fdir, lf, usebl){
  htimes=c("untreat","8h")#,"12h")
  conditions=c("HSV")
  
  condition=conditions[1]
  htime=htimes[1]
  
  for (condition in conditions){
    for (htime in htimes){
      print(htime)
      print(condition)
      
      
      #print(lf)
      #filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|Z22).*quant.sf"),lf)]
      #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|Z22).*quant.genes"),lf)]
      lf
      filenames=lf[grep(pattern = paste0(".*",htime,".*low.*(?:_input_|Z22).*quant.genes.sf"),lf)]
      #filenames=lf[grep(pattern = paste0(".*",htime,".*low.*(?:_input_|Z22).*quant.genes.sf"),lf)]
      
      filenames
      colData=makedata(filenames)
      colData
      ss=which(colData$virus==condition & colData$time==htime& colData$protein!="Z22")
      sscolData=colData
      txi = tximport(filenames, txIn = FALSE, type="salmon", geneIdCol = "Name")
      dim(txi$abundance)
      #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
      ss=grep(paste0(condition,".",htime,".*(?:IgG|Z22)"),rownames(colData))
      
      dds <- DESeqDataSetFromTximport(txi,
                                      colData = sscolData, #design = ~virus+protein)
                                      design = ~protein)
      #KRABS1######################################################
      
      KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
      KRABS=unique(KRABS[,1])
      
      #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
      #KRABS[,2]
      #KRABS=KRABS[,2]
      KRABS=gsub(pattern = "\\s",replacement="",KRABS)
      
      counts <- counts(dds, normalized = FALSE)
      acounts=rownames(counts)
      acounts=gsub(pattern = "\\..*",replacement="",acounts)
      mcounts=match(KRABS,acounts)
      #counts[na.omit(m),]
      
      
      #################################
      if (usebl){
        
        blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
        blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
        blackliist=unique(blackliist)
        
        bl2=match(blackliist,rownames(counts(dds)))
        bl2
        bl2=na.omit(bl2)
        dds <- dds[-na.omit(bl2),]
      }
      #################################################
      
      nm <- assays(dds)[["avgTxLength"]]
      #sf <- estimateSizeFactorsForMatrix(counts(dds))
      dim(counts(dds))
      #blackliist
      #rownames(counts(dds))
      
      #dds
      keep <- rowSums(counts(dds)) > 10
      #keep <- rowMeans(counts(dds)) > 10
      dds <- dds[keep,]
      dds$protein <- relevel(dds$protein, ref = "input")
      dds <- DESeq(dds)
      #plotDispEsts(dds)
      #title(paste0(condition,"  ",htime))
      
      res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
      ####################################
      a=rownames(res)
      a=gsub(pattern = "\\..*",replacement="",a)
      m=match(KRABS,a)
      #m
      #res[na.omit(m),]
      KRABres=res[na.omit(m),c(2,6)]
      
      
      a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
      dim(counts[na.omit(mcounts),])
      colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
      ##################################  
      summary(res)
      #plotMA(res, ylim=c(-10,10))
      res=res[order(res$padj),]
      ensemblsIDS=row.names(res)
      #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
      res$symbols=getgenename(ensemblsIDS)
      res 
      # counts=counts[order(counts[,1]),]
      #df=sort(counts[,1],decreasing = T)
      #getgenename(names(df)[1:20])
      #names(df)[1:20]
      # colSums(counts)
      # tail(counts)
      counts <- counts(dds, normalized = FALSE)
      counts=data.frame(counts)
      #df[with(df, order(1)), ]
      #
      #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
      #mitcounts$names=getgenename(rownames(mitcounts))
      #mitcounts
      counts <- counts(dds, normalized = TRUE)
      counts=data.frame(counts)
      #KRABS2##########################
      counts <- counts(dds, normalized = T)
      
      acounts=rownames(counts)
      acounts=gsub(pattern = "\\..*",replacement="",acounts)
      mcounts=match(KRABS,acounts)
      #counts[na.omit(mcounts),]
      row.names(a)=a$Row.names
      a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
      a$Row.names=NULL
      dim(counts[na.omit(mcounts),])
      row.names(a)=a$Row.names
      a$Row.names=getgenename(a$Row.names)
      #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
      #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
      ################################
      mit_genes
      mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
      mitcounts=data.frame(counts[mc,])
      mitcounts$names=getgenename(rownames(mitcounts))
      colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
      print(mitcounts)
      
      
      mitres=res[mit_genes$Gene.stable.ID,]
      mitres=na.omit(mitres)
      mitres$log2FoldChange=round(mitres$log2FoldChange,2)
      
      #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
      mitres
      #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.Z22_R1","IAV.Z22_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
      
      write.csv(mitres,paste0(fdir,"/",htime,"mittable005.csv"),row.names = F)
      #############
      library(gprofiler2)
      lg=rownames(res[which(res$padj<0.05),])
      lg=gsub(lg, pattern = "\\..*", replacement="")
      
      # counts
      # summary(counts)
      # counts[c(which(counts[,3]==max(counts[,3])),which(counts[,3]==max(counts[,3]))+1),]
      # row.names(counts[which(counts[,3]==max(counts[,3])),])
      #geneCounts<-plotCounts(dds_lrt_time, gene=rownames(res)[1], intgroup = c("protein","time"),returnData = T)
      
      #dim(res)
      #res
      # counts=round(counts[rownames(res),],2)
      # colnames(counts)=gsub(pattern = "HSV.|IAV.|1_",replacement = "",colnames(counts))
      # colnames(counts)=gsub(pattern = "mock",replacement = "",colnames(counts))
      #res=cbind(res[,c(2,5,6,7)],counts[rownames(res),])
      
      write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,"/",htime,"table005.csv"))
      write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
    }
  }
}
#lf=list.files("E:/DATA/mm10/input/mm39/newdata/", recursive = T, full.names = T)


fdir="E:/DATA/mm10/reseq/results/inputvsZ22/"
lf=list.files("E:/DATA/mm10/reseq/", recursive = T, full.names = T,".sf")
lf
WriteInputvsZ22ForEveryGroup(fdir=fdir, lf=lf, usebl = F)
###########################################################

WriteInputvsFLAGForEveryGroup(fdir=fdir, lf=lf, usebl = F)

WriteInputvsFLAGForEveryGroup<-function(fdir, lf, usebl){
  htimes=c("untreat","8h")#,"12h")
  conditions=c("HSV")
  
  condition=conditions[1]
  htime=htimes[1]
  
  for (condition in conditions){
    for (htime in htimes){
      print(htime)
      print(condition)
      
      
      #print(lf)
      #filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|FLAG).*quant.sf"),lf)]
      #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|FLAG).*quant.genes"),lf)]
      lf
      filenames=lf[grep(pattern = paste0(".*",htime,".*low.*(?:_input_|FLAG).*quant.genes.sf"),lf)]
      filenames
      colData=makedata(filenames)
      colData
      ss=which(colData$virus==condition & colData$time==htime& colData$protein!="FLAG")
      sscolData=colData
      txi = tximport(filenames, txIn = FALSE, type="salmon", geneIdCol = "Name")
      dim(txi$abundance)
      #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
      ss=grep(paste0(condition,".",htime,".*(?:IgG|FLAG)"),rownames(colData))
      
      dds <- DESeqDataSetFromTximport(txi,
                                      colData = sscolData, #design = ~virus+protein)
                                      design = ~protein)
      #KRABS1######################################################
      
      KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
      KRABS=unique(KRABS[,1])
      
      #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
      #KRABS[,2]
      #KRABS=KRABS[,2]
      KRABS=gsub(pattern = "\\s",replacement="",KRABS)
      
      counts <- counts(dds, normalized = FALSE)
      acounts=rownames(counts)
      acounts=gsub(pattern = "\\..*",replacement="",acounts)
      mcounts=match(KRABS,acounts)
      #counts[na.omit(m),]
      
      
      #################################
      if (usebl){
        
        blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
        blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
        blackliist=unique(blackliist)
        
        bl2=match(blackliist,rownames(counts(dds)))
        bl2
        bl2=na.omit(bl2)
        dds <- dds[-na.omit(bl2),]
      }
      #################################################
      
      nm <- assays(dds)[["avgTxLength"]]
      #sf <- estimateSizeFactorsForMatrix(counts(dds))
      dim(counts(dds))
      #blackliist
      #rownames(counts(dds))
      
      #dds
      keep <- rowSums(counts(dds)) > 10
      #keep <- rowMeans(counts(dds)) > 10
      dds <- dds[keep,]
      dds$protein <- relevel(dds$protein, ref = "input")
      dds <- DESeq(dds)
      #plotDispEsts(dds)
      #title(paste0(condition,"  ",htime))
      
      res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
      ####################################
      a=rownames(res)
      a=gsub(pattern = "\\..*",replacement="",a)
      m=match(KRABS,a)
      #m
      #res[na.omit(m),]
      KRABres=res[na.omit(m),c(2,6)]
      
      
      a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
      dim(counts[na.omit(mcounts),])
      colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
      ##################################  
      summary(res)
      #plotMA(res, ylim=c(-10,10))
      res=res[order(res$padj),]
      ensemblsIDS=row.names(res)
      #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
      res$symbols=getgenename(ensemblsIDS)
      res 
      # counts=counts[order(counts[,1]),]
      #df=sort(counts[,1],decreasing = T)
      #getgenename(names(df)[1:20])
      #names(df)[1:20]
      # colSums(counts)
      # tail(counts)
      counts <- counts(dds, normalized = FALSE)
      counts=data.frame(counts)
      #df[with(df, order(1)), ]
      #
      #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
      #mitcounts$names=getgenename(rownames(mitcounts))
      #mitcounts
      counts <- counts(dds, normalized = TRUE)
      counts=data.frame(counts)
      #KRABS2##########################
      counts <- counts(dds, normalized = T)
      
      acounts=rownames(counts)
      acounts=gsub(pattern = "\\..*",replacement="",acounts)
      mcounts=match(KRABS,acounts)
      #counts[na.omit(mcounts),]
      row.names(a)=a$Row.names
      a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
      a$Row.names=NULL
      dim(counts[na.omit(mcounts),])
      row.names(a)=a$Row.names
      a$Row.names=getgenename(a$Row.names)
      #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
      #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
      ################################
      mit_genes
      mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
      mitcounts=data.frame(counts[mc,])
      mitcounts$names=getgenename(rownames(mitcounts))
      colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
      print(mitcounts)
      
      
      mitres=res[mit_genes$Gene.stable.ID,]
      mitres=na.omit(mitres)
      mitres$log2FoldChange=round(mitres$log2FoldChange,2)
      
      #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
      mitres
      #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.FLAG_R1","IAV.FLAG_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
      
      write.csv(mitres,paste0(fdir,"/",htime,"mittable005.csv"),row.names = F)
      #############
      library(gprofiler2)
      lg=rownames(res[which(res$padj<0.05),])
      lg=gsub(lg, pattern = "\\..*", replacement="")
      
      # counts
      # summary(counts)
      # counts[c(which(counts[,3]==max(counts[,3])),which(counts[,3]==max(counts[,3]))+1),]
      # row.names(counts[which(counts[,3]==max(counts[,3])),])
      #geneCounts<-plotCounts(dds_lrt_time, gene=rownames(res)[1], intgroup = c("protein","time"),returnData = T)
      
      #dim(res)
      #res
      # counts=round(counts[rownames(res),],2)
      # colnames(counts)=gsub(pattern = "HSV.|IAV.|1_",replacement = "",colnames(counts))
      # colnames(counts)=gsub(pattern = "mock",replacement = "",colnames(counts))
      #res=cbind(res[,c(2,5,6,7)],counts[rownames(res),])
      
      write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,"/",htime,"table005.csv"))
      write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
    }
  }
}
lf=list.files("E:/DATA/mm10/reseq/HSV/", recursive = T, full.names = T)
fdir="E:/DATA/mm10/reseq/results/inputvsFLAG/"
fdir
WriteInputvsFLAGForEveryGroup(fdir=fdir, lf=lf, usebl = F)
#0vs12####################################################
Write0vs12<-function(fdir, lf, usebl){
  htimes=c("mock","8h","12h")
  conditions=c("IAV")
  
  condition=conditions[1]
  htime=htimes[3]
  
  
  for (condition in conditions){
    #print(lf)
    filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|Z22).*quant.sf"),lf)]
    #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|Z22).*quant.genes"),lf)]
    filenames=lf[grep(pattern = paste0(".*(?:mock|12h).*Z22.*quant.genes"),lf)]
    filenames
    colData=makedata(filenames)
    colData
    ss=which(colData$virus==condition & colData$time==htime& colData$protein!="Z22")
    sscolData=colData
    txi = tximport(filenames, txIn = FALSE, type="salmon", geneIdCol = "Name")
    #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
    ss=grep(paste0(condition,".",htime,".*(?:IgG|Z22)"),rownames(colData))
    
    
    
    
    dds <- DESeqDataSetFromTximport(txi,
                                    colData = sscolData, #design = ~virus+protein)
                                    design = ~time)
    #KRABS1######################################################
    
    KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
    KRABS=unique(KRABS[,1])
    
    #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
    #KRABS[,2]
    #KRABS=KRABS[,2]
    KRABS=gsub(pattern = "\\s",replacement="",KRABS)
    
    counts <- counts(dds, normalized = FALSE)
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(m),]
    
    
    #################################
    if (usebl){
      
      blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
      blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
      blackliist=unique(blackliist)
      
      bl2=match(blackliist,rownames(counts(dds)))
      bl2=na.omit(bl2)
      dds <- dds[-na.omit(bl2),]
    }
    #################################################
    
    nm <- assays(dds)[["avgTxLength"]]
    #sf <- estimateSizeFactorsForMatrix(counts(dds))
    dim(counts(dds))
    #blackliist
    #rownames(counts(dds))
    
    #dds
    keep <- rowSums(counts(dds)) > 10
    #keep <- rowMeans(counts(dds)) > 10
    dds <- dds[keep,]
    dds$time <- relevel(dds$time, ref = "0h")
    dds <- DESeq(dds)
    #plotDispEsts(dds)
    #title(paste0(condition,"  ",htime))
    
    res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
    ####################################
    a=rownames(res)
    a=gsub(pattern = "\\..*",replacement="",a)
    m=match(KRABS,a)
    #m
    #res[na.omit(m),]
    KRABres=res[na.omit(m),c(2,6)]
    
    
    a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
    dim(counts[na.omit(mcounts),])
    colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
    ##################################  
    summary(res)
    #plotMA(res, ylim=c(-10,10))
    res=res[order(res$padj),]
    ensemblsIDS=row.names(res)
    #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
    res$symbols=getgenename(ensemblsIDS)
    res 
    # counts=counts[order(counts[,1]),]
    #df=sort(counts[,1],decreasing = T)
    #getgenename(names(df)[1:20])
    #names(df)[1:20]
    # colSums(counts)
    # tail(counts)
    counts <- counts(dds, normalized = FALSE)
    counts=data.frame(counts)
    #df[with(df, order(1)), ]
    #
    #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
    #mitcounts$names=getgenename(rownames(mitcounts))
    #mitcounts
    counts <- counts(dds, normalized = TRUE)
    counts=data.frame(counts)
    #KRABS2##########################
    counts <- counts(dds, normalized = T)
    
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(mcounts),]
    row.names(a)=a$Row.names
    a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
    a$Row.names=NULL
    dim(counts[na.omit(mcounts),])
    row.names(a)=a$Row.names
    a$Row.names=getgenename(a$Row.names)
    #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
    #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
    ################################
    mit_genes
    mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
    mitcounts=data.frame(counts[mc,])
    mitcounts$names=getgenename(rownames(mitcounts))
    colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
    print(mitcounts)
    
    
    mitres=res[mit_genes$Gene.stable.ID,]
    mitres=na.omit(mitres)
    mitres$log2FoldChange=round(mitres$log2FoldChange,2)
    
    #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
    mitres
    #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.Z22_R1","IAV.Z22_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
    
    write.csv(mitres,paste0(fdir,condition,"/",htime,"mittable005.csv"),row.names = F)
    #############
    library(gprofiler2)
    lg=rownames(res[which(res$padj<0.05),])
    lg=gsub(lg, pattern = "\\..*", replacement="")
    
    # counts
    # summary(counts)
    # counts[c(which(counts[,3]==max(counts[,3])),which(counts[,3]==max(counts[,3]))+1),]
    # row.names(counts[which(counts[,3]==max(counts[,3])),])
    #geneCounts<-plotCounts(dds_lrt_time, gene=rownames(res)[1], intgroup = c("protein","time"),returnData = T)
    
    #dim(res)
    #res
    # counts=round(counts[rownames(res),],2)
    # colnames(counts)=gsub(pattern = "HSV.|IAV.|1_",replacement = "",colnames(counts))
    # colnames(counts)=gsub(pattern = "mock",replacement = "",colnames(counts))
    #res=cbind(res[,c(2,5,6,7)],counts[rownames(res),])
    
    write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,condition,"/",htime,"table005.csv"))
    write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,condition,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
  }
  
}

Write0vs12(fdir=fdir, lf=lf, usebl = F)
#################################################
Write0vs8<-function(fdir, lf, usebl){
  htimes=c("untreat","8h")
  conditions=c("HSV")
  
  condition=conditions[1]
  htime=htimes[1]
  
  
  for (condition in conditions){
    #print(lf)
    filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|Z22).*quant.sf"),lf)]
    #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|Z22).*quant.genes"),lf)]
    filenames=lf[grep(pattern = paste0(".*(?:mock|8h).*Z22.*quant.genes"),lf)]
    filenames
    ss=which(colData$virus==condition & colData$time==htime& colData$protein!="Z22")
    sscolData=colData
    fnames1=lf[grep(pattern = paste0(".*8h.*low.*Z22.*quant.genes"),lf)]
    fnames1
    fnames2=lf[grep(pattern = paste0(".*untreat.*low.*Z22.*quant.genes"),lf)]
    fnames2
    filenames=c(fnames1,fnames2)
    filenames
    colData=makedata(c(fnames1,fnames2))
    colData
    #match(rownames(txi2$abundance),rownames(txi1$abundance))
    txi1 = tximport(fnames1, txIn = FALSE, type="salmon", geneIdCol = "Name")
    txi2 = tximport(fnames2, txIn = FALSE, type="salmon", geneIdCol = "Name")
    rn=(rownames(txi1$abundance) %in% rownames(txi2$abundance))
    txi1$abundance=txi1$abundance[rn,]
    txi1$length=txi1$length[rn,]
    txi1$counts=txi1$counts[rn,]
    txi1$abundance=cbind(txi1$abundance,txi2$abundance)
    txi1$counts=cbind(txi1$counts,txi2$counts)
    txi1$length=cbind(txi1$length,txi2$length)
    #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
    ss=grep(paste0(condition,".",htime,".*(?:IgG|Z22)"),rownames(colData))
    
    dds <- DESeqDataSetFromTximport(txi1,
                                    colData = colData, #design = ~virus+protein)
                                    design = ~time)
    #KRABS1######################################################
    
    KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
    KRABS=unique(KRABS[,1])
    
    #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
    #KRABS[,2]
    #KRABS=KRABS[,2]
    KRABS=gsub(pattern = "\\s",replacement="",KRABS)
    
    counts <- counts(dds, normalized = FALSE)
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(m),]
    
    
    #################################
    if (usebl){
      
      blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
      blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
      blackliist=unique(blackliist)
      
      bl2=match(blackliist,rownames(counts(dds)))
      bl2=na.omit(bl2)
      dds <- dds[-na.omit(bl2),]
    }
    #################################################
    
    nm <- assays(dds)[["avgTxLength"]]
    #sf <- estimateSizeFactorsForMatrix(counts(dds))
    dim(counts(dds))
    #blackliist
    #rownames(counts(dds))
    
    #dds
    keep <- rowSums(counts(dds)) > 10
    #keep <- rowMeans(counts(dds)) > 10
    dds <- dds[keep,]
    dds$time <- relevel(dds$time, ref = "0h")
    dds <- DESeq(dds)
    #plotDispEsts(dds)
    #title(paste0(condition,"  ",htime))
    res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
    dim(res)
    #res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greaterAbs")
    #dim(res)
    
    ####################################
    a=rownames(res)
    a=gsub(pattern = "\\..*",replacement="",a)
    m=match(KRABS,a)
    #m
    #res[na.omit(m),]
    KRABres=res[na.omit(m),c(2,6)]
    
    
    a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
    dim(counts[na.omit(mcounts),])
    colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
    ##################################  
    summary(res)
    #plotMA(res, ylim=c(-10,10))
    res=res[order(res$padj),]
    ensemblsIDS=row.names(res)
    #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
    res$symbols=getgenename(ensemblsIDS)
    res 
    # counts=counts[order(counts[,1]),]
    #df=sort(counts[,1],decreasing = T)
    #getgenename(names(df)[1:20])
    #names(df)[1:20]
    # colSums(counts)
    # tail(counts)
    counts <- counts(dds, normalized = FALSE)
    counts=data.frame(counts)
    #df[with(df, order(1)), ]
    #
    #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
    #mitcounts$names=getgenename(rownames(mitcounts))
    #mitcounts
    counts <- counts(dds, normalized = TRUE)
    counts=data.frame(counts)
    #KRABS2##########################
    counts <- counts(dds, normalized = T)
    
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(mcounts),]
    row.names(a)=a$Row.names
    a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
    a$Row.names=NULL
    dim(counts[na.omit(mcounts),])
    row.names(a)=a$Row.names
    a$Row.names=getgenename(a$Row.names)
    #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
    #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
    ################################
    mit_genes
    mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
    mitcounts=data.frame(counts[mc,])
    mitcounts$names=getgenename(rownames(mitcounts))
    colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
    print(mitcounts)
    
    
    mitres=res[mit_genes$Gene.stable.ID,]
    mitres=na.omit(mitres)
    mitres$log2FoldChange=round(mitres$log2FoldChange,2)
    
    #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
    mitres
    #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.Z22_R1","IAV.Z22_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
    
    write.csv(mitres,paste0(fdir,"/",htime,"mittable005.csv"),row.names = F)
    #############
    library(gprofiler2)
    lg=rownames(res[which(res$padj<0.05),])
    lg=gsub(lg, pattern = "\\..*", replacement="")
    
    # counts
    # summary(counts)
    # counts[c(which(counts[,3]==max(counts[,3])),which(counts[,3]==max(counts[,3]))+1),]
    # row.names(counts[which(counts[,3]==max(counts[,3])),])
    #geneCounts<-plotCounts(dds_lrt_time, gene=rownames(res)[1], intgroup = c("protein","time"),returnData = T)
    
    #dim(res)
    #res
    # counts=round(counts[rownames(res),],2)
    # colnames(counts)=gsub(pattern = "HSV.|IAV.|1_",replacement = "",colnames(counts))
    # colnames(counts)=gsub(pattern = "mock",replacement = "",colnames(counts))
    #res=cbind(res[,c(2,5,6,7)],counts[rownames(res),])
    
    write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,"/",htime,"table005.csv"))
    write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
  }
  
}
lf=list.files("E:/DATA/mm10/reseq/", recursive = T, full.names = T)
fdir="E:/DATA/mm10/reseq/results/Z22_0hvsxh/"
Write0vs8(fdir=fdir, lf=lf, usebl = F)
#flag##################################################
WriteFLAG0vs8<-function(fdir, lf, usebl){
  htimes=c("untreat","8h")
  conditions=c("HSV")
  
  condition=conditions[1]
  htime=htimes[1]
  
  
  for (condition in conditions){
    #print(lf)
    filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|FLAG).*quant.sf"),lf)]
    #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|FLAG).*quant.genes"),lf)]
    filenames=lf[grep(pattern = paste0(".*(?:mock|8h).*FLAG.*quant.genes"),lf)]
    filenames
    ss=which(colData$virus==condition & colData$time==htime& colData$protein!="FLAG")
    sscolData=colData
    fnames1=lf[grep(pattern = paste0(".*8h.*low.*FLAG.*quant.genes"),lf)]
    fnames1
    fnames2=lf[grep(pattern = paste0(".*untreat.*low.*FLAG.*quant.genes"),lf)]
    fnames2
    filenames=c(fnames1,fnames2)
    filenames
    colData=makedata(c(fnames1,fnames2))
    colData
    #match(rownames(txi2$abundance),rownames(txi1$abundance))
    txi1 = tximport(fnames1, txIn = FALSE, type="salmon", geneIdCol = "Name")
    txi2 = tximport(fnames2, txIn = FALSE, type="salmon", geneIdCol = "Name")
    rn=(rownames(txi1$abundance) %in% rownames(txi2$abundance))
    txi1$abundance=txi1$abundance[rn,]
    txi1$length=txi1$length[rn,]
    txi1$counts=txi1$counts[rn,]
    txi1$abundance=cbind(txi1$abundance,txi2$abundance)
    txi1$counts=cbind(txi1$counts,txi2$counts)
    txi1$length=cbind(txi1$length,txi2$length)
    #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
    ss=grep(paste0(condition,".",htime,".*(?:IgG|Z22)"),rownames(colData))
    
    dds <- DESeqDataSetFromTximport(txi1,
                                    colData = colData, #design = ~virus+protein)
                                    design = ~time)
    #KRABS1######################################################
    
    KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
    KRABS=unique(KRABS[,1])
    
    #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
    #KRABS[,2]
    #KRABS=KRABS[,2]
    KRABS=gsub(pattern = "\\s",replacement="",KRABS)
    
    counts <- counts(dds, normalized = FALSE)
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(m),]
    
    
    #################################
    if (usebl){
      
      blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
      blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
      blackliist=unique(blackliist)
      
      bl2=match(blackliist,rownames(counts(dds)))
      bl2=na.omit(bl2)
      dds <- dds[-na.omit(bl2),]
    }
    #################################################
    
    nm <- assays(dds)[["avgTxLength"]]
    #sf <- estimateSizeFactorsForMatrix(counts(dds))
    dim(counts(dds))
    #blackliist
    #rownames(counts(dds))
    
    #dds
    keep <- rowSums(counts(dds)) > 10
    #keep <- rowMeans(counts(dds)) > 10
    dds <- dds[keep,]
    dds$time <- relevel(dds$time, ref = "0h")
    dds <- DESeq(dds)
    #plotDispEsts(dds)
    #title(paste0(condition,"  ",htime))
    res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
    dim(res)
    #res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greaterAbs")
    #dim(res)
    
    ####################################
    a=rownames(res)
    a=gsub(pattern = "\\..*",replacement="",a)
    m=match(KRABS,a)
    #m
    #res[na.omit(m),]
    KRABres=res[na.omit(m),c(2,6)]
    
    
    a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
    dim(counts[na.omit(mcounts),])
    colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
    ##################################  
    summary(res)
    #plotMA(res, ylim=c(-10,10))
    res=res[order(res$padj),]
    ensemblsIDS=row.names(res)
    #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
    res$symbols=getgenename(ensemblsIDS)
    res 
    # counts=counts[order(counts[,1]),]
    #df=sort(counts[,1],decreasing = T)
    #getgenename(names(df)[1:20])
    #names(df)[1:20]
    # colSums(counts)
    # tail(counts)
    counts <- counts(dds, normalized = FALSE)
    counts=data.frame(counts)
    #df[with(df, order(1)), ]
    #
    #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
    #mitcounts$names=getgenename(rownames(mitcounts))
    #mitcounts
    counts <- counts(dds, normalized = TRUE)
    counts=data.frame(counts)
    #KRABS2##########################
    counts <- counts(dds, normalized = T)
    
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(mcounts),]
    row.names(a)=a$Row.names
    a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
    a$Row.names=NULL
    dim(counts[na.omit(mcounts),])
    row.names(a)=a$Row.names
    a$Row.names=getgenename(a$Row.names)
    #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
    #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
    ################################
    mit_genes
    mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
    mitcounts=data.frame(counts[mc,])
    mitcounts$names=getgenename(rownames(mitcounts))
    colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
    print(mitcounts)
    
    
    mitres=res[mit_genes$Gene.stable.ID,]
    mitres=na.omit(mitres)
    mitres$log2FoldChange=round(mitres$log2FoldChange,2)
    
    #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
    mitres
    #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.Z22_R1","IAV.Z22_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
    
    write.csv(mitres,paste0(fdir,"/",htime,"mittable005.csv"),row.names = F)
    #############
    library(gprofiler2)
    lg=rownames(res[which(res$padj<0.05),])
    lg=gsub(lg, pattern = "\\..*", replacement="")
    
    # counts
    # summary(counts)
    # counts[c(which(counts[,3]==max(counts[,3])),which(counts[,3]==max(counts[,3]))+1),]
    # row.names(counts[which(counts[,3]==max(counts[,3])),])
    #geneCounts<-plotCounts(dds_lrt_time, gene=rownames(res)[1], intgroup = c("protein","time"),returnData = T)
    
    #dim(res)
    #res
    # counts=round(counts[rownames(res),],2)
    # colnames(counts)=gsub(pattern = "HSV.|IAV.|1_",replacement = "",colnames(counts))
    # colnames(counts)=gsub(pattern = "mock",replacement = "",colnames(counts))
    #res=cbind(res[,c(2,5,6,7)],counts[rownames(res),])
    
    write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,"/",htime,"table005.csv"))
    write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
  }
  
}
lf=list.files("E:/DATA/mm10/reseq/", recursive = T, full.names = T)
fdir="E:/DATA/mm10/reseq/results/FLAG_0hvsxh/"
WriteFLAG0vs8(fdir=fdir, lf=lf, usebl = F)
#input 0vs12#############################
####
WriteInput0vs12<-function(fdir, lf, usebl){
  htimes=c("mock","8h","12h")
  conditions=c("IAV")
  
  condition=conditions[1]
  htime=htimes[3]
  
  
  for (condition in conditions){
    #print(lf)
    filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|Z22).*quant.sf"),lf)]
    #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|Z22).*quant.genes"),lf)]
    filenames=lf[grep(pattern = paste0(".*(?:mock|12h).*_input_.*quant.genes"),lf)]
    filenames
    colData=makedata(filenames)
    colData
    ss=which(colData$virus==condition & colData$time==htime& colData$protein!="Z22")
    sscolData=colData
    txi = tximport(filenames, txIn = FALSE, type="salmon", geneIdCol = "Name")
    #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
    ss=grep(paste0(condition,".",htime,".*(?:IgG|Z22)"),rownames(colData))
    
    
    
    
    dds <- DESeqDataSetFromTximport(txi,
                                    colData = sscolData, #design = ~virus+protein)
                                    design = ~time)
    #KRABS1######################################################
    
    KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
    KRABS=unique(KRABS[,1])
    
    #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
    #KRABS[,2]
    #KRABS=KRABS[,2]
    KRABS=gsub(pattern = "\\s",replacement="",KRABS)
    
    counts <- counts(dds, normalized = FALSE)
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(m),]
    
    
    #################################
    if (usebl){
      
      blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
      blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
      blackliist=unique(blackliist)
      
      bl2=match(blackliist,rownames(counts(dds)))
      bl2=na.omit(bl2)
      dds <- dds[-na.omit(bl2),]
    }
    #################################################
    
    nm <- assays(dds)[["avgTxLength"]]
    #sf <- estimateSizeFactorsForMatrix(counts(dds))
    dim(counts(dds))
    #blackliist
    #rownames(counts(dds))
    
    #dds
    keep <- rowSums(counts(dds)) > 10
    #keep <- rowMeans(counts(dds)) > 10
    dds <- dds[keep,]
    dds$time <- relevel(dds$time, ref = "0h")
    dds <- DESeq(dds)
    #plotDispEsts(dds)
    #title(paste0(condition,"  ",htime))
    
    res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
    ####################################
    a=rownames(res)
    a=gsub(pattern = "\\..*",replacement="",a)
    m=match(KRABS,a)
    #m
    #res[na.omit(m),]
    KRABres=res[na.omit(m),c(2,6)]
    
    
    a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
    dim(counts[na.omit(mcounts),])
    colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
    ##################################  
    summary(res)
    #plotMA(res, ylim=c(-10,10))
    res=res[order(res$padj),]
    ensemblsIDS=row.names(res)
    #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
    res$symbols=getgenename(ensemblsIDS)
    res 
    # counts=counts[order(counts[,1]),]
    #df=sort(counts[,1],decreasing = T)
    #getgenename(names(df)[1:20])
    #names(df)[1:20]
    # colSums(counts)
    # tail(counts)
    counts <- counts(dds, normalized = FALSE)
    counts=data.frame(counts)
    #df[with(df, order(1)), ]
    #
    #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
    #mitcounts$names=getgenename(rownames(mitcounts))
    #mitcounts
    counts <- counts(dds, normalized = TRUE)
    counts=data.frame(counts)
    #KRABS2##########################
    counts <- counts(dds, normalized = T)
    
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(mcounts),]
    row.names(a)=a$Row.names
    a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
    a$Row.names=NULL
    dim(counts[na.omit(mcounts),])
    row.names(a)=a$Row.names
    a$Row.names=getgenename(a$Row.names)
    #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
    #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
    ################################
    mit_genes
    mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
    mitcounts=data.frame(counts[mc,])
    mitcounts$names=getgenename(rownames(mitcounts))
    colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
    print(mitcounts)
    
    
    mitres=res[mit_genes$Gene.stable.ID,]
    mitres=na.omit(mitres)
    mitres$log2FoldChange=round(mitres$log2FoldChange,2)
    
    #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
    mitres
    #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.Z22_R1","IAV.Z22_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
    
    write.csv(mitres,paste0(fdir,condition,"/",htime,"mittable005.csv"),row.names = F)
    #############
    library(gprofiler2)
    lg=rownames(res[which(res$padj<0.05),])
    lg=gsub(lg, pattern = "\\..*", replacement="")
    
    
    write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,condition,"/",htime,"table005.csv"))
    write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,condition,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
  }
  
}
lf=list.files("E:/DATA/mm10/input/mm39/salmon/", recursive = T, full.names = T)
fdir="E:/DATA/mm10/input/mm39/newdata/IAV/input_hvsxh/"
WriteInput0vs12(fdir=fdir, lf=lf, usebl = F)
#################################################
WriteInput0vs8<-function(fdir, lf, usebl){
  htimes=c("mock","8h","8h")
  conditions=c("IAV")
  
  condition=conditions[1]
  htime=htimes[3]
  
  
  for (condition in conditions){
    #print(lf)
    filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|Z22).*quant.sf"),lf)]
    #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|Z22).*quant.genes"),lf)]
    filenames=lf[grep(pattern = paste0(".*low.*_input_.*quant.genes"),lf)]
    filenames
    
    fnames1=lf[grep(pattern = paste0(".*8h.*low.*_input_.*quant.genes"),lf)]
    fnames1
    fnames2=lf[grep(pattern = paste0(".*untreat.*low.*_input_.*quant.genes"),lf)]
    fnames2
    filenames=c(fnames1,fnames2)
    filenames
    colData=makedata(c(fnames1,fnames2))
    colData
    #match(rownames(txi2$abundance),rownames(txi1$abundance))
    txi1 = tximport(fnames1, txIn = FALSE, type="salmon", geneIdCol = "Name")
    txi2 = tximport(fnames2, txIn = FALSE, type="salmon", geneIdCol = "Name")
    rn=(rownames(txi1$abundance) %in% rownames(txi2$abundance))
    txi1$abundance=txi1$abundance[rn,]
    txi1$length=txi1$length[rn,]
    txi1$counts=txi1$counts[rn,]
    txi1$abundance=cbind(txi1$abundance,txi2$abundance)
    txi1$counts=cbind(txi1$counts,txi2$counts)
    txi1$length=cbind(txi1$length,txi2$length)
    colData
    dds <- DESeqDataSetFromTximport(txi1,
                                    colData = colData, #design = ~virus+protein)
                                    design = ~time)
    #KRABS1######################################################
    
    KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
    KRABS=unique(KRABS[,1])
    
    #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
    #KRABS[,2]
    #KRABS=KRABS[,2]
    KRABS=gsub(pattern = "\\s",replacement="",KRABS)
    
    counts <- counts(dds, normalized = FALSE)
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(m),]
    
    
    #################################
    if (usebl){
      
      blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
      blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
      blackliist=unique(blackliist)
      
      bl2=match(blackliist,rownames(counts(dds)))
      bl2=na.omit(bl2)
      dds <- dds[-na.omit(bl2),]
    }
    #################################################
    
    nm <- assays(dds)[["avgTxLength"]]
    #sf <- estimateSizeFactorsForMatrix(counts(dds))
    dim(counts(dds))
    #blackliist
    #rownames(counts(dds))
    
    #dds
    keep <- rowSums(counts(dds)) > 10
    #keep <- rowMeans(counts(dds)) > 10
    dds <- dds[keep,]
    dds$time <- relevel(dds$time, ref = "0h")
    dds <- DESeq(dds)
    #plotDispEsts(dds)
    #title(paste0(condition,"  ",htime))
    
    res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
    ####################################
    a=rownames(res)
    a=gsub(pattern = "\\..*",replacement="",a)
    m=match(KRABS,a)
    #m
    #res[na.omit(m),]
    KRABres=res[na.omit(m),c(2,6)]
    
    
    a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
    dim(counts[na.omit(mcounts),])
    colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
    ##################################  
    summary(res)
    #plotMA(res, ylim=c(-10,10))
    res=res[order(res$padj),]
    ensemblsIDS=row.names(res)
    #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
    res$symbols=getgenename(ensemblsIDS)
    res 
    # counts=counts[order(counts[,1]),]
    #df=sort(counts[,1],decreasing = T)
    #getgenename(names(df)[1:20])
    #names(df)[1:20]
    # colSums(counts)
    # tail(counts)
    counts <- counts(dds, normalized = FALSE)
    counts=data.frame(counts)
    #df[with(df, order(1)), ]
    #
    #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
    #mitcounts$names=getgenename(rownames(mitcounts))
    #mitcounts
    counts <- counts(dds, normalized = TRUE)
    counts=data.frame(counts)
    #KRABS2##########################
    counts <- counts(dds, normalized = T)
    
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(mcounts),]
    row.names(a)=a$Row.names
    a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
    a$Row.names=NULL
    dim(counts[na.omit(mcounts),])
    row.names(a)=a$Row.names
    a$Row.names=getgenename(a$Row.names)
    #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
    #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
    ################################
    mit_genes
    mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
    mitcounts=data.frame(counts[mc,])
    mitcounts$names=getgenename(rownames(mitcounts))
    colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
    print(mitcounts)
    
    
    mitres=res[mit_genes$Gene.stable.ID,]
    mitres=na.omit(mitres)
    mitres$log2FoldChange=round(mitres$log2FoldChange,2)
    
    #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
    mitres
    #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.Z22_R1","IAV.Z22_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
    
    write.csv(mitres,paste0(fdir,"/",htime,"mittable005.csv"),row.names = F)
    #############
    library(gprofiler2)
    lg=rownames(res[which(res$padj<0.05),])
    lg=gsub(lg, pattern = "\\..*", replacement="")
    
    write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,"/",htime,"table005.csv"))
    write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
  }
  
}
lf=list.files("E:/DATA/mm10/input/mm39/newdata/", recursive = T, full.names = T)
WriteInput0vs8(fdir=fdir, lf=lf, usebl = F)
#0vs12####################################################
WriteJ20vs12<-function(fdir, lf, usebl){
  htimes=c("mock","8h","12h")
  conditions=c("IAV")
  
  condition=conditions[1]
  htime=htimes[3]
  
  
  for (condition in conditions){
    #print(lf)
    filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|J2).*quant.sf"),lf)]
    #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|J2).*quant.genes"),lf)]
    filenames=lf[grep(pattern = paste0(".*(?:mock|12h).*J2.*quant.genes"),lf)]
    filenames
    colData=makedata(filenames)
    colData
    ss=which(colData$virus==condition & colData$time==htime& colData$protein!="J2")
    sscolData=colData
    txi = tximport(filenames, txIn = FALSE, type="salmon", geneIdCol = "Name")
    #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
    ss=grep(paste0(condition,".",htime,".*(?:IgG|J2)"),rownames(colData))
    
    
    
    
    dds <- DESeqDataSetFromTximport(txi,
                                    colData = sscolData, #design = ~virus+protein)
                                    design = ~time)
    #KRABS1######################################################
    
    KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
    KRABS=unique(KRABS[,1])
    
    #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
    #KRABS[,2]
    #KRABS=KRABS[,2]
    KRABS=gsub(pattern = "\\s",replacement="",KRABS)
    
    counts <- counts(dds, normalized = FALSE)
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(m),]
    
    
    #################################
    if (usebl){
      
      blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
      blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
      blackliist=unique(blackliist)
      
      bl2=match(blackliist,rownames(counts(dds)))
      bl2=na.omit(bl2)
      dds <- dds[-na.omit(bl2),]
    }
    #################################################
    
    nm <- assays(dds)[["avgTxLength"]]
    #sf <- estimateSizeFactorsForMatrix(counts(dds))
    dim(counts(dds))
    #blackliist
    #rownames(counts(dds))
    
    #dds
    keep <- rowSums(counts(dds)) > 10
    #keep <- rowMeans(counts(dds)) > 10
    dds <- dds[keep,]
    dds$time <- relevel(dds$time, ref = "0h")
    dds <- DESeq(dds)
    #plotDispEsts(dds)
    #title(paste0(condition,"  ",htime))
    
    res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
    ####################################
    a=rownames(res)
    a=gsub(pattern = "\\..*",replacement="",a)
    m=match(KRABS,a)
    #m
    #res[na.omit(m),]
    KRABres=res[na.omit(m),c(2,6)]
    
    
    a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
    dim(counts[na.omit(mcounts),])
    colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
    ##################################  
    summary(res)
    #plotMA(res, ylim=c(-10,10))
    res=res[order(res$padj),]
    ensemblsIDS=row.names(res)
    #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
    res$symbols=getgenename(ensemblsIDS)
    res 
    # counts=counts[order(counts[,1]),]
    #df=sort(counts[,1],decreasing = T)
    #getgenename(names(df)[1:20])
    #names(df)[1:20]
    # colSums(counts)
    # tail(counts)
    counts <- counts(dds, normalized = FALSE)
    counts=data.frame(counts)
    #df[with(df, order(1)), ]
    #
    #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
    #mitcounts$names=getgenename(rownames(mitcounts))
    #mitcounts
    counts <- counts(dds, normalized = TRUE)
    counts=data.frame(counts)
    #KRABS2##########################
    counts <- counts(dds, normalized = T)
    
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(mcounts),]
    row.names(a)=a$Row.names
    a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
    a$Row.names=NULL
    dim(counts[na.omit(mcounts),])
    row.names(a)=a$Row.names
    a$Row.names=getgenename(a$Row.names)
    #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
    #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
    ################################
    mit_genes
    mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
    mitcounts=data.frame(counts[mc,])
    mitcounts$names=getgenename(rownames(mitcounts))
    colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
    print(mitcounts)
    
    
    mitres=res[mit_genes$Gene.stable.ID,]
    mitres=na.omit(mitres)
    mitres$log2FoldChange=round(mitres$log2FoldChange,2)
    
    #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
    mitres
    #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.J2_R1","IAV.J2_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
    
    write.csv(mitres,paste0(fdir,condition,"/",htime,"mittable005.csv"),row.names = F)
    #############
    library(gprofiler2)
    lg=rownames(res[which(res$padj<0.05),])
    lg=gsub(lg, pattern = "\\..*", replacement="")
    
    # counts
    # summary(counts)
    # counts[c(which(counts[,3]==max(counts[,3])),which(counts[,3]==max(counts[,3]))+1),]
    # row.names(counts[which(counts[,3]==max(counts[,3])),])
    #geneCounts<-plotCounts(dds_lrt_time, gene=rownames(res)[1], intgroup = c("protein","time"),returnData = T)
    
    #dim(res)
    #res
    # counts=round(counts[rownames(res),],2)
    # colnames(counts)=gsub(pattern = "HSV.|IAV.|1_",replacement = "",colnames(counts))
    # colnames(counts)=gsub(pattern = "mock",replacement = "",colnames(counts))
    #res=cbind(res[,c(2,5,6,7)],counts[rownames(res),])
    
    write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,condition,"/",htime,"table005.csv"))
    write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,condition,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
  }
  
}
lf=list.files("E:/DATA/mm10/input/mm39/newdata/", recursive = T, full.names = T)
fdir="E:/DATA/mm10/input/mm39/deseq2/J2_0hvsxh/"
WriteJ20vs12(fdir=fdir, lf=lf, usebl = F)
#################################################
WriteJ20vs8<-function(fdir, lf, usebl){
  htimes=c("mock","8h","8h")
  conditions=c("IAV")
  
  condition=conditions[1]
  htime=htimes[3]
  
  
  for (condition in conditions){
    #print(lf)
    filenames=lf[grep(pattern = paste0(condition,"-",htime,".*(?:IgG|J2).*quant.sf"),lf)]
    #filenames=lf[grep(pattern = paste0(".*(?:mock|",condition,"-",htime,").*(?:IgG|J2).*quant.genes"),lf)]
    filenames=lf[grep(pattern = paste0(".*(?:untreat|8h).*J2.*quant.genes"),lf)]
    filenames
    colData=makedata(filenames)
    colData
    ss=which(colData$virus==condition & colData$time==htime& colData$protein!="J2")
    sscolData=colData
    txi = tximport(filenames, txIn = FALSE, type="salmon", geneIdCol = "Name")
    #txi = tximport(filenames, type="salmon",tx2gene = read.csv("E:/DATA/mm10/gencode.vM25.t2g.v.tsv",sep="\t")[,c(1,2)])#, ignoreTxVersion = T)
    ss=grep(paste0(condition,".",htime,".*(?:IgG|J2)"),rownames(colData))
    
    
    dds <- DESeqDataSetFromTximport(txi,
                                    colData = sscolData, #design = ~virus+protein)
                                    design = ~time)
    #KRABS1######################################################
    
    KRABS=read.csv("E:/DATA/mm10/krabs_no_ver.bed")
    KRABS=unique(KRABS[,1])
    
    #KRABS=read.csv("E:/DATA/mm10/KRABSconvert.csv")
    #KRABS[,2]
    #KRABS=KRABS[,2]
    KRABS=gsub(pattern = "\\s",replacement="",KRABS)
    
    counts <- counts(dds, normalized = FALSE)
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(m),]
    
    
    #################################
    if (usebl){
      
      blackliist=read.csv("E:/DATA/mm10/blgenes2.bed",header = F)
      blackliist=gsub(pattern=" ",replacement = "", blackliist[,1])
      blackliist=unique(blackliist)
      
      bl2=match(blackliist,rownames(counts(dds)))
      bl2=na.omit(bl2)
      dds <- dds[-na.omit(bl2),]
    }
    #################################################
    
    nm <- assays(dds)[["avgTxLength"]]
    #sf <- estimateSizeFactorsForMatrix(counts(dds))
    dim(counts(dds))
    #blackliist
    #rownames(counts(dds))
    
    #dds
    keep <- rowSums(counts(dds)) > 10
    #keep <- rowMeans(counts(dds)) > 10
    dds <- dds[keep,]
    dds$time <- relevel(dds$time, ref = "0h")
    dds <- DESeq(dds)
    #plotDispEsts(dds)
    #title(paste0(condition,"  ",htime))
    
    res=DESeq2::results(dds,alpha=0.1, lfcThreshold = 0.585, altHypothesis="greater")
    ####################################
    a=rownames(res)
    a=gsub(pattern = "\\..*",replacement="",a)
    m=match(KRABS,a)
    #m
    #res[na.omit(m),]
    KRABres=res[na.omit(m),c(2,6)]
    
    
    a=merge.data.frame(KRABres,counts[na.omit(mcounts),],by=0, all=TRUE)
    dim(counts[na.omit(mcounts),])
    colnames(a)=gsub(pattern = "Ali.*",replacement = "",colnames(a))
    ##################################  
    summary(res)
    #plotMA(res, ylim=c(-10,10))
    res=res[order(res$padj),]
    ensemblsIDS=row.names(res)
    #symbols <- mapIds(org.Mm.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")
    res$symbols=getgenename(ensemblsIDS)
    res 
    # counts=counts[order(counts[,1]),]
    #df=sort(counts[,1],decreasing = T)
    #getgenename(names(df)[1:20])
    #names(df)[1:20]
    # colSums(counts)
    # tail(counts)
    counts <- counts(dds, normalized = FALSE)
    counts=data.frame(counts)
    #df[with(df, order(1)), ]
    #
    #mitcounts=data.frame(counts[mit_genes$Gene.stable.ID,])
    #mitcounts$names=getgenename(rownames(mitcounts))
    #mitcounts
    counts <- counts(dds, normalized = TRUE)
    counts=data.frame(counts)
    #KRABS2##########################
    counts <- counts(dds, normalized = T)
    
    acounts=rownames(counts)
    acounts=gsub(pattern = "\\..*",replacement="",acounts)
    mcounts=match(KRABS,acounts)
    #counts[na.omit(mcounts),]
    row.names(a)=a$Row.names
    a=merge.data.frame(a,counts[na.omit(mcounts),],by=0, all=TRUE)
    a$Row.names=NULL
    dim(counts[na.omit(mcounts),])
    row.names(a)=a$Row.names
    a$Row.names=getgenename(a$Row.names)
    #if (!usebl){write.csv(a,paste0(fdir,condition,"/",htime,"KRABS.csv"))
    #  write.csv(a[which(a$padj<0.05),c(1,2,3)],paste0(fdir,condition,"/",htime,"KRABS005.csv"))}
    ################################
    mit_genes
    mc=na.omit(match(mit_genes$Gene.stable.ID.version,rownames(counts)))
    mitcounts=data.frame(counts[mc,])
    mitcounts$names=getgenename(rownames(mitcounts))
    colnames(mitcounts)=gsub(pattern = "Ali.*",replacement = "",colnames(mitcounts))
    print(mitcounts)
    
    
    mitres=res[mit_genes$Gene.stable.ID,]
    mitres=na.omit(mitres)
    mitres$log2FoldChange=round(mitres$log2FoldChange,2)
    
    #mitres=merge(mitcounts,data.frame(mitres),by=0, all=TRUE)
    mitres
    #colnames(mitres)=c("genes","IAV.IgG_R1","IAV.IgG_R2","IAV.J2_R1","IAV.J2_R2","names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","symbols")
    
    write.csv(mitres,paste0(fdir,condition,"/",htime,"mittable005.csv"),row.names = F)
    #############
    library(gprofiler2)
    lg=rownames(res[which(res$padj<0.05),])
    lg=gsub(lg, pattern = "\\..*", replacement="")
    
    # counts
    # summary(counts)
    # counts[c(which(counts[,3]==max(counts[,3])),which(counts[,3]==max(counts[,3]))+1),]
    # row.names(counts[which(counts[,3]==max(counts[,3])),])
    #geneCounts<-plotCounts(dds_lrt_time, gene=rownames(res)[1], intgroup = c("protein","time"),returnData = T)
    
    #dim(res)
    #res
    # counts=round(counts[rownames(res),],2)
    # colnames(counts)=gsub(pattern = "HSV.|IAV.|1_",replacement = "",colnames(counts))
    # colnames(counts)=gsub(pattern = "mock",replacement = "",colnames(counts))
    #res=cbind(res[,c(2,5,6,7)],counts[rownames(res),])
    
    write.csv(res[which(res$padj<0.05),c(2,6,7)],paste0(fdir,condition,"/",htime,"table005.csv"))
    write.csv(rownames(res[which(res$padj<0.05),]),paste0(fdir,condition,"/",htime,"genes005.txt"),row.names = F, col.names = F,quote = F)
  }
  
}
lf=list.files("E:/DATA/mm10/input/mm39/salmon/", recursive = T, full.names = T)
fdir="E:/DATA/mm10/input/mm39/deseq2/J2_0hvsxh/"
WriteJ20vs8(fdir=fdir, lf=lf, usebl = F)
################################################33

input0vs8=read.csv("E:/DATA/mm10/input/mm39/deseq2/input_0hvsxh/IAV/8htable005.csv")
input0vs12=read.csv("E:/DATA/mm10/input/mm39/deseq2/input_0hvsxh/IAV/12htable005.csv")

Z220vs8=read.csv("E:/DATA/mm10/input/mm39/deseq2/Z22_0hvsxh/IAV/8htable005.csv")
Z220vs12=read.csv("E:/DATA/mm10/input/mm39/deseq2/Z22_0hvsxh/IAV/12htable005.csv")

Z22vsIgG0h=read.csv("E:/DATA/mm10/input/mm39/deseq2/IgGvsZ22/IAV/mocktable005.csv")
Z22vsIgG12h=read.csv("E:/DATA/mm10/input/mm39/deseq2/IgGvsZ22/IAV/12htable005.csv")
Z22vsIgG8h=read.csv("E:/DATA/mm10/input/mm39/deseq2/IgGvsZ22/IAV/8htable005.csv")

Z22vsinput12h=read.csv("E:/DATA/mm10/input/mm39/deseq2/inputvsZ22/IAV/12htable005.csv")
Z22vsinput8h=read.csv("E:/DATA/mm10/input/mm39/deseq2/inputvsZ22/IAV/8htable005.csv")
Z22vsinput0h=read.csv("E:/DATA/mm10/input/mm39/deseq2/inputvsZ22/IAV/mocktable005.csv")

J20vs8=read.csv("E:/DATA/mm10/input/mm39/deseq2/J2_0hvsxh/IAV/8htable005.csv")
J20vs12=read.csv("E:/DATA/mm10/input/mm39/deseq2/J2_0hvsxh/IAV/12htable005.csv")

library("VennDiagram")
parts=get.venn.partitions(list(input0vs12$X,input0vs8$X,Z220vs8$X,Z220vs12$X), keep.elements = T, force.unique = T)
parts
library(RColorBrewer)
coul <- brewer.pal(4, "Pastel2")
venn.diagram(list(input0vs12$X,Z220vs8$X,input0vs8$X,Z220vs12$X),filename = 'E:/DATA/mm10/Lifs.png',
             category.names = c("input 0hvs12h","Z22 0hvs8h","input 0hvs8h","Z22 0hvs12h"),
             fill = coul,
             output=T)
lg=parts$..values..$`1`
getgenename(lg)
write.table(data.frame(lg,getgenename(lg)),"E:/DATA/mm10/input/mm39/deseq2/Lifs.csv",quote = F,col.names = F,row.names = F)

IgG0h=read.csv("E:/DATA/mm10/input/mm39/deseq2/IgGvsZ22/IAV/mocktable005.csv")
IgG12h=read.csv("E:/DATA/mm10/input/mm39/deseq2/IgGvsZ22/IAV/12htable005.csv")
IgG8h=read.csv("E:/DATA/mm10/input/mm39/deseq2/IgGvsZ22/IAV/8htable005.csv")

parts2=get.venn.partitions(list(input0vs12$X,input0vs8$X,Z220vs8$X,Z220vs12$X,IgG8h$X), keep.elements = T, force.unique = T)
coul3 <- brewer.pal(5, "Pastel2")
venn.diagram(list(input0vs12$X,Z220vs8$X,input0vs8$X,Z220vs12$X,IgG8h$X),filename = 'E:/DATA/mm10/LifsIgG.png',
             category.names = c("input 0hvs12h","Z22 0hvs8h","input 0hvs8h","Z22 0hvs12h","Z22vsIgG 8h"),
             fill = coul3,
             output=T)
IgG8h$X
list(input0vs12$X,input0vs8$X,Z220vs8$X,Z220vs12$X,IgG8h$X,IgG12h$X)
lg2=parts2$..values..$`1`
getgenename(lg2)

partsJ2=get.venn.partitions(list(input0vs12$X,input0vs8$X,Z220vs8$X,Z220vs12$X,J20vs12$X,J20vs8$X), keep.elements = T, force.unique = T)
lgJ2=partsJ2$..values..$`1`
write.table(data.frame(lgJ2,getgenename(lgJ2)),"E:/DATA/mm10/input/mm39/deseq2/LifsJ2.csv",quote = F,col.names = F,row.names = F)

length(lgJ2)

#################
lg
