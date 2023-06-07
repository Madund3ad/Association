#!/usr/bin/Rscript
#Rscript E:/projects/AI/main_assoc.R -i E:/projects/AI/input -g E:/DATA/mm10/Nazar/hg19.ensGene.gtf -c E:/hg19.chrom.sizes

suppressMessages(library(optparse)) # 1 7 3
option_list = list(
  make_option(c("-i", "--dir"), type="character", default=NULL, 
              help="input dir", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output dir [default= %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL, 
              help="path to genome file", metavar="character"),
  make_option(c("-c", "--chrom"), type="character", default=NULL, 
              help="chrom sizes file", metavar="character"),
  make_option(c("-v", "--verbose"), type="character", default=0, 
              help="verbosity level", metavar="character"),
  make_option(c("-p", "--simulate"), type="character", default=0, 
              help="enable simulated p-value", metavar="character"),
  make_option(c("-n", "--nsim"), type="character", default=1000, 
              help="number of simulations [default= %default]", metavar="character"),
  make_option(c("-m", "--mode"), type="character", default="region", 
              help="set scoring mode/n
              region - use number of overlapping regions \n
              bp - use overlap size \n
              weighted_bp use score column", metavar="character"),
  make_option(c("-s", "--strand"), type="logical", default=T, 
              help="ignore strands", metavar="character"),
  make_option(c("-t", "--target"), type="logical", default=T, 
              help="path to target file in one vs many mode", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print_log<-function(text,log_con){
  cat(paste(Sys.time(),text),file = log_con, sep="\n")
}

print(opt$genome)

if (is.null(opt$out)){
  output_dir=paste0(getwd(),"/out")
  print(paste0("output path is not provided writing to ", output_dir))
  dir.create(output_dir,showWarnings = FALSE)
}else(
  output_dir=opt$out
)
if (is.null(opt$genome)){
  print_help(opt_parser)
  stop("genome file must be supplied", call.=FALSE)
}

if (is.null(opt$chrom)){
  print_help(opt_parser)
  stop("chrom sizes file must be supplied", call.=FALSE)
}

#install.packages(c("GenomicFeatures","GenomicRanges","rtracklayer"."regioneR","curl"),)
suppressMessages(library(GenomicFeatures)) #1.46.5
suppressMessages(library(GenomicRanges))#1,46.1
suppressMessages(library(rtracklayer))#1.54
suppressMessages(library(regioneR)) #1.26.1
suppressMessages(library(curl)) #4.3.2
################################
one_vs_many<-function(markdir="E:/DATA/mm10/Nazar/tf/hg19/",target_ranges,verbose=T){
  log_con <- file(paste0(output_dir,"/log.txt"), open="w")
  log_con <- file(paste0(output_dir,"/log.txt"), open="a")
  i=1
  marks=list.files(markdir,full.names = T)
  marks=marks[!file.info(marks)$isdir]
  
  df=data.frame(matrix(NA, nrow = length(marks), ncol = 10))
  rownames(df)=gsub(marks,pattern = ".*/|.bed",replacement = "")
  colnames(df)=c("hits_from_mark","n_mark","hits_from_mark/n_mark",
                 "hits_from_target","n_target","hits_from_target/n_target",
                 "intersect_width","target_width","mark_width","IoU")
  for (mark in marks[2:4]){
    if (verbose){
      print(paste0(Sys.time()," Reading file ¹",i,"/",length(marks)," ", mark))
    }
    print_log(paste("Reading file ¹",i,"/",length(marks)," ", mark),log_con)
    if (!file.exists(mark)){
      print(paste0("File ",mark ,"does not exist!"))
    }
    tryCatch( { mark_ranges=import(mark,format="bed");}
              , error = function(e) {print(paste("Can't import: ",mark));break})
    mark_ranges=import(mark,format="bed")
    n_mark=length(mark_ranges)
    n_target=length(target_ranges)
    seqlevelsStyle(target_ranges)<-"UCSC"
    seqlevelsStyle(mark_ranges)<-"UCSC"
    hits=findOverlaps(target_ranges,mark_ranges,ignore.strand=T)
    hits_from_target=length(unique(hits@from))
    hits_from_mark=length(unique(hits@to))
    intersect_width=sum(width(intersect(target_ranges,mark_ranges,ignore.strand=T)))
    target_width=sum(width(target_ranges))
    mark_width=sum(width(mark_ranges))
    IoU=intersect_width/sum(width(union(target_ranges,mark_ranges,ignore.strand=T)))
    res=c(hits_from_mark,n_mark,hits_from_mark/n_mark,
          hits_from_target,n_target,hits_from_target/n_target,intersect_width,target_width,mark_width,IoU)
    df[i,]=res
    i=i+1
  }
  print(paste(Sys.time(),"Saving results in", output_dir))
  print_log(paste("Saving results in", output_dir),log_con)
  print(paste(Sys.time(),"Job finished."))
  print_log("Job finished.",log_con)
  close(log_con)
  return(df)
}
########################################
opt$genome <-"E:/DATA/mm10/Nazar/hg38.ensGene.gtf"
txdb <- makeTxDbFromGFF(opt$genome, format="gtf")
##############
opt$dir="E:/DATA/mm10/Nazar/tf/hg19/"
########################################
if (is.null(opt$target)){
  print(paste0(Sys.time(), "Starting in many vs many mode"))
  many_vs_many(markdir=opt$dir,verbose=opt$verbose,opt$mode,s_mode=opt$mode,p=opt$simulate)
  }else{
    print(paste0(Sys.time(), "Starting in one vs many mode"))
    one_vs_many(markdir=opt$dir,target_ranges=opt$target,verbose=opt$verbose,s_mode=opt$mode,p=opt$simulate)
}

