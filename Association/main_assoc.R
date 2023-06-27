#!/usr/bin/Rscript
#example Rscript E:/projects/AI/main_assoc.R -i E:/projects/AI/test -g E:/DATA/mm10/Nazar/hg19.ensGene.gtf -c E:/hg19.chrom.sizes
#source("E:/projects/AI/funcs.R")
suppressMessages(library(optparse)) # 1 7 3
option_list = list(
  make_option(c("-i", "--dir"), type="character", default=NULL, 
              help="path to input folder with bed files", metavar="<input dir>"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output files directory [default= %default]", metavar="<path>"),
  make_option(c("-g", "--genome"), type="character", default=NULL, 
              help="path to genome annotation file", metavar="<genome.gtf>"),
  make_option(c("-c", "--chrom"), type="character", default=NULL, 
              help="path to chrom sizes file", metavar="chrom.sizes"),
  make_option(c("-v", "--verbose"), type="logical", default=FALSE, 
              help="verbosity level [default= %default]", metavar="", action = "store_true"),
  make_option(c("-p", "--simulate"), type="logical", default=FALSE, 
              help="enable simulated p-value", metavar="", action = "store_true"),
  make_option(c("-n", "--nsim"), type="integer", default=1000, 
              help="number of simulations [default= %default]", metavar="numeric"),
  make_option(c("-m", "--mode"), type="character", default="region", 
              help="set scoring mode
              region - use number of overlapping regions 
              bp - use overlap size 
              weighted_bp use score column", metavar="character"),
  make_option(c("-s", "--strand"), type="logical", default=FALSE, 
              help="take strand into consideration", metavar="character", action = "store_true"),
  make_option(c("-t", "--target"), type="character", default=NULL, 
              help="path to target file in one vs many mode", metavar="character"),
  make_option(c("-k", "--mask"), type="character", default=NULL, 
              help="path to mask file", metavar="mask.bed"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$out)){
  output_dir=paste0(getwd(),"/storage")
  cat(paste0("output path is not provided writing to ", output_dir,"\n"))
  dir.create(output_dir,showWarnings = FALSE)
}else{
  if (!dir.exists(opt$out)){
  stop(paste("Output directory",opt$out,"does not exist!"), call.=FALSE)
  }
  output_dir=opt$out
}
b=file.create(paste0(output_dir,"/log.txt"),showWarnings = FALSE)
log_con <- file(paste0(output_dir,"/log.txt",open = "a"))


print_log<-function(text,log_con=log_con){
  cat(paste(Sys.time(),text),file = paste0(output_dir,"/log.txt"), sep="\n",append = T)
}
if (is.null(opt$genome)){
  print_help(opt_parser)
  print_log("genome file must be supplied")
  cat("genome file is not supplied")
}else{
  if (!file.exists(opt$genome)){
    print_log("genome file does not exist")
    stop("genome file does not exist", call.=FALSE)
  }
  
}
if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("input dir path must be supplied", call.=FALSE)
}
if (!dir.exists(opt$dir)){
  stop(paste("Input directory",opt$dir,"does not exist!"), call.=FALSE)
}


if (is.null(opt$chrom)){
  print_help(opt_parser)
  stop("chrom sizes file must be supplied", call.=FALSE)
}




cat("Using parameters:\n")
for (i in 1:length(opt)){
  cat(paste(names(opt[i]),"=", opt[i],"\n"))
  if (opt$verbose){
    print_log(paste(names(opt[i]),"=", opt[i],"\n"),log_con)}
}

suppressMessages(library(GenomicFeatures)) #1.46.5
suppressMessages(library(GenomicRanges))#1,46.1
suppressMessages(library(rtracklayer))#1.54
suppressMessages(library(regioneR)) #1.26.1
suppressMessages(library(curl)) #4.3.2
suppressMessages(library(doParallel))
####################################
mscore<-function(dat_j,dat_i,s_mode){
  if (s_mode=="region"){
    dfji=length(subsetByOverlaps(dat_j,dat_i,ignore.strand=opt$strand))
    dfji1=dfji/length(dat_j)
  }
  if (s_mode=="bp"){
    dfji=sum(width(intersect(dat_j,dat_i,ignore.strand=opt$strand)))
    dfji1=dfji/sum(width(dat_j))
  }
  if (s_mode=="weighted_bp"){
    weighted_mark=sum(width(dat_i)*dat_i$score/1000)
    hits=GenomicRanges::findOverlaps(dat_j,dat_i)
    overlaps <- pintersect(dat_i[subjectHits(hits)],dat_j[queryHits(hits)])
    dfji==sum(width(overlaps)*overlaps$score/1000)
    dfji1=dfji/weighted_mark
  }
  return(c(dfji,dfji1))
}
sim_score<-function(track1,track2,ss,opt=opt){
  s_mode=opt$mode
  no_cores <- detectCores() - 1  
  cl <- makeCluster(no_cores)  
  registerDoParallel(cl)
  genome=get_genome_ranges(chr_sizes =opt$chrom)
  if (s_mode=="region"){
    nov=foreach(i = 1:opt$nsim,.packages=c("GenomicRanges","regioneR"),.combine=c) %dopar% {
      length(subsetByOverlaps(randomizeRegions(track1, mask = opt$mask, genome=genome,non.overlapping=T),track2))
    }
    obs=length(subsetByOverlaps(track1,track2))
    ops_pr=obs/length(track1)
  }
  if (s_mode=="bp"){
    nov=foreach(i = 1:opt$nsim,.packages=c("GenomicRanges","regioneR"),.combine=c) %dopar% {
      sum(width(intersect(randomizeRegions(track1, mask = opt$mask, genome=genome,non.overlapping=T),track2)))
    }
    obs=sum(width(intersect(track1,track2)))
    ops_pr=obs/sum(width(track1))
  }
  if (s_mode=="weighted_bp"){
    weighted_mark=sum(width(track2)*track2$score/1000)
    hits=GenomicRanges::findOverlaps(track1,track2)
    overlaps <- pintersect(track1[hits@from],track2[hits@to])
    obs=sum(width(overlaps)*overlaps$score/1000)
    obs_pr=obs/weighted_mark
    
    rr=shuffle_n(track1,genome=get_genome_ranges(chr_sizes = opt$chrom),n=opt$nsim,mask = opt$mask )
    nov=vector()
    for (i in 1:opt$nsim){
      rr_s=rr[(length(track1)*(i-1)+1):(length(track1)*i)]
      hits=GenomicRanges::findOverlaps(rr_s,track2)
      overlaps <- pintersect(track2[hits@to],rr_s[hits@from])
      nov[i]=sum(width(overlaps)*overlaps$score/1000)
    }
  }
  
  m=mean(nov)
  p_lower=(sum(nov<=obs)+1)/(opt$nsim+1)
  p_upper=(sum(nov>=obs)+1)/(opt$nsim+1)
  alt=c("lower", "greater")[which.min(c(p_lower,p_upper))]
  p=min(p_lower,p_upper)
  return(c(m,p,alt))
}
many_vs_many<-function(markdir=opt$dir,verbose=opt$verbose,s_mode=opt$mode,p_calc=opt$simulate,ig_strand=!opt$strand,ss=NULL){
  #dats=lapply(markdir,FUN=read_bed3)
  fs=list.files(markdir,pattern = "*.bed")
  datnames=gsub(x=fs,pattern = ".*/|.bed",replacement = "")
  df=data.frame(matrix(nrow =length(fs) ,ncol = length(fs)*2))
  for (i in 1:length(fs)){
    dat_i=import(paste0(markdir,"/",fs[i]),format = "bed")
    if (verbose){
      cat(paste0(Sys.time(), " Running: ",fs[i],"\n"))
      print_log(text=paste0(Sys.time(), " Running: ",fs[i],"\n"),log_con)
    }
    if (!is.null(ss)){
      dat_i=intersect(ss,dat_i,ignore.strand=T)
    }
    
    for (j in i:length(fs)){
      if (verbose){
        cat(paste0(Sys.time(), " : ",fs[j],"\n"))
        print_log(text=paste0(Sys.time(), " to: ",fs[j],"\n"),log_con)
      }
      dat_j=import(paste0(markdir,"/",fs[j]))
      
      if (!is.null(ss)){
        dat_j=intersect(ss,dat_j,ignore.strand=T)
      }
      if (opt$simulate){
        sim_ij=sim_score(dat_i,dat_j,ss,opt=opt)
        sim_ji=sim_score(dat_j,dat_i,ss,opt=opt)
        p=sim_ij[2]
        r=mscore(dat_i,dat_j,s_mode=opt$mode)
        df[i,j*2-1]=r[1]
        df[i,j*2]=paste0(round(r[2]*100,2),"%"," p-val_",sim_ij[3],"=",p)
        p=sim_ji[2]
        r=mscore(dat_j,dat_i,s_mode=opt$mode)
        df[j,i*2-1]=r[1]
        df[j,i*2]=paste0(round(r[2]*100,2),"%"," p-val_",sim_ji[3],"=",p)
      }else{
      p=chisq_width(dat_i,dat_j,bglen =get_genome_len(opt$chrom) )
      r=mscore(dat_i,dat_j,s_mode=opt$mode)
      df[i,j*2-1]=r[1]
      df[i,j*2]=paste(round(r[2]*100,2),"%"," p-val=",p[[1]]$p.value)
      p=chisq_width(dat_j,dat_i,get_genome_len(opt$chrom))
      r=mscore(dat_j,dat_i,s_mode=opt$mode)
      df[j,i*2-1]=r[1]
      df[j,i*2]=paste(round(r[2]*100,2),"%"," p-val=",p[[1]]$p.value)
      }
    }
    
  }
  colnames(df)=rep(datnames,each=2)
  rownames(df)=datnames
  return(df)
}
shuffle_n<-function(x,genome,n=1000,mask=NULL){
  res=foreach(i = 1:n,.packages=c("GenomicRanges","regioneR"),.combine=c) %dopar% {
    randomizeRegions(x, mask = mask, genome=genome,non.overlapping=T)
  }
  return(res)
}
one_vs_many<-function(markdir=opt$dir,target_ranges=opt$target,verbose=opt$verbose,s_mode=opt$mode,p=opt$simulate,ss=NULL){
  
  i=1
  marks=list.files(markdir,full.names = T)
  marks=marks[!file.info(marks)$isdir]
  target_ranges=import(target_ranges,format="bed")
  if (!is.null(ss)){
    target_ranges=intersect(ss,target_ranges,ignore.strand=T)
  }
  df=data.frame(matrix(NA, nrow = length(marks), ncol = 10))
  rownames(df)=gsub(marks,pattern = ".*/|.bed",replacement = "")
  colnames(df)=c("hits_from_mark","n_mark","hits_from_mark/n_mark",
                 "hits_from_target","n_target","hits_from_target/n_target",
                 "intersect_width","target_width","mark_width","IoU")
  for (mark in marks){
    if (verbose){
      cat(paste0(Sys.time()," Reading file №",i,"/",length(marks)," ", mark,"\n"))
    }
    print_log(paste("Reading file №",i,"/",length(marks)," ", mark),log_con)
    if (!file.exists(mark)){
      cat(paste0("File ",mark ,"does not exist!","\n"))
    }
    tryCatch( { mark_ranges=import(mark,format="bed");}
              , error = function(e) {print(paste("Can't import: ",mark));break})
    mark_ranges=import(mark,format="bed")
    if (!is.null(ss)){
      mark_ranges=intersect(ss,mark_ranges,ignore.strand=T)
    }
    n_mark=length(mark_ranges)
    n_target=length(target_ranges)
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

get_genome_len<-function(chr_sizes){
  df=read.delim(chr_sizes, sep = "\t",header = F)
  return(sum(df[,2]))
}
get_genome_ranges<-function(chr_sizes=opt$chrom){
  df=read.delim(chr_sizes, sep = "\t",header = F)
  genome_ranges=GenomicRanges::GRanges(seqnames = df[,1],ranges=IRanges(start=0,end=df[,2]))
  return(genome_ranges)
}
chisq_width<-function(a,b,bglen=3088269832,simulate.p.value=F){
  a11=sum(width(intersect(a,b,ignore.strand=T)))
  a12=sum(width(a))-a11
  a21=sum(width(b))-a11
  a22=bglen-a21-a12+a11
  m=matrix(c(a11,a12,a21,a22),ncol = 2)
  return(list(chisq.test(m,simulate.p.value=simulate.p.value),m))
}
get_distr<-function(target,s){
  d=vector()
  for (i in 1:length(s)){
    d[i]=sum(width(GenomicRanges::intersect(import(target,format = "bed"),s[[i]],ignore.strand=T)))
  }
  return(d)
}
get_bg<-function(s){
  bg=lapply(s, function(x) sum(width(x)))
  bg=unlist(bg)
  return(bg/sum(bg))
}

################################
chisq_width<-function(a,b,bglen=3088269832,simulate.p.value=F){
  a11=sum(width(intersect(a,b,ignore.strand=T)))
  a12=sum(width(a))-a11
  a21=sum(width(b))-a11
  a22=bglen-a21-a12+a11
  m=matrix(c(a11,a12,a21,a22),ncol = 2)
  return(list(chisq.test(m,simulate.p.value=simulate.p.value),m))
}
chisq_int<-function(track1,track2,genome=opt$genome){
  #bedtools like
  overlapCounts=length(findOverlaps(track1,track2,ignore.strand=opt$strand))
  hits=findOverlaps(track1,track1)
  qq=length(track1)
  sq=length(track2)
  meanInt=mean(width(track1))+2+mean(width(track2))
  
  gl=sum(width(get_genome_ranges(chr_sizes = opt$chrom)))
  gl/meanInt
  a12 = max(0,  qq-overlapCounts )
  a21 = max(0, sq - overlapCounts)
  a22_full = max(a21 + a12 + a11,gl/meanInt)
  a22 = max(0, as.integer(n22_full - n12 - n21 - n11))
  
  a11=length(subsetByOverlaps(a,b,ignore.strand=T))
  a12=sum(width(a))-a11
  a21=sum(width(b))-a11
  a22=bglen-a21-a12+a11
  m=matrix(c(a11,a12,a21,a22),ncol = 2)
  return(list(chisq.test(m,simulate.p.value=F),m))
}
score_genomic_regions<-function(genome=opt$genome,pri_order=c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                              "Downstream", "Intergenic"),
                                up=1500, down=500){
  txdb <- makeTxDbFromGFF(opt$genome, format="gtf")  
  UTR5=fiveUTRsByTranscript(txdb)
  UTR5=unlist(UTR5)
  
  UTR3=threeUTRsByTranscript(txdb)
  UTR3=unlist(UTR3)  
  introns <- intronsByTranscript(txdb)
  introns <- unlist(introns)
  ex=exons(txdb)
  proms=promoters(txdb,upstream = up,downstream = down)
  downstream=flank(transcripts(txdb),300,start=F)
  genes=genes(txdb)
  
  intergenicRegions=setdiff(get_genome_ranges(chr_sizes = opt$chrom),genes,ignore.strand=TRUE)
  
  l=list()
  for (i in 1:7){
    l[[i]]=switch(pri_order[i],"Promoter"=proms,"5UTR"=UTR5, "3UTR"=UTR3, "Exon"=ex, "Intron"=introns,
                  "Downstream"=downstream, "Intergenic"=intergenicRegions)
  }
  s=list()
  s[[1]]=l[[1]]
  for (i in 2:7){
    u=l[[1]]
    for (j in 1:(i-1)){
      u=union(u,l[[j]],ignore.strand=TRUE)
    }
    s[[i]]=GenomicRanges::reduce(GenomicRanges::setdiff(l[[i]],u,ignore.strand=TRUE))
  }
  names(s)=pri_order
  for (i in 1:7){
    if (is.null(opt$target)){
      print(paste0(Sys.time()," Starting in many vs many mode for ",names(s)[i]))
      res=many_vs_many(markdir=opt$dir,verbose=opt$verbose,s_mode=opt$mode,p_calc =opt$simulate,opt$strand,ss = s[[i]])
    }else{
      print(paste0(Sys.time()," Starting in one vs many mode for ",names(s)[i]))
      res=one_vs_many(markdir=opt$dir,target_ranges=opt$target,verbose=opt$verbose,s_mode=opt$mode,p=opt$simulate,ss = s[[i]])
    }
    write.table(x = cbind(from_to = rownames(res), res), file = paste0(output_dir,"/results_",names(s)[i],".csv"),sep = "\t",row.names = F)
  }
}
###########################
if (is.null(opt$target)){
  print(paste0(Sys.time(), " Starting in many vs many mode"))
  res=many_vs_many(markdir=opt$dir,verbose=opt$verbose,s_mode=opt$mode,p=opt$simulate,opt$strand)
  }else{
    print(paste0(Sys.time(), " Starting in one vs many mode"))
  res=one_vs_many(markdir=opt$dir,target_ranges=opt$target,verbose=opt$verbose,s_mode=opt$mode,p=opt$simulate)
  }
write.table(x = cbind(file_name = rownames(res), res), file = paste0(output_dir,"/results.csv"),sep = "\t",row.names = F,quote=F)
################################
if (!is.null(opt$genome)){
  res=score_genomic_regions()
  #write.table(x = cbind(from_to = rownames(res), res), file = paste0(output_dir,"/results_genomic_distribution.csv"),sep = "\t")
}
close(log_con)

