thrshg19=c("E:/DATA/mm10/Nazar/fixed/prc_prediction_hg19_5.bed",
           "E:/DATA/mm10/Nazar/fixed/prc_prediction_hg19_4.bed",
           "E:/DATA/mm10/Nazar/fixed/prc_prediction_hg19_3.bed",
           "E:/DATA/mm10/Nazar/fixed/prc_prediction_hg19_2.bed",
           "E:/DATA/mm10/Nazar/fixed/prc_prediction_hg19_1.bed")
thrsmm9=c("E:/DATA/mm10/Nazar/fixed/prc_prediction_mm9_5.bed",
          "E:/DATA/mm10/Nazar/fixed/prc_prediction_mm9_4.bed",
          "E:/DATA/mm10/Nazar/fixed/prc_prediction_mm9_3.bed",
          "E:/DATA/mm10/Nazar/fixed/prc_prediction_mm9_2.bed",
          "E:/DATA/mm10/Nazar/fixed/prc_prediction_mm9_1.bed")
mkrandommm9<-function(deepZ){
  mm9_chr <- c(1:19,'X','Y')
  n=length(deepZ)
  mm9_chr=paste0("chr",mm9_chr)
  mm9_chr_size<-list(197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296,15902555)
  names(mm9_chr_size)=mm9_chr
  my_random_start  <- vector()
  my_random_end    <- vector()
  my_random_chr    <- vector()
  my_random_chr <- sample(x=mm9_chr,size=n,replace = T,prob=unlist(mm9_chr_size)/sum(unlist(mm9_chr_size)))
  my_max <- mm9_chr_size[my_random_chr]
  my_random_width <- sample(deepZ@ranges@width,size = n,replace = T)
  my_random_start<- runif(n, min=1, max=unlist(my_max)-my_random_width)
  return(GenomicRanges::GRanges(seqnames = my_random_chr,ranges = IRanges(start = my_random_start,width = my_random_width)))
}
mkrandom<-function(deepZ){
  n=length(deepZ)
  my_chr <- c(1:22,'X','Y')
  my_chr=paste0("chr",my_chr)
  my_chr_size <- list(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,  50818468, 156040895,57227415)
  mm9_chr <- c(1:19,'X','Y')
  mm9_chr=paste0("chr",mm9_chr)
  mm9_chr_size<-list(197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296,15902555)
  names(mm9_chr_size)=mm9_chr
  names(my_chr_size)=my_chr
  my_random_start  <- vector()
  my_random_end    <- vector()
  my_random_chr    <- vector()
  my_random_chr <- sample(x=my_chr,size=n,replace = T,prob=unlist(my_chr_size)/sum(unlist(my_chr_size)))
  my_max <- my_chr_size[my_random_chr]
  my_random_width <- sample(deepZ@ranges@width,size = n,replace = T)
  my_random_start<- runif(n, min=1, max=unlist(my_max)-my_random_width)
  return(GenomicRanges::GRanges(seqnames = my_random_chr,ranges = IRanges(start = my_random_start,width = my_random_width)))
}
########################################
df = data.frame(matrix(NA, ncol=10, nrow=8))
library(doParallel)
library(regioneR)
library(rtracklayer)
#run##################################################################
thr=5
for (thr in 5:1){
  deepZ=import(thrsmm9[6-thr],format = "bed")
  f2=list.files("E:/DATA/mm10/Nazar/histones/mm9",full.names = T)
  f=list.files("E:/DATA/mm10/Nazar/tfs/mm9",full.names = T)
  f=c(f,f2)
  f12=f[grep(ignore.case=TRUE,f,pattern = "(:?EP300|BMI1|BRD4|BRD9|CHD7|CTCFL|ERCC3|ERG|ESR1|EZH2|FOXA1|GABPA|GFP|HIF1A|KDM5B|MAX|MED12|MYC|POU5F1|RAD21|RNF2|RUNX1|SMARCA4|SMC1A|SPI1|TRIM28).bed")]
  f22=f2[grep(f2,pattern="(:?H2A.Z|H3|H3.3|H3ac|H3K18ac|H3K27me3|H3K36me3|H3K4me1|H3K4me2|H3K4me3|H3K9ac|H3K9K14ac|H3K9me3)$")]
  f=c(f12,f22)
  p_upper=vector()
  obs=vector()
  perc=vector()
  m=vector()
  s=vector()
  m_perc=vector()
  s_perc=vector()
  p_lower=vector()
  p_upper=vector()
  p_lower_hg=vector()
  p_upper_hg=vector()
  m_hg=vector()
  s_hg=vector()
  obs_hg=vector()
  perc_hg=vector()
  m_perc_hg=vector()
  s_perc_hg=vector()
  #f=list.files("E:/DATA/mm10/Nazar/tf/mm9",full.names = T)
  #f2=list.files("E:/DATA/mm10/Nazar/histones/mm9",full.names = T)
  #f=c(f,f2)
  f=f[!file.info(f)$isdir]
  f=f[grep(f,pattern=".*.bed$")]
  for (a in f[c(13,16,28,39,57,180,300,230)]){
    print(a)
    feature=import(a,format = "bed")
    #noov=sum(s==0)
    #print(c(a,(length(s)-noov)/length(s)))
    #perc=c(perc,(length(s)-noov)/length(s))
    repeats=100
    nov=vector()
    n=length(deepZ)
    cl=parallel::makeCluster(4)
    doParallel::registerDoParallel(cl)
    nov=foreach(i = 1:repeats,.packages="GenomicRanges",.combine=c) %dopar% {
      length(subsetByOverlaps(mkrandommm9(deepZ),feature))
      #length(subsetByOverlaps(randomizeRegions(deepZ, genome="mm9",non.overlapping=T),feature))
    }
    parallel::stopCluster(cl)
    length(feature)
    consZ=length(subsetByOverlaps(deepZ,feature))
    m=c(m,mean(nov))
    s=c(s,sd(nov))
    obs=c(obs,length(subsetByOverlaps(deepZ,feature)))
    p_lower=c(p_lower,(sum(nov<=consZ)+1)/(repeats+1))
    p_upper=c(p_upper,(sum(nov>=consZ)+1)/(repeats+1))
    
    #m=matrix(c(length(s)-noov,noov,length(deepZ)-noov,noov)
    #fisher.test()
  }
  perc=obs/length(deepZ)
  m_perc=m/length(deepZ)
  s_perc=s/length(deepZ)
  deepZ=import(thrshg19[6-thr],format = "bed")
  #f=list.files("E:/DATA/mm10/Nazar/tf/hg19",full.names = T)
  #f2=list.files("E:/DATA/mm10/Nazar/histones/hg19",full.names = T)
  #f=c(f,f2)
  f2=list.files("E:/DATA/mm10/Nazar/histones/hg19",full.names = T)
  f=list.files("E:/DATA/mm10/Nazar/tfs/hg19",full.names = T)
  f12=f[-grep(ignore.case=TRUE,f,pattern = "(:?EP300|BMI1|BRD4|BRD9|CHD7|CTCFL|ERCC3|ERG|ESR1|EZH2|FOXA1|GABPA|GFP|HIF1A|KDM5B|MAX|MED12|MYC|POU5F1|RAD21|RNF2|RUNX1|SMARCA4|SMC1A|SPI1|TRIM28).bed")]
  f22=f2[-grep(f2,pattern="(:?H2A.Z|H3|H3.3|H3ac|H3K18ac|H3K27me3|H3K36me3|H3K4me1|H3K4me2|H3K4me3|H3K9ac|H3K9K14ac|H3K9me3)$")]
  f=c(f12,f22)
  f=f[!file.info(f)$isdir]
  f=f[grep(f,pattern=".*.bed$")]
  f
  ####change back!!!!!##############################################################################
  for (a in f[c(13,16,28,39,57,180,300,230)]){
    print(a)
    feature=import(a,format = "bed")
    #noov=sum(s==0)
    #print(c(a,(length(s)-noov)/length(s)))
    #perc=c(perc,(length(s)-noov)/length(s))
    noov=length(subsetByOverlaps(deepZ,feature))
    obs_hg=c(obs_hg,noov)
    #print(c(a,(length(t)-noov)/length(t)))
    perc_hg=c(perc_hg,noov/length(deepZ))
    repeats=100
    nov=vector()
    n=length(deepZ)
    cl=parallel::makeCluster(6)
    doParallel::registerDoParallel(cl)
    nov=foreach(i = 1:repeats,.packages=c("GenomicRanges","regioneR"),.combine=c) %dopar% {
      length(subsetByOverlaps(mkrandom(deepZ),feature))
      #length(subsetByOverlaps(randomizeRegions(deepZ, genome="hg19",non.overlapping=T),feature))
    }
    parallel::stopCluster(cl)
    #length(feature)
    consZ=length(subsetByOverlaps(deepZ,feature))
    m_hg=c(m_hg,mean(nov))
    s_hg=c(s_hg,sd(nov))
    p_lower_hg=c(p_lower_hg,(sum(nov<=consZ)+1)/(repeats+1))
    p_upper_hg=c(p_upper_hg,(sum(nov>=consZ)+1)/(repeats+1))
    #obs_hg=c(obs_hg,length(subsetByOverlaps(deepZ,feature)))
    
  }
  perc_hg=obs_hg/length(deepZ)
  m_perc_hg=m_hg/length(deepZ)
  s_perc_hg=s_hg/length(deepZ)
  #p_lower_hg=c(p_lower_hg,p_lower_hg2)
  #p_lower_hg=c(p_upper_hg,p_upper_hg2)
  #thrn=gsub(thr,pattern = "(:?.*thr|.bed)",replacement = "")
  #thrn
  tfs=gsub(f,pattern = ".*/",replacement = "")
  df[,thr*2]=perc_hg
  df[,thr*2-1]=perc
  write.csv(data.frame(p_lower_hg,p_upper_hg,obs_hg,m_hg,s_hg,perc_hg,m_perc_hg,s_perc_hg,p_lower,p_upper,obs,m,s,perc,m_perc,s_perc), paste0("E:/DATA/mm10/Nazar/overlaps/final/random_fixed_thr",thr,".csv"),quote = F, row.names = F)
  #write.csv(data.frame(tfs,p_lower_hg,p_upper_hg,obs_hg,m_hg,s_hg,perc_hg,m_perc_hg,s_perc_hg,p_lower,p_upper,obs,m,s,perc,m_perc,s_perc), paste0("E:/DATA/mm10/Nazar/overlaps/final/random_fixed_thr",thr,".csv"),quote = F, row.names = F)
  
  #write.csv(data.frame(tfs,p_lower_hg,p_upper_hg,m_hg,s_hg,p_lower,p_upper,m,s), paste0("E:/DATA/mm10/Nazar/overlaps/pval3",thr,".csv"),quote = F, row.names = F)
}

#regioner version#####################
df = data.frame(matrix(NA, ncol=10, nrow=8))
library(doParallel)
library(regioneR)
library(rtracklayer)
#run##################################################################
thr=5
library(regioneR)
df = data.frame(matrix(NA, ncol=20, nrow=8))
for (thr in 2:1){
  print(thr)
  deepZ=import(thrsmm9[6-thr],format = "bed")
  f2=list.files("E:/DATA/mm10/Nazar/histones/mm9",full.names = T)
  f=list.files("E:/DATA/mm10/Nazar/tfs/mm9",full.names = T)
  f=c(f,f2)
  f12=f[-grep(ignore.case=TRUE,f,pattern = "(:?EP300|BMI1|BRD4|BRD9|CHD7|CTCFL|ERCC3|ERG|ESR1|EZH2|FOXA1|GABPA|GFP|HIF1A|KDM5B|MAX|MED12|MYC|POU5F1|RAD21|RNF2|RUNX1|SMARCA4|SMC1A|SPI1|TRIM28).bed")]
  f12
  f22=f2[-grep(f2,pattern="(:?H2A.Z|H3|H3.3|H3ac|H3K18ac|H3K27me3|H3K36me3|H3K4me1|H3K4me2|H3K4me3|H3K9ac|H3K9K14ac|H3K9me3)$")]
  f=c(f12,f22)
  f
  p=vector()
  alt=vector()
  f=f[!file.info(f)$isdir]
  f=f[grep(f,pattern=".*.bed$")]
  i=1
  for (a in f[c(13,16,28,39,57,180,300,230)]){
    print(a)
    feature=import(a,format = "bed")
    #noov=sum(s==0)
    #print(c(a,(length(s)-noov)/length(s)))
    #perc=c(perc,(length(s)-noov)/length(s))
    t=permTest(A=deepZ, B=feature, ntimes=100,
               randomize.function=randomizeRegions,
               evaluate.function=numOverlaps, count.once=TRUE,
               genome="mm9", mc.set.seed=FALSE, force.parallel = T)
  df[i,thr*4]=t$numOverlaps$pval
  df[i,thr*4-1]=t$numOverlaps$alternative
  i=i+1
  }
  i=1
  deepZ=import(thrshg19[6-thr],format = "bed")

  f2=list.files("E:/DATA/mm10/Nazar/histones/hg19",full.names = T)
  f=list.files("E:/DATA/mm10/Nazar/tfs/hg19",full.names = T)
  f12=f[-grep(ignore.case=TRUE,f,pattern = "(:?EP300|BMI1|BRD4|BRD9|CHD7|CTCFL|ERCC3|ERG|ESR1|EZH2|FOXA1|GABPA|GFP|HIF1A|KDM5B|MAX|MED12|MYC|POU5F1|RAD21|RNF2|RUNX1|SMARCA4|SMC1A|SPI1|TRIM28).bed")]
  f22=f2[-grep(f2,pattern="(:?H2A.Z|H3|H3.3|H3ac|H3K18ac|H3K27me3|H3K36me3|H3K4me1|H3K4me2|H3K4me3|H3K9ac|H3K9K14ac|H3K9me3)$")]
  f=c(f12,f22)
  f=f[!file.info(f)$isdir]
  f=f[grep(f,pattern=".*.bed$")]
  f
  ####change back!!!!!##############################################################################
  for (a in f[c(13,16,28,39,57,180,300,230)]){
    print(a)
    feature=import(a,format = "bed")
    t=permTest(A=deepZ, B=feature, ntimes=100,
               randomize.function=randomizeRegions,
               evaluate.function=numOverlaps, count.once=TRUE,
               genome="hg19", mc.set.seed=FALSE, force.parallel = T)
    df[i,thr*4-2]=t$numOverlaps$pval
    df[i,thr*4-3]=t$numOverlaps$alternative
    i=i+1
  }
  write.csv(df, paste0("E:/DATA/mm10/Nazar/overlaps/final/random_regioner_fixed_thr",thr,".csv"),quote = F, row.names = F)
  #write.csv(data.frame(tfs,p_lower_hg,p_upper_hg,m_hg,s_hg,p_lower,p_upper,m,s), paste0("E:/DATA/mm10/Nazar/overlaps/pval3",thr,".csv"),quote = F, row.names = F)
}
