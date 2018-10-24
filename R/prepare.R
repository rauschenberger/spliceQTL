

#' @name spliceQTL-package
#' @md
#' @aliases spliceQTL
#' 
#' @title
#' 
#' Alternative Splicing
#' 
#' @description
#' 
#' This R package includes various functions
#' for applying the global test of alternative splicing.
#' Some functions only work on the virtual machine (see below).
#' 
#' @seealso 
#' 
#' Prepare BBMRI and Geuvadis data:
#' * \code{\link{get.snps.geuvadis}} (not VM)
#' * \code{\link{get.snps.bbmri}} (only VM)
#' * \code{\link{get.exons.geuvadis}} (only VM)
#' * \code{\link{get.exons.bbmri}} (only VM)
#' 
#' Process samples and covariates:
#' * \code{\link{match.samples}}
#' * \code{\link{adjust.samples}}
#' * \code{\link{adjust.variables}}
#' 
#' Search for exons and SNPs:
#' * \code{\link{map.genes}}
#' * \code{\link{map.exons}}
#' * \code{\link{map.snps}}
#' * \code{\link{drop.trivial}}
#' 
#' Test for alternative splicing:
#' * \code{\link{test.single}}
#' * \code{\link{test.multiple}}
#'
#' @keywords documentation
#' @docType package
#' 
NULL


#' @export
#' @title
#' Get SNP data (Geuvadis)
#' 
#' @description
#' This function transforms SNP data (local machine):
#' downloads missing genotype data from ArrayExpress,
#' transforms variant call format to binary files,
#' removes SNPs with low minor allele frequency,
#' labels SNPs in the format "chromosome:position",
#' changes sample identifiers
#'  
#' @param chr
#' chromosome: integer \eqn{1-22}
#' 
#' @param data
#' local directory for VCF (variant call format)
#' and SDRF (sample and data relationship format) files;
#' 
#' @param path
#' local directory for output
#' 
#' @examples
#' path <- "C:/Users/a.rauschenbe/Desktop/spliceQTL/data"
#' 
get.snps.geuvadis <- function(chr,data=NULL,path=getwd()){
  
  if(is.null(data)){
    data <- path
    # download SNP data
    file <- paste0("GEUVADIS.chr",chr,".PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz")
    url <- paste0("http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/genotypes/",file)
    destfile <- file.path(data,file)
    if(!file.exists(destfile)){
      utils::download.file(url=url,destfile=destfile,method="auto")
    }
    # transform with PLINK
    setwd(data)
    system(paste0("plink --vcf GEUVADIS.chr",chr,".PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz",
                  " --maf 0.05 --geno 0 --make-bed --out snps",chr),invisible=FALSE)
    # obtain identifiers
    file <- "E-GEUV-1.sdrf.txt"
    url <- paste("http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/",file,sep="")
    destfile <- file.path(data,file)
    if(!file.exists(destfile)){
      utils::download.file(url=url,destfile=destfile,method="auto")
    }
  }
  
  # read into R
  bed <- file.path(data,paste("snps",chr,".bed",sep=""))
  bim <- file.path(data,paste("snps",chr,".bim",sep=""))
  fam <- file.path(data,paste("snps",chr,".fam",sep=""))
  X <- snpStats::read.plink(bed=bed,bim=bim,fam=fam)
  X$fam <- NULL; all(diff(X$map$position) > 0)
  
  # filter MAF
  maf <- snpStats::col.summary(X$genotypes)$MAF
  cond <- maf >= 0.05
  X$genotypes <- X$genotypes[,cond]
  X$map <- X$map[cond,]
  
  # format
  colnames(X$genotypes) <- paste0(X$map$chromosome,":",X$map$position)
  snps <- methods::as(object=X$genotypes,Class="numeric")
  class(snps) <- "integer"
  
  # change identifiers
  samples <- utils::read.delim(file=file.path(data,"E-GEUV-1.sdrf.txt"))
  match <- match(rownames(snps),samples$Source.Name)
  rownames(snps) <- samples$Comment.ENA_RUN.[match]
  snps <- snps[!is.na(rownames(snps)),]
  
  save(object=snps,file=file.path(path,paste0("Geuvadis.chr",chr,".RData")))
}


#' @export
#' @title
#' Get SNP data (BBMRI)
#' 
#' @description
#' This function transforms SNP data (virtual machine):
#' limits analysis to specified biobanks,
#' reads in genotype data in chunks,
#' removes SNPs with missing values (multiple biobanks/technologies),
#' removes SNPs with low minor allele frequency,
#' fuses data from multiple biobanks/technologies
#' 
#' @param chr
#' chromosome: integer \eqn{1-22}
#' 
#' @param biobank
#' character "CODAM", "LL", "LLS", "NTR", "PAN", "RS", or NULL (all)
#' 
#' @param path
#' data directory
#' 
#' @param size
#' maximum number of SNPs to read in at once;
#' trade-off between memory usage (low) and speed (high)
#' 
#' @examples
#' path <- "/virdir/Scratch/arauschenberger/trial"
#'
get.snps.bbmri <- function(chr,biobank=NULL,path=getwd(),size=500*10^3){
  
  start <- Sys.time()
  message(rep("-",times=20)," chromosome ",chr," ",rep("-",times=20))
  
  p <- 5*10^6 # (maximum number of SNPs per chromosome, before filtering)
  skip <- seq(from=0,to=p,by=size)
  if(is.null(biobank)){
    study <- c("CODAM","LL","LLS0","LLS1","NTR0","NTR1","PAN","RS")
  } else if(biobank=="LLS"){
    study <- c("LLS0","LLS1")
  } else if(biobank=="NTR"){
    study <- c("NTR0","NTR1")
  } else if(biobank %in% c("CODAM","LL","PAN","RS")){
    study <- biobank
  } else{
    stop("Invalid biobank.",call.=FALSE)
  }
  collect <- matrix(list(),nrow=length(skip),ncol=length(study))
  colnames(collect) <- study
  
  for(i in seq_along(skip)){
    message("\n","chunk ",i,": ",appendLF=FALSE)
    for(j in seq_along(study)){
      message(study[j],"  ",appendLF=FALSE)
      
      # Locating files on virtual machine.
      dir <- study[j]
      if(study[j]=="LLS0"){dir <- "LLS/660Q"}
      if(study[j]=="LLS1"){dir <- "LLS/OmniExpr"}
      if(study[j]=="NTR0"){dir <- "NTR/Affy6"}
      if(study[j]=="NTR1"){dir <- "NTR/GoNL"}
      path0 <- file.path("/mnt/virdir/Backup/RP3_data/HRCv1.1_Imputation",dir)
      path1 <- path
      file0 <- paste0("chr",chr,".dose.vcf.gz")
      file1 <- paste0(study[j],".chr",chr,".dose.vcf.gz")
      file2 <- paste0(study[j],".chr",chr,".dose.vcf")
      
      # Reading in files.
      #vcf <- vcfR::read.vcfR(file=file.path(path1,file2),skip=skip[i],nrows=size,verbose=FALSE)
      vcf <- vcfR::read.vcfR(file=file.path(path0,file0),skip=skip[i],nrows=size,verbose=FALSE)
      vcf <- vcf[vcf@fix[,"CHROM"]!="",] # bug fix
      vcf@fix[,"ID"] <- paste0(vcf@fix[,"ID"],"_",seq_len(dim(vcf)["variants"]))
      collect[i,j][[1]] <- vcf
      stop <- dim(vcf)["variants"]==0
      final <- dim(vcf)["variants"]<size
      if(stop){break}
    }
    print(utils::object.size(collect),units="Gb")
    end <- Sys.time()
    if(stop){break}
    
    # only retain SNPs measured in all studies
    position <- lapply(seq_along(study),function(j) collect[i,j][[1]]@fix[,"POS"])
    common <- Reduce(f=intersect,x=position)
    for(j in seq_along(study)){
      cond <- match(x=common,table=position[[j]])
      collect[i,j][[1]] <- collect[i,j][[1]][cond,]
    }
    
    # Calculating minor allele frequency.
    num <- numeric(); maf <- list()
    for(j in seq_along(study)){
      num[j] <- dim(collect[i,j][[1]])["gt_cols"] # replace by adjusted sample sizes?
      maf[[j]] <- num[j]*vcfR::maf(collect[i,j][[1]])[,"Frequency"]
    }
    cond <- rowSums(do.call(what="cbind",args=maf))/sum(num)>0.05
    if(sum(cond)==0){if(final){break}else{next}}
    
    # Extracting genotypes.
    for(j in seq_along(study)){
      gt <- vcfR::extract.gt(collect[i,j][[1]][cond,])
      gt[gt=="0|0"] <- 0
      gt[gt=="0|1"|gt=="1|0"] <- 1
      gt[gt=="1|1"] <- 2
      storage.mode(gt) <- "integer"
      collect[i,j][[1]] <- gt
    }
    
    if(final){break}
  }
  
  # Removing empty rows.
  cond <- apply(collect,1,function(x) all(sapply(x,length)==0))
  collect <- collect[!cond,,drop=FALSE]
  
  # Fusing all matrices.
  snps <- NULL
  for(i in seq_len(nrow(collect))){
    inner <- NULL
    for(j in seq_len(ncol(collect))){
      add <- collect[i,j][[1]]
      colnames(add) <- paste0(colnames(collect)[j],":",colnames(add))
      inner <- cbind(inner,add)
    }
    snps <- rbind(snps,inner)
  }
  attributes(snps)$time <- end-start
  rownames(snps) <- sapply(strsplit(x=rownames(snps),split="_"),function(x) x[[1]])
  snps <- t(snps)
  
  # Filter samples.
  rownames(snps) <- sub(x=rownames(snps),pattern="LLS0|LLS1",replacement="LLS")
  rownames(snps) <- sub(x=rownames(snps),pattern="NTR0|NTR1",replacement="NTR")
  
  if(is.null(biobank)){
    save(object=snps,file=file.path(path1,paste0("BBMRI.chr",chr,".RData")))
  } else {
    save(object=snps,file=file.path(path1,paste0(biobank,".chr",chr,".RData")))
  }
}


#' @export
#' @title
#' Get exon data (Geuvadis)
#' 
#' @description
#' This function transforms exon data (virtual machine):
#' retains exons on the autosomes,
#' labels exons in the format "chromosome_start_end",
#' extracts corresponding gene names
#' 
#' @param path
#' data directory 
#' 
#' @examples
#' path <- "/virdir/Scratch/arauschenberger/trial"
#' 
get.exons.geuvadis <- function(path=getwd()){
  
  nrows <- 303544
  file <-"/virdir/Scratch/rmenezes/data_counts.txt"
  exons <- utils::read.table(file=file,header=TRUE,nrows=nrows)
  exons <- exons[exons[,"chr"] %in% 1:22,] # autosomes
  rownames(exons) <- exon_id <- paste0(exons[,"chr"],"_",exons[,"start"],"_",exons[,"end"])
  gene_id <- as.character(exons[,4])
  exons <- t(exons[,-c(1:4)])
  
  save(list=c("exons","exon_id","gene_id"),file=file.path(path,"Geuvadis.exons.RData"))
}


#' @export
#' @title
#' Get exon data (BBMRI)
#' 
#' @description
#' This function transforms exon data (virtual machine):
#' loads quality controlled gene expression data,
#' extracts sample identifiers,
#' removes samples without SNP data,
#' loads exon expression data,
#' extracts sample identifiers,
#' retains samples that passed quality control,
#' retains exons on the autosomes
#' 
#' @param path
#' data directory 
#' 
#' @examples
#' path <- "/virdir/Scratch/arauschenberger/trial"
#' 
get.exons.bbmri <- function(path=getwd()){
  
  # sample identifiers:
  # (1) loading quality controlled gene expression data 
  # (2) extracting sample identifiers
  # (3) removing identifiers without SNP data
  # (4) translating identifiers
  utils::data(rnaSeqData_ReadCounts_BIOS_cleaned,package="BBMRIomics") # (1)
  cd <- SummarizedExperiment::colData(counts)[,c("biobank_id","imputation_id","run_id")] # (2)
  counts <- NULL
  names(cd) <- substr(names(cd),start=1,stop=3) # abbreviate names
  cd <- cd[!is.na(cd$imp),] # (3)
  cd$id <- NA # (4)
  cd$id[cd$bio=="CODAM"] <- sapply(strsplit(x=cd$imp[cd$bio=="CODAM"],split="_"),function(x) x[[2]])
  cd$id[cd$bio=="LL"] <- sub(pattern="1_LLDeep_",replacement="",x=cd$imp[cd$bio=="LL"])
  cd$id[cd$bio=="LLS"] <- sapply(strsplit(x=cd$imp[cd$bio=="LLS"],split="_"),function(x) x[[2]])
  cd$id[cd$bio=="NTR"] <- sapply(strsplit(x=cd$imp[cd$bio=="NTR"],split="_"),function(x) x[[2]])
  cd$id[cd$bio=="PAN"] <- cd$imp[cd$bio=="PAN"]
  cd$id[cd$bio=="RS"] <- sub(pattern="RS1_|RS2_|RS3_",replacement="",x=cd$imp[cd$bio=="RS"])
  
  # Identify individual not with "id" but with "bio:id".
  any(duplicated(cd$id)) # TRUE
  sapply(unique(cd$bio),function(x) any(duplicated(cd$id[x]))) # FALSE
  
  # exon data:
  # (1) loading exon expression data
  # (2) extracting sample identifiers
  # (3) retaining autosomes
  # (4) retaining samples from above
  load("/virdir/Backup/RP3_data/RNASeq/v2.1.3/exon_base/exon_base_counts.RData") # (1)
  colnames(counts) <- sub(pattern=".exon.base.count.gz",replacement="",x=colnames(counts)) # (2)
  autosomes <- sapply(strsplit(x=rownames(counts),split="_"),function(x) x[[1]] %in% 1:22) # (3)
  exons <- counts[autosomes,cd$run] # (3) and (4)
  exon_id <- exon_id[autosomes] # (3)
  gene_id <- gene_id[autosomes] # (3)
  colnames(exons) <- paste0(cd$bio,":",cd$id)
  exons <- t(exons)
  
  save(list=c("exons","exon_id","gene_id"),file=file.path(path,"BBMRI.exons.RData"))
}


#' @export
#' @title
#' Prepare data matrices
#' 
#' @description
#' This function removes duplicate samples from each matrix,
#' only retains samples appearing in all matrices,
#' and brings samples into the same order.
#' 
#' @param ...
#' matrices with samples in the rows and variables in the columns,
#' with sample identifiers as rows names
#' 
#' @param message
#' display messages\strong{:} logical
#' 
#' @examples
#' X <- matrix(rnorm(6),nrow=3,ncol=2,dimnames=list(c("A","B","C")))
#' Z <- matrix(rnorm(9),nrow=3,ncol=3,dimnames=list(c("B","A","B")))
#' match.samples(X,Z)
#' 
match.samples <- function(...,message=TRUE){
  
  list <- list(...)
  if(length(list)==1 & is.list(list[[1]])){list <- list[[1]]}
  if(is.null(names(list))){
    names(list) <- sapply(substitute(list(...))[-1],deparse)
  }
  names <- names(list)
  
  # check input
  cond <- sapply(list,function(x) !is.matrix(x))
  if(any(cond)){
    stop("Provide matrices!",call.=FALSE)
  }
  cond <- sapply(list,function(x) is.null(rownames(x)))
  if(any(cond)){
    stop("Provide row names!",call.=FALSE)
  }
  
  # remove duplicated samples
  duplic <- lapply(list,function(x) duplicated(rownames(x)))
  for(i in seq_along(list)){
    number <- round(100*mean(duplic[[i]]))
    if(message){message(number," duplicates in \"",names[i],"\"")}
    list[[i]] <- list[[i]][!duplic[[i]],,drop=FALSE]
  }
  
  # retain overlapping samples
  all <- Reduce(f=intersect,x=lapply(list,rownames))
  for(i in seq_along(list)){
    percent <- round(100*mean(rownames(list[[i]]) %in% all))
    if(message){message(percent,"% overlap in \"",names[i],"\"")}
    list[[i]] <- list[[i]][all,,drop=FALSE]
  }
  
  # check output
  cond <- sapply(list,function(x) any(duplicated(rownames(x))))
  if(any(cond)){
    stop("Duplicate samples!",call.=FALSE)
  }
  cond <- sapply(list,function(x) nrow(x)!=nrow(list[[1]]))
  if(any(cond)){
    stop("Different sample sizes!",call.=FALSE)
  }
  cond <- sapply(list,function(x) any(rownames(x)!=rownames(list[[1]])))
  if(any(cond)){
    stop("Different sample names!",call.=FALSE)
  }
  
  return(list)
}

#' @export
#' @title
#' Adjust library sizes
#' 
#' @description
#' This function adjusts RNA-seq expression data for different library sizes.
#' 
#' @param x
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (variables)
#' 
#' @examples
#' n <- 5; p <- 10
#' x <- matrix(rnbinom(n=n*p,mu=5,size=1/0.5),nrow=n,ncol=p)
#' x[1,] <- 10*x[1,]
#' adjust.samples(x)
#' 
adjust.samples <- function(x){
  if(!is.matrix(x)){
    stop("no matrix argument",call.=FALSE)
  }
  if(!is.numeric(x)){
    stop("no numeric argument",call.=FALSE)
  }
  if(!is.integer(x)&&any(round(x)!=x)){
    warning("non-integer values",call.=FALSE)
  }
  if(any(x<0)){
    warning("negative values",call.=FALSE)
  }
  n <- nrow(x); p <- ncol(x)
  lib.size <- rowSums(x)
  norm.factors <- edgeR::calcNormFactors(object=t(x),lib.size=lib.size)
  gamma <- norm.factors*lib.size/mean(lib.size)
  gamma <- matrix(gamma,nrow=n,ncol=p,byrow=FALSE)
  x <- x/gamma
  return(x)
}

#' @export
#' @title
#' Adjust exon length
#' 
#' @description
#' This function adjusts exon expression data for different exon lengths.
#' 
#' @param x
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (exons)
#' 
#' @param offset
#' exon length\strong{:} vector of length \eqn{p}
#' 
#' @param group
#' gene names\strong{:} vector of length \eqn{p}
#' 
#' @examples
#' NA
#' 
adjust.variables <- function(x,offset,group){
  if(!is.numeric(x)|!is.matrix(x)){
    stop("Argument \"x\" is no numeric matrix.",call.=FALSE)
  }
  if(!is.numeric(offset)|!is.vector(offset)){
    stop("Argument \"offset\" is no numeric vector.",call.=FALSE)
  }
  if(any(offset<0)){
    stop("Argument \"offset\" takes negative values",call.=FALSE)   
  }
  if(!is.character(group)|!is.vector(group)){
    stop("Argument \"group\" is no character vector.",call.=FALSE)
  }
  if(ncol(x)!=length(group)|ncol(x)!=length(offset)){
    stop("Contradictory dimensions.",call.=FALSE)
  }
  n <- nrow(x); p <- ncol(x); names <- dimnames(x)
  x <- as.numeric(x)
  offset <- rep(offset,each=n)
  group <- strsplit(group,split=",")
  group <- sapply(group,function(x) x[[1]][1])
  group <- rep(group,each=n)
  lmer <- lme4::lmer(x ~ offset + (1|group)); gc()
  x <- matrix(stats::residuals(lmer),nrow=n,ncol=p,dimnames=names)
  x <- x-min(x)
  return(x)
}

#' @export
#' @title
#' Search for genes
#' 
#' @description
#' This function retrieves all protein-coding genes on a chromosome.
#' 
#' @param chr
#' chromosome\strong{:} integer 1-22
#' 
#' @param path
#' path to gene transfer format files (.gtf)
#' 
#' @param release
#' character "NCBI36", "GRCh37", or "GRCh38"
#' 
#' @param build
#' integer 49-91
#' 
#' @examples
#' NA
#' 
map.genes <- function(chr,path=getwd(),release="GRCh37",build=71){
  
  # check input
  if((!is.null(chr)) && (!chr %in% 1:22)){
    stop("Invalid argument \"chr\".",call.=FALSE)
  }
  if(!release %in% c("NCBI36","GRCh37","GRCh38")){
    stop("Invalid argument \"release\".",call.=FALSE)
  }
  if(!build %in% 49:91){
    stop("Invalid argument \"build\".",call.=FALSE)
  }
  
  file <- paste0("Homo_sapiens.",release,".",build,".gtf")
  if(!file.exists(file.path(path,file))){
    url <- paste0("ftp://ftp.ensembl.org/pub/release-",build,
                  "/gtf/homo_sapiens/",file,".gz")
    destfile <- file.path(path,paste0(file,".gz"))
    utils::download.file(url=url,destfile=destfile,method="auto")
    R.utils::gunzip(filename=destfile,remove=FALSE,overwrite=TRUE)
  }
  object <- refGenome::ensemblGenome()
  refGenome::basedir(object) <- path
  refGenome::read.gtf(object,filename=file)
  x <- refGenome::getGenePositions(object=object,by="gene_id")
  if(is.null(chr)){
    x <- x[x$gene_biotype=="protein_coding",]
  } else {
    x <- x[x$seqid==chr & x$gene_biotype=="protein_coding",]
  }
  x <- x[,c("gene_id","seqid","start","end","gene_name")] # added "gene_name"
  rownames(x) <- NULL
  colnames(x)[colnames(x)=="seqid"] <- "chr"
  return(x)
}

#' @export
#' @title
#' Search for exons
#' 
#' @description
#' This function attributes exons to genes.
#' 
#' @param gene
#' gene names\strong{:} vector with one entry per gene,
#' including the gene names
#' 
#' @param exon
#' exon names\strong{:} vector with one entry per exon,
#' including the corresponding \emph{gene} names
#' (separated by comma if multiple gene names)
#' 
#' @details
#' The exon names should contain the gene names. For each gene, this function
#' returns the indices of the exons.
#' 
#' @examples
#' gene <- c("A","B","C")
#' exon <- c("A","A,B","B","B,C","C")
#' map.exons(gene,exon)
#'
map.exons <- function(gene,exon){
  p <- length(gene)
  x <- list()
  pb <- utils::txtProgressBar(min=0,max=p,style=3)
  for(i in seq_len(p)){
    utils::setTxtProgressBar(pb=pb,value=i)
    which <- as.integer(grep(pattern=gene[i],x=exon))
    x[[i]] <- which
  }
  close(con=pb)
  names(x) <- gene
  return(x)
}


#' @export
#' @title
#' Search for SNPs
#' 
#' @description
#' This function attributes SNPs to genes.
#' 
#' @param gene.chr
#' chromosome\strong{:}
#' numeric vector with one entry per gene
#' 
#' @param gene.start
#' start position\strong{:}
#' numeric vector with one entry per gene
#' 
#' @param gene.end
#' end position\strong{:}
#' numeric vector with one entry per gene
#' 
#' @param snp.chr
#' integer 1-22
#' 
#' @param snp.pos
#' chromosomal position of SNPs\strong{:}
#' numeric vector with one entry per SNP
#' 
#' @param dist
#' number of base pairs before start position\strong{:}
#' integer
#' 
#' @examples
#' gene.chr <- rep(1,times=5)
#' gene.start <- 1:5
#' gene.end <- 2:6
#'
#' snp.chr <- rep(1,times=100)
#' snp.pos <- seq(from=1,to=4.9,length.out=100)
#' 
#' map.snps(gene.chr,gene.start,gene.end,snp.chr,snp.pos,dist=0)
#'
map.snps <- function(gene.chr,gene.start,gene.end,snp.chr,snp.pos,dist=10^3){
  if(length(gene.chr)!=length(gene.start)|length(gene.chr)!=length(gene.end)){
    stop("Invalid.",call.=FALSE)
  }
  if(!is.numeric(snp.chr)|!is.numeric(snp.pos)){
    stop("Invalid.",call.=FALSE)
  }
  p <- length(gene.start)
  x <- data.frame(from=integer(length=p),to=integer(length=p))
  pb <- utils::txtProgressBar(min=0,max=p,style=3)
  for(i in seq_len(p)){ # 
    utils::setTxtProgressBar(pb=pb,value=i)
    chr <- snp.chr == gene.chr[i]
    if(!any(chr)){next}
    start <- snp.pos >= (gene.start[i] - dist)
    end <- snp.pos <= gene.end[i] + 0
    which <- as.integer(which(chr & start & end))
    if(length(which)==0){next}
    x$from[i] <- min(which)
    x$to[i] <- max(which)
    if(length(which)==1){next}
    if(!all(diff(which)==1)){stop("SNPs are in wrong order!")}
  }
  close(con=pb)
  return(x)
}


#' @export
#' @title
#' Drop trivial tests
#' 
#' @description
#' This function trops trivial tests.
#' 
#' @param map
#' list with names "genes", "exons", and "snps"
#' (output from \code{\link{map.genes}}, \code{\link{map.exons}},
#' and \code{\link{map.snps}})
#' 
#' @details
#' This functions drops tests for genes without SNPs or with a single exon.
#' 
#' @examples
#' NA
#' 
drop.trivial <- function(map){
  
  # check input
  if(length(map)!=3){
    stop("Unexpected argument length.",call.=FALSE)
  }
  if(any(names(map)!=c("genes","exons","snps"))){
    stop("Unexpected argument names.",call.=FALSE)
  }
  
  # search
  p <- nrow(map$genes)
  pass <- rep(NA,times=p)
  pb <- utils::txtProgressBar(min=0,max=p,style=3)
  for(i in seq_len(p)){
    utils::setTxtProgressBar(pb=pb,value=i)
    ys <- map$exons[[i]]
    check <- logical()
    # Exclude genes without SNPs:
    check[1] <- map$snps$from[i] > 0
    check[2] <- map$snps$to[i] > 0
    # Exclude genes with single exon:
    check[3] <- length(ys) > 1
    pass[i] <- all(check)
  }
  close(con=pb)
  
  # check output
  if(any(pass[map$snps$to==0 & map$snps$from==0])){
    stop("Genes without any SNPs.",call.=FALSE)
  }
  if(any(pass[sapply(map$exons,length)<2])){
    stop("Genes without multiple exons.",call.=FALSE)
  }
  
  map$genes <- map$genes[pass,]
  map$exons <- map$exons[pass]
  map$snps <- map$snps[pass,]
  return(map)
}

#' @export
#' @title
#' Plot SNP-exon correlations
#' 
#' @description
#' This function plot the absolute Spearman correlation coefficients
#' between SNP and exons.
#' 
#' @param Y
#' exon expression\strong{:}
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (exons)
#' 
#' @param X
#' SNP genotype\strong{:}
#' matrix with \eqn{n} rows (samples) and \eqn{q} columns (SNPs)
#' 
#' @param map
#' list with names "genes", "exons", and "snps"
#' (output from \code{map.genes}, \code{map.exons}, and \code{map.snps})
#'
#' @param i
#' gene index\strong{:}
#' integer between \eqn{1} and \code{nrow(map$genes)}
#' 
#' @examples
#' # see vignette
#' 
visualise <- function(Y,X,map,i){
  
  # correlation
  ys <- map$exons[[i]]
  y <- as.matrix(Y[,ys,drop=FALSE])
  xs <- map$snps$from[i]:map$snps$to[i]
  x <- X[,xs,drop=FALSE]
  cor <- matrix(NA,nrow=length(ys),ncol=length(xs))
  for(j in seq_along(ys)){
    for(k in seq_along(xs)){
      cor[j,k] <- abs(cor(y[,j],x[,k],method="spearman"))
    }
  }
  
  # plot image
  graphics::par(mar=c(2,2,2,1))
  k <- 9
  inc <- 1.2 # Set to 1 for equal spacing.
  d <- 1/sum(inc^(0:k))
  breaks <- c(0,cumsum(inc^(0:k)*d))
  colour <- grDevices::colorRampPalette(c("white","darkblue"))(k+1)
  graphics::image(cor,xlab="",ylab="",breaks=breaks,col=colour,axes=FALSE)
  graphics::title(main=map$genes$gene_id[i],xlab="exons",ylab="SNPs",line=1)
  graphics::box()
  
  # indicate gene location (redundant if window approx. gene)
  #snp.loc <- as.numeric(sapply(strsplit(x=colnames(X),split=":"),function(x) x[[2]]))
  #local <- which(map$genes$start[i] <= snp.loc[xs] &  snp.loc[xs] <= map$genes$end[i])
  #dx <- 1/(length(xs) - 1)
  #dy <- 1/(length(ys) - 1)
  #y0 <- min(local)/length(xs) - 0.5*dx
  #y1 <- max(local)/length(xs) + 0.5*dx
  #col <- "black"; lwd <- 1
  #segments(x0=-0.5*dy,y0=y0,x1=-0.5*dy,y1=y1,col=col,lwd=4)
  #segments(x0=1+0.5*dy,y0=y0,x1=1+0.5*dy,y1=y1,col=col,lwd=4)
  #segments(x0=-0.5*dy,y0=y0,x1=1+0.5*dy,y1=y0,col=col,lwd=lwd,lty=2)
  #segments(x0=-0.5*dy,y0=y1,x1=1+0.5*dy,y1=y1,col=col,lwd=lwd,lty=2)
  
}


#' @export
#' @title
#' Plot grid
#' 
#' @description
#' This function combines a scatter plot with a density plot.
#' 
#' @param x
#' vector
#' 
#' @param y
#' vector
#' 
#' @param n
#' grid resolution
#' 
#' @param ...
#' graphical parameters
#' 
#' @examples
#' x <- stats::rbeta(n=100,shape1=0.3,shape2=0.5)
#' y <- stats::rbeta(n=100,shape1=0.3,shape2=0.5)
#' grid(x,y)
#' 
grid <- function(x,y,n=10,...){
  #default values
  par <- list(...)
  if(is.null(par$xlim)){par$xlim <- range(x)}
  if(is.null(par$ylim)){par$ylim <- range(y)}
  if(is.null(par$pch)){par$pch <- 16}
  if(is.null(par$cex)){par$cex <- 0.5}
  if(is.null(par$xlab)){par$xlab <- "x"}
  if(is.null(par$ylab)){par$ylab <- "y"}
  
  # open plot
  graphics::plot.new()
  graphics::plot.window(xlim=par$xlim,ylim=par$ylim)
  graphics::box()
  graphics::axis(side=1)
  graphics::axis(side=2)
  
  # density
  xc <- seq(from=par$xlim[1],to=par$xlim[2],length.out=n+1)
  yc <- seq(from=par$ylim[2],to=par$ylim[1],length.out=n+1)
  M <- matrix(integer(),nrow=n,ncol=n)
  for(i in seq_len(n)){
    for(j in seq_len(n)){
      M[i,j] <- sum(x >= xc[i] & x <= xc[i+1] & y <= yc[j] & y >= yc[j+1])
    }
  }
  M <- M/(1.25*max(M))
  
  # fill plot
  for(i in seq_len(n)){
    for(j in seq_len(n)){
      graphics::polygon(x=c(xc[i],xc[i],xc[i+1],xc[i+1]),
                        y=c(yc[j],yc[j+1],yc[j+1],yc[j]),
                        col=grDevices::gray(level=1-M[i,j]),border=NA)
    }
  }
  graphics::segments(x0=xc,y0=par$ylim[1],y1=par$ylim[2],col="white")
  graphics::segments(x0=par$xlim[1],x1=par$xlim[2],y0=yc,col="white")
  graphics::points(x=x,y=y,pch=par$pch,cex=par$cex,col="black")
  graphics::title(xlab=par$xlab,ylab=par$ylab)
}




# test.trial <- function(y,x,limit=NULL,steps=NULL,rho=c(0,0.5,1)){
#     
#     if(is.null(limit)){limit <- 5}
#     if(is.null(steps)){steps <- c(10,20,20,50)}
#     
#     # check input
#     if(!is.numeric(limit)){
#         stop("Argument \"limit\" is not numeric.",call.=FALSE)
#     }
#     if(limit<1){
#         stop("Argument \"limit\" is below one.",call.=FALSE)
#     }
#     if(!is.numeric(steps)|!is.vector(steps)){
#         stop("Argument \"steps\" is no numeric vector.",call.=FALSE)
#     }
#     if(sum(steps)<2){
#         stop("Too few permutations \"sum(steps)\".",call.=FALSE)
#     }
#     
#     # test effects
#     pvalue <- rep(x=NA,times=length(rho))
#     for(j in seq_along(rho)){
#         tstat <- spliceQTL:::G2.multin(
#             dep.data=y,indep.data=x,nperm=steps[1]-1,rho=rho[j])$Sg
#         for(nperm in steps[-1]){
#             tstat <- c(tstat,spliceQTL:::G2.multin(
#                 dep.data=y,indep.data=x,nperm=nperm,rho=rho[j])$Sg[-1])
#             if(sum(tstat >= tstat[1]) >= limit){break}
#         }
#         pvalue[j] <- mean(tstat >= tstat[1],na.rm=TRUE)
#     }
#     
#     return(pvalue)
# }
