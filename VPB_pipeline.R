#!/usr/bin/env Rscript

# VANAPRABHAVA pipeline, only UNIX compatible
# v1.0 created on Dec. 30th 2015
# v1.1 Dec. 31th 2015 -- now can deal with simple conditions between parameters
# remi (dot) tournebize (at) gmail (dot) com
# command-line in UNIX: Rscript path_to_script/VPB.R path_to_the_parameter_file
# make sure that VPB.R script is executable (chmod +750)
argv <- commandArgs(trailingOnly=TRUE)

inputFile <- argv[1]
if (!file.exists(inputFile)) stop("ERROR. Input file does not exist! Exiting.")

# specific interest, so theta = 2 * No * mu (there is no consideration for TWO copies of
# genes per individual)
coalescent.ploidy <- "haploid" # either "haploid" or "diploid"
# please note that for demophylogenetical models, we are working with individuals
# which are the replicators in consideration and there is one copy of each of
# these replicators, so set ploidy to "haploid"
# (even if those are pseudo-haploid individuals in reality)

#___ outDir/
#    |___ temp/
#    |___ phyloTrees/{input_prm_file}.newick
#    |___ SFS/{input_prm_file}.sfs
#    |___ {input_prm_file}.priors

#####################################
#####################################
#####################################
#####################################

replace <- function(X, vals) {
  if (!is.null(X)) {
    if (is.null(dim(X))) X <- matrix(X, nrow=1)
    X <- apply(X, c(1,2), function(x) {
      w <- which(vals[,1]==x)
      return(ifelse(length(w)>0, as.numeric(vals[w,2]), as.numeric(x)))
    })
  }
  if (length(X)==0) X <- NULL
  return(X)
}

generate.ms.command <- function(ms_path, inputFile, firstRun, outDir, coalescent.ploidy) {
  
  if (firstRun==TRUE) {
    con <- file(inputFile, open="r")
    
    n <- Ne <- ej <- en <- ep <- priors <-  c()
    F <- FALSE; read <- c("mu"=F, "n"=F, "Ne"=F, "ej"=F, "en"=F, "ep"=F, "priors"=F)
    while( length(l <- readLines(con, n=1, warn=FALSE)) > 0 ) {
      
      if ( grepl("^\\[", l, perl=T) ) read[read] <- FALSE
      if ( grepl("\\[MUTATION RATE", l) ) read["mu"] <- TRUE
      if ( grepl("\\[DEME SAMPLE SIZES", l) ) read["n"] <- TRUE
      if ( grepl("\\[PRESENT DEME EFFICIENT SIZES", l) ) read["Ne"] <- TRUE
      if ( grepl("\\[DIVERGENCE EVENTS", l) ) read["ej"] <- TRUE
      if ( grepl("\\[PUNCTUAL EVENTS", l) ) read["en"] <- TRUE
      if ( grepl("\\[PERIODIC EVENTS", l) ) read["ep"] <- TRUE
      if ( grepl("\\[PRIORS", l) ) read["priors"] <- TRUE
      
      if ( !grepl("^#", l, perl=T) && !grepl("^\\[", l, perl=T) && nchar(l) >= 1 ) {
        if (read["mu"]) mu <- l
        if (read["n"]) n <- c(n, l)
        if (read["Ne"]) Ne <- c(Ne, l)
        if (read["ej"]) ej <- rbind(ej, unlist(strsplit(l, " ")))
        if (read["en"]) en <- rbind(en, unlist(strsplit(l, " ")))
        if (read["ep"]) ep <- rbind(ep, unlist(strsplit(l, " ")))
        if (read["priors"]) priors <- rbind(priors, unlist(strsplit(l, " ")))
      }
    }
    close(con)
    LPRM <- list(mu=mu, n=n, Ne=Ne, ej=ej, en=en, ep=ep, priors=priors)
    save(LPRM, file=paste(outDir,"/temp/ListPRM.bin",sep=""))
  } else {
    load(paste(outDir,"/temp/ListPRM.bin",sep=""))
    for (l in seq_along(LPRM)) {
      assign(names(LPRM)[l], LPRM[[l]])
    }
  }
  
  # replace all data by their values
  if (length(priors)>0) {
    vals <- as.data.frame(cbind(priors[,1], NA), stringsAsFactors=F)
    names(vals) <- c("param","val"); vals$val <- as.numeric(vals$val)
    for ( i in 1:nrow(priors) ) {
      if (priors[i,3] == "unif") {
        vals[i,2] <- runif(1, as.numeric(priors[i,4]), as.numeric(priors[i,5]))
      } else if (priors[i,3] == "logunif") {
        logval <- runif(1, log10(as.numeric(priors[i,4])), log10(as.numeric(priors[i,5])))
        vals[i,2] <- 10**logval
      } else {
        stop("Unrecognized distribution law...")
      }
      if (priors[i,2]=="1") vals[i,2] <- round(vals[i,2])
    }
  } else { vals <- NULL }
  
  # replace values
  en <- replace(en, vals)
  Ne <- replace(Ne, vals)
  mu <- replace(mu, vals)
  ej <- replace(ej, vals)
  ep <- replace(ep, vals)
  
  # set the ploidy factor
  if (coalescent.ploidy=="haploid") {
    Ne <- Ne / 2
  }
  
  # get No
  No <- as.numeric(Ne[1])
  
  get.past.events.command <- function(ej, en, eno, No) {
    
    EE <- matrix(c(0, as.numeric(Ne)/No), nrow=1)
    if (!is.null(en)) {
      # set the values of -en as a fraction of No for each island
      for (i in 1:nrow(en)) {
        en[i,-1] <- as.numeric(Ne) * as.numeric(en[i,-1]) / No
      }
      # rbind
      EE <- rbind(EE, en)
    }
    EE <- as.data.frame(EE, stringsAsFactors=F)
    names(EE)[1] <- "time"; EE$time <- as.numeric(EE$time)
    # scale times
    EE$time <- EE$time / ( 4 * No )
    if (!is.null(ej)) {
      ej <- as.data.frame(ej, stringsAsFactors=F)
      names(ej)[1] <- "time"; ej$time <- as.numeric(ej$time)
      # scale times
      ej$time <- ej$time / ( 4 * No )
    }
    
    # sort chronologically
    EE <- EE[order(as.numeric(EE$time)),]
    
    if (!is.null(ej)) {
      for ( i in 1:nrow(ej) ) {
        # rescale if TRUE and add the new corresponding -en event
        if ( ej[i,4] == 1 ) {
          prev.en <- rev(which(EE$time <= ej$time[i]))[1]
          new.en <- EE[prev.en,]
          new.en[c(1,as.integer(ej[i,3])+1)] <- c(ej$time[i],
                                                  as.numeric(EE[prev.en,as.integer(ej[i,2])+1]) + as.numeric(EE[prev.en,as.integer(ej[i,3])+1])) 
          # insert the new -en event
          EE <- rbind(EE, new.en)
          # always re-sort
          EE <- EE[order(as.numeric(EE$time)),]
        }
      }
    }
    
    # print
    rownames(EE) <- NULL
    colnames(EE)[-1] <- paste("ISL_",2:ncol(EE)-1,sep="")
    print("Past historical events:")
    print(EE)
    
    # remove stable Ne events on a time-sorted EE matrix
    NEE <- EE
    if ( nrow(EE) > 1) {
      for ( i in 2:nrow(EE) ) {
        dif <- as.numeric(EE[i,-1]) - as.numeric(EE[i-1,-1])
        NEE[i, c(FALSE, dif==0)] <- NA
      }
    }
    NEE[1,2] <- NA # because it is No
    
    if (!is.null(ej)) {
      for ( i in 1:nrow(ej) ) {
        post.en <- which(EE$time > ej$time[i])
        # set to NE previous ones
        if (length(post.en)>0) NEE[post.en, as.integer(ej[i,2])+1] <- NA
      }
    }
    
    # get the commands for ej and en
    # get command for -en
    c.en <- c.en.nama <- c()
    for ( i in 1:nrow(NEE) ) for ( j in 2:ncol(NEE) ) {
      if (!is.na(NEE[i,j])) {
        c.en <- c(c.en, paste("-en ",NEE[i,1]," ",j-1," ",NEE[i,j],sep="")) 
        c.en.nama <- c(c.en.nama, NEE[i,1])  
      }
    }
    names(c.en) <- c.en.nama
    
    # get command for -ej
    if (!is.null(ej)) {
      c.ej <- c.ej.nama <- c()
      for ( i in 1:nrow(ej) ) {
        c.ej <- c(c.ej, paste("-ej ",ej[i,1]," ",ej[i,2]," ",ej[i,3],sep=""))
        c.ej.nama <- c(c.ej.nama, ej[i,1])
      }
      names(c.ej) <- c.ej.nama
    }
    
    if (!is.null(ej)) EE <- c(c.en, c.ej) else EE <- c.en
    EE <- EE[order(as.numeric(names(EE)), EE)]
    EE <- paste(EE, collapse=" ")
    rownames(EE) <- NULL
    
    return(EE)
  }
  
  # get fixed command part
  ms <- paste(ms_path," ",sum(as.numeric(n))," 1 -t ",4 * No * mu,sep="")
  
  # get command for -I
  if ( length(Ne) > 1 ) ms <- paste(ms, " -I ",length(Ne)," ",paste(n, collapse=" "),sep="")
  
  # get past events MS command
  ## sort the historical events and Ne-div-rescale if requested by the user
  EE.ms <- get.past.events.command(ej, en, eno, No)
  
  # get command tail
  tail.ms <- paste("-T -s 1 | grep \"(\" > ",outDir,"/temp/ms.txt",sep="")
  
  # assemble the MS command
  ms <- paste(c(ms, EE.ms, tail.ms), collapse=" ")
  
  return(list(ms=ms, vals=vals, mu=c(mu), No=No, n=n))
}

#####################################
#####################################
#####################################
#####################################

# GET [ADMIN] PARAM VALUES
# default values
verbose <- 1
drawTrees <- 0
forceUltrametricity <- 1
empty_previous_files <- 0
# read [ADMIN] parameters
con <- file(inputFile, open="r")
read <- FALSE
while( length(l <- readLines(con, n=1, warn=FALSE)) > 0 ) {
  if ( grepl("^\\[ADMIN", l, perl=T) ) { read <- TRUE; next }
  if ( read == TRUE && !grepl("^#", l, perl=T) && !grepl("^\\[", l, perl=T) && nchar(l) >= 1 ) {
    l <- gsub(" = ","=",l)
    l <- unlist(strsplit(l, "="))
    if (length(l)==2) {
      assign(l[1], l[2])
    } else {
      stop("ERROR. The value of an [ADMIN] parameter is left missing! Exiting.")
    }
  } else if ( grepl("^\\[", l, perl=T) ) read <- FALSE
}
close(con)
# set types
simul_start <- as.integer(simul_start)
simul_end <- as.integer(simul_end)
verbose <- as.integer(verbose)
drawTrees <- as.integer(drawTrees)
forceUltrametricity <- as.integer(forceUltrametricity)
empty_previous_files <- as.integer(empty_previous_files)

# check the absence of spaces in paths
if ( grepl(" ",ms_path) || !file.exists(ms_path) ) stop("ERROR. Found a space in ms_path or file does not exist. Exiting.")
if ( grepl(" ",VPB_path) || !file.exists(VPB_path) ) stop("ERROR. Found a space in VPB_path or file does not exist. Exiting.")
if ( grepl(" ",outDir) ) stop("ERROR. Found a space in outDir path. Exiting.")

# create the outDir
if (!dir.exists(outDir)) { dir.create(paste(outDir,"/temp",sep=""), recursive=TRUE); print(paste("Created dir: ",outDir)) }
if (!dir.exists(paste(outDir,"/PHYLO",sep=""))) dir.create(paste(outDir,"/PHYLO",sep=""))
if (!dir.exists(paste(outDir,"/SFS",sep=""))) dir.create(paste(outDir,"/SFS",sep=""))
if (empty_previous_files==1) {
  f <- paste(outDir,"/PHYLO/",list.files(paste(outDir,"/PHYLO",sep="")),sep="")
  if (length(f)>0) x=sapply(f, function(x) if (file.exists(x)) unlink(x, recursive=FALSE))
  f <- paste(outDir,"/SFS/",list.files(paste(outDir,"/SFS",sep="")),sep="")
  if (length(f)>0) x=sapply(f, function(x) if (file.exists(x)) unlink(x, recursive=FALSE))
  if (file.exists(paste(outDir,"/temp/ms.txt",sep=""))) unlink(paste(outDir,"/temp/ms.txt",sep=""))
  if (file.exists(paste(outDir,"/temp/ListPRM.bin",sep=""))) unlink(paste(outDir,"/temp/ListPRM.bin",sep=""))
  print("Safely removed previous files from the outDir!")
  if (!is.null(warnings())) print(warnings())
}

# ad ./ to the script/exe paths if absent
if (!grepl("^\\.\\/", ms_path, perl=T) && !grepl("^\\/", ms_path, perl=T)) ms_path <- paste("./",ms_path,sep="")
if (!grepl("^\\.\\/", VPB_path, perl=T) && !grepl("^\\/", ms_path, perl=T)) VPB_path <- paste("./",VPB_path,sep="")
# set the permissions to the scripts
Sys.chmod(ms_path, "0750"); Sys.chmod(VPB_path, "0750")

# get the radical (= model strict name)
RADICAL <- rev(unlist(strsplit(inputFile, "/")))[1]
RADICAL <- unlist(strsplit(RADICAL, "\\."))
RADICAL <- paste(RADICAL[-length(RADICAL)], collapse=".")
print(RADICAL)

#####################################
#####################################
#####################################
#####################################

# set directory
if (!dir.exists(outDir)) stop("outDir does not exist. Exiting.")
inputFile <- normalizePath(inputFile)
setwd(outDir)
print(paste("Set directory to:",outDir))

print("Started at: ")
print(proc.time())

# sink to log file
log <- file(paste(outDir,"/",RADICAL,".log",sep=""), open = "wt")
sink(log, type="output", split=T)

for ( run in simul_start:simul_end ) {
  
  cat("\n\n")
  print("|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
  print("|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
  print(paste("[RUN]",run))
  print("|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
  print("|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||")
  
  firstRun <- ifelse(run==simul_start, TRUE, FALSE)
  
  # get the priors
  LOUT <- generate.ms.command(ms_path, inputFile, firstRun, outDir, coalescent.ploidy)
  for (l in seq_along(LOUT)) assign(names(LOUT)[l], LOUT[[l]])
  
  # export the priors
  vvals <- as.numeric(vals[,2]); names(vvals) <- vals[,1]
  vals <- vvals
  if (empty_previous_files==1 && firstRun==TRUE) {
    write.table(t(vals), paste(outDir,"/",RADICAL,".priors",sep=""), row.names=F, quote=F, sep="\t")
  } else {
    write.table(t(vals), paste(outDir,"/",RADICAL,".priors",sep=""), col.names=F, row.names=F, sep="\t", append=T, quote=F)
  }
  print("Parameter values:")
  print(t(vals))
  
  # get the VPB-python command
  python <- paste(VPB_path," ",
                  paste(outDir,"/temp/ms.txt",sep="")," ",
                  "-m ",sprintf("%.20f", mu)," ",
                  "-s ",4 * No," ",
                  "-I ",paste(n,collapse=" ")," ",
                  "-ophylo ",outDir,"/PHYLO/",RADICAL,"_",run,".newick ",
                  "-osfs ",outDir,"/SFS/",RADICAL,"_",run,".sfs",sep="")
  if (forceUltrametricity==1) python <- paste(python, "--forceUltrametric")
  if (drawTrees==1) python <- paste(python, "--drawTrees")
  if (verbose==0) python <- paste(python, "--quiet")
  
  # systemize
  ## MS
  cat("\n")
  print(ms)
  unlink(paste(outDir,"/temp/ms.txt",sep=""))
  MS.STDOUT <- system(ms, intern=TRUE)
  attr <- attributes(MS.STDOUT)
  if (!is.null(attr) && attr$status!=0) stop("ERROR on exit status of VPB.py")
    #print(MS.STDOUT)
  ## PYTHON
  cat("\n")
  print(python)
  VPB.STDOUT <- system(python, intern=TRUE)
  attr <- attributes(VPB.STDOUT)
  if (!is.null(attr) && attr$status!=0) stop("ERROR on exit status of VPB.py")
  print(VPB.STDOUT)
  
}

# end sinking to log file
sink()

print("Ended at: ")
print(proc.time())
print("ALL DONE")

