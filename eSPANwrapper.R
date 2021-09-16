eSPANwrapper <- function(Input_R1, Input_R2, ChIP_R1, ChIP_R2,
                         BrDU_R1, BrDU_R2, eSPAN_R1, eSPAN_R2,
                         
                         ExpTitle = "None",
                         slidingWindow = "YES", 
                         bin = "None",
                         stepSize = "None", 
                         PlotIPprofile = "Global",
                         PlotAveProfile = "all") {
  
  packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
                "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw",
                "BSgenome.Scerevisiae.UCSC.sacCer3")
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  #
  useDef <- function(a,d) ifelse(isTruthy(a), a,d)
  
  ExpTitle = useDef(ExpTitle, "None")
  bin = useDef(bin, "None")
  stepSize = useDef(stepSize, "None") 
  #
  E_Ori <- read.table("E_Rep.bed", header = TRUE, quote = "/t")
  L_Ori <- read.table("L_Rep.bed", header = TRUE, quote = "/t")
  All_Ori <- read.table("OriginList_Full.bed", header = TRUE, quote = "/t")
  All_Ori_Link <- "OriginList_Full.bed"
  image <- load.image("ForkPic_1.jpeg")
  #
  
  
  if(ExpTitle == "None"){
    Pro_1 <- unlist(strsplit(basename(Input_R1), split='_', fixed=TRUE))[[1]]
  } else {
    Pro_1 <- ExpTitle
  }
  
  
  message(paste0("Experiment: ", Pro_1))
  
  suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1)))   #create directory named with the protein in the Desktop
  
  #Quality check of fastqs'
  
  message("Running QC ...")
  
  if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "QR", ".html"))){
    
    fls = c(Input_R1, Input_R2, ChIP_R1, ChIP_R2, BrDU_R1, BrDU_R2, eSPAN_R1, eSPAN_R2)
    
    names(fls) = sub(".fastq", "", basename(fls))
    
    qas = lapply(seq_along(fls),
                 function(i, fls) qa(readFastq(fls[i]), names(fls)[i]),
                 fls)
    qa = do.call(rbind, qas)
    rpt = report(qa, dest = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "QR", ".html"))
    
  }
  
  ##
  
  message("Running alignments ...")
  message("Reference yeast genome : S288C")
  
  if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", "Bam"))){
    
    RunAlign <- function(File_R1, File_R2, SampName){
      
      tempdir(check = TRUE)
      
      Sam <- tempfile(fileext = ".sam")
      Bam <- tempfile(fileext = ".bam")
      nmCollate <- tempfile(fileext = ".bam")
      fixMat <- tempfile(fileext = ".bam")
      SrtBam <- tempfile(fileext = ".bam")
      
      
      Pro_1 <- Pro_1
      Pro_2 <- SampName
      
      suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Bam")))
      
      AlnLog <- paste0("~/Desktop/", Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".log")
      SFBam <- paste0("~/Desktop/", Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".bam")
      
      #read the indexed reference genome for the alignment of sequenced data
      ref_index <- "bowtie2-2.4.4-macos-x86_64/indexes/S288C_Ref"
      
      #following commands will run the alignemnt, check quality, sort, filter and index the resultant bam file 
      
      system(sprintf("(bowtie2-2.4.4-macos-x86_64/bowtie2 -p 8  --no-discordant --fr -x %s -1 %s -2 %s -S %s) 2> %s", 
                     ref_index, File_R1, File_R2, Sam, AlnLog))
      
      system(sprintf("samtools-1.13/samtools view -bS -@ 15 -q 30 -f 2 %s > %s", Sam, Bam))
      
      system(sprintf("samtools-1.13/samtools collate -@ 15 -o %s %s", nmCollate, Bam))
      
      system(sprintf("samtools-1.13/samtools fixmate -@ 15 -m %s %s", nmCollate, fixMat))
      
      system(sprintf("samtools-1.13/samtools sort -l 9 -@ 15 -m 1024M  -O bam -o %s %s", SrtBam, fixMat))
      
      system(sprintf("samtools-1.13/samtools markdup -@ 15 %s %s", SrtBam, SFBam))
      
      system(sprintf("samtools-1.13/samtools index -@ 15 %s", SFBam))
      
      unlink(c(Sam, Bam, nmCollate, fixMat, SrtBam), recursive = T, force = T)
      
    }
    
    RunAlign(Input_R1, Input_R2, "Input")
    RunAlign(ChIP_R1, ChIP_R2, "ChIP")
    RunAlign(BrDU_R1, BrDU_R2, "BrDU")
    RunAlign(eSPAN_R1, eSPAN_R2, "eSPAN")
    
  }
  
  # binSize
  
  bamFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Bam", "/", "*", ".bam"))
  
  if(bin == "None"){
    mFragS <- c()
    for(i in 1:length(bamFiles)) {
      binMFS <- round_any(mean(c(mFragS, mean(getPESizes(bamFiles[i], param=readParam(pe="both"))$sizes))), 100)
    }
    binSize <- binMFS
  } else {
    binSize <- as.numeric(bin)
  }
  
  message(paste0("Mean fragment size = ", binMFS, " bps"))
  
  message(paste0("Selected bin length = ", binSize, " bps"))
  
  # stepSize
  
  if(slidingWindow == "NO"){
    stepSize <- binSize
  } else {
    if(stepSize == "None"){
      stepSize <- 10
    } else {
      stepSize <- as.numeric(stepSize)
    }
  }
  
  #3 - Calculate genome-wide binned coverage
  
  message("Calculating read coverage ...")
  
  if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", "Coverage"))){
    
    bamFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Bam", "/", "*", ".bam"))
    
    BamCoverage <- function(bamFile, binSize = binSize, stepSize = 10, slidingWindow = "YES"){
      
      Pro_1 <- unlist(strsplit(tools::file_path_sans_ext(basename(bamFile)), split='_', fixed=TRUE))[[1]] #extract protein name
      Pro_2 <- unlist(strsplit(tools::file_path_sans_ext(basename(bamFile)), split='_', fixed=TRUE))[[2]] #extract sample name
      
      suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Coverage"))) 
      
      
      tempdir(check = TRUE)
      
      GenomFile <- tempfile(fileext = ".txt")
      binFile <- tempfile(fileext = ".bed")
      
      command_1 <- "samtools-1.13/samtools idxstats %s | awk 'BEGIN {OFS=\"\\t\"} {if ($2>0) print ($1,$2)}' >  %s"
      system(sprintf(command_1, bamFile, GenomFile))
      
      if(slidingWindow=="YES"){
        command_2 <- "bedtools2/bin/bedtools makewindows -g %s -w %s -s %s > %s"
        system(sprintf(command_2, GenomFile, binSize, stepSize, binFile))
      } else {
        command_2 <- "bedtools2/bin/bedtools makewindows -g %s -w %s > %s"
        system(sprintf(command_2, GenomFile, binSize, binFile))
      }
      
      pncFiles_watson <- tempfile(fileext = ".bed")
      pncFiles_crick <- tempfile(fileext = ".bed")
      
      command_3 <- "samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s | grep -v XS:i: | samtools-1.13/samtools view -@ 8 -b - | bedtools2/bin/bedtools genomecov -5 -d -ibam stdin -strand + | awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' > %s"
      command_4 <- "samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s | grep -v XS:i: | samtools-1.13/samtools view -@ 8 -b - | bedtools2/bin/bedtools genomecov -5 -d -ibam stdin -strand - | awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' > %s"
      
      system(sprintf(command_3, binFile, bamFile, paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"), pncFiles_watson))
      system(sprintf(command_4, binFile, bamFile, paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"), pncFiles_crick))
      
      finFiles_watson <- paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "watson.bed")
      finFiles_crick <- paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "crick.bed")
      
      command_5 <- "bedtools2/bin/bedtools map -a %s -b %s -null 0 -o sum | awk 'BEGIN {OFS=\"\\t\"} {if ($4>=0) print $1,$2,$3,\"%s\",$4}' > %s"
      
      system(sprintf(command_5, binFile, pncFiles_watson, paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"), finFiles_watson))
      system(sprintf(command_5, binFile, pncFiles_crick, paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"), finFiles_crick))
      
      unlink(c(GenomFile, binFile, pncFiles_watson, pncFiles_crick), recursive = T, force = T)
      
    }
    
    for(i in 1:length(bamFiles)){
      BamCoverage(bamFile = bamFiles[i], binSize = binSize)
    }
    
  }
  
  #4 - Calculate Ratio
  
  message("Calculating enrichment ratios ...")
  
  if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", "Ratios"))){
    
    CoverageFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", "*", ".bed"))
    
    CalculateRatio <- function(IP_coverage, Input_coverage){
      
      IP.df = read.table(IP_coverage, header = F)
      In.df = read.table(Input_coverage, header = F)
      
      suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Ratios"))) 
      
      
      IP_Sum <- sum(as.numeric(IP.df[,5]))
      In_Sum <- sum(as.numeric(In.df[,5]))
      corrFactor <- IP_Sum/In_Sum
      Ratio <- round(IP.df[,5]/In.df[,5]/corrFactor, 4)
      In.score.norm <- round(In.df[,5]*corrFactor)
      Ratio[!is.finite(Ratio)] <- 0
      
      strand <- cbind.data.frame(IP.df[,1], IP.df[,2], IP.df[,3], IP.df[,4], IP.df[,5], In.score.norm, Ratio)
      
      chroms <- unique(strand[,1])
      s_strand <- NULL
      for(i in 1:length(chroms)){
        chr <- strand[strand[,1]==chroms[i], ]
        x <- chr[,2]
        y <- chr$Ratio
        splineObject <- smooth.spline(x, y)
        chr$splineSmooth <- round(as.numeric(splineObject$y), 3)
        s_strand <- rbind.data.frame(s_strand, chr)
      }
      
      s_strand$ppois <- ppois(  q=s_strand[,5] - 1, 
                                lambda=s_strand$`In.score.norm`, 
                                lower.tail=FALSE, log=FALSE      )
      
      
      colnames(s_strand) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'ip.score', 'in.score', 'ratio', 'smooth', 'pvalue')
      
      ###
      
      write.table(s_strand, paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", s_strand$name[1], ".bed"),
                  quote=FALSE, row.names=FALSE, sep="\t")
      
    }
    
    CalculateRatio(CoverageFiles[6], CoverageFiles[6])
    CalculateRatio(CoverageFiles[4], CoverageFiles[6])
    CalculateRatio(CoverageFiles[2], CoverageFiles[6])
    CalculateRatio(CoverageFiles[8], CoverageFiles[2])
    
    CalculateRatio(CoverageFiles[5], CoverageFiles[5])
    CalculateRatio(CoverageFiles[3], CoverageFiles[5])
    CalculateRatio(CoverageFiles[1], CoverageFiles[5])
    CalculateRatio(CoverageFiles[7], CoverageFiles[1])
    
  }
  
  #5 - Define and process peaks
  
  message("Processing peaks ...")
  
  if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", "Peaks"))){
    
    bamFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Bam", "/", "*", ".bam"))
    
    ProcessPeaks <- function(bamFiles){
      
      PeakFinder <- function(IPBam, InBam){
        
        Input_combined <- read.table(pipe(sprintf("bedtools2/bin/bedtools bamtobed -i %s", InBam)) )
        Input_Plus <- read.table(pipe(sprintf("bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"+\") print $0}'", InBam)) )
        Input_Minus <- read.table(pipe(sprintf("bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"-\") print $0}'", InBam)) )
        
        IP_combined <- read.table(pipe(sprintf("bedtools2/bin/bedtools bamtobed -i %s", IPBam)) )
        IP_Plus <- read.table(pipe(sprintf("bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"+\") print $0}'", IPBam)) )
        IP_Minus <- read.table(pipe(sprintf("bedtools2/bin/bedtools bamtobed -i %s | awk '{OFS=\"\\t\"} {if($6==\"-\") print $0}'", IPBam)) )
        
        OriginPeaks <- function(IP_DF, Input_DF){
          
          tempdir(check = TRUE)
          IP <- tempfile(fileext=".bed")
          Input <- tempfile(fileext=".bed")
          outDir <- tempdir()
          peakFile <- tempfile(fileext = ".bed")
          ColHeads <- "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\""
          
          write.table(IP_DF, file = IP, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
          write.table(Input_DF, file = Input, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
          
          system(sprintf("macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n %s --outdir %s 2> /dev/null", IP, Input, "Peak", outDir))
          allPeaks <- read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#")
          
          #allPeaks <- suppressMessages(read.delim2(callpeak(tfile = IP, cfile = Input, gsize = 12157105, format = "BED", pvalue = 10e-6, nomodel = F, 
          #                                                  outdir = outDir, name = Pro_1)$outputs[3], comment.char="#"))
          
          write.table(allPeaks, file = peakFile, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
          
          Peaks_at_Origins <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'", 
                                                      peakFile, All_Ori_Link, ColHeads)), header = TRUE ) 
          Peaks_at_Origins <- Peaks_at_Origins[!duplicated(Peaks_at_Origins$oriName), ]
          
          Peaks_at_Origins$oriCenter <- round((Peaks_at_Origins$oriStart + Peaks_at_Origins$oriEnd)/2)
          
          return(Peaks_at_Origins)
          
          unlink(c(IP, Input, peakFile))
          unlink(outDir, recursive = TRUE)
          
        }
        
        Watson <- OriginPeaks(IP_DF=IP_Plus, Input_DF=Input_Plus)
        Crick <- OriginPeaks(IP_DF=IP_Minus, Input_DF=Input_Minus)
        Combo <- OriginPeaks(IP_DF=IP_combined, Input_DF=Input_combined)
        
        ResPeaks <- list(Watson, Crick, Combo)
        names(ResPeaks) <- c("watsonPeaks", "crickPeaks", "comboPeaks")
        
        return(ResPeaks)
        
      }
      
      BrDU_Input <- PeakFinder(IPBam=bamFiles[1], InBam=bamFiles[3])
      
      ChIP_Input <- PeakFinder(IPBam=bamFiles[2], InBam=bamFiles[3])
      
      eSPAN_ChIP <- PeakFinder(IPBam=bamFiles[4], InBam=bamFiles[2])
      
      #Define overlapping peaks
      
      tempdir(check = TRUE)
      
      pF_1 <- tempfile(fileext = ".bed")
      pF_2 <- tempfile(fileext = ".bed")
      
      BrDU_ColHeads <- "\"chrom\\tBWpStart\\tBWpEnd\\tBWpSummit\\tBCpStart\\tBCpEnd\\tBCpSummit\\toriName\\toriCenter\""
      ChIP_ColHeads <- "\"chrom\\tCWpStart\\tCWpEnd\\tCWpSummit\\tCCpStart\\tCCpEnd\\tCCpSummit\\toriName\\toriCenter\""
      eSPAN_ColHeads <- "\"chrom\\tEWpStart\\tEWpEnd\\tEWpSummit\\tECpStart\\tECpEnd\\tECpSummit\\toriName\\toriCenter\""
      
      #Overlaps between both watson and crick strands
      
      #BrDU
      write.table(BrDU_Input$watsonPeaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      write.table(BrDU_Input$crickPeaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      Overlapping_BrDU_Peaks <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$11,$12,$14,$15,$18}'", 
                                                        pF_1, pF_2, BrDU_ColHeads)), header = TRUE ) 
      Overlapping_BrDU_Peaks <- Overlapping_BrDU_Peaks[!duplicated(Overlapping_BrDU_Peaks$oriName), ]
      
      #ChIP
      write.table(ChIP_Input$watsonPeaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      write.table(ChIP_Input$crickPeaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      Overlapping_ChIP_Peaks <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$11,$12,$14,$15,$18}'", 
                                                        pF_1, pF_2, ChIP_ColHeads)), header = TRUE ) 
      Overlapping_ChIP_Peaks <- Overlapping_ChIP_Peaks[!duplicated(Overlapping_ChIP_Peaks$oriName), ]
      
      #eSPAN
      write.table(eSPAN_ChIP$watsonPeaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      write.table(eSPAN_ChIP$crickPeaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      Overlapping_eSPAN_Peaks <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$11,$12,$14,$15,$18}'", 
                                                         pF_1, pF_2, eSPAN_ColHeads)), header = TRUE ) 
      Overlapping_eSPAN_Peaks <- Overlapping_eSPAN_Peaks[!duplicated(Overlapping_eSPAN_Peaks$oriName), ]
      
      #Overlaps between both BrDU and ChIP samples
      
      BrDU_ChIP_ColHeads <- paste0("\"chrom\\tCWpStart\\tCWpEnd\\tCWpSummit\\tCCpStart\\tCCpEnd\\tCCpSummit\\tBWpStart\\tBWpEnd\\tBWpSummit\\",
                                   "tBCpStart\\tBCpEnd\\tBCpSummit\\toriName\\toriCenter\"")
      
      write.table(Overlapping_ChIP_Peaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      write.table(Overlapping_BrDU_Peaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      BrDU_ChIP_Overlapping_Peaks <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$6,$7,$11,$12,$13,$14,$15,$16,$17,$18}'", 
                                                             pF_1, pF_2, BrDU_ChIP_ColHeads)), header = TRUE ) 
      BrDU_ChIP_Overlapping_Peaks <- BrDU_ChIP_Overlapping_Peaks[!duplicated(BrDU_ChIP_Overlapping_Peaks$oriName), ]
      
      
      #Overlaps between both BrDU, ChIP and eSPAN samples
      
      BrDU_ChIP_eSPAN_ColHeads <- paste0("\"chrom\\tCWpStart\\tCWpEnd\\tCWpSummit\\tCCpStart\\tCCpEnd\\tCCpSummit\\",
                                         "tBWpStart\\tBWpEnd\\tBWpSummit\\tBCpStart\\tBCpEnd\\tBCpSummit\\",
                                         "tEWpStart\\tEWpEnd\\tEWpSummit\\tECpStart\\tECpEnd\\tECpSummit\\toriName\\toriCenter\"")
      
      write.table(BrDU_ChIP_Overlapping_Peaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      write.table(Overlapping_eSPAN_Peaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
      
      BrDU_ChIP_eSPAN_Peaks <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$17,$18,$19,$20,$21,$22,$23,$24}'", 
                                                       pF_1, pF_2, BrDU_ChIP_eSPAN_ColHeads)), header = TRUE ) 
      BrDU_ChIP_eSPAN_Peaks <- BrDU_ChIP_eSPAN_Peaks[!duplicated(BrDU_ChIP_eSPAN_Peaks$oriName), ]
      
      unlink(c(pF_1, pF_2))
      
      #add name column to peakLists
      if(dim(Overlapping_BrDU_Peaks)[1]>0){
        Overlapping_BrDU_Peaks$name <- "BrDU"
      }
      
      if(dim(BrDU_ChIP_Overlapping_Peaks)[1]>0){
        BrDU_ChIP_Overlapping_Peaks$name <- "BrDU_ChIP"
      }
      
      if(dim(BrDU_ChIP_eSPAN_Peaks)[1]>0){
        BrDU_ChIP_eSPAN_Peaks$name <- "BrDU_ChIP_eSPAN"
      }
      
      ###Calculate Fork Locations from BrDU data
      SignalRanges <- function(PeakList){
        
        
        Chroms <- paste0("chr", as.roman(1:16))
        
        ForkPos <- NULL
        
        
        for(i in 1:16){
          
          PeakPos <- NULL
          
          i <- i
          
          Chr_Peaks <- PeakList[PeakList$chrom == Chroms[i], ]
          
          
          if(length(Chr_Peaks$chrom)==0) next 
          
          
          for(y in 1:length(Chr_Peaks$chrom)){
            
            y <- y
            
            if(Chr_Peaks$name[1]=="BrDU"){
              WatsonLeft <- round((Chr_Peaks$BWpSummit[y] - Chr_Peaks$BWpStart[y]))
              WatsonRight <- round((Chr_Peaks$BWpEnd[y] - Chr_Peaks$BWpSummit[y]))
              
              CrickLeft <- round((Chr_Peaks$BCpSummit[y] - Chr_Peaks$BCpStart[y]))
              CrickRight <- round((Chr_Peaks$BCpEnd[y] - Chr_Peaks$BCpSummit[y]))
              
              Peaks <- cbind.data.frame(WatsonLeft=WatsonLeft, WatsonRight=WatsonRight, 
                                        CrickLeft=CrickLeft, CrickRight=CrickRight, oriName=Chr_Peaks$oriName[y])
            }
            
            if(Chr_Peaks$name[1]=="BrDU_ChIP"){
              WatsonLeft <- round((Chr_Peaks$BWpSummit[y] - Chr_Peaks$CWpStart[y]))
              WatsonRight <- round((Chr_Peaks$CWpEnd[y] - Chr_Peaks$BWpSummit[y]))
              
              CrickLeft <- round((Chr_Peaks$BCpSummit[y] - Chr_Peaks$CCpStart[y]))
              CrickRight <- round((Chr_Peaks$CCpEnd[y] - Chr_Peaks$BCpSummit[y]))
              
              Peaks <- cbind.data.frame(WatsonLeft=WatsonLeft, WatsonRight=WatsonRight, 
                                        CrickLeft=CrickLeft, CrickRight=CrickRight, oriName=Chr_Peaks$oriName[y])
            }
            
            if(Chr_Peaks$name[1]=="BrDU_ChIP_eSPAN"){
              WatsonLeft <- round((Chr_Peaks$BWpSummit[y] - Chr_Peaks$EWpStart[y]))
              WatsonRight <- round((Chr_Peaks$EWpEnd[y] - Chr_Peaks$BWpSummit[y]))
              
              CrickLeft <- round((Chr_Peaks$BCpSummit[y] - Chr_Peaks$ECpStart[y]))
              CrickRight <- round((Chr_Peaks$ECpEnd[y] - Chr_Peaks$BCpSummit[y]))
              
              Peaks <- cbind.data.frame(WatsonLeft=WatsonLeft, WatsonRight=WatsonRight, 
                                        CrickLeft=CrickLeft, CrickRight=CrickRight, oriName=Chr_Peaks$oriName[y])
            }
            
            PeakPos <- rbind.data.frame(PeakPos, Peaks)
            PeakPos[PeakPos < 0] <- 0
            
          }
          
          ForkPos <- rbind.data.frame(ForkPos, PeakPos)
          
        }
        return(ForkPos)
      }
      
      #BrDU
      if(dim(Overlapping_BrDU_Peaks)[1]>0){
        ##Lagging and Leading strand lengths synthesised 
        ForkPos <- SignalRanges(Overlapping_BrDU_Peaks)
        
        LaggLeadSynthesis <- cbind.data.frame(lagging=ForkPos$CrickLeft+ForkPos$WatsonRight, 
                                              leading=ForkPos$WatsonLeft+ForkPos$CrickRight, 
                                              oriName=ForkPos$oriName)
        
        
        LeadingAverage <- round(mean(LaggLeadSynthesis$leading))
        LaggingAverage <- round(mean(LaggLeadSynthesis$lagging))
        LeadingSd <- round(sd(LaggLeadSynthesis$leading))
        LaggingSd <- round(sd(LaggLeadSynthesis$lagging))
        
        ###Estimate Averaging Window
        
        #Remove extreme outliers from the data
        Left <- c(ForkPos$WatsonLeft, ForkPos$CrickLeft)
        Right <- c(ForkPos$WatsonRight, ForkPos$CrickRight)
        
        Q <- quantile(Left, probs=0.99, na.rm = T)
        I <- IQR(Left)
        up  <-  Q + 1.5*I # Upper Range
        NewLeft <- Left[which(Left < up)]
        
        Q <- quantile(Right, probs=0.99, na.rm = T)
        I <- IQR(Right)
        up  <-  Q + 1.5*I # Upper Range
        NewRight <- Right[which(Right < up)]
        
        #Range for Averaging Window
        LeftLim <- max(NewLeft) #Left pan
        RightLim <- max(NewRight) #Right pan
        
        #Round the window
        LeftSide <- round(LeftLim/1000+0.5)*1000
        RightSide <- round(RightLim/1000+0.5)*1000
        
        if(LeftSide==RightSide){
          AveragingWindow <- LeftSide
        } else {
          AveragingWindow <- max(c(LeftSide, RightSide))
        } 
      }
      
      #BrDU_ChIP
      if(dim(BrDU_ChIP_Overlapping_Peaks)[1]>0){
        
        ChIPSignal <- SignalRanges(BrDU_ChIP_Overlapping_Peaks)
        
        LaggLeadChIP <- cbind.data.frame(lagging=ChIPSignal$CrickLeft+ChIPSignal$WatsonRight,
                                         leading=ChIPSignal$WatsonLeft+ChIPSignal$CrickRight,
                                         oriName=ChIPSignal$oriName)
        
        LeadingChIPAverage <- round(mean(LaggLeadChIP$leading))
        LaggingChIPAverage <- round(mean(LaggLeadChIP$lagging))
        LeadingChIPSd <- round(sd(LaggLeadChIP$leading))
        LaggingChIPSd <- round(sd(LaggLeadChIP$lagging))
        
        LeftChIP <- c(ChIPSignal$WatsonLeft, ChIPSignal$CrickLeft)
        RightChIP <- c(ChIPSignal$WatsonRight, ChIPSignal$CrickRight)
        
        Q <- quantile(LeftChIP, probs=0.99, na.rm = T)
        I <- IQR(LeftChIP)
        up  <-  Q + 1.5*I # Upper Range
        NewLeftChIP <- LeftChIP[which(LeftChIP < up)]
        
        Q <- quantile(RightChIP, probs=0.99, na.rm = T)
        I <- IQR(RightChIP)
        up  <-  Q + 1.5*I # Upper Range
        NewRightChIP <- RightChIP[which(RightChIP < up)]
        
        LeftLimChIP <- max(NewLeftChIP) #Left pan
        RightLimChIP <- max(NewRightChIP) #Right pan
      }
      
      #BrDU_ChIP_eSPAN
      if(dim(BrDU_ChIP_eSPAN_Peaks)[1]>0){
        
        eSPANSignal <- SignalRanges(BrDU_ChIP_eSPAN_Peaks)
        
        LaggLeadeSPAN <- cbind.data.frame(lagging=eSPANSignal$CrickLeft+eSPANSignal$WatsonRight,
                                          leading=eSPANSignal$WatsonLeft+eSPANSignal$CrickRight,
                                          oriName=eSPANSignal$oriName)
        
        LeadingeSPANAverage <- round(mean(LaggLeadeSPAN$leading))
        LaggingeSPANAverage <- round(mean(LaggLeadeSPAN$lagging))
        LeadingeSPANSd <- round(sd(LaggLeadeSPAN$leading))
        LaggingeSPANSd <- round(sd(LaggLeadeSPAN$lagging))
        
        LefteSPAN <- c(eSPANSignal$WatsonLeft, eSPANSignal$CrickLeft)
        RighteSPAN <- c(eSPANSignal$WatsonRight, eSPANSignal$CrickRight)
        
        Q <- quantile(LefteSPAN, probs=0.99, na.rm = T)
        I <- IQR(LefteSPAN)
        up  <-  Q + 1.5*I # Upper Range
        NewLefteSPAN <- LefteSPAN[which(LefteSPAN < up)]
        
        Q <- quantile(RighteSPAN, probs=0.99, na.rm = T)
        I <- IQR(RighteSPAN)
        up  <-  Q + 1.5*I # Upper Range
        NewRighteSPAN <- RighteSPAN[which(RighteSPAN < up)]
        
        LeftLimeSPAN <- max(NewLefteSPAN) #Left pan
        RightLimeSPAN <- max(NewRighteSPAN) #Right pan
      }
      
      #
      StatDat <- rbind.data.frame(LeftLim, RightLim, paste0(LeadingAverage, "±", LeadingSd), paste0(LaggingAverage, "±", LaggingSd), 
                                  AveragingWindow, binSize, stepSize)
      
      rownames(StatDat) <- c('Left', 'Right', "LeadSynthesis", "LaggSynthesis", "AveragingWindow", 'bin', 'slide')
      colnames(StatDat) <- " "
      
      ##
      dir.create(paste0("~/Desktop/", Pro_1, "/", "Peaks"), showWarnings = FALSE)
      
      write.table(BrDU_Input$comboPeaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Primary_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(Overlapping_BrDU_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "BrDU_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(Overlapping_ChIP_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "ChIP_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(Overlapping_eSPAN_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "eSPAN_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(BrDU_ChIP_Overlapping_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "BrDU_ChIP_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(BrDU_ChIP_eSPAN_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "BrDU_ChIP_eSPAN_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(StatDat, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Synthesis.bed"), quote=FALSE, sep="\t")
      
    }
    
    ProcessPeaks(bamFiles)
    
    
  }
  
  ###
  
  peakFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", "*", ".bed"))
  watsonFiles <-  Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", "*", "watson.bed"))
  crickFiles <-  Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", "*", "crick.bed"))
  
  ###
  
  #10- Plot the global profile
  
  message("Plotting IP profiles ...")
  
  if(!file.exists(paste0("~/Desktop/", Pro_1, "/", Pro_1, "_global_profiles.pdf"))){
    
    Plot_function <- function(CoverageFile, peakFile, DataType, SamplType){
      
      CoverageFile = CoverageFile
      peakFile = peakFile
      
      scoreVals <- c(CoverageFile$wat.score, CoverageFile$cri.score)
      
      Q <- quantile(scoreVals, probs=c(.01, .99), na.rm = FALSE)
      I <- IQR(scoreVals)
      up  <-  Q[2]+1.5*I # Upper Range  
      low <- Q[1]-1.5*I # Lower Range
      Scores <- scoreVals[which(scoreVals < up & scoreVals > low)]
      thd <- round(mean(Scores)+(12*sd(Scores)))
      Ylim_sc <- c(-thd, thd)
      
      wcVals <- log2(CoverageFile$wat.score/CoverageFile$cri.score)
      wcVals[!is.finite(wcVals)] <- 0
      
      Q <- quantile(wcVals, probs=c(.01, .99), na.rm = FALSE)
      I <- IQR(wcVals)
      up  <-  Q[2]+1.5*I # Upper Range  
      low <- Q[1]-1.5*I # Lower Range
      wc <- wcVals[which(wcVals < up & wcVals > low)]
      wcs <- round(mean(wc)+(3*sd(wc)))
      Ylim_wc <- c(-wcs, wcs)
      
      yAxis_reads <- c(-0.3*thd, -0.6*thd, -0.9*thd, 0*thd, 0.3*thd, 0.6*thd, 0.9*thd)
      yAxis_wc <- c(-0.3*wcs, -0.6*wcs, -0.9*wcs, 0*wcs, 0.3*wcs, 0.6*wcs, 0.9*wcs)
      
      Coverage_chr <- CoverageFile[CoverageFile$chrom == seqnames(Scerevisiae)[k], ]
      All_Ori_chr <- All_Ori[All_Ori$chrom == seqnames(Scerevisiae)[k], ]
      Peaks2Plot_chr <- peakFile[peakFile$chrom == seqnames(Scerevisiae)[k], ]
      
      Coverage <-  Coverage_chr[Coverage_chr$chromStart>=S & Coverage_chr$chromStart<=E, ]
      Ori_chr <- All_Ori_chr[All_Ori_chr$chromStart>=S & All_Ori_chr$chromStart<=E, ]
      PeakReg <- Peaks2Plot_chr[Peaks2Plot_chr$peakStart>=S & Peaks2Plot_chr$peakStart<=E, ]
      
      CovWat <- Coverage$wat.score
      CovCri <- Coverage$cri.score
      
      ###
      steps <- 10
      
      if(stepSize == steps){
        CovWat <- CovWat
        CovCri <- CovCri
      }
      if(stepSize > steps){
        s <- stepSize/steps
        CovWat <- as.vector(sapply(CovWat, function (x) rep(x,s)))
        CovCri <- as.vector(sapply(CovCri, function (x) rep(x,s)))
      }
      if(stepSize < steps){
        s <- round(steps/stepSize)
        CovWat <- round(as.vector(tapply(CovWat, gl(length(CovWat)/s, s), mean)))
        CovCri <- round(as.vector(tapply(CovCri, gl(length(CovCri)/s, s), mean)))
      }
      
      ### 
      
      WCrat <- log2(CovWat/CovCri)
      WCrat[!is.finite(WCrat)] <- 0
      
      Cov.wat <- round(as.numeric(smooth.spline(1:length(CovWat), CovWat)$y), 3)
      Cov.cri <- round(as.numeric(smooth.spline(1:length(CovCri), CovCri)$y), 3)
      
      WC.rat <- round(as.numeric(smooth.spline(1:length(WCrat), WCrat)$y), 3)
      
      if(length(Cov.wat) < (Length_per_Row/steps - 1)){
        Cov.wat <- c(Cov.wat, rep(NA, (Length_per_Row/steps - 1)-length(Cov.wat)) )
      }
      
      if(length(Cov.cri) < (Length_per_Row/steps - 1)){
        Cov.cri <- c(Cov.cri, rep(NA, (Length_per_Row/steps - 1)-length(Cov.cri)) )
      }
      
      if(length(WC.rat) < (Length_per_Row/steps - 1)){
        WC.rat <- c(WC.rat, rep(NA, (Length_per_Row/steps - 1)-length(WC.rat)) )
      }
      
      #plot
      par(mar = c(0,0,0,0))
      suppressWarnings(
        plot(Cov.wat, type='h', ylim=Ylim_sc, col =  rgb(100,0,0,alpha=180, maxColorValue=255), ylab=' ', xlab=' ', xaxt='n', yaxt='n', lwd=0.07, bty = 'n',  cex.lab=1, las = 2, xaxs='i')
      ) 
      lines(Cov.cri*(-1), type='h', col =  rgb(0,100,0,alpha=180, maxColorValue=255), lwd=0.07)
      
      segments(x0 = c(Coverage$chromStart - S)/steps, x1 = c(Coverage$chromEnd - S)/steps, y0 = 0, y1= 0, lwd = 0.5)
      
      #plot left y axis
      axis(side = 2, at = yAxis_reads, labels = round(yAxis_reads), line = 0, tick = TRUE, lwd.ticks = 1.5, las = 2, cex.axis = 0.8)
      abline(h=yAxis_reads, lwd=0.05, col =  rgb(112,128,144,alpha=225, maxColorValue=255))
      
      #draw the peaks
      if(DataType=="strandedRatio"){
        if(length(PeakReg$peakStart) > 0){
          for(i in 1:length(PeakReg$peakStart)){
            segments(x0 = c(PeakReg$peakStart[i] - S)/steps, x1 = c(PeakReg$peakEnd[i] - S)/steps, 
                     y0 = par('usr')[4]-(thd*0.1), y1= par('usr')[4]-(thd*0.1), lwd = 4, col = "red", xpd = TRUE)
            draw.circle(x = (PeakReg$peakSummit[i] - S)/steps, y = par('usr')[4]-(thd*0.1), radius = 50, border = "yellow", lwd=2, col="blue")
          }
          
          if((PeakReg$peakSummit[length(PeakReg$peakStart)] - S)/steps<9000){
            mtext(side=3, line=-0.70, at=9900, adj=1, cex=0.5, "macs2-peaks", col = 'red')
          } 
        }
      }
      
      #draw Replication origins
      
      if(length(Ori_chr$chromStart) > 0){
        for(i in 1:length(Ori_chr$chromStart)){
          draw.circle(x = ((Ori_chr$chromStart[i]+Ori_chr$chromEnd[i])/2 - S)/steps, y = 0, radius = 60, border = "purple", lwd=2, col="yellow")
        }
      }
      
      #put chromosome name
      if(SamplType=="Input" & DataType=="readCoverage"){
        title(main = paste("Chromosome", gsub("[[:punct:]]*chr[[:punct:]]*", "", seqnames(Scerevisiae)[k])), col="gray", adj = 0, cex.main=1.5, line = 0, outer = TRUE)
      }
      
      #put sample name
      if(DataType=="readCoverage"){
        mtext(side=3, line=0.75, at=-50, adj=1, cex=1, SamplType)
      }
      
      #put the datatype
      if(DataType=="readCoverage"){
        mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, DataType)
      }
      if(DataType=="strandedRatio"){
        if(length(PeakReg$peakStart) == 0){
          mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, DataType)
        } 
        if(length(PeakReg$peakStart) > 0){
          if((PeakReg$peakSummit[1] - S)/steps<1500){
            mtext(side=1, line=-1.25, at=100, adj=0, cex=0.85, DataType)
          } else {
            mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, DataType)
          }
        } 
      }
      
      #put origin names
      
      if(length(Ori_chr$chromStart) > 0){
        
        colors_ars <- rep("gray1", length(Ori_chr$name))
        colors_ars[which(Ori_chr$stat=='early')] <- "red"
        colors_ars[which(Ori_chr$stat=='late')] <- "blue"
        
        Ori_ticks <- round((c(Ori_chr$chromStart+Ori_chr$chromEnd)/2 - S)/steps)
        
        GetDistances <- function(x){
          x <- sort(x)
          Distances <- c()
          for(i in 1:length(x)-1){
            Distances <- c(Distances, x[i + 1] - x[i])
          }
          return(Distances)
        }
        
        Ori_Dists <- GetDistances(Ori_ticks)
        
        DistIndex <- which(Ori_Dists<500)
        
        CloseOriIndices <- unique(sort(c(DistIndex, DistIndex+1)))
        
        if(SamplType=="Input" & DataType=="readCoverage"){
          
          if(length(Ori_ticks)>0){
            
            if(length(CloseOriIndices)==0){
              Ori_ticks_distal <- Ori_ticks
              Ori_name_distal <- Ori_chr$name
              color_distal <- colors_ars
            } else {
              Ori_ticks_distal <- Ori_ticks[-CloseOriIndices]
              Ori_name_distal <- Ori_chr$name[-CloseOriIndices]
              color_distal <- colors_ars[-CloseOriIndices]
            }
            
            if(Ori_ticks_distal[1]>300){
              axis(side = 3, at = Ori_ticks_distal, labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              axis(side = 3, at = Ori_ticks_distal, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              
              mtext(rep(intToUtf8(9650), length(Ori_ticks_distal)), at = Ori_ticks_distal, side = 3, line=0, col = "purple", cex = 1)
              mtext(rep(intToUtf8(9650), length(Ori_ticks_distal)), at = Ori_ticks_distal, side = 3, line=0.075, col = "yellow", cex = 0.65)
              
              text(x = Ori_ticks_distal, y = grconvertY(0.925, from = "ndc"), labels = Ori_name_distal, xpd = NA, srt = 0, col = color_distal, cex = 0.9 )
            } else {
              if(length(Ori_ticks_distal)>1){
                axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0, col = "purple", cex = 1)
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0.075, col = "yellow", cex = 0.65)
                
                text(x = Ori_ticks_distal[-1], y = grconvertY(0.925, from = "ndc"), labels = Ori_name_distal[-1], xpd = NA, srt = 0, col = color_distal[-1], cex = 0.9 )
                ###
                axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
                
                text(x = Ori_ticks_distal[1], y = grconvertY(0.925+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
              } else {
                axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
                
                text(x = Ori_ticks_distal[1], y = grconvertY(0.925+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
              }
              
            }
            
          }
          
          if(length(CloseOriIndices)>1){
            Consec_Oris <- split(CloseOriIndices, cumsum(c(1, diff(CloseOriIndices) != 1)))
            
            if(length(Consec_Oris)>=1){
              for(j in 1:length(Consec_Oris)){
                for(s in 0:(length(Consec_Oris[[j]])-1)){
                  
                  Pos <- Ori_ticks[Consec_Oris[[j]]][s+1]
                  Nam <- Ori_chr$name[Consec_Oris[[j]]][s+1]
                  Col <- colors_ars[Consec_Oris[[j]]][s+1]
                  
                  tckL <- seq(0.05, 0.75, length.out = 5)
                  
                  axis(side = 3, at = Pos, labels = F, line = 0, tick = TRUE, tck=-(tckL[s+1]), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  axis(side = 3, at = Pos, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  
                  lineL <- seq(0, 6.20, length.out = 5)
                  
                  mtext(rep(intToUtf8(9650), length(Pos)), at = Pos, side = 3, line=(lineL[s+1])+0, col = "purple", cex = 0.90)
                  mtext(rep(intToUtf8(9650), length(Pos)), at = Pos, side = 3, line=(lineL[s+1])+0.075, col = "yellow", cex = 0.585)
                  
                  texL <- seq(0, 0.0725, length.out = 5)
                  
                  text(x = Pos, y = grconvertY(texL[s+1]+0.9225, from = "ndc"), labels = Nam, xpd = NA, srt = 0, col = Col, cex = 0.80 )
                }
              }
            }  
          }
          
          
        } else {
          
          if(DataType=="readCoverage"){
            axis(side = 3, at = Ori_ticks, labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
            axis(side = 3, at = Ori_ticks, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
            
            mtext(rep(intToUtf8(9650), length(Ori_ticks)), at = Ori_ticks, side = 3, line=0, col = "purple", cex = 1)
            mtext(rep(intToUtf8(9650), length(Ori_ticks)), at = Ori_ticks, side = 3, line=0.075, col = "yellow", cex = 0.65)
          }
        }
        
      }
      
      
      #plot x axis
      if(E>seqlengths(Scerevisiae)[k]){
        axis_ticks <-  seq(0, (seqlengths(Scerevisiae)[k] - S)/steps, ((seqlengths(Scerevisiae)[k] - S)/(steps*(seqlengths(Scerevisiae)[k] - S)/10000))) 
      } else {
        axis_ticks <-  seq(0, (E - S)/steps, ((E - S)/(steps*(Length_per_Row/10000))))
      }
      axis_ticks[1] <- 1
      axis(side = 1, at = axis_ticks, labels = rep(NA, 11), line = 0, tick = TRUE, tck=0.03, lwd.ticks = 1.5)
      axis(side = 1, at = axis_ticks, labels = rep(NA, 11), line = 0, tick = TRUE, tck=-0.03, lwd.ticks = 1.5)
      title(ylab='readDensity', col="gray", cex.lab=1.25, line = 0, outer = T)
      
      if(DataType=="strandedRatio" & SamplType=="eSPAN"){
        axis(side = 1, at = axis_ticks, labels = round((S/1000)+axis_ticks*(steps/1000)), line = 0, tick = F)
        title(xlab="Chromosomal Coordinates (Kbp)", col="gray", cex.lab=1.25, line = 0, outer = T)
      }
      
      #plot the watson to crick ratio in a new plot   
      par(new=TRUE)
      plot(WC.rat, ylim=Ylim_wc, col =  rgb(0,0,205,alpha=90, maxColorValue=255), ann=FALSE, axes=FALSE, type='l', lty = 1, lwd = 1.5, xaxs='i')
      mtext("log2(watson/crick)", side=4, line=0, outer = T, cex = 0.75)
      
      #plot right y axis
      axis(side = 4, at = yAxis_wc, labels = round(yAxis_wc), line = 0, tick = TRUE, lwd.ticks = 1.5, las = 2, cex.axis = 0.8)
      
      #
      box("figure", col="forestgreen") 
    }
    
    BrDU_watson <- read.table(watsonFiles[1], header = T)
    BrDU_crick <- read.table(crickFiles[1], header = T)
    ChIP_watson <- read.table(watsonFiles[2], header = T)
    ChIP_crick <- read.table(crickFiles[2], header = T)
    Input_watson <- read.table(watsonFiles[3], header = T)
    Input_crick <- read.table(crickFiles[3], header = T)
    eSPAN_B_watson <- read.table(watsonFiles[2], header = T)
    eSPAN_B_crick <- read.table(crickFiles[2], header = T)
    
    Input_coverage <- cbind.data.frame(chrom=Input_watson$chrom, chromStart=Input_watson$chromStart, chromEnd=Input_watson$chromEnd, wat.score=Input_watson$ip.score, cri.score=Input_crick$ip.score)
    ChIP_coverage <- cbind.data.frame(chrom=ChIP_watson$chrom, chromStart=ChIP_watson$chromStart, chromEnd=ChIP_watson$chromEnd, wat.score=ChIP_watson$ip.score, cri.score=ChIP_crick$ip.score)
    BrDU_coverage <- cbind.data.frame(chrom=BrDU_watson$chrom, chromStart=BrDU_watson$chromStart, chromEnd=BrDU_watson$chromEnd, wat.score=BrDU_watson$ip.score, cri.score=BrDU_crick$ip.score)
    eSPAN_coverage <- cbind.data.frame(chrom=eSPAN_B_watson$chrom, chromStart=eSPAN_B_watson$chromStart, chromEnd=eSPAN_B_watson$chromEnd, wat.score=eSPAN_B_watson$ip.score, cri.score=eSPAN_B_crick$ip.score)
    
    ChIP_ratio <- cbind.data.frame(chrom=ChIP_watson$chrom, chromStart=ChIP_watson$chromStart, chromEnd=ChIP_watson$chromEnd, wat.score=ChIP_watson$ratio, cri.score=ChIP_crick$ratio)
    BrDU_ratio <- cbind.data.frame(chrom=BrDU_watson$chrom, chromStart=BrDU_watson$chromStart, chromEnd=BrDU_watson$chromEnd, wat.score=BrDU_watson$ratio, cri.score=BrDU_crick$ratio)
    eSPAN_ratio <- cbind.data.frame(chrom=eSPAN_B_watson$chrom, chromStart=eSPAN_B_watson$chromStart, chromEnd=eSPAN_B_watson$chromEnd, wat.score=eSPAN_B_watson$ratio, cri.score=eSPAN_B_crick$ratio)
    
    BrDU_Peaks <- read.table(peakFiles[3], header = T)
    ChIP_Peaks <- read.table(peakFiles[4], header = T)
    eSPAN_Peaks <- read.table(peakFiles[7], header = T)
    
    raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_global_profiles.pdf"), width = 9, height = 11, units = "in", res = 400)
    par(oma=c(2,2,2,2))
    for(k in 1:16){
      
      #define plot layout for global profiles
      LayOut.Dims <- function(x, y){
        xx <- c()
        y <- y
        for(i in 1:length(x)){
          xx <- c(xx, rep(x[i], y))
        }
        return(xx)
      }
      Plot.Nums <- c(1:7)
      Cols <- 19
      
      mat <- matrix(LayOut.Dims(Plot.Nums, Cols),
                    length(Plot.Nums),Cols,byrow=TRUE)
      vec <- c(1,3,6,9,12)
      new_mat <- matrix(0,nrow=length(Plot.Nums)+5,ncol=Cols)
      new_mat[-vec,] <- mat  
      
      fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(2,3,1,3,3,1,3,3,1,3,3,1)), ]
      #fin_mat <- new_mat[rep(1:ncol(fin_mat), times = c(1,1,1,1,1,5)), ]
      
      layout(fin_mat, c(1,1), c(1,1), TRUE)
      #layout.show(n)
      
      
      k <- k
      
      #define plot intervals by chromosomes
      Length_per_Row <- 100000
      Plotting_Rows <- c(seq(0, seqlengths(Scerevisiae)[k]+Length_per_Row, Length_per_Row))
      Starts <- Plotting_Rows[-length(Plotting_Rows)]
      Ends <- Plotting_Rows[-1]
      
      #plot global profiles
      for(i in 1:length(Plotting_Rows[-length(Plotting_Rows)])){
        
        steps <- stepSize
        i <- i
        S <- Starts[i]
        E <- Ends[i]
        
        Plot_function(CoverageFile = Input_coverage, peakFile = ChIP_Peaks, DataType="readCoverage", SamplType="Input")
        
        Plot_function(CoverageFile = ChIP_coverage, peakFile = ChIP_Peaks, DataType="readCoverage", SamplType="ChIP")
        Plot_function(CoverageFile = ChIP_ratio, peakFile = ChIP_Peaks, DataType="strandedRatio", SamplType="ChIP")
        
        Plot_function(CoverageFile = BrDU_coverage, peakFile = BrDU_Peaks, DataType="readCoverage", SamplType="BrDU")
        Plot_function(CoverageFile = BrDU_ratio, peakFile = BrDU_Peaks, DataType="strandedRatio", SamplType="BrDU")
        
        Plot_function(CoverageFile = eSPAN_coverage, peakFile = eSPAN_Peaks, DataType="readCoverage", SamplType="eSPAN")
        Plot_function(CoverageFile = eSPAN_ratio, peakFile = eSPAN_Peaks, DataType="strandedRatio", SamplType="eSPAN")
        
        
      }
    }
    dev.off()
    
  }
  
  #####
  
  message("Plotting Results ...")
  
  # plot average enrichment profile around the replication origins
  
  if(PlotAveProfile == "all"){
    
    Average_Enrichment <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = TRUE){
      
      Window <- as.numeric(read.table(peakFiles[6], row.names = 1)["AveragingWindow", ])
      BrDU_Input_Peaks <- read.table(peakFiles[5], header = T)
      
      watsonRatio <- read.table(watsonRatio, header = T)
      crickRatio <- read.table(crickRatio, header = T)
      PeakList <- read.table(PeakList, header = T)
      
      
      if(Normalise2Input == TRUE){
        V <- 7
      } else {
        V <- 5
      } 
      
      ##
      
      
      RemList <- anti_join(BrDU_Input_Peaks, PeakList, by = c("chrom", "oriName", "oriCenter"))
      NamList <- anti_join(BrDU_Input_Peaks, RemList, by = c("chrom", "oriName", "oriCenter"))
      
      NamList <- NamList[!duplicated(NamList$oriName), ]
      
      PeakList$BrDUSummit <- NamList$peakSummit
      ###
      
      PeakList$AvBstart <- PeakList$BrDUSummit - Window
      PeakList$AvBend <- PeakList$BrDUSummit + Window
      
      chrS <- paste0("chr", as.roman(1:16))
      
      #watson
      IP_R <- watsonRatio
      IP_W <- NULL
      for(i in 1:length(chrS)){
        
        OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
        IP_ROs <- IP_R[IP_R$chrom == chrS[i], ]
        
        if(length(OriList_ROs$chrom)==0) next 
        
        IP_C <- NULL
        for(y in 1:length(OriList_ROs$chrom)){
          
          IP_Z <- IP_ROs[IP_ROs$chromStart>=OriList_ROs$AvBstart[y] & IP_ROs$chromStart<=OriList_ROs$AvBend[y], ]
          
          if(length(IP_Z[,V])==round(2*Window/stepSize)){
            Ratios <- IP_Z[,V]
          } 
          
          if(length(IP_Z[,V]) < round(2*Window/stepSize)){
            if(length(1:(length(IP_Z[,V])/2)) < round(2*Window/stepSize)/2){
              Ratios <- c(rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)), IP_Z[,V] )
            }
            if(length((length(IP_Z[,V])/2+1):length(IP_Z[,V])) < round(2*Window/stepSize)/2){
              Ratios <- c(IP_Z[,V], rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)) )
            }
          } 
          
          if(length(IP_Z[,V]) > round(2*Window/stepSize)){
            Ratios <- IP_Z[,V][-c((round(2*Window/stepSize)+1):length(IP_Z[,V]))]
          }
          
          Rat <- matrix(Ratios, ncol = 1)
          
          IP_C <- cbind(IP_C, Rat)
        }
        IP_W <- cbind(IP_W, IP_C)
      }
      IP_W <- as.data.frame(IP_W)
      colnames(IP_W) <- c(1:length(PeakList$chrom))
      
      #crick
      IP_R <- crickRatio
      IP_Cr <- NULL
      for(i in 1:length(chrS)){
        OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
        IP_ROs <- IP_R[IP_R$chrom == chrS[i], ]
        
        if(length(OriList_ROs$chrom)==0) next 
        
        IP_C <- NULL
        for(y in 1:length(OriList_ROs$chrom)){
          
          IP_Z <- IP_ROs[IP_ROs$chromStart>=OriList_ROs$AvBstart[y] & IP_ROs$chromStart<=OriList_ROs$AvBend[y], ]
          
          if(length(IP_Z[,V])==round(2*Window/stepSize)){
            Ratios <- IP_Z[,V]
          } 
          
          if(length(IP_Z[,V]) < round(2*Window/stepSize)){
            if(length(1:(length(IP_Z[,V])/2)) < round(2*Window/stepSize)/2){
              Ratios <- c(rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)), IP_Z[,V] )
            }
            if(length((length(IP_Z[,V])/2+1):length(IP_Z[,V])) < round(2*Window/stepSize)/2){
              Ratios <- c(IP_Z[,V], rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)) )
            }
          } 
          
          if(length(IP_Z[,V]) > round(2*Window/stepSize)){
            Ratios <- IP_Z[,V][-c((round(2*Window/stepSize)+1):length(IP_Z[,V]))]
          }
          
          Rat <- matrix(Ratios, ncol = 1)
          
          IP_C <- cbind(IP_C, Rat)
        }
        IP_Cr <- cbind(IP_Cr, IP_C)
      }
      IP_Cr <- as.data.frame(IP_Cr)
      colnames(IP_Cr) <- c(1:length(PeakList$chrom))
      
      ###
      rowStat <- function(DF){
        
        Quantiles <- apply(DF, 1, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
        Ses <- apply(DF, 1, std.error)
        Means <- apply(DF, 1, mean)
        Sds <- apply(DF, 1, sd)
        
        return(list(q.25 = Quantiles["25%",],
                    Median = Quantiles["50%",],
                    q.75 = Quantiles["75%",],
                    Std.err = Ses,
                    Mean = Means,
                    Sd = Sds))
      }
      
      watsonAverage <- cbind.data.frame(watson.q25 = rowStat(IP_W)$q.25, watson.median = rowStat(IP_W)$Median, watson.q75 = rowStat(IP_W)$q.75, watson.se = rowStat(IP_W)$Std.err, watson.mean = rowStat(IP_W)$Mean, watson.sd = rowStat(IP_W)$Sd)
      crickAverage <- cbind.data.frame(crick.q25 = rowStat(IP_Cr)$q.25, crick.median = rowStat(IP_Cr)$Median, crick.q75 = rowStat(IP_Cr)$q.75, crick.se = rowStat(IP_Cr)$Std.err, crick.mean = rowStat(IP_Cr)$Mean, crick.sd = rowStat(IP_Cr)$Sd)
      
      TwoStrands <- cbind.data.frame(watsonAverage, crickAverage)
      
      rownames(TwoStrands) <- paste0("bin", 1:length(TwoStrands[,1]))
      
      return(TwoStrands)
    }
    WatsonOverCrick_Average <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = TRUE){
      
      Window <- as.numeric(read.table(peakFiles[6], row.names = 1)["AveragingWindow", ])
      BrDU_Input_Peaks <- read.table(peakFiles[5], header = T)
      
      
      Watson <- read.table(watsonRatio, header = T)
      Crick <- read.table(crickRatio, header = T)
      PeakList <- read.table(PeakList, header = T)
      
      
      if(Normalise2Input == TRUE){
        V <- 7
      } else {
        V <- 5
      } 
      
      ##
      
      
      RemList <- anti_join(BrDU_Input_Peaks, PeakList, by = c("chrom", "oriName", "oriCenter"))
      NamList <- anti_join(BrDU_Input_Peaks, RemList, by = c("chrom", "oriName", "oriCenter"))
      
      NamList <- NamList[!duplicated(NamList$oriName), ]
      
      PeakList$BrDUSummit <- NamList$peakSummit
      ###
      
      PeakList$AvBstart <- PeakList$BrDUSummit - Window
      PeakList$AvBend <- PeakList$BrDUSummit + Window
      
      chrS <- paste0("chr", as.roman(1:16))
      
      IP_T <- NULL
      for(i in 1:length(chrS)){
        
        i <- i
        
        OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
        Crick_ROs <- Crick[Crick$chrom == chrS[i], ]
        Watson_ROs <- Watson[Watson$chrom == chrS[i], ]
        
        if(length(OriList_ROs$chrom)==0) next 
        
        IP_C <- NULL
        
        for(y in 1:length(OriList_ROs$chrom)){
          
          y <- y
          
          Crick_S <- Crick_ROs[Crick_ROs$chromStart>=OriList_ROs$AvBstart[y] & Crick_ROs$chromStart<=OriList_ROs$AvBend[y], ]
          Watson_S <- Watson_ROs[Watson_ROs$chromStart>=OriList_ROs$AvBstart[y] & Watson_ROs$chromStart<=OriList_ROs$AvBend[y], ]
          
          IP_Z <- log2(Watson_S[,V] / Crick_S[,V]); IP_Z[!is.finite(IP_Z)] <- 0
          
          if(length(IP_Z)==round(2*Window/stepSize)){
            Ratios <- IP_Z
          } 
          
          if(length(IP_Z) < round(2*Window/stepSize)){
            if(length(1:(length(IP_Z)/2)) < round(2*Window/stepSize)/2){
              Ratios <- c(rep(0, (round(2*Window/stepSize))-length(Crick_S$chromStart)), IP_Z )
            }
            if(length((length(IP_Z)/2+1):length(IP_Z)) < round(2*Window/stepSize)/2){
              Ratios <- c(IP_Z, rep(0, (round(2*Window/stepSize))-length(Crick_S$chromStart)) )
            }
          } 
          
          if(length(IP_Z) > round(2*Window/stepSize)){
            Ratios <- IP_Z[-c((round(2*Window/stepSize)+1):length(IP_Z))]
          }
          
          Rat <- matrix(Ratios, ncol = 1)
          
          IP_C <- cbind(IP_C, Rat)
          
        }
        
        IP_T <- cbind(IP_T, IP_C)
        IP_T <- as.data.frame(IP_T)
      }
      colnames(IP_T) <- c(1:length(PeakList$chrom))
      
      ##extract bin stats
      rowStat <- function(DF){
        
        
        Quantiles <- apply(DF, 1, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
        Ses <- apply(DF, 1, std.error)
        Means <- apply(DF, 1, mean)
        Sds <- apply(DF, 1, sd)
        
        return(list(q.25 = Quantiles["25%",],
                    Median = Quantiles["50%",],
                    q.75 = Quantiles["75%",],
                    Std.err = Ses,
                    Mean = Means,
                    Sd = Sds))
      }
      
      watsonOvercrickAvg <- cbind.data.frame(q25 = rowStat(IP_T)$q.25, median = rowStat(IP_T)$Median, q75 = rowStat(IP_T)$q.75, se = rowStat(IP_T)$Std.err, mean = rowStat(IP_T)$Mean, sd = rowStat(IP_T)$Sd)
      
      rownames(watsonOvercrickAvg) <- paste0("bin", 1:length(watsonOvercrickAvg[,1]))
      
      return(watsonOvercrickAvg)
      
    }
    Bias_at_individual_peaks <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = F){
      
      watsonRatio <- read.table(watsonRatio, header = T)
      crickRatio <- read.table(crickRatio, header = T)
      PeakList <- read.table(PeakList, header = T)
      
      
      chrS <- paste0("chr", as.roman(1:16))
      
      IP_T <- NULL
      
      for(i in 1:length(chrS)){
        
        i <- i
        
        IP_ROs_pos <- watsonRatio[watsonRatio$chrom == chrS[i], ]
        IP_ROs_neg <- crickRatio[crickRatio$chrom == chrS[i], ]
        Chr_Peaks <- PeakList[PeakList$chrom == chrS[i], ]
        
        if(length(Chr_Peaks$chrom)==0) next 
        
        ###
        
        IP_C <- NULL
        
        for(y in 1:length(Chr_Peaks$chrom)){
          
          y <- y
          
          if(PeakList$name[1]=="BrDU"){
            IP_Z_pos_left <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$BWpStart[y] & IP_ROs_pos$chromStart<=Chr_Peaks$BWpSummit[y], ]
            IP_Z_pos_right <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$BWpSummit[y] & IP_ROs_pos$chromStart<=Chr_Peaks$BWpEnd[y], ]
            
            IP_Z_neg_left <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$BCpStart[y] & IP_ROs_neg$chromStart<=Chr_Peaks$BCpSummit[y], ]
            IP_Z_neg_right <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$BCpSummit[y] & IP_ROs_neg$chromStart<=Chr_Peaks$BCpEnd[y], ]
          }
          
          if(PeakList$name[1]=="BrDU_ChIP"){
            #BrDU or Input
            if(IP_ROs_pos$name[1]==paste0(Pro_1, "_BrDU_watson") || IP_ROs_pos$name[1]==paste0(Pro_1, "_Input_watson") ){
              IP_Z_pos_left <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$BWpStart[y] & IP_ROs_pos$chromStart<=Chr_Peaks$BWpSummit[y], ]
              IP_Z_pos_right <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$BWpSummit[y] & IP_ROs_pos$chromStart<=Chr_Peaks$BWpEnd[y], ]
            }
            if(IP_ROs_neg$name[1]==paste0(Pro_1, "_BrDU_crick") || IP_ROs_neg$name[1]==paste0(Pro_1, "_Input_crick") ) {
              IP_Z_neg_left <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$BCpStart[y] & IP_ROs_neg$chromStart<=Chr_Peaks$BCpSummit[y], ]
              IP_Z_neg_right <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$BCpSummit[y] & IP_ROs_neg$chromStart<=Chr_Peaks$BCpEnd[y], ]
            }
            #eSPAN
            if(IP_ROs_pos$name[1]==paste0(Pro_1, "_eSPAN_watson") ){
              IP_Z_pos_left <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$CWpStart[y] & IP_ROs_pos$chromStart<=Chr_Peaks$BWpSummit[y], ]
              IP_Z_pos_right <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$BWpSummit[y] & IP_ROs_pos$chromStart<=Chr_Peaks$CWpEnd[y], ]
            }
            if(IP_ROs_neg$name[1]==paste0(Pro_1, "_eSPAN_crick") ){
              IP_Z_neg_left <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$CCpStart[y] & IP_ROs_neg$chromStart<=Chr_Peaks$BCpSummit[y], ]
              IP_Z_neg_right <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$BCpSummit[y] & IP_ROs_neg$chromStart<=Chr_Peaks$CCpEnd[y], ]
            }
          }
          
          if(PeakList$name[1]=="BrDU_ChIP_eSPAN"){
            #BrDU or Input
            if(IP_ROs_pos$name[1]==paste0(Pro_1, "_BrDU_watson") || IP_ROs_pos$name[1]==paste0(Pro_1, "_Input_watson") ){
              IP_Z_pos_left <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$BWpStart[y] & IP_ROs_pos$chromStart<=Chr_Peaks$BWpSummit[y], ]
              IP_Z_pos_right <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$BWpSummit[y] & IP_ROs_pos$chromStart<=Chr_Peaks$BWpEnd[y], ]
            }
            if(IP_ROs_neg$name[1]==paste0(Pro_1, "_BrDU_crick") || IP_ROs_neg$name[1]==paste0(Pro_1, "_Input_crick")){
              IP_Z_neg_left <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$BCpStart[y] & IP_ROs_neg$chromStart<=Chr_Peaks$BCpSummit[y], ]
              IP_Z_neg_right <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$BCpSummit[y] & IP_ROs_neg$chromStart<=Chr_Peaks$BCpEnd[y], ]
            }
            #eSPAN
            if(IP_ROs_pos$name[1]==paste0(Pro_1, "_eSPAN_watson")){
              IP_Z_pos_left <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$EWpStart[y] & IP_ROs_pos$chromStart<=Chr_Peaks$BWpSummit[y], ]
              IP_Z_pos_right <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$BWpSummit[y] & IP_ROs_pos$chromStart<=Chr_Peaks$EWpEnd[y], ]
            }
            if(IP_ROs_neg$name[1]==paste0(Pro_1, "_eSPAN_crick")){
              IP_Z_neg_left <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$ECpStart[y] & IP_ROs_neg$chromStart<=Chr_Peaks$BCpSummit[y], ]
              IP_Z_neg_right <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$BCpSummit[y] & IP_ROs_neg$chromStart<=Chr_Peaks$ECpEnd[y], ]
            }
          }
          
          
          IP_Z_CL <- sum(IP_Z_pos_left[,5], na.rm = TRUE)
          IP_Z_CR <- sum(IP_Z_pos_right[,5], na.rm = TRUE)
          IP_Z_WL <- sum(IP_Z_neg_left[,5], na.rm = TRUE)
          IP_Z_WR <- sum(IP_Z_neg_right[,5], na.rm = TRUE)
          
          In_Z_CL <- sum(IP_Z_pos_left[,6], na.rm = TRUE)
          In_Z_CR <- sum(IP_Z_pos_right[,6], na.rm = TRUE)
          In_Z_WL <- sum(IP_Z_neg_left[,6], na.rm = TRUE)
          In_Z_WR <- sum(IP_Z_neg_right[,6], na.rm = TRUE)
          
          IP_Z_CL_WR <- IP_Z_CL + IP_Z_WR; IP_Z_CL_WR[is.na(IP_Z_CL_WR)]<-0
          IP_Z_WL_CR <- IP_Z_WL + IP_Z_CR; IP_Z_WL_CR[is.na(IP_Z_WL_CR)]<-0
          
          IP_Z_Pvalue <- binom.test(round(c(IP_Z_CL_WR, IP_Z_WL_CR)))$p.value
          
          In_Z_CL_WR <- In_Z_CL + In_Z_WR; In_Z_CL_WR[is.na(In_Z_CL_WR)]<-0
          In_Z_WL_CR <- In_Z_WL + In_Z_CR; In_Z_WL_CR[is.na(In_Z_WL_CR)]<-0
          
          In_Z_LaLe <- In_Z_CL_WR/In_Z_WL_CR; In_Z_LaLe[!is.finite(In_Z_LaLe)] <- 0
          IP_Z_LaLe <- IP_Z_CL_WR/IP_Z_WL_CR; IP_Z_LaLe[!is.finite(IP_Z_LaLe)] <- 0
          
          if(Normalise2Input==FALSE){
            IP_Z_LagLead <- log2(IP_Z_LaLe); IP_Z_LagLead[!is.finite(IP_Z_LagLead)] <- 0
          } else {
            IP_Z_LagLead <- log2(IP_Z_LaLe/In_Z_LaLe); IP_Z_LagLead[!is.finite(IP_Z_LagLead)] <- 0
          }
          
          IP_Z <- cbind.data.frame(round(IP_Z_CL_WR), round(IP_Z_WL_CR), IP_Z_Pvalue, IP_Z_LagLead)
          IP_C <- rbind.data.frame(IP_C, IP_Z)
        }
        IP_C <- cbind.data.frame(Chr_Peaks$chrom, Chr_Peaks$oriName, IP_C)
        colnames(IP_C) <- c('chrom', 'name', "Lagg.sum", "Lead.sum", "p_value", "Bias")
        
        status <- rep('null', length(Chr_Peaks$chrom))
        signif <- rep('null', length(Chr_Peaks$chrom))
        status[which(IP_C$Bias>=0)] <- 'lagging_bias'; status[which(IP_C$Bias<0)] <- 'leading_bias'
        signif[which(IP_C$p_value<=10e-6)] <- 'significant'; signif[which(IP_C$p_value>10e-6)] <- 'not_signif'
        
        IP_C$description <- paste0(signif, "_", status)
        
        IP_T <- rbind(IP_T, IP_C)
        IP_T <- as.data.frame(IP_T)
      }
      return(IP_T)
    }
    
    ##
    LeftLim <- as.numeric(read.table(peakFiles[6], row.names = 1)["Left", ])
    RightLim <- as.numeric(read.table(peakFiles[6], row.names = 1)["Right", ])
    LeadingAverage <-  read.table(peakFiles[6], row.names = 1)["LeadSynthesis", ]
    LaggingAverage <-  read.table(peakFiles[6], row.names = 1)["LaggSynthesis", ]
    binSize <- as.numeric(read.table(peakFiles[6], row.names = 1)["bin", ])
    stepSize <- as.numeric(read.table(peakFiles[6], row.names = 1)["slide", ])
    AveragingWindow <- as.numeric(read.table(peakFiles[6], row.names = 1)["AveragingWindow", ])
    
    PlotEnrichments <- function(DataFile, PlotHeader){
      
      Watson <- round(as.numeric(smooth.spline(1:length(DataFile$watson.median), DataFile$watson.median)$y), 2)
      Crick <- round(as.numeric(smooth.spline(1:length(DataFile$crick.median), DataFile$crick.median)$y), 2)
      
      Wat25 <- round(as.numeric(smooth.spline(1:length(DataFile$watson.q25), DataFile$watson.q25)$y), 2)
      Wat75 <- round(as.numeric(smooth.spline(1:length(DataFile$watson.q75), DataFile$watson.q75)$y), 2)
      
      Cri25 <- round(as.numeric(smooth.spline(1:length(DataFile$crick.q25), DataFile$crick.q25)$y), 2)
      Cri75 <- round(as.numeric(smooth.spline(1:length(DataFile$crick.q75), DataFile$crick.q75)$y), 2)
      
      Y <- max(round(abs(range(c(Wat75, Cri75*(-1))))+0.5))
      
      plot(Watson,
           ylim = c(-Y, Y),
           main = PlotHeader,
           ylab = "Average Enrichment", cex.main=0.8, xlab = "Distance from BrDU peakSummit (Kbp)", 
           xaxt = "n", col = 'brown3', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs = 'i', yaxs = 'i')
      
      lines(Wat25, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
      lines(Wat75, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
      
      polygon(x = c(1:length(Watson), rev(1:length(Watson))), 
              y = c(Wat25, rev(Wat75)), 
              col = adjustcolor("red", alpha.f = 0.2), border = NA)
      
      lines((Crick)*(-1), lwd=2, col = 'cornflowerblue', type = 'l')
      lines((Cri25)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
      lines((Cri75)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
      
      polygon(x = c(1:length(Crick*(-1)), rev(1:length(Crick*(-1)))), 
              y = c(Cri25*(-1), rev(Cri75*(-1))), 
              col = adjustcolor("cornflowerblue", alpha.f = 0.2), border = NA)
      
      rect(0, -Y, (length(Watson))/2 - (LeftLim/stepSize), Y, col = adjustcolor("grey", alpha.f = 0.65), border = NA)
      rect((length(Watson))/2 + (RightLim/stepSize), -Y, length(Watson), Y, 
           col = adjustcolor("grey", alpha.f = 0.65), border = NA)
      
      text((AveragingWindow)*2-500, Y-1, labels = "Watson", cex = 0.9, col = 'brown3')
      text((AveragingWindow)*2-500, -Y+1, labels = "Crick", cex = 0.9, col = 'cornflowerblue')
      
      abline(h=0,lwd=0.4); abline(v=(length(Watson))/2,lwd=0.4)
      axisLabels <- seq(-AveragingWindow,
                        +AveragingWindow,
                        length.out = 9)
      axisLabels[c(2,4,6,8)] <- NA
      At <- (AveragingWindow/stepSize)*seq(0,2,0.25); At[1] <- 1
      axis(1, at=At, labels = signif(axisLabels/1000, 2))
      
    }
    PlotAverages <- function(DataFile, PlotHeader){
      
      
      Med <- round(as.numeric(smooth.spline(1:length(DataFile$median), DataFile$median)$y), 2)
      q25 <- round(as.numeric(smooth.spline(1:length(DataFile$q25), DataFile$q25)$y), 2)
      q75 <- round(as.numeric(smooth.spline(1:length(DataFile$q75), DataFile$q75)$y), 2)
      
      Y <- max(round(abs(range(q75))+0.5))
      
      plot(Med,
           ylim = c(-Y, +Y),
           main = PlotHeader,
           ylab = "log2 watson/crick", cex.main=0.8, xlab = "Distance from BrDU peakSummit (Kbp)", 
           xaxt = "n", col = 'blue', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs='i', yaxs='i')
      
      polygon(x = c(1:length(Med), rev(1:length(Med))), 
              y = c(q25, rev(q75)), 
              col = adjustcolor("blue", alpha.f = 0.2), border = NA)
      
      rect(0, -Y, (length(Med))/2 - (LeftLim/stepSize), Y, col = adjustcolor("grey", alpha.f = 0.65), border = NA)
      rect((length(Med))/2 + (RightLim/stepSize), -Y, length(Med), Y, 
           col = adjustcolor("grey", alpha.f = 0.65), border = NA)
      
      abline(h=0, lwd=0.4); abline(v=(length(Med))/2,lwd=0.4)
      axisLabels <- seq(-AveragingWindow,
                        +AveragingWindow,
                        length.out = 9)
      axisLabels[c(2,4,6,8)] <- NA
      At <- (AveragingWindow/stepSize)*seq(0,2,0.25); At[1] <- 1
      axis(1, at=At, labels = signif(axisLabels/1000, 2))
      
    }
    PlotIdBias <- function(DataFile, PlotHeader){
      
      DataFile$Bias[!is.finite(DataFile$Bias)] <- 0
      DataFile$Bias[is.na(DataFile$Bias)] <- 0
      
      boxplot(DataFile$Bias, ylim = c(-2,2), las = 2, ylab = 'log2 lagging/leading', 
              cex.main=0.8, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
      
      if(length(DataFile$Bias[which(DataFile$p_value > 10e-6)])>0){
        spreadPoints(values=DataFile$Bias[which(DataFile$p_value > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
      }
      if(length(DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)])>0){
        spreadPoints(values=DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)], position=1.0, pointCex=0.65, col="red", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
      }
      if(length(DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)])>0){
        spreadPoints(values=DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)], position=1.0, pointCex=0.65, col="green", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
      }
      
      
      Lagg <- length(DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)])
      Lead <- length(DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)])
      Indt <- length(DataFile$Bias[which(DataFile$p_value > 10e-6)])
      
      legend("bottomright", legend = c(paste0("lagg", " (", Lagg, ")"), 
                                       paste0("lead", " (", Lead, ")"), 
                                       paste0("inde", " (", Indt, ")")), 
             col = c("red", "green", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
      
      
      
      #extract input file's name
      #library(pryr)
      
      #adr <- vapply(ls(.GlobalEnv), function(x) {
      # a <- as.symbol(x)
      # setNames(eval(bquote(address(.(a))),.GlobalEnv), x)
      # }, FUN.VALUE = character(1))
      
      #InName <- names(adr[adr == address(DataFile)])[2]
      
      InName <- PlotHeader
      
      #calculate p value
      #Decision Tree
      biasP <- binom.test(round(c(Lagg+Lead, Indt)))$p.value
      
      
      if(InName=="eSPAN_BrDU" || InName=="eSPAN"){
        
        if(biasP <= 10e-3 & (Lagg+Lead) > Indt){
          
          biasQ <- binom.test(round(c(Lagg, Lead)))$p.value
          
          if(biasQ <= 10e-6){
            if(Lagg > Lead){
              conc <- paste0("Strong lagging bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
            if(Lead > Lagg){
              conc <- paste0("Strong leading bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          }
          
          if(biasQ > 10e-6 & biasQ <= 10e-4){
            if(Lagg > Lead){
              conc <- paste0("Weak lagging bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
            if(Lead > Lagg){
              conc <- paste0("Weak leading bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          }
          
          if(biasQ > 10e-4  & biasQ <= 10e-2){
            if(Lagg > Lead){
              conc <- paste0("Very weak lagging bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
            if(Lead > Lagg){
              conc <- paste0("Very weak leading bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          }
          
          if(biasQ > 10e-2){
            conc <- paste0("Factor may bind one or the other strand", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
          }
        } 
        else 
        {
          conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
        }
      }
      
      if(InName=="BrDU" || InName=="Input"){
        
        if(biasP <= 10e-3 & (Lagg+Lead) > Indt){
          
          biasQ <- binom.test(round(c(Lagg, Lead)))$p.value
          
          if(biasQ <= 10e-6){
            if(Lagg > Lead){
              conc <- paste0("Strong bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
            if(Lead > Lagg){
              conc <- paste0("Strong bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          }
          
          if(biasQ > 10e-6 & biasQ <= 10e-4){
            if(Lagg > Lead){
              conc <- paste0("Weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
            if(Lead > Lagg){
              conc <- paste0("Weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          }
          
          if(biasQ > 10e-4  & biasQ <= 10e-2){
            if(Lagg > Lead){
              conc <- paste0("Very weak bias for lagging synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
            if(Lead > Lagg){
              conc <- paste0("Very weak bias for leading synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
            }
          }
          
          if(biasQ > 10e-2){
            conc <- paste0("No significant strandedness in DNA synthesis", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
          }
        } 
        else 
        {
          conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
        }
      }
      
      
      mtext(conc, side = 1, line = 2, cex = 0.65)
    }
    
    pdf(paste0("~/Desktop/", Pro_1, "/", Pro_1, "_Results.pdf"), width = 10, height = 12)
    
    par(oma=c(0,0,0,0))
    
    PlotMat <- {matrix(c(     0,0,1,1,1,1,1,1,1,0,0,2,2,2,2,2,2,0,0,0,
                              0,0,1,1,1,1,1,1,1,0,0,2,2,2,2,2,2,0,0,0,
                              0,0,1,1,1,1,1,1,1,0,0,2,2,2,2,2,2,0,0,0,
                              0,0,1,1,1,1,1,1,1,0,0,2,2,2,2,2,2,0,0,0,
                              0,0,1,1,1,1,1,1,1,0,0,2,2,2,2,2,2,0,0,0,
                              
                              3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,
                              3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,
                              3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,
                              3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,
                              3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,
                              
                              7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,
                              7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,
                              7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,
                              7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,
                              7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10,
                              
                              11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,
                              11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,
                              11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,
                              11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,
                              11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,14,14,14,
                              
                              15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
                              15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
                              15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
                              15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15), 
                       
                       24,20,byrow=TRUE)}
    
    if(dim(read.table(peakFiles[3], header = T))[1]>0){
      
      BrDU_Input_AvE_Bp <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[3], Normalise2Input = T)
      ChIP_Input_AvE_Bp <- Average_Enrichment(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = peakFiles[3], Normalise2Input = T)
      eSPAN_BrDU_AvE_Bp <- Average_Enrichment(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[3], Normalise2Input = T)
      eSPAN_AvE_Bp <- Average_Enrichment(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[3], Normalise2Input = F)
      
      BrDU_Input_WoC_Bp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[3], Normalise2Input = T)
      BrDU_WoC_Bp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[3], Normalise2Input = F)
      eSPAN_WoC_Bp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[3], Normalise2Input = F)
      eSPAN_BrDU_WoC_Bp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[3], Normalise2Input = T)
      
      Input_Bp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = peakFiles[3], Normalise2Input = F)
      BrDU_Bp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[3], Normalise2Input = F)
      eSPAN_Bp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[3], Normalise2Input = F)
      eSPAN_BrDU_Bp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[3], Normalise2Input = T)
      
      layout(PlotMat, c(1,1), c(1,1), TRUE)
      plot(image, axes = FALSE)
      plot_venn <- {
        
        x <- list(
          BrDU = read.table(peakFiles[3], header = T)$oriName
        )
        
        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        VennObject <- venn.diagram(x, filename = NULL, main = paste0(Pro_1, " Origin Firing Events"),
                                   main.fontface = "bold", main.cex = 0.6, 
                                   main.pos = c(0.5, 1.25),
                                   category.names = c("BrDU"),
                                   # Circles
                                   lwd = 2,
                                   lty = 'blank',
                                   fill = c("#999999"),
                                   # Numbers
                                   cex = .6,
                                   fontface = "italic",
                                   # Set names
                                   cat.cex = 0.5,
                                   cat.fontface = "bold",
                                   cat.default.pos = "outer",
                                   cat.dist = c(0.1))
        
        # Grid regions of current base plot (ie from frame)
        plot.new()  
        vps <- baseViewports()
        pushViewport(vps$inner, vps$figure, vps$plot)
        grid.draw(VennObject)
        popViewport(3)
        
        
      }
      PlotEnrichments(BrDU_Input_AvE_Bp, "BrDU_Input"); PlotEnrichments(ChIP_Input_AvE_Bp, "ChIP_Input"); PlotEnrichments(eSPAN_AvE_Bp, "eSPAN"); PlotEnrichments(eSPAN_BrDU_AvE_Bp, "eSPAN_BrDU")
      PlotAverages(BrDU_Input_WoC_Bp, "BrDU_Input"); PlotAverages(BrDU_WoC_Bp, "BrDU"); PlotAverages(eSPAN_WoC_Bp, "eSPAN"); PlotAverages(eSPAN_BrDU_WoC_Bp, "eSPAN_BrDU")
      PlotIdBias(Input_Bp, "Input"); PlotIdBias(BrDU_Bp, "BrDU"); PlotIdBias(eSPAN_Bp, "eSPAN"); PlotIdBias(eSPAN_BrDU_Bp, "eSPAN_BrDU")
      
      plot(NULL, xlim=c(0,2), ylim=c(0,2), ylab=" ", xlab=" ", yaxt="n", xaxt="n")
      txt <- paste0("99% BrDU Signal are within ", LeftLim, " (left) and ", RightLim, " (right)", " bps from stranded peakSummit", "\n",
                    "Leading synthesis = ", LeadingAverage, " bps;", " Lagging synthesis = ", LaggingAverage, " bps", "\n",
                    "OriginList = All BrDU Peaks", ", ", "binSize = ", binSize, " bps", ", ", "slide = ", stepSize, " bps")
      text(1, 1, labels = txt, cex = 1)
    }
    
    if(dim(read.table(peakFiles[1], header = T))[1]>0){
      
      BrDU_Input_AvE_BCp <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[1], Normalise2Input = T)
      ChIP_Input_AvE_BCp <- Average_Enrichment(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = peakFiles[1], Normalise2Input = T)
      eSPAN_BrDU_AvE_BCp <- Average_Enrichment(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[1], Normalise2Input = T)
      eSPAN_AvE_BCp <- Average_Enrichment(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[1], Normalise2Input = F)
      
      BrDU_Input_WoC_BCp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[1], Normalise2Input = T)
      BrDU_WoC_BCp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[1], Normalise2Input = F)
      eSPAN_WoC_BCp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[1], Normalise2Input = F)
      eSPAN_BrDU_WoC_BCp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[1], Normalise2Input = T)
      
      Input_BCp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = peakFiles[1], Normalise2Input = F)
      BrDU_BCp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[1], Normalise2Input = F)
      eSPAN_BCp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[1], Normalise2Input = F)
      eSPAN_BrDU_BCp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[1], Normalise2Input = T)
      
      layout(PlotMat, c(1,1), c(1,1), TRUE)
      plot(image, axes = FALSE)
      plot_venn <- {
        
        
        x <- list(
          BrDU = read.table(peakFiles[3], header = T)$oriName, 
          ChIP = read.table(peakFiles[4], header = T)$oriName
        )
        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        VennObject <- venn.diagram(x, filename = NULL, main = paste0(Pro_1, " Origin Firing Events"),
                                   main.fontface = "bold", main.cex = 0.6, 
                                   main.pos = c(0.5, 1.25),
                                   category.names = c("BrDU" , "ChIP"),
                                   # Circles
                                   lwd = 2,
                                   lty = 'blank',
                                   fill = c("#999999", "#E69F00"),
                                   # Numbers
                                   cex = .6,
                                   fontface = "italic",
                                   # Set names
                                   cat.cex = 0.5,
                                   cat.fontface = "bold",
                                   cat.default.pos = "outer",
                                   cat.dist = c(0.1, 0.1))
        
        # Grid regions of current base plot (ie from frame)
        plot.new()  
        vps <- baseViewports()
        pushViewport(vps$inner, vps$figure, vps$plot)
        grid.draw(VennObject)
        popViewport(3)
        
        
      }
      PlotEnrichments(BrDU_Input_AvE_BCp, "BrDU_Input"); PlotEnrichments(ChIP_Input_AvE_BCp, "ChIP_Input"); PlotEnrichments(eSPAN_AvE_BCp, "eSPAN"); PlotEnrichments(eSPAN_BrDU_AvE_BCp, "eSPAN_BrDU")
      PlotAverages(BrDU_WoC_BCp, "BrDU_Input"); PlotAverages(BrDU_WoC_BCp, "BrDU"); PlotAverages(eSPAN_WoC_BCp, "eSPAN"); PlotAverages(eSPAN_BrDU_WoC_BCp, "eSPAN_BrDU")
      PlotIdBias(Input_BCp, "Input"); PlotIdBias(BrDU_BCp, "BrDU"); PlotIdBias(eSPAN_BCp, "eSPAN"); PlotIdBias(eSPAN_BrDU_BCp, "eSPAN_BrDU")
      
    }
    
    if(dim(read.table(peakFiles[2], header = T))[1]>0){
      
      BrDU_Input_AvE_BCEp <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[2], Normalise2Input = T)
      ChIP_Input_AvE_BCEp <- Average_Enrichment(watsonRatio=watsonFiles[2], crickRatio=crickFiles[2], PeakList = peakFiles[2], Normalise2Input = T)
      eSPAN_BrDU_AvE_BCEp <- Average_Enrichment(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[2], Normalise2Input = T)
      eSPAN_AvE_BCEp <- Average_Enrichment(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[2], Normalise2Input = F)
      
      BrDU_Input_WoC_BCEp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[2], Normalise2Input = T)
      BrDU_WoC_BCEp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[2], Normalise2Input = F)
      eSPAN_WoC_BCEp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[2], Normalise2Input = F)
      eSPAN_BrDU_WoC_BCEp <- WatsonOverCrick_Average(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[2], Normalise2Input = T)
      
      Input_BCEp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[3], crickRatio=crickFiles[3], PeakList = peakFiles[2], Normalise2Input = F)
      BrDU_BCEp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[2], Normalise2Input = F)
      eSPAN_BCEp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[2], Normalise2Input = F)
      eSPAN_BrDU_BCEp <- Bias_at_individual_peaks(watsonRatio=watsonFiles[4], crickRatio=crickFiles[4], PeakList = peakFiles[2], Normalise2Input = T)
      
      
      layout(PlotMat, c(1,1), c(1,1), TRUE)
      
      plot(image, axes = FALSE)
      plot_venn <- {
        
        x <- list(
          BrDU = read.table(peakFiles[3], header = T)$oriName, 
          ChIP = read.table(peakFiles[4], header = T)$oriName,
          eSPAN = read.table(peakFiles[7], header = T)$oriName
        )
        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        VennObject <- venn.diagram(x, filename = NULL, main = paste0(Pro_1, " Origin Firing Events"),
                                   main.fontface = "bold", main.cex = 0.6, 
                                   main.pos = c(0.5, 1.25),
                                   category.names = c("BrDU" , "ChIP", "eSPAN"),
                                   # Circles
                                   lwd = 2,
                                   lty = 'blank',
                                   fill = c("#999999", "#E69F00", "#56B4E9"),
                                   # Numbers
                                   cex = .6,
                                   fontface = "italic",
                                   # Set names
                                   cat.cex = 0.5,
                                   cat.fontface = "bold",
                                   cat.default.pos = "outer",
                                   cat.dist = c(0.1, 0.1, 0.1))
        
        # Grid regions of current base plot (ie from frame)
        plot.new()  
        vps <- baseViewports()
        pushViewport(vps$inner, vps$figure, vps$plot)
        grid.draw(VennObject)
        popViewport(3)
        
        
      }
      PlotEnrichments(BrDU_Input_AvE_BCEp, "BrDU_Input"); PlotEnrichments(ChIP_Input_AvE_BCEp, "ChIP_Input"); PlotEnrichments(eSPAN_AvE_BCEp, "eSPAN"); PlotEnrichments(eSPAN_BrDU_AvE_BCEp, "eSPAN_BrDU")
      PlotAverages(BrDU_Input_WoC_BCEp, "BrDU_Input"); PlotAverages(BrDU_WoC_BCEp, "BrDU"); PlotAverages(eSPAN_WoC_BCEp, "eSPAN"); PlotAverages(eSPAN_BrDU_WoC_BCEp, "eSPAN_BrDU")
      PlotIdBias(Input_BCEp, "Input"); PlotIdBias(BrDU_BCEp, "BrDU"); PlotIdBias(eSPAN_BCEp, "eSPAN"); PlotIdBias(eSPAN_BrDU_BCEp, "eSPAN_BrDU")
      
    }
    
    dev.off()
    
    
  }
  
  ###
  
  message("Analysis complete.")
  
  message(paste0("Check results at Desktop folder - ", Pro_1))
  
  
  
  rm(list=ls())
  gc()
  
}