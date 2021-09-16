
ChIPwrapper <- function(Input_R1, 
                        Input_R2, 
                        IP_R1, 
                        IP_R2,
                        
                        ExpTitle = "None",
                        SamplType = "ChIP", 
                        SeqType = "paired",
                        
                        slidingWindow = "YES", 
                        bin = "None",
                        stepSize = "None", 
                        
                        PlotIPprofile = "Global",
                        ChromCoords = "4:50000-150000",
                        
                        peakSet = "RO",
                        
                        AveragingPan = "3000:3000", 
                        Normalise = "YES"
                        ) {
  
  packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
                "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw",
                "BSgenome.Scerevisiae.UCSC.sacCer3")
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  #
  useDef <- function(a,d) ifelse(isTruthy(a), a,d)
  
  Input_R1 = useDef(Input_R1, NULL)
  Input_R2 = useDef(Input_R2, NULL)
  IP_R1 = useDef(IP_R1, NULL)
  IP_R2 = useDef(IP_R2, NULL)

  ExpTitle = useDef(ExpTitle, "None")
  bin = useDef(bin, "None")
  stepSize = useDef(stepSize, "None")
  ChromCoords = useDef(ChromCoords, "4:50000-150000")
  AveragingPan = useDef(AveragingPan, "3000:3000") 

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
    Pro_1 <- paste0(ExpTitle, "-", SamplType)
  }
  
  
  message(paste0("Experiment: ", Pro_1))
  
  suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1)))   #create directory named with the protein in the Desktop
  
  #Quality check of fastqs'
  
  message("Running QC ...")
  
  if(!dir.exists(paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "QR", ".html"))){
    
    fls = c(Input_R1, Input_R2, IP_R1, IP_R2)
    
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
    
    RunAlign_paired <- function(File_R1, File_R2, SampName){
      
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
      
      system(sprintf("(bowtie2-2.4.4-macos-x86_64/bowtie2 -p 16  --no-discordant --fr -x %s -1 %s -2 %s -S %s) 2> %s", 
                     ref_index, File_R1, File_R2, Sam, AlnLog))
      
      system(sprintf("samtools-1.13/samtools view -bS -@ 15 -q 30 -f 2 %s > %s", Sam, Bam))
      
      system(sprintf("samtools-1.13/samtools collate -@ 15 -o %s %s", nmCollate, Bam))
      
      system(sprintf("samtools-1.13/samtools fixmate -@ 15 -m %s %s", nmCollate, fixMat))
      
      system(sprintf("samtools-1.13/samtools sort -l 9 -@ 15 -m 1024M  -O bam -o %s %s", SrtBam, fixMat))
      
      system(sprintf("samtools-1.13/samtools markdup -@ 15 %s %s", SrtBam, SFBam))
      
      system(sprintf("samtools-1.13/samtools index -@ 15 %s", SFBam))
      
      unlink(c(Sam, Bam, nmCollate, fixMat, SrtBam), recursive = T, force = T)
      
    }
    
    RunAlign_single <- function(File_R1, SampName){
      
      tempdir(check = TRUE)
      
      Sam <- tempfile(fileext = ".sam")
      Bam <- tempfile(fileext = ".bam")
      SrtBam <- tempfile(fileext = ".bam")
      
      
      Pro_1 <- Pro_1
      Pro_2 <- SampName
      
      suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Bam")))
      
      AlnLog <- paste0("~/Desktop/", Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".log")
      SrtBam <- paste0("~/Desktop/", Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".bam")
      
      #read the indexed reference genome for the alignment of sequenced data
      ref_index <- "bowtie2-2.4.4-macos-x86_64/indexes/S288C_Ref"
      
      #following commands will run the alignemnt, check quality, sort, filter and index the resultant bam file 
      
      system(sprintf("(bowtie2-2.4.4-macos-x86_64/bowtie2 -p 16  --no-unal -x %s -U %s -S %s) 2> %s", 
                     ref_index, File_R1, Sam, AlnLog))
      
      system(sprintf("samtools-1.13/samtools view -bS -@ 15 -q 30 %s > %s", Sam, Bam))
      
      system(sprintf("samtools-1.13/samtools sort -l 9 -@ 15 -m 1024M  -O bam -o %s %s", SrtBam, Bam))
      
      system(sprintf("samtools-1.13/samtools index -@ 15 %s", SrtBam))
      
      unlink(c(Sam, Bam), recursive = T, force = T)
      
    }
    
    if(SeqType == "paired"){
      
      RunAlign_paired(Input_R1, Input_R2, "Input")
      
      if(SamplType == "ChIP"){
        RunAlign_paired(IP_R1, IP_R2, "ChIP")
      }
      if(SamplType == "BrDU"){
        RunAlign_paired(IP_R1, IP_R2, "BrDU")
      }
      
    }
    
    if(SeqType == "single"){
      
      RunAlign_single(Input_R1, "Input")
      
      if(SamplType == "ChIP"){
        RunAlign_single(IP_R1, "ChIP")
      }
      if(SamplType == "BrDU"){
        RunAlign_single(IP_R1, "BrDU")
      }
      
    }
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
    
    BamCoverage <- function(bamFile, binSize = binSize, stepSize = stepSize, slidingWindow = "YES"){
      
      Pro_1 <- unlist(strsplit(tools::file_path_sans_ext(basename(bamFile)), split='_', fixed=TRUE))[[1]] #extract protein name
      Pro_2 <- unlist(strsplit(tools::file_path_sans_ext(basename(bamFile)), split='_', fixed=TRUE))[[2]] #extract sample name
      
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
      
      suppressWarnings(dir.create(paste0("~/Desktop/", Pro_1, "/", "Coverage"))) 
      
      finFiles_watson <- paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "watson.bed")
      finFiles_crick <- paste0("~/Desktop/", Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "crick.bed")
      
      command_5 <- "bedtools2/bin/bedtools map -a %s -b %s -null 0 -o sum | awk 'BEGIN {OFS=\"\\t\"} {if ($4>=0) print $1,$2,$3,\"%s\",$4}' > %s"
      
      system(sprintf(command_5, binFile, pncFiles_watson, paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"), finFiles_watson))
      system(sprintf(command_5, binFile, pncFiles_crick, paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"), finFiles_crick))
      
      unlink(c(GenomFile, binFile, pncFiles_watson, pncFiles_crick), recursive = T, force = T)
      
    }
    
    for(i in 1:length(bamFiles)){
      BamCoverage(bamFile = bamFiles[i], binSize = binSize, stepSize = stepSize, slidingWindow = slidingWindow)
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
    
    CalculateRatio(CoverageFiles[1], CoverageFiles[3])
    CalculateRatio(CoverageFiles[2], CoverageFiles[4])
    
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
          
          allPeaks <- allPeaks[,c(1,2,3,4,5)]
          
          names(allPeaks) <- c("chrom", 'peakStart', 'peakEnd', 'peakLength', 'peakSummit')
          
          ALLbutROs <- anti_join(allPeaks, Peaks_at_Origins, by = c("chrom", 'peakStart', 'peakEnd'))
          
          PeakList <- list(Peaks_at_Origins, allPeaks, ALLbutROs)
          
          names(PeakList) <- c("ROpeaks", "ALLpeaks", "ALLbutROs")
          
          return(PeakList)
          
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
      
      IP_Input <- PeakFinder(IPBam=bamFiles[1], InBam=bamFiles[2])
      
      #Define overlapping peaks
      
      tempdir(check = TRUE)
      
      pF_1 <- tempfile(fileext = ".bed")
      pF_2 <- tempfile(fileext = ".bed")
      
      #Overlaps between both watson and crick strands
      
        IP_ColHeads <- "\"chrom\\tCWpStart\\tCWpEnd\\tCWpSummit\\tCCpStart\\tCCpEnd\\tCCpSummit\\toriName\\toriCenter\""
        
        write.table(IP_Input$watsonPeaks$ROpeaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        write.table(IP_Input$crickPeaks$ROpeaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        Overlapping_RO_Peaks <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$11,$12,$14,$15,$18}'",
                                                        pF_1, pF_2, IP_ColHeads)), header = TRUE )
        Overlapping_RO_Peaks <- Overlapping_RO_Peaks[!duplicated(Overlapping_RO_Peaks$oriName), ]
        
        Primary_RO_Peaks <- IP_Input$comboPeaks$ROpeaks
        
      
      ##
        IP_ColHeads <- "\"chrom\\tCWpStart\\tCWpEnd\\tCWpSummit\\tCCpStart\\tCCpEnd\\tCCpSummit\""
        
        write.table(IP_Input$watsonPeaks$ALLpeaks, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        write.table(IP_Input$crickPeaks$ALLpeaks, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        Overlapping_ALL_Peaks <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$7,$8,$10}'",
                                                        pF_1, pF_2, IP_ColHeads)), header = TRUE )
        
        Primary_ALL_Peaks <- IP_Input$comboPeaks$ALLpeaks
        
      
        ##
        IP_ColHeads <- "\"chrom\\tCWpStart\\tCWpEnd\\tCWpSummit\\tCCpStart\\tCCpEnd\\tCCpSummit\""
        
        write.table(IP_Input$watsonPeaks$ALLbutROs, file = pF_1, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        write.table(IP_Input$crickPeaks$ALLbutROs, file = pF_2, quote=FALSE, row.names=FALSE, sep="\t", col.names = FALSE)
        Overlapping_ALLminusRO_Peaks <- read.table(pipe(sprintf("bedtools2/bin/bedtools intersect -wa -wb -a %s -b %s | awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$5,$7,$8,$10}'",
                                                        pF_1, pF_2, IP_ColHeads)), header = TRUE )
        
        Primary_ALLminusRO_Peaks <- IP_Input$comboPeaks$ALLbutROs
        
      
      unlink(c(pF_1, pF_2))
      
      #
      
      dir.create(paste0("~/Desktop/", Pro_1, "/", "Peaks"), showWarnings = FALSE)
      
      write.table(Primary_RO_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Primary_RO_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(Overlapping_RO_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Overlapping_RO_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      
      write.table(Primary_ALL_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Primary_ALL_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(Overlapping_ALL_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Overlapping_ALL_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      
      write.table(Primary_ALLminusRO_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Primary_ALLbutOri_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      write.table(Overlapping_ALLminusRO_Peaks, file = paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", Pro_1, "_", "Overlapping_ALLbutOri_Peaks.bed"), quote=FALSE, row.names=FALSE, sep="\t")
      
    }
    
    ProcessPeaks(bamFiles)
    
  }
    
    
  ###
    
    if(peakSet == "RO"){
      peakFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", "*", "_RO_Peaks.bed"))
    }
    
    if(peakSet == "ALL"){
      peakFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", "*", "_ALL_Peaks.bed"))
    }
    
    if(peakSet == "ALL-RO"){
      peakFiles <- Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Peaks", "/", "*", "_ALLbutOri_Peaks.bed"))
    }
  
  watsonFiles <-  Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", "*", "watson.bed"))
  crickFiles <-  Sys.glob(paste0("~/Desktop/", Pro_1, "/", "Ratios", "/", "*", "crick.bed"))
  
  
  ChIP_watson <- read.table(watsonFiles[1], header = T)
  ChIP_crick <- read.table(crickFiles[1], header = T)
  
  ChIP_coverage <- cbind.data.frame(chrom=ChIP_watson$chrom, chromStart=ChIP_watson$chromStart, chromEnd=ChIP_watson$chromEnd, wat.score=ChIP_watson$ip.score, cri.score=ChIP_crick$ip.score)
  ChIP_ratio <- cbind.data.frame(chrom=ChIP_watson$chrom, chromStart=ChIP_watson$chromStart, chromEnd=ChIP_watson$chromEnd, wat.score=ChIP_watson$ratio, cri.score=ChIP_crick$ratio)
  
  ChIP_Peaks <- read.table(peakFiles[2], header = T)
  
  #10- Plot IP profiles
  
  message("Plotting IP profiles ...")
  
  ChIP_Plot_function <- function(CoverageFile, peakFile, DataType, SamplType){
    
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
      plot(Cov.wat, type='h', ylim=Ylim_sc, col =  rgb(100,0,0,alpha=180, maxColorValue=255), ylab=' ', 
           xlab=' ', xaxt='n', yaxt='n', lwd=0.07, bty = 'n',  cex.lab=1, las = 2, xaxs='i')
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
    
    if(SamplType=="ChIP" & DataType=="readCoverage"){
      title(main = paste(Pro_1, " ", SamplType, " Chromosome", gsub("[[:punct:]]*chr[[:punct:]]*", "", seqnames(Scerevisiae)[k])), col="gray", adj = 0, cex.main=1.5, line = 0, outer = TRUE)
    }
    
    #put the datatype
    if(DataType=="readCoverage"){
      mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "read-Coverage")
    }
    if(DataType=="strandedRatio"){
      if(length(PeakReg$peakStart) == 0){
        mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
      } 
      if(length(PeakReg$peakStart) > 0){
        if((PeakReg$peakSummit[1] - S)/steps<1500){
          mtext(side=1, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
        } else {
          mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
        }
      } 
    }
    
    #put origin names
    
    S1 <- Starts[seq(1, length(Starts), 3)]
    S2 <- Starts[seq(2, length(Starts), 3)]
    S3 <- Starts[seq(3, length(Starts), 3)]
    
    if(length(Ori_chr$chromStart) > 0){
      
      if(DataType=="readCoverage"){
        
        Draw_name <- function(Ss, x){
          
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
          
          for(i in 1:length(Ss)){
            if(S == Ss[i]){
              
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
                  
                  text(x = Ori_ticks_distal, y = grconvertY(x, from = "ndc"), labels = Ori_name_distal, xpd = NA, srt = 0, col = color_distal, cex = 0.9 )
                } else {
                  if(length(Ori_ticks_distal)>1){
                    axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    
                    mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0, col = "purple", cex = 1)
                    mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0.075, col = "yellow", cex = 0.65)
                    
                    text(x = Ori_ticks_distal[-1], y = grconvertY(x, from = "ndc"), labels = Ori_name_distal[-1], xpd = NA, srt = 0, col = color_distal[-1], cex = 0.9 )
                    ###
                    axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    
                    mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
                    mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
                    
                    text(x = Ori_ticks_distal[1], y = grconvertY(x+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
                  } else {
                    axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    
                    mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
                    mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
                    
                    text(x = Ori_ticks_distal[1], y = grconvertY(x+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
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
                      
                      text(x = Pos, y = grconvertY(texL[s+1]+x, from = "ndc"), labels = Nam, xpd = NA, srt = 0, col = Col, cex = 0.80 )
                    }
                  }
                }
              }
              
            }
          }
          
        }
        
        Draw_name(S1, 0.92)
        
        Draw_name(S2, 0.615)
        
        Draw_name(S3, 0.31)
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
    
    if(DataType=="strandedRatio" & SamplType=="ChIP"){
      axis(side = 1, at = axis_ticks, labels = round((S/1000)+axis_ticks*(steps/1000)), line = 0, tick = F)
      title(xlab="Chromosomal Coordinates (Kbp)", col="gray", cex.lab=1.25, line = 0, outer = T)
    }
    
    # #plot the watson to crick ratio in a new plot   
    # par(new=TRUE)
    # plot(WC.rat, ylim=Ylim_wc, col =  rgb(0,0,205,alpha=90, maxColorValue=255), ann=FALSE, axes=FALSE, type='l', lty = 1, lwd = 1.5, xaxs='i')
    # mtext("log2(watson/crick)", side=4, line=0, outer = T, cex = 0.75)
    # 
    # #plot right y axis
    # axis(side = 4, at = yAxis_wc, labels = round(yAxis_wc), line = 0, tick = TRUE, lwd.ticks = 1.5, las = 2, cex.axis = 0.8)
    # 
    #
    box("figure", col="forestgreen") 
  }
  
  Local_Profile <- function(CoverageFile, peakFile, DataType, SamplType, ChromCoords = ChromCoords){
    
    Coverage_chr <- CoverageFile[CoverageFile$chrom == seqnames(Scerevisiae)[k], ]
    All_Ori_chr <- All_Ori[All_Ori$chrom == seqnames(Scerevisiae)[k], ]
    Peaks2Plot_chr <- peakFile[peakFile$chrom == seqnames(Scerevisiae)[k], ]
    
    Coverage <-  Coverage_chr[Coverage_chr$chromStart>=S & Coverage_chr$chromStart<=E, ]
    Ori_chr <- All_Ori_chr[All_Ori_chr$chromStart>=S & All_Ori_chr$chromStart<=E, ]
    PeakReg <- Peaks2Plot_chr[Peaks2Plot_chr$peakStart>=S & Peaks2Plot_chr$peakStart<=E, ]
    
    CovWat <- Coverage$wat.score
    CovCri <- Coverage$cri.score
    
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
    
    
    scoreVals <- c(CovWat, CovCri)
    
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
    
    WCrat <- log2(CovWat/CovCri)
    WCrat[!is.finite(WCrat)] <- 0
    
    Cov.wat <- round(as.numeric(smooth.spline(1:length(CovWat), CovWat)$y), 3)
    Cov.cri <- round(as.numeric(smooth.spline(1:length(CovCri), CovCri)$y), 3)
    
    WC.rat <- round(as.numeric(smooth.spline(1:length(WCrat), WCrat)$y), 3)
    
    Length_per_Row <- E - S
    
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
      plot(Cov.wat, type='h', ylim=Ylim_sc, col =  rgb(100,0,0,alpha=180, maxColorValue=255), ylab=' ', 
           xlab=' ', xaxt='n', yaxt='n', lwd=0.07, bty = 'n',  cex.lab=1, las = 2, xaxs='i')
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
    
    if(SamplType=="ChIP" & DataType=="readCoverage"){
      title(main = paste0(Pro_1, "_", SamplType, "_Chromosome", gsub("[[:punct:]]*chr[[:punct:]]*", "", seqnames(Scerevisiae)[k]), ":", S, "-", E), col="gray", adj = 0, cex.main=1.5, line = 0, outer = TRUE)
    }
    
    
    # #put sample name
    # if(DataType=="readCoverage"){
    #   mtext(side=3, line=0.75, at=-50, adj=1, cex=1, SamplType)
    # }
    
    #put the datatype
    if(DataType=="readCoverage"){
      mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "read-Coverage")
    }
    if(DataType=="strandedRatio"){
      if(length(PeakReg$peakStart) == 0){
        mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
      } 
      if(length(PeakReg$peakStart) > 0){
        if((PeakReg$peakSummit[1] - S)/steps<1500){
          mtext(side=1, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
        } else {
          mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
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
      
      
      if(DataType=="readCoverage"){
        
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
            
            text(x = Ori_ticks_distal, y = grconvertY(0.715, from = "ndc"), labels = Ori_name_distal, xpd = NA, srt = 0, col = color_distal, cex = 0.9 )
          } else {
            if(length(Ori_ticks_distal)>1){
              axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              
              mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0, col = "purple", cex = 1)
              mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0.075, col = "yellow", cex = 0.65)
              
              text(x = Ori_ticks_distal[-1], y = grconvertY(0.715, from = "ndc"), labels = Ori_name_distal[-1], xpd = NA, srt = 0, col = color_distal[-1], cex = 0.9 )
              ###
              axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              
              mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
              mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
              
              text(x = Ori_ticks_distal[1], y = grconvertY(0.715+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
            } else {
              axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
              
              mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
              mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
              
              text(x = Ori_ticks_distal[1], y = grconvertY(0.715+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
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
                
                text(x = Pos, y = grconvertY(texL[s+1]+0.715, from = "ndc"), labels = Nam, xpd = NA, srt = 0, col = Col, cex = 0.80 )
              }
            }
          }
        }
        
        
        
        
      } else {
        
        if(DataType=="strandedRatio"){
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
    
    if(DataType=="strandedRatio"){
      axis(side = 1, at = axis_ticks, labels = round((S/1000)+axis_ticks*(steps/1000)), line = 0, tick = F)
      title(xlab="Chromosomal Coordinates (Kbp)", col="gray", cex.lab=1.25, line = 0, outer = T)
    }
    
    # #plot the watson to crick ratio in a new plot   
    # par(new=TRUE)
    # plot(WC.rat, ylim=Ylim_wc, col =  rgb(0,0,205,alpha=90, maxColorValue=255), ann=FALSE, axes=FALSE, type='l', lty = 1, lwd = 1.5, xaxs='i')
    # mtext("log2(watson/crick)", side=4, line=0, outer = T, cex = 0.75)
    # 
    # #plot right y axis
    # axis(side = 4, at = yAxis_wc, labels = round(yAxis_wc), line = 0, tick = TRUE, lwd.ticks = 1.5, las = 2, cex.axis = 0.8)
    # 
    #
    box("figure", col="forestgreen") 
  }
  
  if(PlotIPprofile == "Global"){
    
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
      Plot.Nums <- c(1:6)
      Cols <- 18
      
      mat <- matrix(LayOut.Dims(Plot.Nums, Cols),
                    length(Plot.Nums),Cols,byrow=TRUE)
      
      vec <- c(1,4,7,10)
      new_mat <- matrix(0,nrow=length(Plot.Nums)+4,ncol=Cols)
      new_mat[-vec,] <- mat  
      
      fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(2,3,3,2,3,3,2,3,3,1)), ]
      
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
        
        S <- Starts[i]
        E <- Ends[i]
        
        ChIP_Plot_function(CoverageFile = ChIP_coverage, peakFile = ChIP_Peaks, DataType="readCoverage", SamplType)
        ChIP_Plot_function(CoverageFile = ChIP_ratio, peakFile = ChIP_Peaks, DataType="strandedRatio", SamplType)
        
      }
    }
    dev.off()
    
  }
  
  if(PlotIPprofile == "Local"){
      
      k <- as.numeric(unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[1]])
      SE <- unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[2]]
      S <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[1]])
      E <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[2]])
      
      raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "Chrom", as.roman(k), "_",  as.character(S), "-", as.character(E), ".pdf"), width = 9, height = 11, units = "in", res = 400)
      
      par(oma=c(2,2,2,2))
      LayOut.Dims <- function(x, y){
        xx <- c()
        y <- y
        for(i in 1:length(x)){
          xx <- c(xx, rep(x[i], y))
        }
        return(xx)
      }
      Plot.Nums <- c(1:2); Cols <- 12
      mat <- matrix(LayOut.Dims(Plot.Nums, Cols), length(Plot.Nums),Cols,byrow=TRUE)
      vec <- c(1,3,5)
      new_mat <- matrix(0,nrow=length(Plot.Nums)+3,ncol=Cols)
      new_mat[-vec,] <- mat  
      fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(5,3,1,3,5)), ]
      layout(fin_mat, c(1,1), c(1,1), TRUE)
      
      Local_Profile(CoverageFile = ChIP_coverage, peakFile = ChIP_Peaks, DataType="readCoverage", SamplType, ChromCoords)
      Local_Profile(CoverageFile = ChIP_ratio, peakFile = ChIP_Peaks, DataType="strandedRatio", SamplType, ChromCoords)
      
      dev.off()
      
  }
  
  if(PlotIPprofile == "Both"){
    
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
      Plot.Nums <- c(1:6)
      Cols <- 18
      
      mat <- matrix(LayOut.Dims(Plot.Nums, Cols),
                    length(Plot.Nums),Cols,byrow=TRUE)
      
      vec <- c(1,4,7,10)
      new_mat <- matrix(0,nrow=length(Plot.Nums)+4,ncol=Cols)
      new_mat[-vec,] <- mat  
      
      fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(2,3,3,2,3,3,2,3,3,1)), ]
      
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
        
        S <- Starts[i]
        E <- Ends[i]
        
        ChIP_Plot_function(CoverageFile = ChIP_coverage, peakFile = ChIP_Peaks, DataType="readCoverage", SamplType)
        ChIP_Plot_function(CoverageFile = ChIP_ratio, peakFile = ChIP_Peaks, DataType="strandedRatio", SamplType)
        
      }
    }
    dev.off()
  
    # plot local profile
    
    k <- as.numeric(unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[1]])
    SE <- unlist(strsplit(basename(ChromCoords), split=':', fixed=TRUE))[[2]]
    S <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[1]])
    E <- as.numeric(unlist(strsplit(basename(SE), split='-', fixed=TRUE))[[2]])
    
    raster_pdf(file = paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", "Chrom", as.roman(k), "_",  S, "-", E, ".pdf"), width = 9, height = 11, units = "in", res = 400)
    
    par(oma=c(2,2,2,2))
    LayOut.Dims <- function(x, y){
      xx <- c()
      y <- y
      for(i in 1:length(x)){
        xx <- c(xx, rep(x[i], y))
      }
      return(xx)
    }
    Plot.Nums <- c(1:2); Cols <- 12
    mat <- matrix(LayOut.Dims(Plot.Nums, Cols), length(Plot.Nums),Cols,byrow=TRUE)
    vec <- c(1,3,5)
    new_mat <- matrix(0,nrow=length(Plot.Nums)+3,ncol=Cols)
    new_mat[-vec,] <- mat  
    fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(5,3,1,3,5)), ]
    layout(fin_mat, c(1,1), c(1,1), TRUE)
    
    Local_Profile(CoverageFile = ChIP_coverage, peakFile = ChIP_Peaks, DataType="readCoverage", SamplType, ChromCoords)
    Local_Profile(CoverageFile = ChIP_ratio, peakFile = ChIP_Peaks, DataType="strandedRatio", SamplType, ChromCoords)
    
    dev.off()
    
  }
    
    # plot average enrichment profile around the replication origins
  
  upstream <- as.numeric(unlist(strsplit(basename(AveragingPan), split=':', fixed=TRUE))[[1]])
  downstream <- as.numeric(unlist(strsplit(basename(AveragingPan), split=':', fixed=TRUE))[[2]])
  
  if(upstream == downstream){
    AveragingWindow <- upstream
  } else {
    AveragingWindow <- 3000
  }
  
  
    Average_Enrichment <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = "YES"){
      
      watsonRatio <- read.table(watsonRatio, header = T)
      crickRatio <- read.table(crickRatio, header = T)
      PeakList <- read.table(PeakList, header = T)
      
      Window = AveragingWindow
      
      if(Normalise2Input == "YES"){
        V <- 7
      } else {
        V <- 5
      } 
      
      ###
      
      if(peakSet == "RO"){
        PeakList$AvBstart <- PeakList$oriCenter - Window
        PeakList$AvBend <- PeakList$oriCenter + Window
      } else {
        PeakList$AvBstart <- PeakList$peakSummit - Window
        PeakList$AvBend <- PeakList$peakSummit + Window
      }
      
      
      
      chrS <- paste0("chr", as.roman(1:16))
      
      #watson
      IP_R <- watsonRatio
      IP_W <- NULL
      for(i in 1:length(chrS)){
        
        i <- i
        
        OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
        IP_ROs <- IP_R[IP_R$chrom == chrS[i], ]
        
        if(length(OriList_ROs$chrom)==0) next 
        
        IP_C <- NULL
        for(y in 1:length(OriList_ROs$chrom)){
          
          y <- y
          
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
    WatsonOverCrick_Average <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = "YES"){
      
      Watson <- read.table(watsonRatio, header = T)
      Crick <- read.table(crickRatio, header = T)
      PeakList <- read.table(PeakList, header = T)
      
      Window = AveragingWindow
      
      if(Normalise2Input == "YES"){
        V <- 7
      } else {
        V <- 5
      } 
      
      ###
      
      if(peakSet == "RO"){
        PeakList$AvBstart <- PeakList$oriCenter - Window
        PeakList$AvBend <- PeakList$oriCenter + Window
      } else {
        PeakList$AvBstart <- PeakList$peakSummit - Window
        PeakList$AvBend <- PeakList$peakSummit + Window
      }
      
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
    
    PlotEnrichments <- function(DataFile){
      
      Window = AveragingWindow
      
      Watson <- round(as.numeric(smooth.spline(1:length(DataFile$watson.median), DataFile$watson.median)$y), 2)
      Crick <- round(as.numeric(smooth.spline(1:length(DataFile$crick.median), DataFile$crick.median)$y), 2)
      
      Wat25 <- round(as.numeric(smooth.spline(1:length(DataFile$watson.q25), DataFile$watson.q25)$y), 2)
      Wat75 <- round(as.numeric(smooth.spline(1:length(DataFile$watson.q75), DataFile$watson.q75)$y), 2)
      
      Cri25 <- round(as.numeric(smooth.spline(1:length(DataFile$crick.q25), DataFile$crick.q25)$y), 2)
      Cri75 <- round(as.numeric(smooth.spline(1:length(DataFile$crick.q75), DataFile$crick.q75)$y), 2)
      
      Y <- max(round(abs(range(c(Wat75, Cri75*(-1))))+0.5))
      
      if(peakSet == "RO"){
        txt <- "Distance from RO center (Kbp)"
      } else {
        txt <- "Distance from peak summit (Kbp)"
      }
      
      plot(Watson,
           main = paste0(Pro_1, " Average Enrichment"),
           ylim = c(-Y, Y),
           ylab = "Average Enrichment", xlab = txt, 
           xaxt = "n", col = 'brown3', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs = 'i', yaxs = 'i', 
           cex.main = 1.5, cex.lab = 1.5)
      
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
      
      text(100, Y-0.5, labels = "Watson", cex = 1, col = 'brown3')
      text(100, -Y+0.5, labels = "Crick", cex = 1, col = 'cornflowerblue')
      
      abline(h=0,lwd=0.4); abline(v=(length(Watson))/2,lwd=0.4)
      axisLabels <- seq(-Window,
                        +Window,
                        length.out = 9)
      axisLabels[c(2,4,6,8)] <- NA
      At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
      axis(1, at=At, labels = signif(axisLabels/1000, 2))
      
    }
    PlotAverages <- function(DataFile){
      
      
      Med <- round(as.numeric(smooth.spline(1:length(DataFile$median), DataFile$median)$y), 2)
      q25 <- round(as.numeric(smooth.spline(1:length(DataFile$q25), DataFile$q25)$y), 2)
      q75 <- round(as.numeric(smooth.spline(1:length(DataFile$q75), DataFile$q75)$y), 2)
      
      Y <- max(round(abs(range(q75))+0.5))
      
      if(peakSet == "RO"){
        txt <- "Distance from RO center (Kbp)"
      } else {
        txt <- "Distance from peak summit (Kbp)"
      }
      
      plot(Med,
           ylim = c(-Y, +Y),
           main = paste0(Pro_1, " watson over crick average"),
           ylab = "log2 watson/crick", cex.main = 1.5, cex.lab = 1.5, xlab = txt, 
           xaxt = "n", col = 'blue', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs='i', yaxs='i')
      
      polygon(x = c(1:length(Med), rev(1:length(Med))), 
              y = c(q25, rev(q75)), 
              col = adjustcolor("blue", alpha.f = 0.2), border = NA)
      
      text(100, Y-0.5, labels = "Watson", cex = 1, col = 'brown3')
      text(100, -Y+0.5, labels = "Crick", cex = 1, col = 'cornflowerblue')
      
      abline(h=0, lwd=0.4); abline(v=(length(Med))/2,lwd=0.4)
      axisLabels <- seq(-AveragingWindow,
                        +AveragingWindow,
                        length.out = 9)
      axisLabels[c(2,4,6,8)] <- NA
      At <- (AveragingWindow/stepSize)*seq(0,2,0.25); At[1] <- 1
      axis(1, at=At, labels = signif(axisLabels/1000, 2))
      
    }
  
  

  message("Plotting average profiles ...")
  
  ChIP_Input_AvE <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[2], Normalise2Input = Normalise)
  ChIP_Input_WoC <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = peakFiles[2], Normalise2Input = Normalise)
  
  if(peakSet == "RO"){Ptxt <- "ROpeaks"}
  if(peakSet == "ALL"){Ptxt <- "ALLpeaks"}
  if(peakSet == "ALL-RO"){Ptxt <- "ALL-ROpeaks"}
    
    pdf(paste0("~/Desktop/", Pro_1, "/", Pro_1, "_", Ptxt, "_Profile.pdf"), width = 7, height = 6)
    PlotMat <- {matrix(c(     0,0,1,1,1,1,1,1,1,0,0,
                              0,0,1,1,1,1,1,1,1,0,0,
                              0,0,1,1,1,1,1,1,1,0,0,
                              0,0,1,1,1,1,1,1,1,0,0,
                              0,0,1,1,1,1,1,1,1,0,0,
                              0,0,1,1,1,1,1,1,1,0,0
                              
    ), 
    
    6,11,byrow=TRUE)}
    layout(PlotMat, c(1,1), c(1,1), TRUE)
    PlotEnrichments(ChIP_Input_AvE)
    PlotAverages(ChIP_Input_WoC)
    dev.off()
    

    
    
  message("Analysis complete.")
  
  message(paste0("Check results at Desktop folder - ", Pro_1))
  
  rm(list=ls())
  gc()
  
}


