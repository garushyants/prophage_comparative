#!/usr/bin/env Rscript

library(ggplot2)
#library(ggpubr)
#library(ggrepel)
library(gggenes)
library(ggnewscale)
library(dplyr)
library(tidyr)
library(zoo)
library(stringr)
library(optparse)

###Do input arguments
option_list = list(
  # make_option(c("-d", "--deviationincoverage"), type="double", default=0.3,
  #             help="Coverage percentile cutoff to search for regions of interest [default= %default]", metavar="character"),
  make_option(c("-u", "--minimumderivative"), type="double", default=0.2,
              help="Coverage percentile cutoff to search for regions of interest [default= %default]", metavar="character"),
  make_option(c("-s", "--skipshort"), type="integer", default=200, 
              help="Maximum length of short regions to skip when finding ends on nucleotide level [default= %default]", metavar="character"),
  make_option(c("-g", "--genomeid"), type="character", default="ecor61",#NULL, 
              help="file with genome coverage by BLAST hits", metavar="character"),
  make_option(c("-b", "--blast-file"), type="character", default=NULL,
              help="file with genome coverage by BLAST hits", metavar="character"),
  make_option(c("-l", "--flanks"), type="integer", default=20000,
              help="The length of flanks (nt) to add to the genomad output for blast vizualization [default= %default]", metavar="character"),
  make_option(c("-w", "--windowSizeLarge"), type="integer", default=4000,
              help="Window size for prophage borders prediction from BLAST output in the first pass [default= %default]", metavar="character"),
  make_option(c("-m", "--windowSizeSmall"), type="integer", default=500,
              help="Window size for prophage borders prediction from BLAST output in the second pass [default= %default]", metavar="character"),
  make_option(c("-e", "--genomad-folder"), type="character", default=NULL,
              help="Path to the folder with geNomad results", metavar="character"),
  make_option(c("-c", "--checkv-folder"), type="character", default=NULL,
              help="Path to the folder with CheckV output", metavar="character"),
  make_option(c("-o", "--outfolder"), type="character", default=NULL,
              help="Path to the outfolder", metavar="character"),
  make_option(c("-a", "--phage"), type="character", default="../20240321_ECOR_Laub_map_phages/ECOR_L_phages_coordinates_20240425.tsv",#NULL,
              help="Mapping of phages to ECOLI genomes", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
##
#Get missing parameters
if (is.null(opt$`blast-file`)) {
  opt$`blast-file`<-paste0(opt$genomeid,"_nucl_cov_megablast.txt")
} 

if (is.null(opt$`genomad-folder`)) {
  opt$`genomad-folder`<-paste0(opt$genomeid,"_cons_fdr/")
} 

if (is.null(opt$`checkv-folder`)) {
  opt$`checkv-folder`<-paste0(opt$genomeid,"_cons_fdr_checkv/")
} 
#####

getpathparts<-strsplit(opt$`genomad-folder`, split = "/")
PreoutfilePrefix<-getpathparts[[1]]
genomeid<-opt$genomeid#str_replace(PreoutfilePrefix,"_cons_fdr","")
outfilePrefix<-str_split(genomeid,"\\.")[[1]][1]

if(is.null(opt$outfolder)) {
  opt$outfolder<-paste0("./",outfilePrefix,"_results_with_phage")
}

folderForResults<-opt$outfolder
if (!dir.exists(folderForResults)){
  dir.create(folderForResults)
}else{
  print("dir exists")
}
#tsv out
outtsvfolder<-"prophage_prediction_tsvs_20240502/"
if (!dir.exists(outtsvfolder)){
  dir.create(outtsvfolder)
}else{
  print("dir exists")
}

#### New 20240321 Read phage data
PhageMapping<-read.csv(opt$phage, header = F, sep="\t")
#get normal start and end
PhageMapping$start<-ifelse(PhageMapping$V9 < PhageMapping$V10,
                           PhageMapping$V9,
                           PhageMapping$V10)

PhageMapping$end<-ifelse(PhageMapping$V9 > PhageMapping$V10,
                         PhageMapping$V9,
                         PhageMapping$V10)

#####

####Read geNomad results
genomadMain<-read.csv(paste0(opt$`genomad-folder`,"/",genomeid,
                             "_summary/",genomeid,
                             "_virus_summary.tsv"),
                      sep="\t",
                      header=T)
genomadMain<-separate(data = genomadMain, col = coordinates, into = c("start", "end"), sep = "-")
genomadMain$start<-as.integer(genomadMain$start)
genomadMain$end<-as.integer(genomadMain$end)

##filter only prophages
#it is important, because otherwise I can pick up all the contamination from some short contigs
genomadMainFiltered<-subset(genomadMain, genomadMain$topology == "Provirus")

#get information about genes
genomadGenes<-read.csv(paste0(opt$`genomad-folder`,"/",genomeid,
                              "_annotate/",genomeid,
                              "_genes.tsv"),
                       sep="\t",
                       header=T)
genespldf<-data.frame(do.call('rbind',strsplit(genomadGenes$gene, "_\\s*(?=[^_]+$)", perl=TRUE)))
genomadGenes$contig<-genespldf$X1

####Read checkv results
checkvMain<-read.csv(paste0(opt$`checkv-folder`,
                            "/quality_summary.tsv"),
                     sep="\t",
                     header=T)
if("Yes" %in% checkvMain$provirus)
{
  checkvContamination<-read.csv(paste0(opt$`checkv-folder`,
                                       "/contamination.tsv"),
                                sep="\t",
                                header=T)
  contList<-str_split(checkvContamination$region_types,",")
  coordList<-str_split(checkvContamination$region_coords_bp,",")
  #virpos<-vector()
  checkv_viral_region_coords<-vector()
  for (j in c(1:length(contList)))
  {
    #j<-2
    if(is.na(contList[[j]][1]))
    {
      #virpos<-c(virpos,NA)
      checkv_viral_region_coords<-c(checkv_viral_region_coords,NA)
      
    }else{
      ###it has to be a little bit more complicated,
      ###because there are situation where there are viral,host,viral
      posinrow<-which(contList[[j]]=="viral")
      if(length(posinrow) == 1)
      {
        #virpos<-c(virpos,posinrow)
        checkv_viral_region_coords<-c(checkv_viral_region_coords,coordList[[j]][posinrow])
      }
      else{
        ##if it is a complicated situation with mixture of viral and host
        ##I amnot trying to update coordinates
        checkv_viral_region_coords<-c(checkv_viral_region_coords,NA)
      }
    }
  }
  #tmpdf<-data.frame(do.call("bind_rows",str_split(checkvContamination$region_coords_bp,",")))
  checkvMain$newcoords<-checkv_viral_region_coords
  
}else{
  ###resolving situation where no proviruses in checkv as with ECOR08
  checkvMain$newcoords<-rep(NA, length(checkvMain$contig_id))
  
}

#####
#Merge geNomad and checkv outputs
GenomadCheckvMain<-merge(genomadMainFiltered, checkvMain, by.x="seq_name", by.y = "contig_id")
GenomadCheckvMain<-separate(data = GenomadCheckvMain, col = newcoords, into = c("checkv_nstart", "checkv_nend"), sep = "-")

GenomadCheckvMain$checkv_nstart<-as.integer(GenomadCheckvMain$checkv_nstart)
GenomadCheckvMain$checkv_nend<-as.integer(GenomadCheckvMain$checkv_nend)

GenomadCheckvMain$checkv_adjstart<-ifelse(is.na(GenomadCheckvMain$checkv_nstart),GenomadCheckvMain$start, 
                                          GenomadCheckvMain$start+
                                            GenomadCheckvMain$checkv_nstart-1)
GenomadCheckvMain$checkv_adjend<-ifelse(is.na(GenomadCheckvMain$checkv_nend),GenomadCheckvMain$end, 
                                        GenomadCheckvMain$start+
                                          GenomadCheckvMain$checkv_nend-1)

####Read blast coverage and do some estimates
BlastCovQuery<-read.csv(opt$`blast-file`, sep="\t", header=F)
names(BlastCovQuery)<-c("QAcc","Pos","Cov")
BlastCovQuery$Cov<-as.integer(BlastCovQuery$Cov)

#learn more about genome
ChromosomesInfo<-BlastCovQuery%>% group_by(QAcc) %>%
  summarize(max=max(Pos))

###################################
findBorders<-function(row)
{
  #row<-1
  flank<-as.integer(opt$flanks)
  #genomad
  RowDf<-GenomadCheckvMain[row,]
  prophageid<-RowDf[,"seq_name"]
  seqname<-unlist(str_split(prophageid,"\\|"))[1]
  
  start<-RowDf[,"checkv_adjstart"]
  end<-RowDf[,"checkv_adjend"]
  
  #blast
  startfl<-ifelse((start -flank)> 0, 
                  (start -flank),
                  1)
  endfl<-ifelse((end +flank)< as.integer(ChromosomesInfo[ChromosomesInfo$QAcc==seqname,2]), 
                (end +flank),
                as.integer(ChromosomesInfo[ChromosomesInfo$QAcc==seqname,2]))
  
  ####get status of prediction
  status<-ifelse((start -flank)< 0 & (end +flank)> as.integer(ChromosomesInfo[ChromosomesInfo$QAcc==seqname,2]), 
                 "ShortFlanks",
                 ifelse((start -flank)< 0, "ShortLeftFl",
                        ifelse((end +flank)> as.integer(ChromosomesInfo[ChromosomesInfo$QAcc==seqname,2]), "ShortRightFl","FlanksPresent")))
  
  ###get genes
  genomadGenesSub<-subset(genomadGenes,
                          genomadGenes$contig == seqname &
                            genomadGenes$start>= startfl &
                            genomadGenes$end <= endfl)
  genomadGenesSub$forward<-ifelse(genomadGenesSub$strand == 1, T,F)
  
  ##
  
  BlastCovSub<-subset(BlastCovQuery,
                      BlastCovQuery$QAcc == seqname &
                        BlastCovQuery$Pos >= startfl &
                        BlastCovQuery$Pos <= endfl)
  #####New 20240321 get quantiles of coverage for flanks only
  BlastCovFlanks<-subset(BlastCovQuery,
                         BlastCovQuery$QAcc == seqname &
                           ((BlastCovQuery$Pos >= startfl & BlastCovQuery$Pos < start) |
                              (BlastCovQuery$Pos < endfl & BlastCovQuery$Pos >= end)))
  MedianCoverageInFlanks<-median(BlastCovFlanks$Cov)
  
  ##Median coverage per window
  ##Windiow size of 5000 nt allows to easily skip regions with common genes
  
  getCoverageInWindow<-function(Df,windowsize){
    IntDf<-Df
    MedianPerWindowVec<-c()
    WindowVec<-c()
    numBreaks<-nrow(IntDf)%/% windowsize +1
    for (i in seq(numBreaks)){
      if (((i-1)*windowsize+1) < nrow(IntDf)){
        partCov<- IntDf[((i-1)*windowsize+1):(min(nrow(IntDf), i*windowsize)), ]
        MedianCovPerWindow<-median(partCov$Cov)
        MedianPerWindowVec<-c(MedianPerWindowVec, rep(MedianCovPerWindow, nrow(partCov)))
        WindowVec<-c(WindowVec,rep(i,nrow(partCov)))
      }
    }
    IntDf$CovMedianWindow<-MedianPerWindowVec
    IntDf$WindowID<-WindowVec
    return(IntDf)
  }
  
  BlastCovSub<-getCoverageInWindow(BlastCovSub,opt$windowSizeLarge)
  
  ############
  #find start and end in large windows
  
  #select large windows and prefilter them based on derivative
  LargeWindowDf<-BlastCovSub
  LargeWindowDf$Deriv<-(LargeWindowDf$CovMedianWindow - lag(LargeWindowDf$CovMedianWindow, n=1))/
    (LargeWindowDf$Pos - lag(LargeWindowDf$Pos,n=1))
  LargeWindowDfNoSmall<-subset(LargeWindowDf,
                               abs(LargeWindowDf$Deriv) > MedianCoverageInFlanks*opt$minimumderivative)
  
  ##Add backup start or end if left or right flank is missing
  if(status == "ShortLeftFl")
  {
    if(LargeWindowDf[1,]$CovMedianWindow<MedianCoverageInFlanks)
    {
      MissingStart<-LargeWindowDf[1,]
      MissingStart$Deriv<-c(-10000)
      LargeWindowDfNoSmall<-rbind(MissingStart,LargeWindowDfNoSmall)
    }
  }else if(status == "ShortRightFl"){
    if(LargeWindowDf[1,]$CovMedianWindow < MedianCoverageInFlanks)
    {
      MissingEnd<-LargeWindowDf[nrow(LargeWindowDf),]
      MissingEnd$Deriv<-c(10000)#this is here for the flanks
      LargeWindowDfNoSmall<-rbind(LargeWindowDfNoSmall,MissingEnd)
    }
  }
  
  
  
  #####This block is obviously the wrong way to filter things, and that is why I am getting artifacts,
  #####Because window border can be within the pit or on the other hand on the cliff, and it affects
  #####the median window coverage
  # if(nrow(LargeWindowDfNoSmallEnds)>1)
  # {
  #   LargeWindowDfNoSmallEnds<-subset(LargeWindowDfNoSmallEnds, LargeWindowDfNoSmallEnds$CovMedianWindow > MedianCoverageInFlanks*(1-opt$deviationincoverage) &
  #                                     LargeWindowDfNoSmallEnds$CovMedianWindow < MedianCoverageInFlanks*(1 + opt$deviationincoverage))
  # 
  # }
  
  
  ####The best way as I see it now (2024/04/24) is to find longest fragment 
  ####And theb extend it if outside derivatives of larger magnitude are available
  LargeWindowDfNoSmall$length<-lead(LargeWindowDfNoSmall$Pos)-LargeWindowDfNoSmall$Pos
  LargeWindowDfNoSmall$leadPos<-lead(LargeWindowDfNoSmall$Pos)
  LargeWindowDfNoSmall$leadDeriv<-lead(LargeWindowDfNoSmall$Deriv)
  LargeWindowDfNoSmall$leadWinID<-lead(LargeWindowDfNoSmall$WindowID)
  
  ##Subset only Pits
  LargePitsDf<-subset(LargeWindowDfNoSmall, LargeWindowDfNoSmall$Deriv < 0 & LargeWindowDfNoSmall$leadDeriv >0)
  #Check whether just take the longest of extend
  LongestPit<-LargePitsDf[order(-LargePitsDf$length),][1,]
  
  #introduce borders before if statement
  LeftBlastBorder<-NA
  RightBlastBorder<-NA
  LengthBlastVsCheckV<-NA
  AgreementCheckV<-NA
  
  if(nrow(LargePitsDf)>0)
  {
    
    #get start
    LargeWindowStart <- LongestPit$Pos
    #skip small segments of elevated coverage
    MinRow <- LargePitsDf[LargePitsDf$Deriv == min(LargePitsDf$Deriv), ]
    if ((LongestPit$WindowID - MinRow$leadWinID) <3 &
        (LongestPit$WindowID - MinRow$leadWinID) >0)
    {
      LargeWindowStart <- MinRow$Pos
    }
    #get end
    LargeWindowEnd<-LongestPit$leadPos
    MaxRow <- LargePitsDf[LargePitsDf$leadDeriv == max(LargePitsDf$leadDeriv), ]
    #upd 2024/07/01 In case if MaxRow longer than one line
    MaxRowOne<-subset(MaxRow, abs(MaxRow$leadPos-LargeWindowEnd) == min(abs(
      MaxRow$leadPos-LargeWindowEnd)))
    if ((MaxRowOne$WindowID - LongestPit$leadWinID) < 3 &
        (MaxRowOne$WindowID - LongestPit$leadWinID) >0)
    {
      LargeWindowEnd <- MaxRow$leadPos
    }
    
    #check if there are additional segment up or down is needed
    
    #do two cycles of improvement just in case
    CheckSegSorted<-LargeWindowDfNoSmall[order(LargeWindowDfNoSmall$Pos),]
    #end
    for(j in 1:2)
    {
      if( ! is.na(match(LargeWindowEnd,CheckSegSorted$Pos)))
      {
        if(! is.na(CheckSegSorted[match(LargeWindowEnd,CheckSegSorted$Pos),]$leadDeriv)){
          if(CheckSegSorted[match(LargeWindowEnd,CheckSegSorted$Pos),]$leadDeriv > CheckSegSorted[match(LargeWindowEnd,CheckSegSorted$leadPos),]$leadDeriv & 
             CheckSegSorted[match(LargeWindowEnd,CheckSegSorted$Pos),]$Deriv > 0  & 
             CheckSegSorted[match(LargeWindowEnd,CheckSegSorted$Pos),]$leadDeriv > 0)
          {
            LargeWindowEnd<-CheckSegSorted[match(LargeWindowEnd,CheckSegSorted$Pos),]$leadPos
          }
        }
      }
      
      #start
      if( ! is.na(match(LargeWindowStart,CheckSegSorted$leadPos)))
      {
        if(CheckSegSorted[match(LargeWindowStart,CheckSegSorted$leadPos),]$Deriv < CheckSegSorted[match(LargeWindowStart,CheckSegSorted$Pos),]$Deriv & 
           CheckSegSorted[match(LargeWindowStart,CheckSegSorted$leadPos),]$Deriv < 0  & 
           CheckSegSorted[match(LargeWindowStart,CheckSegSorted$Pos),]$Deriv < 0)
        {
          LargeWindowStart<-CheckSegSorted[match(LargeWindowStart,CheckSegSorted$leadPos),]$Pos
        }
      }
    }
    
    
    
    BordersDfLongest<-BlastCovSub[BlastCovSub$Pos == LargeWindowStart |  BlastCovSub$Pos == LargeWindowEnd,]
    
    #############
    ###deal with the fact that fragments have no or very short left or right flanks
    
    
    #####Do a second pass with a smaller window

    
    ###If the region is found the idea here is to go with smaller windows to smooth the things and get more accurate borders
    getstartendDf<-function(Df,rownum)
    {
      SubsetDf<-subset(Df,
                       Df$WindowID == (BordersDfLongest[rownum,c("WindowID")]-1) |
                         Df$WindowID == BordersDfLongest[rownum,c("WindowID")])
      return(SubsetDf)
    }
    StartWinSecondDf<-getstartendDf(BlastCovSub,1)
    StartWinSecondDf<-getCoverageInWindow(StartWinSecondDf,opt$windowSizeSmall)
    EndWinSecondDf<-getstartendDf(BlastCovSub,2)
    EndWinSecondDf<-getCoverageInWindow(EndWinSecondDf,opt$windowSizeSmall)
    
    #####New April 02, Calculate derivative to get exact border
    derivwindowsize<-50
    
    getDerivStartEndSimple<-function(Df,ColName,derivwindow,start=F)
    {
      # Df<-EndWinSecondDf
      # derivwindow<-1
      # ColName<-"CovMedianWindow"
      IntDf<-Df
      IntDf$Deriv<-(lead(Df[[ColName]], n=derivwindow) - Df[[ColName]])/
        (lead(Df$Pos,n=derivwindow)- Df$Pos)
      EndToReturn<--1000
      if(start)
      {
        EDf<-head(IntDf[order(IntDf$Deriv,-IntDf$Pos),],n=1)
        if(EDf$Deriv < -1)
        {
          EndToReturn<-EDf$Pos
        }
      }else{
        EDf<-head(IntDf[order(-IntDf$Deriv,IntDf$Pos),],n=1)
        if(EDf$Deriv > 1)
        {
          EndToReturn<-EDf$Pos
        }
        else{
          ###if it does not work return thelast position in window [upd 2024/07/01]
          EndToReturn<-max(IntDf$Pos)
        }
      }
      return(EndToReturn)
    }
    
    SmallWindowStart<-getDerivStartEndSimple(StartWinSecondDf,"CovMedianWindow",1,start=T)
    SmallWindowEnd<-getDerivStartEndSimple(EndWinSecondDf,"CovMedianWindow",1)
    
    ####Now get the exact point without windows
    
    BlastCovSub$Deriv<-(lead(BlastCovSub$Cov, n=derivwindowsize) - BlastCovSub$Cov)/(lead(BlastCovSub$Pos,n=derivwindowsize)- BlastCovSub$Pos)
    
    StartNuclSpace<-subset(BlastCovSub, (BlastCovSub$Pos > SmallWindowStart-opt$windowSizeSmall-1) &
                             (BlastCovSub$Pos < SmallWindowStart + opt$windowSizeSmall))
    
    EndNuclSpace<-subset(BlastCovSub, (BlastCovSub$Pos > SmallWindowEnd-opt$windowSizeSmall-1) &
                           (BlastCovSub$Pos < SmallWindowEnd + opt$windowSizeSmall))
    
    #get start
    BordersDfStart<-head(StartNuclSpace[order(StartNuclSpace$Deriv,-StartNuclSpace$Pos),],n=1)
    
    BackupStart<-subset(StartNuclSpace,StartNuclSpace$Pos == SmallWindowStart)
    
    BordersDfStart<-rbind(BordersDfStart, BackupStart)
    
    #######get end
    BordersDfEnd<- head(EndNuclSpace[order(-EndNuclSpace$Deriv,EndNuclSpace$Pos),],n=1)
    #Get backup on the edge of window
    BackupEnd<-subset(EndNuclSpace,EndNuclSpace$Pos == SmallWindowEnd)
    BordersDfEnd<-rbind(BordersDfEnd, BackupEnd)
    
    #####Finally get borders
    #remove them where there are no flanks
    
    LeftBlastBorder<-ifelse(BordersDfLongest[1,]$Pos == 1, 1,
                            ifelse((nrow(BordersDfStart)==0 |
                                      BordersDfStart[1,]$Pos == -1000), NA, BordersDfStart[1,]$Pos))
    RightBlastBorder<-ifelse((nrow(BordersDfEnd)==0 |
                                BordersDfEnd[1,]$Pos == -1000), NA, BordersDfEnd[1,]$Pos)
    
    ###Get intersection between checkv and blast
    LengthBlastVsCheckV<-(RightBlastBorder-LeftBlastBorder+1)*100/(end - start+1)
    AgreementCheckV<-ifelse((start>=LeftBlastBorder) & (end<=RightBlastBorder), 100,
                            ifelse(start< LeftBlastBorder, (end - LeftBlastBorder+1)*100/(end - start+1),
                                   (RightBlastBorder - start+1)*100/(end - start+1)))
    
    #####additional criteria if blast borders are way longer on one side, then use genomad border
    ########
    if(abs(LengthBlastVsCheckV -100) >= 20 & abs(LengthBlastVsCheckV -100) < 40){
      if((RightBlastBorder -RowDf$end)> 0.1*(RightBlastBorder-LeftBlastBorder+1))
      {
        RightBlastBorder<-RowDf$end
      }
      if((RowDf$start - LeftBlastBorder)> 0.1*(RightBlastBorder-LeftBlastBorder+1))
      {
        LeftBlastBorder<-RowDf$start
      }
    }
    
  }
  #####Subsampling and sliding window for vizualization
  #I use it only for vizualization purposes
  BlastCovSubLowDens<-sample_n(BlastCovSub, 10000) # I can adjust this parameter later on
  
  ###Do the plot for the predicted prophage region with all information combined
  #It is done like that because otherwise ggarrows cannot be aligned
  
  
  #annotationborder
  AnnotStart<-max(RightBlastBorder, end, RowDf$end, na.rm = T)
  
  Maxcov<-max(BlastCovSubLowDens$Cov)
  
  genomadcheckvplot<-ggplot()+
    #genes
    geom_hline(yintercept = Maxcov+Maxcov*0.150,color = "#969696")+
    geom_gene_arrow(data=genomadGenesSub,
                    aes(xmin=as.integer(start),
                        xmax=as.integer(end),
                        forward=forward,
                        y=Maxcov+Maxcov*0.150,
                        fill=as.factor(virus_hallmark)))+
    scale_fill_manual(values=c("#bdbdbd","#f03b20"),
                      labels=c("not hallmark","viral hallmark"),
                      name="Genes annotation")+
    geom_text(data=genomadGenesSub,
              aes(x=start+(end-start)/2,
                  y=Maxcov+Maxcov*0.180, label=annotation_description),
              size=2.5,
              angle=45,
              hjust=0, vjust=1,
              na.rm=T)+
    #predictions
    geom_rect(data=RowDf, aes(xmin=start,xmax=end,ymin=Maxcov+Maxcov*0.070,ymax=Maxcov+Maxcov*0.085),
              fill="#0570b0")+
    annotate(geom="text",y=Maxcov+Maxcov*0.080,x=AnnotStart+150,label="geNomad",hjust = 0)+
    #checkv prediction
    new_scale_fill()+
    # #get cutoff level
    # geom_hline(yintercept = CoverageCutOff,
    #            color = "#800026", size=0.5)+
    geom_rect(data=RowDf, aes(xmin=checkv_adjstart,xmax=checkv_adjend,
                              fill=completeness,
                              ymin=Maxcov+Maxcov*0.045,ymax=Maxcov+Maxcov*0.060))+#,
    #fill="#88419d")+
    scale_fill_gradient(low="#f7fbff",high="#08306b", name = "CheckV completeness",
                        limits=c(0,100))+
    annotate(geom="text",y=Maxcov+Maxcov*0.055,x=AnnotStart+150,label="checkv", hjust= 0)
  if(! (is.na(LeftBlastBorder) | is.na(RightBlastBorder)))
  {
    #coverage
    genomadcheckvplot<- genomadcheckvplot + geom_rect(aes(xmin=LeftBlastBorder,xmax=RightBlastBorder,ymin=Maxcov+Maxcov*0.020,ymax=Maxcov+Maxcov*0.035),
                                                      fill="#31a354")+
      annotate(geom="text",y=Maxcov+Maxcov*0.035,x=AnnotStart+150,label="blast", hjust= 0)
  }
  
  genomadcheckvplot<-genomadcheckvplot+
    geom_line(data=BlastCovSubLowDens,
              aes(x=Pos,y=CovMedianWindow),
              size=0.5,
              color = "#31a354")+
    geom_point(data=BlastCovSubLowDens,
               aes(x=Pos,y=Cov),
               size=0.5,
               shape=16,
               color = "#004529")+
    #scale_y_continuous(limits=c(0,Maxcov+Maxcov*0.5), breaks=seq(0,Maxcov,100), name="# of BLAST hits")+
    scale_y_continuous(breaks=seq(0,Maxcov,100), name="# of BLAST hits",
                       limits = c(NA,Maxcov*1.6))+
    scale_x_continuous(limits = c(startfl,endfl+700), breaks=seq(0,10000000,10000), 
                       labels=seq(0,10000,10), name="genome position (kb)")+
    ggtitle(prophageid)+
    theme_minimal()+
    theme(legend.position = "bottom")
  genomadcheckvplot
  
  ####New 20240321 add info about phage is present
  
  GenomeOfInt<-str_replace(opt$`blast-file`,"_nucl_cov_megablast.txt","")
  PhageSubsetDf<-subset(PhageMapping, PhageMapping$V1 == GenomeOfInt &
                          PhageMapping$start > startfl &
                          PhageMapping$end < endfl)
  
  if (length(PhageSubsetDf$V1)> 0)
  {
    PhageSubsetDf$ymin<-seq(-Maxcov*0.03,length(PhageSubsetDf$V1)*(-Maxcov*0.03), -Maxcov*0.03)
    genomadcheckvplot<-genomadcheckvplot +new_scale_fill()+
      
      geom_rect(data=PhageSubsetDf,
                aes(xmin=start, xmax=end, ymin=ymin, ymax=ymin+Maxcov*0.02, fill = V2))
    #fill = "#cc4c02")
  }
  
  
  
  
  
  returnList<-list("graph"=genomadcheckvplot,
                   "blastborders"=c(LeftBlastBorder,RightBlastBorder),
                   "status"=status,
                   "IntersectionWithCheckV"=AgreementCheckV,
                   "LengthBlastVsCheckV"=LengthBlastVsCheckV)
  
  return(returnList)
}
###############

BlastBordersDf<-data.frame(seq_name=vector(),
                           blast_l=vector(),
                           blast_r=vector(),
                           status=vector())

for (j in c(1:nrow(GenomadCheckvMain)))
{
  #j<-5
  rowout<-findBorders(j)
  
  tmpdf<-data.frame(seq_name=GenomadCheckvMain[j,"seq_name"],
                    blast_l=rowout[["blastborders"]][1],
                    blast_r=rowout[["blastborders"]][2],
                    status=rowout[["status"]],
                    IntersectionWithCheckV=rowout[["IntersectionWithCheckV"]],
                    LengthBlastVsCheckV=rowout[["LengthBlastVsCheckV"]])
  BlastBordersDf<-rbind(BlastBordersDf,tmpdf)
  
  #Save plots to results folder
  ggsave(paste0(outfilePrefix,"_",str_replace_all(GenomadCheckvMain[j,"seq_name"], c(`\\|`="_", `\\.`="_")),"_plot.pdf"),
         plot=rowout[["graph"]],
         path=folderForResults,
         width=55, height =25, units="cm",dpi=300,
         bg = "white")
}

#Save merged table with predictions for the genome
AddBlastDf<-merge(genomadMain,BlastBordersDf, by="seq_name")

checkvToSave<-GenomadCheckvMain[,c("seq_name","provirus","checkv_quality",
                                   "completeness",
                                   "completeness_method",
                                   "checkv_adjstart",
                                   "checkv_adjend")]
names(checkvToSave)<-c("seq_name","checkv_provirus","checkv_quality",
                       "checkv_completeness",
                       "checkv_completeness_method",
                       "checkv_adjstart",
                       "checkv_adjend")

OutputDfToSave<-merge(AddBlastDf,checkvToSave,
                      by="seq_name")
OutputDfToSave$blastLength<-OutputDfToSave$blast_r - OutputDfToSave$blast_l+1
OutputDfToSave$genome<-rep(genomeid,nrow(OutputDfToSave))

write.table(OutputDfToSave,file=paste0(folderForResults,"/",genomeid,"_prophage_prediction_output.tsv"),
            sep="\t", quote = F, row.names = F)

























