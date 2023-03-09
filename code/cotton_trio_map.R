#### Cotton map ####
library(readxl)
source("/mnt/repo/trio_genetic_map/code/trio_mapping_functions.R")


## step0. load pedigree data, marker map and genotype data
dataFolder = "/mnt/cotton/data/"
genoData = read.table(paste0(dataFolder,"Trios_unique_lines_FP_calls_from_MarkerCallSetsV2_aggregated_02272023.txt"),
                      header=T,check.names = F)
row.names(genoData) = genoData$marker_name
genoData = genoData[,-1]
pedData= read_excel(paste0(dataFolder,"Trios_list_and_subject_Id_for_Tao.xlsx"),sheet=1)
pedData = as.data.frame(pedData %>% select(lineCode_subject_ID,parent1_subject_ID,parent2_subject_ID))
colnames(pedData) = c("pedigree","prt1","prt2")
markerData = read.delim(paste0(dataFolder,"Physical_genetic_map_FP_markers_with_DP393_09292022.txt"),header=T)
markerData = markerData %>% arrange(DP393_Chr_Andrei,DP393_V2_Pos_Andrei)
row.names(markerData) = markerData$MARKER_ROOT_NAME


# step1. informative markers per line
tmpFolder = "/mnt/cotton/tmpData"
infoMarkerByLine(pedData,genoData,nCores = 30, tmpFolder, genos = c("AA","TT","CC","GG"))
linePrtInfo = collectAllInfoMarkers(tmpFolder)
linePrtInfo$p1Pct = round(linePrtInfo$p1Count / (linePrtInfo$p1Count + linePrtInfo$p2Count),2)
hist(linePrtInfo$p1Pct)
qstLines = subset(linePrtInfo,p1Pct < 0.3 | p1Pct > 0.7)$line # lines are questionable, with excessive p1 or p2 ratio
file.remove(paste0(tmpFolder,"/",paste0(qstLines,".csv"))) # remove the data from questionable lines

# step3. collect and get rec Data
generateRecData(tmpFolder,row.names(genoData),outFileFolder = "/mnt/cotton/recResults",verbose = TRUE)

# step4. 

totalTrioCount <- readRDS(paste0(outFileFolder,"/totalTrioCount_allChr.rds"))
totalRecRate <- readRDS(paste0(outFileFolder,"/recRate_allChr.rds"))

phyPos = markerData %>% filter(str_detect(DP393_Chr_Andrei,"Gh_DP393")) %>%select(MARKER_ROOT_NAME,DP393_Chr_Andrei,DP393_V2_Pos_Andrei)
colnames(phyPos) = c("Q_ID", "chr", "Sbegin")

mkLinkageInfo = getLinkedMarkerInfo(phyPos,colnames(totalRecRate),totalTrioCount,totalRecRate,minCount = 100,maxRecRate = 0.1)
mkLinkageInfo$chrRatio = round(mkLinkageInfo$totalMkLinkedOwnChr / mkLinkageInfo$totalMkLinked,2)
pbMkData = mkLinkageInfo %>% filter((totalMkLinkedOwnChr == 1 & numLinkedChr >= 2) | chrRatio <= 0.5) # a singleton or less linked with own chr 
pbMarkers = pbMkData$Q_ID

mkLinkageInfo2 = getLinkedMarkerInfo(phyPos,mkLinkageInfo$Q_ID[!mkLinkageInfo$Q_ID %in% pbMarkers], totalTrioCount, totalRecRate,minCount=100, maxRecRate=0.1)
pbMkData2 = mkLinkageInfo2 %>% filter((totalMkLinkedOwnChr == 1 & numLinkedChr >= 2) | totalMkLinkedOwnChr/totalMkLinked <= 0.5) # a singleton or less linked with own chr 
pbMarkers = c(pbMarkers,pbMkData2$Q_ID)

mkLinkageInfo2 = mkLinkageInfo2 %>% select(-c("LDMarkers"))
timestamp()


### step5.
tmpOutFolder = "/mnt/cotton/tmpData/tmpAnchorData"
for(chrN in unique(mkLinkageInfo2$chr)){
  print(paste(chrN, Sys.time()))
  findAnchorsByChr(chrN,mkLinkageInfo2,totalTrioCount,totalRecRate,tmpOutFolder,minCount = 100,minLinkedMk = 5,
                   mkPerGroup = 5,minDataPoint = 20, minDataPointPerMk = 2, maxDist = 0.5,numOfDataBwGroups = 5)
}


anchorResAllChrs = data.frame()
for(testChr in unique(mkLinkageInfo2$chr)){
  tmp = read.csv(paste0(tmpOutFolder,"/anchorRes_chr",testChr,".csv"))
  tmp$chr = testChr
  anchorResAllChrs = bind_rows(anchorResAllChrs,tmp)
}

anchorResAllChrs$maxAnchorDist = sapply(anchorResAllChrs$distVec,function(x){max(as.numeric(unlist(strsplit(x,"_"))),na.rm=T)})

ggplot(anchorResAllChrs,aes(as.character(chr),totalGenLen,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),MAE,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),MAE_pct,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),maxAnchorAE,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),maxAnchorDist,col=as.character(chr))) + geom_boxplot() + theme_bw()


#### step6. select anchor path ####
# request at least 10 groups, with chr5 as a exception as the max is 9; 
selAnchorPath = data.frame()
for(testChr in sort(unique(mkLinkageInfo2$chr))){
  tmp = read.csv(paste0(tmpOutFolder,"/anchorRes_chr",testChr,".csv"))
  tmp$chr = testChr
  tmp = tmp %>% filter(numOfGroups >= 10 | numOfGroups == max(numOfGroups)) %>% arrange(MAE_pct)
  selAnchorPath = bind_rows(selAnchorPath,tmp[1,])
}

selAnchorPath$maxAnchorDist = sapply(selAnchorPath$distVec,function(x){max(as.numeric(unlist(strsplit(x,"_"))),na.rm=T)})


