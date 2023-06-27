#### Cotton map ####
library(readxl)
source("/mnt/repo/trio_genetic_map/code/trio_mapping_functions.R")
 
 
#### step0. load pedigree data, marker map and genotype data ####
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
 
 
#### step1. informative markers per line ####
tmpFolder = "/mnt/cotton/tmpData/tmp"
if(!dir.exists(tmpFolder)){dir.create(tmpFolder)}
infoMarkerByLine(pedData,genoData,nCores = 30, tmpFolder, genos = c("AA","TT","CC","GG"))
linePrtInfo = collectAllInfoMarkers(tmpFolder)
linePrtInfo$p1Pct = round(linePrtInfo$p1Count / (linePrtInfo$p1Count + linePrtInfo$p2Count),2)
hist(linePrtInfo$p1Pct)
qstLines = subset(linePrtInfo,p1Pct < 0.3 | p1Pct > 0.7)$line # lines are questionable, with excessive p1 or p2 ratio
file.remove(paste0(tmpFolder,"/",paste0(qstLines,".csv"))) # remove the data from questionable lines
 
#### step2. collect and get rec Data ####
#generateRecData(tmpFolder,row.names(genoData),outFileFolder = "/mnt/cotton/recResults",verbose = TRUE)
 
#### step3.  ####
outFileFolder = "/mnt/cotton/recResults"
totalTrioCount <- readRDS(paste0(outFileFolder,"/totalTrioCount_allChr.rds"))
totalRecRate <- readRDS(paste0(outFileFolder,"/recRate_allChr.rds"))
 
phyPos1 = read.delim("/mnt/cotton/data/physical_position_DP393_HIFI_GMAP.txt",header=T)
phyPos1$SNPpos = as.numeric(phyPos1$SNPpos)
phyPos1$SNPposInferred = as.numeric(phyPos1$SNPposInferred)
 
phyPos2 = read.delim("/mnt/cotton/data/physical_position_DP393_HIFI_GMAP_part2.txt",header=T)
phyPosAll = bind_rows(phyPos1,phyPos2)
 
phyPos = phyPosAll %>% 
  mutate(chr = as.numeric(str_replace_all(Subj,"GOSHI_DP393_v2_0_0_chr",""))) %>% 
  filter(Coverage >= 90, PercentID >= 90, chr <= 1000000) %>% 
  arrange(chr,Sbegin) %>% 
  dplyr::select(Q_ID,chr,Sbegin,SNPpos)
row.names(phyPos) = phyPos$Q_ID
 
 
#identify markers that are questionable (not consistent beteen physical and genetic chrs)
testMarkers = phyPos$Q_ID[phyPos$Q_ID %in% colnames(totalTrioCount)]
mkLinkageInfo = getLinkedMarkerInfo(phyPos,testMarkers,totalTrioCount,totalRecRate,minCount = 100,maxRecRate = 0.1)
mkLinkageInfo$chrRatio = round(mkLinkageInfo$totalMkLinkedOwnChr / mkLinkageInfo$totalMkLinked,2)
pbMkData = mkLinkageInfo %>% filter((totalMkLinkedOwnChr == 1 & numLinkedChr >= 2) | chrRatio <= 0.5) # a singleton or less linked with own chr 
pbMarkers = pbMkData$Q_ID
 
mkLinkageInfo2 = getLinkedMarkerInfo(phyPos,mkLinkageInfo$Q_ID[!mkLinkageInfo$Q_ID %in% pbMarkers], totalTrioCount, totalRecRate,minCount=100, maxRecRate=0.1)
pbMkData2 = mkLinkageInfo2 %>% filter((totalMkLinkedOwnChr == 1 & numLinkedChr >= 2) | totalMkLinkedOwnChr/totalMkLinked <= 0.5) # a singleton or less linked with own chr 
pbMarkers = c(pbMarkers,pbMkData2$Q_ID)
 
mkLinkageInfo2 = mkLinkageInfo2 %>% select(-c("LDMarkers"))
mkLinkageInfo2 = mkLinkageInfo2 %>% arrange(chr,Sbegin)
timestamp()
 
 
### step4. define anchors and assemble by chr ####
tmpOutFolder = "/mnt/cotton/tmpData/tmpAnchorData" 
 
for(chrN in sort(unique(phyPos$chr))){
  print(paste(chrN, Sys.time()))
  # findAnchorsByChr(chrN,mkLinkageInfo2,totalTrioCount,totalRecRate,tmpOutFolder,minCount = 100,minGroupDist = 5,maxGroupDist = 30,minGroupDistEdge = 1,
  #                  minLinkedMk = 5, mkPerGroup = 1,minDataPoint = 1, minDataPointPerMk = 1, maxDist = 0.5,numOfDataBwGroups = 1) # only one marker in a group
  
  findAnchorsByChr(chrN,mkLinkageInfo2,totalTrioCount,totalRecRate,tmpOutFolder,minCount = 100,minGroupDist = 5,maxGroupDist = 30,minGroupDistEdge = 1,
                   minLinkedMk = 5, mkPerGroup= 2, minDataPoint = 4, minDataPointPerMk = 2, maxDist = 0.5,numOfDataBwGroups = 2,maxRangeBwGroups = 4) # only two markers in a group
  
}
 
 
anchorResAllChrs = data.frame()
for(testChr in sort(unique(phyPos$chr))){
  tmp = read.csv(paste0(tmpOutFolder,"/anchorRes_chr",testChr,".csv"))
  tmp$chr = testChr
  anchorResAllChrs = bind_rows(anchorResAllChrs,tmp)
}
 
anchorResAllChrs$maxAnchorDist = sapply(anchorResAllChrs$distVec,function(x){max(as.numeric(unlist(strsplit(x,"_"))),na.rm=T)})
write.csv(anchorResAllChrs,paste0(tmpFolder,"/anchorResALlChrs.csv"))
 
ggplot(anchorResAllChrs,aes(x = chr,totalGenLen,group=chr,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),MAE,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),MAE_pct,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),maxAnchorAE,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),maxAnchorDist,col=as.character(chr))) + geom_boxplot() + theme_bw()
 
 
#### step5. select anchor path ####
# request at least 10 groups, with chr5 as a exception as the max is 9; 
 
qCut = 0.15 # quantitle cutoff
 
selAnchorPath = data.frame()
for(testChr in sort(unique(phyPos$chr))){
  tmp = read.csv(paste0(tmpOutFolder,"/anchorRes_chr",testChr,".csv"))
  tmp$chr = testChr
  tmp$maxAnchorDist = sapply(tmp$distVec,function(x){if(!is.na(x)){max(as.numeric(unlist(strsplit(x,"_"))),na.rm=T)}else(NA)})
  tmp$beginIndex = sapply(tmp$selIndex,function(x){if(!is.na(x)){return(min(as.numeric(unlist(strsplit(x,"_"))),na.rm=T))}else(NA)})
  tmp$endIndex = sapply(tmp$selIndex,function(x){if(!is.na(x)){return(max(as.numeric(unlist(strsplit(x,"_"))),na.rm=T))}else(NA)})
  qtGroup = quantile(tmp$numOfGroups,na.rm=T,probs = c(qCut,1-qCut))
  qtLen = quantile(tmp$totalGenLen,na.rm=T,probs = c(qCut,1-qCut))
  tmp = tmp %>% filter(numOfGroups >= max(8,qtGroup[1]), totalGenLen >= max(qtLen[1],30), beginIndex <= 30,maxAnchorAE <= 10) %>% arrange(MAE_pct)
  selAnchorPath = bind_rows(selAnchorPath,tmp[1,])
}
 
row.names(selAnchorPath) = selAnchorPath$chr
 
#### step6. all group distances  ####
# step1: map high quality anchors (groups) onto the map #selGroupSm = groupSm %>% filter(numOfData >= 80, minDataPerMk >= 5, max <= 0.5)
# step2: map sub-high quality anchors onto the map (particually those present in large gaps)
# step3: map individual markers onto the map (because of the slightly fluctuation, decide to replace this step by finding anchor markers at every 1 cM distance)
 
markerMapRes = data.frame()
 
largestGaps = rep(NA,nrow(selAnchorPath))
names(largestGaps) = selAnchorPath$chr
 
mkPerGroup = 2 # numOfMarkers per Group
minDataWithAnchor = 1 # the minumum data point between a marker and a anchor (10 markers); 
 
oldNewChrDict = 1:26
names(oldNewChrDict) = selAnchorPath$chr
 
markerMap = markerData
colnames(markerMap)[1:3] = c("markerName","oldMapChr","oldMapPos") # change names to fit for pipeline
 
for(testChr in selAnchorPath$chr){
  print(paste(testChr,Sys.time()))
  selAnchors = as.numeric(unlist(strsplit(selAnchorPath[testChr,"selIndex"],"_")))
  distVect = as.numeric(unlist(strsplit(selAnchorPath[testChr,"distVec"],"_")))
  pwGroupSm = readRDS(paste0(tmpOutFolder,"/pwGroupSm_chr",testChr,".rds"))
  anchorRes = read.csv(paste0(tmpOutFolder,"/anchorRes_chr",testChr,".csv"))
  withGroupSm = read.csv(paste0(tmpOutFolder,"/groupSm_chr",testChr,".csv"))
  goodAnchors = withGroupSm %>% filter(numOfData >= 4, minDataPerMk >= 2, maxD <= 0.5) #
  chrMapDist = readRDS(paste0(tmpOutFolder,"/subChrMapDist_chr",testChr,".rds"))
  # get all anchor positions
  allGroupDist  = goodAnchorDist(goodAnchors,selAnchors,anchorRes,pwGroupSm)
  goodAnchors$groupPos = round(allGroupDist,3)
  goodAnchors = goodAnchors %>% filter(!is.na(groupPos))
  distGap = goodAnchors$groupPos[2:nrow(goodAnchors)] - goodAnchors$groupPos[1:(nrow(goodAnchors) - 1)]
  
  # remove anchors that are in wrong oder, unless there is a big gap
  # if(min(distGap) < -2.5){
  #   print(warning(paste0("Warning!!! chr",testChr," is skipped due to large mis-gaps:",min(distGap))))
  # }
  # 
  # while(any(distGap < 0) & min(distGap) > -2.5){
  #   goodAnchors = goodAnchors[-(which(distGap < 0) + 1),]
  #   distGap = goodAnchors$groupPos[2:nrow(goodAnchors)] - goodAnchors$groupPos[1:(nrow(goodAnchors) - 1)]
  # }
  
  while(any(distGap < 0)){
    goodAnchors = goodAnchors[-(which(distGap < 0) + 1),]
    distGap = goodAnchors$groupPos[2:nrow(goodAnchors)] - goodAnchors$groupPos[1:(nrow(goodAnchors) - 1)]
  }
  
  # interpolate markers onto the map
  subChrData = read.csv(paste0(tmpOutFolder,"/subChrData_chr",testChr,".csv"))
  row.names(subChrData) = subChrData$Q_ID
  newMarkerRes = getMarkerPos(subChrData,goodAnchors,mkPerGroup)
  
  # re-adjust marker position at big gaps (centromere regions are usually wrongly interpolated, eg. chr2 15-30 Mb
  bigGapIndex = which(distGap >= 5)
  for(index in bigGapIndex){
    a1Index = goodAnchors$gIndex[index] + mkPerGroup
    a2Index = goodAnchors$gIndex[index+1] - 1
    tmpPos = reCalMkDist(goodAnchors,index,chrMapDist,newMarkerRes,mkPerGroup,minDataWithAnchor = 1)
    if(!all(is.na(tmpPos))){
      newMarkerRes$newMapPos[a1Index:a2Index] = NA
      newMarkerRes$newMapPos[as.numeric(names(tmpPos))] = tmpPos
      newGaps = diff(c(goodAnchors$groupPos[index],tmpPos,goodAnchors$groupPos[index+1]))
      #if(any(newGaps >= 5)){print(c(goodAnchors$groupPos[index],tmpPos,goodAnchors$groupPos[index+1]))}
    }
  }
  
  # interpolate all other markers with physical positions
  subPhyPos = phyPos %>% filter(chr == testChr) %>% arrange(Sbegin)
  subPhyPos = subPhyPos[!(duplicated(subPhyPos$Q_ID)),]
  row.names(subPhyPos) = subPhyPos$Q_ID
  subPhyPos$newMapPos = newMarkerRes[subPhyPos$Q_ID,"newMapPos"]
  intPosPhy = interpolate(subPhyPos$Sbegin,subPhyPos$newMapPos)
  subPhyPos$newMapPos = intPosPhy[,3]
  
  # interpolation based old gen map for markers without physical positions
  subMarkerMap = markerMap %>% 
    filter(oldMapChr == oldNewChrDict[testChr])
  subMarkerMap$newMapPos = subPhyPos[subMarkerMap$markerName,"newMapPos"]
  subMarkerMap$Sbegin = subPhyPos[subMarkerMap$markerName,"Sbegin"]
  subMarkerMap = subMarkerMap %>% arrange(oldMapPos)
  intIndex = which(is.na(subMarkerMap$Sbegin),is.na(subMarkerMap$newMapPos)) # only interpolate markers that having no physical positions
  intPosGen = interpolate(subMarkerMap$oldMapPos,subMarkerMap$newMapPos)
  subMarkerMap$newMapPos[intIndex] = intPosGen[intIndex,3] # only interpolate markers that having no physical positions
  extraMarkers = subMarkerMap$markerName[!subMarkerMap$markerName %in% phyPos$Q_ID]
  
  # combine phyPos,oldGenMap together
  subPhyPos[extraMarkers,] = NA
  subPhyPos[extraMarkers,"Q_ID"] = extraMarkers
  subPhyPos$oldMapChr = markerMap[subPhyPos$Q_ID,"oldMapChr"]
  subPhyPos$oldMapPos = markerMap[subPhyPos$Q_ID,"oldMapPos"]
  subPhyPos[extraMarkers,"newMapPos"] = subMarkerMap[extraMarkers,"newMapPos"]
  subPhyPos  = subPhyPos %>% arrange(newMapPos,Sbegin)
  
  # combin all chr data together
  markerMapRes = bind_rows(markerMapRes,subPhyPos)
  #print(table(subPhyPos$newMapPos[2:nrow(subPhyPos)] - subPhyPos$newMapPos[1:(nrow(subPhyPos)-1)] < 0))
  largestGaps[testChr] = max(subPhyPos$newMapPos[2:nrow(subPhyPos)] - subPhyPos$newMapPos[1:(nrow(subPhyPos)-1)],na.rm=T)
  
}
 
 
markerMapRes$chr_clr = markerMap[markerMapRes$Q_ID,"DP393_Chr_Andrei"]
markerMapRes$pos_clr = markerMap[markerMapRes$Q_ID,"DP393_V2_Pos_Andrei"]
markerMapRes = markerMapRes %>% arrange(chr,Sbegin)
write.csv(markerMapRes,"markerMapRes_hifi.csv",row.names=F)
saveRDS(markerMapRes,"markerMapRes_hifi.rds")
 
#### step7 manully add missing positions ####
## the beginning and end of each chr are usually missed from interpolate
# rules by priority: 1. linkage data; 2. slope; 3. old map
table(is.na(markerMapRes$newMapPos),markerMapRes$chr) # with phyPos but not new genetic pos
table(is.na(markerMapRes$newMapPos),markerMapRes$oldMapChr) # with old genetic pos but not new genetic pos
 
allChrs = 1:26
names(allChrs) = sort(unique(markerMapRes$chr))
 
#chr1
# the start of chr1
testChr = 1
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
stPhyPos = 0; edPhyPos = 1 * 10^6 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
 
addDist = 10 # add to this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#the end of chr1 
# had some issue of calling positions at the end of chr1
testChr = 1
stPhyPos = 117000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
testData$newMapPos = NA
testData$newMapPos[1] = 114.469
testData$newMapPos[nrow(testData)] = 114.469 + 16.3 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
# chr2
testChr = 2
 
# start of chr2
# 2 NA at the begining and 1 at the end
markerMapRes[c("NGHIR009568920"),"newMapPos"] = 0
 
# end of chr2
# 3 NA without physical position but with old genetic pos
# "NGHIR009405770" is not convinced to place on the new map
stPhyPos = 107600000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 1.23 # add cM to the end;
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
# chr3
#the end of chr
testChr = 3
stPhyPos = 112000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 1.5 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#chr4
testChr = 4
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr4
stPhyPos = 0; edPhyPos = 1500000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 5 # add to this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#the end of chr4
stPhyPos = 86000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 9 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#chr5
testChr = 5
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr5
stPhyPos = 0; edPhyPos = 1700000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 6 # add this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#the end of chr5
stPhyPos = 108000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
 
# the last at least 4 markers should on homeologous chrs D05, but not A05; CLR mapping is correct; These markers also present in the pbMarkers
# NG0210332, NGHIR009570877, NGHIR009574612, NGHIR009539591
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 8.5 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# Would not be able to place their positions NGHIR009578999; Since they are at the terminal in the old map,set their as the end
markerMapRes[c("NGHIR009578999"),"newMapPos"] = 214.2
 
 
#chr6
testChr = 6
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr6
stPhyPos = 0; edPhyPos = 1900000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 17.8 # add this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#the end of chr6
stPhyPos = 126300000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 5 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# another 9 markers without physicla position. 
subMarkerMap = markerMapRes %>% 
  filter(oldMapChr == testChr) %>% 
  arrange(oldMapChr,oldMapPos)
 
intPos = interpolate(subMarkerMap$oldMapPos,subMarkerMap$newMapPos)[,3]
table(is.na(subMarkerMap$newMapPos),intPos == subMarkerMap$newMapPos)
markerMapRes[subMarkerMap$Q_ID,"newMapPos"] = interpolate(subMarkerMap$oldMapPos,subMarkerMap$newMapPos)[,3]
 
#chr7
testChr = 7
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr7
addDist = 1.1 # add this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
markerMapRes["NGHIR009554094","newMapPos"] = 0
 
# end of chr7; limited evidence. 
stPhyPos = 97800000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 4.6 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
###chr8
testChr = 8
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr8
stPhyPos = 0; edPhyPos = 1500000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 3.8 # add this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#the end of chr8
stPhyPos = 123900000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 12 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# another 6 markers without physicla position. 
subMarkerMap = markerMapRes %>% 
  filter(oldMapChr == testChr) %>% 
  arrange(oldMapChr,oldMapPos)
 
intPos = interpolate(subMarkerMap$oldMapPos,subMarkerMap$newMapPos)[,3]
table(is.na(subMarkerMap$newMapPos),intPos == subMarkerMap$newMapPos)
markerMapRes[subMarkerMap$Q_ID,"newMapPos"] =  intPos
 
 
 
 
###chr9
testChr = 9
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr9
stPhyPos = 0; edPhyPos = 2000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 3.5 # add this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#the end of chr9
stPhyPos = 89000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 0.5 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
###chr10
testChr = 10
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr10
stPhyPos = 0; edPhyPos = 500000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 1.6 # add this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
#the end of chr10
stPhyPos = 117000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 4.5 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# set this based on old map position. 
markerMapRes["NGHIR009559066","newMapPos"] = 0
 
 
###chr11; problem at the 38- 50 
testChr = 11
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr11; 2 markers set as 0
markerMapRes[c("NGHIR009560556","NGHIR009561428"),"newMapPos"] = 0
 
#the end of chr11
markerMapRes[c("NGHIR009575887"),"newMapPos"] = 151.7
 
 
 
###chr12
testChr = 12
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr12
addDist = 1 # add this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
 
# NG0204400 & NG0206798
markerMapRes[c("NG0204400","NG0206798"),"newMapPos"] = 0
 
 
# end of chr12; 
stPhyPos = 108800000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
 
markerMapRes[c("NGHIR009549524","NGHIR009550040","NGHIR009550227"),"newMapPos"] = 165.9 # set them as 164.9
 
 
###chr13
testChr = 13
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr13; A few chr13 markers should be linked to chr26, and these markers messed up the linkage calls
stPhyPos = 0; edPhyPos = 1500000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = -50.54 + 15.6 # subtract this from the chr; A few markers that are questionable make the incorrect calls; Adjusted. 
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = NA
testData$newMapPos[1] = 0
testData$newMapPos[nrow(testData)] = 15.6
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
#the end of chr13
stPhyPos = 111800000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
markerMapRes[c("NGHIR009557434","NG0203817","NGHIR009556892","NGHIR009556570","NG0204842"),"newMapPos"] = 147.5 # set them as 0
 
### chr14
testChr = 14
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr14
stPhyPos = 0; edPhyPos = 500000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 0.1 # add this chr based on old map
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# end of chr14
markerMapRes[c("NGHIR009535134"),"newMapPos"] = 117.0 # set them based on oldmap
 
 
### chr15
testChr = 15
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr15
stPhyPos = 0; edPhyPos = 1100000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 11.4 # add this chr based on old map
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
### chr16
testChr = 16
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr16
markerMapRes["NGHIR009581578","newMapPos"] = 0 # set the last one 
 
# end of chr16
stPhyPos = 54000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
 
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 2 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
 
### chr17
testChr = 17
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr17
stPhyPos = 0; edPhyPos = 700000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 0.2 # add this chr based on old map
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# end of chr17
markerMapRes[c("NGHIR009580391","NGHIR009538581","NGHIR009537689"),"newMapPos"] = 122.6# set the end of 124.37
 
 
### chr18
testChr = 18
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr18
markerMapRes[c("NGHIR009580580","NGHIR009565869","NGHIR009538049"),"newMapPos"] = 0
 
# end of chr18
stPhyPos = 66200000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
 
# no data and use the old map distance
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 14 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
### chr19
testChr = 19
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr19 
stPhyPos = 0; edPhyPos = 500000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 0.3 # add this chr
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# end of chr19
stPhyPos = 65200000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 4.5 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# another 1 markers without physicla position. 
 
subMarkerMap = markerMapRes %>% 
  filter(oldMapChr == testChr) %>% 
  arrange(oldMapChr,oldMapPos)
 
intPos = interpolate(subMarkerMap$oldMapPos,subMarkerMap$newMapPos)[,3]
table(is.na(subMarkerMap$newMapPos))
table(intPos == subMarkerMap$newMapPos)
markerMapRes[subMarkerMap$Q_ID,"newMapPos"] = intPos
 
### chr20
testChr = 20
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr20
markerMapRes[c("NGHIR009539758"),"newMapPos"]  = 0
 
 
# end of chr20
stPhyPos = 55000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 42 # add cM to the end
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# another 2 markers without physicla position. 
 
subMarkerMap = markerMapRes %>% 
  filter(oldMapChr == testChr) %>% 
  arrange(oldMapChr,oldMapPos)
 
intPos = interpolate(subMarkerMap$oldMapPos,subMarkerMap$newMapPos)[,3]
table(is.na(subMarkerMap$newMapPos))
table(intPos == subMarkerMap$newMapPos)
markerMapRes[subMarkerMap$Q_ID,"newMapPos"] = intPos
 
 
### chr21
testChr = 21
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr21
markerMapRes[c("NGHIR009570442"),"newMapPos"]  = 0
 
# end of chr21
stPhyPos = 71000000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 1 # add cM to the end based old map
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
### chr22
testChr = 22
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr22
stPhyPos = 0; edPhyPos = 3000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 25 # add this chr based on old map
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# end of chr22
stPhyPos = 53600000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 6 # add cM to the end based on old map
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
### chr23
testChr = 23
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr23
stPhyPos = 0; edPhyPos = 2500000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 8 # add this chr based on both
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# end of chr23
stPhyPos = 68900000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
 
testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 2.2 # add cM to the end, based on old map
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# another 1 markers without physicla position. 
 
subMarkerMap = markerMapRes %>% 
  filter(oldMapChr == testChr) %>% 
  arrange(oldMapChr,oldMapPos)
 
intPos = interpolate(subMarkerMap$oldMapPos,subMarkerMap$newMapPos)[,3]
table(is.na(subMarkerMap$newMapPos))
table(intPos == subMarkerMap$newMapPos)
markerMapRes[subMarkerMap$Q_ID,"newMapPos"] = intPos
 
 
### chr24
testChr = 24
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr24
stPhyPos = 0; edPhyPos = 500000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 1 # add this chr based on proportion
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# end of chr24
stPhyPos = 73300000; edPhyPos = 200000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist)
#write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))
 
markerMapRes[c("NG0204897","NG0210714","NGHIR009544836"),"newMapPos"] =  max(testData$newMapPos,na.rm=T) + 0.2 # add cM to the end
 
### chr25
testChr = 25
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
# start of chr25
stPhyPos = 0; edPhyPos = 4000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 31 # add this chr based on proportion
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
 
 
### chr26
testChr = 26
chrMks = subset(markerMapRes,chr == testChr)$Q_ID
 
#start of chr26
stPhyPos = 0; edPhyPos = 4000000 # check their linkage with more markers
testData = markerMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist = bind_cols(as.data.frame(markerMapRes[colnames(testMkPwDist),"newMapPos"]), testMkPwDist) 
#write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
 
addDist = 12.5 # add this chr based on proportion
markerMapRes[chrMks,"newMapPos"] = markerMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
markerMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
 
# another markers without physicla position. 
 
subMarkerMap = markerMapRes %>% 
  filter(oldMapChr == testChr) %>% 
  arrange(oldMapChr,oldMapPos)
 
intPos = interpolate(subMarkerMap$oldMapPos,subMarkerMap$newMapPos)[,3]
table(is.na(subMarkerMap$newMapPos))
table(intPos == subMarkerMap$newMapPos)
markerMapRes[subMarkerMap$Q_ID,"newMapPos"] = intPos
 
 
## add oldChr to markers without phy position
table(markerMapRes$chr,markerMapRes$oldMapChr)
markerMapRes$phyChr = markerMapRes$chr
mkIndex = which(is.na(markerMapRes$chr))
markerMapRes$chr[mkIndex] = markerMapRes$oldMapChr[mkIndex]
 
 
# add clr positions
genMap_clr = read.csv("/mnt/cotton/cottonMarkerMapRes_clr.csv")
row.names(genMap_clr) = genMap_clr$Q_ID
 
markerMapRes$genMapChr_clr = genMap_clr[markerMapRes$Q_ID,"genMapChr"]
markerMapRes$genMapPos_clr = genMap_clr[markerMapRes$Q_ID,"newMapPos"]
 
markerMapRes$pbMks = NA
 
# markers that have suspected linkage info; 1) oldMapChr != chr 2) from validation data
extraPbMks = c("NG0207810","NGHIR009559655","NGHIR009574723","NG0203344",	"NGHIR009568859","NGHIR009581329","NG0203423",
               "NG0210277","NGHIR009582419","NGHIR009570036","NGHIR009570362","NGHIR009582420", 
               "NGHIR009547662", "NGHIR009548704", "NGHIR009554065", "NGHIR009556969", "NGHIR009541628", "NGHIR009559230",
               "NGHIR009557628", "NGHIR009569076", "NGHIR009564449","NGHIR009564031", "NGHIR009564879", "NGHIR009548352",
               "NGHIR009575436", "NGHIR009567438", "NGHIR009579092", "NGHIR009576235", "NGHIR009570572", "NGHIR009542065")
 
markerMapRes[c(pbMarkers,extraPbMks),"pbMks"] = "Yes"
 
 
# chr11 middle has some issue
# incorrect call at stPhyPos = 6730002; edPhyPos = 112600000 # check their linkage with more markers
# 113 and 112 has large gap, remove it. 
 
testChr = 11
testData = markerMapRes %>% filter(chr == testChr) %>% arrange(chr,newMapPos)
tmpMapPos = testData$genMapPos_clr 
table(tmpMapPos[2:nrow(testData)] - tmpMapPos[1:(nrow(testData)-1)] >= 0)
stPos = 0
for(i in 1:(nrow(testData)-1)){
  if(is.na(tmpMapPos[i])){next}
  if(tmpMapPos[i] < stPos){
    tmpMapPos[i] <- NA
  }else{
    if(!is.na(tmpMapPos[i+1])){
      if(tmpMapPos[i] > tmpMapPos[i+1]){
        tmpMapPos[i] <- NA
      }else{
        stPos = tmpMapPos[i]
      }
    }
    
  }
}
 
indexInt = which(!is.na(testData$Sbegin))
testData$newMapPos[indexInt] = interpolate(testData$Sbegin[indexInt],tmpMapPos[indexInt])[,3]
markerMapRes[testData$Q_ID,"newMapPos"] = testData$newMapPos
 
#### check ####
newMapSummary = as.data.frame(markerMapRes %>% filter(!is.na(markerMapRes$chr)) %>% group_by(chr) %>% 
                                summarise(maxLenNew = max(newMapPos,na.rm=T),maxLenOld = max(oldMapPos,na.rm=T),
                                          minLenNew = min(newMapPos,na.rm=T),minLenOld = min(oldMapPos,na.rm=T)))
chrStPos = newMapSummary$minLenNew
 
# set chr start to 0 if not yet
for(i in 1:nrow(newMapSummary)){
  testChr = newMapSummary$chr[i]
  minChrLen = newMapSummary$minLenNew[i]
  if(minChrLen != 0){
    chrMkIndex = which(markerMapRes$chr == testChr)
    markerMapRes$newMapPos[chrMkIndex] = markerMapRes$newMapPos[chrMkIndex] - minChrLen
  }
}
 
# round mapPos to 0.1 & put 1/10^6 into the same bin
markerMapRes = markerMapRes %>% arrange(chr,newMapPos,Sbegin)
markerMapRes$newMapPos2 = round(markerMapRes$newMapPos,1)
 
for(chrN in names(allChrs)){
  tmp  = markerMapRes %>% filter(chr == chrN)
  for(tmpInt in unique(tmp$newMapPos2)){
    tmpIndex = which(tmp$newMapPos2 == tmpInt)
    if(length(tmpIndex) > 0){
      tmp$newMapPos2[tmpIndex] = tmpInt + seq(1,length(tmpIndex),1) * 10^(-6)
    }
  }
  
  markerMapRes[row.names(tmp),"genMapPos_2023"] = tmp$newMapPos2
  
}
 
markerMapRes$genMapChr = allChrs[markerMapRes$chr]
 
 
table(markerMapRes$genMapPos_2023[2:nrow(markerMapRes)] - markerMapRes$genMapPos_2023[1:(nrow(markerMapRes)-1)] < 0)
 
table(markerMapRes$chr)
 
# add another three missed based on markerData
 
cropFolder = "/mnt/cotton"
write.csv(markerMapRes,paste0(cropFolder,"/cottonMarkerMapRes_hifi.csv"),row.names=F)
 
ggplot(markerMapRes %>% filter(!is.na(chr)),aes(Sbegin/10^6,newMapPos,col=as.character(oldMapChr))) + geom_point() + facet_wrap(vars(chr),scale="free") + theme_bw() +
  ylab("newGenMap (2023_HIFI)") + xlab("phyPos (DP393_HIFI)") + ggtitle("CottonGeneticMap (2023_HIFI) vs PhyMap(DP393_HIFI)")
 
ggplot(markerMapRes %>% filter(!is.na(chr)),aes(newMapPos,Sbegin/10^6,col=as.character(oldMapChr))) + geom_point() + facet_wrap(vars(chr),scale="free") + theme_bw() +
  xlab("newGenMap (2023_HIFI)") + ylab("phyPos (DP393_HIFI)") + ggtitle("CottonGeneticMap (2023_HIFI) vs PhyMap(DP393_HIFI)")
 
ggplot(markerMapRes %>% filter(!is.na(chr)),aes(oldMapPos,Sbegin/10^6,col=as.character(oldMapChr))) + geom_point() + facet_wrap(vars(chr),scale="free") + theme_bw() +
  xlab("oldGenMap (xxxx)") + ylab("phyPos (DP393_HIFI)") + ggtitle("CottonGeneticMap (xxxx) vs PhyMap(DP393_HIFI)")
 
ggplot(markerMapRes %>% filter(!is.na(chr)),aes(newMapPos,oldMapPos,col=as.character(oldMapChr))) + geom_point() + facet_wrap(vars(chr),scale="free") + theme_bw() +
  ylab("oldGenMap (xxxx)") + xlab("newGenMap (2023_HIFI)") + ggtitle("CottonGeneticMap (xxxx) vs CottonGeneticMap(2023_HIFI)")
 
ggplot(markerMapRes %>% filter(!(is.na(chr_clr) | chr_clr == "")),aes(genMapPos_2023,genMapPos_clr,col=chr_clr)) + geom_point() + facet_wrap(vars(chr),scale="free") + theme_bw() +
  xlab("genMapPos (HIFI)") + ylab("genMapPos (CLR)") + ggtitle("HIFI vs CLR (DP393 genMap)")
 
tmp = markerMapRes %>% group_by(chr) %>% summarise(chrLenOld = max(oldMapPos,na.rm=T), chrLenNew = max(newMapPos,na.rm=T))
 
 
 
 
 
 
#### validate - dist check ####
tmpMarkers = colnames(totalRecRate)
mkLinkageInfo2$newMapPos = markerMapRes[mkLinkageInfo2$Q_ID,"newMapPos"]
 
for(chrN in sort(unique(markerMapRes$chr))){
  for(pos in seq(10,200,by=20)){
    if(pos >= max(subset(mkLinkageInfo2,chr == chrN)$newMapPos,na.rm=T) - 10){next}
    subDf = mkLinkageInfo2 %>% filter(chr == chrN, newMapPos >= pos - 2, newMapPos <= pos + 2,!is.na(totalMkLinkedOwnChr))
    if(nrow(subDf) >= 1){
      mk = subDf$Q_ID[which.max(subDf$totalMkLinkedOwnChr)]
      recData2 = plotSingleMarkerRecRateByNewMap(which(tmpMarkers == mk),chrs = chrN,mkPosMap = markerMapRes[tmpMarkers,])
      # Customizing the output
      pdf(paste0(chrN,"_",pos,"cM.pdf"),   width = 12, height = 7)
      print(recData2)
      dev.off() 
    }
    
  }
}
 
 
 
###
