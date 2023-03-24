
library(parallel)
library(MASS)
library(foreach)
library(doParallel)
library(tidyverse)
library(lattice)

source("/mnt/repo/trio_genetic_map/code/trio_mapping_functions.R")

#### Soy Genetic Map ####
# The purpose of this project is to denovo construct a soy consensus genetic map, by leveraging HiFi A3555 genome and Trio mapping;
# 1. Markers were mapped to A3555 genome to locate their positions and orders (Done by Alex Brohammer (GMAP) and Joe Zhou)
# 2. Trio data was collected from Danny's Soy_FP_Parental_Haplotypes; Criterial was applied to filter Trios
# 3. FpString was collected by Danny and saved at https://domino.science-at-scale.io/u/dchal/pullFpGenos/browse/Soybean__MON_2010_v2r1


#### 1. fpString data ####
# read fpString (fpGeno) data (FpString was collected by Danny and saved at https://domino.science-at-scale.io/u/dchal/pullFpGenos/browse/Soybean__MON_2010_v2r1)
fpString = read.delim(gzfile("/domino/datasets/Soybean__MON_2010_v2r1/Soybean__MON_2010_v2r1_extramarkers__withData_LshFpGenos.gz"),header=F)
fpLines = fpString$V1
row.names(fpString) = fpString$V1


#### 2. Soy trio data ####
# get soybean trio information
# ls /domino/datasets/Soy_FP_Parental_Haplotypes/metaData/ >  /mnt/data/metaData_soyFP.csv
# cat * > /mnt/data/metaData_allLines.txt 

df <- read.delim("/mnt/data/metaData_allLines.txt",sep=",")
index <- seq(1, nrow(df),by=2)
df <- df[index,] # the even number rows are just header. Removed; 
row.names(df) <- df$pedigree
df$numOfPrt <- sapply(df$parents,function(x){return(length(unlist(strsplit(x,";"))))})
df$allPrtFp <- sapply(df$parents,function(x){return(all(unlist(strsplit(x,";")) %in% fpLines))})
table(df$allPrtFp,df$numOfPrt)
df$rawPctContribPrt1 <- as.numeric(sapply(df$rawPctContrib,function(x){return(unlist(strsplit(x,";"))[1])}))
df$rawPctContribPrt2 <- as.numeric(sapply(df$rawPctContrib,function(x){return(unlist(strsplit(x,";"))[2])}))
df <- df %>% filter(numOfPrt == 2) # only use bi-cross 
df$prt1 <- sapply(df$parents,function(x){unlist(strsplit(x,";"))[1]})
df$prt2 <- sapply(df$parents,function(x){unlist(strsplit(x,";"))[2]})
table(df$rawPctContribPrt1 > 35, df$rawPctContribPrt2 > 35)
ggplot(data=df %>% filter(rawPctOther < 1),aes(rawPctContribPrt1,rawPctContribPrt2,col=as.numeric(rawPctOther))) + geom_point() + theme_bw()
table(str_detect(df$originUsed,"\\*"))

# remove backgross & lines with excessive rawPctOthers
df <- df %>% filter(!str_detect(originUsed,"\\*"), rawPctOther <= 1)

# only keep lines with both parents FP'ed
df <- df %>% filter(allPrtFp)

# read marker data
markerMap <- read.table("/domino/datasets/Soybean__MON_2010_v2r1/Soybean__MON_2010_v2r1_extramarkers__Map.txt",header=T,stringsAsFactors = F)
row.names(markerMap) <- markerMap$markerName

##### 2. pull fpString ####
# fpDataVec; row as marker and col as line
# takes about 3 mins to generate the fpMatrix; 67717 markers x 15634 lines


fpDataVect = convertFpStringToGeno(df,fpString,markerMap)
#saveRDS(fpDataVec,"/mnt/intermediateData/fpDataVec_fpGeno.rds") # decide not save since it is faster to generate data

#### 3. calculate pairsise recommbination data #####
# get informative marker status by each line

tmpFolder = "/mnt/tmpData/fpGeno"
infoMk = infoMarkerByLine(df,fpDataVec,nCores,tmpFolder)

# get all markers that have recFraction data
# for any pair of markers, count #totalTros and #recombinedTrios
# Output are two matrix, each one is a large matrix with the dimension as #markers x #markers (around 65k x 65k)

outFileFolder = "/mnt/results/20221024"
runRecRate = generateRecData(tmpFolder,markerMap,outFileFolder,verbose = TRUE)


#### 4. check linked markers ####
# 1. OffChr markers
# 2. OffPos markers
# 3. marker clusters (use 10 markers/window); Get the 5 cM as accurate as possible
#1. Calculate map distance within clusters
#2. Calculate map distance between clusters
#3. Define boarders of every 5 cM (staring from teh first clusters)
outFileFolder = "/mnt/results/20221024"
totalTrioCount <- readRDS(paste0(outFileFolder,"/totalTrioCount_allChr.rds"))
totalRecRate <- readRDS(paste0(outFileFolder,"/trioRecRate_allChr.rds"))

# physical position (A3555)
# 55,292 Infinium markers against A3555 v3 (Done by Alex Brohammer)

phyPos1 = read.delim("/mnt/data/phyPos/gmap_parsed_A3555_part1.txt",header=T) %>% select(-c("NumUnknowns","ConsDiffs"))
phyPos2 = read.delim("/mnt/data/phyPos/Soy_SNP2018_gmap_parsed.txt",header=T) %>% select(-c("NumUnknowns","ConsDiffs"))
phyPos2$Desc = NA # to make this consistent with phyPos1$Desc
phyPos = unique(bind_rows(phyPos1,phyPos2))

phyPos = phyPos %>% 
  mutate(chr = as.numeric(str_replace_all(Subj,"GLYMA_A3555_v3_0_0_chr",""))) %>% 
  filter(Coverage >= 90, PercentID >= 90, chr <= 1000000) %>% 
  arrange(chr,Sbegin)

table(phyPos$chr)

# non-unique mapped markers
mulPosMk = phyPos$Q_ID[duplicated(phyPos$Q_ID)]

# overlapped markers
ovlpMk = phyPos$Q_ID[phyPos$Q_ID %in% row.names(totalTrioCount) & (!phyPos$Q_ID %in% mulPosMk)]

# read marker data
tmpPhyPos = phyPos %>% filter(!(Q_ID %in% mulPosMk))
row.names(tmpPhyPos) = tmpPhyPos$Q_ID
markerMap$phyChr = tmpPhyPos[markerMap$markerName,"chr"]

totalTrioCount = totalTrioCount[ovlpMk,ovlpMk]
totalRecRate = totalRecRate[ovlpMk,ovlpMk]

#
## first check linked markers 

#res = getLinkedMarkerInfo(phyPos,ovlpMk, totalTrioCount,totalRecRate, markerMap,minCount=100, maxRecRate=0.1)
#res$chrRatio = round(res$totalMkLinkedOwnChr / res$totalMkLinked,2)
#saveRDS(res,paste0("/mnt/results/20221024/marker_linkage_distribution_",nrow(res),".rds"))

res = readRDS("/mnt/results/20221024/marker_linkage_distribution_50340.rds")

# find markers that are suspected; likely linked into a wrong group

pbMkData = res %>% filter((totalMkLinkedOwnChr == 1 & numLinkedChr >= 2) | chrRatio < 0.5) # a singleton or less linked with own chr 
pbMarkers = pbMkData$Q_ID
#res_2 = getLinkedMarkerInfo(phyPos,ovlpMk[!ovlpMk %in% pbMarkers], totalTrioCount, totalRecRate, markerMap,minCount=100, maxRecRate=0.1)
#write.csv(res_2[,c(1:8,10:15)],paste0("/mnt/results/20221024/marker_linkage_distribution_",nrow(res_2),".csv"))
res_2 = read.csv("/mnt/results/20221024/marker_linkage_distribution_50284.csv")
sort(table(res_2$chrDist),decreasing = T)

#### 4.1 run functions to get anchor distance ####
# for(chrN in 1:20){
#   print(paste(chrN, Sys.time()))
#   findAnchorsByChr(chrN,res_2,totalTrioCount,totalRecRate)
# }

rm(totalTrioCount,totalRecRate)
gc()

# plot 
anchorResAllChrs = data.frame()
for(testChr in 1:20){
  tmp = read.csv(paste0("/mnt/tmpData/tmpAnchorData/anchorRes_chr",testChr,".csv"))
  tmp$chr = testChr
  anchorResAllChrs = bind_rows(anchorResAllChrs,tmp)
}

anchorResAllChrs$maxAnchorDist = sapply(anchorResAllChrs$distVec,function(x){max(as.numeric(unlist(strsplit(x,"_"))),na.rm=T)})

ggplot(anchorResAllChrs,aes(as.character(chr),totalGenLen,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),MAE,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),MAE_pct,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),maxAnchorAE,col=as.character(chr))) + geom_boxplot() + theme_bw()
ggplot(anchorResAllChrs,aes(as.character(chr),maxAnchorDist,col=as.character(chr))) + geom_boxplot() + theme_bw()



#### 5. select anchor path ####
# request at least 10 groups, with chr5 as a exception as the max is 9; 
selAnchorPath = data.frame()
for(testChr in 1:20){
  tmp = read.csv(paste0("/mnt/tmpData/tmpAnchorData/anchorRes_chr",testChr,".csv"))
  tmp$chr = testChr
  tmp = tmp %>% filter(numOfGroups >= 10 | numOfGroups == max(numOfGroups)) %>% arrange(MAE_pct)
  selAnchorPath = bind_rows(selAnchorPath,tmp[1,])
}

selAnchorPath$maxAnchorDist = sapply(selAnchorPath$distVec,function(x){max(as.numeric(unlist(strsplit(x,"_"))),na.rm=T)})


#### 6. all group distances  ####
# step1: map high quality anchors (groups) onto the map #selGroupSm = groupSm %>% filter(numOfData >= 80, minDataPerMk >= 5, max <= 0.5)
# step2: map sub-high quality anchors onto the map (particually those present in large gaps)
# step3: map individual markers onto the map (because of the slightly fluctuation, decide to replace this step by finding anchor markers at every 1 cM distance)



markerMapRes = data.frame()
minDataWithAnchor = 5 # the minumum data point between a marker and a anchor (10 markers); 

oldNewChrDict = c(8,16,6,18,7,9,11,3,19,17,4,10,13,20,2,1,12,5,14,15) # old map chr vs phy map chr

markerMapRes = data.frame()
largestGaps = rep(NA,20)

for(testChr in 1:20){
  print(paste(testChr,Sys.time()))
  selAnchors = as.numeric(unlist(strsplit(selAnchorPath$selIndex[testChr],"_")))
  distVect = as.numeric(unlist(strsplit(selAnchorPath$distVec[testChr],"_")))
  pwGroupSm = readRDS(paste0("/mnt/tmpData/tmpAnchorData/pwGroupSm_chr",testChr,".rds"))
  anchorRes = read.csv(paste0("/mnt/tmpData/tmpAnchorData/anchorRes_chr",testChr,".csv"))
  withGroupSm = read.csv(paste0("/mnt/tmpData/tmpAnchorData/groupSm_chr",testChr,".csv"))
  goodAnchors = withGroupSm %>% filter(numOfData >= 80, minDataPerMk >= 5, max <= 0.5) # adjust the threshold to 60_3_0.5 does not reduce big gaps,eg. chr2
  chrMapDist = readRDS(paste0("/mnt/tmpData/tmpAnchorData/subChrMapDist_chr",testChr,".rds"))
  # get all anchor positions
  allGroupDist  = goodAnchorDist(goodAnchors,selAnchors,anchorRes,pwGroupSm)
  goodAnchors$groupPos = round(allGroupDist,3)
  goodAnchors = goodAnchors %>% filter(!is.na(groupPos))
  distGap = goodAnchors$groupPos[2:nrow(goodAnchors)] - goodAnchors$groupPos[1:(nrow(goodAnchors) - 1)]
  
  # remove anchors that are in wrong oder, unless there is a big gap
  # as checked, the distGap in all 20 chrs are less than 0.66; after adjusted, no rong order
  # so it is ok to just remove the anchors.
  
  while(any(distGap < 0) & min(distGap) > -1){
    goodAnchors = goodAnchors[-(which(distGap < 0) + 1),]
    distGap = goodAnchors$groupPos[2:nrow(goodAnchors)] - goodAnchors$groupPos[1:(nrow(goodAnchors) - 1)]
  }
  
  # interpolate markers onto the map
  subChrData = read.csv(paste0("/mnt/tmpData/tmpAnchorData/subChrData_chr",testChr,".csv"))
  row.names(subChrData) = subChrData$Q_ID
  newMarkerRes = getMarkerPos(subChrData,goodAnchors)
  
  # re-adjust marker position at big gaps (centromere regions are usually wrongly interpolated, eg. chr2 15-30 Mb
  bigGapIndex = which(distGap >= 5)
  for(index in bigGapIndex){
    a1Index = goodAnchors$gIndex[index] + 10
    a2Index = goodAnchors$gIndex[index+1] - 1
    tmpPos = reCalMkDist(goodAnchors,index,chrMapDist,newMarkerRes)
    if(!all(is.na(tmpPos))){
      newMarkerRes$newMapPos[a1Index:a2Index] = NA
      newMarkerRes$newMapPos[as.numeric(names(tmpPos))] = tmpPos
      newGaps = diff(c(goodAnchors$groupPos[index],tmpPos,goodAnchors$groupPos[index+1]))
      if(any(newGaps >= 5)){print(c(goodAnchors$groupPos[index],tmpPos,goodAnchors$groupPos[index+1]))}
    }
  }
  
  # interpolate all other markers with A3555 physical positions
  subPhyPos = phyPos %>% filter(chr == testChr) %>% arrange(Sbegin)
  subPhyPos = subPhyPos[!(duplicated(subPhyPos$Q_ID)),]
  row.names(subPhyPos) = subPhyPos$Q_ID
  subPhyPos$newMapPos = newMarkerRes[subPhyPos$Q_ID,"newMapPos"]
  intPosPhy = interpolate(subPhyPos$Sbegin,subPhyPos$newMapPos)
  subPhyPos$newMapPos = intPosPhy[,3]
  
  # interpolation based old gen map for markers without A3555 positions
  subMarkerMap = markerMap %>% 
    filter(genMapChr == oldNewChrDict[testChr]) %>% 
    select(markerName, genMapChr,genMapPos) %>% 
    arrange(genMapChr,genMapPos)
  subMarkerMap$newMapPos = subPhyPos[subMarkerMap$markerName,"newMapPos"]
  intPosGen = interpolate(subMarkerMap$genMapPos,subMarkerMap$newMapPos)
  subMarkerMap$newMapPos = intPosGen[,3]
  extraMarkers = subMarkerMap$markerName[!subMarkerMap$markerName %in% phyPos$Q_ID]
  
  # combine phyPos,oldGenMap together
  subPhyPos[extraMarkers,] = NA
  subPhyPos[extraMarkers,"Q_ID"] = extraMarkers
  subPhyPos$oldMapChr = subMarkerMap[subPhyPos$Q_ID,"genMapChr"]
  subPhyPos$oldMapPos = subMarkerMap[subPhyPos$Q_ID,"genMapPos"]
  subPhyPos[extraMarkers,"newMapPos"] = subMarkerMap[extraMarkers,"newMapPos"]
  subPhyPos  = subPhyPos %>% arrange(newMapPos)
  
  # combin all chr data together
  markerMapRes = bind_rows(markerMapRes,subPhyPos)
  #print(table(subPhyPos$newMapPos[2:nrow(subPhyPos)] - subPhyPos$newMapPos[1:(nrow(subPhyPos)-1)] < 0))
  largestGaps[testChr] = max(subPhyPos$newMapPos[2:nrow(subPhyPos)] - subPhyPos$newMapPos[1:(nrow(subPhyPos)-1)],na.rm=T)
  
}

markerMapRes$oldMapChr = markerMap[markerMapRes$Q_ID,"genMapChr"]
markerMapRes$oldMapPos = markerMap[markerMapRes$Q_ID,"genMapPos"]

write.csv(markerMapRes,"/mnt/soybean/results/markerMapRes_v2.csv",row.names=F)

ggplot(markerMapRes %>% filter(!is.na(chr)),aes(Sbegin/10^6,newMapPos,col=as.character(chr))) + geom_point() + facet_wrap(vars(chr),scale="free") + theme_bw() +
  ylab("newGenMap (2023)") + xlab("phyPos (A3555 v3)") + ggtitle("SoyGeneticMap (2023) vs PhyMap(A3555)")

ggplot(markerMapRes %>% filter(!is.na(chr)),aes(newMapPos,Sbegin/10^6,col=as.character(chr))) + geom_point() + facet_wrap(vars(chr),scale="free") + theme_bw() +
  xlab("newGenMap (2023)") + ylab("phyPos (A3555 v3)") + ggtitle("SoyGeneticMap (2023) vs PhyMap(A3555)")

ggplot(markerMapRes %>%  filter(!is.na(chr)),aes(newMapPos,oldMapPos,col=as.character(oldMapChr))) + geom_point() + 
  facet_wrap(vars(chr),scale="free") + theme_bw() + xlab("newGenMap (2023)") + ylab(" oldGenMap (2010)") + 
  ggtitle("SoyGeneticMap (2023) vs SoyGeneticMap (2010)")

ggplot(markerMapRes %>% filter(!is.na(chr)),aes(Sbegin/10^6,oldMapPos,col=as.character(oldMapChr))) + geom_point() + 
  facet_wrap(vars(chr),scale="free") + theme_bw() + xlab("phyPos (A3555 v3)") + ylab(" oldGenMap (2010)")



plot(subChrData$newMapPos,subChrData$Sbegin)
plot(subChrData$Sbegin,subChrData$newMapPos)

#### 6.1 manual fill in NA ####

oldNewChrDict = c(8,16,6,18,7,9,11,3,19,17,4,10,13,20,2,1,12,5,14,15) # old map chr vs phy map chr

newMapRes = read.csv("/mnt/soybean/results/markerMapRes_v2.csv")
row.names(newMapRes) = newMapRes$Q_ID

res_2$newMapPos = newMapRes[res_2$Q_ID,"newMapPos"]



# chr is NA; map oldMapChr to newMapChr; These marker positions are interpolated based on old map; 
chrNaIndex = which(is.na(newMapRes$chr) & !(is.na(newMapRes$newMapPos)) & !(is.na(newMapRes$oldMapChr)))
newOldChr = 1:20
names(newOldChr) = oldNewChrDict
newMapRes[chrNaIndex,"chr"] = newOldChr[as.character(newMapRes[chrNaIndex,"oldMapChr"])]




#chr1; five with phyPos but no genPos
testChr = 1;
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))
stPhyPos = 0; edPhyPos = 300000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = getPwKsbDist(testData,totalTrioCount, totalRecRate)
max(testMkPwDist,na.rm=T)
#since the max dist is 0, so set everyone to 0
newMapRes[testMk,"newMapPos"] = 0

#chr2;  
testChr = 2; 
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr2
stPhyPos = 0; edPhyPos = 500000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] # about 0.75 from NA to 0

newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + 0.75 # add 0.75 to every chr2 markers
testData$newMapPos = testData$newMapPos + 0.75
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr2
stPhyPos = 54700000; edPhyPos = 55000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] # about 0.2 to the end
testData$newMapPos[nrow(testData)] = 105.96 + 0.2
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr3;  
testChr = 3; 
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr3
stPhyPos = 0; edPhyPos = 300000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 1.4 # add 1.4 to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr3
stPhyPos = 47600000; edPhyPos = 49000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 3 # add 3 cM
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr4;  
testChr = 4
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))
# 4 markers only and only 2 have recRate data. Manually assign
newMapRes["NS0099895","newMapPos"] = 0
newMapRes[testMk[2:4],"newMapPos"] = 83.095
###

#chr5;  
testChr = 5; 
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr5
stPhyPos = 0; edPhyPos = 1000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 3.17 # add 3.17 to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr5
stPhyPos = 42990000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 6.7 # add 6.7 cM to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr6;  
testChr = 6; 
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr6
stPhyPos = 0; edPhyPos = 500000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 2 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr6
stPhyPos = 52450000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 3.7 # add 6.7 cM to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr7;  
testChr = 7
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr7
stPhyPos = 0; edPhyPos = 2710000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 11 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr7
stPhyPos = 48300000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 0.35 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr8;  
testChr = 8
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr8
stPhyPos = 0; edPhyPos = 1500000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 6.6 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr8
stPhyPos = 50100000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 0.2 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr9;  
testChr = 9
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr9
stPhyPos = 0; edPhyPos = 500000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 2.4 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr9
stPhyPos = 50800000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 0.4 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr10  
testChr = 10
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr10
stPhyPos = 0; edPhyPos = 200000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 0.2 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr10
stPhyPos = 53700000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 3.4 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]


###

#chr11  
testChr = 11
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))
newMapRes[testMk[1:5],"newMapPos"] = 0
newMapRes[testMk[6:length(testMk)],"newMapPos"] = 87.949 # 0 recombination
###

#chr12
testChr = 12
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

#end of chr12
stPhyPos = 42400000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 6.8 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr13  
testChr = 13
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr13
stPhyPos = 0; edPhyPos = 3000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 0.5 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr13
stPhyPos = 45700000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 0.2 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr14
testChr = 14
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr14
stPhyPos = 0; edPhyPos = 800000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 2.2 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr14
stPhyPos = 52300000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 1.5 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr15
testChr = 15
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr15
stPhyPos = 0; edPhyPos = 250000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 0.15 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr15
newMapRes[testMk[16:length(testMk)],"newMapPos"] = 77.915 + addDist # 0 recombination

###
#chr16
testChr = 16
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr16
stPhyPos = 0; edPhyPos = 200000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 0.2 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr17
testChr = 17
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr17
stPhyPos = 0; edPhyPos = 300000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 0.3 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr17
stPhyPos = 43350000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 0.5 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr18
testChr = 18
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr18
newMapRes["NGMAX007935888","newMapPos"] = 0

#end of chr18
stPhyPos = 60270000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 2.2 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr19
testChr = 19
chrMks = subset(newMapRes,chr == testChr)$Q_ID
testMk = subset(newMapRes,chr == testChr & is.na(newMapPos))$Q_ID
table(testMk %in% colnames(totalTrioCount))

# start of chr19
stPhyPos = 0; edPhyPos = 100000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_st.csv"))
addDist = 0 # add this amount to this chr
newMapRes[chrMks,"newMapPos"] = newMapRes[chrMks,"newMapPos"] + addDist
testData$newMapPos = testData$newMapPos + addDist
testData$newMapPos[1] = 0
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]

#end of chr19
stPhyPos = 53000000; edPhyPos = 100000000 # check their linkage with more markers
testData = newMapRes %>% filter(chr == testChr, Sbegin >= stPhyPos, Sbegin <= edPhyPos) %>% arrange(chr,Sbegin)
testMkPwDist = as.data.frame(getPwKsbDist(testData,totalTrioCount, totalRecRate))
max(testMkPwDist,na.rm=T)
testMkPwDist[,"newMapPos"] = newMapRes[colnames(testMkPwDist),"newMapPos"] 
write.csv(testMkPwDist,paste0("chr",testChr,"_ed.csv"))

testData$newMapPos[nrow(testData)] = max(testData$newMapPos,na.rm=T) + 4 # add this much to the end
newMapRes[testData$Q_ID,"newMapPos"] = interpolate(testData$Sbegin,testData$newMapPos)[,3]
###

#chr20; No missing



# interpolation based old gen map for last 161 markers
for(chrN in 1:20){
  tmp = newMapRes %>% 
    filter(oldMapChr == chrN, chr == newOldChr[as.character(chrN)] | is.na(chr)) %>% 
    select(newMapPos,oldMapPos,chr,oldMapChr) %>% arrange(oldMapPos)
  if(is.na(tmp$newMapPos[1])){tmp$newMapPos[1] = min(tmp$newMapPos,na.rm=T)}
  if(is.na(tmp$newMapPos[nrow(tmp)])){tmp$newMapPos[nrow(tmp)] = max(tmp$newMapPos,na.rm=T)}
  naIndex = which(is.na(tmp$newMapPos))
  if(length(naIndex) > 0){
    tmp$newMapPos = interpolate(tmp$oldMapPos,tmp$newMapPos)[,3]
  }
  newMapRes[row.names(tmp),"chr"] = newOldChr[as.character(chrN)]
  newMapRes[row.names(tmp),"newMapPos"] = tmp$newMapPos
}

# round mapPos to 0.1 & put 1/10^6 into the same bin
newMapRes = newMapRes %>% arrange(chr,newMapPos,Sbegin)
newMapRes$newMapPos2 = round(newMapRes$newMapPos,1)

for(chrN in 1:20){
  tmp  = newMapRes %>% filter(chr == chrN)
  for(tmpInt in unique(tmp$newMapPos2)){
    tmpIndex = which(tmp$newMapPos2 == tmpInt)
    if(length(tmpIndex) > 0){
      tmp$newMapPos2[tmpIndex] = tmpInt + seq(1,length(tmpIndex),1) * 10^(-6)
    }
  }
  
  newMapRes[row.names(tmp),"newMapPos3"] = tmp$newMapPos2
  
}

table(newMapRes$newMapPos3[2:nrow(newMapRes)] - newMapRes$newMapPos3[1:(nrow(newMapRes)-1)] < 0)

#### 6.2 adjust chr4 and chr19 regions ####
# chr4; shorten old 11.16-16.22 from 5 cM to 1.1 cM; change the positon after 16.26  
chr4Data = newMapRes %>% filter(chr == 4)
index1 = which(chr4Data$newMapPos >= 11.155 & chr4Data$newMapPos <= 16.22)
mapPos = chr4Data$newMapPos[index1]
chr4Data$newMapPos[index1] = NA 
chr4Data$newMapPos[index1[1]] = 11.15502
chr4Data$newMapPos[index1[length(index1)]] = 11.15502 + 1.1
mapPosAdj = as.numeric(round(interpolate(mapPos,chr4Data$newMapPos[index1])[,3],1))

for(tmpInt in unique(mapPosAdj)){
  tmpIndex = which(mapPosAdj == tmpInt)
  if(length(tmpIndex) > 0){
    mapPosAdj[tmpIndex] = tmpInt + seq(1,length(tmpIndex),1) * 10^(-6)
  }
}

chr4Data$newMapPos3[index1] = mapPosAdj
chr4Data$newMapPos3[(index1[length(index1)] + 1) : nrow(chr4Data)] =chr4Data$newMapPos3[(index1[length(index1)] + 1) : nrow(chr4Data)]  - 3.9
newMapRes[chr4Data$Q_ID,"newMapPos3"] = chr4Data$newMapPos3

# chr19; shorten the distance; from 48.98 - 59.95 (11 cM) to 5 cM

chr19Data = newMapRes %>% filter(chr == 19)
index2 = which(chr19Data$newMapPos >= 48.98 & chr19Data$newMapPos <= 60.02)
mapPos = chr19Data$newMapPos[index2]
chr19Data$newMapPos[index2] = NA 
chr19Data$newMapPos[index2[1]] = 48.98566
chr19Data$newMapPos[index2[length(index2)]] = 48.98566 + 5
mapPosAdj = as.numeric(round(interpolate(mapPos,chr19Data$newMapPos[index2])[,3],1))

for(tmpInt in unique(mapPosAdj)){
  tmpIndex = which(mapPosAdj == tmpInt)
  if(length(tmpIndex) > 0){
    mapPosAdj[tmpIndex] = tmpInt + seq(1,length(tmpIndex),1) * 10^(-6)
  }
}

chr19Data$newMapPos3[index2] = mapPosAdj
chr19Data$newMapPos3[(index2[length(index2)] + 1) : nrow(chr19Data)] =chr19Data$newMapPos3[(index2[length(index2)] + 1) : nrow(chr19Data)]  - 6
newMapRes[chr19Data$Q_ID,"newMapPos3"] = chr19Data$newMapPos3

newMapRes$inChrLinkage = NA
newMapRes[pbMarkers,"inChrLinkage"] = "FALSE"

write.csv(newMapRes,"/mnt/soybean/results/markerMapRes_v3.csv",row.names=F)

table(newMapRes$newMapPos3[2:nrow(newMapRes)] - newMapRes$newMapPos3[1:(nrow(newMapRes)-1)] < 0)

ggplot(newMapRes,aes(newMapPos3,oldMapPos))

#### 6.3 adjust the GM_W82_CRxx issue ####
newMapRes = read.csv("/mnt/soybean/results/markerMapRes_v3.csv")
row.names(newMapRes) = newMapRes$Q_ID
GmW82Index = which(str_detect(newMapRes$Q_ID,"Gm_W82"))
GmW82Mk = newMapRes$Q_ID[GmW82Index]
GmW82Base = str_replace(GmW82Mk,"Gm_W82_CR[0-9]{2}","")
table(GmW82Mk %in% phyPos$Q_ID,GmW82Mk %in% markerMap$markerName)
#       FALSE
# TRUE  5088
table(GmW82Base %in% phyPos$Q_ID, GmW82Base %in% markerMap$markerName)
#         FALSE TRUE
# FALSE    90 2320
# TRUE      0 2678

for(i in 1:length(GmW82Mk)){
  mk = GmW82Mk[i]
  mkBase = GmW82Base[i]
  if(mkBase %in% phyPos$Q_ID){ # both have phyPos; just remove the Gm_W82 version; 2678 markers
    newMapRes[mk,"Q_ID"] = NA
  }else if (mkBase %in% markerMap$markerName){ # base does not have phyPos, but in the map; Remove it and change Gm_W82 to base; 2320 markers
    newMapRes[mk,"Q_ID"] = mkBase
    newMapRes[mkBase,"Q_ID"] = NA
    newMapRes[mk,c("oldMapChr","oldMapPos")] = newMapRes[mkBase,c("oldMapChr","oldMapPos")]
  }else{ # base not include in the map; 90 markers
    newMapRes[mk,"Q_ID"] = mkBase
  }
  
}

newMapRes = newMapRes %>% filter(!is.na(Q_ID))
row.names(newMapRes) = newMapRes$Q_ID

# round mapPos to 0.1 & put 1/10^6 into the same bin

newMapRes$newMapPos2 = round(newMapRes$newMapPos3,1)
newMapRes = newMapRes %>% arrange(chr,newMapPos2,Sbegin)
tmpInts = as.numeric()
for(chrN in 1:20){
  tmp  = newMapRes %>% filter(chr == chrN)
  for(tmpInt in unique(tmp$newMapPos2)){
    tmpIndex = which(tmp$newMapPos2 == tmpInt)
    if(length(tmpIndex) > 0){
      tmpInts = c(tmpInts,tmpInt)
      tmp$newMapPos2[tmpIndex] = tmpInt + seq(1,length(tmpIndex),1) * 10^(-6)
    }
  }
  
  newMapRes[row.names(tmp),"newMapPos3"] = tmp$newMapPos2
  
}

table(newMapRes$newMapPos3[2:nrow(newMapRes)] - newMapRes$newMapPos3[1:(nrow(newMapRes)-1)] < 0) # check if in oder

newMapRes$W82_2010_chr = markerMap[newMapRes$Q_ID,"scaffoldName"]
newMapRes$W82_2010_pos = markerMap[newMapRes$Q_ID,"physPos"]

table(newMapRes$chr,newMapRes$W82_2010_chr)

ggplot(newMapRes,aes(newMapPos3,Sbegin/10^6,col=as.character(W82_2010_chr)))  + geom_point() + facet_wrap(vars(chr),scale="free")
ggplot(newMapRes,aes(Sbegin/10^6,newMapPos3,col=as.character(W82_2010_chr)))  + geom_point() + facet_wrap(vars(chr),scale="free")
ggplot(newMapRes,aes(oldMapPos,newMapPos3,col=as.character(W82_2010_chr)))  + geom_point() + facet_wrap(vars(chr),scale="free")

write.csv(newMapRes,"/mnt/soybean/results/soyGenMap_2023_final.csv")

#### 7. check newGenMap vs Kosambi Dist ####
#
# new map

tmpMarkers = colnames(totalRecRate)
res_2$newMapPos = newMapRes[res_2$Q_ID,"newMapPos3"]
newMapRes$newMapPos = newMapRes$newMapPos3
for(chrN in 1:20){
  for(pos in c(10,30,50,70,90)){
    subDf = res_2 %>% filter(chr == chrN, newMapPos >= pos - 1, newMapPos <= pos + 1,!is.na(totalMkLinkedOwnChr))
    if(nrow(subDf) >= 1){
      mk = subDf$Q_ID[which.max(subDf$totalMkLinkedOwnChr)]
      recData2 = plotSingleMarkerRecRateByNewMap(which(tmpMarkers == mk),chrs = chrN)
      png(file=paste0("chr",chrN,"_",pos,"cM.png"),width=1200, height=650)
      print(recData2)
      dev.off()
    }
    
  }
}

tmp = plotSingleMarkerRecRateByNewMap(which(tmpMarkers == "NGMAX006209071"),chrs = 4)
tmp = plotSingleMarkerRecRateByNewMap(which(tmpMarkers == "NGMAX006212546"),chrs = 4)
tmp = plotSingleMarkerRecRateByNewMap(which(tmpMarkers == "NGMAX006214195"),chrs = 4)

# problematic markers (on wrong chr)
tmp = plotSingleMarkerRecRateByNewMap(which(tmpMarkers == "NGMAX009439305"),chrs = c(1,2))
tmp = plotSingleMarkerRecRateByNewMap(which(tmpMarkers == "NGMAX009487630"),chrs = c(18,19))
tmp = plotSingleMarkerRecRateByNewMap(which(tmpMarkers == "NGMAX009441052"),chrs = c(11,2))

#### 8. validate distance with ####
# use the soy GBS to validate 
library(tidyverse)
df = read.table(gzfile("/mnt/data/soy_hackathon_raw_genos_wide_v2.tsv.gz"),sep = '\t',header=T)
row.names(df) = df$pedigree
prtGeno = read.table(gzfile("/mnt/data/soy_hackathon_raw_parent_genos_wide.tsv.gz"),sep = '\t',header=T)
row.names(prtGeno) = prtGeno$pedigree
prtGeno = prtGeno %>% filter(active_marker_count > 100)

origin = "MK0214C2-C0DNN/900Y71"
prts = unlist(strsplit(origin,"/"))
if(all(prts %in% prtGeno$pedigree)){
  subPrtGeno = prtGeno[prts,8:ncol(prtGeno)]
  validMks = apply(subPrtGeno,2,function(x){if(x[1] != x[2] & x[1] %in% c("A/A","C/C","G/G","T/T") & 
                                               x[2] %in% c("A/A","C/C","G/G","T/T")){return(TRUE)}else{return(FALSE)}})
  subPrtGeno = subPrtGeno[,which(validMks)]
  subDf = df[which(df$origin == origin),colnames(subPrtGeno)[colnames(subPrtGeno) %in% colnames(df)]]
  
  validIndex = as.numeric()
  for(j in 1:ncol(subDf)){
    mkName = colnames(subDf)[j]
    pGenos = subPrtGeno[,mkName]
    index1 = which(subDf[,j] == pGenos[1])
    index2 = which(subDf[,j] == pGenos[2])
    subDf[,j] = "Missing"
    subDf[index1,j] = "P1"
    subDf[index2,j] = "P2"
    if(length(index1) > 0 & length(index2) > 0){
      validIndex = c(validIndex,j)
    }
  }
  
  subDf = subDf[,validIndex]
}


len = ncol(subDf)
minPop = 100

for(count in 1:len){
  R = rep(NA,len)
  for(i in 1:len){
    tmp = data.frame(subDf[,i],subDf[,count])
    colnames(tmp) = c("pos1","pos2")
    tmp = tmp %>% filter(pos1 %in% c("P1","P2"), pos2 %in% c("P1","P2"))
    if(nrow(tmp) >= minPop){
      tableData = table(tmp[,1] == tmp[,2])
      R[i] = 1 - tableData["TRUE"] / sum(tableData) 
    }
  }
  
  subRecData = markerMapRes[colnames(subDf),c(1,24:27)]
  subRecData$recRate = R
  ggplot(subRecData,aes(newMapPos,recRate,col=as.character(chr))) + geom_point() + facet_wrap(vars(chr))
  
  R = ifelse(R >= 0.5, 0.49, R)
  r = R/(2-2*R)
  mapDist = as.numeric(round(100*0.25*log((1+2*r)/(1-2*r)),1))
  mapDist = ifelse(mapDist >= 50,50, mapDist) 
  
} 





# use the soy GBS to validate 
library(tidyverse)
df = read.table(gzfile("/mnt/data/soy_hackathon_raw_genos_wide_v2.tsv.gz"),sep = '\t',header=T)
mkInfo = apply(df[,9:53155],2,function(x){y = table(x);if(any(names(y) == "")){return(length(y)-1)}else{return(length(y))}})

popData = df[,c(1:8,8+which(mkInfo > 1))]

majorAF = apply(popData[,9:32900],2,function(x){y = sort(table(x),decreasing = T);naIndex = which(names(y) == "");if(naIndex == 1){return(y[2]/sum(y[2:length(y)]))}
else{return(y[1]/(sum(y)-y[naIndex]))}})
hist(majorAF)
testMkNames = colnames(popData)[9:32900]
testMkNames = testMkNames[testMkNames %in% markerMapRes$Q_ID]

# 
fpGeno = read.delim(gzfile("/domino/datasets/Soybean__MON_2010_v2r1/Soybean__MON_2010_v2r1_extramarkers__withData_LshFpGenos.gz"),header=F)
fpLines = fpGeno$V1
row.names(fpGeno) = fpGeno$V1

lines <- unique(unlist(strsplit(unique(df$origin),"/")))
lines = lines[lines %in% fpLines]
numOfLines <- length(lines)

fpDataVec <- data.frame(matrix("-",ncol=numOfLines,nrow=nrow(markerMap)))
colnames(fpDataVec) <- lines
row.names(fpDataVec) <- markerMap$markerName
count = 0
for(line in lines){
  if(count %% 2000 == 0){print(paste(count,Sys.time()))}
  count = count + 1
  str <- str_to_upper(fpString[line,"V2"])
  fpDataVec[,line] <- unlist(str_split(str,""))
}


# 
allOrigins = unique(popData$origin)
originGenoExist = sapply(allOrigins,function(x){all(unlist(strsplit(x,"/")) %in% lines)})
newOrigins = which(originGenoExist)


for(newOrigin in names(newOrigins)[2]){
  print(paste(newOrigin,Sys.time()))
  subPopData = popData %>% filter(origin %in% newOrigin)
  prts = unlist(strsplit(newOrigin,"/"))
  for(j in 9:ncol(subPopData)){
    mkName = colnames(subPopData)[j]
    if(mkName %in% row.names(fpDataVec)){
      pGenos = fpDataVec[mkName,prts]
      if(all(pGenos %in% c("A","T","C","G")) & pGenos[1] != pGenos[2]){
        pGenos = paste0(pGenos,"/",pGenos)
        index1 = which(subPopData[,j] == pGenos[1])
        index2 = which(subPopData[,j] == pGenos[2])
        subPopData[,j] = "Missing"
        subPopData[index1,j] = "P1"
        subPopData[index2,j] = "P2"
      }else{
        #subPopData[,j] = "Missing"
      }
      
    }
  }
  
}

subPopData = subPopData[,9:ncol(subPopData)]
validMks = apply(subPopData,2,function(x){if(all(c("P1","P2") %in% unique(x))){return(TRUE)}else{return(FALSE)}})
subPopData = subPopData[,which(validMks)]
subPopData = subPopData[,markerMap$markerName[markerMap$markerName %in% colnames(subPopData)]]

af = apply(subPopData,2,function(x){return(length(which(x == "P1")) + length(which(x == "P2")))})

subNewPopData = subPopData[,which(af >= 150)]

allTestMarkers = colnames(subNewPopData)

numCores <- detectCores()
numCores
registerDoParallel(30)

len = length(allTestMarkers)
minPop = 100

foreach(count = 1718:len,.combine=rbind) %dopar% {
  
  # R = apply(chrData,2,function(x){tmp = data.frame(x,chrData[,count]); 
  #     colnames(tmp) = c("pos1","pos2"); tmp = tmp %>% filter(pos1 %in% allPrts, pos2 %in% allPrts); if(nrow(tmp) == 0){return(NA)}else{
  #     tableData = table(tmp[,1] == tmp[,2]); return(1 - tableData["TRUE"] / sum(tableData)) }})
  
  R = rep(NA,len)
  for(i in 1:len){
    tmp = data.frame(subNewPopData[,i],subNewPopData[,count])
    colnames(tmp) = c("pos1","pos2")
    tmp = tmp %>% filter(pos1 %in% c("P1","P2"), pos2 %in% c("P1","P2"))
    if(nrow(tmp) >= minPop){
      tableData = table(tmp[,1] == tmp[,2])
      R[i] = 1 - tableData["TRUE"] / sum(tableData) 
    }
  }
  
  
  
  R = ifelse(R >= 0.5, 0.49, R)
  r = R/(2-2*R)
  mapDist = as.numeric(round(100*0.25*log((1+2*r)/(1-2*r)),1))
  mapDist = ifelse(mapDist >= 50,50, mapDist) 
  #pwDist[count,] = pwDist[,count] = mapDist
  saveRDS(mapDist,paste0("/mnt/tmpData/validPop/mk_",count,".rds"))
}

mkDistMat = matrix(NA,nrow=length(allTestMarkers),ncol=length(allTestMarkers),dimnames = list(allTestMarkers,allTestMarkers))
for(i in 1:len){
  tmp = readRDS(paste0("/mnt/tmpData/validPop/mk_",count,".rds"))
  mkDistMat[i,] = tmp
}
row.names(markerMapRes) = markerMapRes$Q_ID 
subMarkerMapRes  = markerMapRes[allTestMarkers,]
count = 5
test = readRDS(paste0("/mnt/tmpData/validPop/mk_",count,".rds"))
names(test) = allTestMarkers
testMk = allTestMarkers[count]
testMkChr = markerMapRes[testMk,"chr"]
hist(test)
subMarkerMapRes$mkDist = test
subMarkerMapRes[count,]
ggplot(subMarkerMapRes %>% filter(chr == testMkChr),aes(newMapPos,mkDist,col=as.character(chr))) + geom_point() 




tmp = markerMapRes %>% filter(Q_ID %in% allTestMarkers,chr == testMkChr)
tmp$distWithTestMk = test[tmp$Q_ID]
tmp$distInNewMap = abs(tmp$newMapPos - tmp[testMk,"newMapPos"])
plot(tmp$distInNewMap,tmp$distWithTestMk)






#### 9. validation by parental haplotypes ####

soyBinMap = read.csv("/mnt/data/soyNewPredictionMap.csv")

#
newPrtData <- read.csv(gzfile("/domino/datasets/2022_NA_Soybean_Hackathon_ParentHaplotypes_training_version2/results4/allMeta.csv.gz")) 
newPrtData <- newPrtData[,c("pedigree",	"midasGermplasmID",	"originUsed",	"parents",	"includedAncestors")]
row.names(newPrtData) <- newPrtData$pedigree

# take the top 20 origins 
selOrigins = head(sort(table(newPrtData$originUsed),decreasing = T),20)
selOrigins = selOrigins[!str_detect(names(selOrigins),"RM1514A5-C0YNN")] # RM1514A5-C0YNN  does not have LSH data

# pedigree-filename 
pedDict = read.table("/mnt/data/filename_ped_key.txt",header = T)

# pull parental haplotype data

pullParentalHaploData <- function(selOrigins){
  for(origin in names(selOrigins)){
    prts = unlist(strsplit( str_replace_all(origin,"\\[|\\]",""),"/"))
    progLines = subset(newPrtData,originUsed == origin)$pedigree
    tmpPop = as.data.frame(matrix(NA,ncol=nrow(soyBinMap),nrow=length(progLines)))
    row.names(tmpPop) = progLines
    
    for(line in progLines){
      fileName <- paste0("/domino/datasets/2022_NA_Soybean_Hackathon_ParentHaplotypes_training_version2/results4/maxHaploProbs/",pedDict[line,"file_prefix"],".csv")
      maxProbData <- read.csv(fileName, check.names = F)
      tmpPop[line,] = maxProbData$ProcessedParentalHap
    }
    saveRDS(tmpPop,paste0("/mnt/tmpData/popParentalData/",paste(prts,collapse = "___"),".rds"))
    
  }
  
}

#pullParentalHaploData(selOrigins)

# collect all data together
allPopPrtData = data.frame()
for(origin in names(selOrigins)){
  prts = unlist(strsplit( str_replace_all(origin,"\\[|\\]",""),"/"))
  tmpData = readRDS(paste0("/mnt/tmpData/popParentalData/",paste(prts,collapse = "___"),".rds"))
  tmpData$origin = origin
  allPopPrtData = bind_rows(allPopPrtData,tmpData)
}

# pairwise prt LSH data
lshData = read.csv("/mnt/tmpData/popParentalData/pairwise_LSH-2022-12-13.csv")
sharedBins = matrix(0,ncol=nrow(soyBinMap),nrow=nrow(allPopPrtData))  
for(origin in names(selOrigins)){
  prts = unlist(strsplit( str_replace_all(origin,"\\[|\\]",""),"/"))
  tmpLsh = lshData %>% filter(Line1 %in% prts, Line2 %in% prts)
  rowIndex = which(allPopPrtData$origin == origin)
  for(i in 1:nrow(tmpLsh)){
    stIndex = subset(soyBinMap, genMapChr == tmpLsh$Chr[i] & genMapPos == tmpLsh$StartGenPos[i])$index
    edIndex = subset(soyBinMap, genMapChr == tmpLsh$Chr[i] & genMapPos == tmpLsh$EndGenPos[i])$index
    sharedBins[rowIndex,stIndex:edIndex] = 1
  }
}

# calculate pw distance for each chr
chrLength = c(0, cumsum(table(soyBinMap$genMapChr)))
# 0  1487  3132  5535  7526  9167 10792 12117 14139 16125 17870 19606 21437 23776 25644 27446 29397 31364 33258 34933 36457
allPrts = unique(as.character(unlist(sapply(names(selOrigins),function(x){(strsplit(str_replace_all(x,"\\[|\\]",""),"/"))}))))

getPwDistByChr <- function(allPopPrtData,chrLength,allPrts,chr,sharedBins,minPop = 200){
  numCores <- detectCores()
  numCores
  registerDoParallel(numCores - 2)
  
  st = chrLength[chr] + 1
  ed = chrLength[chr+1]
  len = ed - st + 1
  chrData = allPopPrtData[,st:ed]
  print(paste(chr,st,ed,len,sep=";"))
  pwDist = matrix(NA,nrow=len,ncol=len)
  chrSharedBins = sharedBins[,st:ed]
  
  foreach(count = 1:len,.combine=rbind) %dopar% {
    
    # R = apply(chrData,2,function(x){tmp = data.frame(x,chrData[,count]); 
    #     colnames(tmp) = c("pos1","pos2"); tmp = tmp %>% filter(pos1 %in% allPrts, pos2 %in% allPrts); if(nrow(tmp) == 0){return(NA)}else{
    #     tableData = table(tmp[,1] == tmp[,2]); return(1 - tableData["TRUE"] / sum(tableData)) }})
    
    R = rep(NA,len)
    for(i in 1:len){
      tmp = data.frame(chrData[,i],chrData[,count])
      colnames(tmp) = c("pos1","pos2")
      tmpSharedBins = chrSharedBins[,c(i,count)]
      valIndex = which(tmpSharedBins[,1] == 0 & tmpSharedBins[,2] == 0)
      tmp = tmp[valIndex,] %>% filter(pos1 %in% allPrts, pos2 %in% allPrts)
      if(nrow(tmp) >= minPop){
        tableData = table(tmp[,1] == tmp[,2])
        R[i] = 1 - tableData["TRUE"] / sum(tableData) 
      }
    }
    
    
    
    R = ifelse(R >= 0.5, 0.49, R)
    r = R/(2-2*R)
    mapDist = as.numeric(round(100*0.25*log((1+2*r)/(1-2*r)),1))
    mapDist = ifelse(mapDist >= 50,50, mapDist) 
    #pwDist[count,] = pwDist[,count] = mapDist
    saveRDS(mapDist,paste0("/mnt/tmpData/popParentalData/tmpPWDist/pwDist_parentalHaplo_chr",chr,"_",count,".rds"))
  }
  
  # collect all data together
  pwDist = matrix(NA,nrow=len,ncol=len)
  for(count in 1:len){
    tmp = readRDS(paste0("/mnt/tmpData/popParentalData/tmpPWDist/pwDist_parentalHaplo_chr",chr,"_",count,".rds"))
    pwDist[count,] = pwDist[,count] = tmp
  }
  
  saveRDS(pwDist,paste0("/mnt/tmpData/popParentalData/chrPWDist/pwDist_parentalHaplo_chr",chr,".rds"))
}


# for(chr in 1:20){
#   print(paste(chr,Sys.time()))
#   getPwDistByChr(allPopPrtData,chrLength,allPrts,chr,sharedBins)
# }

for(chr in 1:20){
  pwDist = readRDS(paste0("/mnt/tmpData/popParentalData/chrPWDist/pwDist_parentalHaplo_chr",chr,".rds"))
  print(levelplot(pwDist, main=paste0("PwGroupSm - ParentalHaplo: Chr",chr), xlab="", ylab=""))
}



# validate map distances of selected anchors

phyGenKey = 1:20
names(phyGenKey) = c(16,15,8,11,18,3,5,1,6,12,7,17,13,19,20,2,10,4,9,14)

tmpChrRes = data.frame(phyChr=as.numeric(),genChr = as.numeric(),anchor1 = as.numeric(), anchor2 = as.numeric(), mkIndex1=as.numeric(), 
                       mkIndex2=as.numeric(),mkBinPos1=as.numeric(),mkBinPos2=as.numeric(),ancDist=as.numeric(),parHapDist=as.numeric())

for(phyChr in 1:20){
  selPaths = as.numeric(unlist(strsplit(selAnchorPath$selIndex[phyChr],"_")))
  anchorDist = as.numeric(unlist(strsplit(selAnchorPath$distVec[phyChr],"_")))
  gmapData = read.csv(paste0("/mnt/tmpData/tmpAnchorData/subChrData_chr",phyChr,".csv"))
  anchorRes = read.csv(paste0("/mnt/tmpData/tmpAnchorData/anchorRes_chr",phyChr,".csv"))
  genChr = as.numeric(phyGenKey[as.character(phyChr)])
  genPwDist = readRDS(paste0("/mnt/tmpData/popParentalData/chrPWDist/pwDist_parentalHaplo_chr",genChr,".rds"))
  
  for(st in 1:(length(selPaths)-1)){
    mkIndexSt1 = anchorRes$gIndex[selPaths[st]]
    for(ed in (st+1) : length(selPaths)){
      if(sum(anchorDist[st:(ed-1)]) <= 45){
        ancDist = sum(anchorDist[st:(ed-1)])
        mkIndexEd1 = anchorRes$gIndex[selPaths[ed]]
        for(mkSt in mkIndexSt1 : (mkIndexSt1 + 9)){
          for(mkEd in mkIndexEd1 : (mkIndexEd1 + 9)){
            mkBin1 = floor(gmapData$genMapPos[mkSt]*10) + 1
            mkBin2 = floor(gmapData$genMapPos[mkEd]*10) + 1
            
            parHapDist = genPwDist[mkBin1,mkBin2]
            tmpChrRes = bind_rows(tmpChrRes,data.frame(phyChr,genChr,anchor1 = st, anchor2 = ed, mkIndex1 = mkSt, mkIndex2 = mkEd,
                                                       mkBinPos1 = mkBin1, mkBinPos2 = mkBin2, ancDist = ancDist, parHapDist))
          }
          
        }
        
      }
      
    }
    
  }
  
}


#### 10. validation based on Xuehui's POP ####
popData = read.csv("/mnt/data/21GJT_SG2756.csv")
testMarkers = markerMapRes$Q_ID[markerMapRes$Q_ID %in% colnames(popData)]
popData = popData[,testMarkers]


numCores <- detectCores()
numCores
registerDoParallel(30)

len = length(testMarkers)
minPop = 100

foreach(count = 1:len,.combine=rbind) %dopar% {
  
  # R = apply(chrData,2,function(x){tmp = data.frame(x,chrData[,count]); 
  #     colnames(tmp) = c("pos1","pos2"); tmp = tmp %>% filter(pos1 %in% allPrts, pos2 %in% allPrts); if(nrow(tmp) == 0){return(NA)}else{
  #     tableData = table(tmp[,1] == tmp[,2]); return(1 - tableData["TRUE"] / sum(tableData)) }})
  
  R = rep(NA,len)
  for(i in 1:len){
    tmp = data.frame(popData[,i],popData[,count])
    colnames(tmp) = c("pos1","pos2")
    tmp = tmp %>% filter(pos1 %in% c("A","B"), pos2 %in% c("A","B"))
    if(nrow(tmp) >= minPop){
      tableData = table(tmp[,1] == tmp[,2])
      R[i] = 1 - tableData["TRUE"] / sum(tableData) 
    }
  }
  
  
  
  R = ifelse(R >= 0.5, 0.49, R)
  r = R/(2-2*R)
  mapDist = as.numeric(round(100*0.25*log((1+2*r)/(1-2*r)),1))
  mapDist = ifelse(mapDist >= 50,50, mapDist) 
  #pwDist[count,] = pwDist[,count] = mapDist
  saveRDS(mapDist,paste0("/mnt/tmpData/validPop/mk_",count,".rds"))
}

diffRes = data.frame()
mkDistMat = matrix(NA,nrow=length(testMarkers),ncol=length(testMarkers),dimnames = list(testMarkers,testMarkers))
for(count in 1:len){
  if(count %% 1000 == 0){print(paste(count,Sys.time()))}
  if(!file.exists(paste0("/mnt/tmpData/validPop/mk_",count,".rds"))){next;}
  tmp = readRDS(paste0("/mnt/tmpData/validPop/mk_",count,".rds"))
  if(length(tmp) != length(testMarkers)){next;}
  names(tmp) = testMarkers
  mk = testMarkers[count]
  mkNewPos = markerMapRes[mk,"newMapPos"]
  mkOldPos = markerMapRes[mk,"oldMapPos"]
  mkChr = markerMapRes[mk,"chr"]
  subMkData = markerMapRes %>% filter(chr == mkChr, abs(newMapPos - mkNewPos) < 40)
  if(nrow(subMkData) < 1){next;}
  newMapDist = abs(subMkData$newMapPos - mkNewPos) 
  oldMapDist = abs(subMkData$oldMapPos - mkOldPos) 
  diffRes = bind_rows(diffRes,data.frame(makerName = mk,popDist = tmp[subMkData$Q_ID], newMapDist = newMapDist, oldMapDist = oldMapDist ))
  #mkDistMat[i,] = tmp
}

ggplot(diffRes,aes(popDist,newMapDist)) + geom_point() + theme_bw() + xlab("mapDistance_validationPop") + ylab("mapDistance (2023)")
ggplot(diffRes,aes(popDist,oldMapDist)) + geom_point() + theme_bw() + xlab("mapDistance_validationPop") + ylab("mapDistance (2010)")

ggplot(diffRes,aes(popDist - newMapDist)) + geom_histogram() + theme_bw() + ggtitle("DiffDistance (cM) - 2023 map vs validationPop")
ggplot(diffRes,aes(popDist - oldMapDist)) + geom_histogram() + theme_bw() + ggtitle("DiffDistance (cM) - 2010 map vs validationPop")

ggplot(diffRes,aes(popDist - newMapDist)) + geom_density() + theme_bw() + ggtitle("DiffDistance (cM) - 2023 map vs validationPop")
ggplot(diffRes,aes(popDist - oldMapDist)) + geom_density() + theme_bw() + ggtitle("DiffDistance (cM) - 2010 map vs validationPop")

ggplot(diffRes,aes(abs(popDist - oldMapDist)/abs(popDist - newMapDist) )) + geom_density() + theme_bw() + ggtitle("DiffDistance (cM) - 2023 map vs validationPop") +
  xlim(-5,5)

ggplot(diffRes,aes(abs(popDist - oldMapDist) - abs(popDist - newMapDist) )) + geom_density() + theme_bw() + ggtitle("DiffDistance (cM) - 2023 map vs validationPop")

ggplot(diffRes,aes((popDist - newMapDist), (popDist - oldMapDist))) + geom_point() + theme_bw() + ggtitle("DiffDistance (cM) - 2023 map vs validationPop")

ggplot(diffRes,aes((popDist - oldMapDist), (popDist - newMapDist))) + geom_point() + theme_bw() + ggtitle("DiffDistance (cM) - 2023 map vs validationPop")


####


tmpChrRes$genMapDist = (tmpChrRes$mkBinPos2 - tmpChrRes$mkBinPos1) / 10 
tmpChrRes$neighborAnc = ifelse(tmpChrRes$anchor2 - tmpChrRes$anchor1 == 1, "Yes","No")

saveRDS(tmpChrRes,"/mnt/results/20221024/parHapDist_anchorDist.rds")


plot(tmpChrRes$ancDist,tmpChrRes$parHapDist)
cor(tmpChrRes$ancDist,tmpChrRes$parHapDist,use="na.or.complete")
cor(tmpChrRes$genMapDist,tmpChrRes$parHapDist,use="na.or.complete")

ggplot(data=tmpChrRes,aes(parHapDist,ancDist,col=as.character(phyChr))) + geom_point() + facet_wrap(vars(phyChr))+ theme_bw()
ggplot(data=tmpChrRes,aes(parHapDist - ancDist)) + geom_histogram() + theme_bw()
ggplot(data=tmpChrRes,aes(parHapDist - genMapDist)) + geom_histogram() + theme_bw()
ggplot(data=tmpChrRes,aes(ancDist - genMapDist)) + geom_histogram() + theme_bw()

ggplot(data=tmpChrRes %>% filter(genMapDist > 0),aes(parHapDist,genMapDist,col=as.character(phyChr))) + geom_point() + facet_wrap(vars(phyChr)) + theme_bw()


testRes = tmpChrRes %>% group_by(phyChr,anchor1, anchor2) %>% summarise(ancDist = mean(ancDist),avgParHapDist = mean(parHapDist))
testRes$neighborAnc = ifelse(testRes$anchor2 - testRes$anchor1 == 1, "Yes","No")
cor(testRes$avgParHapDist,testRes$ancDist,use="na.or.complete")

ggplot(data=testRes,aes(avgParHapDist,ancDist,col=as.character(phyChr))) + geom_point() + facet_wrap(vars(phyChr))+ theme_bw()
ggplot(data=testRes %>% filter(neighborAnc == "Yes"),aes(avgParHapDist,ancDist,col=as.character(phyChr))) + 
  geom_point() + facet_wrap(vars(phyChr),scale="free")+ theme_bw() + ggtitle("Neigthbor anchor distances")

testRes %>% filter(neighborAnc == "Yes") %>%  group_by(phyChr) %>% summarise(totalAncLen = sum(ancDist), totalParHapLen = sum(avgParHapDist,na.rm=T))
testRes2 = testRes %>% filter(neighborAnc == "Yes")
cor(testRes2$avgParHapDist,testRes2$ancDist,use="na.or.complete")

#

## lineResRec data
convertBack2PBH <- function(cndStr){
  # cndStr is a condense string from condensStr
  # eg. A::4;T::3;C::1;A::2;T::2 => AAAATTTCAATT
  tmp <- unlist(lapply(unlist(strsplit(cndStr,";")),
                       function(x){y <- unlist(strsplit(x,"::")); z <- as.numeric(y[2]);
                       names(z) <- y[1];return(z)}))
  tmp1 <- data.frame(name=names(tmp),repTimes = as.numeric(tmp),stringsAsFactors = F)
  str1 <- as.character(unlist(apply(tmp1,1,function(x){return(rep(x[1],x[2]))})))
  return(str1)
}

cbh = readRDS("/mnt/data/strCondensed_AncestralComboHaplotypes_20220825.rds")


lineRecRes = readRDS("/mnt/results/20221024/markerRecByLine.rds")
mk1 = "NGMAX005840819"
mk2 = "NGMAX009505606"
totalRecRate[mk1,mk2]
totalTrioCount[mk1,mk2]
tmpRecRes = lineRecRes %>% filter(markers %in% c(mk1,mk2))
test = tmpRecRes %>% select(markers,parent,line) %>% spread(markers,parent)
colnames(test)[2:3] = c("mk1","mk2")
test = test %>% filter(!is.na(mk1),!is.na(mk2))

testLines = c(test$line,"A3555")
testLines = testLines[testLines %in% cbh$lineName] 

subPbh <- sapply(cbh[testLines,"cbh"],convertBack2PBH)
colnames(subPbh) <- testLines

markerMap[c(mk1,mk2),]
soyBinMap %>% filter(genMapChr == 8, genMapPos == 3.7)
soyBinMap %>% filter(genMapChr == 7, genMapPos == 120.3)
head(sort(table(subPbh[12155,]),decreasing = T))
head(sort(table(subPbh[11996,]),decreasing = T))