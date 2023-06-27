#### functions ####
# These are functions used for trio mapping

## load all required packages
lapply(c('tidyverse','parallel','MASS','foreach','doParallel','lattice'), require, character.only = TRUE)

### 1. get genotype data (table)

#'  This function is to convert fpString to a table format,with row as the number of markers and column as the number of lines;
#'  takes about 3 mins to generate the fpMatrix for soybean: 67717 markers x 15634 lines
#'  
#' @param pedInfo data.frame the file with pedigree information; containing at least three columns: "pedigree", "prt1", "prt2"
#' @param fpGenoe data.frame a table format, with row as line, and column V2 as the fpString for each line. Cloolected by Danny and team
#' @param markerMap data.frame a table format; The genetic map for all the markers used in fpGeno
#' 
#' @return data.frame a table format, with row as the number of markers and column as the number of lines
convertFpStringToGeno <- function(pedInfo,fpGeno,markerMap){
  lines <- unique(c(pedInfo$pedigree,pedInfo$prt1,pedInfo$prt2))
  fpDataVec <- data.frame(matrix("-",ncol=length(lines),nrow=nrow(markerMap),dimnames = list(markerMap$markerName, lines)))
  for(line in lines){
    if(line %in% row.names(fpGeno)){
      str <- str_to_upper(fpGeno[line,"V2"])
      fpDataVec[,line] <- unlist(str_split(str,""))
    }else{
      warning(paste0("Warning!!! ",line," is not included in the fpGeno dataset"))
    }
  }
  return(fpDataVec)
}

### 2. recombination rate
#### 2..1 informative markers by trio; Including matching status, P1 or P2;

#'  This function is to find informative markers for each trio and then return its matching status with parents: P1 or P2
#'
#' @param pedInfo data.frame the file with pedigree information; containing at least three columns: "pedigree", "prt1", "prt2"
#' @param fpDataVec data.frame a table format, with row as the number of markers and column as the number of lines
#' @param nCores integer the number of cores to be used
#' @param tmpFolder string a tmp folder to save the output
#' @param genos vector<string> homozygous genotypes that the fpDataVec would include
#' 
#' @return null output is saved to tmpFolder in a csv format with one marker per row and following columns, O, P1, P2, markers, parent
infoMarkerByLine <- function(pedInfo,fpDataVec,nCores,tmpFolder,genos = c("AA","TT","CC","GG")){
  if(!dir.exists(tmpFolder)){dir.create(tmpFolder)}
  numCores <- detectCores()
  if(nCores < numCores){
    registerDoParallel(nCores)
    print(paste("Using",nCores,"Cores out of the computer's maximum of", numCores,"cores"))
  }else{
    nCores = numCores
    registerDoParallel(numCores)
    print(paste("Using",nCores,"Cores out of the computer's maximum of", numCores,"cores"))
  }
  
  linesPerRun <- ceiling(nrow(pedInfo) / nCores)
  
  foreach(batchIndex = 1:nCores,.combine=rbind) %dopar% {
    st <- (batchIndex-1) * linesPerRun + 1
    ed <- min( batchIndex*linesPerRun,nrow(pedInfo))
    
    subPedInfo <- pedInfo[st:ed,]
    
    for(k in 1:nrow(subPedInfo)){
      qLine <- subPedInfo$pedigree[k]
      prts <- c(subPedInfo$prt1[k],subPedInfo$prt2[k])
      if(all(c(qLine,prts) %in% colnames(fpDataVec))){
        outData <- fpDataVec[,as.character(c(qLine,prts))]
        colnames(outData) <- c("O","P1","P2")
        
        # 1) skip if any missing data or heterozygous presents
        outData <- outData %>% 
                        filter(O %in% genos,P1 %in% genos, P2 %in% genos,P1 != P2)
        outData$markers = row.names(outData)
        outData$parent = 0 # not matching any parent
        outData[outData$O == outData$P1,"parent"] = 1
        outData[outData$O == outData$P2,"parent"] = 2
        outData = outData %>% filter(parent %in% c(1,2))
        
        # 2) skip if no data or only one parental genotype
        if(nrow(outData) <= 1){next;}
        if(length(unique(outData$parent)) == 1){next;}
        
        write.csv(outData,paste0(tmpFolder,"/",qLine,".csv"),row.names=F)
      }
    }
  }
}


#### 2.2 recombination data for all trios; #total Trios and #recombined Trios

#'  This function is to count #totalTros and #recombinedTrios for each pair of markers
#'
#' @param tmpFolder string a tmp folder from each trio; Informative markers and Parental Status;
#' @param allMarkers vector<string> a list of all markers from the genetic map
#' @param outFileFolder, string a folder to store all the output data
#' @param verbose string whether or not to print information
#' 
#' @return null outputs are saved to outFileFolder;  Output are two matrix, each one is a large matrix with the dimension as #markers x #markers (around 65k x 65k)

generateRecData = function(tmpFolder,allMarkers,outFileFolder,verbose = TRUE){
  # This function is to generate recbinate rate for any pair of markers
  # For every trio, collect the informative markers and found which pair of markers are recombined (saved at tmpFolder)
  # collect all infoMarkerByLine results, and then count recombined lines and non-recombined lines;
  
  if(!dir.exists(outFileFolder)){dir.create(outFileFolder)}
  
  totalTrioCount <- data.frame(matrix(0,nrow=length(allMarkers),ncol=length(allMarkers)))
  row.names(totalTrioCount) <- allMarkers
  colnames(totalTrioCount) <- allMarkers
  
  recTrioCount <- data.frame(matrix(0,nrow=length(allMarkers),ncol=length(allMarkers)))
  row.names(recTrioCount) <- allMarkers
  colnames(recTrioCount) <- allMarkers
  
  allFiles = list.files(tmpFolder)
  testedMarkers= as.character()
  count = 0
  
  for(fileName in allFiles){
    count = count + 1
    if((count-1) %% 1000 == 0 & verbose){print(paste("Processed",count,"lines;",Sys.time()))}
    tmpData = read.csv(paste0(tmpFolder,"/",fileName))
    qLine = unlist(strsplit(fileName,"\\."))[1]

    testedMarkers = unique(c(testedMarkers,tmpData$markers))
    
    tmpDataMarkers = tmpData$markers
    
    totalTrioCount[tmpDataMarkers,tmpDataMarkers] <- totalTrioCount[tmpDataMarkers,tmpDataMarkers] + 1 # sum recombined and none-recombined trios
    
    p1Index <- which(tmpData$parent == 1)
    p2Index <- which(tmpData$parent == 2)
    
    # recombination happened between markers having different parent
    recTrioCount[tmpDataMarkers[p1Index],tmpDataMarkers[p2Index]] <- recTrioCount[tmpDataMarkers[p1Index],tmpDataMarkers[p2Index]] + 1 # only recombined trios
    recTrioCount[tmpDataMarkers[p2Index],tmpDataMarkers[p1Index]] <- recTrioCount[tmpDataMarkers[p2Index],tmpDataMarkers[p1Index]] + 1 # only recombined trios
  }
  
  recTrioCount = recTrioCount[testedMarkers,testedMarkers]
  totalTrioCount = totalTrioCount[testedMarkers,testedMarkers]
  
  cat(paste("... Saving Data ....;", "\n",Sys.time()))
  saveRDS(totalTrioCount,paste0(outFileFolder,"/totalTrioCount_allChr.rds"))
  saveRDS(recTrioCount,paste0(outFileFolder,"/recTrioCount_allChr.rds"))
  testRecRate <- round(recTrioCount/totalTrioCount,4) # recombination rate for each markers (#recTrios / #totalTrios)
  saveRDS(testRecRate,paste0(outFileFolder,"/recRate_allChr.rds"))
  cat(paste("... Rec Data Done!!! ....;", "\n",Sys.time()))
}


### 3. define groups and determine anchors

#'  This function is to find closely linked markers based on maxRecRate; This could tell if a marker is off chr/pos (bwteen genetic and phyical pos)
#'
#' @param phyPos data.frame a table format that lists the physical positon of each marker; Includes at least three columns Q_ID/markerName, chr, Sbegin/pos
#' @param testMarkers vector<string> a list of markers
#' @param totalTrioCount data.frame a table format that includes all the trios for each pair of marker (row and column)
#' @param totalRecRate data.frame a table format that includes all the recombination rate for each pair of marker (row and column)
#' @param minCount integer the minimum number of trio required
#' @param maxRecRate numeric the maximum recombination rate required
#' 
#' @return data.frame to includes all the information for each marker. For each marker, find other markers that linked with this marker and summarize the linkage information;

getLinkedMarkerInfo <- function(phyPos,testMarkers, totalTrioCount, totalRecRate,minCount, maxRecRate){
  
  res <- phyPos %>% filter(Q_ID %in% testMarkers) %>% dplyr::select(Q_ID, chr, Sbegin)
  row.names(res) = res$Q_ID
  res = res[testMarkers,]
  res$chrDist <- NA
  res$numPerChr <- NA
  res$totalMkLinkedOwnChr <- NA
  res$totalMkLinked <- NA
  res$LDMarkers <- NA
  
  for(rowIndex in 1:nrow(res)){
    tmp <- data.frame(chr = res$chr,totalCount = totalTrioCount[testMarkers,res$Q_ID[rowIndex]],
                      recRate = totalRecRate[testMarkers,res$Q_ID[rowIndex]]) 
    tmp$marker = testMarkers
    tmp <- tmp %>% filter(totalCount >= minCount, recRate <=  maxRecRate)
    tableTmp <- table(tmp$chr) 
    chrDist <- paste(names(tableTmp),collapse =";")
    chrNum <- paste(tableTmp,collapse =";")
    res$chrDist[rowIndex] <- chrDist
    res$numPerChr[rowIndex] <- chrNum
    res$totalMkLinkedOwnChr[rowIndex] = tableTmp[as.character(res$chr[rowIndex])]
    res$totalMkLinked[rowIndex] = sum(tableTmp)
    res$LDMarkers[rowIndex] = paste(tmp$marker,collapse =";")
    
  }
  
  res$numLinkedChr = sapply(res$chrDist,function(x){return(length(unlist(strsplit(x,";"))))})
  res$index = 1:nrow(res)
  
  return(res)
}


#'  This function is to find anchor to each chr; All tmporary data are saved at tmpOutFolder
#'
#' @param testChr integer the query chr
#' @param mkLinkageInfo data.frame a table format that includes marker linkage information
#' @param totalTrioCount data.frame a table format that includes all the trios for each pair of marker (row and column)
#' @param totalRecRate data.frame a table format that includes all the recombination rate for each pair of marker (row and column)
#' @param tmpOutFolder string a tmp folder to save the output
#' @param minCount integer the minimum number of trio required
#' @param minLinkedMk integer the minimum number of linked markers
#' @param mkPerGroup integer the required number of markers per group
#' @param minDataPoint integer the minimum number of data point per group
#' @param minDataPointPerMk integer the minimum number of data point per marker
#' @param maxDist numeric the maximum distance
#' @param numOfDataBwGroups integer the minimum number of data between two groups

#' 
#' @return null

findAnchorsByChr = function(testChr, mkLinkageInfo, totalTrioCount, totalRecRate,tmpOutFolder,minCount = 200,minLinkedMk = 5,
                            mkPerGroup = 10,minDataPoint = 80, minDataPointPerMk = 5, maxDist = 0.5,numOfDataBwGroups = 20) {
  if(!dir.exists(tmpOutFolder)){dir.create(tmpOutFolder)}
  testChrData = rmkLinkageInfo %>% filter(chr == testChr, chrDist == testChr, totalMkLinkedOwnChr >= minLinkedMk) %>% arrange(Sbegin)
  testChrMk = testChrData$Q_ID
  chrRecCount = totalTrioCount[testChrMk,testChrMk]
  chrRecRate = totalRecRate[testChrMk,testChrMk]
  chrRecCount = ifelse(chrRecCount >= minCount,1,NA)
  chrRecRate = as.matrix(chrRecRate * chrRecCount)
  chrRecRate = ifelse(chrRecRate >= 0.5,0.49,chrRecRate)
  chrRecRate = chrRecRate / (2 - 2 * chrRecRate)
  chrMapDist = 0.25 * log((1+2*chrRecRate) / (1 - 2*chrRecRate)) * 100 #kosambi distance; need to set
  chrMapDist = ifelse(chrMapDist >= 50,50,round(chrMapDist,2)) #kosambi distance; need to set 
  
  
  plotMapDist = chrMapDist
  row.names(plotMapDist) = NULL
  colnames(plotMapDist) = NULL
  
  # pdf(file=paste0("/mnt/results/20221024/heatmaps/pwMarkerDist_Chr",testChr,".pdf"),width=12,height=12)
  # print(levelplot(plotMapDist, main=paste0("pwMarkerDist: Chr",testChr), xlab="", ylab=""))
  # dev.off()
  
  
  rm(plotMapDist)
  
  offPosMkData = findOffPosMk(chrMapDist)
  offPosMks = subset(offPosMkData,pbMks=="T")$mkName
  
  # summarize within group map distance; Ignore markers that are questionable
  
  subChrData = testChrData %>% filter(!Q_ID %in% offPosMks) %>% arrange(Sbegin)
  write.csv(subChrData,paste0(tmpOutFolder,"/subChrData_chr",testChr,".csv"))
  
  subChrMapDist = chrMapDist[subChrData$Q_ID,subChrData$Q_ID]
  saveRDS(subChrMapDist,paste0(tmpOutFolder,"/subChrMapDist_chr",testChr,".rds"))
  
  numOfGroups = nrow(subChrData) - (mkPerGroup - 1)
  groupSm = data.frame(gIndex = 1:numOfGroups,numOfData = rep(NA,numOfGroups),minDataPerMk = rep(NA,numOfGroups),minD=rep(NA,numOfGroups),
                       maxD = rep(NA,numOfGroups), meanD = rep(NA,numOfGroups), medianD = rep(NA,numOfGroups),stdD = rep(NA,numOfGroups))
  for(i in 1:numOfGroups){
    st = i
    ed = i + (mkPerGroup - 1)
    posIndex = st:ed
    subMapDist = subChrMapDist[posIndex,posIndex]
    groupSm[i,] = c(i,length(which(subMapDist >= 0)),min(apply(subMapDist,1,function(x){length(which(x >=0))})),min(subMapDist,na.rm = T),
                    max(subMapDist,na.rm = T),mean(as.matrix(subMapDist),na.rm = T),median(as.matrix(subMapDist),na.rm = T),sd(as.matrix(subMapDist),na.rm = T))
    
  }
  
  write.csv(groupSm,paste0(tmpOutFolder,"/groupSm_chr",testChr,".csv"))
  
  pwGroupSm = matrix(NA,nrow=numOfGroups,ncol=numOfGroups)
  
  for(i in 1:numOfGroups){
    st1 = i
    ed1 = i + (mkPerGroup - 1)
    for(j in i:numOfGroups){
      st2 = j
      ed2 = j + (mkPerGroup - 1)
      if(i == j){
        pwGroupSm[i,j] = 0
        pwGroupSm[j,i] = 0
      }else{
        subMapDist = subChrMapDist[st1:ed1,st2:ed2]
        if(length(which(subMapDist >= 0)) >= numOfDataBwGroups){
          pwGroupSm[i,j] = mean(subMapDist,na.rm=T)
          pwGroupSm[j,i] = mean(subMapDist,na.rm=T)
        }
      }
    }
    
  }
  
  #print(levelplot(pwGroupSm, main=paste0("pwGroupSm: Chr",testChr, " (>= 20% data)"), xlab="", ylab=""))
  saveRDS(pwGroupSm,paste0(tmpOutFolder,"/pwGroupSm_chr",testChr,".rds"))
  
  # pick groups that are most reliable
  selGroupSm = groupSm %>% filter(numOfData >= minDataPoint, minDataPerMk >= minDataPointPerMk, maxD <= maxDist)
  subPwGroupSm = pwGroupSm[selGroupSm$gIndex,selGroupSm$gIndex]
  
  # pdf(file=paste0("/mnt/results/20221024/heatmaps/PwGroupSm_Chr",testChr,".pdf"),width=12,height=12)
  # print(levelplot(subPwGroupSm, main=paste0("PwGroupSm - anchor: Chr",testChr, " (10 mk per group)"), xlab="", ylab=""))
  # dev.off()
  
  print(dim(subPwGroupSm))
  
  #1. from every position, go both left and right. Sample many times, get summary statistics
  #2. For every sampled data, correlate newly calculated distance versus the pair-wise data
  
  anchorRes = data.frame(stIndex = 1:nrow(subPwGroupSm),gIndex = selGroupSm$gIndex,totalGenLen = rep(NA,nrow(subPwGroupSm)),MAE = rep(NA,nrow(subPwGroupSm)),
                         MAE_pct = rep(NA,nrow(subPwGroupSm)),distVec = rep(NA,nrow(subPwGroupSm)), selIndex = rep(NA,nrow(subPwGroupSm)),
                         numOfGroups = rep(NA,nrow(subPwGroupSm)),maxAnchorAE = rep(NA,nrow(subPwGroupSm)),maxAnchorAEGroups = rep(NA,nrow(subPwGroupSm)))
  
  for(i in 1:nrow(subPwGroupSm)){
    anchors = findLinkedGroups(i,pwDist = subPwGroupSm)
    if(length(anchors$selIndex) <= 1){next}
    tmpDiff = getGroupAnchorDiff(anchors,pwDist = subPwGroupSm)
    tmpDiff$sumAE = apply(tmpDiff[,c(1,3,5)],1,function(x){sum(x,na.rm=T)})
    tmpDiff$sumCount = apply(tmpDiff[,c(2,4,6)],1,function(x){sum(x,na.rm=T)})
    mae = sum(tmpDiff[,c(1,3,5)],na.rm = T) / sum(tmpDiff[,c(2,4,6)], na.rm=T)
    maePct = sum(tmpDiff$sumAE / tmpDiff$sumCount / anchors$distVec * tmpDiff$sumCount,na.rm = T) / sum(tmpDiff$sumCount,na.rm = T )
    anchorRes[i,c(3:5,8)] = c(sum(anchors$distVec,na.rm=T),mae,maePct,length(anchors$selIndex))
    anchorRes[i,6:7] = c(paste(round(anchors$distVec,2),collapse = "_"),paste(anchors$selIndex,collapse = "_"))
    aeByAnchor = apply(tmpDiff[,c(1,3,5)],1,function(x){return(sum(x,na.rm=T))})
    numByAnchor = apply(tmpDiff[,c(2,4,6)],1,function(x){return(sum(x,na.rm=T))})
    anchorRes[i,9:10] = c(max(aeByAnchor/numByAnchor,na.rm=T),numByAnchor[which.max(aeByAnchor/numByAnchor)])
  }
  
  write.csv(anchorRes,paste0(tmpOutFolder,"/anchorRes_chr",testChr,".csv"))
}


#'  This function is to find informative markers for each trio and then return its matching status with parents: P1 or P2
#'
#' @param chrMapDist data.frame pairwise distance
#' @param minDist numeric the minimum distance 
#' @param maxDist numeric the maximum distance
#' 
#' @return data.frame a table format; 

getNeigboringGroups = function(chrMapDist,minDist, maxDist){
  
  tmpRes = data.frame(mk = 1:ncol(chrMapDist),lMks=NA,rMks=NA)
 
  for(i in 1:ncol(chrMapDist)){
    tmpDist = chrMapDist[,i]
    mkIndex = which(tmpDist >= minDist & tmpDist <= maxDist)
    rightIndex = mkIndex[which(mkIndex > i)]
    leftIndex = mkIndex[which(mkIndex < i)]
    if(length(leftIndex) > 0){
      tmpRes$lMks[i] = paste(leftIndex,collapse = "_")
    }
    
    if(length(rightIndex) > 0){
      tmpRes$rMks[i] = paste(rightIndex,collapse = "_")
    }
    
  }
  return(tmpRes)
}


#'  This function is to calcuate Pairwise marker MAE; using the markers between any given pairs to calculate MAE
#'
#' @param chrMapDist data.frame pairwise distance
#' @return data.frame a table format; 

getPwMkMAE = function(chrMapDist){
  tmpRes = matrix(NA,nrow=nrow(chrMapDist),ncol=ncol(chrMapDist),dimnames = list(colnames(chrMapDist),colnames(chrMapDist)))
  for(i in 1:(ncol(chrMapDist)-2)){
    for(j in (i+2) : ncol(chrMapDist)){
      if(is.na(chrMapDist[i,j])){next;}
      tmp = chrMapDist[i:j,c(i,j)]
      rowSum = rowSums(tmp)
      rowSum = rowSum[which(!is.na(rowSum))]
      if(length(rowSum) > 2) {# having at least 1 extra value
        tmpRes[i,j] = tmpRes[j,i] = round(mean(abs(rowSum - rowSum[1])[2:(length(rowSum)-1)]),2)   
      }
    }
  }
  return(tmpRes)
}


#'  This function is to find closely linked anchors from both left and right; Anchor path
#'
#' @param initialIndex integer the query index
#' @param pwDist data.frame the pairsie group similary data
#' @param minMapDist integer the minimum distance between two anchors/groups
#' @param minMapDistEdge integer the minimum distance between two anchors/groups at the end of a chromosome
#' 
#' @return list

findLinkedGroups = function(initialIndex,pwDist = subPwGroupSm,minMapDist = 5,minMapDistEdge = 1){
  if(initialIndex == 1){
    rightSel = findLinkedGroupFromRight(initialIndex, pwDist,minMapDist,minMapDistEdge)
    return(list(distVec=rightSel$distVec, selIndex = c(initialIndex,rightSel$selIndex)))
  }else if(initialIndex == nrow(pwDist)){
    leftSel = findLinkedGroupFromLeft(initialIndex, pwDist,minMapDist,minMapDistEdge)
    return(list(distVec=leftSel$distVec, selIndex = leftSel$selIndex))
  }else{
    rightSel = findLinkedGroupFromRight(initialIndex, pwDist,minMapDist,minMapDistEdge)
    leftSel = findLinkedGroupFromLeft(initialIndex, pwDist,minMapDist,minMapDistEdge)
    return(list(distVec=c(leftSel$distVec,rightSel$distVec), selIndex = c(leftSel$selIndex,rightSel$selIndex)))
  }
}

#'  This function is to find closely linked anchors from right; Anchor path
#'
#' @param initialIndex integer the query index
#' @param pwDist data.frame the pairsie group similary data
#' @param minMapDist integer the minimum distance between two anchors/groups
#' @param minMapDistEdge integer the minimum distance between two anchors/groups at the end of a chromosome
#' 
#' @return list

findLinkedGroupFromRight = function(initialIndex,pwDist = subPwGroupSm,minMapDist,minMapDistEdge){
  # 
  if(initialIndex >= nrow(pwDist) | initialIndex <= 0){
    print("Error!! Wrong initial index")
    stop()
  }
  distVec = as.numeric()
  selGroupIndex = as.numeric()
  stIndex = initialIndex
  
  while(stIndex < nrow(pwDist)){
    qIndex = stIndex:nrow(pwDist)
    sIndex = which(pwDist[stIndex,qIndex] >= minMapDist)
    if(length(sIndex) < 1){
      if(max(pwDist[stIndex, qIndex],na.rm=T) >= minMapDistEdge){# require at least 1 cM length
        distVec = c(distVec,max(pwDist[stIndex, qIndex],na.rm=T))
        selGroupIndex = c(selGroupIndex,qIndex[which.max(pwDist[stIndex, qIndex])])
      } 
      stIndex = nrow(pwDist)
    }else{
      tmpIndex = min(sIndex)
      tmpIndex = qIndex[tmpIndex]
      selGroupIndex = c(selGroupIndex,tmpIndex)
      distVec = c(distVec,pwDist[stIndex, tmpIndex])
      stIndex = tmpIndex
    }
  }
  return(list(distVec = distVec,selIndex = selGroupIndex))
  
}

#'  This function is to find closely linked anchors from left; Anchor path
#'
#' @param initialIndex integer the query index
#' @param pwDist data.frame the pairsie group similary data
#' @param minMapDist integer the minimum distance between two anchors/groups
#' @param minMapDistEdge integer the minimum distance between two anchors/groups at the end of a chromosome
#' 
#' @return list
findLinkedGroupFromLeft = function(initialIndex,pwDist = subPwGroupSm,minMapDist,minMapDistEdge){
  # 
  if(initialIndex > nrow(pwDist) | initialIndex <= 1){
    print("Error!! Wrong initial index")
    stop()
  }
  distVec = as.numeric()
  selGroupIndex = stIndex = initialIndex
  while(stIndex > 1 ){
    qIndex = 1: stIndex
    sIndex = which(pwDist[stIndex,qIndex] >= minMapDist)
    if(length(sIndex) < 1){
      if(max(pwDist[stIndex, qIndex],na.rm = T) >= minMapDistEdge){ # request at least 1cM distance
        distVec = c(distVec,max(pwDist[stIndex, qIndex],na.rm = T))
        selGroupIndex = c(selGroupIndex,qIndex[which.max(pwDist[stIndex, qIndex])])
      }
      
      stIndex = 1
    }else{
      tmpIndex = max(sIndex)
      tmpIndex = qIndex[tmpIndex]
      selGroupIndex = c(selGroupIndex,tmpIndex)
      distVec = c(distVec,pwDist[stIndex, tmpIndex])
      stIndex = tmpIndex
    }
  }
  return(list(distVec = rev(distVec),selIndex = rev(selGroupIndex)))
  
}

#'  This function is to evaluate anchors based on MAE; For any pair of anchor, MAE is calculated from groups from their left (after last anchor), middle and right (before last anchor)

#'
#' @param anchors list a output from the function of "findLinkedGroups", with two components: distVec and selIndex; 
#' @param pwDist data.frame the pairsie group similary data
#' 
#' @return data.frame

getGroupAnchorDiff <- function(anchors, pwDist = subPwGroupSm){
  distVec = anchors$distVec
  selIndex = anchors$selIndex
  
  res = data.frame(totalLeft = rep(NA,length(selIndex)-1), numLeft = rep(NA,length(selIndex)-1),
                   totalMid = rep(NA,length(selIndex)-1), numMid = rep(NA,length(selIndex)-1),
                   totalRight = rep(NA,length(selIndex)-1), numRight = rep(NA,length(selIndex)-1))
  
  for(i in 1:(length(selIndex) - 1)){
    a1 = selIndex[i]
    a2 = selIndex[i+1]
    leftIndex ='if'(i > 1, (selIndex[i-1]+1):(a1-1),1:(a1-1))
    rightIndex = 'if'(i < length(selIndex) - 1, (a2+1):(selIndex[i+2]-1), (a2+1):ncol(pwDist))
    midIndex = (a1 + 1) : (a2 - 1)
    
    # diff from left
    if(!(length(leftIndex) == 2 & leftIndex[2] < leftIndex[1])){
      tmp = abs(pwDist[a2,leftIndex] - pwDist[a1,leftIndex] - rep(distVec[i],length(leftIndex)))
      res[i,1:2] = c(sum(tmp,na.rm = T),length(which(!is.na(tmp))))
    }
    
    # diff from middle 
    if(!(length(midIndex) == 2 & midIndex[2] < midIndex[1])){
      tmp = abs(pwDist[a2,midIndex] + pwDist[a1,midIndex] - rep(distVec[i],length(midIndex)))
      res[i,3:4] = c(sum(tmp,na.rm = T),length(which(!is.na(tmp))))
    }
    
    # diff from right
    if(!(length(rightIndex) == 2 & rightIndex[2] < rightIndex[1])){
      tmp = abs(pwDist[a1,rightIndex] - pwDist[a2,rightIndex] - rep(distVec[i],length(rightIndex)))
      res[i,5:6] =  c(sum(tmp,na.rm = T),length(which(!is.na(tmp))))
    }
    
  }
  
  return(res)
  
}




### 4. anchor distances and marker positions

#'  This function is to coonect all anchors from a given anchor
#'
#' @param goodAnchors data.frame selected anchors that are in good quality in terms of data point, maxDistance, etc
#' @param selAnchorPath string the path for a select anchor 
#' @param anchorRes data.frame a table format including anchor information starting from a given anchor
#' @param pwGroupSm data.frame the pairsie group similary data
#' 
#' @return vector<string> the distance between groups in the given anchor path

goodAnchorDist = function(goodAnchors,selAnchorPath,anchorRes,pwGroupSm){
  # Starting from the selected anchor path, map all the other anchors to the map
  stDist = 0
  allGroupDist = as.numeric()
  lastDist = 0
  
  for(i in 0:length(selAnchorPath)){
    if(i == 0){
      if(selAnchorPath[1] > 1){
        gIndex1 = anchorRes$gIndex[1]
        gIndex2 = anchorRes$gIndex[selAnchorPath[1]]
        
      }else{next;}
    }else if(i == length(selAnchorPath)){
      if(selAnchorPath[i] < nrow(anchorRes)){
        gIndex1 = anchorRes$gIndex[selAnchorPath[i]]
        gIndex2 = anchorRes$gIndex[nrow(anchorRes)]
      }else{next;}
    }else{
      gIndex1 = anchorRes$gIndex[selAnchorPath[i]]
      gIndex2 = anchorRes$gIndex[selAnchorPath[i+1]]
    }
    
    goodAnchorIndex = subset(goodAnchors,gIndex >= gIndex1 & gIndex <= gIndex2)$gIndex
    tmpGroupSm = pwGroupSm[goodAnchorIndex, goodAnchorIndex]
    
    distToAnchors = tmpGroupSm[1,] + tmpGroupSm[ncol(tmpGroupSm),]
    stDist = stDist + lastDist
    lastDist = tmpGroupSm[1,ncol(tmpGroupSm)]
    tmpGroupSm = tmpGroupSm / (distToAnchors) * tmpGroupSm[1,ncol(tmpGroupSm)]
    groupDist = tmpGroupSm[,1]
    groupDist[which(distToAnchors == 0)] = 0
    groupDist = groupDist + stDist
    if(i == 0 | (i==1 & selAnchorPath[1] == 1)){
      allGroupDist = c(allGroupDist,groupDist)
    }else{
      allGroupDist = c(allGroupDist,groupDist[2:length(groupDist)])
    }
    
  }
  return(allGroupDist)
}


#'  This function is to interpolate position
#'
#' @param x vector<numeric>  a numeric vector with no missing values (e.g. physical positions on a chr)
#' @param y vector<numeric> a numeric vector the same length as x (e.g. genetic positions on a chr)
#' @param tol numeric a defalut value
#' 
#' @return data.frame with interpolated values

interpolate <- function(x, y, tol=1e-5) {
  #x is a numeric vector with no missing values (e.g. physical positions on a chr)
  #y is a numeric vector the same length as x (e.g. genetic positions on a chr)
  #
  #the result is, after sorting y in the same order as x, interpolating the missing values of y 
  #based on the nearest non-missing values of y and their corresponding values of x
  if (!is.vector(x) || !is.vector(y) || length(x) != length(y) || !is.numeric(x) || !is.numeric(y) || any(is.na(x))) 
    stop("x and y need to be numeric vectors of the same length without any missing values in x")
  xsorted <- sort(x)
  xorder <- order(x)
  ysorted <- y[xorder]
  yna <- is.na(ysorted)
  if (sum(yna) == 0) return(cbind(x=xsorted,y.in=ysorted,y.out=ysorted))
  ynai <- which(yna)
  nyna <- which(!yna)
  starts <- sapply(ynai, function(z) max(nyna[nyna < z]))
  stops  <- sapply(ynai, function(z) min(nyna[nyna > z]))
  diffsy <- ysorted[stops] - ysorted[starts]
  diffsx <- xsorted[stops] - xsorted[starts]
  ratios <- ifelse(abs(diffsx) > tol, diffsy / diffsx, 0)
  interp <- sapply(1:sum(yna), function(z) ysorted[starts[z]] + ratios[z]*(xsorted[ynai[z]] - xsorted[starts[z]]))
  yout <- ysorted
  yout[ynai] <- interp
  #if(y[min(nyna)] == 0 & min(nyna) > 1) {yout[1:(min(nyna)-1)] <- 0} # the first overlap is 0, then assigned any position before the first match to 0
  out <- cbind(x=xsorted,y.in=ysorted,y.out=yout)
  out
}

#'  This function is to interpolate marker positions based on goodAnchors' position
#'
#' @param subChrData data.frame
#' @param goodAnchors data.frame selected anchors that are in good quality in terms of data point, maxDistance, etc
#' 
#' @return data.frame

getMarkerPos = function(subChrData,goodAnchors){
  # interpolate marker positions based on goodAnchors' position
  subChrData$newMapPos = NA
  currMapLen = 0
  
  for(i in 1:nrow(goodAnchors)){
    stIndex = goodAnchors$gIndex[i]
    anchorMkIndex = stIndex:(stIndex + 9)
    
    if(i == 1){
      # the first 5 has half of the max distance; assume linear and the anchor performs a single unit
      subChrData$newMapPos[anchorMkIndex[1:5]] = seq(0,goodAnchors$max[i]/2,length.out = 5)
      currMapLen = currMapLen + goodAnchors$max[i]/2
    }else if(i == nrow(goodAnchors)){
      # the last 5 has half of the max distance for the last anchor; 
      subChrData$newMapPos[anchorMkIndex[6:10]] = seq(0,goodAnchors$max[i]/2,length.out = 5) + currMapLen
    }else{
      stIndex0 = goodAnchors$gIndex[i-1]
      anchorMkIndex0 = stIndex0:(stIndex0 + 9)
      
      newMkIndexSt = anchorMkIndex0[5]
      newMkIndexEd = anchorMkIndex[5]
      
      anchorDist = goodAnchors$groupPos[i] - goodAnchors$groupPos[i-1]
      
      
      if(goodAnchors$max[i-1] == 0){
        subChrData$newMapPos[anchorMkIndex0[6:10]] = 0 + currMapLen
        newMkIndexSt = max(anchorMkIndex0)
      }
      
      if(goodAnchors$max[i] == 0){
        subChrData$newMapPos[anchorMkIndex[1:5]] = 0 + currMapLen + anchorDist
        newMkIndexEd = min(anchorMkIndex)
      }
      
      if(newMkIndexSt < newMkIndexEd){
        x = subChrData$Sbegin[newMkIndexSt:newMkIndexEd]
        y = rep(NA,length(x)); y[1] = 0;y[length(x)] = anchorDist
        interpolatePos = interpolate(x,y)
        subChrData$newMapPos[newMkIndexSt:newMkIndexEd] = interpolatePos[,3] + currMapLen
      }else{
        #do nothing; 
      }
      currMapLen = currMapLen + anchorDist
      
    }
    
  }
  return(subChrData)
}


#'  This function is to calculate marker distance
#'
#' @param goodAnchors data.frame selected anchors that are in good quality in terms of data point, maxDistance, etc.
#' @param index integer
#' @param chrMapDist data.frame pairwise distance
#' @param newMarkerRes data.frame marker information with genetic map information
#' @param minDataWithAnchor integer minimum data for a given anchor
#' 
#' @return data.frame

reCalMkDist = function(goodAnchors,index,chrMapDist,newMarkerRes,minDataWithAnchor = 3){
  # for the adjacent anchors with big gaps,instead of interpolate by linear, recalculate the marker distance; 
  #
  index1 = goodAnchors$gIndex[index]
  index2 = goodAnchors$gIndex[index+1]
  
  anchorMkIndex = index1:(index1 + 9)
  anchorMkIndex2 = index2:(index2 + 9)
  
  a1Pos = newMarkerRes$newMapPos[index1+9]
  a2Pos = newMarkerRes$newMapPos[index2]
  
  if(index2 - index1 > 20){ # at least have 10 markers between two anchors
    tmpMkDist = as.numeric()
    for(markerIndex in (index1 + 10):(index2 - 1)){
      tmpDist1 = ifelse(length(which(!is.na(chrMapDist[markerIndex,anchorMkIndex]))) < minDataWithAnchor, NA, mean(chrMapDist[markerIndex,anchorMkIndex],na.rm=T))
      tmpDist2 = ifelse(length(which(!is.na(chrMapDist[markerIndex,anchorMkIndex2]))) < minDataWithAnchor, NA, mean(chrMapDist[markerIndex,anchorMkIndex2],na.rm=T))
      tmpDist = ifelse(tmpDist1 + tmpDist2 == 0, a1Pos, a1Pos + tmpDist1 / (tmpDist1 + tmpDist2) * (a2Pos - a1Pos))
      tmpMkDist = c(tmpMkDist,tmpDist)
    }
    names(tmpMkDist) = (index1 + 10):(index2 - 1)
    index3 = which(!is.na(tmpMkDist))
    tmpMkDist2 = tmpMkDist[index3]
    if(length(tmpMkDist2) < 1){
      return(NA)
    }else if(length(tmpMkDist2) > 1){
      distGap = tmpMkDist2[2:length(tmpMkDist2)] - tmpMkDist2[1:(length(tmpMkDist2) - 1)]
      
      while(any(distGap < 0)){
        if(length(tmpMkDist2) == 2){
          tmpSwap = names(tmpMkDist2)
          tmpMkDist2 = c(tmpMkDist2[2],tmpMkDist2[1])
          names(tmpMkDist2) = tmpSwap
          break
        }
        tmpMkDist2 = tmpMkDist2[-(which(distGap < 0) + 1)]
        distGap = tmpMkDist2[2:length(tmpMkDist2)] - tmpMkDist2[1:(length(tmpMkDist2) - 1)]
      }
      return(tmpMkDist2)
      
      
    }else{
      return(tmpMkDist2)
    }
    
  }else{
    return(NA)
  }
  
}

### 5. pair-wise marker kosambi distances

#'  This function is to calculate pairwise kosambi distance for the trio data
#'
#' @param testData data.frame having at least one column named "Q_ID"
#' @param totalTrioCount data.frame a table format that includes all the trios for each pair of marker (row and column)
#' @param totalRecRate data.frame a table format that includes all the recombination rate for each pair of marker (row and column)
#' @param minCount integer minimum trio data required to calculate kosambi distance
#'
#' @return data.frame

getPwKsbDist = function(testData,totalTrioCount, totalRecRate,minCount = 100){
  testMarkers  = testData$Q_ID
  testMarkers = testMarkers[testMarkers %in% colnames(totalTrioCount)]
  subRecRate = totalRecRate[testMarkers,testMarkers]
  subTotalCount = totalTrioCount[testMarkers,testMarkers]
  subTotalCount = ifelse(subTotalCount >= minCount,1,NA)
  subRecRate = as.matrix(subRecRate * subTotalCount)
  subRecRate = ifelse(subRecRate >= 0.5,0.49,subRecRate)
  subRecRate = subRecRate / (2 - 2 * subRecRate)
  chrMapDist = 0.25 * log((1+2*subRecRate) / (1 - 2*subRecRate)) * 100 #kosambi distance; need to set
  chrMapDist = ifelse(chrMapDist >= 50,50,round(chrMapDist,2)) #kosambi distance; need to set 
  return(chrMapDist)
}


### 6. others

#'  This function is to find the offPos markers
#'
#' @param chrMapDist data.frame pairwise distance
#' @param minCM numeric minimum genetic distance
#' @param maxDist numeirc maximum genetic distance
#'
#' @return data.frame

findOffPosMk = function(chrMapDist,minCM=0.2,maxDist = 150){
  # offPosMk: markers that are tightly linked with other markers that are farther way and/or 
  # markers that have no closely linked markers
  N = nrow(chrMapDist)
  pbMks = data.frame(mkName = row.names(chrMapDist), maxDistMk = rep(NA,N), pbMks = rep("F",N))
  
  for(i in 1:N){
    posIndex = which(as.numeric(chrMapDist[i,]) <= minCM)
    if(length(posIndex) < 1){
      pbMks[i,3] = "T"
    }else{
      pbMks[i,2] = max(abs(posIndex - i))
      if(max(abs(posIndex - i)) >= maxDist){
        pbMks[i,3] = "T"
      }
    }
  }
  return(pbMks)
}

#'  This function is to get all markers that have recFraction data
#'
#' @param tmpFolder string the path that stores all the tmp data for each trio
#'
#' @return data.frame

collectAllInfoMarkers <- function(tmpFolder){
  # collect all infoMarkerByLine results
  
  allFiles = list.files(tmpFolder)
  lineRes = data.frame(line=as.character(),p1Count = as.numeric(),p2Count = as.numeric())

  for(fileName in allFiles){
    tmpData = read.csv(paste0(tmpFolder,"/",fileName))
    qLine = unlist(strsplit(fileName,"\\."))[1]
    tmpData$line = qLine
    lineRes = bind_rows(lineRes,data.frame(line = qLine,p1Count = table(tmpData$parent)[1],p2Count = table(tmpData$parent)[2]))
  }
  return(lineRes)
}
