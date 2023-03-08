#### functions ####
# These are functions used for trio mapping

### 1. genotype data (table): Convert fpString to table format
convertFpStringToGeno <- function(df,fpGeno,markerMap){
  # this function is to convert fpString to a table format, with row as the number of markers and column as the number of lines;
  # df is the file with pedigree information; containing at least three columns: "pedigree", "prt1", "prt2"
  # fpGeno is the fpString collected by Danny and team
  # markerMap is the genetic map for all the markers used in fpGeno
  # takes about 3 mins to generate the fpMatrix for soybean: 67717 markers x 15634 lines

  lines <- unique(c(df$pedigree,df$prt1,df$prt2))
  fpDataVec <- data.frame(matrix("-",ncol=length(lines),nrow=nrow(markerMap)))
  colnames(fpDataVec) <- lines
  row.names(fpDataVec) <- markerMap$markerName
  count = 0
  for(line in lines){
    count = count + 1
    str <- str_to_upper(fpGeno[line,"V2"])
    fpDataVec[,line] <- unlist(str_split(str,""))
  }
  return(fpDataVec)
}

### 2. recombination rate
infoMarkerByLine <- function(df,fpDataVec,nCores,tmpFolder){
  # This function is to find informative markers for each trio and then return its status by line: P1 or P2
  # data will be saved by each line
  # missing data and het call in any parent will not be considered as informative

  numCores <- detectCores()
  if(nCores < numCores){
    registerDoParallel(nCores)
    print(paste("Using",nCores,"Cores out of the computer's maximum of", numCores,"cores"))
  }else{
    nCores = numCores
    registerDoParallel(numCores)
    print(paste("Using",nCores,"Cores out of the computer's maximum of", numCores,"cores"))
  }
  
  linesPerRun <- ceiling(nrow(df) / nCores)
  
  foreach(count = 1:nCores,.combine=rbind) %dopar% {
    
    st <- (count-1) * linesPerRun + 1
    ed <- min(count*linesPerRun,nrow(df))
    
    subDf <- df[st:ed,]
    
    for(k in 1:nrow(subDf)){
      qLine <- subDf$pedigree[k]
      prts <- unlist(str_split(subDf$parents[k],";"))
      if(all(c(qLine,prts) %in% colnames(fpDataVec))){
        tmpData <- fpDataVec[,c(qLine,prts)]
        colnames(tmpData) <- c("O","P1","P2")
        
        # 1) skip if any missing data or heterozygous presents
        tmpData <- subset(tmpData,!str_detect(O,"\\.|[hHnN]") & !str_detect(P1,"\\.|[hHnN]") & !str_detect(P2,"\\.|[hHnN]") & P1 != P2)
        tmpDataMarkers <- row.names(tmpData)
        tmpData$markers = tmpDataMarkers
        tmpData$parent <- ifelse(tmpData$O == tmpData$P1,1,2)
        
        # 2) skip if no data or only one parental genotype
        if(nrow(tmpData) <= 1){next;}
        if(length(unique(tmpData$parent)) == 1){next;}
        
        write.csv(tmpData,paste0(tmpFolder,"/",qLine,".csv"),row.names=F)
      }
    }
  }
}

# get all markers that have recFraction data
collectAllInfoMarkers <- function(tmpFolder){
  # collect all infoMarkerByLine results
  
  allFiles = list.files(tmpFolder)
  lineRes = data.frame()

  for(fileName in allFiles){
    tmpData = read.csv(fileName)
    qLine = unlist(strsplit(fileName,"\\."))[1]
    tmpData$line = qLine
    lineRes = bind_rows(lineRes,tmpData)
  }
  return(lineRes)
}

# for any pair of markers, count #totalTros and #recombinedTrios
# Output are two matrix, each one is a large matrix with the dimension as #markers x #markers (around 65k x 65k)

generateRecData = function(tmpFolder,markerMap,outFileFolder,verbose = TRUE){
  # This function is to generate recbinate rate for any pair of markers
  # For every trio, collect the informative markers and found which pair of markers are recombined (saved at tmpFolder)
  # collect all infoMarkerByLine results, and then count recombined lines and non-recombined lines;
  
  allMarkers = markerMap$markerName
  
  testResEffCount <- data.frame(matrix(0,nrow=length(allMarkers),ncol=length(allMarkers)))
  row.names(testResEffCount) <- allMarkers
  colnames(testResEffCount) <- allMarkers
  
  testResRecCount <- data.frame(matrix(0,nrow=length(allMarkers),ncol=length(allMarkers)))
  row.names(testResRecCount) <- allMarkers
  colnames(testResRecCount) <- allMarkers
  
  allFiles = list.files(tmpFolder)
  testedMarkers= as.character()
  count = 0
  
  for(fileName in allFiles){
    count = count + 1
    if((count-1) %% 1000 == 0 & verbose){print(paste("Processed",count,"lines;",Sys.time()))}
    tmpData = read.csv(fileName)
    qLine = unlist(strsplit(fileName,"\\."))[1]

    testedMarkers = unique(c(testedMarkers,tmpData$markers))
    
    tmpDataMarkers = tmpData$markers
    
    testResEffCount[tmpDataMarkers,tmpDataMarkers] <- testResEffCount[tmpDataMarkers,tmpDataMarkers] + 1 # sum recombined and none-recombined trios
    
    p1Index <- which(tmpData$parent == 1)
    p2Index <- which(tmpData$parent == 2)
    
    # recombination happened between markers having different parent
    testResRecCount[tmpDataMarkers[p1Index],tmpDataMarkers[p2Index]] <- testResRecCount[tmpDataMarkers[p1Index],tmpDataMarkers[p2Index]] + 1 # only recombined trios
    testResRecCount[tmpDataMarkers[p2Index],tmpDataMarkers[p1Index]] <- testResRecCount[tmpDataMarkers[p2Index],tmpDataMarkers[p1Index]] + 1 # only recombined trios
  }
  
  testResRecCount = testResRecCount[testedMarkers,testedMarkers]
  testResEffCount = testResEffCount[testedMarkers,testedMarkers]
  
  cat(paste("... Saving Data ....;", "\n",Sys.time()))
  saveRDS(testResEffCount,paste0(outFileFolder,"/totalTrioCount_allChr.rds"))
  saveRDS(testResRecCount,paste0(outFileFolder,"/recTrioCount_allChr.rds"))
  testRecRate <- round(testResRecCount/testResEffCount,4) # recombination rate for each markers (#recTrios / #totalTrios)
  saveRDS(testRecRate,paste0(outFileFolder,"/recRate_allChr.rds"))
  cat(paste("... Rec Data Done!!! ....;", "\n",Sys.time()))
}


### 3. define groups and determine anchors
# find closely linked markers
getLinkedMarkerInfo <- function(phyPos,testMarkers, testTrioCount, totalRecRate, markerMap,minCount, maxRecRate){
  # For each marker, find other markers that linked with this marker and summarize the linkage information;
  # return summarized information
  
  res <- phyPos %>% filter(Q_ID %in% testMarkers) %>% select(Q_ID, chr, Sbegin, Send)
  row.names(res) = res$Q_ID
  res = res[testMarkers,]
  res$chrDist <- NA
  res$numPerChr <- NA
  res$totalMkLinkedOwnChr <- NA
  res$totalMkLinked <- NA
  res$LDMarkers <- NA
  
  for(pos in 1:nrow(res)){
    tmp <- data.frame(chr = res$chr,pos = res$Sbegin,totalCount = testTrioCount[testMarkers,res$Q_ID[pos]],
                      recRate = totalRecRate[testMarkers,res$Q_ID[pos]]) 
    tmp$marker = testMarkers
    tmp <- tmp %>% filter(totalCount >= minCount, recRate <=  maxRecRate)
    tableTmp <- table(tmp$chr) 
    chrDist <- paste(names(tableTmp),collapse =";")
    chrNum <- paste(tableTmp,collapse =";")
    res$chrDist[pos] <- chrDist
    res$numPerChr[pos] <- chrNum
    res$totalMkLinkedOwnChr[pos] = tableTmp[as.character(res$chr[pos])]
    res$totalMkLinked[pos] = sum(tableTmp)
    res$LDMarkers[pos] = paste(tmp$marker,collapse =";")
    
  }
  
  res$numLinkedChr = sapply(res$chrDist,function(x){return(length(unlist(strsplit(x,";"))))})
  res$index = 1:nrow(res)
  
  res = left_join(res,markerMap[,c(2:6)],by=c("Q_ID" = "markerName"))
  return(res)
}

# find the offPos markers
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

getGroupAnchorDiff <- function(anchors, pwDist = subPwGroupSm){
  # this function is to evaluate anchors based on MAE; 
  # For any pair of anchor, MAE is calculated from groups from their left (after last anchor), middle and right (before last anchor)
  # anchors is a output from the function of "findLinkedGroups"
  # its a list with two components: distVec and selIndex; 
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

findAnchorsByChr = function(testChr, res_2, testTrioCount, totalRecRate) {
  
  minCount = 200
  testChrData = res_2 %>% filter(chr == testChr, chrDist == testChr, totalMkLinkedOwnChr >= 5) %>% arrange(Sbegin)
  testChrMk = testChrData$Q_ID
  chrRecCount = testTrioCount[testChrMk,testChrMk]
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
  write.csv(subChrData,paste0("/mnt/tmpData/tmpAnchorData/subChrData_chr",testChr,".csv"))
  
  subChrMapDist = chrMapDist[subChrData$Q_ID,subChrData$Q_ID]
  saveRDS(subChrMapDist,paste0("/mnt/tmpData/tmpAnchorData/subChrMapDist_chr",testChr,".rds"))
  
  numOfGroups = nrow(subChrData) - 9
  groupSm = data.frame(gIndex = 1:numOfGroups,numOfData = rep(NA,numOfGroups),minDataPerMk = rep(NA,numOfGroups),min=rep(NA,numOfGroups),
                       max = rep(NA,numOfGroups), mean = rep(NA,numOfGroups), median = rep(NA,numOfGroups),std = rep(NA,numOfGroups))
  for(i in 1:numOfGroups){
    st = i
    ed = i + 9
    posIndex = st:ed
    subMapDist = subChrMapDist[posIndex,posIndex]
    groupSm[i,] = c(i,length(which(subMapDist >= 0)),min(apply(subMapDist,1,function(x){length(which(x >=0))})),min(subMapDist,na.rm = T),
                    max(subMapDist,na.rm = T),mean(as.matrix(subMapDist),na.rm = T),median(as.matrix(subMapDist),na.rm = T),sd(as.matrix(subMapDist),na.rm = T))
    
  }
  
  write.csv(groupSm,paste0("/mnt/tmpData/tmpAnchorData/groupSm_chr",testChr,".csv"))
  
  pwGroupSm = matrix(NA,nrow=numOfGroups,ncol=numOfGroups)
  
  for(i in 1:numOfGroups){
    st1 = i
    ed1 = i + 9
    for(j in i:numOfGroups){
      st2 = j
      ed2 = j + 9
      if(i == j){
        pwGroupSm[i,j] = 0
        pwGroupSm[j,i] = 0
      }else{
        subMapDist = subChrMapDist[st1:ed1,st2:ed2]
        if(length(which(subMapDist >= 0)) >= 20){
          pwGroupSm[i,j] = mean(subMapDist,na.rm=T)
          pwGroupSm[j,i] = mean(subMapDist,na.rm=T)
        }
      }
    }
    
  }
  
  #print(levelplot(pwGroupSm, main=paste0("pwGroupSm: Chr",testChr, " (>= 20% data)"), xlab="", ylab=""))
  saveRDS(pwGroupSm,paste0("/mnt/tmpData/tmpAnchorData/pwGroupSm_chr",testChr,".rds"))
  
  # pick groups that are most reliable
  selGroupSm = groupSm %>% filter(numOfData >= 80, minDataPerMk >= 5, max <= 0.5)
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
  
  write.csv(anchorRes,paste0("/mnt/tmpData/tmpAnchorData/anchorRes_chr",testChr,".csv"))
}


### 4. anchor distances and marker positions

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

getPwKsbDist = function(testData,testTrioCount, totalRecRate,minCount = 100){
  testMarkers  = testData$Q_ID
  testMarkers = testMarkers[testMarkers %in% colnames(testTrioCount)]
  subRecRate = totalRecRate[testMarkers,testMarkers]
  subTotalCount = testTrioCount[testMarkers,testMarkers]
  subTotalCount = ifelse(subTotalCount >= minCount,1,NA)
  subRecRate = as.matrix(subRecRate * subTotalCount)
  subRecRate = ifelse(subRecRate >= 0.5,0.49,subRecRate)
  subRecRate = subRecRate / (2 - 2 * subRecRate)
  chrMapDist = 0.25 * log((1+2*subRecRate) / (1 - 2*subRecRate)) * 100 #kosambi distance; need to set
  chrMapDist = ifelse(chrMapDist >= 50,50,round(chrMapDist,2)) #kosambi distance; need to set 
  return(chrMapDist)
}



