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