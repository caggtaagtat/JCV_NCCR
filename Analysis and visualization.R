setwd("C:/Users/joh/Desktop/HomeOffice/JC")

##Load distance 
library(VarCon)
library(usedist)
library(clValid)
library(cluster)
library(DECIPHER)
library(seqinr)
library(ips)

library(ggtreeExtra)
library(gridExtra)
library(ggplot2)
library(ape)
library(msa)
library(ggtree)


files <- list.files(path="C:/Users/joh/Desktop/HomeOffice/JC/results/", pattern=NULL, all.files=FALSE, 
                    full.names=T)

res <- list()

for(i in files){
   
   # i <- files[1]
   curFile <- read.csv2(i)
   curFileName <- gsub("C:/Users/joh/Desktop/HomeOffice/JC/results/","", i)
   curFileName <- strsplit(curFileName, "-")[[1]]
   curFileName <- paste0(curFileName[1],"-", curFileName[2])
   curFile$length_consensus_seq <- nchar(curFile$consensus_seq)
   curFile$sample <- curFileName
   
   curFile$sample_totalPercentage <- sum(curFile$read_percentage)
   res[[curFileName]] <- curFile
}


library(data.table)
resAll <- rbindlist(res)

resAll$ClusterID <- paste0(resAll$Cluster, "_", resAll$sample)

## Create DNAStringSet from cluster consensus sequences
resAll <- resAll[order(resAll$length_consensus_seq, decreasing = T),]
ttSeq <- resAll$consensus_seq
names(ttSeq) <- resAll$ClusterID
ttDNAstringset <- DNAStringSet(ttSeq)

ttDNAstringsetaln <- msa(ttDNAstringset, method="Muscle")

ttDNAstringsetalnD <- DNAStringSet(ttDNAstringsetaln)

ttdist <- DistanceMatrix(ttDNAstringsetalnD, includeTerminalGaps=T, type ="dist")

resAll$length <- nchar(resAll$consensus_seq)


save(resAll, file="resAll")

### Pairwise alignemnt for changeing cluster on HPC

library(parallel)
library(data.table)
library(ips)
library(tidyverse)
library(msa)
library(DECIPHER)


load("clusterSave")

cl <- makeCluster(50)


clusterExport(cl, c("resAll", "msa", "DNAStringSet", "DistanceMatrix", "msaMuscle"))

allAlign <- parLapply(cl, resAll$ClusterID, function(clusterID){
   
   #clusterID <- resAll$ClusterID[1]
   curAl <- data.frame(id=clusterID, partner=resAll$ClusterID[resAll$ClusterID != clusterID], dist= 999)
   
   for(i in 1:nrow(curAl)){
      
      ttDNAstringset <- DNAStringSet(c(resAll$consensus_seq[resAll$ClusterID %in% c(curAl$id[i], curAl$partner[i])]))
      ttDNAstringsetaln <- msa(ttDNAstringset, method="Muscle")
      ttDNAstringsetalnD <- DNAStringSet(ttDNAstringsetaln)
      
      ttdist <- DistanceMatrix(ttDNAstringsetalnD, type="dist", includeTerminalGaps=T)
      curAl$dist[i] <- ttdist[[1]]
      
   }
   
   sort(c(clusterID, curAl$partner[curAl$dist == 0]))
   
   
})

save(allAlign, file="allAlign")







## Continue on mashien
load("allAlign")
load("resAll")
allAlign <- unique(allAlign)

changeClusters <- list()

for(i in 1:length(allAlign)){
   
   curGlob <- data.frame(cluster=i, name=allAlign[[i]])
   changeClusters[[i]] <- curGlob
   
}

changeClusters <- rbindlist(changeClusters)


clusterTable <- data.frame(table(changeClusters$cluster))
clusterTable$Var1 <- as.character(clusterTable$Var1)
clusterTable$oldClusterID <- ""

names(clusterTable) <- c("ClusterID", "frequency", "oldClusterID")
clusterTable$oldClustersPercentage <- ""

for(i in 1:nrow(clusterTable)){
   
   curCluster <- clusterTable$ClusterID[i]
   changedClusters <- changeClusters$name[changeClusters$cluster == curCluster]
   
   clusterTable$oldClusterID[i] <- paste(changedClusters, collapse =", ")
   clusterTable$oldClustersPercentage[i] <- paste0(changedClusters, "(", resAll$read_percentage[match(changedClusters, resAll$ClusterID) ], "%)", collapse =", ")
   
}

resAll$globalCluster <- changeClusters$cluster[match(resAll$ClusterID, changeClusters$name)]

##mark archetype NCCR sequence 
resAll$globalCluster <- as.character(resAll$globalCluster)
resAll$globalCluster[resAll$consensus_seq == "GCCTCGGCCTCCTGTATATATAAAAAAAAGGGAAGGTAGGGAGGAGCTGGCTAAAACTGGATGGCTGCCAGCCAAGCATGAGCTCATACCTAGGGAGCCAACCAGCTGACAGCCAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCCCTCCCATAGTCCTTAATCACAAGTAAACAAAGCACAAGGGGAAGTGGAAAGCAGCCAAGGGAACATGTTTTGCGAGCCAGAGCTGTTTTGGCTTGTCACCAGCTGGCC"] <- "Archetype"

#clusterID <- resAll$ClusterID[1]
curAl <- data.frame(id="ar", partner= resAll$ClusterID, dist= 9)

for(i in 1:nrow(curAl)){
   
   parseq <- resAll$consensus_seq[resAll$ClusterID == curAl$partner[i]]
   if(nchar(parseq) > 200){
      ttDNAstringset <- DNAStringSet(c(ar, parseq))
      ttDNAstringsetaln <- msa(ttDNAstringset, method="Muscle")
      ttDNAstringsetalnD <- DNAStringSet(ttDNAstringsetaln)
      
      ttdist <- DistanceMatrix(ttDNAstringsetalnD, type="dist", includeTerminalGaps=T)
      curAl$dist[i] <- ttdist[[1]]
   }
}


resAll$globalCluster[resAll$ClusterID %in% curAl$partner[curAl$dist <= 0.01]] <-  "Archetype"

library("xlsx")
# Write the first data set in a new workbook
write.xlsx(resAll, file = "JC_Haplotypes2.xlsx",
           sheetName = "Sample Clusters", append = FALSE)

# Add a second data set in a new worksheet
write.xlsx(clusterTable, file = "JC_Haplotypes2.xlsx", 
           sheetName="Global Cluster Table", append=TRUE)

#save(resAll, file="resAll2")









## Test 
load("resAll2")
samples <- read.csv("C:/Users/joh/Desktop/HomeOffice/JC/JC NCCR für Johannes Ptok.csv", sep=";")
samples$Donor <- as.character(samples$Donor)
samples$DonorID <- as.character(lapply(samples$Donor, function(x){ strsplit(x, "\\.")[[1]][[1]]}))
samples$Sample <- gsub("D1", "D-1", samples$Sample)
samples$Sample <- gsub(" \\(Test\\)", "", samples$Sample)


resAll$DonorID <- samples$DonorID[match(resAll$sample, samples$Sample)]
resAll$DonorTime <- samples$Date[match(resAll$sample, samples$Sample)]
resAll$sampleType <- samples$Type[match(resAll$sample, samples$Sample)]

## mache datum zu datum


library(tidyr)
pacman::p_load(lubridate)

resAll$DonorTime2 <- resAll$DonorTime
resAll$DonorTime2 <- gsub("2020", "20", resAll$DonorTime2)
resAll$DonorTime2 <- gsub("2021", "21", resAll$DonorTime2)
resAll$DonorTime2 <- gsub("2022", "22", resAll$DonorTime2)
resAll$DonorTime2 <- gsub("2023", "23", resAll$DonorTime2)
resAll$DonorTime2 <- gsub("2024", "24", resAll$DonorTime2)

resAll$DonorTime <- as.character(gsub("\\.", "-",resAll$DonorTime))
resAll$DonorTime <- dmy(resAll$DonorTime)

resAll <- resAll[order(resAll$DonorTime, decreasing = F),]


resAll <- resAll[which(resAll$sample %in% samples$Sample),]
samples <- samples[which(samples$DonorID %in% resAll$DonorID),]

resAllSave <- resAll


##Function to visualize NCCR segment composition

get_table <- function(inputSequence){
   
   #inputSequence <- haplotypes$sequence[haplotypes$name ==  "Haplotype_13_maxFreq_1.5%"]
   #inputSequence <- readClipboard()
   
   ori <- "GCCTCGGCCTCC"
   a <- "CCTGTATATATAAAAAAAAGGGAAG"
   b <- "GTAGGGAGGAGCTGGCTAAAACT"
   c <- "GGATGGCTGCCAGCCAAGCATGAGCTCATACCTAGGGAGCCAACCAGCTGACAGC"
   d <- "CAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCCCTCCCATAGTCCTTAATCAC"
   e <- "AAGTAAACAAAGCACAAG"
   f <- "GGGAAGTGGAAAGCAGCCAAGGGAACATGTTTTGCGAGCCAGAGCTGTTTTGGCTTGTCACCAGCTGGCC"
   
   if(as.numeric(regexpr("GAAACCACTCCCA" , inputSequence)) != -1){
      
      f <- "GGGAAGTGGAAAGCAGCCAGGGGAACATGTTTTGCGAGCCAGAGCTGTTTTGGCTTGTCACCAGCTGGCC"
      d <- "CAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCACTCCCATAGTCCTTAATCAC"
   }
   
   getOverlappingVectorsFromVector2 <- function(largeVector, subvectorLength, subvectorOverlap) 
   {
      startPositions = seq(1, length(largeVector), by = subvectorLength - 
                              subvectorOverlap)
      end_positions = startPositions + subvectorLength - 1
      end_positions[end_positions > length(largeVector)] = length(largeVector)
      lapply(seq_len(length(startPositions)), function(x) largeVector[startPositions[x]:end_positions[x]])
   }
   
   
   
   creatSubreadVec <- function(vec, len=10){
      
      vec <- strsplit(vec, "")[[1]]
      vec <- getOverlappingVectorsFromVector2(vec, 10, 9)
      vec <- as.character(lapply(vec, function(x){paste0(x, collapse="")} ))
      vec <- vec[nchar(vec) == 10]
      return(vec)
   }
   
   oriV <- creatSubreadVec(ori)
   aV <- creatSubreadVec(a)
   bV <- creatSubreadVec(b)
   cV <- creatSubreadVec(c)
   dV <- creatSubreadVec(d)
   eV <- creatSubreadVec(e)
   fV <- creatSubreadVec(f)
   
   nameVec <- c(rep("ori", length(oriV)),
                rep("a", length(aV)),
                rep("b", length(bV)),
                rep("c", length(cV)),
                rep("d", length(dV)),
                rep("e", length(eV)),
                rep("f", length(fV)))
   
   nameVecPos <- c(seq(1, length(oriV), 1),
                   seq(1, length(aV), 1),
                   seq(1, length(bV), 1),
                   seq(1, length(cV), 1),
                   seq(1, length(dV), 1),
                   seq(1, length(eV), 1),
                   seq(1, length(fV), 1))
   
   multimer <- data.frame(name= nameVec, pos=nameVecPos, 
                          multimer=c(oriV, aV, bV, cV, dV, eV, fV))
   
   
   ori2 <- substr(ori,nchar(ori)-8 ,nchar(ori))
   a1 <-  substr(a, 1, 9)
   b1 <-  substr(b, 1, 9)
   c1 <-  substr(c, 1, 9)
   d1 <-  substr(d, 1, 9)
   e1 <-  substr(e, 1, 9)
   f1 <-  substr(f, 1, 9)
   
   a2 <-  substr(a,nchar(a)-8 ,nchar(a))
   b2 <-  substr(b,nchar(b)-8 ,nchar(b))
   c2 <-  substr(c,nchar(c)-8 ,nchar(c))
   d2 <-  substr(d,nchar(d)-8 ,nchar(d))
   e2 <-  substr(e,nchar(e)-8 ,nchar(e))
   
   oa <- paste0(ori2,a1)
   ab <- paste0(a2,b1)
   bc <- paste0(b2,c1)
   cd <- paste0(c2,d1)
   de <- paste0(d2,e1)
   ef <- paste0(e2,f1)
   
   oaV <- creatSubreadVec(oa)
   abV <- creatSubreadVec(ab)
   bcV <- creatSubreadVec(bc)
   cdV <- creatSubreadVec(cd)
   deV <- creatSubreadVec(de)
   efV <- creatSubreadVec(ef)
   
   
   nameVec <- c(rep("ori a", length(oaV)),
                rep("a b", length(abV)),
                rep("b c", length(bcV)),
                rep("c d", length(cdV)),
                rep("d e", length(deV)),
                rep("e f", length(efV)))
   
   
   multimerBorders  <- data.frame(name= nameVec, pos= NA, 
                                  multimer=c(oaV, abV, bcV, cdV, deV, efV))
   
   creatPlottingData <- function(promotorSeq){
      
      # promotorSeq <- inputSequence
      
      if(nchar(promotorSeq)<20 )
         stop("ERROR during setting of variable 'promotorSeq'. The entered sequence", 
              " must be at least 20nt long.")  
      if (!all(strsplit(promotorSeq, "")[[1]] %in% c("a", "c", "g", "t", 
                                                     "G", "C", "T", "A"))) 
         stop("ERROR during setting of variable 'promotorSeq'. The entered sequence", 
              " must be a character string of A C G and T.")    
      
      promotorSeq <- toupper(promotorSeq)
      wtV <- creatSubreadVec(promotorSeq)
      
      wtDF <- data.frame(pos= 0, seq=wtV)
      
      wtDF$segment <- multimer$name[match(wtDF$seq, multimer$multimer)]
      wtDF$segmentPos <- multimer$pos[match(wtDF$seq, multimer$multimer)]
      
      wtDF$pos <- 1:nrow(wtDF)
      
      i <- 1
      while(i < nrow(wtDF)){
         
         while(!is.na(wtDF$segment[i]) & i < nrow(wtDF)) i <- i+1
         
         counter <- list()
         
         while( is.na(wtDF$segment[i]) & ((i+1) <= nrow(wtDF))){
            
            counter[[i]] <- i
            i <- i+1
         }
         
         counter <- unlist(counter)
         
         if(length(counter) > 9 ){
            wtDF$segment[counter] <-  multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)] ##test
            
            if(any(!is.na(multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)]))){
               
               pos <- which(!is.na(multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)]))
               
               before10mer <- wtDF$seq[(sort(counter)[1]-1)] 
               after10mer <- wtDF$seq[(sort(counter)[length(counter)]+1)]  
               
               if(before10mer == oriV[length(oriV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "ori a"] <- "a"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == aV[length(aV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "a b"] <- "b"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == bV[length(bV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "b c"] <- "c"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == cV[length(cV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "c d"] <- "d"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == dV[length(dV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "d e"] <- "e"
                  wtDF$segment[counter[1:9]] <- test
               }   
               
               if(before10mer == eV[length(eV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "e f"] <- "f"
                  wtDF$segment[counter[1:9]] <- test
               }  
               
               lengthcounter <- length(wtDF$segment[counter])
               
               if(after10mer == aV[1]){
                  test <- wtDF$segment[counter]
                  test <- test[(lengthcounter-8):lengthcounter]
                  test[test == "ori a"] <- "ori"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }    
               
               if(after10mer == bV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "a b"] <- "a"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == cV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "b c"] <- "b"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == dV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "c d"] <- "c"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == eV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "d e"] <- "d"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == fV[1]){
                  test <- wtDF$segment[counter]
                  test <- test[(lengthcounter-8):lengthcounter]
                  test[test == "e f"] <- "e"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               wtDF$segmentPos[sort(counter)[pos]] <- 999  
            }else{
               
               wtDF$segment[counter] <- "others"
               wtDF$segmentPos[sort(counter)] <- 1:length(wtDF$segment[counter] )
            }
            
            
         }
         print(i)
         i <- i+1
      }
      
      ## Second time, in case the other_sequence is at the border of a sequence segment
      i <- 1
      while(i < nrow(wtDF)){
         
         while(!is.na(wtDF$segment[i]) & i < nrow(wtDF)) i <- i+1
         
         counter <- list()
         
         while( is.na(wtDF$segment[i]) & ((i+1) <= nrow(wtDF))){
            
            counter[[i]] <- i
            i <- i+1
         }
         
         counter <- unlist(counter)
         #print(counter)
         if(length(counter) > 9 ){
            wtDF$segment[counter] <-  multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)] ##test
            
            if(any(!is.na(multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)]))){
               
               pos <- which(!is.na(multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)]))
               
               before10mer <- wtDF$seq[(sort(counter)[1]-1)]  
               after10mer <- wtDF$seq[(sort(counter)[length(counter)]+1)]  
               
               if(before10mer == oriV[length(oriV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "ori a"] <- "a"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == aV[length(aV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "a b"] <- "b"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == bV[length(bV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "b c"] <- "c"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == cV[length(cV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "c d"] <- "d"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == dV[length(dV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "d e"] <- "e"
                  wtDF$segment[counter[1:9]] <- test
               }   
               
               if(before10mer == eV[length(eV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "e f"] <- "f"
                  wtDF$segment[counter[1:9]] <- test
               }  
               
               lengthcounter <- length(wtDF$segment[counter])
               
               if(after10mer == aV[1]){
                  test <- wtDF$segment[counter]
                  test <- test[(lengthcounter-8):lengthcounter]
                  test[test == "ori a"] <- "ori"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }    
               
               if(after10mer == bV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "a b"] <- "a"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == cV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "b c"] <- "b"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == dV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "c d"] <- "c"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == eV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "d e"] <- "d"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == fV[1]){
                  test <- wtDF$segment[counter]
                  test <- test[(lengthcounter-8):lengthcounter]
                  test[test == "e f"] <- "e"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               wtDF$segmentPos[sort(counter)[pos]] <- 999  
            }else{
               
               wtDF$segment[counter] <- "others"
               wtDF$segmentPos[sort(counter)] <- 1:length(wtDF$segment[counter] )
            }
            
            
         }
         
      }
      
      
      wtDF <- na.omit(wtDF)
      wtDF$pos <- 1:nrow(wtDF)
      
      wtDF$segmentPos[wtDF$segmentPos == 999] <- 0
      
      wtDF$segmentColor[wtDF$segment == "ori"] <- "#66C2A5"
      wtDF$segmentColor[wtDF$segment == "a"] <- "#FC8D62"
      wtDF$segmentColor[wtDF$segment == "b"] <- "#8DA0CB"
      wtDF$segmentColor[wtDF$segment == "c"] <- "#E78AC3"
      wtDF$segmentColor[wtDF$segment == "d"] <- "#A6D854"
      wtDF$segmentColor[wtDF$segment == "e"] <- "#FFD92F"
      wtDF$segmentColor[wtDF$segment == "f"] <- "#E5C494"
      wtDF$segmentColor[wtDF$segment == "others"] <- "#B3B3B3"
      
      
      
      library(ggplot2)
      
      rectangeDF <- list()
      
      i <- 1
      while(i  < nrow(wtDF)){
         
         posi <- wtDF$pos[i]
         lab <- wtDF$segment[i]
         labPos <- wtDF$segmentPos[i]
         
         segmentColor <- wtDF$segmentColor[i]
         
         while(lab == wtDF$segment[i+1] & ((wtDF$segmentPos[i] + wtDF$segmentPos[i+1] == 0 ) | (wtDF$segmentPos[i+1] == wtDF$segmentPos[i]+1 ))   ){
            i <- i+1
            if(i == nrow(wtDF)) break
         }
         
         #posis <- wtDF$segmentPos[labPosStart:labPosEnd]
         
         posi2 <- wtDF$pos[i]
         labPos2 <- wtDF$segmentPos[i]
         mid <- posi+((posi2-posi)/2)
         labPos3 <- labPos2+10-1
         lapPosLength <- posi2-posi+10
         #labPos <- min(posis)
         #labPos3 <- max(posis)+10-1
         labIn <- paste0(lab,"\n", labPos, "-", labPos3)
         labIn3 <- paste0(lab," (", labPos, "-", labPos3,")")
         
         if( wtDF$segmentPos[i] == 0 ) lapPosLength <- posi2-posi+1
         if( wtDF$segmentPos[i] == 0 ) labIn <- paste0(lab,"\n", lapPosLength,"nt")
         if( wtDF$segmentPos[i] == 0 ) labIn3 <- paste0(lab," (", lapPosLength,"nt)")
         
         rectangeDF[[i]] <- data.frame(seg= lab , xmin=posi, xmax=posi2, ymin=4, ymax=6, segmentColor=segmentColor, mid=mid, labIn=labIn, length=lapPosLength, labIn3=labIn3)  
         rectangeDF[[i]]
         if(i != nrow(wtDF)) i <- i + 1
      }
      
      rectangeDF <- rbindlist(rectangeDF)
      
      ## correct insertion length
      rectangeDF$length[rectangeDF$seg == "others"] <- rectangeDF$length[rectangeDF$seg == "others"]+1-10-9
      
      ## correct xmin and xmax by actual sequence length
      rectangeDF$xmin2 <- 0
      rectangeDF$xmax2 <- 0
      
      rectangeDF$xmin2[1] <- 1
      rectangeDF$xmax2[1] <- 1+rectangeDF$length[1]-1
      
      i <- 2
      while(i  <= nrow(rectangeDF)){
         rectangeDF$xmin2[i] <- rectangeDF$xmax2[i-1]+1
         rectangeDF$xmax2[i] <- rectangeDF$xmin2[i] + rectangeDF$length[i]-1
         rectangeDF$mid[i] <- rectangeDF$xmin2[i]+((rectangeDF$xmax2[i]-rectangeDF$xmin2[i])/2)
         i <- i+1
      }
      
      return(rectangeDF)
      
   }
   
   plottingData <- creatPlottingData(inputSequence)
   
   if(all(plottingData$seg == "others")){
      
      message("No JCV NCCR segments found")
      
   }else{
      
      return(plottingData)
      
   }
   
}



## Figure 3

## Figure 3 A
## Serum samples

singleDonors <- c("62", "116", "125", "6" )
singleDonorsSample <- c("JT-46", "JT-61", "JT-17", "JT-6" )

donorList <- list()

##Generate data for plot
for(d in singleDonors){
   
   # d <- 62
   resT <- resAllSave[resAllSave$DonorID == d,]
   resT <- na.omit(resT)
   
   fillup <- unique(resT[,c(6,7)])
   fillup$unassigned <- 100-fillup$sample_totalPercentage
   fillup2 <- resT[c(1:length(unique(resT$sample))),]
   fillup2$sample <- fillup$sample
   fillup2$read_percentage <- fillup$unassigned
   fillup2$globalCluster  <- as.character(fillup2$globalCluster )
   fillup2$globalCluster <- "Others"
   fillup2$DonorTime <- resT$DonorTime[match(fillup2$sample, resT$sample)]
   fillup2$sampleType <- resT$sampleType[match(fillup2$sample, resT$sample)]
   
   resT$globalCluster <- as.character(resT$globalCluster)
   resT <- resT[order(resT$DonorTime,-resT$read_percentage, decreasing = F),]
   resT <- rbind(resT, fillup2)
   
   resT$cluster <- factor(resT$globalCluster)
   
   resT$cluster <- factor(resT$globalCluster, levels = unique(resT$globalCluster[order(resT$read_percentage, decreasing = T)]))
   
   sample_sum_reads <- data.frame( tapply(resT$n_reads, resT$sample, sum) )
   sample_sum_reads$sample <- row.names(sample_sum_reads)
   names(sample_sum_reads) <- c("sum", "sample")
   
   resT$sampleSumHaplotypeReads <-  sample_sum_reads$sum[match(resT$sample, sample_sum_reads$sample)]
   
   #resT$percentage <- round(resT$n_reads/resT$sampleSumHaplotypeReads*100,1)
   
   resT$percentage <- resT$read_percentage
   
   resT$sample <- as.character(resT$sample)
   
   resT$DonorTime_min <- min(resT$DonorTime)
   resT$DonorTime_diff <- resT$DonorTime -   resT$DonorTime_min
   resT$DonorTime_diff <- paste0("Day ", resT$DonorTime_diff)
   
   resT <- resT[resT$sample %in% singleDonorsSample,]
   
   resT$cluster <- as.character(resT$cluster )
   factorLevels <- unique(resT$globalCluster[order(resT$read_percentage, decreasing = T)])
   factorLevels <- factorLevels[factorLevels != "Others"]
   factorLevels <- c(factorLevels, "Others")
   resT$cluster <- factor(resT$globalCluster, levels = factorLevels)
   
   donorList[[d]] <- resT
   
}

library(data.table)
donorTbl <- rbindlist(donorList)

resT <- donorTbl

## Replace internal donor ID
resT$DonorID2 <- testD$newID[match(resT$DonorID, testD$internID)]
resT$DonorID2 <- paste0("P", resT$DonorID2)


resT$DonorID2 <- factor(resT$DonorID2, levels= c(  "P13" ,"P18", "P19", "P2"))

# Sort cluster levels by percentage, ensuring "Others" is always at the bottom
resT$cluster <- factor(
   resT$globalCluster,
   levels = c(
      unique(resT$globalCluster[resT$globalCluster != "Others"][order(-resT$percentage[resT$globalCluster != "Others"])]),
      "Others"
   )
)

# Add cluster length to cluster names
clusters <- levels(resT$cluster)
clustersNamesNew <- sapply(clusters, function(x) {
   if (x != "Others") {
      paste0(x, " (", round(mean(resT$length_consensus_seq[resT$globalCluster == x], na.rm = TRUE)), "nt)")
   } else {
      x
   }
})

resT$sampleType[resT$sampleType == "Liquor"] <- "CSF"
#resT$sampleType[resT$sampleType == "Serum"] <- "Blood"
resT$sampleType <- factor(resT$sampleType, levels= c("CSF", "Serum"))
library(ggplot2)
library(scales)

# Generate the plot
g <- ggplot(data = resT, aes(x = sampleType, y = percentage, fill = cluster)) +
   geom_bar(stat = "identity", position = "fill") + 
   scale_y_continuous(labels = scales::percent_format()) +
   xlab("") +
   ylab("Quasispecies fraction") +
   theme(text = element_text(size = 15))

g <- g + facet_grid(. ~  resT$DonorID2, scales = "free", space = "free") +
   #ggtitle(paste0("Donor ", d)) +
   theme(plot.title = element_text(hjust = 0.5))

g <- g + guides(fill = guide_legend(title = "Cluster ID"))


# Add a color palette for the clusters, setting "Others" to grey
clusters2 <- as.numeric(clusters)
clusters2 <- sort(clusters2)
if(any(clusters == "Archetype")) clusters2 <- c("Archetype", clusters2)
clusters2 <- c(clusters2, "Others")

clustersNamesNew2 <- sapply(clusters2, function(x) {
   if (x != "Others") {
      paste0(x, " (", round(mean(resT$length_consensus_seq[resT$globalCluster == x], na.rm = TRUE)), "nt)")
   } else {
      x
   }
})


g <- g + scale_fill_manual(
   name = "NCCR haplotype ID",
   values = c(hue_pal()(length(clusters) - 1), "grey"),
   breaks =  clusters2 ,
   labels = clustersNamesNew2
)

g <- g+ guides(fill = guide_legend(ncol = 2))
g

#ggsave(g, filename=paste0("Single sample haplotype examplary composition Blood.svg"), width=7.8, height = 4)
ggsave(g, filename=paste0("Figure 3 A.svg"), width=7.8, height = 4)



## Composition of exemplary CSF samples
## Figure 3 B

## CSF first
singleDonors <- c("90", "19", "81", "6" )
singleDonorsSample <- c("JT-60", "D-1617", "JT-58", "JT-7" )

donorList <- list()

##Generate data for plot
for(d in singleDonors){
   
   # d <- 6
   resT <- resAllSave[resAllSave$DonorID == d,]
   resT <- na.omit(resT)
   
   fillup <- unique(resT[,c(6,7)])
   fillup$unassigned <- 100-fillup$sample_totalPercentage
   fillup2 <- resT[c(1:length(unique(resT$sample))),]
   fillup2$sample <- fillup$sample
   fillup2$read_percentage <- fillup$unassigned
   fillup2$globalCluster  <- as.character(fillup2$globalCluster )
   fillup2$globalCluster <- "Others"
   fillup2$DonorTime <- resT$DonorTime[match(fillup2$sample, resT$sample)]
   fillup2$sampleType <- resT$sampleType[match(fillup2$sample, resT$sample)]
   
   resT$globalCluster <- as.character(resT$globalCluster)
   resT <- resT[order(resT$DonorTime,-resT$read_percentage, decreasing = F),]
   resT <- rbind(resT, fillup2)
   
   resT$cluster <- factor(resT$globalCluster)
   
   resT$cluster <- factor(resT$globalCluster, levels = unique(resT$globalCluster[order(resT$read_percentage, decreasing = T)]))
   
   sample_sum_reads <- data.frame( tapply(resT$n_reads, resT$sample, sum) )
   sample_sum_reads$sample <- row.names(sample_sum_reads)
   names(sample_sum_reads) <- c("sum", "sample")
   
   resT$sampleSumHaplotypeReads <-  sample_sum_reads$sum[match(resT$sample, sample_sum_reads$sample)]
   
   #resT$percentage <- round(resT$n_reads/resT$sampleSumHaplotypeReads*100,1)
   
   resT$percentage <- resT$read_percentage
   
   resT$sample <- as.character(resT$sample)
   
   resT$DonorTime_min <- min(resT$DonorTime)
   resT$DonorTime_diff <- resT$DonorTime -   resT$DonorTime_min
   resT$DonorTime_diff <- paste0("Day ", resT$DonorTime_diff)
   
   resT <- resT[resT$sample %in% singleDonorsSample,]
   
   resT$cluster <- as.character(resT$cluster )
   factorLevels <- unique(resT$globalCluster[order(resT$read_percentage, decreasing = T)])
   factorLevels <- factorLevels[factorLevels != "Others"]
   factorLevels <- c(factorLevels, "Others")
   resT$cluster <- factor(resT$globalCluster, levels = factorLevels)
   
   donorList[[d]] <- resT
   
}

library(data.table)
donorTbl <- rbindlist(donorList)

resT <- donorTbl

## Replace internal donor ID
resT$DonorID2 <- testD$newID[match(resT$DonorID, testD$internID)]
resT$DonorID2 <- paste0("P", resT$DonorID2)


resT$DonorID2 <- factor(resT$DonorID2, levels= c(  "P16" ,"P4", "P15", "P2"))


# Sort cluster levels by percentage, ensuring "Others" is always at the bottom
resT$cluster <- factor(
   resT$globalCluster,
   levels = c(
      unique(resT$globalCluster[resT$globalCluster != "Others"][order(-resT$percentage[resT$globalCluster != "Others"])]),
      "Others"
   )
)

# Add cluster length to cluster names
clusters <- levels(resT$cluster)
clustersNamesNew <- sapply(clusters, function(x) {
   if (x != "Others") {
      paste0(x, " (", round(mean(resT$length_consensus_seq[resT$globalCluster == x], na.rm = TRUE)), "nt)")
   } else {
      x
   }
})

resT$sampleType[resT$sampleType == "Liquor"] <- "CSF"
#resT$sampleType[resT$sampleType == "Serum"] <- "Blood"
resT$sampleType <- factor(resT$sampleType, levels= c("CSF", "Serum"))

# Generate the plot
g <- ggplot(data = resT, aes(x = sampleType  , y = percentage, fill = cluster)) +
   geom_bar(stat = "identity", position = "fill") + 
   scale_y_continuous(labels = scales::percent_format()) +
   xlab("") +
   ylab("Quasispecies fraction") +
   theme(text = element_text(size = 15))

g <- g + facet_grid(. ~  resT$DonorID2, scales = "free", space = "free") +
   #ggtitle(paste0("Donor ", d)) +
   theme(plot.title = element_text(hjust = 0.5))

g <- g + guides(fill = guide_legend(title = "Cluster ID"))

# Add a color palette for the clusters, setting "Others" to grey
clusters2 <- as.numeric(clusters)
clusters2 <- sort(clusters2)
if(any(clusters == "Archetype")) clusters2 <- c("Archetype", clusters2)
clusters2 <- c(clusters2, "Others")

clustersNamesNew2 <- sapply(clusters2, function(x) {
   if (x != "Others") {
      paste0(x, " (", round(mean(resT$length_consensus_seq[resT$globalCluster == x], na.rm = TRUE)), "nt)")
   } else {
      x
   }
})

g <- g + scale_fill_manual(
   name = "NCCR haplotype ID",
   values = c(hue_pal()(length(clusters) - 1), "grey"),
   breaks = clusters2,
   labels = clustersNamesNew2
)



g <- g+ guides(fill = guide_legend(ncol = 2))

g

#ggsave(g, filename=paste0("Single sample haplotype examplary composition CSF.png"), width=7.8, height = 4)
ggsave(g, filename=paste0("Figure 3 B.svg"), width=7.8, height = 4)




## Number of distinct haplotyes per sample
## Figure 3 C
clinical <- read.csv("C:/Users/joh/Desktop/HomeOffice/JC/Viruslasten.csv", sep=";")
clinical$Seq.ID <- gsub(" \\(Test\\)","",clinical$Seq.ID)

resAllSave$c <- 1
co <- data.frame(tapply(resAllSave$c, resAllSave$sample, sum))
co$sample <- row.names(co)
names(co) <- c("count", "sample")
co$type <- resAllSave$sampleType[match(co$sample, resAllSave$sample)]
co$c <- 1

co$donor <- resAllSave$DonorID[ match(co$sample, resAllSave$sample)]
co$date <- resAllSave$DonorTime2[ match(co$sample, resAllSave$sample)]

co$date2 <- as.character(lapply(co$date , function(x){
   
   y <- strsplit(x, "\\.")[[1]]
   paste0("20", y[[3]],"-", y[[2]],"-", y[[1]])
   
}))

co$SampleIdentificationID <- paste(co$donor, co$date2, co$type)

co$titer <- clinical$Mittelwert.Viruslast[match(co$sample, clinical$Seq.ID)]
co <- co[co$donor != "23",]

co$titer <- as.numeric(co$titer)


co$earliestDate <- lapply(c(1:nrow(co)), function(x){
   
   curD <- co$date2[co$donor==co$donor[x]]
   
   as.numeric(as.Date(co$date2[x]) - as.Date(min(curD)))
   
   
   
})

co2 <- co[co$earliestDate == 0,]

coco <- data.frame(table( co2$type, co2$count))

coco$Var1 <- as.character(coco$Var1)
coco$Var1[coco$Var1 == "Liquor"] <- "CSF"
#coco$Var1[coco$Var1 == "Serum"] <- "Blood"

g <- ggplot(coco, aes(x = Var2, y = Freq, fill = Var1)) +
   geom_bar(stat = "identity", position = "dodge") +
   labs(x = "Number of detected haplotypes per sample", y = "Number of samples") 

g <- g +
   scale_y_continuous(
      breaks = seq(0, 10, by = 2),  # Set breaks at intervals of 10
      labels = seq(0, 10, by = 2)# Optional: Add custom labels
   )+
   theme(legend.title=element_blank())
g
#ggsave(g, file="Frequency of haplotypes per sample.png", width=5, height = 3)
ggsave(g, file="Figure 3 C.svg", width=5, height = 3)


## Figure 3D 
## Main variant coverage per sample
mvu <- data.frame(tapply(resAllSave$read_percentage, resAllSave$sample, max))

mvu$sample <- row.names(mvu)
mvu$sample_type <- resAllSave$sampleType[match(mvu$sample, resAllSave$sample)]
names(mvu) <- c("coverage", "sample", "sampleType")
mvu$sampleType[mvu$sampleType == "Liquor"] <- "CSF"
#mvu$sampleType[mvu$sampleType == "Serum"] <- "Blood"
mvu <- mvu[mvu$sample %in% co2$sample,]

library(ggplot2)
g <- ggplot(mvu, aes(x=coverage, group=sampleType, fill=sampleType)) +
   geom_density(adjust=1.5, alpha=.4)+xlab("Main variant coverage [%]") +
   theme(legend.title=element_blank())


##ggsave(g, file="Distribution of main variant abundance per sample.png", width=5, height = 3)
ggsave(g, file="Figure 3 D.svg", width=5, height = 3)



## Figure 4

## Figure 4 A
## sample  pairs where quasispecies is similar

# Donor 45
# JT-37 JT-36
# 
# Donor 125
# JT-20 JT-18
# 
# Donor 24
# JT-16 JT-15
# 
# Donor 32
# JT-25 JT-24
# 
# Donor 62 
# JT-44 JT-45 JT-46 JT-47

# singleDonors <- c("24", "32", "45", "62")
# singleSamples <- c("JT-37", "JT-36",  "JT-16", "JT-15", "JT-25", "JT-24", "JT-44", "JT-45" )

singleDonors <- c("24", "32", "45", "62", "125")
singleSamples <- c("JT-37", "JT-36",  "JT-16", "JT-15", "JT-25", "JT-24", "JT-44", "JT-45","JT-46", "JT-47","JT-20", "JT-18" )


donorList <- list()

##Generate data for plot
for(d in singleDonors){
   
   # d <-24
   resT <- resAllSave[resAllSave$DonorID == d,]
   resT <- na.omit(resT)
   
   fillup <- unique(resT[,c(6,7)])
   fillup$unassigned <- 100-fillup$sample_totalPercentage
   fillup2 <- resT[c(1:length(unique(resT$sample))),]
   fillup2$sample <- fillup$sample
   fillup2$read_percentage <- fillup$unassigned
   fillup2$globalCluster  <- as.character(fillup2$globalCluster )
   fillup2$globalCluster <- "Others"
   fillup2$DonorTime <- resT$DonorTime[match(fillup2$sample, resT$sample)]
   fillup2$sampleType <- resT$sampleType[match(fillup2$sample, resT$sample)]
   
   resT$globalCluster <- as.character(resT$globalCluster)
   resT <- resT[order(resT$DonorTime,-resT$read_percentage, decreasing = F),]
   resT <- rbind(resT, fillup2)
   
   resT$cluster <- factor(resT$globalCluster)
   
   resT$cluster <- factor(resT$globalCluster, levels = unique(resT$globalCluster[order(resT$read_percentage, decreasing = T)]))
   
   sample_sum_reads <- data.frame( tapply(resT$n_reads, resT$sample, sum) )
   sample_sum_reads$sample <- row.names(sample_sum_reads)
   names(sample_sum_reads) <- c("sum", "sample")
   
   resT$sampleSumHaplotypeReads <-  sample_sum_reads$sum[match(resT$sample, sample_sum_reads$sample)]
   
   #resT$percentage <- round(resT$n_reads/resT$sampleSumHaplotypeReads*100,1)
   
   resT$percentage <- resT$read_percentage
   
   resT$sample <- as.character(resT$sample)
   
   resT$DonorTime_min <- min(resT$DonorTime)
   resT$DonorTime_diff <- resT$DonorTime -   resT$DonorTime_min
   resT$DonorTime_diff <- paste0("Day ", resT$DonorTime_diff)
   
   resT <- resT[resT$sample %in% singleSamples,]
   
   resT$cluster <- as.character(resT$cluster )
   factorLevels <- unique(resT$globalCluster[order(resT$read_percentage, decreasing = T)])
   factorLevels <- factorLevels[factorLevels != "Others"]
   factorLevels <- c(factorLevels, "Others")
   resT$cluster <- factor(resT$globalCluster, levels = factorLevels)
   
   donorList[[d]] <- resT
   
}

library(data.table)
donorTbl <- rbindlist(donorList)



resT <- donorTbl

resT$DonorID2 <- testD$newID[match(resT$DonorID, testD$internID)]
resT$DonorID2 <- paste0("P", resT$DonorID2)

#resT$DonorID2 <- factor(resT$DonorID2, levels= c(  "D90" ,"D42", "D105", "D19", "D142",  "D133" ))



# Sort cluster levels by percentage, ensuring "Others" is always at the bottom
resT$cluster <- factor(
   resT$globalCluster,
   levels = c(
      unique(resT$globalCluster[resT$globalCluster != "Others"][order(-resT$percentage[resT$globalCluster != "Others"])]),
      "Others"
   )
)

# Add cluster length to cluster names
clusters <- levels(resT$cluster)
clustersNamesNew <- sapply(clusters, function(x) {
   if (x != "Others") {
      paste0(x, " (", round(mean(resT$length_consensus_seq[resT$globalCluster == x], na.rm = TRUE)), "nt)")
   } else {
      x
   }
})

resT$sampleType[resT$sampleType == "Liquor"] <- "CSF"
#resT$sampleType[resT$sampleType == "Serum"] <- "Blood"
resT$sampleType <- factor(resT$sampleType, levels= c("CSF","Serum"))

resT$DonorID3 <- resT$DonorID2

resT$DonorID3[resT$DonorID3 == "P19"] <- "P19\npair 1"
resT$DonorID3[resT$sample %in% c("JT-44","JT-45")] <- "P13\npair 1"
resT$DonorID3[resT$sample %in% c("JT-46","JT-47")] <- "P13\npair 2"


# Donor 62 
# JT-44 JT-45 JT-46 JT-47

# Generate the plot
g <- ggplot(data = resT, aes(x = sampleType, y = percentage, fill = cluster)) +
   geom_bar(stat = "identity", position = "fill") + 
   scale_y_continuous(labels = scales::percent_format()) +
   xlab("") +
   ylab("Quasispecies fraction") +
   theme(text = element_text(size = 15))

# Facet by DonorID, grouping CSF and Blood beside each other within each Donor
g <- g + facet_grid(. ~ DonorID3, scales = "free_x", space = "free_x") +
   theme(
      plot.title = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "grey85"),
      strip.text.x = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1, "lines")#,
      # axis.text.x = element_text(angle = 45, hjust = 1)
   )

# Add a color palette for the clusters, setting "Others" to grey
clusters2 <- as.numeric(clusters)
clusters2 <- sort(clusters2)
if(any(clusters == "Archetype")) clusters2 <- c("Archetype", clusters2)
clusters2 <- c(clusters2, "Others")

clustersNamesNew2 <- sapply(clusters2, function(x) {
   if (x != "Others") {
      paste0(x, " (", round(mean(resT$length_consensus_seq[resT$globalCluster == x], na.rm = TRUE)), "nt)")
   } else {
      x
   }
})

g <- g + scale_fill_manual(
   name = "NCCR haplotype ID",
   values = c(hue_pal()(length(clusters) - 1), "grey"),
   breaks = clusters2,
   labels = clustersNamesNew2
)


g

##ggsave(g, filename=paste0("Samples with CSF and Blood haplotype composition no change.png"),width=13, height = 4.6)
ggsave(g, filename=paste0("Figure 4 A.svg"),width=13, height = 4.6)




## Figure 4 B

singleDonors <- c("125", "6","13", "81", "116")
singleSamples <- c("JT-17", "JT-19", "JT-5", "JT-4", "JT-7", "JT-6", "JT-10", "JT-11" ,"JT-52", "JT-51", "JT-55", "JT-61")


donorList <- list()

##Generate data for plot
for(d in singleDonors){
   
   # d <- 6
   resT <- resAllSave[resAllSave$DonorID == d ,]
   resT <- na.omit(resT)
   
   fillup <- unique(resT[,c(6,7)])
   fillup$unassigned <- 100-fillup$sample_totalPercentage
   fillup2 <- resT[c(1:length(unique(resT$sample))),]
   fillup2$sample <- fillup$sample
   fillup2$read_percentage <- fillup$unassigned
   fillup2$globalCluster  <- as.character(fillup2$globalCluster )
   fillup2$globalCluster <- "Others"
   fillup2$DonorTime <- resT$DonorTime[match(fillup2$sample, resT$sample)]
   fillup2$sampleType <- resT$sampleType[match(fillup2$sample, resT$sample)]
   
   # singleSamples <- c("JT-17", "JT-19", "JT-5", "JT-4", "JT-52", "JT-51", "JT-55", "JT-61")
   resT$globalCluster <- as.character(resT$globalCluster)
   resT <- resT[order(resT$DonorTime,-resT$read_percentage, decreasing = F),]
   resT <- rbind(resT, fillup2)
   
   resT$cluster <- factor(resT$globalCluster)
   
   resT$cluster <- factor(resT$globalCluster, levels = unique(resT$globalCluster[order(resT$read_percentage, decreasing = T)]))
   
   sample_sum_reads <- data.frame( tapply(resT$n_reads, resT$sample, sum) )
   sample_sum_reads$sample <- row.names(sample_sum_reads)
   names(sample_sum_reads) <- c("sum", "sample")
   
   resT$sampleSumHaplotypeReads <-  sample_sum_reads$sum[match(resT$sample, sample_sum_reads$sample)]
   
   #resT$percentage <- round(resT$n_reads/resT$sampleSumHaplotypeReads*100,1)
   
   resT$percentage <- resT$read_percentage
   
   resT$sample <- as.character(resT$sample)
   
   resT$DonorTime_min <- min(resT$DonorTime)
   resT$DonorTime_diff <- resT$DonorTime -   resT$DonorTime_min
   resT$DonorTime_diff <- paste0("Day ", resT$DonorTime_diff)
   
   resT <- resT[resT$sample %in% singleSamples,]
   
   resT$cluster <- as.character(resT$cluster )
   factorLevels <- unique(resT$globalCluster[order(resT$read_percentage, decreasing = T)])
   factorLevels <- factorLevels[factorLevels != "Others"]
   factorLevels <- c(factorLevels, "Others")
   resT$cluster <- factor(resT$globalCluster, levels = factorLevels)
   
   donorList[[d]] <- resT
   
}

library(data.table)
donorTbl <- rbindlist(donorList)



resT <- donorTbl


resT$DonorID2 <- testD$newID[match(resT$DonorID, testD$internID)]
resT$DonorID2 <- paste0("P", resT$DonorID2)


# Sort cluster levels by percentage, ensuring "Others" is always at the bottom
resT$cluster <- factor(
   resT$globalCluster,
   levels = c(
      unique(resT$globalCluster[resT$globalCluster != "Others"][order(-resT$percentage[resT$globalCluster != "Others"])]),
      "Others"
   )
)

# Add cluster length to cluster names
clusters <- levels(resT$cluster)
clustersNamesNew <- sapply(clusters, function(x) {
   if (x != "Others") {
      paste0(x, " (", round(mean(resT$length_consensus_seq[resT$globalCluster == x], na.rm = TRUE)), "nt)")
   } else {
      x
   }
})

resT$sampleType[resT$sampleType == "Liquor"] <- "CSF"
#resT$sampleType[resT$sampleType == "Serum"] <- "Blood"
resT$sampleType <- factor(resT$sampleType, levels= c("CSF","Serum"))

resT$DonorID3 <- resT$DonorID2

resT$DonorID3[resT$sample %in% c("JT-5", "JT-4")] <- "P2\npair 1" 
resT$DonorID3[resT$sample %in% c("JT-6", "JT-7")] <- "P2\npair 2" 
resT$DonorID3[resT$DonorID3 == "P19"] <- "P19\npair 2"

# Generate the plot
g <- ggplot(data = resT, aes(x =  sampleType, y = percentage, fill = cluster)) +
   geom_bar(stat = "identity", position = "fill") + 
   scale_y_continuous(labels = scales::percent_format()) +
   xlab("") +
   ylab("Quasispecies fraction") +
   theme(text = element_text(size = 15))

# Facet by DonorID, grouping CSF and Blood beside each other within each Donor
g <- g + facet_grid(. ~ DonorID3, scales = "free_x", space = "free_x") +
   theme(
      plot.title = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "grey85"),
      strip.text.x = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1, "lines")#,
      # axis.text.x = element_text(angle = 45, hjust = 1)
   )

# Add a color palette for the clusters, setting "Others" to grey
clusters2 <- as.numeric(clusters)
clusters2 <- sort(clusters2)
if(any(clusters == "Archetype")) clusters2 <- c("Archetype", clusters2)
clusters2 <- c(clusters2, "Others")

clustersNamesNew2 <- sapply(clusters2, function(x) {
   if (x != "Others") {
      paste0(x, " (", round(mean(resT$length_consensus_seq[resT$globalCluster == x], na.rm = TRUE)), "nt)")
   } else {
      x
   }
})

library(scales)
g <- g + scale_fill_manual(
   name = "NCCR haplotype ID",
   values = c(hue_pal()(length(clusters) - 1), "grey"),
   breaks = clusters2,
   labels = clustersNamesNew2
)

g

##ggsave(g, filename=paste0("Samples with CSF and Serum haplotype composition.png"), width=14, height = 4.6)
ggsave(g, filename=paste0("Figure 4 B.svg"), width=14, height = 4.6)


## Figure 5

## Figure 5 A
library(VarCon)
library(DECIPHER)
library(ModCon)
library(ape)
library(MASS)
library(ggtreeExtra)
library(gridExtra)
library(ggplot2)
library(ape)
library(msa)
library(ggtree)

## load tree
aln <- prepareReferenceFasta("Haplotye_consensus.aln.fasta")
aln <- aln[names(aln) != "artificialARC"]


aln_names <- names(aln)
aln_names <- as.character(lapply(aln_names, function(x){ 
   if(x != "ARC") {
      strsplit(x, "_")[[1]][[2]]
   }else{ "ARC"}
}))

names(aln) <- aln_names
aln <- aln[names(aln) != "26"]


ttdist <- DistanceMatrix(aln, includeTerminalGaps=T, type ="dist", penalizeGapLetterMatches=T)
my_nj <- ape::nj(ttdist)

my_nj <- root(my_nj, "ARC", resolve.root = TRUE)

p <- ggtree(my_nj, options(ignore.negative.edge=TRUE) ) +
   geom_tiplab(size=3, align=TRUE, linesize=.5)+ggtitle("original")
p

haplotype_donor <- data.frame(haplotype= sort(unique(resAllSave$globalCluster)))

resAllSave$DonorID2 <- testD$newID[match(resAllSave$DonorID, testD$internID)]
haplotype_donor$DonorID <- as.character(lapply(haplotype_donor$haplotype, function(x){
   y <- unique(resAllSave$DonorID2[resAllSave$globalCluster == x])
   y <- y[! y  %in% c("23")]
   paste(y, collapse=",")
   
}))

haplotype_donor <- haplotype_donor[!haplotype_donor$haplotype %in% c("108","113","117","26","27"),]



haplotype_donor$haplotype[haplotype_donor$haplotype == "Archetype"] <- "ARC"
haplotype_donor2 <- haplotype_donor
row.names(haplotype_donor2) <- haplotype_donor2$haplotype
haplotype_donor2$haplotype <- NULL
haplotype_donor2$DonorID[haplotype_donor2$DonorID == "NA,3"] <- "3"

library(ape)

write.csv2(haplotype_donor2, file="annotation_new.csv", row.names = T)

library(metafolio)
n = 22
cols = gg_color_hue(n)

haplotype_donor3 <- data.frame(donor=unique(haplotype_donor2$DonorID),type="label", color=cols, font= "bold 2")

haplotype_donor4 <-  haplotype_donor2
haplotype_donor4$color <- haplotype_donor3$color[match(haplotype_donor4$DonorID,haplotype_donor3$donor )]
haplotype_donor4$font <-  "bold 2"
haplotype_donor4$DonorID <- NULL
haplotype_donor4$lab <- "label"
haplotype_donor4 <- haplotype_donor4[,c(1,4,2,3)]

write.csv2(haplotype_donor4, file="annotationColor_new.csv", row.names = T)


## Do Tree creation in iTOL with "MyTree" which will be saved in the following



## Figure 5 B
#### Segment composition heatmap
library(VarCon)
haplotypes <- prepareReferenceFasta("C:\\Users\\joh\\Desktop\\HomeOffice\\JC\\Haplotye_consensus.fasta")

haplotypes <- data.frame(name= names(haplotypes), sequence = haplotypes)

haplotypes$length <- nchar(haplotypes$sequence)
haplotypes <- haplotypes[order(haplotypes$length, decreasing = T),]
haplotypes$length <- NULL


## add WT to haplotype dataframe
wt <- data.frame(name="ARC", sequence= "GCCTCGGCCTCCTGTATATATAAAAAAAAGGGAAGGTAGGGAGGAGCTGGCTAAAACTGGATGGCTGCCAGCCAAGCATGAGCTCATACCTAGGGAGCCAACCAGCTGACAGCCAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCCCTCCCATAGTCCTTAATCACAAGTAAACAAAGCACAAGGGGAAGTGGAAAGCAGCCAAGGGAACATGTTTTGCGAGCCAGAGCTGTTTTGGCTTGTCACCAGCTGGCC")

haplotypes <- rbind(wt, haplotypes)
haplotypes$cl <- lapply(haplotypes$name, function(x){
   
   if(x != "ARC"){
      strsplit(x, "_")[[1]][[2]]
   } else{ "ARC"}
})


haplotypes$cl[haplotypes$name == "ARC"] <- "ARC"
haplotypes <- haplotypes[haplotypes$cl != "Archetype",]




## Continue analysis

resAll <- list()
yStart <- 100

library(data.table)
for(i in 1:nrow(haplotypes)){
   
   res <- get_table(haplotypes$sequence[i])
   
   res$ymin <- yStart-0.5
   res$ymax <- yStart+0.5
   res$ymid <- yStart
   
   res$allMid <- max(res$xmax)/2
   res$allLab <- haplotypes$name[i]
   
   yStart <- yStart-1
   
   resAll[[i]] <- res
}

resAll <- rbindlist(resAll)

resAll$labIn[resAll$seg %in% "ori"] <- ""
resAll$labIn[resAll$seg %in% "others"] <- ""

resAll$labIn3[resAll$seg %in% "ori"] <- ""
resAll$labIn3[resAll$seg %in% "others"] <- ""


resAll$ylabminus <- -20

segs <- c("ori", letters[1:6])
segs2 <- paste0(segs, " (", c(12,25,23,55,66,18,70), "nt)")
segs2 <- c(segs2, "others")
segs <- c(segs, "others")


resAll$hapLab <- as.character(lapply(1:nrow(resAll), function(x){
   
   fileName <- resAll$allLab[x]
   
   if(fileName != "ARC"){
      
      fileName <- strsplit(fileName, "_")[[1]]
      
      paste0(fileName[[2]], " (maxFreq ", fileName[[4]], ")")
   }else{
      "ARC"
   }
}))

resAll$hapLabSmall <- lapply(resAll$hapLab, function(x){
   
   strsplit(x, " ")[[1]][[1]]
})

resAll <- resAll[!resAll$hapLabSmall %in% c("108","113","117","26","27"),]


## Plot heatmap phylotree with adjusted segmentlength coloring
resAll$hapLabSmall <- as.character( resAll$hapLabSmall)

my_nj2 <- keep.tip(my_nj, unique(resAll$hapLabSmall))

write.tree(my_nj2, file="Mytree.nwk")
p <- ggtree(my_nj2)  + geom_tiplab(size=4, align=TRUE, linesize=.5)

tree2 <- ladderize(my_nj2, right = FALSE)

is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
ordered_tips <- tree2$edge[is_tip, 2]
clusterOrder <- my_nj2$tip.label[ordered_tips]
clusterOrder2 <- rev(clusterOrder)

resT2 <- data.frame(globalCluster = clusterOrder)

resT2$globalClusterFactor <- factor(resT2$globalCluster, levels=clusterOrder)

##Create data
resT2$label <-  resT2$globalClusterFactor 
resT2 <- data.frame(label= unique(resT2$label))

resT2$a <- as.numeric(lapply(resT2$label, function(haplo){
   
   curTab <- resAll[resAll$hapLabSmall == haplo,]
   curTab <- curTab[curTab$seg == "a",]
   sum(curTab$length)
   
}))

resT2$b <- as.numeric(lapply(resT2$label, function(haplo){
   
   curTab <- resAll[resAll$hapLabSmall == haplo,]
   curTab <- curTab[curTab$seg == "b",]
   sum(curTab$length)
   
}))

resT2$c <- as.numeric(lapply(resT2$label, function(haplo){
   
   curTab <- resAll[resAll$hapLabSmall == haplo,]
   curTab <- curTab[curTab$seg == "c",]
   sum(curTab$length)
   
}))

resT2$d <- as.numeric(lapply(resT2$label, function(haplo){
   
   curTab <- resAll[resAll$hapLabSmall == haplo,]
   curTab <- curTab[curTab$seg == "d",]
   sum(curTab$length)
   
}))

resT2$e <- as.numeric(lapply(resT2$label, function(haplo){
   
   curTab <- resAll[resAll$hapLabSmall == haplo,]
   curTab <- curTab[curTab$seg == "e",]
   sum(curTab$length)
   
}))

resT2$f <- as.numeric(lapply(resT2$label, function(haplo){
   
   curTab <- resAll[resAll$hapLabSmall == haplo,]
   curTab <- curTab[curTab$seg == "f",]
   sum(curTab$length)
   
}))

resT2$norm_a <- round(resT2$a/25*100)
resT2$norm_b <- round(resT2$b/23*100)
resT2$norm_c <- round(resT2$c/55*100)
resT2$norm_d <- round(resT2$d/66*100)
resT2$norm_e <- round(resT2$e/18*100)
resT2$norm_f <- round(resT2$f/70*100)

resT2$a <- NULL
resT2$b<- NULL
resT2$c<- NULL
resT2$d<- NULL
resT2$e<- NULL
resT2$f<- NULL

row.names(resT2) <- resT2$label

resT2$label <- NULL

names(resT2) <- c("a", "b", "c", "d", "e", "f")

resTsave <- resT2

resT2[resT2 == 0] <- NA


names(resT2) <- toupper(names(resT2))

pp <-    gheatmap(p, resT2 ,offset=0.03, width=0.5,
                  #colnames_angle=45, colnames_offset_y = 0.25,colnames_offset_x =0.001, colnames=TRUE,
                  colnames_angle=0, colnames_offset_y = .5,colnames_offset_x =0, colnames=TRUE,
                  colnames_position='top',font.size = 2,
                  legend_title = "Number of block\noccurrences")+
   scale_fill_gradientn(
      colors = c("black","darkblue", "lightblue", "gray90", "yellow", "red"),
      values = scales::rescale(c(0, 0.00001,99.999, 100, 101, 350)), # Breakpoints rescaled between 0 and 1
      limits = c(0, 350),  # Ensure scale limits match data range
      na.value="black"
   )




pp <- pp + labs(fill="Difference in segment content\nrelative to ARC [%]")
pp <- pp+ theme(
   legend.position = "bottom",
   legend.direction = "horizontal",
   legend.text = element_text(angle = 0)  # ensure horizontal labels
)
pp

##ggsave(pp, file=paste("Haplotype segment composition heatmap small offset.svg"), width=12, height = 14)
ggsave(pp, file=paste("Figure 5 B.svg"), width=12, height = 14)



### Prepare sequences for Figure 6
## Get first disrupted d segment for alignment

## Redefine getTable function to extract first d segment

get_DSeq <- function(inputSequence){
   
   #inputSequence <- haplotypes$sequence[haplotypes$name ==  "Haplotype_3_maxFreq_2.2%"]
   #inputSequence <- readClipboard()
   ori <- "GCCTCGGCCTCC"
   a <- "CCTGTATATATAAAAAAAAGGGAAG"
   b <- "GTAGGGAGGAGCTGGCTAAAACT"
   c <- "GGATGGCTGCCAGCCAAGCATGAGCTCATACCTAGGGAGCCAACCAGCTGACAGC"
   d <- "CAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCCCTCCCATAGTCCTTAATCAC"
   e <- "AAGTAAACAAAGCACAAG"
   f <- "GGGAAGTGGAAAGCAGCCAAGGGAACATGTTTTGCGAGCCAGAGCTGTTTTGGCTTGTCACCAGCTGGCC"
   
   getOverlappingVectorsFromVector2 <- function(largeVector, subvectorLength, subvectorOverlap) 
   {
      startPositions = seq(1, length(largeVector), by = subvectorLength - 
                              subvectorOverlap)
      end_positions = startPositions + subvectorLength - 1
      end_positions[end_positions > length(largeVector)] = length(largeVector)
      lapply(seq_len(length(startPositions)), function(x) largeVector[startPositions[x]:end_positions[x]])
   }
   
   
   
   creatSubreadVec <- function(vec, len=10){
      
      vec <- strsplit(vec, "")[[1]]
      vec <- getOverlappingVectorsFromVector2(vec, 10, 9)
      vec <- as.character(lapply(vec, function(x){paste0(x, collapse="")} ))
      vec <- vec[nchar(vec) == 10]
      return(vec)
   }
   
   oriV <- creatSubreadVec(ori)
   aV <- creatSubreadVec(a)
   bV <- creatSubreadVec(b)
   cV <- creatSubreadVec(c)
   dV <- creatSubreadVec(d)
   eV <- creatSubreadVec(e)
   fV <- creatSubreadVec(f)
   
   nameVec <- c(rep("ori", length(oriV)),
                rep("a", length(aV)),
                rep("b", length(bV)),
                rep("c", length(cV)),
                rep("d", length(dV)),
                rep("e", length(eV)),
                rep("f", length(fV)))
   
   nameVecPos <- c(seq(1, length(oriV), 1),
                   seq(1, length(aV), 1),
                   seq(1, length(bV), 1),
                   seq(1, length(cV), 1),
                   seq(1, length(dV), 1),
                   seq(1, length(eV), 1),
                   seq(1, length(fV), 1))
   
   multimer <- data.frame(name= nameVec, pos=nameVecPos, 
                          multimer=c(oriV, aV, bV, cV, dV, eV, fV))
   
   
   ori2 <- substr(ori,nchar(ori)-8 ,nchar(ori))
   a1 <-  substr(a, 1, 9)
   b1 <-  substr(b, 1, 9)
   c1 <-  substr(c, 1, 9)
   d1 <-  substr(d, 1, 9)
   e1 <-  substr(e, 1, 9)
   f1 <-  substr(f, 1, 9)
   
   a2 <-  substr(a,nchar(a)-8 ,nchar(a))
   b2 <-  substr(b,nchar(b)-8 ,nchar(b))
   c2 <-  substr(c,nchar(c)-8 ,nchar(c))
   d2 <-  substr(d,nchar(d)-8 ,nchar(d))
   e2 <-  substr(e,nchar(e)-8 ,nchar(e))
   
   oa <- paste0(ori2,a1)
   ab <- paste0(a2,b1)
   bc <- paste0(b2,c1)
   cd <- paste0(c2,d1)
   de <- paste0(d2,e1)
   ef <- paste0(e2,f1)
   
   oaV <- creatSubreadVec(oa)
   abV <- creatSubreadVec(ab)
   bcV <- creatSubreadVec(bc)
   cdV <- creatSubreadVec(cd)
   deV <- creatSubreadVec(de)
   efV <- creatSubreadVec(ef)
   
   
   nameVec <- c(rep("ori a", length(oaV)),
                rep("a b", length(abV)),
                rep("b c", length(bcV)),
                rep("c d", length(cdV)),
                rep("d e", length(deV)),
                rep("e f", length(efV)))
   
   
   multimerBorders  <- data.frame(name= nameVec, pos= NA, 
                                  multimer=c(oaV, abV, bcV, cdV, deV, efV))
   
   creatPlottingData <- function(promotorSeq){
      
      # promotorSeq <- inputSequence
      
      if(nchar(promotorSeq)<20 )
         stop("ERROR during setting of variable 'promotorSeq'. The entered sequence", 
              " must be at least 20nt long.")  
      if (!all(strsplit(promotorSeq, "")[[1]] %in% c("a", "c", "g", "t", 
                                                     "G", "C", "T", "A"))) 
         stop("ERROR during setting of variable 'promotorSeq'. The entered sequence", 
              " must be a character string of A C G and T.")    
      
      promotorSeq <- toupper(promotorSeq)
      wtV <- creatSubreadVec(promotorSeq)
      
      wtDF <- data.frame(pos= 0, seq=wtV)
      
      wtDF$segment <- multimer$name[match(wtDF$seq, multimer$multimer)]
      wtDF$segmentPos <- multimer$pos[match(wtDF$seq, multimer$multimer)]
      
      wtDF$pos <- 1:nrow(wtDF)
      
      i <- 1
      while(i < nrow(wtDF)){
         
         while(!is.na(wtDF$segment[i]) & i < nrow(wtDF)) i <- i+1
         
         counter <- list()
         
         while( is.na(wtDF$segment[i]) & ((i+1) <= nrow(wtDF))){
            
            counter[[i]] <- i
            i <- i+1
         }
         
         counter <- unlist(counter)
         
         if(length(counter) > 9 ){
            wtDF$segment[counter] <-  multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)] ##test
            
            if(any(!is.na(multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)]))){
               
               pos <- which(!is.na(multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)]))
               
               before10mer <- wtDF$seq[(sort(counter)[1]-1)] 
               after10mer <- wtDF$seq[(sort(counter)[length(counter)]+1)]  
               
               if(before10mer == oriV[length(oriV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "ori a"] <- "a"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == aV[length(aV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "a b"] <- "b"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == bV[length(bV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "b c"] <- "c"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               # if(before10mer == cV[length(cV)]){
               #    test <- wtDF$segment[counter]
               #    test <- test[1:9]
               #    test[test == "c d"] <- "d"
               #    wtDF$segment[counter[1:9]] <- test
               # } 
               
               if(before10mer == dV[length(dV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "d e"] <- "e"
                  wtDF$segment[counter[1:9]] <- test
               }   
               
               if(before10mer == eV[length(eV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "e f"] <- "f"
                  wtDF$segment[counter[1:9]] <- test
               }  
               
               lengthcounter <- length(wtDF$segment[counter])
               
               if(after10mer == aV[1]){
                  test <- wtDF$segment[counter]
                  test <- test[(lengthcounter-8):lengthcounter]
                  test[test == "ori a"] <- "ori"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }    
               
               if(after10mer == bV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "a b"] <- "a"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == cV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "b c"] <- "b"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == dV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "c d"] <- "c"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == eV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "d e"] <- "d"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == fV[1]){
                  test <- wtDF$segment[counter]
                  test <- test[(lengthcounter-8):lengthcounter]
                  test[test == "e f"] <- "e"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               wtDF$segmentPos[sort(counter)[pos]] <- 999  
            }else{
               
               wtDF$segment[counter] <- "others"
               wtDF$segmentPos[sort(counter)] <- 1:length(wtDF$segment[counter] )
            }
            
            
         }
         #print(i)
         i <- i+1
      }
      
      
      wtDF <- na.omit(wtDF)
      wtDF$row <- 1:nrow(wtDF)
      
      
      dStart <- 1
      while(wtDF$segment[dStart] != "d"){
         
         if(dStart == nrow(wtDF)) break
         
         dStart <- dStart +1
      }  
      
      dEnd <- dStart
      while(wtDF$segment[dEnd] == "d"){
         
         if(wtDF$segment[dEnd+1] != "d") break
         dEnd <- dEnd +1
         
      }      
      
      if(dStart != nrow(wtDF)){
         
         return(substr(inputSequence, wtDF$pos[dStart] , wtDF$pos[dEnd]+9))
      }else{
         return("")
      }
      
   }
   
   dSegSeq <- creatPlottingData(inputSequence)
   
   return(dSegSeq)
   
}


dSegList <- list()

for(i in 1:nrow(haplotypes)){
   
   print(i)
   res <- get_DSeq(haplotypes$sequence[i])
   
   dSegList[[i]] <- res
}

dSegList <- unlist(dSegList)

sSegDF <- data.frame(name= unlist(haplotypes$cl), dSeg=dSegList )


sSegDF <- sSegDF[!sSegDF$name %in% c("108", "113", "117", "26", "27"),]


sSegDF2 <- sSegDF[sSegDF$dSeg != "",]


i <-1
sSegDF2$alnSeq <- ""

for(i in 1:nrow(sSegDF2)){
   
   curDSeq <- sSegDF2$dSeg[i]
   names(curDSeq) <- sSegDF2$name[i]
   Dseq <- "CAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCCCTCCCATAGTCCTTAATCACAAGTAAACAAAGCACAAG"
   names(Dseq) <- "ARC_dSeg"
   
   ttDNAstringset <- DNAStringSet(c(curDSeq,Dseq ))
   ttDNAstringsetaln <- msa(ttDNAstringset, method="Muscle")
   ttDNAstringsetalnD <- DNAStringSet(ttDNAstringsetaln)
   ttDNAstringsetalnD <- as.character(ttDNAstringsetalnD)
   ttDNAstringsetalnD <- ttDNAstringsetalnD[names(ttDNAstringsetalnD) ==  sSegDF2$name[i]]
   ttDNAstringsetalnD <- substr(ttDNAstringsetalnD, 1, 66)
   
   sSegDF2$alnSeq[i] <- ttDNAstringsetalnD
   
}


library(seqinr)


sequences <- sSegDF2$alnSeq
namesIn <- sSegDF2$name

for(i in seq_len(length(sequences))){
   
   write.fasta(sequences=sequences[i], names=namesIn[i], 
               as.string=FALSE, file.out="Haplotypes_DSegment_aln2.fasta",
               open = "a",nbchar = 60)
   
}



### Figure 7

h2 <- haplotypes[as.character(haplotypes$cl) %in% row.names(haplotype_donor4),]

hall <- data.frame(ID=unlist(h2$cl), sequence=h2$sequence, repeatEncoding="", structuralClassification="")

i <- "ARC"
##get repetition coding
curRes <- resAll[resAll$hapLabSmall == i]
curRes <- curRes[curRes$seg != "ori",]

curCode <- curRes$seg
curCode <- toupper(curCode)
curCodeDF <- data.frame(seg=curCode, len=curRes$length)
curCodeDF_arc <- curCodeDF

find_largest_repeating_substring <- function(s) {
   n <- nchar(s)
   max_substr <- ""
   
   # Loop through possible substrings
   for (len in n:1) {  # Start from the longest possible length
      for (iiii in 1:(n - len + 1)) {
         substring_candidate <- substr(s, iiii, iiii + len - 1)
         
         # Check if the substring repeats elsewhere in the string
         if (length(gregexpr(substring_candidate, s)[[1]]) > 1) {
            if (nchar(substring_candidate) > nchar(max_substr)) {
               max_substr <- substring_candidate
            }
         }
      }
      # If we've found a repeating substring of this length, stop
      if (nchar(max_substr) > 0) {
         break
      }
   }
   return(max_substr)
}

hall$largestRepeatSeq <- ""
hall$largestRepeatSeqCount <- 0
hall$structuralClassification <- ""
hall$structuralClassification2 <- ""

hall <- hall[hall$ID != "27",]


for(i in hall$ID){
   
   # i <- "59"
   
   ##get repetition coding
   curRes <- resAll[resAll$hapLabSmall == i]
   curRes <- curRes[curRes$seg != "ori",]
   
   curCode <- curRes$seg
   curCode <- toupper(curCode)
   curCodeDF <- data.frame(seg=curCode, len=curRes$length)
   curCodeDF$arcLen <- curCodeDF_arc$len[match(curCodeDF$seg, curCodeDF_arc$seg)]
   curCodeDF$del <- 0
   curCodeDF$del[curCodeDF$len != curCodeDF$arcLen] <- 1
   
   curCodeDF$seg2 <- ""
   curCodeDF$seg2 <- as.character(lapply(1:nrow(curCodeDF), function(x){
      
      if(curCodeDF$del[x] == 1){
         curCodeDF$seg2[x] <- tolower(curCodeDF$seg[x])
      } else{
         curCodeDF$seg2[x] <- curCodeDF$seg[x]
      }
      
   }))
   
   curCodeDF <- curCodeDF[curCodeDF$seg != "OTHERS",]
   
   ## To Be corrected!
   curCodeDF <- curCodeDF[curCodeDF$seg != "E F",]
   
   y <- paste0(curCodeDF$seg2, collapse = "")
   
   x <-  paste0(curCodeDF$seg, collapse = "")
   x <- gsub("AAA", "A", x)
   x <- gsub("BBB", "B", x)
   x <- gsub("CCC", "C", x)
   x <- gsub("DDD", "D", x)
   x <- gsub("EEE", "E", x)
   x <- gsub("FFF", "F", x)
   x <- gsub("AA", "A", x)
   x <- gsub("BB", "B", x)
   x <- gsub("CC", "C", x)
   x <- gsub("DD", "D", x)
   x <- gsub("EE", "E", x)
   x <- gsub("FF", "F", x)
   rm(countRep)
   rm(largestRepetition)
   
   largestRepetition <- find_largest_repeating_substring(x)
   hall$repeatEncoding[hall$ID == i] <- y
   
   if(largestRepetition !=""){ 
      countRep <- unlist(gregexpr2(largestRepetition, x))
      countRep <- countRep[countRep != -1]
      hall$largestRepeatSeq[hall$ID == i] <- largestRepetition
      hall$largestRepeatSeqCount[hall$ID == i] <- length(countRep)
   }
   
   x <- strsplit(x, "")[[1]]
   hall$structuralClassification[hall$ID == i] <- "Type I"
   if(any(x %in% c("B", "D"))) hall$structuralClassification[hall$ID == i] <- "Type II"
   
   hall$structuralClassification2[hall$ID == i] <- "S"
   if(any(duplicated(x))) hall$structuralClassification2[hall$ID == i] <- "R"
   
}

hall$structuralClassificationCombined <- paste0(hall$structuralClassification, "-", hall$structuralClassification2)

hall$donor <- as.character(lapply(hall$ID, function(x){
   
   curDo <- resAllSave$DonorID[resAllSave$globalCluster == x]
   paste(unique(curDo), collapse=",")
   
}))

hall$donor[hall$ID == "ARC"] <- "13"

hall$sampleType <- as.character(lapply(hall$ID, function(x){
   
   curDo <- resAllSave$sampleType[resAllSave$globalCluster == x]
   paste(unique(curDo), collapse=",")
   
}))

hall$sampleType[hall$ID == "ARC"] <- "Serum"

hall$sampleDate <- as.character(lapply(hall$ID, function(x){
   
   curDo <- resAllSave$DonorTime2[resAllSave$globalCluster == x]
   paste(unique(curDo), collapse=",")
   
}))

hall$sampleDate[hall$ID == "ARC"] <- "28.12.21"

#save.image("C:/Users/joh/Desktop/HomeOffice/JC/finaldata3.RData")


hall$structuralClassification <- NULL
hall$structuralClassification2 <- NULL


# Find the largest repeating substring
largest_substring <- find_largest_repeating_substring(s)
print(largest_substring)


## Test bar chart 5A
bars <- unique(data.frame(rep=hall$largestRepeatSeq, donor=hall$donor))
bars2 <- data.frame(rep="", donor=c(32,6,105,24))

bars <- bars[bars$donor != "32,6,105,24",]

bars <- rbind(bars, bars2)

bars <- unique(bars)
bars$c <- 1

bars3 <- data.frame(table(bars$rep))

bars3$len <- lapply(as.character(bars3$Var1), nchar)
bars3$len <- as.numeric(bars3$len)

bars3 <- bars3[order(  bars3$len,bars3$Var1, decreasing = T),]

bars3$Var1 <- as.character(bars3$Var1)

bars3$Var1[bars3$Var1 == ""] <- "None"

bars3$Var1Fac <- factor(bars3$Var1, levels= c(bars3$Var1))

bars$c <- 1



bars$len <- as.numeric(lapply(bars$rep, nchar))
bars <- bars[order(  bars$len,bars$rep, decreasing = T),]

bars$rep <- as.character(bars$rep)

bars$rep[bars$rep == ""] <- "None"

bars$Var1Fac <- factor(bars$rep, levels= unique(bars$rep))
bars$Donor <- bars$donor

bars$Donor2 <- as.character(testD$newID[match(bars$Donor, testD$internID)])
bars$Donor2 <- factor(bars$Donor2, levels=sort(unique(as.numeric(bars$Donor2))))

bars <- bars[bars$rep != "None",]
g <- ggplot(bars, aes(x = rep, y = c, fill=Donor2)) +
   geom_bar(stat = "identity") +
   labs(x = "Longest duplication found in at least one NCCR haplotype", y = "Number of patients with\nduplications",
        fill= "Patient ID") +
   scale_y_continuous(breaks = seq(0, max(bars3$Freq), by = 2)) +
   theme(
      text = element_text(size = 12),
      legend.position = "top",
      legend.justification = "center",
      #legend.title = element_text("Patient ID"), # optional, remove if you want to keep title
      legend.text = element_text(size = 10),
      legend.key.width = unit(.7, "cm"),
      legend.key.height =   unit(.4, "cm")# increase width as needed
   ) +
   guides(fill = guide_legend(nrow = 3)) # make legend items spread horizontally

g

ggsave(g, filename="Figure 7 A.svg", width=11, height=3)


## FIgure 7B

bars2 <- data.frame(donor = unique(bars$Donor))
bars2$A <- 0
bars2$B <- 0
bars2$C <- 0
bars2$D <- 0
bars2$E <- 0
bars2$F <- 0

for(i in bars2$donor){
   
   curRep <- bars$rep[bars$Donor == i]
   curRep <- strsplit(curRep, "")
   
   curRepTest <- lapply(curRep, function(x){ any(x == "A")})
   bars2$A[bars2$donor == i] <- as.numeric(any(unlist(curRepTest)))
   
   curRepTest <- lapply(curRep, function(x){ any(x == "B")})
   bars2$B[bars2$donor == i] <- as.numeric(any(unlist(curRepTest)))
   
   curRepTest <- lapply(curRep, function(x){ any(x == "C")})
   bars2$C[bars2$donor == i] <- as.numeric(any(unlist(curRepTest)))
   
   curRepTest <- lapply(curRep, function(x){ any(x == "D")})
   bars2$D[bars2$donor == i] <- as.numeric(any(unlist(curRepTest)))
   
   curRepTest <- lapply(curRep, function(x){ any(x == "E")})
   bars2$E[bars2$donor == i] <- as.numeric(any(unlist(curRepTest)))
   
   curRepTest <- lapply(curRep, function(x){ any(x == "F")})
   bars2$F[bars2$donor == i] <- as.numeric(any(unlist(curRepTest)))
   
   
}

bars2$outcome <- ""

bars2A <- data.frame(donor=bars2$donor, value=bars2$A, outcome=bars2$outcome, type="A")
bars2B <- data.frame(donor=bars2$donor, value=bars2$B, outcome=bars2$outcome, type="B")
bars2C <- data.frame(donor=bars2$donor, value=bars2$C, outcome=bars2$outcome, type="C")
bars2D <- data.frame(donor=bars2$donor, value=bars2$D, outcome=bars2$outcome, type="D")
bars2E <- data.frame(donor=bars2$donor, value=bars2$E, outcome=bars2$outcome, type="E")
bars2F <- data.frame(donor=bars2$donor, value=bars2$F, outcome=bars2$outcome, type="F")

bars3 <- rbind(bars2A, bars2B, bars2C,bars2D,bars2E,bars2F)

bars4 <- data.frame(tapply(bars3$value, bars3$type, sum))

bars4$Seg <- row.names(bars4) 


names(bars4) <- c("count", "Seg")

bars4$percentage <- round(bars4$count/19*100,1)
 
g <- ggplot(bars4, aes(x = Seg, y = count, fill = Seg)) +
   geom_bar(stat = "identity", width = 0.7) +
   geom_text(aes(label = paste0(percentage, "%")), vjust = 1.1) +
   scale_fill_brewer(palette = "Pastel1") +
   labs(
      x = "Segment fully or partially duplicated in at least on haplotype",
      y = "Number of patients with\nduplications"
   ) +
   #theme_minimal(base_size = 12)  +
   theme(
      #axis.text.x = element_text(angle = 15, hjust = 1),
      legend.position = "none",
      text = element_text(size=15)
   )+
   scale_y_continuous(
      breaks = seq(0, max(bars4$count, na.rm = TRUE), by = 2),  # y-axis ticks every 2
      expand = expansion(mult = c(0, 0.05))  # optional: avoid top clipping
   )
g

ggsave(g, filename="Figure 7 B.svg", width = 6.5,height = 3)



### Figure 7 C
balken <- data.frame(name=c("No duplications", "All duplications\nbuild sequentially",
                            "Combinations of\ntwo haplotypes",
                            "Duplications from\ndistinct indels"),
                     count=c(2, 13, 3, 3),
                     percentage=0)

balken$percentage <- round(balken$count/21*100,1)

balken$nameFakt <- factor(balken$name, levels=c("All duplications\nbuild sequentially",
                                                "Combinations of\ntwo haplotypes",
                                                "Duplications from\ndistinct indels", "No duplications"))

# Plot
g <- ggplot(balken, aes(x = name, y = count, fill = name)) +
   geom_bar(stat = "identity", width = 0.7) +
   geom_text(aes(label = paste0(percentage, "%")), vjust = 1.1) +
   scale_fill_brewer(palette = "Pastel1") +
   labs(
      x = "",
      y = "Number of patients"
   ) +
   #theme_minimal(base_size = 12)  +
   theme(
      #axis.text.x = element_text(angle = 15, hjust = 1),
      legend.position = "none",
      text = element_text(size=13)
   )+
   scale_y_continuous(
      breaks = seq(0, max(bars4$count, na.rm = TRUE), by = 2),  # y-axis ticks every 2
      expand = expansion(mult = c(0, 0.05))  # optional: avoid top clipping
   )
g
ggsave(g, filename="Figure 7 C.svg", width = 5.6, height = 3)



##Figure 7 D-F taken from following supplementary Figure 2
## Supplementary Figure 2
## save plot data for manual plot generation
resAll_list <- list()


for(d in unique(samples$DonorID)){
   
   gC <- data.frame(globalCluster= resAllSave$globalCluster[resAllSave$DonorID == d], 
                    sequence = resAllSave$consensus_seq[resAllSave$DonorID == d])
   
   gCcons <- data.frame(globalCluster=unique(gC$globalCluster), consensus="")
   
   for(g in gCcons$globalCluster){
      
      #g <- 41
      geno <- gC$sequence[gC$globalCluster == g]
      genoTbl <- table(geno)
      gCcons$consensus[gCcons$globalCluster== g] <- names(genoTbl[genoTbl == max(genoTbl)][1])
   }
   names(gCcons)  <- c("name", "sequence")
   haplotypes <- gCcons
   
   haplotypes <- haplotypes[!haplotypes$name %in% c("Others", "108", "113", "117", "26", "Archetype"),]
   
   haplotypes$length <- nchar(haplotypes$sequence)
   haplotypes <- haplotypes[order(haplotypes$length, decreasing = F),]
   haplotypes$length <- NULL
   
   ##Add arc 
   wt <- data.frame(name="ARC", sequence= "GCCTCGGCCTCCTGTATATATAAAAAAAAGGGAAGGTAGGGAGGAGCTGGCTAAAACTGGATGGCTGCCAGCCAAGCATGAGCTCATACCTAGGGAGCCAACCAGCTGACAGCCAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCCCTCCCATAGTCCTTAATCACAAGTAAACAAAGCACAAGGGGAAGTGGAAAGCAGCCAAGGGAACATGTTTTGCGAGCCAGAGCTGTTTTGGCTTGTCACCAGCTGGCC")
   haplotypes <- rbind(haplotypes, wt)
   
   
   resAll <- list()
   yStart <- 100
   
   for(i in 1:nrow(haplotypes)){
      
      res <- get_table(haplotypes$sequence[i])
      
      if(i == nrow(haplotypes)) yStart <- yStart-.5
      
      res$ymin <- yStart-0.5
      res$ymax <- yStart+0.5
      res$ymid <- yStart
      
      res$allMid <- max(res$xmax)/2
      res$allLab <- haplotypes$name[i]
      
      yStart <- yStart-1
      
      resAll[[i]] <- res
   }
   
   resAll <- rbindlist(resAll)
   
   resAll$labIn[resAll$seg %in% "ori"] <- ""
   resAll$labIn[resAll$seg %in% "others"] <- ""
   
   resAll$labIn3[resAll$seg %in% "ori"] <- ""
   resAll$labIn3[resAll$seg %in% "others"] <- ""
   
   
   resAll$ylabminus <- -20
   
   segs <- c("ori", letters[1:6])
   segs2 <- paste0(segs, " (", c(12,25,23,55,66,18,70), "nt)")
   segs2 <- c(segs2, "others")
   segs <- c(segs, "others")
   
   resAll$labIn[resAll$length <=5] <- ""
   
   ### correct genotype related "e f" notation
   resAll$seg[resAll$seg == "e f"] <- "e"
   resAll$donor <- d
   
   resAll_list[[d]] <- resAll
}


names(resAll_list)

resAll <- resAll_list[["116"]]

g <- ggplot(resAll)
g <- g +  geom_rect(aes(xmin = xmin2, xmax = xmax2, ymin = ymin, ymax = ymax, fill = seg)) +
   annotate(geom="text", x=resAll$ylabminus, y=resAll$ymid, label=resAll$allLab, size=6)

g <- g + scale_fill_manual(name = "Segment",
                           values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
                           breaks = segs, labels = segs2)+
   theme(text = element_text(size=20))+xlab("Length (nt)")+
   theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank())

g <- g + annotate(geom="text", x=resAll$mid, y=resAll$ymid, label=resAll$labIn, size=6)
g <- g + ggtitle(paste0("Donor ", unique(resAll$donor)))

height <- 5


height <- (length(unique(resAll$allLab))-1)*.83
height
height <- 6.9
height <- 8.2
ggsave(g, file=paste("Donor ",unique(resAll$donor), " Segment mapping.svg"), width=20, height = height)






## Supplementary Figure 3

## Get first disrupted d segment for alignment

## Redefine getTable function to extract first d segment

get_CSeq <- function(inputSequence){
   
   #inputSequence <- haplotypes$sequence[haplotypes$name ==  "Haplotype_97_maxFreq_1.7%"]
   #inputSequence <- haplotypes$sequence[haplotypes$name ==  "ARC"] 
   #inputSequence <- readClipboard()
   ori <- "GCCTCGGCCTCC"
   a <- "CCTGTATATATAAAAAAAAGGGAAG"
   b <- "GTAGGGAGGAGCTGGCTAAAACT"
   c <- "GGATGGCTGCCAGCCAAGCATGAGCTCATACCTAGGGAGCCAACCAGCTGACAGC"
   d <- "CAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCCCTCCCATAGTCCTTAATCAC"
   e <- "AAGTAAACAAAGCACAAG"
   f <- "GGGAAGTGGAAAGCAGCCAAGGGAACATGTTTTGCGAGCCAGAGCTGTTTTGGCTTGTCACCAGCTGGCC"
   
   getOverlappingVectorsFromVector2 <- function(largeVector, subvectorLength, subvectorOverlap) 
   {
      startPositions = seq(1, length(largeVector), by = subvectorLength - 
                              subvectorOverlap)
      end_positions = startPositions + subvectorLength - 1
      end_positions[end_positions > length(largeVector)] = length(largeVector)
      lapply(seq_len(length(startPositions)), function(x) largeVector[startPositions[x]:end_positions[x]])
   }
   
   
   
   creatSubreadVec <- function(vec, len=10){
      
      vec <- strsplit(vec, "")[[1]]
      vec <- getOverlappingVectorsFromVector2(vec, 10, 9)
      vec <- as.character(lapply(vec, function(x){paste0(x, collapse="")} ))
      vec <- vec[nchar(vec) == 10]
      return(vec)
   }
   
   oriV <- creatSubreadVec(ori)
   aV <- creatSubreadVec(a)
   bV <- creatSubreadVec(b)
   cV <- creatSubreadVec(c)
   dV <- creatSubreadVec(d)
   eV <- creatSubreadVec(e)
   fV <- creatSubreadVec(f)
   
   nameVec <- c(rep("ori", length(oriV)),
                rep("a", length(aV)),
                rep("b", length(bV)),
                rep("c", length(cV)),
                rep("d", length(dV)),
                rep("e", length(eV)),
                rep("f", length(fV)))
   
   nameVecPos <- c(seq(1, length(oriV), 1),
                   seq(1, length(aV), 1),
                   seq(1, length(bV), 1),
                   seq(1, length(cV), 1),
                   seq(1, length(dV), 1),
                   seq(1, length(eV), 1),
                   seq(1, length(fV), 1))
   
   multimer <- data.frame(name= nameVec, pos=nameVecPos, 
                          multimer=c(oriV, aV, bV, cV, dV, eV, fV))
   
   
   ori2 <- substr(ori,nchar(ori)-8 ,nchar(ori))
   a1 <-  substr(a, 1, 9)
   b1 <-  substr(b, 1, 9)
   c1 <-  substr(c, 1, 9)
   d1 <-  substr(d, 1, 9)
   e1 <-  substr(e, 1, 9)
   f1 <-  substr(f, 1, 9)
   
   a2 <-  substr(a,nchar(a)-8 ,nchar(a))
   b2 <-  substr(b,nchar(b)-8 ,nchar(b))
   c2 <-  substr(c,nchar(c)-8 ,nchar(c))
   d2 <-  substr(d,nchar(d)-8 ,nchar(d))
   e2 <-  substr(e,nchar(e)-8 ,nchar(e))
   
   oa <- paste0(ori2,a1)
   ab <- paste0(a2,b1)
   bc <- paste0(b2,c1)
   cd <- paste0(c2,d1)
   de <- paste0(d2,e1)
   ef <- paste0(e2,f1)
   
   oaV <- creatSubreadVec(oa)
   abV <- creatSubreadVec(ab)
   bcV <- creatSubreadVec(bc)
   cdV <- creatSubreadVec(cd)
   deV <- creatSubreadVec(de)
   efV <- creatSubreadVec(ef)
   
   
   nameVec <- c(rep("ori a", length(oaV)),
                rep("a b", length(abV)),
                rep("b c", length(bcV)),
                rep("c d", length(cdV)),
                rep("d e", length(deV)),
                rep("e f", length(efV)))
   
   
   multimerBorders  <- data.frame(name= nameVec, pos= NA, 
                                  multimer=c(oaV, abV, bcV, cdV, deV, efV))
   
   creatPlottingData <- function(promotorSeq){
      
      # promotorSeq <- inputSequence
      
      if(nchar(promotorSeq)<20 )
         stop("ERROR during setting of variable 'promotorSeq'. The entered sequence", 
              " must be at least 20nt long.")  
      if (!all(strsplit(promotorSeq, "")[[1]] %in% c("a", "c", "g", "t", 
                                                     "G", "C", "T", "A"))) 
         stop("ERROR during setting of variable 'promotorSeq'. The entered sequence", 
              " must be a character string of A C G and T.")    
      
      promotorSeq <- toupper(promotorSeq)
      wtV <- creatSubreadVec(promotorSeq)
      
      wtDF <- data.frame(pos= 0, seq=wtV)
      
      wtDF$segment <- multimer$name[match(wtDF$seq, multimer$multimer)]
      wtDF$segmentPos <- multimer$pos[match(wtDF$seq, multimer$multimer)]
      
      wtDF$pos <- 1:nrow(wtDF)
      
      i <- 1
      while(i < nrow(wtDF)){
         
         while(!is.na(wtDF$segment[i]) & i < nrow(wtDF)) i <- i+1
         
         counter <- list()
         
         while( is.na(wtDF$segment[i]) & ((i+1) <= nrow(wtDF))){
            
            counter[[i]] <- i
            i <- i+1
         }
         
         counter <- unlist(counter)
         
         if(length(counter) > 9 ){
            wtDF$segment[counter] <-  multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)] ##test
            
            if(any(!is.na(multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)]))){
               
               pos <- which(!is.na(multimerBorders$name[match(wtDF$seq[counter], multimerBorders$multimer)]))
               
               before10mer <- wtDF$seq[(sort(counter)[1]-1)] 
               after10mer <- wtDF$seq[(sort(counter)[length(counter)]+1)]  
               
               if(before10mer == oriV[length(oriV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "ori a"] <- "a"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               if(before10mer == aV[length(aV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "a b"] <- "b"
                  wtDF$segment[counter[1:9]] <- test
               } 
               
               # if(before10mer == bV[length(bV)]){
               #    test <- wtDF$segment[counter]
               #    test <- test[1:9]
               #    test[test == "b c"] <- "c"
               #    wtDF$segment[counter[1:9]] <- test
               # } 
               
               if(before10mer == cV[length(cV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "c d"] <- "d"
                  wtDF$segment[counter[1:9]] <- test
               }
               
               if(before10mer == dV[length(dV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "d e"] <- "e"
                  wtDF$segment[counter[1:9]] <- test
               }   
               
               if(before10mer == eV[length(eV)]){
                  test <- wtDF$segment[counter]
                  test <- test[1:9]
                  test[test == "e f"] <- "f"
                  wtDF$segment[counter[1:9]] <- test
               }  
               
               lengthcounter <- length(wtDF$segment[counter])
               
               if(after10mer == aV[1]){
                  test <- wtDF$segment[counter]
                  test <- test[(lengthcounter-8):lengthcounter]
                  test[test == "ori a"] <- "ori"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }    
               
               if(after10mer == bV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "a b"] <- "a"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == cV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "b c"] <- "b"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == dV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "c d"] <- "c"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == eV[1]){
                  test <- wtDF$segment[counter]
                  test <-  test[(lengthcounter-8):lengthcounter]
                  test[test == "d e"] <- "d"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               if(after10mer == fV[1]){
                  test <- wtDF$segment[counter]
                  test <- test[(lengthcounter-8):lengthcounter]
                  test[test == "e f"] <- "e"
                  wtDF$segment[counter[(lengthcounter-8):lengthcounter]] <- test
               }  
               
               wtDF$segmentPos[sort(counter)[pos]] <- 999  
            }else{
               
               wtDF$segment[counter] <- "others"
               wtDF$segmentPos[sort(counter)] <- 1:length(wtDF$segment[counter] )
            }
            
            
         }
         #print(i)
         i <- i+1
      }
      
      
      wtDF <- na.omit(wtDF)
      wtDF$row <- 1:nrow(wtDF)
      
      
      dStart <- 1
      while(wtDF$segment[dStart] != "c"){
         
         if(dStart == nrow(wtDF)) break
         
         dStart <- dStart +1
      }  
      
      dEnd <- dStart
      while(wtDF$segment[dEnd] == "c"){
         
         if(wtDF$segment[dEnd+1] != "c") break
         dEnd <- dEnd +1
         
      }     
      
      row.names(wtDF) <- NULL
      
      ## If C is intact but there is a not intact c downstream, take the not intact c seg downstream
      if((dEnd-dStart == 45) & any(wtDF$segment[(dEnd+1):(nrow(wtDF)+1)] == "c", na.rm=T)){
         
         dStart <- dEnd+1
         while(wtDF$segment[dStart] != "c"){
            
            if(dStart == nrow(wtDF)) break
            
            dStart <- dStart +1
         }  
         
         dEnd <- dStart
         while(wtDF$segment[dEnd] == "c"){
            
            if(wtDF$segment[dEnd+1] != "c") break
            dEnd <- dEnd +1
            
         } 
         
      }
      
      if(dStart != nrow(wtDF)){
         
         return(substr(inputSequence, wtDF$pos[dStart] , wtDF$pos[dEnd]+9))
      }else{
         return("")
      }
      
   }
   
   dSegSeq <- creatPlottingData(inputSequence)
   
   return(dSegSeq)
   
}


CSegList <- list()

load("C:/Users/joh/Desktop/HomeOffice/JC/finaldata3.RData")

for(i in 1:nrow(haplotypes)){
   
   print(i)
   res <- get_CSeq(haplotypes$sequence[i])
   
   CSegList[[i]] <- res
}

CSegList <- unlist(CSegList)

sSegDF <- data.frame(name= unlist(haplotypes$cl), dSeg=CSegList )


sSegDF <- sSegDF[!sSegDF$name %in% c("108", "113", "117", "26", "27"),]


sSegDF2 <- sSegDF[sSegDF$dSeg != "",]


i <-1
sSegDF2$alnSeq <- ""

for(i in 1:nrow(sSegDF2)){
   
   curCSeq <- sSegDF2$dSeg[i]
   names(curCSeq) <- sSegDF2$name[i]
   Cseq <- "GGATGGCTGCCAGCCAAGCATGAGCTCATACCTAGGGAGCCAACCAGCTGACAGCCAGAGGGAGCCCTGGCTGCATGCCACTGGCAGTTATAGTGAAACCCCTCCCATAGTCCTTAATCAC"
   names(Cseq) <- "ARC_cSeg"
   
   ttDNAstringset <- DNAStringSet(c(curCSeq,Cseq ))
   ttDNAstringsetaln <- msa(ttDNAstringset, method="Muscle")
   ttDNAstringsetalnD <- DNAStringSet(ttDNAstringsetaln)
   ttDNAstringsetalnD <- as.character(ttDNAstringsetalnD)
   ttDNAstringsetalnD <- ttDNAstringsetalnD[names(ttDNAstringsetalnD) ==  sSegDF2$name[i]]
   ttDNAstringsetalnD <- substr(ttDNAstringsetalnD, 1, 55)
   
   sSegDF2$alnSeq[i] <- ttDNAstringsetalnD
   
}


library(seqinr)


sSegDF2 <- sSegDF2[nchar(sSegDF2$dSeg) != 55,]
sequences <- sSegDF2$alnSeq
namesALT <- sSegDF2$name

for(i in seq_len(length(sequences))){
   
   write.fasta(sequences=sequences[i], names=namesALT[i], 
               as.string=FALSE, file.out="Haplotypes_CSegment_aln_withDels.fasta",
               open = "a",nbchar = 60)
   
}



