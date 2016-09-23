################################# Our old way to calculate the slopes #################################

plotCols=grep("H7.8", colnames(Trans)):ncol(Trans)
Hrs = as.numeric(gsub("H", "", colnames(Trans))[plotCols])
TransSlopes=data.frame(matrix(nrow=nrow(choices), ncol=ncol(Trans)-1))



for (i in 1:nrow(choices)){
  
  
  sub              =  Trans[Trans$Feature == choices$Name[i],]
  if(nrow(sub)==0){
    sub              =  Trans[Trans$Feature.Systematic.Name == choices$Name[i],]
  }
  
  if(nrow(sub)>0){
    
    diffs            =  diff(as.numeric(sub[,plotCols]))
    HrsDiffs         =  diff(Hrs)
    slopes11         =  diffs/HrsDiffs
    
    TransSlopes[i,1:(plotCols[1]-1)]      =  sub[,-plotCols] 
    TransSlopes[i,-(1:(plotCols[1]-1))]   =  slopes11
  } 
}
colnames(TransSlopes)    =  colnames(Trans)[-ncol(Trans)]
TransSlopes    =  TransSlopes[!is.na(TransSlopes$H7.8),]
rownames(TransSlopes) = 1:nrow(TransSlopes)
write.csv(TransSlopes, "11 Slopes for all choices with Trans data.csv", row.names=F)

C1<-TransSlopes[,c(2,3,5,6,12:ncol(TransSlopes))]

head(C1)
ColN = colnames(C1)<-c("pro","anti", "Feature", "Feature.Systematic.Name", "I1C1", "I2C1", "I3C1", "I4C1", "I5C1", "I6C1", "I7C1", "I8C1", "I9C1", "I10C1", "I11C1")



C1[1,-(1:4)]

plotColsC2=grep("H7.8", colnames(Trans)):ncol(Trans)
Hrs = as.numeric(gsub("H", "", colnames(Trans)[plotCols]))
TransSlopes=data.frame(matrix(nrow=nrow(choices), ncol=ncol(Trans)-1))
for (i in 1:nrow(choices)){
  
  max(Trans[,12:23])
  
  
  plotCols2 = plotCols[-length(plotCols)]
  TransAbsSlopes              =  TransSlopes
  TransAbsSlopes[,plotCols2]   =  abs(TransSlopes[,plotCols2])
  
  write.csv(TransAbsSlopes, "Abs 11 Slopes for all choices with Trans data.csv", row.names=F)
  
  
  write.csv(TransSlopes[TransSlopes$pro == 1, ], "11 Slopes for Pro Longevity genes.csv", row.names=F)
  write.csv(TransAbsSlopes[TransAbsSlopes$pro == 1 ], "Abs for 11 Slopes for Pro longevity Genes.csv", row.names=F)
  write.csv(TransSlopes[TransSlopes$anti == 1, ], "11  for Anti Longevity genes.csv", row.names=F)
  write.csv(TransAbsSlopes[TransAbsSlopes$anti == 1, ], "Abs Slopes anti.csv", row.names = F)
  
}