#Plotting the average expression for pro- and anti longevity genes

# selects the columns from which the values to be plotted for levels will be taken from
plotCols = grep("H7.8", colnames(Trans)):ncol(Trans)

# removes the "H" from the header and uses its numeric values for the X axis
Hrs = as.numeric(gsub("H", "", colnames(Trans)[plotCols]))

ProLevels = colMeans(Trans[Trans$pro == 1, plotCols],na.rm = T)
AntiLevels = colMeans(Trans[Trans$anti ==1, plotCols],na.rm = T)

pdf("R_Output/AveragePDF.pdf", width = 12, height = 9)

# sets the margins of PDF
par(mar=c(5.1, 4.5, 6.0, 5.1))

plot(Hrs, ProLevels, type="l", lwd=50, col="blue", 
     cex.axis=2, cex.lab=2, yaxt="n", xaxt="n",
     panel.first=grid(col="lightgray", lwd=0.1, nx=NA, ny=NULL), ylab="", xlab="Hours")
abline(v=Hrs, col="lightgray", lwd=0.1, lty=1)
points(Hrs, ProLevels, pch=24, cex=2, col="red", bg="blue")
axis(side=2, cex.axis=2, col.axis="blue", font=2)
axis(side=1, cex.axis=0.8, col.axis="red", font=2, at=Hrs, labels = Hrs)
mtext("Pro-Longevity mRNA",side=2,line=3, col = "blue", font=2, cex=1.3)



par(new=TRUE)
plot(Hrs, AntiLevels, type="l", lwd=50, col="green", 
     ylab="", xlab="", yaxt="n", xaxt="n", lty=1)
points(Hrs, AntiLevels, pch=23, cex=2, col="red", bg="green")
axis(side=4, cex.axis=2, col.axis="green", font=2)
mtext("AntiLongevity mRNA",side=4,line=3, col = "green", font=2, cex=1.3)


if(FALSE){
  main=c(paste(idxTrans[i], subTrans$Feature[i],sep = ": "), 
         subTrans$Feature.Systematic.Name[i], subTrans$Feature.Type[i])
  title(main=main[1], line=3.5, cex.main=2, font=2)
  title(main=main[2], line=2, cex.main=2, font=2)
  title(main=main[3], line=0.5, cex.main=2, font=2)
}
dev.off()

for(f in 1:length(fName)){
OutFileName = paste0("PathwayPlots/",gsub(".csv", ".pdf" ,pName[f]))
if(length(TransChoiceIdx)==0 & length(ProtChoiceIdx) == 0){
  noPlot = c(noPlot, pName[f])
} else {
  plotLevel(idxTrans = TransChoiceIdx, idxProt = ProtChoiceIdx, makePDF=T, PDFname=OutFileName, nr=5, nc=5)
}
}

write.csv(noPlot, "PathwayPlots/noPlot.csv", row.names = F)


# csv output
if(TRUE){
  TransResults = GetChoices(dat=Trans, datname="Trans", verbose = T)
  write.csv(TransResults$NotFound, "R_Output/Found_in_choices_BUT_NOT_in_Trans.csv", row.names=F)
  write.csv(TransResults$Selected, "R_Output/Found_in_choices_AND_Trans.csv", row.names=F)
  
  ProtResults = GetChoices(dat=Prot, datname="Prot", verbose = T)
  write.csv(ProtResults$NotFound, "R_Output/Found in choices BUT NOT in Prot.csv", row.names=F)
  write.csv(ProtResults$Selected, "R_Output/Found in choices AND in Prot.csv", row.names=F)
  
  NeitherResults = unique(c(TransResults$NotFound, ProtResults$NotFound))
  BothResults = TransResults$Selected[TransResults$Selected %in% ProtResults$Selected]
  At_Least_one_plot = unique(c(TransResults$Selected, ProtResults$Selected))
  write.csv(NeitherResults, "R_Output/Found in choices BUT neither in Trans or Prot.csv", row.names=F)
  write.csv(BothResults, "R_Output/Found in choices AND Trans AND Prot.csv", row.names=F)
  write.csv(At_Least_one_plot, "R_Output/At least one plot was made for these genes.csv", row.names=F)
  
}



xxxxx = GetChoices(dat=Trans, datname="Trans", verbose = F)
yyyyy = GetChoices(dat=Prot, datname="Prot", verbose = F)





################## Feature selection part: Gettting Trans Slopes 1 to 11 Finite difference is finding all the derivatives at once.  If we want to see a particular derivative, we have to type changes[[i]]
plotCols                     =  grep("H7.8", colnames(Trans)):ncol(Trans)
Hrs                          =  as.numeric(gsub("H", "", colnames(Trans)[plotCols]))
HrsDiffs                     =  diff(Hrs)
GeneInfo                     =  Trans[,1:(plotCols[1]-1)]
changes                      =  list()
changes[[1]]                 =  t(sapply(1:nrow(Trans), function(x) diff(as.numeric(Trans[x,plotCols]),differences=1)/HrsDiffs))


for(i in 1:(length(Hrs)-2)){
  DifMat                       =  t(sapply(1:nrow(Trans), function(x) diff(changes[[1]][x,],differences=i)))
  changes[[i+1]]               =  cbind(GeneInfo, 'if'(i==10,t(DifMat),DifMat)) 
  colnames(changes[[i+1]])     =  c(colnames(GeneInfo),paste0("T",1:(ncol(changes[[i+1]])-ncol(GeneInfo))))
}
changes[[1]]                   =  cbind(GeneInfo, changes[[1]])
colnames(changes[[1]])         =  c(colnames(GeneInfo),paste0("T",1:(ncol(changes[[1]])-ncol(GeneInfo))))

Ord_highest = order(changes[[11]]$T1, decreasing = T)
Sorted_10th_order_diff = changes[[11]][Ord_highest,]
write.csv(Sorted_10th_order_diff, "sorted gene expression similarity predictor.csv", row.names=F)

U_shaped  =  order(-changes[[10]]$T1, changes[[10]]$T2,decreasing=T)
Sorted_by_the_U_shaped = changes[[10]][U_shaped,]
write.csv(Sorted_by_the_U_shaped, "sorted by U shaped.csv", row.names  = F)

########################## How many tines does the sign of slope switches? ###################

slopes = changes[[1]]
slopeCols = which(colnames(slopes) == "T1"):ncol(slopes)
vals = as.numeric(slopes[1,slopeCols])
gr0 = as.numeric(vals>0)
transition2 = diff(gr0)
transitions = t(apply(slopes[,slopeCols], 1, function(x) diff(as.numeric(x>0))))
nTransitions = cbind(GeneInfo, N=rowSums(abs(transitions)))
proRow   = which(GeneInfo$pro > 0.5)
antiRow  = which(GeneInfo$anti > 0.5)
nTransitions$N[proRow]
#GeneInfo[proRow[which.max(nTransitions$N[proRow])],]
nTransitions$N[antiRow]
Nrange=range(nTransitions$N)
pdf("hist2.pdf")
par(mfrow=c(3,1))
hist(nTransitions$N, col="blue", breaks=seq(Nrange[1],Nrange[2]+1,by=1), main="all", right=F)
hist(nTransitions$N[proRow], col="green", breaks=seq(Nrange[1],Nrange[2]+1,by=1), main="Pro", right=F)
hist(nTransitions$N[antiRow], col="red", breaks=seq(Nrange[1],Nrange[2]+1,by=1), main="Anti", right=F)
dev.off()


########### Working with Ana

write.csv(changes[[1]], "slopes1.csv", row.names = F)
#sum([,slopeCols])
#sum(slopes[,slopeCols])
