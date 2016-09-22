setwd("D:/Tang/12_time_point_dataset")

################## beginning of reading in the GO term table with Surendra Gupta started on Tuesday evening, Sept. 20th, 2016 ##################################

goYeast = read.csv("SGD/gene_association.sgd", sep = "\t", stringsAsFactors = F)
goLookUp = read.csv("Databases/GO/All_GO_Terms_go.obo", sep = "\t", stringsAsFactors = F)







########### end  of reading in the GO term table with Surendra Gupta started on Tuesday evening, Sept. 20th, 2016 ##################################

###### Reading in the chemical pathway script from Stephen from Tuesday, September 20th 2016 ###

setwd("D:/Tang/12_time_point_dataset")

pathways = read.csv("SGD/biochemical_pathways.tab", sep = "\t", stringsAsFactors = F)


pathways = pathways[,c("Pathway", "Gene")]
which(duplicated(pathways))
pathways = pathways[!duplicated(pathways),]
pathways = pathways[!pathways$Gene == "",]
pTable = as.data.frame(sort(table(pathways$Pathway),decreasing = T), stringsAsFactors = F)
gTable = as.data.frame(sort(table(pathways$Gene),decreasing = T), stringsAsFactors = F)

colnames(pTable)[1]="Pathway"
colnames(gTable)[1]="Gene"

clean <- function(x){ tolower(gsub("([^[:alnum:] ])", "", x)) }


for(i in 1:nrow(pTable)){
  
  pName = pTable$Pathway[i]
  # pName = "xylulose degradation"
  choices = data.frame(Name = pathways$Gene[pathways$Pathway == pName], stringsAsFactors = F)
  
  fName = paste0("Pathway_choices/", clean(pName),"_choices.csv")
  write.csv(choices, fName, row.names = F)
  
  ###### End of  the chemical pathway script from Stephen from Tuesday, September 20th 2016 ###
  
  
  ####################################################################################################################################################################################################################################################### reading in choices ######################################################################################################################################################################################
  
  ###### Beginning of the regular plotting script by Stephen Nawara as of Sept. 20th, 2016 ###
  
  setwd("D:/Tang/12_time_point_dataset")
  
  pName = list.files("Pathway_choices/")
  fName = paste0("Pathway_choices/", pNames)
  noPlot = NULL
  for(f in 1:length(fName)){
    
    
    #choices = read.csv(file.choose(), stringsAsFactors = F)
    choices  = read.csv(fName[f], stringsAsFactors = F)
    
    #choices=data.frame(Name=antiGenes)
    AllChr   = read.csv("S288c_gene_info/AllChr.csv", stringsAsFactors = F)
    Wuttke   = read.csv("Wuttke/Complete_list_of_lifespan_changing_genes_in_yeast_by_Wuttke.tsv", sep = "\t", stringsAsFactors = F)
    Trans    = read.csv("Table_Fig2_S5_final_transcriptome_after_having_deleted_the_two_genes_which_are_not_in_AllChr.csv", 
                        stringsAsFactors = F)
    Prot     = read.csv("Table_Fig2_S4_Proteome_without_semicolon_1443_genes.csv", 
                        stringsAsFactors = F)
    
    # Deleting not needed columns from Wuttke (long text strings)
    toDel=names(which(sapply(Wuttke,function(x) length(unique(x))) == 1))
    rem = colnames(Wuttke) %in% c(toDel, "Phenotype.Description", "Bibliographic.reference", "Gene.Name") 
    Wuttke = Wuttke[,!rem]
    
    
    
    
    
    # converting characters of the columns below into lower case and remove white space
    choices$Name                   = tolower(gsub(" ", "", choices$Name))
    AllChr$Feature                 = tolower(gsub(" ", "", AllChr$Feature))
    AllChr$Feature.Systematic.Name = tolower(gsub(" ", "", AllChr$Feature.Systematic.Name))
    Trans$ORF                      = tolower(gsub(" ", "", Trans$ORF))
    Prot$ORF                       = tolower(gsub(" ", "", Prot$ORF))
    Wuttke$Gene.Symbol             = tolower(gsub(" ", "", Wuttke$Gene.Symbol))
    
    
    
    # Remove duplications from Wuttke
    WGenes = unique(Wuttke$Gene.Symbol)
    probs  = data.frame(matrix(nrow=length(WGenes), ncol=4))
    
    for(i in 1:length(WGenes)){
      sub    = Wuttke[Wuttke$Gene.Symbol==WGenes[i], ]
      P = sapply(c("pro", "anti", "fitness"),function(x) sum(sub$Longevity.Influence == x)/nrow(sub))
      probs[i,] = c(WGenes[i], P)
    }
    probs[,2:4] = sapply(probs[,2:4],as.numeric)
    colnames(probs) = c("Common_Name", "pro", "anti", "fitness")
    
    # removing the lines with a semicolon in column A from Prot
    Prot = Prot[!grepl(";",Prot$ORF),]
    
    
    
    # Matching the rows of Trans with AllChr
    mTrans = match(Trans$ORF, AllChr$Feature.Systematic.Name)
    # Matching the rows of Prot with AllChr
    mProt = match(Prot$ORF, AllChr$Feature.Systematic.Name)
    
    
    # Adding all columns of AllChr to Trans
    Trans  = cbind(AllChr[mTrans,], Trans)
    # Adding all columns of AllChr to Prot
    Prot   = cbind(AllChr[mProt,], Prot)
    
    mProbTrans = match(Trans$Feature, probs$Common_Name)
    mProbProt = match(Prot$Feature, probs$Common_Name)
    
    
    Trans = cbind(probs[mProbTrans,], Trans)
    Prot = cbind(probs[mProbProt,], Prot)
    
    
    # Giving it the row names of Trans
    rownames(Trans)=1:nrow(Trans)
    # Giving it the row names of Prot
    rownames(Prot)=1:nrow(Prot)
    
    
    
    
    # Finding mis-matches between Trans and AllChr
    NoMatchTrans=which(is.na(mTrans))
    if(length(NoMatchTrans)>0){
      MissingTrans=Trans$ORF[NoMatchTrans]
      cat(paste(length(MissingTrans), "are found in Trans BUT NOT in AllChr.\n These are", 
                paste(MissingTrans, collapse=", ")))
    }
    # Finding mis-matches between Prot and AllChr
    NoMatchProt=which(is.na(mProt))
    if(length(NoMatchProt)>0){
      MissingProt=Prot$ORF[NoMatchProt]
      cat(paste(length(MissingProt), "are found in Prot but not in AllChr.\n These are", 
                paste(MissingProt, collapse=", ")))
    }
    
    # Get Choices
    
    GetChoices <- function(dat, datname, verbose = T){
      
      # Matches the string of each row in choices with the common name column in Trans
      CommonIdx = match(choices$Name, dat$Feature)
      
      # Matches the string of each row in choices with the systematic name column of Trans 
      SystIdx = match(choices$Name, dat$Feature.Systematic.Name)
      
      # Combines the row numbers of Trans of CommonIdx and SystIdx into choiceIdx and removes the NAs
      choiceIdx = c(CommonIdx, SystIdx)
      choiceIdx = unique(choiceIdx[!is.na(choiceIdx)])
      
      
      ########### What is below is not needed for plot selection but only for error messages ##########
      Common = dat$Feature[CommonIdx]
      Syst = dat$Feature.Systematic.Name[SystIdx]
      
      Categories = cbind(choices, Common, Syst, 
                         Missing = as.numeric(is.na(Common) & is.na(Syst)))
      
      NotFound=Categories$Name[Categories$Missing==1]
      Selected=Categories$Name[Categories$Missing==0]
      
      if(length(NotFound)>0 & verbose == T){
        cat(paste0(length(NotFound), " are found in choices BUT NOT in ", datname,".\n These are ",
                   paste(NotFound, collapse = ", "),"\n\n"))
      }
      if(length(Selected)>0 & verbose == T){  
        cat(paste0(length(Selected), " are found in both, choices AND ",datname,".\n These are ", 
                   paste(Selected, collapse=", "),"\n\n"))
      }
      ########### What is above is not needed for plot selection but only for error messages ##########
      return(list(choiceIdx = choiceIdx, NotFound = NotFound, Selected = Selected))
    }
    
    
    TransChoiceIdx = GetChoices(dat=Trans, datname="Trans", verbose = T)$choiceIdx
    ProtChoiceIdx  = GetChoices(dat=Prot, datname="Prot", verbose = T)$choiceIdx
    
    
    
    # Plot function
    plotLevel=function(idxTrans, idxProt, makePDF=T, PDFname="if running plotLevel alone.pdf", nr=3, nc=3){
      # selects the rows from Trans specified by idx for plotting
      subTrans      = Trans[idxTrans,]
      subProt       = Prot[idxProt,]
      
      if(nrow(subTrans)==0 & nrow(subProt)==0){
        #stop("Nothing to plot", call.=F)
        print("Nothing to plot")
        break  
      }
      
      # selects the columns from which the values to be plotted for levels will be taken from
      plotCols = grep("H7.8", colnames(subTrans)):ncol(subTrans)
      
      # removes the "H" from the header and uses its numeric values for the X axis
      Hrs = as.numeric(gsub("H", "", colnames(subTrans)[plotCols]))
      
      # Initializes the PDF and its size
      if(makePDF){pdf(PDFname, width = 120, height = 90)}
      
      # tells the rows and columns of graphs per printed PDF page
      
      par(mfrow=c(nr,nc))
      
      # sets the margins of PDF
      par(mar=c(5.1, 4.5, 6.0, 5.1))
      
      # goes through all the steps for printing the graph for one gene
      if(nrow(subTrans)>0){
        for(i in 1:nrow(subTrans)){
          
          # determines the Y values of the plots
          level = subTrans[i,plotCols]
          PlotProt = which(subProt$ORF == subTrans$ORF[i])
          if(length(PlotProt)>0){
            level2 = subProt[PlotProt,plotCols]
          }
          
          
          # Makes the plots in the following order, x values, y values, draws the line, line width, line color
          plot(Hrs, level, type="l", lwd=50, col="blue", 
               yaxt="n", xaxt="n",
               panel.first=grid(col="lightgray", lwd=0.1, nx=NA, ny=NULL), ylab="", xlab="Hours")
          abline(v=Hrs, col="lightgray", lwd=0.1, lty=2)
          points(Hrs, level, pch=24, cex=2, col="red", bg="blue")
          axis(side=2, cex.axis=2, col.axis="blue", font=2)
          axis(side=1, cex.axis=0.8, col.axis="red", font=2, at=Hrs, labels = Hrs)
          mtext("mRNA Levels",side=2,line=3, col = "blue", font=2, cex=1.3)
          
          
          if(length(PlotProt)>0){
            par(new=TRUE)
            plot(Hrs, level2, type="l", lwd=50, col="darkgreen", 
                 ylab="", xlab="", yaxt="n", xaxt="n", lty=1)
            points(Hrs, level2, pch=23, cex=2, col="red", bg="green")
            axis(side=4, cex.axis=2, col.axis="green", font=2)
            mtext("Protein Levels",side=4,line=3, col = "green", font=2, cex=1.3)
            
          }
          
          main=c(paste(idxTrans[i], subTrans$Feature[i],sep = ": "), 
                 subTrans$Feature.Systematic.Name[i], subTrans$Feature.Type[i])
          title(main=main[1], line=3.5, cex.main=2, font=2)
          title(main=main[2], line=2, cex.main=2, font=2)
          title(main=main[3], line=0.5, cex.main=2, font=2)
        }
      }
      
      
      
      if(nrow(subProt)>0){
        for(i in 1:nrow(subProt)){
          PlotProt = which(subTrans$ORF == subProt$ORF[i])
          if(length(PlotProt)==0){
            level2 = subProt[i,plotCols]
            
            
            
            plot(Hrs, level2, type="l", lwd=50, col="green", 
                 cex.axis=2, cex.lab=2, yaxt="n", lty=1, xaxt = "n",
                 panel.first=grid(col="lightgray", lwd=0.1, nx=NA, ny=NULL), ylab="", xlab="Hours")
            points(Hrs, level2, pch=23, cex=2, col="red", bg="green")
            axis(side=4, cex.axis=2, col.axis="green", font=2)
            axis(side=1, cex.axis=0.8, col.axis="red", font=2, at=Hrs, labels = Hrs)
            mtext("Protein Levels",side=4,line=3, col = "green", font=2, cex=1.3)
            abline(v=Hrs, col="lightgray", lwd=0.1, lty=2)
            
            
            main=c(paste(idxProt[i], subProt$Feature[i],sep = ": "), 
                   subProt$Feature.Systematic.Name[i], subProt$Feature.Type[i])
            title(main=main[1], line=3.5, cex.main=2, font=2)
            title(main=main[2], line=2, cex.main=2, font=2)
            title(main=main[3], line=0.5, cex.main=2, font=2)
          }
        }  
      }
      if(makePDF){dev.off()}
    }
    
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
    
    OutFileName = paste0("PathwayPlots/",gsub(".csv", "Rohit.pdf" ,pName[f]))
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
    write.csv(ProtResults$Not}Found, "R_Output/Found in choices BUT NOT in Prot.csv", row.names=F)
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
sum([,slopeCols])
sum(slopes[,slopeCols])


GeneInfo[proRow[which(nTransitions$N[proRow]==3)],]
nTproRow = nTransitions$N[proRow]
aTprowRow = sum(nTproRow)/length(nTproRow)

H7.8minusH72.3 = Trans$H7.8-Trans$H72.3
order(H7.8minusH72.3)
Ord_highest = order(changes[[11]]$T1, decreasing = T)
Sorted_10th_order_diff = changes[[11]][Ord_highest,]

########################### Calulating net change ####################

Trans
score      = as.numeric(Trans$H72.3 > Trans$H7.8)
proScore   = score[which(GeneInfo$pro > 0.5)]
antiScore  = score[which(GeneInfo$anti > 0.5)]


################### Calculating the proportion of plots whose absolute value is larger in the last than the second last intreval


#score      = as.numeric(abs(changes[[1]]$T11) > abs(changes[[1]]$T10))
#proScore   = score[which(GeneInfo$pro > 0.5)]
#antiScore  = score[which(GeneInfo$anti > 0.5)]
prop <- function(s){ 
  n= length(s)
  c(prop = sum(s)/n, n = n)
}
res = cbind(overall = prop(score), pro = prop(proScore), anti = prop(antiScore))


test.results = list()
# pro - anti
test.results[[1]] = prop.test(x = res[1,2:3]*res[2,2:3], n=res[2,2:3])
names(test.results[[1]]$estimate)=c("pro", "anti")

# pro - Overall
test.results[[2]] = prop.test(x = res[1,c(2,1)]*res[2,c(2,1)], n=res[2,c(2,1)])
names(test.results[[2]]$estimate)=c("pro", "overall")

# anti - Overall
test.results[[3]] = prop.test(x = res[1,c(3,1)]*res[2,c(3,1)], n=res[2,c(3,1)])
names(test.results[[3]]$estimate)=c("anti", "overall")

res[1,2]-res[1,1]

require(binom)
ci.results = binom.confint(x = res[1,]*res[2,], n=res[2,])
ci.results = cbind(group = colnames(res), ci.results)
ci.results = ci.results[order(ci.results$group),]
###

proGenes = Trans$Feature[which(Trans$pro>0.5)]
antiGenes = Trans$Feature[which(Trans$anti>0.5)]

write.csv (proGenes, "25 pro genes.csv", row.names = F)
write.csv (antiGenes, "201 anti genes.csv", row.names = F)



################################# Our old way to calculate the slopes #################################

plotCols=grep("H7.8", colnames(Trans)):ncol(Trans)
Hrs = as.numeric(gsub("H", "", colnames(Trans)[plotCols]))
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
  
  
  
  ##### all slopes will be given in TransSlopes
  
  
  ######  Queries
  
  
  
  