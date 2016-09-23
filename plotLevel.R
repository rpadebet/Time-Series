
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
      title(main=main[1], line=3.5, cex.main=4, font=2)
      title(main=main[2], line=2, cex.main=4, font=2)
      title(main=main[3], line=0.5, cex.main=4, font=2)
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
        title(main=main[1], line=3.5, cex.main=4, font=2)
        title(main=main[2], line=2, cex.main=4, font=2)
        title(main=main[3], line=0.5, cex.main=4, font=2)
      }
    }  
  }
  if(makePDF){dev.off()}
}