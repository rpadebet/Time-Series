###### Reading in the chemical pathway script from Stephen from Tuesday, September 20th 2016 ###

pathways = read.csv("biochemical_pathways.tab", sep = "\t", stringsAsFactors = F)


pathways = pathways[,c("Pathway", "Gene")]
which(duplicated(pathways))
pathways = pathways[!duplicated(pathways),]
pathways = pathways[!pathways$Gene == "",]

pTable = as.data.frame(table(pathways$Pathway),stringsAsFactors = F)
pTable = pTable[order(pTable$Freq,decreasing = TRUE),]
gTable = as.data.frame(table(pathways$Gene), stringsAsFactors = F)
gTable = gTable[order(gTable$Freq,decreasing = TRUE),]

colnames(pTable)[1]="Pathway"
colnames(gTable)[1]="Gene"

clean <- function(x){ tolower(gsub("([^[:alnum:] ])", "", x)) }


for(i in 1:nrow(pTable))
{ 
  pName = pTable$Pathway[i]

# pName = "xylulose degradation"
choices = data.frame(Name = pathways$Gene[pathways$Pathway == pName], stringsAsFactors = F)

fName = paste0("Pathway_choices/", clean(pName),"_choices.csv")
write.csv(choices, fName, row.names = F)
}


###### End of  the chemical pathway script from Stephen from Tuesday, September 20th 2016 ###

pName = list.files("Pathway_choices/")
fName = paste("Pathway_choices/", pName,sep = "")
noPlot = NULL

for(f in 1:length(fName)){
  choices[f,1]  = as.character(read.csv(fName[f], stringsAsFactors = F))
}

####################################################################################################################################################################################################################################################### reading in choices ######################################################################################################################################################################################
