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