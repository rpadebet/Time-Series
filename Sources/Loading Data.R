#setwd("D:/Tang/12_time_point_dataset")
require(binom)
################## beginning of reading in the GO term table with Surendra Gupta started on Tuesday evening, Sept. 20th, 2016 ##################################

goYeast = read.csv("gene_association.sgd", sep = "\t", stringsAsFactors = F)
goLookUp = read.csv("All_GO_Terms_go.obo", sep = "\t", stringsAsFactors = F)
AllChr   = read.csv("AllChr.csv", stringsAsFactors = F)
Wuttke   = read.csv("Complete_list_of_lifespan_changing_genes_in_yeast_by_Wuttke.tsv", sep = "\t", stringsAsFactors = F)
Trans    = read.csv("Table_Fig2_S5_final_transcriptome_after_having_deleted_the_two_genes_which_are_not_in_AllChr.csv", 
                    stringsAsFactors = F)
Prot     = read.csv("Table_Fig2_S4_Proteome_without_semicolon_1443_genes.csv", 
                    stringsAsFactors = F)

########### end  of reading in the GO term table with Surendra Gupta started on Tuesday evening, Sept. 20th, 2016 ##################################



