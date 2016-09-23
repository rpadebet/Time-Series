

#GeneInfo[proRow[which(nTransitions$N[proRow]==3)],]
nTproRow = nTransitions$N[proRow]
aTprowRow = sum(nTproRow)/length(nTproRow)

H7.8minusH72.3 = Trans$H7.8-Trans$H72.3
order(H7.8minusH72.3)
Ord_highest = order(changes[[11]]$T1, decreasing = T)
Sorted_10th_order_diff = changes[[11]][Ord_highest,]

########################### Calulating net change ####################

#Trans
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

#res[1,2]-res[1,1]


ci.results = binom.confint(x = res[1,]*res[2,], n=res[2,])
ci.results = cbind(group = colnames(res), ci.results)
ci.results = ci.results[order(ci.results$group),]
###

proGenes = Trans$Feature[which(Trans$pro>0.5)]
antiGenes = Trans$Feature[which(Trans$anti>0.5)]

write.csv (proGenes, "25 pro genes.csv", row.names = F)
write.csv (antiGenes, "201 anti genes.csv", row.names = F)
