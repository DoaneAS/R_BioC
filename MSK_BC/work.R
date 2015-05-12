close(exprsGSE31519)
rm(exprsGSE31519)
rm(pdataGSE31519)

niceCols <- 
display.brewer.pal(3, "Set1")
display.brewer.pal(3, "Pastel1")
display.brewer.pal(3, "Dark2")
display.brewer.pal(3, "Set3")


###ERP A B analysis (3 way) LIMMA ####
f <- factor(MSK.set$class, levels=c("P","A","B"))
design <- model.matrix(~0+f)

colnames(design) <- c("P","A","B")

fit <- lmFit(MSK.set, design)
contrast.matrix <- makeContrasts(A-P, B-A, B-P, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)
results
vennDiagram(results)
topTableF(fit2, number=30)





myclusters <- kmeans()




##CONSENSUS CLUSTER WORK

d = exprs(MSK.set)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:7000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))


title=getwd()
results = ConsensusClusterPlus(d,maxK=7,reps=1000,pItem=0.8,pFeature=1,verbose = TRUE,
                                 title=title,
                                 clusterAlg="kmdist",
                                 distance="euclidean",
                                 seed=1262118388.71279,
                                 plot="pdf")

icl = calcICL(results,title=title,plot="pdf")
icl[["clusterConsensus"]]


title=getwd()
results = ConsensusClusterPlus(d,maxK=6,reps=1000,pItem=0.8,pFeature=1,verbose = TRUE,
                               title=title,
                               clusterAlg="hc",
                               distance="pearson",
                               seed=1262118388.71279,
                               plot="pdf")

icl = calcICL(results,title=title,plot="pdf")
icl[["clusterConsensus"]]


title=getwd()
results2 = ConsensusClusterPlus(d,maxK=10,reps=1000,pItem=0.8,pFeature=0.8,
                               title=title, 
                               clusterAlg="kmdist",
                               distance="euclidian",
                               seed=1262118388.71279,
                               plot="pdf")

icl = calcICL(results2,title=title,plot="pdf")
icl[["clusterConsensus"]]








##Subsetting ######

for(id in 1:nrow(tab)){
  exp$symbol[exp$Affy %in% tab$Affy[id]] <- tab$symbol[id]
  exp$genename[exp$Affy %in% tab$Affy[id]] <- tab$genename[id]
}


mutate(exp, Affy = rownames(exp), symbol = tab[rownames(exp),]$symbol, genename = as.character(tab[rownames(exp),]$genename))



expAnn <- exp %>%
  mutate(Affy = tab[rownames(exp),]$Affy, symbol = as.character(tab[rownames(exp),]$symbol))

expAffy$symbol

expSym <- exp %>%
  mutate(symbol = as.character(tab[rownames(exp),]$symbol))

expSym$symbol
  
  inner_join(tab)
    
    
    symbol = tab[rownames(exp),]$symbol,
    genename = as.character(tab[rownames(exp),]$genename))
  
exp %>% mutate(symbol = tab[rownames(exp),]$symbol)

exp = data.frame(y[rownames(tab),])

exp_anx <- exp %>%
  mutate(Affy = tab[rownames(exp),]$Affy) %>%
  inner_join(tab)

exp_anx %>% inner_join(tab)







#using dplyr with mutate and inner_join:##########
###WORKING##
###
exp = data.frame(y[rownames(tab),])
exp <- tbl_df(exp)

exp_ann_all <- exp %>%
  mutate(Affy = tab[rownames(exp),]$Affy) %>%
  inner_join(tab) %>%
  mutate(sym_affy = paste(symbol, affy, sep= '_'))

exp_a <- exp_ann_all %>%
  select(sym_affy, starts_with("GSM"))
rownames(exp_a) <- exp_a$sym_affy
exp_a$sym_affy <- NULL

texp = as.data.frame(t(exp_a))
texp = tbl_df(data.frame(t(exp_a)))
head(texp)



x <- dcast(melt(exp_a), ...~sym_affy)

names(res)

nm = names(res)[2:11]
for (i in seq_along(nm)) {
  ggplot(res,aes_string(x = nm[i])) + geom_boxplot(alpha = .5,fill = "dodgerblue")

}

library(ggplot2)

ggplot(res, aes( AGR2_209173_at)) + geom_boxplot()

boxplotplus2(res)
boxplot(res)
boxplot(res[c(1,2)])
