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
topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)
results
vennDiagram(results)
topTableF(fit2, number=30)





myclusters <- kmeans()






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


