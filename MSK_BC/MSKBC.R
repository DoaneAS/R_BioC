
cel.dir= "~/R_BiocProjects.AD/MSK_BC/GSE2603_RAW/"
output.dir = getwd()

process.raw.affy.data(cel.dir, output.dir)

###Normalize arrays
#######################
library(limma)
ds=dd$ds
dsQ=normalizeBetweenArrays(ds,method='quantile')

dd$dsQ <- dsQ
dd$dsR <- ds
dd$ds <- dsQ

pdata <- read.csv("ids.csv", header = TRUE, stringsAsFactors = FALSE)
p <-  read.csv("pheno1.csv", header = TRUE, stringsAsFactors = FALSE)
p <- p[1:99,]

pdata$class[match(p$id, pdata$ids)] <- p$class
row.names(pdata) <- pdata$gsm
pdata$class <- factor(pdata$class, levels= c("P", "A", "B"))
s <- rownames(pdata)
x <- match(s , colnames(dd$ds))
MSKexp <- dd$ds[, x]


#### can also matrch the revserse way?
match(colnames(dd$ds), pdata$gsm) 
i <- match(pdata$gsm, colnames(dd$ds))
exp <- dd$ds[,i]
#####

#####
##Define phenodata ##
#####

##add

pdata$ER[pdata$class == "A" | pdata$class == "B" ] <- "N"
pdata$ER[pdata$class == "P"] <- "P"

pdata$ER <- factor(pdata$ER, levels = c("N", "P"), labels = c("ER+", "ER-"))

#pheno
p <- colnames(pdata)
varInfo <- data.frame(labelDescription= (p), row.names=(p))
phenoMSK <- new("AnnotatedDataFrame", data=pdata, varMetadata = varInfo) #pData for expressionset

#####

MSK.set <- new("ExpressionSet", exprs=MSKexp,
                     phenoData=phenoMSK,
                     annotation="hgu133a.db")



MSK.set <- fAnnot(MSK.set)

save(MSK.set, file="MSKset.Rdata",compress=T)








#########
########

###################

pdf("UnsupclusterAllnorm.pdf", width=32, height=10)
distance <- dist(t(dd$ds),method="euclidian")
clusters <- hclust(distance, method="ward.D2")
par(cex=0.5)

plot(clusters, hang=-1)
dev.off()
clust <- as.dendrogram(clusters)

#####
d <- exprs(MSK.set)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
Index <- as.factor(pData(MSK.set)$class)
ser <- levels(MSK.set$class)
labcol <- ser[Index]
or <- order(ser[Index])
niceCols <- brewer.pal(5, "Set1")
z <- niceCols[Index]
###HEATMAP####
pdf("Cluster.unsupdend.rescale2.pdf",width=36,height=24)
heatmap.2(d,
          #labCol=labcol,
          scale="row", 
          #dendrogram="none", 
          Colv=clust,
          col=my_palette, trace= "none", density="density",
          ColSideColors=z,
          key=FALSE)
dev.off()



Class <- phenoData(MSK.set)$class
classA <- grep("A", Class)
classB <- grep("B", Class)
classP <- grep("P", Class)

ER <- phenoData(MSK.set)$ER
ERp <- grep("P", ER)
ERn <- grep("N", ER)

classAExp <- MSK.set[,classA]  
classBExp <- MSK.set[,classB]  
classPExp <- MSK.set[,classP]
ernExp <- MSK.set[,ERn]
erpExp <- MSK.set[,ERp]

table(pData(MSK.set)$class)
Index <- as.numeric(pData(MSK.set)$ER)
Index2 <- as.numeric(pData(MSK.set)$class)
y = exprs(MSK.set)
library(limma)

ERgroup <- factor(pData(MSK.set)[,4 ] , levels = levels(pData(MSK.set)[,4]))

#design <- model.matrix(~factor(MSK.set$ER))
#or
design <- model.matrix(~ERgroup)
colnames(design) <- c("ER-", "ER+vER-")
fit = lmFit(MSK.set, design)
ebayes = eBayes(fit)
tab <- topTable(ebayes, coef="ER+vER-", adjust="BH", n=50, genelist =  fit$genes$symbol)

labCol <- c("N", "P")[Index]
or <- order(ser[Index])
niceCols <- brewer.pal(5, "Set1")
z <- niceCols[Index]
#labCol <- c("P", "A", "B")[Index2]

 y = exprs(MSK.set)


pdf("Cluster.unsupdend.rescale.pdf",width=36,height=48)
heatmap.2(y[rownames(tab),],
          labRow= tab$ID,
          labCol = labCol,
          col=my_palette, trace= "none",
          density="density",
          ColSideColors=z,
          key=FALSE)
dev.off()



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
tab = topTableF(fit2, number=175)

#heatmap setup
Index2 <- as.numeric(pData(MSK.set)$class)
labCol <- c("P", "A", "B")[Index2]
ser <- levels(MSK.set$class)
or <- order(ser[Index2])
niceCols <- brewer.pal(5, "Set1")
z <- niceCols[Index2]

pdf("Cluster.PvAvBLimma.pdf",width=36,height=48)
heatmap.2(y[rownames(tab),],
          labRow= tab$symbol,
          labCol = labCol,
          col=my_palette, trace= "none",
          density="density",
          ColSideColors=z,
          key=FALSE)
dev.off()



#er n only#

table(pData(ernExp)$class)

xx <- y[tab[,"ID"],]
dim(xx)

#create exp set of sig genes only
selected <- p.adjust(fit$p.value[, 2]) <0.05
esetSel <- MSK.set[selected,]
esetSel


