###COMMON GENE EXPRESSION FUNCTIONS AND SCRITPS###
##################################################


###################
#GCRMA processing
#####################
process.raw.affy.data<-function(cel.dir, output.dir){
  library(gcrma)    
  library(affy)
  
  ## set output file names based on date
  rdata.file = paste(geo.number,"GCRMA.Rdata",sep="_")
  
  cat(c("Reading .CEL files in ",cel.dir,"..."))
  setwd(cel.dir)
  
  ## read cel files in current working directory
  cel.files=ReadAffy(celfile.path=".")
  cat(c("Done.","\n"))
  
  ## run GCRMA on raw data
  cat(c("Running GCRMA...","\n"))
  start=proc.time()
  dd=gcrma(cel.files)
  cat(c("Done.","\n"))
  proc.time()-start
  cat(c("Making detection calls...","\n"))
  ## make detection calls
  dc=mas5calls(cel.files)
  cat(c("Done.","\n"))
  proc.time()-start
  
  cat(c("Creating and saving data matrices..."))
  ## extract normalized values
  ds=exprs(dd)
  dp=exprs(dc)
  
  ## remove .CEL extension from sample names in column headers
  names=colnames(ds)
  pos=regexpr(".cel|.CEL",names)
  colnames(ds)=substr(names,1,pos-1)
  colnames(dp)=substr(names,1,pos-1)
  
  ## remove probes containing "AFFX"
  ii=grep("AFFX",rownames(ds))
  ds=ds[-ii,]
  dp=dp[-ii,]
  
  ## create data frame of normalized expression and P/A/M calls
  dd=NULL
  dd$ds=ds
  dd$dp=dp
  
  setwd(output.dir)
  
  ## save to Rdata file
  save(dd,file=rdata.file,compress=T)
  cat(c("Done.","\n"))
  
  return(dd)
}





load("~/R Bioconductor Data/TNBC/metabList2.Rdata")
library("survival")
library("genefu")
library("rmeta")
library(genefilter)
load("~/R Bioconductor Data/BC.Koutcher/metabList.Rdata")
load("~/R Bioconductor Data/BC.Koutcher/Latest/TNBC.combat.set.Rdata")
#################################
####SURVIVAL ANALYSIS####
################################
surv.analysis<-function(gn.expr,status,stime,method="two_thirds", alt="g")
{
  if(method=="median")
    m1=median(gn.expr)
  if(method=="two_thirds")
  {
    o1=order(gn.expr)
    m1=gn.expr[o1[round(1*length(gn.expr)/3)]]
  }    
  gn.status=rep(NA,length(gn.expr))
  gn.status[which(gn.expr<=m1)]="low"
  gn.status[which(gn.expr>m1)]="high"
  gn.status1=gn.status
  gn.status=as.factor(gn.status)
  
  m.low=median(gn.expr[gn.expr<=m1])
  m.high=median(gn.expr[gn.expr>m1])
  
  dat.list=list(time=stime,status=status,gn.status=gn.status)
  mfit=survfit(Surv(time,status)~gn.status,data=dat.list)
  test.diff=survdiff(Surv(time,status)~gn.status,data=dat.list)
  #print(survdiff(Surv(time,status)~gn.status,data=dat.list))
  #test.res=surv_test(Surv(time, status) ~ gn.status, data = dat.list,alt=alt)
  #print(test.res)
  
  ###calculate p-value ### 
  #pval=getPval(test.diff)
  #print(pval)
  
  #pval2 <- 1 - pchisq(test.diff$chisq, length(test.diff$n) - 1) 
  #print(paste("pval2=",pval2))
  #return(list(mfit=mfit,test.res=test.res,test.diff=test.diff,pval=pval,pval2=pval2))
  return(list(mfit=mfit,test.diff=test.diff,m.low=m.low,m.high=m.high))
}

#####################
####Rescale 
#####################
rescale <- function(x, na.rm=FALSE, q=0.05) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return((x - 0.5) * 2)
}

#######


esetALL <-tnbc.Allset  #TNBC.combat.set
#esetALL.filt <- nsFilter(esetALL, remove.dupEntrez=TRUE, filterByQuantile=TRUE,  var.cutof =0.95)
#esetALL.filt <- esetALL.filt$eset
esetALL.filt <- tnbc.Allset #if no filtering applied here

#  VTMUUSsurvset
#  TN.UnTx.set
#  TNset
#esetALL  <- TNBCsurVset #expression set, this is GSE31519
  
tc= 5*365
surv.data <- censor.time(surv.time=c(pData(esetALL)[ ,"stime"])/12, surv.event=c(pData(esetALL)[ ,"status"]), time.cens=tc / 365)

#Define survival data
status <- surv.data$surv.event.cens
stime <- surv.data$surv.time.cens
surv <- Surv(stime, status)
##setup survival eset
#esetALL <- TNBC.combat.set
#  TNBC.ALLset



#Define gene.list
  #gene list by probesets
gene.list.probe <- metabList$PROBEID
isna <- is.na(gene.list.probe)
gene.list.probe <- gene.list.probe[!isna]
gene.list <- gene.list.probe

##gene list by symbol, if converting
#gene.list <- metabList$SYMBOL
#isna <- is.na(gene.list)
#gene.list <- gene.list[!isna]
##i <- gene.list %in% row.names(exprs(esetALL.filt))
##gene.list <- gene.list[i] #ensures eset contains all genes from gene.list

#create eset filtered for gene list by probesets
ii <- row.names(exprs(esetALL.filt)) %in% gene.list.probe
esetMetab <- esetALL.filt[ii,]
featureNames(esetMetab)
#label rows by symbol
#row.names(esetMetab) <- 
 # paste(c(as.character(fData(esetMetab)$symbol, as.character(fData(esetMetab)$Affy))
featureNames(esetMetab)
#row.names(esetALL.filt) <- fData(esetALL.filt)$symbol

i <- gene.list %in% row.names(esetMetab)
gene.list <- gene.list[i] #ensures eset contains all genes from gene.list
gsList <- gene.list <- featureNames(esetMetab)


######SURV CURVES TWO THRDS/MEDIAN METHODS#######
library(annotate)
library(hgu133a.db)
datALL <-exprs(esetMetab) #gene expression matrix
gene.list <- rownames(datALL)
dataset <- "TNBCset All"
ndata <- as.numeric(nrow(pData(esetMetab)))
pdf("TNBC.All.rescale.pdf",width=12.5,height=9)
for(the.gene in gene.list) { 
  gn.expr=datALL[the.gene,] 
  #res=surv.analysis((rescale(x=gn.expr, q=0.05)-0.5)*2, status=status,stime=stime,method="two_thirds")
  res=surv.analysis(gn.expr,status=status,stime=stime,method="two_thirds")
  mfit=res$mfit
  tt=res$test.diff
  pval=round(1-pchisq(tt$chisq,1),11)
  gn.s <- getSYMBOL(the.gene,"hgu133a")
  plot(mfit,conf.int=FALSE,col=c("red","blue"),lty=1:1,lwd=2.5,xlab="Years",ylab="Relapse-free survival (%)",
       main=paste("Kaplan-Meier: p-value=",pval), sub=paste("dataset:", dataset, " n= ", ndata, sep=""), cex.lab=1.25, cex.main=1.5,)
       legend("topright", c(paste("high ",gn.s,sep=""),paste("low ",gn.s,sep="")), lty=c(1), lwd=2.5, cex=1.5, col=c("red","blue"))
}
dev.off()

########


###COMMON GENE EXPRESSION FUNCTIONS AND SCRITPS###
#####

###################
#rescale fx
###################
rescale <- function(x, na.rm=FALSE, q=0.05) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return((x - 0.5) * 2)
}
ex <- dd$dsQ
ex <- ((rescale(x=ex, q=0.05) - 0.5) * 2)
#or
d <- ex
d <- exprs(tnbc.Allset)
d <- sweep(d,1, apply(d,1,median,na.rm=T))
exprs(esetALL) <- d
ex <- d

esetALL <- tnbc.Allset  #TNBC.combat.set
esetALL.filt <- nsFilter(esetALL, remove.dupEntrez=TRUE, filterByQuantile=TRUE,  var.cutof =0.2)
esetALL.filt<- esetALL.filt$eset

##############################
#ANNOTATE EXPRES SETFUNCTION
##############################

fAnnot <- function(expset) {
  x <- expset
  symbol <- unlist(mget(rownames(x),hgu133aSYMBOL))
  genename <-unlist(mget(rownames(x),hgu133aGENENAME))
  entrezid <- unlist(mget(rownames(x),hgu133aENTREZID))
  fData <- data.frame(Affy=rownames(x), symbol=symbol, genename=genename, entrezid=entrezid);
  fmetadata <- data.frame(labelDescription=c("Affymetrix accession", "Gene symbol", "Gene Name", "Entrez ID"),
                          row.names=c("Affy", "symbol", "genename", "entrezid"))
  
  fd <- new("AnnotatedDataFrame", data=fData, varMetadata = fmetadata)
  featureData(x)<- fd
  expsetann <- x
  return(expsetann)
}

tnbc.Allset <- fAnnot(tnbc.Allset)
TNBC.gpl96RMAonQCset <- fAnnot(TNBC.gpl96RMAonQCset)

########

#GCRMA TO EXPRESSIONSET
eIDs<- sampleNames(dd) <- toupper(sampleNames(dd))
#phenoAll <- phenoData(TNBC.gpl96.RMAset)
p.TNBC <- phenoTNBC[eIDs,]
#TNBC.gpl96QCset <- ExpressionSet(assayData=dd$ds, phenoData=p.gpl96, annotation="hgu133a.db")
tnbc.Allset <- ExpressionSet(assayData=dd$ds, phenoData=p.TNBC, annotation="hgu133a.db")
tnbc.Allsetgene <- ExpressionSet(assayData=dd$ds.gn.lvl, phenoData=p.TNBC, annotation="hgu133a.db")






# biocLite("genefilter")
# boxplot(tnbc.gpl96.RMAset)
# ex <- dd$ds
# ex <- exprs(tnbc.gpl96.RMAset)
# ex <- exprs(esetALL)
# ex.s <- scale(t(ex))
# # Ward Hierarchical Clustering
# d <- dist(ex, method = "euclidean") # distance matrix
# fit <- hclust(d, method="ward") 
# plot(fit) # display dendogram
# groups <- cutree(fit, k=5) # cut tree into 5 clusters
# # draw dendogram with red borders around the 5 clusters 
# rect.hclust(fit, k=5, border="red")

ex <- dd$dsQ
d <- ex
d <- exprs(tnbc.Allset)
d <- sweep(d,1, apply(d,1,median,na.rm=T))
exprs(esetALL) <- d
ex <- d
ii <- row.names(ex) %in% metabList$PROBEID
ex.m <- ex[ii,]



pdf("UnsupclusterAllnorm.pdf", width=32, height=10)
distance <- dist(t(ex),method="euclidian")
clusters <- hclust(distance, method="ward.D2")
par(cex=0.5)

plot(clusters, hang=-1)
dev.off()
clust <- as.dendrogram(clusters)

#####
d <- exprs(esetALL.filt)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
Index <- as.factor(pData(tnbc.Allset)$series_id)
ser <- levels(tnbc.Allset$series_id)
labcol <- ser[Index]
or <- order(ser[Index])
niceCols <- brewer.pal(10, "Spectral")
z <- niceCols[Index]
###HEATMAP####
pdf("Cluster.unsupdend.rescale.pdf",width=36,height=24)
heatmap.2(ex,
          #labCol=labcol,
          scale="row", 
          #dendrogram="none", 
          Colv=clust,
          col=my_palette, trace= "none", density="density",
          ColSideColors=z,
          key=FALSE,)

# par(lend = 1)
# legend("topleft", 
#        legend = ser, 
#        col = niceCols, 
#        lty= 1, lwd = 10, inset = -1.5)
dev.off()



##############################
#QC report
##############################
library("arrayQualityMetrics", lib.loc="/Users/ashleysdoane/Library/R/3.1/library")

arrayQualityMetrics(expressionset = tnbc.gpl96.RMAset,
                    outdir = "Report.tnbcRMA",
                    force = TRUE,
                    intgroup = c("series_id"))



##############################
#Annotation
##############################

library(annotate)
library(hgu133a.db)

###GENE convert from metabolic gene list to probeIDs######
m.genes <- read.table("m.genes.txt", stringsAsFactors=FALSE)
m.genes<- m.genes[,1]
library(org.Hs.eg.db)
# keys <- keys(org.Hs.eg.db, keytype = "ENTREZID"), n = 2)
keys <- keys(hgu133a.db)
#  columns <- c("PFAM", "GO", "SYMBOL")
#  select(org.Hs.eg.db, keys, columns, keytype = "PROBEID")
keytypes(hgu133a.db)
library("hgu133a.db")
dbc=hgu133a_dbconn()
mtgids <- findEGs(dbc, m.genes)
k <- mtgids$gene_id
metabList <- select(hgu133a.db, keys=k, columns = c("SYMBOL","PROBEID", "GENENAME"), keytype="ENTREZID")
save(metabList, file="metabList3.Rdata")
metabProbe <- metabList$PROBEID

