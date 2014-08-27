cel.dir <- getwd()
output.dir <- "/Users/ashleysdoane/R_Bioc.AD/MSK_BC/AR_Cells"

geo.number <- "GSE28789"

RMA.process.raw.affy.data<-function(cel.dir, output.dir){    
  library(affy)
  
  ## set output file names based on date
  rdata.file = paste(geo.number,"RMA.Rdata",sep="_")
  
  cat(c("Reading .CEL files in ",cel.dir,"..."))
  setwd(cel.dir)
  
  ## read cel files in current working directory
  cel.files=ReadAffy(celfile.path=".")
  cat(c("Done.","\n"))
  
  ## run GCRMA on raw data
  cat(c("Running RMA...","\n"))
  start=proc.time()
  dd=rma(cel.files)
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


library(dplyr)
p453 <- read.table("453.txt", header = TRUE, stringsAsFactors = FALSE)

edit(p453)
p453 %>%
  mutate(cond = tx + t)

x <- mapply(c, as.character(p453$tx), as.character(p453$t))
x
(c, p453$tx, p453$t)

RMA.process.raw.affy.data(cel.dir, output.dir)


sampleNames(dd) <- strsplit2(sampleNames(dd), "_")[,1]

sampleNames(dd) %in% pheno453$gsm
sampleNames(pheno453)
sampleNames(pheno453) <- pheno453$gsm

##limit dd$ds to specified columns (samples)###

l <- c("V", "D4", "D16", "D48")

p453$desc <- strsplit2(p453$desc, "_")[,2]
p453$desc <- factor(p453$desc, levels = l)
p453$desc
rownames(p453) <- p453$gsm
p <- colnames(p453)
varInfo <- data.frame(labelDescription= (p), row.names=(p))
pheno453 <- new("AnnotatedDataFrame", data=p453, varMetadata = varInfo) #pData for expressionset


s <- rownames(pheno453)
x <- match(s , colnames(dd$ds))
exprs <- dd$ds[, x]
exprs[1,]
M453.set <- new("ExpressionSet", exprs=exprs,
              phenoData=pheno453,
              annotation="hgu133plus2.db")

fAnnot <- function(expset) {
  x <- expset
  symbol <- unlist(mget(rownames(x),hgu133plus2SYMBOL))
  genename <-unlist(mget(rownames(x),hgu133plus2GENENAME))
  entrezid <- unlist(mget(rownames(x),hgu133plus2ENTREZID))
  fData <- data.frame(Affy=rownames(x), symbol=symbol, genename=genename, entrezid=entrezid);
  fmetadata <- data.frame(labelDescription=c("Affymetrix accession", "Gene symbol", "Gene Name", "Entrez ID"),
                          row.names=c("Affy", "symbol", "genename", "entrezid"))
  
  fd <- new("AnnotatedDataFrame", data=fData, varMetadata = fmetadata)
  featureData(x)<- fd
  expsetann <- x
  return(expsetann)
}
M453.set <- fAnnot(M453.set)
save(M453.set, file="M453.set.Rdata",compress=T)

pData(M453.set)


lev <- levels(M453.set$desc)
f <- M453.set$desc
design <- model.matrix(~0+f)
colnames(design) <- lev
fit <- lmFit(M453.set, design)

cont.D <- makeContrasts("D4-V", "D16-D4", "D48-D16", levels=design)
fit2 <- contrasts.fit(fit, cont.D)
fit2 <- eBayes(fit2)
tab453 <- topTableF(fit2, n=100, adjust="BH")

save(tab453, file="453DHT.Rdata")


