library("affy")
library("gcrma")
library("hgu133plus2.db")
library("annotate")
library("gplots")
library("GEOquery")

source("http://bioconductor.org/biocLite.R")
biocLite("oligo")
biocLite("lumi")
biocLite("affy")
biocLite("hgu133a2frmavecs")
biocLite("u133aaofav2cdf")
biocLite("u133aaofav2probe")
biocLite("maPredictDSC")
#DL GEO data#

geo.number="GSE3494"
geo.number="GSE31519"
[]
geo.number="GSE25066"
geo.number="GSE1456"
geo.number="GSE2990"

getGEOSuppFiles(geo.number)

process.raw.affy.data(w,w)

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

w<-getwd()
w
process.raw.affy.data(w,w)



phenoALL<- read.csv("/Users/ashleysdoane/R Bioconductor Data/KoutcherData/phenoAll.csv", header=TRUE)
head(phenoALL)
keep <- c("VDX", "UPP", "UNT", "CAL", "TRANSBIG", "MAINZ", "EMC2", "DFHCC4", "DFHCC")
keep <- as.factor(keep)
n <- which(phenoALL$dataset %in% "UNT")
phenoF <- phenoDATA[n, ]
phenoF <- data.frame(phenoF)
names(phenoF)

#Merge JNCI paper data and GEO pheno data###
phenoDATA$er <- as.factor(phenoDATA$er)

e <- "STK" #Data set name
e <- which(phenoDATA$dataset %in% e)
phenoSTK <- phenoDATA[e, ]

#text manipulation- remove letters from number ids, ect
STK <- read.table("GSE1456_STK.txt", header=TRUE)
STK$id2 <- as.character(STK$id)
STK$id2 <- gsub("^*[A-z]", "", STK$id2) #remove letters at beginning
STK$id2 <- gsub("*[A-z]$", "", STK$id2) #remove letters at end
STK$id2
STK$id3 <- as.numeric(STK$id2)
STK$id3
STK$id <- as.character(STK$id3)
STK$id2 <- NULL
STK$id3 <- NULL

phenoDFHCC$id<- gsub("^.*?_","",phenoMAINZx$id)
#phenoMAINZx$id <- as.numeric(phenoMAINZx$id)
#write.table(phenoMAINZx, "phenoMainz.txt")

e <- "UPP" #Data set name
e <- which(phenoDATA$dataset %in% e)
phenoUPP <- phenoDATA[e, ]
UPP <- read.table("UPPGSE.txt", header=TRUE, stringsAsFactors=FALSE) #GEO file with geno_accn
UPP$geno_accn <- toupper(UPP$geno_accn)
UPP$affy <- NULL
#MAINZ <- read.table("MainzGSE.txt", header=TRUE, stringsAsFactors=FALSE) #GEO file with geno_accn
phenoUPP$id <- as.character(phenoUPP$id)


UPPp <- merge(UPP, phenoUPP)
UPPp <- UPPp[ ,c(2,1,3:28)] #set geno_accn as first column
rownames(UPPp) <- UPPp[ ,1] #set geno_accn as rownames
UPPp <- UPPp[order(UPPp$geno_accn),] #sort on #geno_accn
p<- colnames(UPPp)
p
View(UPPp)
#write.table(UPPp, "UPPp.txt")
#UPPp <- read.table("UPPp.txt")
varInfo <- data.frame(labelDescription= (p), row.names=(p))
vmd <- new("AnnotatedDataFrame", data=UPPp, varMetadata = varInfo) #pData for expressionset

sampleNames(vmd) #should be geno_accn and match columns of dd$ds



##limit dd$ds to specified columns (samples)###
s <- rownames(UPPp)
x <- match(s , colnames(dd$ds))
exprs <- dd$ds[, x]

UPPSet <- new("ExpressionSet", exprs=exprs,
              phenoData=vmd,
              annotation="hgu133a.db")

save(UPPSet, file="UPPSetv2.Rdata",compress=T)

#er status data for UPPSet

UPPSet[$phenoData]

#######UNT##
UNT <- read.table("UNTGSE.txt", header=TRUE, stringsAsFactors=FALSE) #GEO file with geno_accn
fix(UNT)

un$samplename <- as.character(un$samplename)
unt <- merge(un, phenoUNT)

e <- "UNT" #Data set name
e <- which(phenoDATA$dataset %in% e)
phenoUNT <- phenoDATA[e, ]
phenoUNT <- subset.data.frame (phenoDATA, dataset=="UNT")
phenoUN == phenoUNT #should be true for methods to be same

?expressionset

phenoData(UPPSet)$er
which(UNT[ ,()]  
UNTp <- merge(UNT, phenoUNT)
UNTp <- UNTp[ ,c(2,1,3:28)] #set geno_accn as first column
#low <-  tolower(UNTp$geno_accn)
#UNTp$geno_accn <-  low
rownames(UNTp) <- UNTp[ ,1] #set geno_accn as rownames
UNTp <- UNTp[order(UNTp$geno_accn),] #sort on #geno_accn
p<- colnames(UNTp)
p
View(UNTp)
write.table(UNTp, "UNTp.txt")
#UPPp <- read.table("UPPp.txt")
varInfo <- data.frame(labelDescription= (p), row.names=(p))
vmd <- new("AnnotatedDataFrame", data=UNTp, varMetadata = varInfo) #pData for expressionset

s <- sampleNames(vmd) #should be geno_accn and match columns of dd$ds

ddds <- as.data.frame(dd$ds)
rownames(ddds)
match(s, colnames(ddds))
colnames(ddds)
 <- sampleNames(dd)
y <- sampleNames(vmd)


x <- match(s , colnames(dd$ds))
exprs <- dd$ds[, x]



UNTSet <- new("ExpressionSet", exprs=exprs,
              phenoData=vmd,
              annotation="hgu133a.db")

save(UNTSet, file="UNTSet.Rdata",compress=T)


###STK SET#####

#text manipulation- remove letters from number ids, ect
STK <- read.table("GSE1456_STK2.txt", header=TRUE)
STK$id <- as.character(STK$ID_REF)
STK$id <- gsub("^*[A-z]", "", STK$id) #remove letters at beginning
STK$id <- gsub("*[A-z]*$", "", STK$id) #remove letters at end
STK$id
STK$id <- as.numeric(STK$id)
STK$id3
STK$id <- as.character(STK$id)
STK$id2 <- NULL
STK$id3 <- NULL

STK$id <- as.numeric(STK$id)
STK <- STK[with(STK, order(id)), ]
STK$id
colnames(STK) <- toupper(colnames(STK))
phenoSTK$id <- as.numeric(phenoSTK$id)
phenoSTK$id
phenoSTK <- phenoSTK[with(phenoSTK, order(id)), ]

DAYS <- with(STK, SURV_RELAPSE * 365)
DAYS<- as.integer(DAYS)
STK$DAYS <- DAYS
comp <- data.frame(tot$DAYS, tot$t.rfs)
print(comp)


STKp <- tot[,c(1:27,36)]
STKp$geno_accn <- STKp$GENO_ACCN
STKp <- STKp[, c(29, 1:27)]

rownames(STKp) <- STKp[ ,1] #set geno_accn as rownames
STKp <- STKp[order(STKp$geno_accn),] #sort on #geno_accn
p<- colnames(STKp)
p

p
write.table(STKp, "STKp.txt")
varInfo <- data.frame(labelDescription= (p), row.names=(p))
vmd <- new("AnnotatedDataFrame", data=STKp, varMetadata = varInfo) #pData for expressionset

STKset <- new("ExpressionSet", exprs=dd$ds,
              phenoData=vmd,
              annotation="hgu133a.db")

save(STKset, file="STKset.Rdata",compress=T)




###VDX fRMA GSE2034 and GSE5327

p <- pData(VDTMset)

vdx<- p[which(p$dataset == "VDX"),] #not all vdx sampless are in VDTMset (no GSE5327)
phenoVDX <- phenoDATA[which(phenoDATA$dataset == "VDX"),]

VDX <- read.table("VDXGSE.txt", head=TRUE)

VDXp <- merge(VDX, phenoVDX)
VDXp <- VDXp[ ,c(2,1,3:28)] #set geno_accn as first column
rownames(VDXp) <- VDXp[ ,1] #set geno_accn as rownames
VDXp <- VDXp[order(VDXp$geno_accn),] #sort on #geno_accn

p<- colnames(VDXp)
varInfo <- data.frame(labelDescription= (p), row.names=(p))
vmd <- new("AnnotatedDataFrame", data=VDXp, varMetadata = varInfo) #pData for expressionset

VDXset <- new("ExpressionSet", exprs=dd$ds,
              phenoData=vmd,
              annotation="hgu133a.db")

save(VDXset, file="VDXset2.Rdata",compress=T)



##fRMA Data Sets####

phenoEMC2 <- p[which(p$dataset == "EMC2"),]
phenoEMC2 <- phenoEMC2[order(phenoEMC2$geno_accn),] #sort on #geno_accn
#or from existing
phenoEMC2 <- pData(EMC2Set)

p<- colnames(phenoEMC2)
varInfo <- data.frame(labelDescription= (p), row.names=(p))
vmd <- new("AnnotatedDataFrame", data=phenoEMC2, varMetadata = varInfo) #pData for expressionset

##limit dd$ds to specified columns (samples)###
s <- rownames(phenoEMC2)
x <- match(s , colnames(dd$ds))
exprs <- dd$ds[, x]

EMC2set <- new("ExpressionSet", exprs=exprs,
              phenoData=vmd,
              annotation="hgu133plus2.db")

save(EMC2set, file="EMC2setv2.Rdata",compress=T)


####

i <- grep("GSE2", gsm2gse$series_id)
gsm2gse[i,2]

gsm2gse[(grep("GSE2", gsm2gse$series_id)), 2]




#####GSE21653 Work####

GSE21653 <- read.table("phenoGSE21653.txt", header=TRUE)

GSE21653<- GSE21653[order(GSE21653$geno_accn),]

GSE21653$er.status[GSE21653$er.status == 1] <- "P" 
GSE21653$er.status[GSE21653$er.status == 0] <- "N" 
GSE21653$er.status <- factor(GSE21653$er.status)

i <- complete.cases(GSE21653$er.status)
GSE21653[!i,]  #OK 3 cases NA

GSE21653$her2.status[GSE21653$he2r.status == 1] <- "P" 
GSE21653$her2.status[GSE21653$her2.status == 0] <- "N" 
GSE21653$her2.status <- factor(GSE21653$er.status)

EGSE21653 <- GSE21653[GSE21653$er.status == "P",]

ERP <- GSE21653 %in% 
ERP <-  GSE21653[GSE21653$er.status == "P",]
ERN <-  GSE21653[GSE21653$er.status != "P",]

TN <- ERN[ERN$her2.status != "P",]

i <- complete.cases(TN$geno_accn)

TN <- TN[i,]
TN

i <- complete.cases(TN$er.status)
TN[i,]

TN.GSE21653 <- TN
save(TN.GSE21653, file="pData.TN.GSE21653.Rdata")

GSM <- as.character(TN.GSE21653$geno_accn)

GSM.cel <- paste(GSM, ".CEL", sep="")
GSM.cel

files <- list.files("/Volumes/AD/Data/GSE.Data/GSE21653.cel.files/GSE21653_RAW/")
files
filesn <- gsub("_[012].*CEL", ".CEL", files)

i <- filesn %in% GSM.cel
TNfiles <- files[i]
TNfiles
for(.cel in TNfiles) {
  A = paste("/Volumes/AD/Data/GSE.Data/GSE21653.cel.files/GSE21653_RAW/", .cel, sep="")
  B = paste("/Volumes/AD/Data/GSE.Data/GSE21653.cel.files/TN.cel.files/", basename(.cel), sep="")
  file.copy(A,B)
}
setwd("/Volumes/AD/Data/GSE.Data/GSE21653.cel.files/TN.cel.files/")
w <- getwd()
files

for(.cel in files) {
  A= paste()
  
}
grep("GSM.*CEL", files)
list.files(w)

list.files()
file.rename(list.files(pattern="GSM.*CEL", paste0("water_", 1:700))

#rename files
  A = list.files(w)
  B = gsub("_[012].*CEL", ".CEL", A)
  file.rename(A, B)

file.rename
files <- list.celfiles("/Volumes/AD/Data/GSE.Data/GSE21653.cel.files/TN.cel.files/")
list.celfiles("/Volumes/AD/Data/GSE.Data/GSE21653.cel.files/GSE21653_RAW/")
files



#process fRMA
w <- getwd()
cel.dir <- "/Volumes/AD/Data/GSE.Data/GSE21653.cel.files/TN.cel.files/"
geo.number = "GSE21653"

frma.process.raw.affy.data(cel.dir, w)


#build eset

i <- rownames(ERN) %in% rownames(dd$ds)

ERN <- ERN[i,]
p<- colnames(ERN)
varInfo <- data.frame(labelDescription= (p), row.names=(p))
vmd <- new("AnnotatedDataFrame", data=ERN, varMetadata = varInfo)
sampleNames(vmd) 
rownames(exprsGSE31519surv)

TNBCsurVset <- new("ExpressionSet", exprs=t(exprsGSE31519surv),
                   phenoData=vmd,
                   annotation="hgu133a.db")

save(TNBCsurVset, file="TNBCsurVset.Rdata",compress=T)]


rownames(dd$ds)

i <-
rownames(TN.GSE21653) <- TN.GSE21653[,2]
i <- rownames(TN.GSE21653) %in% rownames(dd$ds)
i
rownames(TN.GSE21653)

ERN <- ERN[i,]
p<- colnames(TN.GSE21653)
varInfo <- data.frame(labelDescription= (p), row.names=(p))
vmd <- new("AnnotatedDataFrame", data=TN.GSE21653, varMetadata = varInfo)
pTNBCsurv
sampleNames(vmd) 
s <- rownames(TN.GSE21653)
x <- match(s , colnames(dd$ds))
exprs <- dd$ds[, x]

TN.GSE21653set <- new("ExpressionSet", exprs=exprs,
                   phenoData=vmd,
                   annotation="hgu133plus2.db")

save(TN.GSE21653set, file="TN.GSE21653set.Rdata",compress=T)

TN.GSE21653set


#####fRMA on GSE31519 by chip##

#GSE31519.gpl96
#AND
#GSE31519.gpl570 #includes array express samples
#GSE31519.AEgpl96

cel.file.list <- list.celfiles(cel.dir)
AE.affybatch <- ReadAffy(celfile.path="/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AE/")
#AE.gpl96 <- convertPlatform(AE.affybatch, "hgu133a"), did nor work, use convertCel below
table(ptn$Array_type)

geo.number="GSE31519.AE.gpl96"
cel.dir <- "/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AEU.133A/"
w <- getwd()
frma.process.raw.affy.data(cel.dir, w)

#readCel("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AE/LBL_POP_W405609.H02.CEL")
#readCel("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/TNBC.cel.files.gpl96/GSM107076.CEL")

cel.file.list<- list.celfiles("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AE.backup/")
for(file in cel.file.list)
  {
  outFile = file.path("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AEU.133A", (file))
  filename = file.path("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AE.backup/", (file))
  #newChipType = "HG-U133A"
  convertCel(filename, outFile, newChipType="HG-U133A")
}
  
  file.path("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AEU.133A/", basename(cel[1]))
  
  convertCel("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AE/LBL_POP_W405609.H02.CEL", "/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AEU.133A/LBL_POP_W405609.H02.CEL", newChipType = "HG-U133A")


frma.process.raw.affy.data(cel.dir, w)



#rename cel files
cel.dir <- "/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AEU.133A"
cel.file.list <- list.celfiles(cel.dir)
cel.file.list <- cel.file.list[order(cel.file.list)]
cel.file.list
AE.key <- TNBC.AEkey[,c("Source.Name", "Array.Data.File")]

AE.key <- AE.key[order(AE.key[,1]),]

AE.key[,2] %in% cel.file.list

AE.key <- AE.key[which(AE.key$Array.Data.File == cel.file.list),]


for(file in AE.key[,2])
{
  outFile = file.path("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AE.rename/", (file))
  filename = file.path("/Volumes/AD/Data/GSE.Data/TNBC.cel.files/AE.backup/", (file))
  file.rename(filename, outFile)
}
file.rename