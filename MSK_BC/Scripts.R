###FILES###

load("~/R_Bioc.AD/MSK_BC/MSKset.Rdata")


library(dplyr)


t(exp)


###scripts##
############

###ERP A B analysis (3 way) LIMMA ####
###return limma results as 'tab' ######

library(limma)
f <- factor(MSK.set$class, levels=c("P","A","B"))
design <- model.matrix(~0+f)

colnames(design) <- c("P","A","B")
fit <- lmFit(MSK.set, design)
contrast.matrix <- makeContrasts(A-P, B-A, B-P, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)
venn <- vennDiagram(results)
tab = topTableF(fit2, number=100)

save(tab, file="tab.Rdata")

#setup df dplyr#######




y = exprs(MSK.set)

###After limma => define gene exp matrix of limma genes;
### Include grouping variables, using dplyr;
#####

pDat <- pData(MSK.set)
pDat <- tbl_df(pDat)
save(pDat, file="pDat.Rdata")
gene = as.character(tab[,2])
affy= as.character(tab[,1])

exp = data.frame(y[rownames(tab),])

#using indexing
exp$Affy <- row.names(exp)
for(id in 1:nrow(tab)){
  exp$symbol[exp$Affy %in% tab$Affy[id]] <- tab$symbol[id]
  exp$genename[exp$Affy %in% tab$Affy[id]] <- tab$genename[id]
}

exp
#or using indexing in dplyr fx:  dplyer mutate
annotate_exp <- function(exp, tab){
  exp_annot <- exp %>%
    mutate(
      Affy = tab[rownames(exp),]$Affy,
      symbol = as.character(tab[rownames(exp),]$symbol),
      genename = as.character(tab[rownames(exp),]$genename))
  return(exp_annot)
}

exp_a2 <- annotate_exp(exp, tab)

#works, but symbol and genename as char not factors
exp_anx <- exp %>%
  mutate(
    Affy = tab[rownames(exp),]$Affy,
    symbol = as.character(tab[rownames(exp),]$symbol),
    genename = as.character(tab[rownames(exp),]$genename))




#using dplyr with mutate and inner_join:##########
###WORKING##
###
exp = data.frame(y[rownames(tab),])
exp <- tbl_df(exp)

exp_ann_all <- exp %>%
  mutate(Affy = tab[rownames(exp),]$Affy) %>%
  inner_join(tab) %>%
  mutate(sym_affy = paste(symbol, Affy, sep= '_'))
save(exp_ann_all, file="exp_ann_all.Rdata")

exp_a <- exp_ann_all %>%
  select(sym_affy, starts_with("GSM"))
rownames(exp_a) <- exp_a$sym_affy
exp_a$sym_affy <- NULL

#transpose
texp = as.data.frame(t(exp_a))
texp2 = tbl_df(data.frame(t(exp)))
texp$gsm <- rownames(texp)
head(texp)

#or with reshape2 (working?)
library(reshape2)

exp_a2 <- exp_ann_all %>%
  select(symbol, Affy, genename, starts_with("GSM"))
x <- tbl_df(dcast(melt(exp_a2), variable~symbol + Affy))


names(x)[1] <- "gsm"
texp <- x
texp
#texp = as.data.frame(t(exp_a))
#texp = tbl_df(data.frame(t(exp_a)))
#texp

###########
########

#texp = as.data.frame(t(exp))

#make this a match or index 
colnames(texp) 
= gene

texp <- tbl_df(texp)



head(texp)



#texp$gsm <- row.names(texp)
#texp <- select(texp, !ids)
#texp$ids <- NULL
texp <- tbl_df(texp)


pDat <- pData(MSK.set)

pDat$gsm %in% texp$gsm

MSKtb <- pDat %>% inner_join(texp)
save(MSKtb, file = "MSKtbl.Rdata")

#byclass <- MSKtb %>% group_by(ER)
res <- MSKtb %>% select(-ids, -ER, -gsm) %>% group_by(class)

head(res)





