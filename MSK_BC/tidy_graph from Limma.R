#from Limma result
library(dplyr)
library(tidyr)
pDat <- tbl_df(pData(MSK.set))

##fix eset 
pData(MSK.set)$gsm <-  as.factor(pData(MSK.set)$gsm)
pDat <- pData(MSK.set)
###

y = exprs(MSK.set)
exp = data.frame(y[rownames(tab),])
exp <- tbl_df(exp)


exp_ann_all <- exp %>%
  mutate(Affy = tab[rownames(exp),]$Affy) %>%
  inner_join(tab) %>%
  mutate(sym_affy = paste(symbol, Affy, sep= '_'))

pDat <- pDat %>%
  select(class, gsm)

expData <- exp_ann_all %>%
  select(sym_affy, starts_with("GSM")) %>%
  gather(gsm, exp, starts_with("GSM")) %>%
  left_join(pDat) %>%
  spread(sym_affy, exp)

annData <- exp_ann_all %>%
  select(Affy, symbol, genename, sym_affy, P.Value, adj.P.Val) %>%
  unique()

MSKtidy <- tbl_df(expData)
MSKtidy
save(MSKtidy, file="MSKtidy.Rdata")
##ggood##, now graph##



library(ggplot2)
library(RColorBrewer)

classd <- MSKtidy$class

MSKtidy[(gs[c(i,2)])]
gs <- annData$sym_affy[1:10]
gs
annData[,2]
gss <- names(MSKtidy)[c(-1,-2)]
#works
ggplot(MSKtidy, aes_string(x = "class", annData[1,4])) + geom_boxplot(width=0.5, outlier.size = 0) + 
  geom_point(aes(colour=class), alpha=0.7, size=2, position = position_jitter(w = 0.1, h = 0.1)) +
  scale_color_brewer(palette="Set1") + ggtitle(paste(annData[1,2], "\n","adj p val=", (signif(annData[1,6],2))))

                                               
paste(annData[1,2], "\n", annData[1,6], digits = 2)

paste(signif(annData[1,6],2))

x <- annData[1,6]
str(x)
paste(x, digits = 3)


annData$symbol <- as.character(annData$symbol)


pdf("gs.pdf")
print(ggplot(MSKtidy, aes_string(x = "class", gs[1])) + geom_boxplot())
dev.off()

boxp

theme_set(theme_bw(base_family="Helvetica", base_size=10))

###working ggplot graph######
pdf("gs.pdf")
for (i in 1:length(annData$sym_affy)) {
  print(ggplot(MSKtidy, aes_string(x = "class", annData[i,4])) + geom_boxplot(width=0.5, outlier.size = 0) + geom_point(aes(colour=class), alpha=0.7, size=2, position = position_jitter(w = 0.1, h = 0.1)) +scale_color_brewer(palette="Set1") + ggtitle(annData[i,2]))}
dev.off()
############

pdf("gsx.pdf")
for (i in 1:length(gs)) {
  print(ggplot(MSKtidy, aes_string(x = "class", gs[i])) + geom_boxplot() + geom_point(position = position_jitter(w = 0.1, h = 0.1)))}
dev.off()



x <- MSKtidy[(gs[c(i,2)])]
boxplotplus2(x)



for (i in 1:length(gs)) {
  ggsave(filename="gs.pdf", MSKtidy, aes_string(x = "class", gs[i])) + geom_boxplot())}



boxplot(gs[1] ~ class, data= MSKtidy)
class <- MSKtidy$class
boxplot(x ~class)

boxplot(MSKtidy[gs[1]] ~MSKtidy$class)

gs <- annData$sym_affy[1:10]
gs[1]
as.symbol(gs[1])
for (i in 1:length(gs)) {
  print(MSKtidy[(gs[c(i,2)])])}

for (i in 1:length(gs)) {
  print(MSKtidy[(gs[c(2,i)])])}

MSKtidy[(2:4)]
MSKtidy[,c(paste(c(gs[1])))]

c(paste(c(gs[1])))

       
MSKtidy[which()]

MSKtidy[c(2,4)]
for (i in 1:length(gs)) {
  print(MSKtidy[(gs[i])])}

for (i in 1:length(gs)) {
  print(ggplot(MSKtidy, aes(class, gs[i])) + geom_boxplot())}

names(MSKtidy)[5]

names(MSKtidy[(gs[1])])
MSKtidy[(gs[1])]

ggplot(MSKtidy, aes(class, MSKtidy[(gs[1])])) + 