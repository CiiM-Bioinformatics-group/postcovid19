## Methylation 
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)
library(ggplot2)
library(wateRmelon)
library(future)
library(gtools)
library(matrixStats)
library(data.table)
library(MASS) 
library(sandwich) 
library(lmtest) 
library(parallel) 
library(R.utils)
library(data.table)
library(readxl)
library(FlowSorted.Blood.EPIC) 
library(DirichletReg) 
library(qqman)
library(tidyr)
library(viridis)
library(ggridges)

targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)
Pd <- pData(RGset)

############ sample QC ####################
agegender <- read_excel("group.xlsx",sheet=2,na="NA")
age <- agegender$age
gender <- agegender$gender
group <- agegender$group
source <- agegender$source
pd <- cbind(Pd, age, gender, group, source)
MSet.raw <- preprocessRaw(RGset) 
qc <- getQC(MSet.raw)
plotQC(qc)

# smaple filter
#gender checking
ratioSet <- ratioConvert(MSet.raw, what = "both", keepCN = TRUE)
gset <- mapToGenome(ratioSet)
predictedsex <- getSex(gset, cutoff = -2)
gset <- addSex(gset, predictedsex)
pd2 <- pData(gset)
plotSex(gset, id = pd2$Sample_Name)
documentedSex <-pd$gender
predictedSex <- predictedsex$predictedSex
gcheck <- function(docsex,predsex,samnam) {
  docsex <- as.character(docsex)
  sex.match <- identical(docsex,predsex)
  if (sex.match == FALSE) {
    sex.mismatch <- samnam[which(docsex != predsex)]
    write.csv(sex.mismatch, file="sexmismatch.csv") }
}
gcheck(documentedSex,predictedSex, Pd$Sample_Name)


det.p <- detectionP(RGset)
failed <- det.p > 0.01
failed.fraction <- colMeans(failed) 
failed.fraction <- data.frame(failed.fraction)
failed.fraction <- cbind(Pd$Sample_Name,failed.fraction)
failed.fraction <- data.frame(failed.fraction[failed.fraction[,2] > 0.01,])
fail.detp<-failed.fraction$Pd.Sample_Name 
write.csv(fail.detp, file="failed_samples_by_probe_frequency.csv")
fail.sample.index <- which(Pd$Sample_Name %in% fail.detp)
filtered_RGset <- RGset[,-fail.sample.index]
colnames(det.p) <- Pd$Sample_Name
barplot(colMeans(det.p), las=2, cex.names=0.8, ylim = c(0,0.01), ylab="Mean detection p-values")

## remove all failed samples from RGset
fail.sample<-unique(c(sex.mismatch,fail.detp))
fail.sample.index<-which(Pd$Sample_Name %in% fail.sample)
filtered_RGset <- RGset[,-fail.sample.index]
pd.flt <- pd[-fail.sample.index,]

############## probes filtering and normaliztion ##############
det.p <- detectionP(filtered_RGset)
bad.probes <- rowMeans(det.p > 0.01) > 0.1
table(bad.probes)
bad.probe.names.detP <- rownames(det.p[bad.probes,])

m.set <- preprocessRaw(filtered_RGset)
ratioSet <- ratioConvert(m.set, what = "both", keepCN = TRUE)
gm.set <- mapToGenome(ratioSet)

## remove the cross-reactive probes, polymorphic probes. 
## "13059_2016_1066_MOESM4_ESM.csv", "13059_2016_1066_MOESM5_ESM.csv", and "13059_2016_1066_MOESM1_ESM.csv" were download from the reference (Pidsley et al., 2016).
snpprobes1 <- read.csv(file=paste(cros, "13059_2016_1066_MOESM4_ESM.csv", sep="/"), sep=",", stringsAsFactors=FALSE)
keep2 <- !(featureNames(gm.set) %in% snpprobes1$PROBE[snpprobes1$EUR_AF>0.05])
bad.probe.names.pm1 <- rownames(gm.set[!keep2,])
table(keep2)

snpprobes2 <- read.csv(file=paste(cros, "13059_2016_1066_MOESM5_ESM.csv", sep="/"), sep=",", stringsAsFactors=FALSE)
keep3 <- !(featureNames(gm.set) %in% snpprobes2$PROBE[snpprobes2$EUR_AF>0.05])
bad.probe.names.pm2 <- rownames(gm.set[!keep3,])
table(keep3)

reactive.probes1 <- read.csv(file=paste(cros, "13059_2016_1066_MOESM1_ESM.csv", sep="/"), sep=",", stringsAsFactors=FALSE)
reactive.probes2 <- unlist(reactive.probes1)
keep1 <- !(featureNames(gm.set) %in% reactive.probes2)
bad.probe.names.cr <- rownames(gm.set[!keep1,])
table(keep1)

## Remove probes on X and Y chr
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
X <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% c("chrX")]
Y <- featureNames(gm.set) %in% ann850k$Name[ann850k$chr %in% c("chrY")]
probe.names.x <- rownames(gm.set[X,])
probe.names.y <- rownames(gm.set[Y,])

## Normalization. we did normalization before removing the probes on X and Y chr. and then remove them. 
remove.probenoxy <- unique(c(bad.probe.names.detP, bad.probe.names.pm1, bad.probe.names.pm2, bad.probe.names.cr))
MSetraw <- preprocessRaw(filtered_RGset)
MSetrp <- MSetraw[!(rownames(MSetraw) %in% remove.probenoxy),]
gmsetrp <- mapToGenome(MSetrp)
MSet.sq <- preprocessQuantile(gmsetrp)
remove.probexy <- unique(c(probe.names.x,probe.names.y))
m.set.flt <- MSet.sq[!(rownames(MSet.sq) %in% remove.probexy),]
save(m.set.flt, file="normalization.Rdata")

## mds plot ##
Mvalue <- getM(m.set.flt)
mdsPlot(Mvalue, numPositions = 1000, sampGroups = pd.flt$group, main="M values")

# density plot before and after normalization 
par(mfrow=c(1,2))
densityPlot(MSet.raw, sampGroups=pd$group, main="Raw", legend=FALSE) 
legend("topleft", legend = levels(factor(pd$group)), text.col=brewer.pal(8,"Dark2")) 
densityPlot(getBeta(m.set.flt), sampGroups=pd$group, main="Normalized", legend=FALSE)
legend("topleft", legend = levels(factor(pd$group)), text.col=brewer.pal(8,"Dark2"))

################### estimate cell proportion  #############
cells <- estimateCellCounts2(filtered_RGset, sex = pd.flt$gender, referencePlatform = "IlluminaHumanMethylationEPIC",cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono"),compositeCellType = "Blood")
cells <- cells$counts
cells[cells < 0] <- 0
cells.box <- data.table(name = gsub("_.*", "", pd.flt$Sample_Name),
                        group = pd.flt$group,
                        cells)
cells2.box <- data.table::melt(cells.box, id = 1:2, measure = colnames(cells.box)[-c(1,2)])
my.dirichreg <- function(celltype = "CD4T", group1 = "healthy-control", group2 = "post-covid19", data = cells, targets = pd.flt){
  data <- data.frame(CpG = pd.flt$Sample_Name,
                     group = pd.flt$group,
                     CD4T = data[,celltype],
                     Other = 1 - data[,celltype])
  data <- subset(data, group %in% c("healthy control", "post-covid19"))
  data$Smp <- DR_data(data[,3:4])
  res <- DirichReg(Smp~group, data, model = "alternative", base = 1, verbosity = 0, control = list(iterlim = 1000)) 
  x <- summary(res)
  result <- data.frame(CellType = celltype, Comparison = paste0(healthy-control, "versus", post-covid19), Z = x$coef.mat[2, 3], P = x$coef.mat[2, 4])}

results <- mclapply(colnames(cells), my.dirichreg, group1 = "healthy-control", group2 = "post-covid19", data = cells, targets = pd3)

## figure 4A 
ggplot(cells2.box, aes(x = variable, y = value)) +
   geom_boxplot(aes(fill = group)) +
   facet_grid(.~variable) +
   geom_jitter(shape=16, size = 0.5, position=position_jitter(0.2)) +
   ggtitle(label = "Cell Proportions ") +
   xlab("group") +
   ylab("Percentage") +
   theme_light() +
   # geom_hline(yintercept = 0, shape = group) +
   theme(axis.text = element_text(size = 12),
         axis.title = element_text(size = 16))

area <- ggplot(cells2.box, aes(x = variable, y = value)) +
   geom_boxplot(aes(fill = group)) +
   ggtitle(label = "Cell Proportions ") +
   xlab("cell type") +
   ylab("Percentage") +
   theme_light() +
   # geom_hline(yintercept = 0, shape = group) +
   theme(axis.text = element_text(size = 12),
         axis.title = element_text(size = 16))
area + scale_fill_manual(values = c("steelblue", "tomato")) + geom_jitter(aes(group=factor(group)),position=position_dodge(width=0.8))

################### EWAS #################################################
# we use limma package
group <- pd.flt$group
age <- as.numeric(pd.flt$age)
gender <- as.factor(pd.flt$gender) 

de <- model.matrix(~group + age + gender)
fit <- lmFit(Mvalue, de)
fit <- eBayes(fit)
Rresult2 <- topTable(fit, adjust.method = 'BH', coef=2, sort.by = "p", number = 794189)

## manhattan plot, supplementary figure 5A
annotdata<-ann850k[c( "chr","pos","UCSC_RefGene_Name","Relation_to_Island")]
result2 <- data.frame(Rresult2)
Ml<-data.frame(merge(result2,annotdata,by="row.names"))
Ml$BP<-as.integer(Ml$pos)
Ml$SNP<-as.character(Ml$Row.names)
Ml$CHR<-as.numeric(gsub("[^0-9]", "", (Ml$chr)))
Ml$P<-Ml$P.Value
Ml1<-Ml[,c("CHR","BP","P","SNP")]
manhattan(Ml1, main = "Manhattan Plot", cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = 4)
## supplementary figure 5B
qq(Ml1$P)
lambda = median(qchisq(as.numeric(as.character(Rresult2$adj.P.Val)),df=1,lower.tail = F),na.rm=T)/qchisq(0.5,1)
## supplementary figure 5C 
ggplot(Rresult2, aes(x = logFC, y = -log(P.Value, 10)))

############################# correlation analysis ##################
Rresult2$site <- rownames(Rresult2)
top50 <- Rresult2_1[1:50,]
top75 <- Rresult2_1[1:75,]
top100 <- Rresult2_2[1:100,]
top150 <- Rresult2[1:150,]

## figure 4B left
## we can get the number of up/down or hypo/hypermethylated sites from the top50, top75, top100, top150 object
x <- rep(c("top50","top75","top100","top150"), each = 2)
y <- rep(c("up","down"), times= 4)
z <- c(17,33,28,47,34,66,56,94)
df1 <- data.frame(x=x, y=y,z=z)
df1$x = factor(x,levels = c("top50","top75","top100","top150"))

ggplot(data = df1, mapping = aes(x = x, y = z, fill = y))+geom_bar(stat = "identity", position = "dodge")+scale_x_discrete(breaks = c("top50","top75","top100","top150"), labels = c("top50","top75","top100","top150"))+scale_fill_manual(values = c("steelblue","tomato"),name="Post-COVID-19 vs Control")+theme_classic()+ylab("Number of differentially methylated sites")+
  ylim(0,100)+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"))
## figure 4B right
x <- rep(c("top50","top75","top100","top150"), each = 2)
y <- rep(c("hypermethylation","hypomethylation"), times= 4)
z <- c(29,21,42,33,53,47,83,67)
df1 <- data.frame(x=x, y=y,z=z)
df1$x = factor(x,levels = c("top50","top75","top100","top150"))
ggplot(data = df1, mapping = aes(x = x, y = z, fill = y))+geom_bar(stat = "identity", position = "dodge")+scale_x_discrete(breaks = c("top50","top75","top100","top150"), labels = c("top50","top75","top100","top150"))+scale_fill_manual(values = c("steelblue","tomato"),name="Post-COVID-19 vs Control")+theme_classic()+ylab("Number of differentially methylated sites")+
  ylim(0,100)+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"))

## figure 4C, 4D, supplementary figure 5D, supplementary 5E
up150 <- subset(top150, logFC > 0) 
down150 <- subset(top150, logFC < 0)
pd <- data.frame(pd.flt)
beta <- getBeta(m.set.flt)
beta <- data.frame(beta)
B.up <- beta[which(rownames(beta) %in% rownames(up150)),]
B.down <- beta[which(rownames(beta) %in% rownames(down150)),]

colnames(B.all) <- pd.flt$Sample_Name
cells <- t(cells)
colnames(cells) <- pd.flt$Sample_Name
sc2 <- cells
rownames(sc2) <- c("CD8", "CD4", "NK", "B", "Mono", "Neu")
cell<-c()
site<-c()
cor_r<-c()
pvalue<-c()

for (i in 1:nrow(sc2)){
  for (r in 1:nrow(B.all)){
    c1=rownames(sc2)[i]
    s2=rownames(B.all)[r]
    c_s=cor(as.numeric(sc2[i,]),as.numeric(B.all[r,]),method="spearman")
    p=cor.test(as.numeric(sc2[i,]),as.numeric(B.all[r,]),method ="spearman")[[3]]
    cell=c(cell,c1)
    site=c(site,s2)
    cor_r=c(cor_r,c_s)
    pvalue=c(pvalue,p)
  }
}
data_cor<-data.frame(cell,site,cor_r,pvalue)
ggplot(data_cor, aes(x = cor_r, y = cell, fill = cell)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")+
  xlab("Correlation coefficient")+
  ylab("Cell type")+
  labs(title = "The correlation between estimated cell proportion and CpG sites", subtitle = "up-regulated CpG sites in top 150 sites")

####################################################################################################
