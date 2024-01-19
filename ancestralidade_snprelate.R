##############
# Libraries  #
##############

library("gdsfmt")
library("SNPRelate")
library("dplyr")
library("scatterplot3d")
library("ca")
library("ggplot2")
library("tidyr")

###################################################

##################
# Open Databases #
##################

setwd('/Users/dasacamila/Desktop/USP-Paper')


bed.fn <- "isamerge1000.bed"
bim.fn <- "isamerge1000.bim"
fam.fn <- "isamerge1000.fam"

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, out.gdsfn="isamerge1000.gds", verbose=TRUE)
genofile <- snpgdsOpen("isamerge1000.gds")

# ISA after LD filter

bed2.fn <- "isa_hg37.prunedld.bed"
bim2.fn <- "isa_hg37.prunedld.bim"
fam2.fn <- "isa_hg37.prunedld.fam"

snpgdsBED2GDS(bed2.fn, fam2.fn, bim2.fn, out.gdsfn="isa.filter.gds", verbose=TRUE)
genofile2 <- snpgdsOpen("isa.filter.gds")

# ISA before LD filter

bed3.fn <- "isa_hg37_filter2.bed"
bim3.fn <- "isa_hg37_filter2.bim"
fam3.fn <- "isa_hg37_filter2.fam"

snpgdsBED2GDS(bed3.fn, fam3.fn, bim3.fn, out.gdsfn="isafilter3.gds", verbose=TRUE)
genofile3 <- snpgdsOpen("isafilter3.gds")

# Related Samples 

bed4.fn <- "isacomm_excluidas_2.bed"
bim4.fn <- "isacomm_excluidas_2.bim"
fam4.fn <- "isacomm_excluidas_2.fam"

snpgdsBED2GDS(bed4.fn, fam4.fn, bim4.fn, out.gdsfn="isafilter4.gds", verbose=TRUE)
genofile4 <- snpgdsOpen("isafilter4.gds")

# Raw data - All SNPs

bed5.fn <- "isa_hg37.bed"
bim5.fn <- "isa_hg37.bim"
fam5.fn <- "isa_hg37.fam"

snpgdsBED2GDS(bed5.fn, fam5.fn, bim5.fn, out.gdsfn="isafilter5.gds", verbose=TRUE)
genofile5 <- snpgdsOpen("isafilter5.gds")

###################################################

##################
# SNP Frequency  #
##################

# After LD filter

snplist <- snpgdsSNPList(genofile2)
snplistchr <- data.frame(snplist$chromosome, snplist$snp.id)
freq <- as.data.frame(table(snplistchr$snplist.chromosome))
color1 <- colorRampPalette(c("darkblue","lightblue"))

sum(freq$Freq)
barplot(freq$Freq,
        col = color1(22), 
        main = "Quantidade de SNPs por cromossomo",
        ylab = "Numero de SNPs",
        xlab = "Cromossomos",
        ylim=c(0,60000),
        names.arg = c(freq$Var1))

p <- ggplot(freq, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Quantidade de SNPs por cromossomo",
       x = "Cromossomos",
       y = "Numero de SNPs") +
  ylim(0, 60000) +
  theme_minimal()

print(p)

# Before LD filter

snplist2 <- snpgdsSNPList(genofile3)
snplistchr2 <- data.frame(snplist2$chromosome, snplist2$snp.id)
freq2 <- as.data.frame(table(snplistchr2$snplist2.chromosome))
color2 <- colorRampPalette(c("darkgreen","lightgreen"))

sum(freq2$Freq)
barplot(freq2$Freq,
        col = color2(22), 
        main = "Quantidade de SNPs por cromossomo sem filtro de LD",
        ylab = "Numero de SNPs",
        xlab = "Cromossomos",
        ylim=c(0,80000),
        names.arg = c(freq2$Var1))

# All SNPs

snplist5 <- snpgdsSNPList(genofile5)
snplistchr5 <- data.frame(snplist5$chromosome, snplist5$snp.id)
freq5 <- as.data.frame(table(snplistchr5$snplist5.chromosome))
color5 <- colorRampPalette(c("darkgreen","lightgreen"))

sum(freq5$Freq)
barplot(freq5$Freq,
        col = "darkred", 
        main = "Quantidade de SNPs por cromossomo - Todos os Marcadores",
        ylab = "Numero de SNPs",
        xlab = "Cromossomos",
        ylim=c(0,80000),
        names.arg = c(freq5$Var1))

p <- ggplot(freq5, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "darkred") +
  labs(title = "Quantidade de SNPs por cromossomo - Todos Marcadores",
       x = "Cromossomos",
       y = "Numero de SNPs") +
  ylim(0, 60000) +
  theme_minimal()

print(p)


###################################################

###################
# Calculating PCA #
###################

# 681 samples and 1000g

CP <- snpgdsPCA(genofile)
head(CP)
summary(CP)

# ISA after LD filter

CP2 <- snpgdsPCA(genofile2)

###################################################

#############################################
# Saving eigenvectors and eigenval on file  #
#############################################

write.table(CP$eigenvec, file="isamerge1000_eigenvec.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

write.table(CP$eigenvect, file="eigenvec.txt", col.names=TRUE,
            row.names=TRUE, sep="\t" , quote=FALSE)

write.table(CP$sample.id, file="samples.txt", col.names=TRUE,
            row.names=TRUE, sep="\t" , quote=FALSE)

write.table(CP$eigenval, file="isamerge1000_eigenval.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

write.table(CP$eigenval, file="isamerge1000_eigenval.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

write.table(CP$snp.id, file="snps.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

write.table(CP2$eigenvec, file="isa_eigenvec.txt", col.names=TRUE,
            row.names=TRUE, sep=",")

write.table(CP2$eigenval, file="isa_eigenval.txt", col.names=TRUE,
            row.names=TRUE, sep=",")

summary(CP$eigenvect)
summary(CP2$eigenvect)
summary(CP$eigenval)
summary(CP2$eigenval)

write.table(CP$eigenval[1:6], file="summary_isamerge_eigenval1-6.txt", col.names=TRUE,
            row.names=TRUE, sep=",")
write.table(CP2$eigenval[1:6], file="summary_isa_eigenval1-6.txt", col.names=TRUE,
            row.names=TRUE, sep=",")
write.table(CP$eigenvect[1:6], file="summary_isamerge_eigenvect1-6.txt", col.names=TRUE,
            row.names=TRUE, sep=",")
write.table(CP2$eigenvect[1:6], file="summary_isa_eigenvect1-6.txt", col.names=TRUE,
            row.names=TRUE, sep=",")

###################################################

#########################################
#  Explanation percentage of components #
#########################################

# 681 samples and 1000g

pc.percent <- CP$eigenval[1:100] / sum(CP$eigenval[1:25])
write.table(pc.percent, file="percent_eigenval_isamerge1000.txt", sep=",", quote = FALSE)
cpplot <- c(pc.percent[1:6])
cp1cp2 <- cpplot[1] + cpplot[2]

barplot(cpplot,
        beside = TRUE, 
        width = 0.7, 
        ylim = c(0,0.6),
        col=1:32,
        legend = c("CP1","CP2","CP3","CP4","CP5","CP6"),
        main = "Explanation percentage of components - Isa and 1000g")

# ISA after LD filter

pc.percent2 <- CP2$eigenval[1:100] / sum(CP2$eigenval[1:25])
write.table(pc.percent2, file="percent_eigenval_isa.txt", sep=",")

cpplot2 <- c(pc.percent2[1:6])
cp1cp22 <- cpplot2[1] + cpplot2[2]

barplot(cpplot2,
        beside = TRUE, 
        width = 0.7, 
        ylim = c(0,0.6),
        col=1:8,
        legend = c("CP1","CP2","CP3","CP4","CP5","CP6"),
        main = "Explanation percentage of components - Isa")

#############################################################

#########################
# Building the PCA Plot #
#########################

par(mai=c(1.3, 1.1, .2, 0.2))
plot(CP$eigenvect[ , 2] ~ CP$eigenvect[ , 1],col=rgb(0,0,150,50,maxColorValue=255),las=1,
     pch=19,xlab="PC1",ylab="PC2", main = "Distribution - ISA and 1000g")

par(mai=c(1.3, 1.1, .2, 0.2))
plot(CP2$eigenvect[ , 2] ~ CP2$eigenvect[ , 1],col=rgb(0,0,155,155,maxColorValue=255),las=1,
     pch=19,xlab="PC1",ylab="PC2", main ="Distribution - ISA" )


##############################################################

#######################################################
# Plot of the relationships of the first 6 components #
#######################################################

names.PCA <- paste("PC", 1:6, "\n", format(pc.percent[1:6]), "%", sep="")
pairs(CP$eigenvec[ , 1:6], labels=names.PCA, 
      main = "Relationships of the first 6 components")

names.PCA2 <- paste("PC", 1:6, "\n", format(pc.percent2[1:6]), "%", sep="")
pairs(CP2$eigenvec[ , 1:6], labels=names.PCA2, 
      main = "Relationships of the first 6 components - ISA")

##############################################################

#############################################
# Correlation between eigenvectors and SNPs #
#     Obtaining the chromosome index        #
#############################################

# 681 and 1000g

chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
CORR <- snpgdsPCACorr(CP, genofile, eig.which=1:6)

# CP1
plot(abs(CORR$snpcorr[1,]), ylim=c(0,1.0), xlab="SNP Index",
     ylab=paste("PC", 1), main="Correlation CP1 - chromosomes 1 a 22", 
     col=chr, pch="+")

# CP2
plot(abs(CORR$snpcorr[2,]), ylim=c(0,1.0), xlab="SNP Index",
     ylab=paste("PC", 2), main="Correlation CP2 - chromosomes 1 a 22", 
     col=chr, pch="+")

# ISA after LD filter

chr2 <- read.gdsn(index.gdsn(genofile2, "snp.chromosome"))
CORR2 <- snpgdsPCACorr(CP2, genofile2, eig.which=1:6)

# CP1
plot(abs(CORR2$snpcorr[1,]), ylim=c(0,1), xlab="SNP Index",
     ylab=paste("PC", 1), main="Correlation CP1 - chromosomes 1 a 22 - ISA", 
     col=chr2, pch="+")

# CP2
plot(abs(CORR2$snpcorr[2,]), ylim=c(0,1), xlab="SNP Index",
     ylab=paste("PC", 2), main="Correlation CP2 - chromosomes 1 a 22 - ISA", 
     col=chr2, pch="+")

##############################################################

#############################
# Plot Ancestry Proportions #
#############################

### 681 and 1000g

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
group <- read.delim("isamerge1000g_2.txt")

tab <- data.frame(sample.id = samp.id, group = factor(group$SuperPop),
                  CP1 = CP$eigenvect[,1],    
                  CP2 = CP$eigenvect[,2],    
                  stringsAsFactors = FALSE)

plot(tab$CP1, tab$CP2, col=as.integer(tab$group),
     xlab="CP 1", ylab="CP 2", pch = 16,
     main = "Ancestry Proportions - SuperPopulations - 681")
legend("topright", legend=levels(tab$group), col=1:6, cex = 0.8, pch = 16)

##############################################################

#############################################
# Plot Ancestry Proportions - Declared Race #
#############################################

### 681 and 1000g

write.table(ISA_Camila_Agosto2023_n_681, file="681.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote=FALSE)

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
group2 <- read.delim("isamerge1000g-race_2.txt")

tab2 <- data.frame(sample.id = samp.id, group = factor(group2$SuperPop),
                   CP1 = CP$eigenvect[,1],    
                   CP2 = CP$eigenvect[,2],    
                   stringsAsFactors = FALSE)
summary(tab2)

# All colors

plot(tab2$CP1, tab2$CP2, col=c("black","deeppink3","blue","royalblue","seagreen",
                               "limegreen","green2","green3","green4", "darkgreen",
                               "darkviolet")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681")
legend("topright", legend=levels(tab2$group), col=c("black","deeppink3",
                                                    "blue","royalblue","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'darkviolet'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))


# ISA colored

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "limegreen","green2","green3","green4", "darkgreen",
                               "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681 - ISA")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'gray'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Without Brancos

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "limegreen","green2","green3","gray", "darkgreen",
                               "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681 - sem Brancos")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "limegreen","green2","green3",
                                                    "gray","darkgreen",'gray'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Brancos

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "gray","gray","gray","green4", "gray",
                                       "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681 - ISA - Brancos")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","gray","gray",
                                                           "green4","gray",'gray'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Amarelos

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "gray","gray","gray","gray", "darkgreen",
                                       "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681 - Amarelos")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","darkgreen",'gray'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Indigenas

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "limegreen","gray","gray","gray", "gray",
                                       "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681 - Indigenas")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "limegreen","gray","gray",
                                                           "gray","gray",'gray'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Outras

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "gray","gray","green3","gray", "gray",
                                       "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681 - Outras")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","gray","green3",
                                                           "gray","gray",'gray'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Parda

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "gray","green2","gray","gray", "gray",
                                       "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681 - Parda")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","green2","gray",
                                                           "gray","gray",'gray'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Preta

plot(tab2$CP1, tab2$CP2, col=c("gray","gray","gray","gray","seagreen",
                                       "gray","gray","gray","gray", "gray",
                                       "gray")[tab2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab2$group],
     main = "Ancestry Proportions - Declared Race - 681 - Preta")
legend("topright", legend=levels(tab2$group), col=c("gray","gray",
                                                           "gray","gray","seagreen",
                                                           "gray","gray","gray",
                                                           "gray","gray",'gray'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,16,19,15,3,16))

##############################################################

##########################################
# Plot Ancestry Proportions - birthplace #
##########################################

### 681 and 1000g

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
group3 <- read.delim("isamerge1000g-nascimento_2.txt")
tail(group3)

tab3 <- data.frame(sample.id = samp.id, group = factor(group3$SuperPop),
                   CP1 = CP$eigenvect[,1],    
                   CP2 = CP$eigenvect[,2],    
                   stringsAsFactors = FALSE)

# All colors

plot(tab3$CP1, tab3$CP2, col=c("black","deeppink3",
                               "blue","royalblue","seagreen",
                               "limegreen","green2","green3",
                               "green4","darkgreen",'darkgreen')[tab3$group],
     xlab="CP 1", ylab="CP 2", pch = c(16,16,16,16,17,8,18,19,15,3,4)[tab3$group],
     main = "Ancestry Proportions - birthplace")
legend("bottomright", legend=levels(tab3$group), col=c("black","deeppink3",
                                                    "blue","royalblue","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'darkgreen'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,18,19,15,3,4))


# ISA colored

plot(tab3$CP1, tab3$CP2, col=c("gray","gray",
                               "gray","gray","seagreen",
                               "limegreen","green2","green3",
                               "green4","darkgreen",'darkgreen')[tab3$group],
     xlab="CP 1", ylab="CP 2", pch = c(16,16,16,16,17,8,18,19,15,3,4)[tab3$group],
     main = "Ancestry Proportions - birthplace")
legend("bottomright", legend=levels(tab3$group), col=c("gray","gray",
                                                       "gray","gray","seagreen",
                                                       "limegreen","green2","green3",
                                                       "green4","darkgreen",'darkgreen'), 
       cex = 0.65, pch = c(16,16,16,16,17,8,18,19,15,3,4))

##############################################################

##########################################
# Plot Ancestry Proportions - Individual #
##########################################

### 681 and 1000g

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
groups <- list(AMR = samp.id[group$SuperPop == "AMR"],
               AFR = samp.id[group$SuperPop == "AFR"],
               EUR = samp.id[group$SuperPop == "EUR"])

prop <- snpgdsAdmixProp(CP, groups=groups,bound=TRUE)
prop_sort <- sort(prop)
propt <- t(prop)

### Save Admix Props on file

write.table(propt, file="admix-props_681samples.txt", col.names=TRUE,
            row.names=TRUE, sep="\t")

## All samples

BRpropt <- t(prop[tab$group=='BR - ISA',])
barplot(BRpropt, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA - 681",space = FALSE,
        ylim=c(0,1))
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Declared Race - Brancos

Branca <- (prop[tab2$group=='ISA - White',])
Brancat <- t(na.omit(Branca))

barplot(Brancat, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race White-ISA", 
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

# Declared Race - Pardos

Parda <- (prop[tab2$group=='ISA - Mixed',])
Pardat <- t(na.omit(Parda))

barplot(Pardat, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Mixed-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

# Declared Race - Negros

Preta <- (prop[tab2$group=='ISA - Black',])
Pretat <- t(na.omit(Preta))

barplot(Pretat, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Black-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

###############################################################################

### Suggestion to organize the order of the figures

## All samples

# Test 1

BRpropt <- t(prop[tab$group=='BR - ISA',])
TBRpropt <- data.frame(t(BRpropt))
T1BRpropt_ord  <- TBRpropt[order(TBRpropt$AFR),]
T3BRpropt_ord  <- T1BRpropt_ord[order(T1BRpropt_ord$AMR),]
BRpropt_ord <- t(T3BRpropt_ord)

barplot(BRpropt_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Test 2

BRpropt <- t(prop[tab$group=='BR - ISA',])
TBRpropt <- data.frame(t(BRpropt))
T1BRpropt_ord  <- TBRpropt[order(TBRpropt$AMR),]
T2BRpropt_ord  <- T1BRpropt_ord[order(T1BRpropt_ord$AFR),]
BRpropt_ord <- t(T2BRpropt_ord)

barplot(BRpropt_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Test 3

#BRpropt <- t(prop[tab$group=='BR - ISA',])
#TBRpropt <- data.frame(t(BRpropt))
#T2BRpropt_ord  <- TBRpropt_ord[order(TBRpropt_ord$EUR),]
#T3BRpropt_ord  <- T2BRpropt_ord[order(T2BRpropt_ord$AMR),]
#BRpropt_ord <- t(T3BRpropt_ord)

#barplot(BRpropt_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        #border=NA,axisnames=FALSE,main="Ancestry of BR-ISA",
        #ylim=c(0,1),space=0)
#legend("bottomright", c("AMR","AFR","EUR"),
       #lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

## Declared Race - Brancos - Test 1

Branca <- (prop[tab2$group=='ISA - White',])
Brancat <- t(na.omit(Branca))
TBrancat <- data.frame(t(Brancat))
T1Brancat_ord  <- TBrancat[order(TBrancat$AFR),]
T2Brancat_ord  <- T1Brancat_ord[order(T1Brancat_ord$AMR),]
Brancapropt_ord <- t(T2Brancat_ord)

barplot(Brancapropt_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race White-ISA", 
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

## Declared Race - Pardos - Test 1

Parda <- (prop[tab2$group=='ISA - Mixed',])
Pardat <- t(na.omit(Parda))
TPardat <- data.frame(t(Pardat))
T1Pardat_ord  <- TPardat[order(TPardat$AFR),]
T2Pardat_ord  <- T1Pardat_ord[order(T1Pardat_ord$AMR),]
Pardat_ord <- t(T2Pardat_ord)

barplot(Pardat_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Mixed-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

## Declared Race - Negros - Test 1

Preta <- (prop[tab2$group=='ISA - Black',])
Pretat <- t(na.omit(Preta))
TPretat <- data.frame(t(Pretat))
T1Pretat_ord  <- TPretat[order(TPretat$AFR),]
T2Pretat_ord  <- T1Pretat_ord[order(T1Pretat_ord$AMR),]
Pretat_ord <- t(T2Pretat_ord)

barplot(Pretat_ord, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of Declared Race Black-ISA",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

###############################################################################

prop_isa <- as.data.frame(prop, stringsAsFactors = FALSE)
rownames(prop_isa) <- rownames(prop)

prop_isa_filt <- prop_isa %>%
  filter(grepl('^a550778', rownames(prop_isa)))
summary(prop_isa_filt)

isa_race <- data.frame(group2$SuperPop[1586:2266])
prop_isa_final <- data.frame(prop_isa_filt,isa_race)

mixed <- subset(prop_isa_final, group2.SuperPop.1586.2266. == "ISA - Mixed")
white <- subset(prop_isa_final, group2.SuperPop.1586.2266. == "ISA - White")
black <- subset(prop_isa_final, group2.SuperPop.1586.2266. == "ISA - Black")
yellow <- subset(prop_isa_final, group2.SuperPop.1586.2266. == "ISA - Yellow")
indigenous <- subset(prop_isa_final, group2.SuperPop.1586.2266. == "ISA - Indigenous")
na <- subset(prop_isa_final, group2.SuperPop.1586.2266. == "ISA - NA")

media_mixed <- colMeans(mixed[, 1:3])
media_white <- colMeans(white[, 1:3])
media_black <- colMeans(black[, 1:3])
media_yellow <- colMeans(yellow[, 1:3])
media_indigenous <- colMeans(indigenous[, 1:3])
media_na <- colMeans(na[, 1:3])

means_prop <- data.frame(media_black,media_indigenous,media_mixed,media_na,media_white,media_yellow)
means_propt <- data.frame(t(means_prop))
rownames(means_propt) = c("Black","Indigenous","Mixed","NA","White","Yellow")
means_propt$types <- rownames(means_propt)
means_propt

dados <- data.frame(
  Categoria = c("Black", "Indigenous", "Mixed", "NA", "White", "Yellow"),
  AMR = c(0.09074315, 0.25409154, 0.12343520, 0.09608049, 0.09119493, 1.00000000),
  AFR = c(0.5343041, 0.2663917, 0.3035574, 0.1757728, 0.1369259, 0.0000000),
  EUR = c(0.3749527, 0.4795167, 0.5730074, 0.7281467, 0.7718792, 0.0000000)
)

dados_long <- dados %>%
  pivot_longer(cols = c(AMR, AFR, EUR), names_to = "Região", values_to = "Valor")

prop_plot <- ggplot(dados_long, aes(x = Categoria, y = Valor, fill = Região)) +
  geom_bar(stat = "identity") +
  labs(title = "Average proportions by declared race",
       x = "Categoria",
       y = "Valor") +
  #geom_text(aes(label = round(Valor, 2)), position = position_stack(vjust = 0.5)) +
  geom_text(aes(label = round(Valor, 2)),
            position = position_stack(vjust = 0.5),
            fontface = "bold") +
  scale_fill_manual(values = c("AMR" = "red", "AFR" = "green", "EUR" = "blue")) +
  theme_minimal() +
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = NULL))+
  xlab(NULL) +
  ylab(NULL)

prop_plot
ggsave(filename = "avg_prop.png", plot = prop_plot)

###############################################################################

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
groups <- list(AMR = samp.id[group$SuperPop == "AMR"],
               AFR = samp.id[group$SuperPop == "AFR"],
               EUR = samp.id[group$SuperPop == "EUR"],
               EAS = samp.id[group$SuperPop == "EAS"])

prop_4 <- snpgdsAdmixProp(CP, groups=groups,bound=TRUE)
prop_isa2 <- as.data.frame(prop_4, stringsAsFactors = FALSE)
rownames(prop_isa2) <- rownames(prop_4)

prop_isa_filt2 <- prop_isa2 %>%
  filter(grepl('^a550778', rownames(prop_isa2)))
summary(prop_isa_filt2)

isa_race <- data.frame(group2$SuperPop[1586:2266])
prop_isa_final2 <- data.frame(prop_isa_filt2,isa_race)

mixed <- subset(prop_isa_final2, group2.SuperPop.1586.2266. == "ISA - Mixed")
white <- subset(prop_isa_final2, group2.SuperPop.1586.2266. == "ISA - White")
black <- subset(prop_isa_final2, group2.SuperPop.1586.2266. == "ISA - Black")
yellow <- subset(prop_isa_final2, group2.SuperPop.1586.2266. == "ISA - Yellow")
indigenous <- subset(prop_isa_final2, group2.SuperPop.1586.2266. == "ISA - Indigenous")
na <- subset(prop_isa_final2, group2.SuperPop.1586.2266. == "ISA - NA")

media_mixed <- colMeans(mixed[, 1:4])
media_white <- colMeans(white[, 1:4])
media_black <- colMeans(black[, 1:4])
media_yellow <- colMeans(yellow[, 1:4])
media_indigenous <- colMeans(indigenous[, 1:4])
media_na <- colMeans(na[, 1:4])

means_prop <- data.frame(media_black,media_indigenous,media_mixed,media_na,media_white,media_yellow)
means_propt <- data.frame(t(means_prop))
rownames(means_propt) = c("Black","Indigenous","Mixed","NA","White","Yellow")
means_propt$types <- rownames(means_propt)
means_propt

dados <- data.frame(
  Categoria = c("Black", "Indigenous", "Mixed", "NA", "White", "Yellow"),
  AMR = c(0.09074315, 0.25409154, 0.12343520, 0.09608049, 0.09119493, 1.00000000),
  AFR = c(0.5343041, 0.2663917, 0.3035574, 0.1757728, 0.1369259, 0.0000000),
  EUR = c(0.3749527, 0.4795167, 0.5730074, 0.7281467, 0.7718792, 0.0000000),
  EAS = c(0.0009600817,0.0022235487,0.0030197046,0.0014159880,0.0073365995,0.0000000000)
)

dados_long <- dados %>%
  pivot_longer(cols = c(AMR, AFR, EUR,EAS), names_to = "Região", values_to = "Valor")

prop_plot <- ggplot(dados_long, aes(x = Categoria, y = Valor, fill = Região)) +
  geom_bar(stat = "identity") +
  labs(title = "Average proportions by declared race",
       x = "Categoria",
       y = "Valor") +
  #geom_text(aes(label = round(Valor, 2)), position = position_stack(vjust = 0.5)) +
  geom_text(aes(label = round(Valor, 2)),
            position = position_stack(vjust = 0.5),
            fontface = "bold") +
  scale_fill_manual(values = c("AMR" = "red", "AFR" = "green", "EUR" = "blue", "EAS" = "yellow")) +
  theme_minimal() +
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = NULL))+
  xlab(NULL) +
  ylab(NULL)

prop_plot
ggsave(filename = "avg_prop_4pop.png", plot = prop_plot)


###############################################################################



################################
# Get the most correlated SNPs #
################################

## 681 ISA samples and 1000genomes

CORR_data_cp1 <- data.frame(CORR$snp.id,CORR$snpcorr[1,])
CORR_data_cp1_ord  <- CORR_data_cp1[order(-CORR_data_cp1$CORR.snpcorr.1...),]
head(CORR_data_cp1_ord)
write.table(CORR_data_cp1_ord, file="CORR_data_cp1_ord_681samplesAND1000g.csv", sep=","
            ,quote = FALSE)

CORR_data_cp2 <- data.frame(CORR$snp.id,CORR$snpcorr[2,])
CORR_data_cp2_ord  <- CORR_data_cp2[order(-CORR_data_cp2$CORR.snpcorr.2...),]
head(CORR_data_cp2_ord)
write.table(CORR_data_cp2_ord, file="CORR_data_cp2_ord_681samplesAND1000g.csv", sep=","
            ,quote = FALSE)

## Only 681 ISA samples

CORR2_data_cp1 <- data.frame(CORR2$snp.id,CORR2$snpcorr[1,])
CORR2_data_cp1_ord  <- CORR2_data_cp1[order(-CORR2_data_cp1$CORR2.snpcorr.1...),]
head(CORR2_data_cp1_ord)
write.table(CORR2_data_cp1_ord, file="CORR2_data_cp1_ord_681samples.csv", sep=","
            ,quote = FALSE)

CORR2_data_cp2 <- data.frame(CORR2$snp.id,CORR2$snpcorr[2,])
CORR2_data_cp2_ord  <- CORR2_data_cp2[order(-CORR2_data_cp2$CORR2.snpcorr.2...),]
head(CORR2_data_cp2_ord)
write.table(CORR2_data_cp2_ord, file="CORR2_data_cp2_ord_681samples.csv", sep=","
            ,quote = FALSE)

###############################################################################

###################################################
#              Joint Related Samples              #
# Calculate sample eigenvectors from SNP loadings #
###################################################

### Loading SNPs from 681 samples + 1000g  

SnpLoad <- snpgdsPCASNPLoading(CP, genofile)
names(SnpLoad)
dim(SnpLoad$snploading)

### PCA Loading - Related Samples

sample.id <- read.gdsn(index.gdsn(genofile4, "sample.id"))
samp_load <- snpgdsPCASampLoading(SnpLoad, genofile4, sample.id=sample.id)

# Check CP

#diff <- CP$eigenvect[1:100,] - samp_load$eigenvect
#summary(c(diff))

### Joint PCA Class

isamerge1000gfull <- list(
  sample.id = c(CP$sample.id, samp_load$sample.id),
  snp.id = CP$snp.id,
  eigenval = c(CP$eigenval, samp_load$eigenval),
  eigenvect = rbind(CP$eigenvect, samp_load$eigenvect),
  varprop = c(CP$varprop, samp_load$varprop),
  TraceXTX = CP$TraceXTX
)
class(isamerge1000gfull) <- "snpgdsPCAClass"

### New Database - All samples

isamerge1000gfull

###############################################################################

###########################################
# Plot Ancestry Proportions - All Samples #
###########################################

### Related Samples Colored

samp.id_all <- isamerge1000gfull$sample.id
group_all <- read.delim("isamerge1000g_all_relacionados.txt")

tab_all <- data.frame(sample.id = samp.id_all, group = factor(group_all$SuperPop),
                      CP1 = isamerge1000gfull$eigenvect[,1],    
                      CP2 = isamerge1000gfull$eigenvect[,2],    
                      stringsAsFactors = FALSE)

plot(tab_all$CP1, tab_all$CP2, col=as.integer(tab_all$group),
     xlab="CP 1", ylab="CP 2", pch = 16,
     main = "Ancestry Proportions - SuperPopulations - All samples + 1000g")
legend("bottomright", legend=levels(tab_all$group), col=1:6, cex = 0.6, pch = 16)

### All Samples

group_all1 <- read.delim("isamerge1000g_all.txt")
tab_all1 <- data.frame(sample.id = samp.id_all, 
                       group = factor(group_all1$SuperPop),
                      CP1 = isamerge1000gfull$eigenvect[,1],    
                      CP2 = isamerge1000gfull$eigenvect[,2],    
                      stringsAsFactors = FALSE)

plot(tab_all1$CP1, tab_all1$CP2, col=as.integer(tab_all1$group),
     xlab="CP 1", ylab="CP 2", pch = 16,
     main = "Ancestry Proportions - SuperPopulations -  All Samples")
legend("bottomright", legend=levels(tab_all1$group), col=1:6, cex = 0.8, pch = 16)

###############################################################################

###########################################################
# Plot Ancestry Proportions - Declared Race - All Samples #
###########################################################

group_all2 <- read.delim("isamerge1000g_all_race.txt")
tail(group_all2)
summary(group_all2$SuperPop)

tab_all2 <- data.frame(sample.id = samp.id_all, group = factor(group_all2$SuperPop),
                      CP1 = isamerge1000gfull$eigenvect[,1],    
                      CP2 = isamerge1000gfull$eigenvect[,2],    
                      stringsAsFactors = FALSE)
summary(tab_all2)
tail(tab_all2)


plot(tab_all2$CP1, tab_all2$CP2, col=c("black","deeppink3","blue","royalblue","seagreen",
                               "limegreen","green2","green3","green4", "darkgreen",
                               "darkviolet")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g")
legend("bottomright", legend=levels(tab_all2$group), col=c("black","deeppink3",
                                                    "blue","royalblue","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'darkviolet'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))
    
# Only ISA colored

plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "limegreen","green2","green3","green4", "darkgreen",
                               "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4","darkgreen",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))


plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "gray","green2","green3","green4", "darkgreen",
                               "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g - Sem Brancos")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "gray","green2","green3",
                                                    "green4","darkgreen",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Brancos

plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "limegreen","gray","gray","gray", "gray",
                                       "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g - Brancos")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "limegreen","gray","gray",
                                                           "gray","gray",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Amarelos

plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","seagreen",
                                       "gray","gray","gray","gray", "gray",
                                       "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g - Amarelos")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                           "gray","gray","seagreen",
                                                           "gray","gray","gray",
                                                           "gray","gray",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Indigenas

plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "gray","green2","gray","gray", "gray",
                                       "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g - Indigenas")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","green2","gray",
                                                           "gray","gray",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Outras

plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "gray","gray","green3","gray", "gray",
                                       "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g - Outras")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","gray","green3",
                                                           "gray","gray",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Parda

plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "gray","gray","gray","green4", "gray",
                                       "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g - Parda")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","gray","gray",
                                                           "green4","gray",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))

# Only Preta

plot(tab_all2$CP1, tab_all2$CP2, col=c("gray","gray","gray","gray","gray",
                                       "gray","gray","gray","gray", "darkgreen",
                                       "gray")[tab_all2$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,17,8,16,19,15,3,20)[tab_all2$group],
     main = "Ancestry Proportions - Declared Race - All samples + 1000g - Preta")
legend("bottomright", legend=levels(tab_all2$group), col=c("gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","gray","gray",
                                                           "gray","darkgreen",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,17,8,16,19,15,3,16))


###############################################################################

########################################################
# Plot Ancestry Proportions - birthplace - All Samples #
########################################################

group_all3 <- read.delim("isamerge1000g_all_nascimento.txt")

tab_all3 <- data.frame(sample.id = samp.id_all, group = factor(group_all3$SuperPop),
                   CP1 = isamerge1000gfull$eigenvect[,1],    
                   CP2 = isamerge1000gfull$eigenvect[,2],    
                   stringsAsFactors = FALSE)

plot(tab_all3$CP1, tab_all3$CP2, col=c("black","deeppink3","blue","royalblue","seagreen",
                               "limegreen","green2","green3","green4", 
                               "darkviolet")[tab_all3$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,8,17,16,15,4,20)[tab_all3$group],
     main = "Ancestry Proportions - Local de Nascimento - All samples + 1000g")
legend("bottomright", legend=levels(tab_all3$group), col=c("black","deeppink3",
                                                    "blue","royalblue","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4",'darkviolet'), 
       cex = 0.60, pch = c(16,16,16,16,8,17,16,15,4,16))

# Only ISA colored

plot(tab_all3$CP1, tab_all3$CP2, col=c("gray","gray","gray","gray","seagreen",
                               "limegreen","green2","green3","green4", 
                               "gray")[tab_all3$group],
     xlab="CP 1", ylab="CP 2", pch = c(20,20,20,20,8,17,16,15,4,20)[tab_all3$group],
     main = "Ancestry Proportions - Local de Nascimento - All samples + 1000g")
legend("bottomright", legend=levels(tab_all3$group), col=c("gray","gray",
                                                    "gray","gray","seagreen",
                                                    "limegreen","green2","green3",
                                                    "green4",'gray'), 
       cex = 0.6, pch = c(16,16,16,16,8,17,16,15,4,16))


###############################################################################

########################################################
# Plot Ancestry Proportions - Individual - All Samples #
########################################################

groups4 <- list(AMR = samp.id_all[group$SuperPop == "AMR"],
               AFR = samp.id_all[group$SuperPop == "AFR"],
               EUR = samp.id_all[group$SuperPop == "EUR"])

prop4 <- snpgdsAdmixProp(isamerge1000gfull, groups=groups4,bound=TRUE)
prop_sort4 <- sort(prop4)
propt4 <- t(prop4)
head(prop4)

### Save Admix Prop - All samples

write.table(prop4, file="admix-props_all.txt", col.names=TRUE,
            row.names=TRUE, sep="\t", quote = FALSE)

### All Samples

BRpropt4 <- t(prop4[tab_all$group=='BR - ISA',])
barplot(BRpropt4, col=c("red","green","blue"),xlab="Individual", 
        ylab="Ancestry", 
        border=NA,axisnames=FALSE,
        main="Ancestry of BR-ISA - All samples + 1000g",space = FALSE,
        ylim=c(0,1))
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)


# Test 1

BRpropt4 <- t(prop4[tab_all$group=='BR - ISA',])
TBRpropt4 <- data.frame(t(BRpropt4))
T1BRpropt_ord4  <- TBRpropt4[order(TBRpropt4$AFR),]
T3BRpropt_ord4  <- T1BRpropt_ord4[order(T1BRpropt_ord4$AMR),]
BRpropt_ord4 <- t(T3BRpropt_ord4)

barplot(BRpropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA - All samples - 1",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Test 2

BRpropt4 <- t(prop4[tab_all$group=='BR - ISA',])
TBRpropt4 <- data.frame(t(BRpropt4))
T1BRpropt_ord4  <- TBRpropt4[order(TBRpropt4$AMR),]
T2BRpropt_ord4  <- T1BRpropt_ord4[order(T1BRpropt_ord4$AFR),]
BRpropt_ord4 <- t(T2BRpropt_ord4)

barplot(BRpropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA - All samples - 2",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Test 3

BRpropt4 <- t(prop4[tab_all$group=='BR - ISA',])
TBRpropt4 <- data.frame(t(BRpropt4))
T1BRpropt_ord4  <- TBRpropt4[order(TBRpropt4$EUR),]
T2BRpropt_ord4  <- T1BRpropt_ord4[order(T1BRpropt_ord4$AMR),]
BRpropt_ord4 <- t(T2BRpropt_ord4)

barplot(BRpropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA - All samples - 3",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)


# Test 4 

BRpropt4 <- t(prop4[tab_all$group=='BR - ISA',])
TBRpropt4 <- data.frame(t(BRpropt4))
T1BRpropt_ord4  <- TBRpropt4[order(TBRpropt4$AMR),]
T2BRpropt_ord4  <- T1BRpropt_ord4[order(T1BRpropt_ord4$EUR),]
BRpropt_ord4 <- t(T2BRpropt_ord4)

barplot(BRpropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA - All samples - 4",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

# Teste 5 

BRpropt4 <- t(prop4[tab_all$group=='BR - ISA',])
TBRpropt4 <- data.frame(t(BRpropt4))
T1BRpropt_ord4  <- TBRpropt4[order(TBRpropt4$AFR),]
T2BRpropt_ord4  <- T1BRpropt_ord4[order(T1BRpropt_ord4$EUR),]
BRpropt_ord4 <- t(T2BRpropt_ord4)

barplot(BRpropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestry of BR-ISA - All samples - 5",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)


### Declared Race - Brancos - All Samples

Branca4 <- (prop4[tab_all2$group=='Isa - Branca',])
Brancat4 <- t(na.omit(Branca4))

barplot(Brancat4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,
        main="Ancestry of Declared Race Brancos-ISA - All samples + 1000g", 
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)


Brancopropt4 <- t(prop4[tab_all2$group=='Isa - Branca',])
TBrancopropt4 <- data.frame(t(Brancopropt4))
T1Brancopropt_ord4  <- TBrancopropt4[order(TBrancopropt4$AFR),]
T2Brancopropt_ord4  <- T1Brancopropt_ord4[order(T1Brancopropt_ord4$EUR),]
Brancopropt_ord4 <- t(T2Brancopropt_ord4)

barplot(Brancopropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestralidade Brancos - All samples - 1",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

Brancopropt4 <- t(prop4[tab_all2$group=='Isa - Branca',])
TBrancopropt4 <- data.frame(t(Brancopropt4))
T1Brancopropt_ord4  <- TBrancopropt4[order(TBrancopropt4$AFR),]
T2Brancopropt_ord4  <- T1Brancopropt_ord4[order(T1Brancopropt_ord4$AMR),]
Brancopropt_ord4 <- t(T2Brancopropt_ord4)

barplot(Brancopropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestralidade Brancos - All samples - 2",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)


### Declared Race - Pardos - All Samples

Parda4 <- (prop4[tab_all2$group=='Isa - Parda',])
Pardat4 <- t(na.omit(Parda4))

barplot(Pardat4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,
        main="Ancestry of Declared Race Pardos-ISA - All samples + 1000g",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

Pardapropt4 <- t(prop4[tab_all2$group=='Isa - Parda',])
TPardapropt4 <- data.frame(t(Pardapropt4))
T1Pardapropt_ord4  <- TPardapropt4[order(TPardapropt4$AFR),]
T2Pardapropt_ord4  <- T1Pardapropt_ord4[order(T1Pardapropt_ord4$EUR),]
Pardapropt_ord4 <- t(T2Pardapropt_ord4)

barplot(Pardapropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestralidade Parda - All samples - 1",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

Pardapropt4 <- t(prop4[tab_all2$group=='Isa - Parda',])
TPardapropt4 <- data.frame(t(Pardapropt4))
T1Pardapropt_ord4  <- TPardapropt4[order(TPardapropt4$AFR),]
T2Pardapropt_ord4  <- T1Pardapropt_ord4[order(T1Pardapropt_ord4$AMR),]
Pardapropt_ord4 <- t(T2Pardapropt_ord4)

barplot(Pardapropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestralidade Parda - All samples - 2",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

### Declared Race - Negros - All Samples

Preta4 <- (prop4[tab_all2$group=='Isa - Preta',])
Pretat4 <- t(na.omit(Preta4))

barplot(Pretat4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,
        main="Ancestry of Declared Race Pretos-ISA  - All samples + 1000g",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),cex=0.5)

Pretapropt4 <- t(prop4[tab_all2$group=='Isa - Preta',])
TPretapropt4 <- data.frame(t(Pretapropt4))
T1Pretapropt_ord4  <- TPretapropt4[order(TPretapropt4$AFR),]
T2Pretapropt_ord4  <- T1Pretapropt_ord4[order(T1Pretapropt_ord4$EUR),]
Pretapropt_ord4 <- t(T2Pretapropt_ord4)

barplot(Pretapropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestralidade Preta - All samples - 1",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)


Pretapropt4 <- t(prop4[tab_all2$group=='Isa - Preta',])
TPretapropt4 <- data.frame(t(Pretapropt4))
T1Pretapropt_ord4  <- TPretapropt4[order(TPretapropt4$AFR),]
T2Pretapropt_ord4  <- T1Pretapropt_ord4[order(T1Pretapropt_ord4$AMR),]
Pretapropt_ord4 <- t(T2Pretapropt_ord4)

barplot(Pretapropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestralidade Preta - All samples - 2",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

Pretapropt4 <- t(prop4[tab_all2$group=='Isa - Preta',])
TPretapropt4 <- data.frame(t(Pretapropt4))
T1Pretapropt_ord4  <- TPretapropt4[order(TPretapropt4$AMR),]
T2Pretapropt_ord4  <- T1Pretapropt_ord4[order(T1Pretapropt_ord4$AFR),]
Pretapropt_ord4 <- t(T2Pretapropt_ord4)

barplot(Pretapropt_ord4, col=c("red","green","blue"),xlab="Individual ", ylab="Ancestry", 
        border=NA,axisnames=FALSE,main="Ancestralidade Preta - All samples - 3",
        ylim=c(0,1),space=0)
legend("bottomright", c("AMR","AFR","EUR"),
       lwd=4, col=c("red","green","blue"),bg="white",cex=0.5)

###############################################################################

###############################
# Plot Simplex Representation #
###############################

### Using Prop4 Data frame - All Proportions 

prop4_frame <- data.frame(prop4)
prop4_frame2 <- tibble::rownames_to_column (prop4_frame, "Samples")
simplex_data <- dplyr::filter(prop4_frame2, 
                              grepl('a550778',prop4_frame2$Samples))

# Statistics 

dim(simplex_data)
colMeans(simplex_data[,2:4])
cov(simplex_data[,2:4])
sqrt(var(simplex_data[,2]))
sqrt(var(simplex_data[,3]))
sqrt(var(simplex_data[,4]))

# 5 trinominal on Simplex Representation

scatterplot3d(simplex_data[,2:4], 
              main ="Representacao das 5 trinomiais no simplex",
              angle = 55,color = "darkblue",pch = 16)

# Color by Group

group_all2_2 <- read.delim("isamerge1000g_all_race_2.txt")
simplex_data_2 <- cbind.data.frame(prop4_frame2,group_all2_2$SuperPop)
head(simplex_data_2)
table(simplex_data_2$`group_all2_2$SuperPop`)

colors <- c("black","deeppink3","blue","royalblue","seagreen","limegreen",
      "green2","green3","green4", "darkgreen","olivedrab4")
colors <- colors[as.factor(simplex_data_2$`group_all2_2$SuperPop`)]
s3d <- scatterplot3d(simplex_data_2[,2:4], color=colors,
                     main = "Representacao Simplex 
                     Proporcoes de Ancestralidade por Raca Declarada - All Samples",
                     type="h", pch=10,  angle = 120)
legend("topright",
       legend=levels(as.factor(simplex_data_2$`group_all2_2$SuperPop`)),
       col = c("black","deeppink3","blue","royalblue","seagreen","limegreen",
               "green2","green3","green4", "darkgreen","olivedrab4"), 
       pch = 16,cex = 0.55,inset = 0.06)

write.csv(simplex_data_2, file = "admix_prop_simplex.csv")

# Correspondence Analysis

fit.ca <- ca(simplex_data_2[,2:4])