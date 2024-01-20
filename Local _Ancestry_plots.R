################################################################################
# Script for:
# Mean ancestry contribution by chromossomes.
# Ideograms and pie plots from infered data.
# Using ".Viterbi.txt" data from runRFmix
################################################################################

library(tidyverse)
library(ggrepel)
#library(devtools)
#devtools::install_github('TickingClock1992/RIdeogram')
library(RIdeogram)
library(ggplot2)


# 01 Mean ancestry contribution by chromosome -----------------------------

# proportions

freqtable <- data.frame()
for (i in seq(22)) {
  print(paste0('chr: ',i))
  
  tempvit <- read.csv(paste0('outs.rfmix/isa.rfmix_chr',i,'.rfmix.1.Viterbi.txt'),
                              sep=' ',header=F, colClasses='numeric')
  
  freqtable <- data.frame(rbind(
    freqtable,
    c(AFR = sum(rowSums(tempvit[,1:1414]==1)),
      AMR = sum(rowSums(tempvit[,1:1414]==2)),
      EAS = sum(rowSums(tempvit[,1:1414]==3)),
      EUR = sum(rowSums(tempvit[,1:1414]==4)) ) ) )
}

colnames(freqtable) <- c("AFR","AMR","EAS","EUR")

freqtable <- data.frame(cbind(
  freqtable, 
  sum.anc=rowSums(freqtable) ))

freqtable <- data.frame(cbind(
  freqtable,
  AFR.p = round((freqtable$AFR/freqtable$sum.anc)*100,3),
  AMR.p = round((freqtable$AMR/freqtable$sum.anc)*100,3),
  EAS.p = round((freqtable$EAS/freqtable$sum.anc)*100,3),
  EUR.p = round((freqtable$EUR/freqtable$sum.anc)*100,3) ))


# read fam (no header line: FID IID FATHER MOTHER SEX PHENO)
brids <- read.table('isa/isa.iid', comment.char = "", header = F, check.names = F)


# fill ancestries by individual ; last column = NA (viterbi)
sdchr <- data.frame()
for (i in seq(22)) {
  print(i)

  tempvit <- read.csv(paste0('outs.rfmix/isa.rfmix_chr',i,'.rfmix.1.Viterbi.txt'),
                              sep=' ',header=F, colClasses='numeric')

  freqind <- data.frame()
  for (j in seq(length(brids[,1])) ) { 
    tempcols <- tempvit[,c(2*(j-1)+1,2*(j-1)+2)]
    freqind <- as.data.frame(rbind(
      freqind,
      cbind(
        AFR = sum(rowSums(tempcols[,]==1)),
        AMR = sum(rowSums(tempcols[,]==2)),
        EAS = sum(rowSums(tempcols[,]==3)),
        EUR = sum(rowSums(tempcols[,]==4)),
        sum.anc = (dim(tempcols)[1])*2 ) ))
  }

  colnames(freqind) <- c("AFR","AMR","EAS","EUR","sum.anc")
  freqind <- as.data.frame(cbind(
    freqind,
    AFR.p = round((freqind$AFR/freqind$sum.anc)*100,3),
    AMR.p = round((freqind$AMR/freqind$sum.anc)*100,3),
    EAS.p = round((freqind$EAS/freqind$sum.anc)*100,3),
    EUR.p = round((freqind$EUR/freqind$sum.anc)*100,3) ))
  sdchr <- as.data.frame(rbind(
    sdchr,  cbind(
      chr = i,
      AFR.sd = round(sd(freqind$AFR.p),2),
      AMR.sd = round(sd(freqind$AMR.p),2),
      EAS.sd = round(sd(freqind$EAS.p),2),
      EUR.sd = round(sd(freqind$EUR.p),2) ) ) )
}


# join percentages with sd
freqtable <- as.data.frame(cbind( freqtable,sdchr ))


# bar all chromosomes
pchr1 <- reshape::melt(freqtable[,c(10,6:9)], id = 'chr')

pchr1$chr <- factor(pchr1$chr, levels=unique(pchr1$chr))
png(paste0('R_results/isa.MeanAncestryByChr.png'), width=300*2, height=300*2)
ggplot(pchr1, aes(y=value, x=chr,fill=variable)) +
  geom_bar(stat="identity")+
  scale_fill_manual(name='', labels=c("AFR","AMR","EAS","EUR"),
                    values=c("#778500","#ffd845","#018ead","#e03c31","#710027"))+
  theme_light()+ ylab('Ancestry contribution') + xlab('')+
  theme(axis.title.y=element_text(color="black",size=18),
        axis.text.x=element_text(size=10,angle=30),axis.text.y=element_text(size=11),
        legend.justification=c('right',"top"),legend.position='top',
        legend.text=element_text(size=12))
dev.off()


# ancestry proportions by individual
brids <- read.table('isa/isa.iid', comment.char = "", header = F, check.names = F)

for (i in seq(length(brids[,1]))) {
  assign(x = paste0('ind.',i), value = data.frame())   }


# fill ancestries by individual
for (i in seq(22)) {
  print(i)
  tempvit <- read.csv(paste0('outs.rfmix/isa.rfmix_chr',i,'.rfmix.1.Viterbi.txt'),
                              sep=' ',header=F, colClasses='numeric')
  tempmap <- read.csv(paste0('outs.rfmix/isa.rfmix_chr',i,'.map'),sep='\t',header=F)

  for (j in seq( length(brids[,1]) ) ) {
    
    tempcols <- tempvit[,c(2*(j-1)+1,2*(j-1)+2)]

    tempind <- data.frame()
    tempind <- as.data.frame(cbind(
      chr = rep(i,times=length(tempmap[,1])),
      bp = tempmap[,1],
      cM = tempmap[,2],
      rsID = tempmap[,3],
      crom.1 = tempcols[,1],
      crom.2 = tempcols[,2]  ))
    tempind <- as.data.frame(rbind(
      get(paste0('ind.',j)), tempind     ))
    assign(paste0('ind.',j), value = tempind)
  }
}



# 02 Ideograms and pie plot -----------------------------------------------

# Size of each chromossome
# https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37

human_karyotype <- as.data.frame(cbind(
  Chr=seq(1,22),
  size=c(249250621, 243199373, 198022430, 191154276, 180915260,
         171115067, 159138663, 146364022, 141213431, 135534747,
         135006516, 133851895, 115169878, 107349540, 102531392,
         90354753,  81195210,  78077248,  59128983,  63025520,
         48129895,  51304566) ))

testcolor <- c("#778500", "#ffd845","#018ead", "#b36c9e")
brids <- read.table('isa/isa.iid', comment.char = "", header = F, check.names = F)

set.seed(565465)
samplepositions <- sample(549066,20000)


# Ideogram
for (i in seq(length(brids[,1])) ) {

  ind0 <- data.frame()
  tempind <- data.frame()
  
  for (j in seq(22)) {
    
    tempvit <- read.csv(paste0('outs.rfmix/isa.rfmix_chr',j,'.rfmix.1.Viterbi.txt'),
                                sep=' ',header=F, colClasses='numeric')
    tempmap <- read.csv(paste0('outs.rfmix/isa.rfmix_chr',j,'.map'),sep='\t',header=F)

    tempcols <- tempvit[,c(2*(i-1)+1,2*(i-1)+2)]
    tempind <- as.data.frame(cbind(
      chr = rep(j,times=length(tempmap[,1])),
      bp = tempmap[,1],
      cM = tempmap[,2],
      rsID = tempmap[,3],
      crom.1 = tempcols[,1],
      crom.2 = tempcols[,2]  ))
    ind0 <- as.data.frame(rbind( ind0, tempind ))
  }
  ind0[,c(1,2,3)] <- apply(ind0[,c(1,2,3)], 2, as.numeric)

  inputkaryo <- reshape::melt(ind0[,c(1,2,5,6)], id = c('chr','bp'))
  inputkaryo$chrsize  <- human_karyotype$size[match(inputkaryo$chr, human_karyotype$Chr)]
  colnames(inputkaryo) <- c('chr', 'bp','cromatide','Ancestry','chrsize')
  x <- sum(inputkaryo$bp>inputkaryo$chrsize)
  print(paste0('SNPs outside chromosome boundaries: ',x))

  # Plot 1
  bp <- inputkaryo %>% ggplot(aes(x=chr, y = chrsize, fill=as.factor(cromatide))) +
    geom_bar(stat = "identity", position = "dodge", color="black", size=0.1) +
    scale_fill_manual(name='',  values=c("grey92", "grey92"))+
    scale_color_manual(values = testcolor, guide='none')+ 
    scale_x_continuous("", labels = as.character(unique(inputkaryo$chr)), 
                       breaks = unique(inputkaryo$chr)) +
    guides(fill="none")+    theme_classic()+   ylab("Genetic position (bp)")+
    theme(axis.line = element_blank(),axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, color = "#666666") ) +
    geom_segment(data=inputkaryo[inputkaryo$cromatide=='crom.1',][samplepositions,],
                 aes(x=chr-0.4,xend=chr-0.05, y=bp,yend=bp, colour = Ancestry ) )+
    geom_segment(data=inputkaryo[inputkaryo$cromatide=='crom.2',][samplepositions,],
                 aes(x=chr+0.05,xend=chr+0.4, y=bp,yend=bp, colour = Ancestry ) )


  # Pie plot
  piedata <- data.frame( round(100*(table(c(ind0[,5],ind0[,6]))/(2*nrow(ind0))),2)   )
  piedata$Var1 <- structure(as.factor(piedata$Var1), .Label = c("AFR","AMR","EAS","EUR"), class = "factor" )
  df2 <- piedata %>%  mutate(cs = rev(cumsum(rev(Freq))), 
           pos = Freq/2 + lead(cs, 1), pos = if_else(is.na(pos), Freq/2, pos))
  piecolor=testcolor[1:length(unique( c(ind0[,5],ind0[,6]) ))]

  # Plot 2
  pie <- ggplot(piedata, aes(x="", y=Freq, fill=as.factor(Var1)))+
    geom_bar(stat="identity", width=1)+
    coord_polar("y", start=0)+
    scale_fill_manual(values=piecolor )+
    geom_label_repel(aes(y = pos, label = paste0(Freq, "%")), fill = "white",
                     data = df2, size=4, show.legend = F, nudge_x = 1) +
    theme(axis.text.x=element_blank())+
    labs(x = NULL, y = NULL, fill = 'Ancestry', )+
    theme_classic()+ 
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), legend.title=element_text(size=12),
          plot.title = element_text(hjust = 0.5, color = "#666666"))
  bp + annotation_custom(  ggplotGrob(pie), xmin=10,xmax=22, ymin=1.25e+08, ymax=2.8e+08 )
}
