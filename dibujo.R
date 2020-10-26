library(gggenes)
library(ggplot2)
library(data.table)
library(dplyr)
library(grid)
library(purrr)
library(patchwork)

setwd("~/figure_synteny")

addi<-fread('gms2_comma.csv')
addi$molecule<-'G11C_GeneMark'
addi$direction<-ifelse(addi$direction=='+',1,-1)
addi$strand<-ifelse(addi$direction=='1','forward','reverse')
addi<-addi %>% mutate(gene = paste("gene",gene,sep=''))
addi2<-addi %>% select(molecule,start,end,direction,strand, gene)

datos<-fread('for_synteny_2.csv',data.table=FALSE)
datos<-rbind(datos, addi2) #adicion caso especial

meta<-fread('gene_info_R.csv',data.table=FALSE)
mg<-merge(x = datos, y = meta, by = "gene", all.x = TRUE)


dummies <- make_alignment_dummies(mg,aes(xmin = start, xmax = end, y = molecule, id = gene),on = "gene18")
dummies['direction']<-1
dummies['cato']<-'Core biosynthetic genes'
dummies['colo']<-'blue'

umi<-c("#E0FFFF","#1954A6","yellow","gray","orange","#34a617","purple")

#paste2 <- function(x, y, sep = "\n") paste(x, y, sep = sep)
#tik1<- meta$annotation[1:19] %>% reduce(paste2)
#tik2<- meta$annotation[20:38] %>% reduce(paste2)

#textframe <- data.frame(x = c(1, 1),labels = c(tik1, tik2),molecule="G11C")

tiff(file='synteny_1.21_comparative.tiff',units="in", width=12, height=6, res=250)
ggplot(mg, aes(xmin = start, xmax = end, y = molecule,fill=cato,forward=direction)) +
  geom_gene_arrow() +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  theme_genes()+
  theme(legend.title = element_blank())+
  geom_text(aes(x=(start+end)/2,y=molecule,label=num,vjust = c(-1.8)),size=1.5)+
  ylab('Strain')+
  xlab('Position (bp)')+
  scale_fill_manual(values=umi)
dev.off()

