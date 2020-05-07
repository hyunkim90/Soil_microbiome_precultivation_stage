### Correlation between OTU abundance and soil physicochemical properties without threshold
ps1.meta.all <- b.meta.all
bac.clean.rel <- microbiome::transform(bac.clean.ss, transform = "compositional")
otu.bac.clean.rel <- otu_table(bac.clean.rel)

otu.bac <- data.frame(otu_table(bac.clean.ss))
high_read=colSums(otu.bac)>=1000

head(otu.bac)

otu.bac$Total<-rowSums(otu.bac)

otu.bac.total <- otu.bac%>%select('Total')

otu.bac.total$RA <- otu.bac.total$Total/colSums(otu.bac.total)

otu.bac.ra <- otu.bac.total%>%select('RA')

tab_ra.b = otu.bac.ra

tab_c.b <- subset(tab_ra.b, RA > 0.0001) #### important cutoff, when is something counted as present ?! orig 0.001
length(tab_c.b$RA)


write.table(otu.bac.clean.rel, "otu.bac.clean.rel.tsv", sep = '\t', quote =F)
write.table(ps1.meta.all, "ps1.meta.all.tsv", sep = '\t', quote =F)

otu.bac.clean.rel<-read.table("otu.bac.clean.rel.tsv", sep = '\t', header =T)
ps1.meta.all<-read.table("ps1.meta.all.tsv", sep = '\t', header =T)

otu.bac.clean.rel.t <- subset(otu.bac.clean.rel, rownames(otu.bac.clean.rel) %in% rownames(tab_c))
head(t(otu.bac.clean.rel))
transpose.otu.bac.clean.rel<-t(otu.bac.clean.rel)
df.transpose.otu.bac.clean.rel <- data.frame(transpose.otu.bac.clean.rel)

ps1.meta.all.OTU <- merge(ps1.meta.all, transpose.otu.bac.clean.rel, by = "row.names")
rownames(ps1.meta.all.OTU) <- ps1.meta.all.OTU$Row.names

ps1.meta.all.OTU.t<-ps1.meta.all.OTU[-c(1:7)]

cor_5.b <- Hmisc::rcorr(as.matrix(ps1.meta.all.OTU.t), type="spearman")
M.b <- cor_5.b$r
p_mat.b <- cor_5.b$P

cor.bac.OTU <- flattenCorrMatrix(M.b,p_mat.b)
head(cor.bac.OTU)

#library(gtools)
#cor.bac.OTU$sig <- stars.pval(cor.bac.OTU$p)

physicochemicals<-colnames(ps1.meta.all.OTU.t[c(1:12)])
otunames.bac<-colnames(ps1.meta.all.OTU.t[-c(1:12)])

cor.bac.OTU.trim <- subset(cor.bac.OTU, row%in%physicochemicals)
cor.bac.OTU.trim <- subset(cor.bac.OTU.trim, column%in%otunames.bac)
head(cor.bac.OTU.trim)

write.table(cor.bac.OTU.trim, "cor.bac.OTU.trim.tsv", sep='\t', quote = F)




flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}




cor.bac.OTU.trim.sig.trans <- dcast(cor.bac.OTU.trim, row ~ column, value.var = "cor")
head(cor.bac.OTU.trim.sig.trans$row)

cor.bac.OTU.trim.sig.trans.p <- dcast(cor.bac.OTU.trim, row ~ column, value.var = "p")

write.xlsx(cor.bac.OTU.trim.sig.trans, "Relative abundance_Spearman correlation between soil and bacterial OTU.xlsx")
write.xlsx(cor.bac.OTU.trim.sig.trans.p, "Relative abundance_Spearman p-values for spearman correlation between soil and bacterial OTU.xlsx")


## Archaea
ps2.meta.all <- b.meta.all
arch.clean.rel <- microbiome::transform(arch.clean.ss, transform = "compositional")
otu.arch.clean.rel <- otu_table(arch.clean.rel)



otu.arch.clean.rel.t <- subset(otu.arch.clean.rel, rownames(otu.arch.clean.rel) %in% rownames(tab_c))
head(t(otu.arch.clean.rel))
transpose.otu.arch.clean.rel<-t(otu.arch.clean.rel)

ps2.meta.all.OTU <- merge(ps2.meta.all, transpose.otu.arch.clean.rel, by = "row.names")
rownames(ps2.meta.all.OTU) <- ps2.meta.all.OTU$Row.names
length(colnames(ps2.meta.all.OTU))

ps2.meta.all.OTU.t<-ps2.meta.all.OTU[-c(1:7)]

sapply(ps2.meta.all.OTU.t, class) ## soil properties are class 'factor'
indx.div <- sapply(ps2.meta.all.OTU.t, is.factor)
ps2.meta.all.OTU.t[indx.div] <- lapply(ps2.meta.all.OTU.t[indx.div], function(x) as.numeric(as.character(x)))
head(ps2.meta.all.OTU.t)
sapply(ps2.meta.all.OTU.t, class)

cor_5.a <- Hmisc::rcorr(as.matrix(ps2.meta.all.OTU.t), type="spearman")
M.a <- cor_5.a$r
p_mat.a <- cor_5.a$P

cor.arch.OTU <- flattenCorrMatrix(M.a,p_mat.a)
head(cor.arch.OTU)

library(gtools)
cor.arch.OTU$sig <- stars.pval(cor.arch.OTU$p)

physicochemicals<-colnames(ps2.meta.all.OTU.t[c(1:12)])
otunames.arch<-colnames(ps2.meta.all.OTU.t[-c(1:12)])

cor.arch.OTU.trim <- subset(cor.arch.OTU, row%in%physicochemicals)
cor.arch.OTU.trim <- subset(cor.arch.OTU.trim, column%in%otunames.arch)
head(cor.arch.OTU.trim)

write.table(cor.arch.OTU.trim, "cor.arch.OTU.trim.tsv", sep='\t', quote = F)

cor.arch.OTU.trim.sig.trans <- dcast(cor.arch.OTU.trim, row ~ column, value.var = "cor")
head(cor.arch.OTU.trim.sig.trans$row)

cor.arch.OTU.trim.sig.trans.p <- dcast(cor.arch.OTU.trim, row ~ column, value.var = "p")

write.xlsx(cor.arch.OTU.trim.sig.trans, "Relative abundance_Spearman correlation between soil and archaeal OTU_all.xlsx")
write.xlsx(cor.arch.OTU.trim.sig.trans.p, "Relative abundance_Spearman p-values for spearman correlation between soil and archaeal OTU_all.xlsx")



## Fungi
ps3.meta.all <- f.meta.all
fun.clean.rel <- microbiome::transform(fun.clean.ss, transform = "compositional")
otu.fun.clean.rel <- otu_table(fun.clean.rel)

otu.fun <- data.frame(otu_table(fun.clean.ss))
high_read=colSums(otu.fun)>=1000

head(otu.fun)

otu.fun$Total<-rowSums(otu.fun)

otu.fun.total <- otu.fun%>%select('Total')

otu.fun.total$RA <- otu.fun.total$Total/colSums(otu.fun.total)

otu.fun.ra <- otu.fun.total%>%select('RA')

tab_ra = otu.fun.ra

tab_c <- subset(tab_ra, RA > 0.0001) #### important cutoff, when is something counted as present ?! orig 0.001
length(tab_c$RA)


otu.fun.clean.rel.t <- subset(otu.fun.clean.rel, rownames(otu.fun.clean.rel) %in% rownames(tab_c))
head(t(otu.fun.clean.rel))
transpose.otu.fun.clean.rel<-t(otu.fun.clean.rel)

ps3.meta.all.OTU <- merge(ps3.meta.all, transpose.otu.fun.clean.rel, by = "row.names")
rownames(ps3.meta.all.OTU) <- ps3.meta.all.OTU$Row.names
length(colnames(ps3.meta.all.OTU))

ps3.meta.all.OTU.t<-ps3.meta.all.OTU[-c(1:7)]

sapply(ps3.meta.all.OTU.t, class) ## soil properties are class 'factor'
indx.div <- sapply(ps3.meta.all.OTU.t, is.factor)
ps3.meta.all.OTU.t[indx.div] <- lapply(ps3.meta.all.OTU.t[indx.div], function(x) as.numeric(as.character(x)))
head(ps3.meta.all.OTU.t)
sapply(ps3.meta.all.OTU.t, class)

cor_5.f <- Hmisc::rcorr(as.matrix(ps3.meta.all.OTU.t), type="spearman")
M.f <- cor_5.f$r
p_mat.f <- cor_5.f$P

cor.fun.OTU <- flattenCorrMatrix(M.f,p_mat.f)
head(cor.fun.OTU)

library(gtools)
cor.fun.OTU$sig <- stars.pval(cor.fun.OTU$p)

physicochemicals<-colnames(ps3.meta.all.OTU.t[c(1:12)])
otunames.fun<-colnames(ps3.meta.all.OTU.t[-c(1:12)])

cor.fun.OTU.trim <- subset(cor.fun.OTU, row%in%physicochemicals)
cor.fun.OTU.trim <- subset(cor.fun.OTU.trim, column%in%otunames.fun)
head(cor.fun.OTU.trim)

write.table(cor.fun.OTU.trim, "cor.fun.OTU.trim.tsv", sep='\t', quote = F)

cor.fun.OTU.trim.sig.trans <- dcast(cor.fun.OTU.trim, row ~ column, value.var = "cor")
head(cor.fun.OTU.trim.sig.trans$row)

cor.fun.OTU.trim.sig.trans.p <- dcast(cor.fun.OTU.trim, row ~ column, value.var = "p")

write.xlsx(cor.fun.OTU.trim.sig.trans, "Relative abundance_Spearman correlation between soil and fungal OTU_all.xlsx")
write.xlsx(cor.fun.OTU.trim.sig.trans.p, "Relative abundance_Spearman p-values for spearman correlation between soil and fungal OTU_all.xlsx")



### Volcano plot
cor.bac.OTU.trim <- read.table("cor.bac.OTU.trim.tsv", sep = '\t', header = T)
cor.arch.OTU.trim <- read.table("cor.arch.OTU.trim.tsv", sep = '\t', header = T)
cor.fun.OTU.trim <- read.table("cor.fun.OTU.trim.tsv", sep = '\t', header = T)

cor.bac.OTU.trim.sig <- subset(cor.bac.OTU.trim, p < 0.05)
cor.bac.OTU.trim.sig.SOM <- subset(cor.bac.OTU.trim.sig, row == "SOM")

length(cor.bac.OTU.trim.sig.SOM$column)

length(cor.bac.OTU.trim.sig.SOM$column[cor.bac.OTU.trim.sig.SOM$cor<0]) #304 #putative_oligotroph #835
length(cor.bac.OTU.trim.sig.SOM$column[cor.bac.OTU.trim.sig.SOM$cor>0]) #238 #putative_copiotroph #775


cor.bac.OTU.trim.sig <- subset(cor.bac.OTU.trim, p < 0.05)
cor.bac.OTU.trim.sig.TN <- subset(cor.bac.OTU.trim.sig, row == "TN")

length(cor.bac.OTU.trim.sig.TN$column)

length(cor.bac.OTU.trim.sig.TN$column[cor.bac.OTU.trim.sig.TN$cor<0]) #277 #putative_oligotroph #976
length(cor.bac.OTU.trim.sig.TN$column[cor.bac.OTU.trim.sig.TN$cor>0]) #255 #putative_copiotroph #690


tab_ra.b
tab_ra.a
tab_ra.f



cor.arch.OTU.trim.sig <- subset(cor.arch.OTU.trim, p < 0.05)
cor.arch.OTU.trim.sig.SOM <- subset(cor.arch.OTU.trim.sig, row == "SOM")

length(cor.arch.OTU.trim.sig.SOM$column)

length(cor.arch.OTU.trim.sig.SOM$column[cor.arch.OTU.trim.sig.SOM$cor<0]) #34 #putative_oligotroph #42
length(cor.arch.OTU.trim.sig.SOM$column[cor.arch.OTU.trim.sig.SOM$cor>0]) #36 #putative_copiotroph #71

cor.arch.OTU.trim.sig <- subset(cor.arch.OTU.trim, p < 0.05)
cor.arch.OTU.trim.sig.TN <- subset(cor.arch.OTU.trim.sig, row == "TN")

length(cor.arch.OTU.trim.sig.TN$column)

length(cor.arch.OTU.trim.sig.TN$column[cor.arch.OTU.trim.sig.TN$cor<0]) #33 #putative_oligotroph#47
length(cor.arch.OTU.trim.sig.TN$column[cor.arch.OTU.trim.sig.TN$cor>0]) #34 #putative_copiotroph#51





cor.fun.OTU.trim.sig <- subset(cor.fun.OTU.trim, p < 0.05)
cor.fun.OTU.trim.sig.SOM <- subset(cor.fun.OTU.trim.sig, row == "SOM")

length(cor.fun.OTU.trim.sig.SOM$column)

length(cor.fun.OTU.trim.sig.SOM$column[cor.fun.OTU.trim.sig.SOM$cor<0]) #97 #putative_oligotroph #495
length(cor.fun.OTU.trim.sig.SOM$column[cor.fun.OTU.trim.sig.SOM$cor>0]) #123 #putative_copiotroph #538

cor.fun.OTU.trim.sig <- subset(cor.fun.OTU.trim, p < 0.05)
cor.fun.OTU.trim.sig.TN <- subset(cor.fun.OTU.trim.sig, row == "TN")

length(cor.fun.OTU.trim.sig.TN$column)

length(cor.fun.OTU.trim.sig.TN$column[cor.fun.OTU.trim.sig.TN$cor<0]) #66 #putative_oligotroph #507
length(cor.fun.OTU.trim.sig.TN$column[cor.fun.OTU.trim.sig.TN$cor>0]) #98 #putative_copiotroph #501



## Distribution of putative oligo and copiotrophs in Fertilized field and non fertilized field
#or oligotrophic and copiotrophic environment

tab_ra.f.SOM.copi <- subset(tab_ra.f, rownames(tab_ra.f)%in% cor.fun.OTU.trim.sig.SOM$column[cor.fun.OTU.trim.sig.SOM$cor>0])
tab_ra.f.TN.copi <- subset(tab_ra.f, rownames(tab_ra.f)%in% cor.fun.OTU.trim.sig.TN$column[cor.fun.OTU.trim.sig.TN$cor>0])
sum(tab_ra.f.SOM.copi$RA) #0.3434138 #0.3507278
sum(tab_ra.f.TN.copi$RA)  #0.153926 #0.1604727

common.fun.copi<- intersect(rownames(tab_ra.f.SOM.copi), rownames(tab_ra.f.TN.copi))
TN_specific.copi <- subset(tab_ra.f.TN.copi, !rownames(tab_ra.f.TN.copi) %in% common.fun.copi)
length(TN_specific.copi$RA)#160
tab_ra.f.SOM.TN.copi <- rbind(tab_ra.f.SOM.copi, TN_specific.copi)
length(tab_ra.f.SOM.TN.copi$RA)
sum(tab_ra.f.SOM.TN.copi$RA) #0.3633143
length(tab_ra.f.SOM.TN.copi$RA) #698

tab_ra.f.SOM.oligo <- subset(tab_ra.f, rownames(tab_ra.f)%in% cor.fun.OTU.trim.sig.SOM$column[cor.fun.OTU.trim.sig.SOM$cor<0])
tab_ra.f.TN.oligo <- subset(tab_ra.f, rownames(tab_ra.f)%in% cor.fun.OTU.trim.sig.TN$column[cor.fun.OTU.trim.sig.TN$cor<0])
sum(tab_ra.f.SOM.oligo$RA) #0.3076559 #0.312305
sum(tab_ra.f.TN.oligo$RA)  #0.1223175 #0.1282304

common.fun.oligo<-intersect(rownames(tab_ra.f.SOM.oligo), rownames(tab_ra.f.TN.oligo))
TN_specific.oligo <- subset(tab_ra.f.TN.oligo, !rownames(tab_ra.f.TN.oligo) %in% common.fun.oligo)
length(TN_specific.oligo$RA)#162
tab_ra.f.SOM.TN.oligo <- rbind(tab_ra.f.SOM.oligo,TN_specific.oligo)
sum(tab_ra.f.SOM.TN.oligo$RA) #0.3258077
length(tab_ra.f.SOM.TN.oligo$RA) #657


##Bacteria
tab_ra.b.SOM.copi <- subset(tab_ra.b, rownames(tab_ra.b)%in% cor.bac.OTU.trim.sig.SOM$column[cor.bac.OTU.trim.sig.SOM$cor>0])
tab_ra.b.TN.copi <- subset(tab_ra.b, rownames(tab_ra.b)%in% cor.bac.OTU.trim.sig.TN$column[cor.bac.OTU.trim.sig.TN$cor>0])
sum(tab_ra.b.SOM.copi$RA) #0.2077602
sum(tab_ra.b.TN.copi$RA)  #0.2408266

common.bac.copi<-intersect(rownames(tab_ra.b.SOM.copi), rownames(tab_ra.b.TN.copi))
TN_specific.copi <- subset(tab_ra.b.TN.copi, !rownames(tab_ra.b.TN.copi) %in% common.bac.copi)
length(TN_specific.copi$RA)#179
tab_ra.b.SOM.TN.copi <- rbind(tab_ra.b.SOM.copi,TN_specific.copi)
sum(tab_ra.b.SOM.TN.copi$RA) # 0.2562043
length(tab_ra.b.SOM.TN.copi$RA) #954

tab_ra.b.SOM.oligo <- subset(tab_ra.b, rownames(tab_ra.b)%in% cor.bac.OTU.trim.sig.SOM$column[cor.bac.OTU.trim.sig.SOM$cor<0])
tab_ra.b.TN.oligo <- subset(tab_ra.b, rownames(tab_ra.b)%in% cor.bac.OTU.trim.sig.TN$column[cor.bac.OTU.trim.sig.TN$cor<0])
sum(tab_ra.b.SOM.oligo$RA) # 0.2237516
sum(tab_ra.b.TN.oligo$RA)  #0.2001819

common.bac.oligo<-intersect(rownames(tab_ra.b.SOM.oligo), rownames(tab_ra.b.TN.oligo))
TN_specific.oligo <- subset(tab_ra.b.TN.oligo, !rownames(tab_ra.b.TN.oligo) %in% common.bac.oligo)
length(TN_specific.oligo$RA)#285
tab_ra.b.SOM.TN.oligo <- rbind(tab_ra.b.SOM.oligo,TN_specific.oligo)
sum(tab_ra.b.SOM.TN.oligo$RA) #0.2410767
length(tab_ra.b.SOM.TN.oligo$RA) #1120

##Archaea
tab_ra.a.SOM.copi <- subset(tab_ra.a, rownames(tab_ra.a)%in% cor.arch.OTU.trim.sig.SOM$column[cor.arch.OTU.trim.sig.SOM$cor>0])
tab_ra.a.TN.copi <- subset(tab_ra.a, rownames(tab_ra.a)%in% cor.arch.OTU.trim.sig.TN$column[cor.arch.OTU.trim.sig.TN$cor>0])
sum(tab_ra.a.SOM.copi$RA) #0.2614325
sum(tab_ra.a.TN.copi$RA)  #0.2694488

common.arch.copi<-intersect(rownames(tab_ra.a.SOM.copi), rownames(tab_ra.a.TN.copi))
TN_specific.copi <- subset(tab_ra.a.TN.copi, !rownames(tab_ra.a.TN.copi) %in% common.arch.copi)
length(TN_specific.copi$RA)#11
tab_ra.a.SOM.TN.copi <- rbind(tab_ra.a.SOM.copi,TN_specific.copi)
sum(tab_ra.a.SOM.TN.copi$RA) #0.3204761
length(tab_ra.a.SOM.TN.copi$RA) #81


tab_ra.a.SOM.oligo <- subset(tab_ra.a, rownames(tab_ra.a)%in% cor.arch.OTU.trim.sig.SOM$column[cor.arch.OTU.trim.sig.SOM$cor<0])
tab_ra.a.TN.oligo <- subset(tab_ra.a, rownames(tab_ra.a)%in% cor.arch.OTU.trim.sig.TN$column[cor.arch.OTU.trim.sig.TN$cor<0])
sum(tab_ra.a.SOM.oligo$RA) # 0.30229
sum(tab_ra.a.TN.oligo$RA)  #0.3187515

common.arch.oligo<-intersect(rownames(tab_ra.a.SOM.oligo), rownames(tab_ra.a.TN.oligo))
TN_specific.oligo <- subset(tab_ra.a.TN.oligo, !rownames(tab_ra.a.TN.oligo) %in% common.arch.oligo)
length(TN_specific.oligo$RA)#11
tab_ra.a.SOM.TN.oligo <- rbind(tab_ra.a.SOM.oligo,TN_specific.oligo)
sum(tab_ra.a.SOM.TN.oligo$RA) #0.3305097
length(tab_ra.a.SOM.TN.oligo$RA) #53



write.table(tab_ra.f.SOM.TN.oligo, "tab_ra.f.SOM.TN.oligo.tsv", sep= '\t', quote=F)
write.table(tab_ra.f.SOM.TN.copi, "tab_ra.f.SOM.TN.copi.tsv", sep= '\t', quote=F)

write.table(tab_ra.a.SOM.TN.oligo, "tab_ra.a.SOM.TN.oligo.tsv", sep= '\t', quote=F)
write.table(tab_ra.a.SOM.TN.copi, "tab_ra.a.SOM.TN.copi.tsv", sep= '\t', quote=F)

write.table(tab_ra.b.SOM.TN.oligo, "tab_ra.b.SOM.TN.oligo.tsv", sep= '\t', quote=F)
write.table(tab_ra.b.SOM.TN.copi, "tab_ra.b.SOM.TN.copi.tsv", sep= '\t', quote=F)


tab_ra.b
cor.bac.OTU.trim.sig.SOM
rownames(cor.bac.OTU.trim.sig.SOM) <- cor.bac.OTU.trim.sig.SOM$column
bac.RA.Corr<-merge(tab_ra.b, cor.bac.OTU.trim.sig.SOM, by = "row.names")
rownames(bac.RA.Corr) <- bac.RA.Corr$Row.names
bac.RA.Corr<- bac.RA.Corr[-c(1)]
bac.RA.Corr$log2_RA <- log2(bac.RA.Corr$RA)
bac.RA.Corr$log10_p <- log10(bac.RA.Corr$p)

## Visualization
theme_set(theme_bw())
library(ggrepel)

ggplot(bac.RA.Corr, aes(x=cor, y=log2_RA)) +
  xlab('\n Correlation coefficient')+
  ylab("log2 mean RA\n") +
  geom_point(aes(size =RA),  alpha=0.6) +
  #scale_color_manual(labels = c('Wild'="Wild enriched\n(FDR Q<0.01 & log2FC>2)     ",'Domesticated'= "Domesticated enriched\n(FDR Q<0.01 & log2FC>2)     ",'Middle'='log2FC<|2|     ','Bottom'='NS\n(FDR Q>0.01)     '), values = c("Wild"= "red3", 'Domesticated'='forestgreen',"Middle"="lightblue","Bottom"= "yellow3"))+
  
  #ggtitle("Volcano Plot \n") +
  #theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #geom_hline(yintercept=h, color="maroon4")+
  #geom_vline(xintercept=0.2, color="maroon4")+
  #geom_vline(xintercept=0, color="black")+
  #geom_vline(xintercept=-0.2, color="maroon4")+
  #theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  #scale_x_continuous(breaks=seq(-20,20,1))+
  #scale_y_continuous(breaks=seq(0,60,10))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
#geom_text_repel(data=subset(sub.mydata, volcano_y > h),aes(label=id), size=3, nudge_x= 1)+
#geom_text_repel(data=subset(sub.mydata, volcano_y <= h),aes(label=id), direction = "both", nudge_y =-2, size=3)



dev.off()


tab_ra.a
cor.arch.OTU.trim.sig.SOM
rownames(cor.arch.OTU.trim.sig.SOM) <- cor.arch.OTU.trim.sig.SOM$column
arch.RA.Corr<-merge(tab_ra.a, cor.arch.OTU.trim.sig.SOM, by = "row.names")
rownames(arch.RA.Corr) <- arch.RA.Corr$Row.names
arch.RA.Corr<- arch.RA.Corr[-c(1)]
arch.RA.Corr$log2_RA <- log2(arch.RA.Corr$RA)

## Visualization
theme_set(theme_bw())
library(ggrepel)

ggplot(arch.RA.Corr, aes(x=cor, y=log2_RA)) +
  xlab('\n Correlation coefficient')+
  ylab("log2 mean RA\n") +
  geom_point(aes(size = RA),  alpha=0.6) +
  #scale_color_manual(labels = c('Wild'="Wild enriched\n(FDR Q<0.01 & log2FC>2)     ",'Domesticated'= "Domesticated enriched\n(FDR Q<0.01 & log2FC>2)     ",'Middle'='log2FC<|2|     ','Bottom'='NS\n(FDR Q>0.01)     '), values = c("Wild"= "red3", 'Domesticated'='forestgreen',"Middle"="lightblue","Bottom"= "yellow3"))+
  
  #ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #geom_hline(yintercept=h, color="maroon4")+
  #geom_vline(xintercept=0.2, color="maroon4")+
  #geom_vline(xintercept=0, color="black")+
  #geom_vline(xintercept=-0.2, color="maroon4")+
  #theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  #scale_x_continuous(breaks=seq(-0.7,0.7,0.1))+
  #scale_y_continuous(breaks=seq(0,60,10))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
#geom_text_repel(data=subset(sub.mydata, volcano_y > h),aes(label=id), size=3, nudge_x= 1)+
#geom_text_repel(data=subset(sub.mydata, volcano_y <= h),aes(label=id), direction = "both", nudge_y =-2, size=3)



dev.off()



tab_ra.f
cor.fun.OTU.trim.sig.SOM
rownames(cor.fun.OTU.trim.sig.SOM) <- cor.fun.OTU.trim.sig.SOM$column
fun.RA.Corr<-merge(tab_ra.f, cor.fun.OTU.trim.sig.SOM, by = "row.names")
rownames(fun.RA.Corr) <- fun.RA.Corr$Row.names
fun.RA.Corr<- fun.RA.Corr[-c(1)]
fun.RA.Corr$log2_RA <- log2(fun.RA.Corr$RA)

## Visualization
theme_set(theme_bw())
library(ggrepel)

ggplot(fun.RA.Corr, aes(x=cor, y=log2_RA)) +
  xlab('\n Correlation coefficient')+
  ylab("log2 mean RA\n") +
  geom_point(aes(size = RA),  alpha=0.6) +
  #scale_color_manual(labels = c('Wild'="Wild enriched\n(FDR Q<0.01 & log2FC>2)     ",'Domesticated'= "Domesticated enriched\n(FDR Q<0.01 & log2FC>2)     ",'Middle'='log2FC<|2|     ','Bottom'='NS\n(FDR Q>0.01)     '), values = c("Wild"= "red3", 'Domesticated'='forestgreen',"Middle"="lightblue","Bottom"= "yellow3"))+
  
  #ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #geom_hline(yintercept=h, color="maroon4")+
  #geom_vline(xintercept=0.2, color="maroon4")+
  #geom_vline(xintercept=0, color="black")+
  #geom_vline(xintercept=-0.2, color="maroon4")+
  #theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.fes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  #scale_x_continuous(breaks=seq(-0.7,0.7,0.1))+
  #scale_y_continuous(breaks=seq(0,60,10))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
#geom_text_repel(data=subset(sub.mydata, volcano_y > h),aes(label=id), size=3, nudge_x= 1)+
#geom_text_repel(data=subset(sub.mydata, volcano_y <= h),aes(label=id), direction = "both", nudge_y =-2, size=3)



dev.off()

bac.clean.ss.field <- merge_samples(bac.clean.ss, "Field")
bac.clean.ss.field.rel <- microbiome::transform(bac.clean.ss.field, transform = "compositional")
sample_data(bac.clean.ss.field.rel)
bac.clean.ss.field.rel.JJ1.UF2 <- subset_samples(bac.clean.ss.field.rel, Field %in% c(17, 25))
bac.clean.ss.field.rel.JJ1.UF2.t <- phyloseq::filter_taxa(bac.clean.ss.field.rel.JJ1.UF2, function(x) sum(x) != 0, TRUE)
otu.bac.clean.ss.field.rel.JJ1.UF2.t <- t(otu_table(bac.clean.ss.field.rel.JJ1.UF2.t))
otu.bac.clean.ss.field.rel.JJ1.UF2.t.SOM <- subset(otu.bac.clean.ss.field.rel.JJ1.UF2.t, rownames(otu.bac.clean.ss.field.rel.JJ1.UF2.t)%in%bac.RA.Corr$column)


otu.bac.clean.ss.field.rel.JJ1.UF2.t.SOM.oligo <- subset(otu.bac.clean.ss.field.rel.JJ1.UF2.t, rownames(otu.bac.clean.ss.field.rel.JJ1.UF2.t)%in%rownames(tab_ra.b.SOM.TN.oligo))
tab_bac.oligo<-otu.bac.clean.ss.field.rel.JJ1.UF2.t.SOM.oligo
df.tab_bac.oligo <- data.frame(tab_bac.oligo)

head(df.tab_bac.oligo)
tab_bac.oligo$OTU <- rownames(tab_bac.oligo)

sum(df.tab_bac.oligo$JJ1) # 0.314945
sum(df.tab_bac.oligo$UF2) # 0.356419



### IS and NJ
sample_data(bac.clean.ss.field.rel)
bac.clean.ss.field.rel.IS1.NJ1 <- subset_samples(bac.clean.ss.field.rel, Field %in% c(13, 21))
bac.clean.ss.field.rel.IS1.NJ1.t <- phyloseq::filter_taxa(bac.clean.ss.field.rel.IS1.NJ1, function(x) sum(x) != 0, TRUE)
otu.bac.clean.ss.field.rel.IS1.NJ1.t <- t(otu_table(bac.clean.ss.field.rel.IS1.NJ1.t))

otu.bac.clean.ss.field.rel.IS1.NJ1.t.SOM.oligo <- subset(otu.bac.clean.ss.field.rel.IS1.NJ1.t, rownames(otu.bac.clean.ss.field.rel.IS1.NJ1.t)%in%rownames(tab_ra.b.SOM.TN.oligo))
tab_bac.oligo.2<-otu.bac.clean.ss.field.rel.IS1.NJ1.t.SOM.oligo
df.tab_bac.oligo.2 <- data.frame(tab_bac.oligo.2)

sum(df.tab_bac.oligo.2$IS1) #0.1139649
sum(df.tab_bac.oligo.2$NJ1) #0.1684452




## 
otu.bac.clean.rel.SOM.oligo <- subset(otu.bac.clean.rel, rownames(otu.bac.clean.rel)%in%rownames(tab_ra.b.SOM.TN.oligo))
df.otu.bac.clean.rel.SOM.oligo <- data.frame(otu.bac.clean.rel.SOM.oligo)

df.otu.bac.clean.rel.SOM.oligo.colsum<-colSums(df.otu.bac.clean.rel.SOM.oligo, na.rm = FALSE, dims = 1)
df.bac.oligo.colsum <- data.frame(df.otu.bac.clean.rel.SOM.oligo.colsum)
names(df.bac.oligo.colsum)[1] <- "bac_RA_oligo"


otu.bac.clean.rel.SOM.copi <- subset(otu.bac.clean.rel, rownames(otu.bac.clean.rel)%in%rownames(tab_ra.b.SOM.TN.copi))
df.otu.bac.clean.rel.SOM.copi <- data.frame(otu.bac.clean.rel.SOM.copi)

df.otu.bac.clean.rel.SOM.copi.colsum<-colSums(df.otu.bac.clean.rel.SOM.copi, na.rm = FALSE, dims = 1)
df.bac.copio.colsum <- data.frame(df.otu.bac.clean.rel.SOM.copi.colsum)
names(df.bac.copio.colsum)[1] <- "bac_RA_copio"


ps1.meta.trophic <- cbind(ps1.meta, df.bac.oligo.colsum, df.bac.copio.colsum)

##Archaea
otu.arch.clean.rel
otu.arch.clean.rel.SOM.oligo <- subset(otu.arch.clean.rel, rownames(otu.arch.clean.rel)%in%rownames(tab_ra.a.SOM.TN.oligo))
df.otu.arch.clean.rel.SOM.oligo <- data.frame(otu.arch.clean.rel.SOM.oligo)

df.otu.arch.clean.rel.SOM.oligo.colsum<-colSums(df.otu.arch.clean.rel.SOM.oligo, na.rm = FALSE, dims = 1)
df.arch.oligo.colsum <- data.frame(df.otu.arch.clean.rel.SOM.oligo.colsum)
names(df.arch.oligo.colsum)[1] <- "arch_RA_oligo"


otu.arch.clean.rel.SOM.copi <- subset(otu.arch.clean.rel, rownames(otu.arch.clean.rel)%in%rownames(tab_ra.a.SOM.TN.copi))
df.otu.arch.clean.rel.SOM.copi <- data.frame(otu.arch.clean.rel.SOM.copi)

df.otu.arch.clean.rel.SOM.copi.colsum<-colSums(df.otu.arch.clean.rel.SOM.copi, na.rm = FALSE, dims = 1)
df.arch.copio.colsum <- data.frame(df.otu.arch.clean.rel.SOM.copi.colsum)
names(df.arch.copio.colsum)[1] <- "arch_RA_copio"

ps2.meta.trophic <- cbind(ps2.meta, df.arch.oligo.colsum, df.arch.copio.colsum)


##Fungi
otu.fun.clean.rel
otu.fun.clean.rel.SOM.oligo <- subset(otu.fun.clean.rel, rownames(otu.fun.clean.rel)%in%rownames(tab_ra.f.SOM.TN.oligo))
df.otu.fun.clean.rel.SOM.oligo <- data.frame(otu.fun.clean.rel.SOM.oligo)

df.otu.fun.clean.rel.SOM.oligo.colsum<-colSums(df.otu.fun.clean.rel.SOM.oligo, na.rm = FALSE, dims = 1)
df.fun.oligo.colsum <- data.frame(df.otu.fun.clean.rel.SOM.oligo.colsum)
names(df.fun.oligo.colsum)[1] <- "fun_RA_oligo"


otu.fun.clean.rel.SOM.copi <- subset(otu.fun.clean.rel, rownames(otu.fun.clean.rel)%in%rownames(tab_ra.f.SOM.TN.copi))
df.otu.fun.clean.rel.SOM.copi <- data.frame(otu.fun.clean.rel.SOM.copi)

df.otu.fun.clean.rel.SOM.copi.colsum<-colSums(df.otu.fun.clean.rel.SOM.copi, na.rm = FALSE, dims = 1)
df.fun.copio.colsum <- data.frame(df.otu.fun.clean.rel.SOM.copi.colsum)
names(df.fun.copio.colsum)[1] <- "fun_RA_copio"

ps3.meta.trophic <- cbind(ps3.meta, df.fun.oligo.colsum, df.fun.copio.colsum)

head(ps3.meta.trophic)
### boxplot

##Fungi

##Oligo
ps3.meta.trophic.CJYSMYNJ <- subset(ps3.meta.trophic, Field %in% c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
ps3.meta.trophic.CJYSMYNJ$Field <- factor(ps3.meta.trophic.CJYSMYNJ$Field, levels = c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
test <- aggregate(ps3.meta.trophic.CJYSMYNJ$fun_RA_oligo, by = list(ps3.meta.trophic.CJYSMYNJ$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(fun_RA_oligo ~ Field, data = ps3.meta.trophic.CJYSMYNJ)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(fun_RA_oligo ~ Field,
              data=ps3.meta.trophic.CJYSMYNJ,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.meta.trophic.CJYSMYNJ, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
ps3.meta.trophic.CJYSMYNJ <- subset(ps3.meta.trophic, Field %in% c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
ps3.meta.trophic.CJYSMYNJ$Field <- factor(ps3.meta.trophic.CJYSMYNJ$Field, levels = c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
test <- aggregate(ps3.meta.trophic.CJYSMYNJ$fun_RA_copio, by = list(ps3.meta.trophic.CJYSMYNJ$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(fun_RA_copio ~ Field, data = ps3.meta.trophic.CJYSMYNJ)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(fun_RA_copio ~ Field,
              data=ps3.meta.trophic.CJYSMYNJ,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.meta.trophic.CJYSMYNJ, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Bacteria

##Oligo
ps1.meta.trophic.CJYSMYNJ <- subset(ps1.meta.trophic, Field %in% c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
ps1.meta.trophic.CJYSMYNJ$Field <- factor(ps1.meta.trophic.CJYSMYNJ$Field, levels = c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
test <- aggregate(ps1.meta.trophic.CJYSMYNJ$bac_RA_oligo, by = list(ps1.meta.trophic.CJYSMYNJ$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(bac_RA_oligo ~ Field, data = ps1.meta.trophic.CJYSMYNJ)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(bac_RA_oligo ~ Field,
              data=ps1.meta.trophic.CJYSMYNJ,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.meta.trophic.CJYSMYNJ, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
ps1.meta.trophic.CJYSMYNJ <- subset(ps1.meta.trophic, Field %in% c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
ps1.meta.trophic.CJYSMYNJ$Field <- factor(ps1.meta.trophic.CJYSMYNJ$Field, levels = c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
test <- aggregate(ps1.meta.trophic.CJYSMYNJ$bac_RA_copio, by = list(ps1.meta.trophic.CJYSMYNJ$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(bac_RA_copio ~ Field, data = ps1.meta.trophic.CJYSMYNJ)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(bac_RA_copio ~ Field,
              data=ps1.meta.trophic.CJYSMYNJ,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.meta.trophic.CJYSMYNJ, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


##Archaea

##Oligo
ps2.meta.trophic.CJYSMYNJ <- subset(ps2.meta.trophic, Field %in% c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
ps2.meta.trophic.CJYSMYNJ$Field <- factor(ps2.meta.trophic.CJYSMYNJ$Field, levels = c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
test <- aggregate(ps2.meta.trophic.CJYSMYNJ$arch_RA_oligo, by = list(ps2.meta.trophic.CJYSMYNJ$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(arch_RA_oligo ~ Field, data = ps2.meta.trophic.CJYSMYNJ)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(arch_RA_oligo ~ Field,
              data=ps2.meta.trophic.CJYSMYNJ,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.meta.trophic.CJYSMYNJ, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
ps2.meta.trophic.CJYSMYNJ <- subset(ps2.meta.trophic, Field %in% c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
ps2.meta.trophic.CJYSMYNJ$Field <- factor(ps2.meta.trophic.CJYSMYNJ$Field, levels = c("CJ1", "CJ2", "YS1", "YS2", "MY1", "MY2", "NJ1", "NJ2"))
test <- aggregate(ps2.meta.trophic.CJYSMYNJ$arch_RA_copio, by = list(ps2.meta.trophic.CJYSMYNJ$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(arch_RA_copio ~ Field, data = ps2.meta.trophic.CJYSMYNJ)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(arch_RA_copio ~ Field,
              data=ps2.meta.trophic.CJYSMYNJ,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.meta.trophic.CJYSMYNJ, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p




### For barpolt
melted.arch.clean.rel <- psmelt(arch.clean.rel)


lb <- c('CC1',
        'CC2',
        'CJ1',
        'CJ2',
        'DG1',
        'DG2',
        'JJ1',
        'JJ2',
        'MY1',
        'MY2',
        'NJ1',
        'NJ2',
        'IS1',
        'IS2',
        'YS1',
        'YS2',
        'UF1',
        'UF2',
        'CC1.18',
        'CC2.18',
        'CJ1.18',
        'CJ2.18',
        'DG1.18',
        'DG2.18',
        'IS1.18',
        'IS2.18',
        'UF1.18',
        'UF2.18')

melted.arch.clean.rel$Field <- factor(melted.arch.clean.rel$Field, levels = lb)

melted.arch.clean.rel.oligo <- subset(melted.arch.clean.rel, OTU %in% rownames(tab_ra.a.SOM.TN.oligo))
melted.arch.clean.rel.oligo$Trophic <- "Oligo"

melted.arch.clean.rel.copi <- subset(melted.arch.clean.rel, OTU %in% rownames(tab_ra.a.SOM.TN.copi))
melted.arch.clean.rel.copi$Trophic <- "Copio"


mer.arch.clean.rel.trophic <- rbind(melted.arch.clean.rel.oligo, melted.arch.clean.rel.copi)
mer.arch.clean.rel.trophic$Abundance



df.trophic.rel.p1 <- ggplot(mer.arch.clean.rel.trophic, aes(x=Field, y = Abundance/9, fill = Trophic)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("#669933", "#CCCC99")) +
  xlab('')+
  ylab("RA\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.trophic.rel.p1



##Fungi
melted.fun.clean.rel <- psmelt(fun.clean.rel)


lb <- c('CC1',
        'CC2',
        'CJ1',
        'CJ2',
        'DG1',
        'DG2',
        'JJ1',
        'JJ2',
        'MY1',
        'MY2',
        'NJ1',
        'NJ2',
        'IS1',
        'IS2',
        'YS1',
        'YS2',
        'UF1',
        'UF2',
        'CC1.18',
        'CC2.18',
        'CJ1.18',
        'CJ2.18',
        'DG1.18',
        'DG2.18',
        'IS1.18',
        'IS2.18',
        'UF1.18',
        'UF2.18')

melted.fun.clean.rel$Field <- factor(melted.fun.clean.rel$Field, levels = lb)

melted.fun.clean.rel.oligo <- subset(melted.fun.clean.rel, OTU %in% rownames(tab_ra.f.SOM.TN.oligo))
melted.fun.clean.rel.oligo$Trophic <- "Oligo"

melted.fun.clean.rel.copi <- subset(melted.fun.clean.rel, OTU %in% rownames(tab_ra.f.SOM.TN.copi))
melted.fun.clean.rel.copi$Trophic <- "Copio"


mer.fun.clean.rel.trophic <- rbind(melted.fun.clean.rel.oligo, melted.fun.clean.rel.copi)
mer.fun.clean.rel.trophic$Abundance



df.trophic.rel.p2 <- ggplot(mer.fun.clean.rel.trophic, aes(x=Field, y = Abundance/9, fill = Trophic)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("#669933", "#CCCC99")) +
  xlab('')+
  ylab("RA\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.trophic.rel.p2


##Bacteria
melted.bac.clean.rel <- psmelt(bac.clean.rel)


lb <- c('CC1',
        'CC2',
        'CJ1',
        'CJ2',
        'DG1',
        'DG2',
        'JJ1',
        'JJ2',
        'MY1',
        'MY2',
        'NJ1',
        'NJ2',
        'IS1',
        'IS2',
        'YS1',
        'YS2',
        'UF1',
        'UF2',
        'CC1.18',
        'CC2.18',
        'CJ1.18',
        'CJ2.18',
        'DG1.18',
        'DG2.18',
        'IS1.18',
        'IS2.18',
        'UF1.18',
        'UF2.18')

melted.bac.clean.rel$Field <- factor(melted.bac.clean.rel$Field, levels = lb)

melted.bac.clean.rel.oligo <- subset(melted.bac.clean.rel, OTU %in% rownames(tab_ra.b.SOM.TN.oligo))
melted.bac.clean.rel.oligo$Trophic <- "Oligo"

melted.bac.clean.rel.copi <- subset(melted.bac.clean.rel, OTU %in% rownames(tab_ra.b.SOM.TN.copi))
melted.bac.clean.rel.copi$Trophic <- "Copio"


mer.bac.clean.rel.trophic <- rbind(melted.bac.clean.rel.oligo, melted.bac.clean.rel.copi)

df.trophic.rel.p3 <- ggplot(mer.bac.clean.rel.trophic, aes(x=Field, y = Abundance/9, fill = Trophic)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("#669933", "#CCCC99")) +
  xlab('')+
  ylab("RA\n") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 1,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())

df.trophic.rel.p3


## relationship between diversity and abundance of copio and oligo OTUs
ps1.meta.trophic


lmMod <- lm(bac_RA_oligo ~ SOM, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.7286 , p-value: < 2.2e-16

lmMod <- lm(bac_RA_copio ~ SOM, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.6753 , p-value: < 2.2e-16

ggplot(ps2.meta.trophic, aes(x=Simpson, y=arch_RA_copio)) +
  xlab('\n Simpson')+
  ylab("archaeal RA copio \n") +
  geom_point(size=3, alpha=0.7) +
  #scale_colour_manual(labels = c('high','medium','low'), values = c("#CC9900", "#0066CC", '#336633'))+
  #scale_shape_manual(values=c(16,15,17))+
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #geom_hline(yintercept=quantile(df.norm.degree$y2018_all_deg,prob=1-1/100) , color="maroon4", linetype='dotted')+
  #geom_vline(xintercept=quantile(df.norm.degree$y2017_all_deg,prob=1-1/100), color="maroon4", linetype='dotted')+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred")


lmMod <- lm(bac_RA_oligo ~ Sand, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.3901  , p-value: < 2.2e-16

lmMod <- lm(bac_RA_copio ~ Sand, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.1939 , p-value: 1.35e-13

lmMod <- lm(bac_RA_oligo ~ pH, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.3962  , p-value: < 2.2e-16

lmMod <- lm(bac_RA_copio ~ pH, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.2822   , p-value: < 2.2e-16


lmMod <- lm(bac_RA_copio ~ Shannon, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.07368 , p-value: 7.376e-06

lmMod <- lm(bac_RA_oligo ~ Shannon, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2:0.01963, p-value: 0.01479

lmMod <- lm(bac_RA_copio ~ Simpson, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.003465 , p-value: 0.1724

lmMod <- lm(bac_RA_oligo ~ Simpson, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.1926 , p-value: 1.646e-13


lmMod <- lm(arch_RA_oligo ~ SOM, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.3814 , p-value: < 2.2e-16

lmMod <- lm(arch_RA_copio ~ SOM, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.4568 , p-value: < 2.2e-16

lmMod <- lm(arch_RA_oligo ~ Sand, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.2132  , p-value: 6.256e-15

lmMod <- lm(arch_RA_copio ~ Sand, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2:0.1644  , p-value: 1.306e-11

lmMod <- lm(arch_RA_copio ~ Shannon, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.0005734, p-value: 0.2858

lmMod <- lm(arch_RA_oligo ~ Shannon, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.01022, p-value: 0.05925

lmMod <- lm(arch_RA_copio ~ Simpson, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.01768 , p-value: 0.01961

lmMod <- lm(arch_RA_oligo ~ Simpson, data=ps1.meta.trophic)  # build the model
summary (lmMod) #R2: 0.118 , p-value: 1.303e-08




lmMod <- lm(fun_RA_oligo ~ SOM, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2:0.2031 , p-value:3.126e-14

lmMod <- lm(fun_RA_copio ~ SOM, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2:0.1766  , p-value: 2.017e-12

lmMod <- lm(fun_RA_oligo ~ Sand, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.2659 , p-value: < 2.2e-16

lmMod <- lm(fun_RA_copio ~ Sand, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.135  , p-value: 1.018e-09

lmMod <- lm(fun_RA_oligo ~ Clay, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.2651   , p-value: < 2.2e-16

lmMod <- lm(fun_RA_copio ~ Clay, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.1536    , p-value: .086e-09

lmMod <- lm(fun_RA_oligo ~ Silt, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.1345  , p-value: 1.163e-09

lmMod <- lm(fun_RA_copio ~ Silt, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.04306   , p-value: 0.0005383

lmMod <- lm(fun_RA_copio ~ Shannon, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2:-0.002851 , p-value: 0.593

lmMod <- lm(fun_RA_oligo ~ Shannon, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.06837, p-value: 1.558e-05

lmMod <- lm(fun_RA_copio ~ Simpson, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.05707  , p-value: 7.591e-05

lmMod <- lm(fun_RA_oligo ~ Simpson, data=ps3.meta.trophic)  # build the model
summary (lmMod) #R2: 0.01766  , p-value: 0.01967



###Comparison of oligo and copiotroph at the field level
## JJ1 UF2

ps3.meta.trophic.JJ1.UF2<- subset(ps3.meta.trophic, Field %in% c("JJ1", "UF2"))

x <- subset(ps3.meta.trophic.JJ1.UF2, Field=='JJ1')$fun_RA_oligo
y <- subset(ps3.meta.trophic.JJ1.UF2, Field=='UF2')$fun_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406


p <- ggplot(data =ps3.meta.trophic.JJ1.UF2, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps3.meta.trophic.JJ1.UF2, Field=='JJ1')$fun_RA_copio
y <- subset(ps3.meta.trophic.JJ1.UF2, Field=='UF2')$fun_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.01875771


p <- ggplot(data =ps3.meta.trophic.JJ1.UF2, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


x <- subset(ps3.meta.trophic.JJ1.UF2, Field=='JJ1')$Shannon
y <- subset(ps3.meta.trophic.JJ1.UF2, Field=='UF2')$Shannon
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.01875771


p <- ggplot(data =ps3.meta.trophic.JJ1.UF2, aes(x=Field, y=Shannon)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("Shannon \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


## JJ1 DG1

ps3.meta.trophic.JJ1.DG1<- subset(ps3.meta.trophic, Field %in% c("JJ1", "DG1"))

x <- subset(ps3.meta.trophic.JJ1.DG1, Field=='JJ1')$fun_RA_oligo
y <- subset(ps3.meta.trophic.JJ1.DG1, Field=='DG1')$fun_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.JJ1.DG1, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps3.meta.trophic.JJ1.DG1, Field=='JJ1')$fun_RA_copio
y <- subset(ps3.meta.trophic.JJ1.DG1, Field=='DG1')$fun_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #8.227067e-05


p <- ggplot(data =ps3.meta.trophic.JJ1.DG1, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



## IS1 NJ1

ps3.meta.trophic.IS1.NJ1<- subset(ps3.meta.trophic, Field %in% c("IS1", "NJ1"))

x <- subset(ps3.meta.trophic.IS1.NJ1, Field=='IS1')$fun_RA_oligo
y <- subset(ps3.meta.trophic.IS1.NJ1, Field=='NJ1')$fun_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.01061292


p <- ggplot(data =ps3.meta.trophic.IS1.NJ1, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps3.meta.trophic.IS1.NJ1, Field=='IS1')$fun_RA_copio
y <- subset(ps3.meta.trophic.IS1.NJ1, Field=='NJ1')$fun_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps3.meta.trophic.IS1.NJ1, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


x <- subset(ps3.meta.trophic.IS1.NJ1, Field=='IS1')$Shannon
y <- subset(ps3.meta.trophic.IS1.NJ1, Field=='NJ1')$Shannon
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps3.meta.trophic.IS1.NJ1, aes(x=Field, y=Shannon)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("Shannon \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

## IS2 YS2

ps3.meta.trophic.IS2.YS2<- subset(ps3.meta.trophic, Field %in% c("IS2", "YS2"))

x <- subset(ps3.meta.trophic.IS2.YS2, Field=='IS2')$fun_RA_oligo
y <- subset(ps3.meta.trophic.IS2.YS2, Field=='YS2')$fun_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.02443439


p <- ggplot(data =ps3.meta.trophic.IS2.YS2, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps3.meta.trophic.IS2.YS2, Field=='IS2')$fun_RA_copio
y <- subset(ps3.meta.trophic.IS2.YS2, Field=='YS2')$fun_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps3.meta.trophic.IS2.YS2, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


## CJ1 MY2

ps3.meta.trophic.CJ1.MY2<- subset(ps3.meta.trophic, Field %in% c("CJ1", "MY2"))

x <- subset(ps3.meta.trophic.CJ1.MY2, Field=='CJ1')$fun_RA_oligo
y <- subset(ps3.meta.trophic.CJ1.MY2, Field=='MY2')$fun_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CJ1.MY2, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps3.meta.trophic.CJ1.MY2, Field=='CJ1')$fun_RA_copio
y <- subset(ps3.meta.trophic.CJ1.MY2, Field=='MY2')$fun_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CJ1.MY2, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



## Bacteria
## JJ1 UF2

ps1.meta.trophic.JJ1.UF2<- subset(ps1.meta.trophic, Field %in% c("JJ1", "UF2"))

x <- subset(ps1.meta.trophic.JJ1.UF2, Field=='JJ1')$bac_RA_oligo
y <- subset(ps1.meta.trophic.JJ1.UF2, Field=='UF2')$bac_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps1.meta.trophic.JJ1.UF2, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.JJ1.UF2, Field=='JJ1')$bac_RA_copio
y <- subset(ps1.meta.trophic.JJ1.UF2, Field=='UF2')$bac_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0007815714


p <- ggplot(data =ps1.meta.trophic.JJ1.UF2, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



## JJ1 DG1

ps1.meta.trophic.JJ1.DG1<- subset(ps1.meta.trophic, Field %in% c("JJ1", "DG1"))

x <- subset(ps1.meta.trophic.JJ1.DG1, Field=='JJ1')$bac_RA_oligo
y <- subset(ps1.meta.trophic.JJ1.DG1, Field=='DG1')$bac_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.JJ1.DG1, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.JJ1.DG1, Field=='JJ1')$bac_RA_copio
y <- subset(ps1.meta.trophic.JJ1.DG1, Field=='DG1')$bac_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.JJ1.DG1, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



## IS1 NJ1

ps1.meta.trophic.IS1.NJ1<- subset(ps1.meta.trophic, Field %in% c("IS1", "NJ1"))

x <- subset(ps1.meta.trophic.IS1.NJ1, Field=='IS1')$bac_RA_oligo
y <- subset(ps1.meta.trophic.IS1.NJ1, Field=='NJ1')$bac_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.IS1.NJ1, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.IS1.NJ1, Field=='IS1')$bac_RA_copio
y <- subset(ps1.meta.trophic.IS1.NJ1, Field=='NJ1')$bac_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps1.meta.trophic.IS1.NJ1, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



x <- subset(ps1.meta.trophic.IS1.NJ1, Field=='IS1')$Shannon
y <- subset(ps1.meta.trophic.IS1.NJ1, Field=='NJ1')$Shannon
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps1.meta.trophic.IS1.NJ1, aes(x=Field, y=Shannon)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

## IS2 YS2

ps1.meta.trophic.IS2.YS2<- subset(ps1.meta.trophic, Field %in% c("IS2", "YS2"))

x <- subset(ps1.meta.trophic.IS2.YS2, Field=='IS2')$bac_RA_oligo
y <- subset(ps1.meta.trophic.IS2.YS2, Field=='YS2')$bac_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.02443439


p <- ggplot(data =ps1.meta.trophic.IS2.YS2, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.IS2.YS2, Field=='IS2')$bac_RA_copio
y <- subset(ps1.meta.trophic.IS2.YS2, Field=='YS2')$bac_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps1.meta.trophic.IS2.YS2, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


## CJ1 MY2

ps1.meta.trophic.CJ1.MY2<- subset(ps1.meta.trophic, Field %in% c("CJ1", "MY2"))

x <- subset(ps1.meta.trophic.CJ1.MY2, Field=='CJ1')$bac_RA_oligo
y <- subset(ps1.meta.trophic.CJ1.MY2, Field=='MY2')$bac_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CJ1.MY2, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CJ1.MY2, Field=='CJ1')$bac_RA_copio
y <- subset(ps1.meta.trophic.CJ1.MY2, Field=='MY2')$bac_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CJ1.MY2, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


## Archaea
## JJ1 UF2

ps1.meta.trophic.JJ1.UF2<- subset(ps1.meta.trophic, Field %in% c("JJ1", "UF2"))

x <- subset(ps1.meta.trophic.JJ1.UF2, Field=='JJ1')$arch_RA_oligo
y <- subset(ps1.meta.trophic.JJ1.UF2, Field=='UF2')$arch_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps1.meta.trophic.JJ1.UF2, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.JJ1.UF2, Field=='JJ1')$arch_RA_copio
y <- subset(ps1.meta.trophic.JJ1.UF2, Field=='UF2')$arch_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0007815714


p <- ggplot(data =ps1.meta.trophic.JJ1.UF2, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



## JJ1 DG1

ps1.meta.trophic.JJ1.DG1<- subset(ps1.meta.trophic, Field %in% c("JJ1", "DG1"))

x <- subset(ps1.meta.trophic.JJ1.DG1, Field=='JJ1')$arch_RA_oligo
y <- subset(ps1.meta.trophic.JJ1.DG1, Field=='DG1')$arch_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.JJ1.DG1, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.JJ1.DG1, Field=='JJ1')$arch_RA_copio
y <- subset(ps1.meta.trophic.JJ1.DG1, Field=='DG1')$arch_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.JJ1.DG1, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



## IS1 NJ1

ps1.meta.trophic.IS1.NJ1<- subset(ps1.meta.trophic, Field %in% c("IS1", "NJ1"))

x <- subset(ps1.meta.trophic.IS1.NJ1, Field=='IS1')$arch_RA_oligo
y <- subset(ps1.meta.trophic.IS1.NJ1, Field=='NJ1')$arch_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.IS1.NJ1, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.IS1.NJ1, Field=='IS1')$arch_RA_copio
y <- subset(ps1.meta.trophic.IS1.NJ1, Field=='NJ1')$arch_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps1.meta.trophic.IS1.NJ1, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


## IS2 YS2

ps1.meta.trophic.IS2.YS2<- subset(ps1.meta.trophic, Field %in% c("IS2", "YS2"))

x <- subset(ps1.meta.trophic.IS2.YS2, Field=='IS2')$arch_RA_oligo
y <- subset(ps1.meta.trophic.IS2.YS2, Field=='YS2')$arch_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.02443439


p <- ggplot(data =ps1.meta.trophic.IS2.YS2, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.IS2.YS2, Field=='IS2')$arch_RA_copio
y <- subset(ps1.meta.trophic.IS2.YS2, Field=='YS2')$arch_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.0002879473


p <- ggplot(data =ps1.meta.trophic.IS2.YS2, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


## CJ1 MY2

ps1.meta.trophic.CJ1.MY2<- subset(ps1.meta.trophic, Field %in% c("CJ1", "MY2"))

x <- subset(ps1.meta.trophic.CJ1.MY2, Field=='CJ1')$arch_RA_oligo
y <- subset(ps1.meta.trophic.CJ1.MY2, Field=='MY2')$arch_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CJ1.MY2, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CJ1.MY2, Field=='CJ1')$arch_RA_copio
y <- subset(ps1.meta.trophic.CJ1.MY2, Field=='MY2')$arch_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CJ1.MY2, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



## Similar texture different nutrients
## CC1 UF1

ps1.meta.trophic.CC1.UF1<- subset(ps1.meta.trophic, Field %in% c("CC1", "UF1"))

x <- subset(ps1.meta.trophic.CC1.UF1, Field=='CC1')$arch_RA_oligo
y <- subset(ps1.meta.trophic.CC1.UF1, Field=='UF1')$arch_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC1.UF1, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CC1.UF1, Field=='CC1')$arch_RA_copio
y <- subset(ps1.meta.trophic.CC1.UF1, Field=='UF1')$arch_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC1.UF1, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p




#Fungi
ps3.meta.trophic.CC1.UF1<- subset(ps3.meta.trophic, Field %in% c("CC1", "UF1"))

x <- subset(ps3.meta.trophic.CC1.UF1, Field=='CC1')$fun_RA_oligo
y <- subset(ps3.meta.trophic.CC1.UF1, Field=='UF1')$fun_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CC1.UF1, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps3.meta.trophic.CC1.UF1, Field=='CC1')$fun_RA_copio
y <- subset(ps3.meta.trophic.CC1.UF1, Field=='UF1')$fun_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CC1.UF1, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


#Bacteria
ps1.meta.trophic.CC1.UF1<- subset(ps1.meta.trophic, Field %in% c("CC1", "UF1"))

x <- subset(ps1.meta.trophic.CC1.UF1, Field=='CC1')$bac_RA_oligo
y <- subset(ps1.meta.trophic.CC1.UF1, Field=='UF1')$bac_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC1.UF1, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CC1.UF1, Field=='CC1')$bac_RA_copio
y <- subset(ps1.meta.trophic.CC1.UF1, Field=='UF1')$bac_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC1.UF1, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CC1.UF1, Field=='CC1')$Simpson
y <- subset(ps1.meta.trophic.CC1.UF1, Field=='UF1')$Simpson
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05

p <- ggplot(data =ps1.meta.trophic.CC1.UF1, aes(x=Field, y=Simpson)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

## CC2 JJ2

ps1.meta.trophic.CC2.JJ2<- subset(ps1.meta.trophic, Field %in% c("CC2", "JJ2"))

x <- subset(ps1.meta.trophic.CC2.JJ2, Field=='CC2')$arch_RA_oligo
y <- subset(ps1.meta.trophic.CC2.JJ2, Field=='JJ2')$arch_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.JJ2, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CC2.JJ2, Field=='CC2')$arch_RA_copio
y <- subset(ps1.meta.trophic.CC2.JJ2, Field=='JJ2')$arch_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.JJ2, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

ps1.meta.trophic.CC2.JJ2<- subset(ps1.meta.trophic, Field %in% c("CC2", "JJ2"))

x <- subset(ps1.meta.trophic.CC2.JJ2, Field=='CC2')$Simpson
y <- subset(ps1.meta.trophic.CC2.JJ2, Field=='JJ2')$Simpson
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.JJ2, aes(x=Field, y=Simpson)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


#Fungi
ps3.meta.trophic.CC2.JJ2<- subset(ps3.meta.trophic, Field %in% c("CC2", "JJ2"))

x <- subset(ps3.meta.trophic.CC2.JJ2, Field=='CC2')$fun_RA_oligo
y <- subset(ps3.meta.trophic.CC2.JJ2, Field=='JJ2')$fun_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CC2.JJ2, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps3.meta.trophic.CC2.JJ2, Field=='CC2')$fun_RA_copio
y <- subset(ps3.meta.trophic.CC2.JJ2, Field=='JJ2')$fun_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CC2.JJ2, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


#Bacteria
ps1.meta.trophic.CC2.JJ2<- subset(ps1.meta.trophic, Field %in% c("CC2", "JJ2"))

x <- subset(ps1.meta.trophic.CC2.JJ2, Field=='CC2')$bac_RA_oligo
y <- subset(ps1.meta.trophic.CC2.JJ2, Field=='JJ2')$bac_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.JJ2, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CC2.JJ2, Field=='CC2')$bac_RA_copio
y <- subset(ps1.meta.trophic.CC2.JJ2, Field=='JJ2')$bac_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.JJ2, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


## CC2.18 DG2.18

ps1.meta.trophic.CC2.18.DG2.18<- subset(ps1.meta.trophic, Field %in% c("CC2.18", "DG2.18"))

x <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='CC2.18')$arch_RA_oligo
y <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='DG2.18')$arch_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.18.DG2.18, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='CC2.18')$arch_RA_copio
y <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='DG2.18')$arch_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.18.DG2.18, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p




#Fungi
ps3.meta.trophic.CC2.18.DG2.18<- subset(ps3.meta.trophic, Field %in% c("CC2.18", "DG2.18"))

x <- subset(ps3.meta.trophic.CC2.18.DG2.18, Field=='CC2.18')$fun_RA_oligo
y <- subset(ps3.meta.trophic.CC2.18.DG2.18, Field=='DG2.18')$fun_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CC2.18.DG2.18, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps3.meta.trophic.CC2.18.DG2.18, Field=='CC2.18')$fun_RA_copio
y <- subset(ps3.meta.trophic.CC2.18.DG2.18, Field=='DG2.18')$fun_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CC2.18.DG2.18, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")

x <- subset(ps3.meta.trophic.CC2.18.DG2.18, Field=='CC2.18')$Shannon
y <- subset(ps3.meta.trophic.CC2.18.DG2.18, Field=='DG2.18')$Shannon
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps3.meta.trophic.CC2.18.DG2.18, aes(x=Field, y=Shannon)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("Shannon \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


x <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='CC2.18')$Shannon
y <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='DG2.18')$Shannon
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05

p <- ggplot(data =ps1.meta.trophic.CC2.18.DG2.18, aes(x=Field, y=Shannon)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("Shannon \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p

#Bacteria
ps1.meta.trophic.CC2.18.DG2.18<- subset(ps1.meta.trophic, Field %in% c("CC2.18", "DG2.18"))

x <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='CC2.18')$bac_RA_oligo
y <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='DG2.18')$bac_RA_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.18.DG2.18, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

x <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='CC2.18')$bac_RA_copio
y <- subset(ps1.meta.trophic.CC2.18.DG2.18, Field=='DG2.18')$bac_RA_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #4.113534e-05


p <- ggplot(data =ps1.meta.trophic.CC2.18.DG2.18, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p




## 
otu.list <- rbind(phy.list, arch.list, fun.list)

cor.bac.OTU.trim.sig.SOM.TN <- subset(cor.bac.OTU.trim.sig, row %in% c("SOM", "TN"))
names(cor.bac.OTU.trim.sig.SOM)[2] <- "OTU"
cor.bac.OTU.trim.sig.SOM.otu.id <- merge(cor.bac.OTU.trim.sig.SOM, otu.list, by = "OTU")
head(cor.bac.OTU.trim.sig.SOM.otu.id)
names(cor.bac.OTU.trim.sig.TN)[2] <- "OTU"
cor.bac.OTU.trim.sig.TN.otu.id <- merge(cor.bac.OTU.trim.sig.TN, otu.list, by = "OTU")
head(cor.bac.OTU.trim.sig.TN.otu.id)


names(cor.fun.OTU.trim.sig.SOM)[2] <- "OTU"
cor.fun.OTU.trim.sig.SOM.otu.id <- merge(cor.fun.OTU.trim.sig.SOM, otu.list, by = "OTU")
head(cor.fun.OTU.trim.sig.SOM.otu.id)
names(cor.fun.OTU.trim.sig.TN)[2] <- "OTU"
cor.fun.OTU.trim.sig.TN.otu.id <- merge(cor.fun.OTU.trim.sig.TN, otu.list, by = "OTU")
head(cor.fun.OTU.trim.sig.TN.otu.id)


names(cor.arch.OTU.trim.sig.SOM)[2] <- "OTU"
cor.arch.OTU.trim.sig.SOM.otu.id <- merge(cor.arch.OTU.trim.sig.SOM, otu.list, by = "OTU")
head(cor.arch.OTU.trim.sig.SOM.otu.id)
names(cor.arch.OTU.trim.sig.TN)[2] <- "OTU"
cor.arch.OTU.trim.sig.TN.otu.id <- merge(cor.arch.OTU.trim.sig.TN, otu.list, by = "OTU")
head(cor.arch.OTU.trim.sig.TN.otu.id)

cor.bac.OTU.trim.sig.SOM.otu.id$Category <- ifelse(cor.bac.OTU.trim.sig.SOM.otu.id$cor > 0, "Copio", "Oligo")
cor.bac.OTU.trim.sig.TN.otu.id$Category <- ifelse(cor.bac.OTU.trim.sig.TN.otu.id$cor > 0, "Copio", "Oligo")

cor.arch.OTU.trim.sig.SOM.otu.id$Category <- ifelse(cor.arch.OTU.trim.sig.SOM.otu.id$cor > 0, "Copio", "Oligo")
cor.arch.OTU.trim.sig.TN.otu.id$Category <- ifelse(cor.arch.OTU.trim.sig.TN.otu.id$cor > 0, "Copio", "Oligo")

cor.fun.OTU.trim.sig.SOM.otu.id$Category <- ifelse(cor.fun.OTU.trim.sig.SOM.otu.id$cor > 0, "Copio", "Oligo")
cor.fun.OTU.trim.sig.TN.otu.id$Category <- ifelse(cor.fun.OTU.trim.sig.TN.otu.id$cor > 0, "Copio", "Oligo")


write.xlsx(cor.bac.OTU.trim.sig.SOM.otu.id, "List of bacterial putative oligo and copiotrophs_SOM.xlsx")
write.xlsx(cor.arch.OTU.trim.sig.SOM.otu.id, "List of archaeal putative oligo and copiotrophs_SOM.xlsx")
write.xlsx(cor.fun.OTU.trim.sig.SOM.otu.id, "List of fungal putative oligo and copiotrophs_SOM.xlsx")


write.xlsx(cor.bac.OTU.trim.sig.TN.otu.id, "List of bacterial putative oligo and copiotrophs_TN.xlsx")
write.xlsx(cor.arch.OTU.trim.sig.TN.otu.id, "List of archaeal putative oligo and copiotrophs_TN.xlsx")
write.xlsx(cor.fun.OTU.trim.sig.TN.otu.id, "List of fungal putative oligo and copiotrophs_TN.xlsx")


## Taxonomic composition of oligo and copiotrophs
list.putative_oligo.bac <- subset(otu.list, OTU %in% rownames(tab_ra.b.SOM.TN.oligo))
list.putative_oligo.arch <- subset(otu.list, OTU %in% rownames(tab_ra.a.SOM.TN.oligo))
list.putative_oligo.fun <- subset(otu.list, OTU %in% rownames(tab_ra.f.SOM.TN.oligo))


list.putative_copio.bac <- subset(otu.list, OTU %in% rownames(tab_ra.b.SOM.TN.copi))
list.putative_copio.arch <- subset(otu.list, OTU %in% rownames(tab_ra.a.SOM.TN.copi))
list.putative_copio.fun <- subset(otu.list, OTU %in% rownames(tab_ra.f.SOM.TN.copi))


### counts of bacteria OTUs
otu.tax.bac <- tax_table(bac.clean.ss)

length(rownames(otu.tax.bac)) #22623
otu.tax.bac <- data.frame(otu.tax.bac, stringsAsFactors = F)
otu.tax.bac$Phylum2 <-otu.tax.bac$Phylum

rownames(otu.tax.bac)[is.na(rownames(otu.tax.bac))]
otu.tax.bac$Phylum2[which(otu.tax.bac$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
otu.tax.bac$Phylum2[which(otu.tax.bac$Class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
otu.tax.bac$Phylum2[which(otu.tax.bac$Class == "Deltaproteobacteria")] <- "Deltaproteobacteria"




bacteria_oligo <- as.data.frame(table(otu.tax.bac[list.putative_oligo.bac$OTU, "Phylum2"] ) )
colnames(bacteria_oligo) <- c("Class", "Oligo")
bacteria_copio <- as.data.frame(table(otu.tax.bac[list.putative_copio.bac$OTU, "Phylum2"] ) )
colnames(bacteria_copio) <- c("Class", "Copio")

bacteria_copi_oli <- merge(bacteria_oligo, bacteria_copio, all=T, by="Class") 
 
bacteria_copi_oli

bacteria_all_OTUs <- as.data.frame(table(otu.tax.bac[, "Phylum2"] ) )
colnames(bacteria_all_OTUs) <- c("Class", "all bOTUs")
bacteria_modules <- merge(bacteria_copi_oli, bacteria_all_OTUs, all=T, by="Class") 
bacteria_modules

bacteria_modules_mat <- bacteria_modules[2:4]
rownames(bacteria_modules_mat) <- bacteria_modules$Class
bacteria_modules_mat[is.na(bacteria_modules_mat)] <- 0
colSums(bacteria_modules_mat)

bacteria_modules_prop <- t(t(bacteria_modules_mat)/colSums(bacteria_modules_mat) ) * 1
bacteria_modules_prop
colSums(bacteria_modules_prop)


bp <- barplot(cbind(bacteria_modules_prop[,1:2], NA, bacteria_modules_prop[,3]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_16s[rownames(bacteria_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(bacteria_modules_mat)[1:4], NA,
                       colSums(bacteria_modules_mat)[5]), xpd=T, cex=.6, las=2)
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )


write.xlsx(bacteria_modules_mat, "Bacteria_oligo and copio_count.xlsx")

### counts of archaea OTUs
otu.tax.arch <- tax_table(arch.clean.ss)

length(rownames(otu.tax.arch)) #1139
otu.tax.arch <- data.frame(otu.tax.arch, stringsAsFactors = F)
otu.tax.arch$Phylum2 <-otu.tax.arch$Phylum

archaea_oligo <- as.data.frame(table(otu.tax.arch[list.putative_oligo.arch$OTU, "Phylum2"] ) )
colnames(archaea_oligo) <- c("Class", "Oligo")
archaea_copio <- as.data.frame(table(otu.tax.arch[list.putative_copio.arch$OTU, "Phylum2"] ) )
colnames(archaea_copio) <- c("Class", "Copio")

archaea_copi_oli <- merge(archaea_oligo, archaea_copio, all=T, by="Class") 

archaea_copi_oli

archaea_all_OTUs <- as.data.frame(table(otu.tax.arch[, "Phylum2"] ) )
colnames(archaea_all_OTUs) <- c("Class", "all aOTUs")
archaea_modules <- merge(archaea_copi_oli, archaea_all_OTUs, all=T, by="Class") 
archaea_modules

archaea_modules_mat <- archaea_modules[2:4]
rownames(archaea_modules_mat) <- archaea_modules$Class
archaea_modules_mat[is.na(archaea_modules_mat)] <- 0
colSums(archaea_modules_mat)

archaea_modules_prop <- t(t(archaea_modules_mat)/colSums(archaea_modules_mat) ) * 1
archaea_modules_prop
colSums(archaea_modules_prop)


ap <- barplot(cbind(archaea_modules_prop[,1:2], NA, archaea_modules_prop[,3]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_arch[rownames(archaea_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(archaea_modules_mat)[1:4], NA,
                       colSums(archaea_modules_mat)[5]), xpd=T, cex=.6, las=2)

plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_arch_legend$labels), 
       fill=rev(PHYLA_label_cols_arch_legend$cols), 
       border=rev(PHYLA_label_cols_arch_legend$cols) )


write.xlsx(archaea_modules_mat, "Archaea_oligo and copio_count.xlsx")
### counts of fungi OTUs
otu.tax.fun <- tax_table(fun.clean.ss)

length(rownames(otu.tax.fun)) #1139
otu.tax.fun <- data.frame(otu.tax.fun, stringsAsFactors = F)
otu.tax.fun$Phylum2 <-otu.tax.fun$Phylum

fungi_oligo <- as.data.frame(table(otu.tax.fun[list.putative_oligo.fun$OTU, "Phylum2"] ) )
colnames(fungi_oligo) <- c("Class", "Oligo")
fungi_copio <- as.data.frame(table(otu.tax.fun[list.putative_copio.fun$OTU, "Phylum2"] ) )
colnames(fungi_copio) <- c("Class", "Copio")

fungi_copi_oli <- merge(fungi_oligo, fungi_copio, all=T, by="Class") 

fungi_copi_oli

fungi_all_OTUs <- as.data.frame(table(otu.tax.fun[, "Phylum2"] ) )
colnames(fungi_all_OTUs) <- c("Class", "all aOTUs")
fungi_modules <- merge(fungi_copi_oli, fungi_all_OTUs, all=T, by="Class") 
fungi_modules

fungi_modules_mat <- fungi_modules[2:4]
rownames(fungi_modules_mat) <- fungi_modules$Class
fungi_modules_mat[is.na(fungi_modules_mat)] <- 0
colSums(fungi_modules_mat)

fungi_modules_prop <- t(t(fungi_modules_mat)/colSums(fungi_modules_mat) ) * 1
fungi_modules_prop
colSums(fungi_modules_prop)

PHYLA_label_cols_its$labels[PHYLA_label_cols_its$labels == "Unassigned"] <- "unidentified"
PHYLA_label_cols_its$cols[PHYLA_label_cols_its$labels == "Glomeromycota"] <- "black"
rownames(PHYLA_label_cols_its) <- PHYLA_label_cols_its$labels

fp <- barplot(cbind(fungi_modules_prop[,1:2], NA, fungi_modules_prop[,3]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_its[rownames(fungi_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(fungi_modules_mat)[1:4], NA,
                       colSums(fungi_modules_mat)[5]), xpd=T, cex=.6, las=2)
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=1, 
       legend=rev(PHYLA_label_cols_its_legend$labels),
       fill=rev(PHYLA_label_cols_its_legend$cols),
       border=rev(PHYLA_label_cols_its_legend$cols) )

write.xlsx(fungi_modules_mat, "Fungi_oligo and copio_count.xlsx")



### list of hub

df.CJ1_all_deg.t <- df.CJ1_all_deg[1]
df.CJ1_all_betweenness.t <-df.CJ1_all_betweenness[1]

df.CJ2_all_deg.t <- df.CJ2_all_deg[1]
df.CJ2_all_betweenness.t <-df.CJ2_all_betweenness[1]


df.CJ1_18_all_deg.t <- df.CJ1_18_all_deg[1]
df.CJ1_18_all_betweenness.t <-df.CJ1_18_all_betweenness[1]

df.CJ2_18_all_deg.t <- df.CJ2_18_all_deg[1]
df.CJ2_18_all_betweenness.t <-df.CJ2_18_all_betweenness[1]

df.MY1_all_deg.t <- df.MY1_all_deg[1]
df.MY1_all_betweenness.t <-df.MY1_all_betweenness[1]

df.MY2_all_deg.t <- df.MY2_all_deg[1]
df.MY2_all_betweenness.t <-df.MY2_all_betweenness[1]


df.NJ1_all_deg.t <- df.NJ1_all_deg[1]
df.NJ1_all_betweenness.t <-df.NJ1_all_betweenness[1]

df.NJ2_all_deg.t <- df.NJ2_all_deg[1]
df.NJ2_all_betweenness.t <-df.NJ2_all_betweenness[1]

df.YS1_all_deg.t <- df.YS1_all_deg[1]
df.YS1_all_betweenness.t <-df.YS1_all_betweenness[1]

df.YS2_all_deg.t <- df.YS2_all_deg[1]
df.YS2_all_betweenness.t <-df.YS2_all_betweenness[1]

CJ1.degree.betw<-cbind(df.CJ1_all_deg.t, df.CJ1_all_betweenness.t)
CJ1.degree.betw$OTU_id <- rownames(CJ1.degree.betw)
CJ1.degree.betw.id <- merge(CJ1.degree.betw, otu.list, by = "OTU_id")

CJ2.degree.betw<-cbind(df.CJ2_all_deg.t, df.CJ2_all_betweenness.t)
CJ2.degree.betw$OTU_id <- rownames(CJ2.degree.betw)
CJ2.degree.betw.id <- merge(CJ2.degree.betw, otu.list, by = "OTU_id")

CJ1_18.degree.betw<-cbind(df.CJ1_18_all_deg.t, df.CJ1_18_all_betweenness.t)
CJ1_18.degree.betw$OTU_id <- rownames(CJ1_18.degree.betw)
CJ1_18.degree.betw.id <- merge(CJ1_18.degree.betw, otu.list, by = "OTU_id")

CJ2_18.degree.betw<-cbind(df.CJ2_18_all_deg.t, df.CJ2_18_all_betweenness.t)
CJ2_18.degree.betw$OTU_id <- rownames(CJ2_18.degree.betw)
CJ2_18.degree.betw.id <- merge(CJ2_18.degree.betw, otu.list, by = "OTU_id")

YS1.degree.betw<-cbind(df.YS1_all_deg.t, df.YS1_all_betweenness.t)
YS1.degree.betw$OTU_id <- rownames(YS1.degree.betw)
YS1.degree.betw.id <- merge(YS1.degree.betw, otu.list, by = "OTU_id")

YS2.degree.betw<-cbind(df.YS2_all_deg.t, df.YS2_all_betweenness.t)
YS2.degree.betw$OTU_id <- rownames(YS2.degree.betw)
YS2.degree.betw.id <- merge(YS2.degree.betw, otu.list, by = "OTU_id")

MY1.degree.betw<-cbind(df.MY1_all_deg.t, df.MY1_all_betweenness.t)
MY1.degree.betw$OTU_id <- rownames(MY1.degree.betw)
MY1.degree.betw.id <- merge(MY1.degree.betw, otu.list, by = "OTU_id")

MY2.degree.betw<-cbind(df.MY2_all_deg.t, df.MY2_all_betweenness.t)
MY2.degree.betw$OTU_id <- rownames(MY2.degree.betw)
MY2.degree.betw.id <- merge(MY2.degree.betw, otu.list, by = "OTU_id")

NJ1.degree.betw<-cbind(df.NJ1_all_deg.t, df.NJ1_all_betweenness.t)
NJ1.degree.betw$OTU_id <- rownames(NJ1.degree.betw)
NJ1.degree.betw.id <- merge(NJ1.degree.betw, otu.list, by = "OTU_id")

NJ2.degree.betw<-cbind(df.NJ2_all_deg.t, df.NJ2_all_betweenness.t)
NJ2.degree.betw$OTU_id <- rownames(NJ2.degree.betw)
NJ2.degree.betw.id <- merge(NJ2.degree.betw, otu.list, by = "OTU_id")


write.xlsx(CJ1.degree.betw.id, "Network properties of CJ1 network and hub.xlsx")
write.xlsx(CJ2.degree.betw.id, "Network properties of CJ2 network and hub.xlsx")

write.xlsx(CJ1_18.degree.betw.id, "Network properties of CJ1_18 network and hub.xlsx")
write.xlsx(CJ2_18.degree.betw.id, "Network properties of CJ2_18 network and hub.xlsx")

write.xlsx(MY1.degree.betw.id, "Network properties of MY1 network and hub.xlsx")
write.xlsx(MY2.degree.betw.id, "Network properties of MY2 network and hub.xlsx")

write.xlsx(NJ1.degree.betw.id, "Network properties of NJ1 network and hub.xlsx")
write.xlsx(NJ2.degree.betw.id, "Network properties of NJ2 network and hub.xlsx")

write.xlsx(YS1.degree.betw.id, "Network properties of YS1 network and hub.xlsx")
write.xlsx(YS2.degree.betw.id, "Network properties of YS2 network and hub.xlsx")
## Hub and oligo, copio
hub.spar.CJ1
intersect(hub.spar.CJ1, list.putative_oligo.bac$OTU_id) #"B213_f_SC-I-84"     "B214_Sulfurifustis" "B1723_o_NA"         "B212_f_SC-I-84"     "B23_uncultured"  
intersect(hub.spar.CJ1, list.putative_oligo.arch$OTU_id)
intersect(hub.spar.CJ1, list.putative_oligo.fun$OTU_id) #"F70_Mortierella"   "F1253_o_NA"        "F645_unidentified"

intersect(hub.spar.CJ1, list.putative_copio.bac$OTU_id)
intersect(hub.spar.CJ1, list.putative_copio.arch$OTU_id)
intersect(hub.spar.CJ1, list.putative_copio.fun$OTU_id) #"F88_unidentified"


intersect(hub.spar.CJ2, list.putative_oligo.bac$OTU_id) #"B95_o_RCP2-54" "B17_uncultured Thermaerobacter sp." "B88_o_NA" 
intersect(hub.spar.CJ2, list.putative_oligo.arch$OTU_id)
intersect(hub.spar.CJ2, list.putative_oligo.fun$OTU_id)

intersect(hub.spar.CJ2, list.putative_copio.bac$OTU_id) #"B116_Reyranella"
intersect(hub.spar.CJ2, list.putative_copio.arch$OTU_id)
intersect(hub.spar.CJ2, list.putative_copio.fun$OTU_id) #"F25_o_Pleosporales" "F62_Nigrospora" 


intersect(hub.spar.YS1, list.putative_oligo.bac$OTU_id) #""B497_uncultured" "B212_f_SC-I-84"
intersect(hub.spar.YS1, list.putative_oligo.arch$OTU_id)
intersect(hub.spar.YS1, list.putative_oligo.fun$OTU_id) #"F277_o_Sordariales" "F29_Tetracladium"   "F707_o_NA"  "F902_Coniochaeta"   "F333_unidentified" 

intersect(hub.spar.YS1, list.putative_copio.bac$OTU_id) #"B103_Porphyrobacter" "B5_Sideroxydans"     "B300_mle1-7"
intersect(hub.spar.YS1, list.putative_copio.arch$OTU_id) #"A3_o_NA"
intersect(hub.spar.YS1, list.putative_copio.fun$OTU_id) #"F291_Byssonectria"


intersect(hub.spar.YS2, list.putative_oligo.bac$OTU_id) #"B163_uncultured Bellilinea sp."
intersect(hub.spar.YS2, list.putative_oligo.arch$OTU_id)
intersect(hub.spar.YS2, list.putative_oligo.fun$OTU_id) #"F3_Solicoccozyma"

intersect(hub.spar.YS2, list.putative_copio.bac$OTU_id) #"B40_f_Nocardioidaceae"   "B45_f_Blastocatellaceae" "B277_Desulfatiglans" 
intersect(hub.spar.YS2, list.putative_copio.arch$OTU_id)
intersect(hub.spar.YS2, list.putative_copio.fun$OTU_id) #"F197_o_NA"


intersect(hub.spar.MY1, list.putative_oligo.bac$OTU_id) #
intersect(hub.spar.MY1, list.putative_oligo.arch$OTU_id)
intersect(hub.spar.MY1, list.putative_oligo.fun$OTU_id) #

intersect(hub.spar.MY1, list.putative_copio.bac$OTU_id) #"B36_o_Sva0485"    "B5_Sideroxydans"  "B15_Thiobacillus"
intersect(hub.spar.MY1, list.putative_copio.arch$OTU_id) #
intersect(hub.spar.MY1, list.putative_copio.fun$OTU_id) #"F48_unidentified"  "F88_unidentified"  "F138_Rhizophydium"


intersect(hub.spar.MY2, list.putative_oligo.bac$OTU_id) #"B143_Rhodanobacter" "B12_GOUTA6"
intersect(hub.spar.MY2, list.putative_oligo.arch$OTU_id)
intersect(hub.spar.MY2, list.putative_oligo.fun$OTU_id) #

intersect(hub.spar.MY2, list.putative_copio.bac$OTU_id) #""B25_f_SC-I-84" "B106_Oryzihumus"  "B105_Defluviicoccus"  "B119_uncultured organism""B15_Thiobacillus"  
intersect(hub.spar.MY2, list.putative_copio.arch$OTU_id)
intersect(hub.spar.MY2, list.putative_copio.fun$OTU_id) #


intersect(hub.spar.NJ1, list.putative_oligo.bac$OTU_id) # "B22_uncultured bacterium" "B458_uncultured Aminicenantes bacterium" "B28_Desulfobacca"     
intersect(hub.spar.NJ1, list.putative_oligo.arch$OTU_id)
intersect(hub.spar.NJ1, list.putative_oligo.fun$OTU_id) #"F54_unidentified"

intersect(hub.spar.NJ1, list.putative_copio.bac$OTU_id) #"B866_o_KI89A clade"            "B19_uncultured soil bacterium" "B154_uncultured" 
intersect(hub.spar.NJ1, list.putative_copio.arch$OTU_id) #
intersect(hub.spar.NJ1, list.putative_copio.fun$OTU_id) #"F27_Neonectria" "F560_o_NA" "F92_Paraphaeosphaeria" "F248_unidentified"  


intersect(hub.spar.NJ2, list.putative_oligo.bac$OTU_id) #
intersect(hub.spar.NJ2, list.putative_oligo.arch$OTU_id)
intersect(hub.spar.NJ2, list.putative_oligo.fun$OTU_id) #

intersect(hub.spar.NJ2, list.putative_copio.bac$OTU_id) #""B25_f_SC-I-84" 
intersect(hub.spar.NJ2, list.putative_copio.arch$OTU_id)
intersect(hub.spar.NJ2, list.putative_copio.fun$OTU_id) #



hub.spar.CJ1_18

intersect(hub.spar.CJ1_18, list.putative_oligo.bac$OTU_id) # "B203_uncultured Sorangiineae bacterium"
intersect(hub.spar.CJ1_18, list.putative_oligo.arch$OTU_id)#"A20_Rice Cluster I"
intersect(hub.spar.CJ1_18, list.putative_oligo.fun$OTU_id) #"F184_o_Sordariales" "F335_unidentified" 

intersect(hub.spar.CJ1_18, list.putative_copio.bac$OTU_id) #"B106_Oryzihumus" "B38_uncultured Sphingobacteriales bacterium"
intersect(hub.spar.CJ1_18, list.putative_copio.arch$OTU_id)
intersect(hub.spar.CJ1_18, list.putative_copio.fun$OTU_id) #

hub.spar.CJ2_18

intersect(hub.spar.CJ2_18, list.putative_oligo.bac$OTU_id) #
intersect(hub.spar.CJ2_18, list.putative_oligo.arch$OTU_id)#
intersect(hub.spar.CJ2_18, list.putative_oligo.fun$OTU_id) #

intersect(hub.spar.CJ2_18, list.putative_copio.bac$OTU_id) #"B19_uncultured soil bacterium" "B15_Thiobacillus"
intersect(hub.spar.CJ2_18, list.putative_copio.arch$OTU_id)
intersect(hub.spar.CJ2_18, list.putative_copio.fun$OTU_id) #"F386_Scutellinia"

## Absolute count of oligo and copiotrophs in each field
tab_ra.b.SOM.TN.oligo

tab_ra.a.SOM.TN.oligo

tab_ra.f.SOM.TN.oligo


meta(bac.clean.ss.field)

bac.clean.ss.field <- merge_samples(bac.clean.ss, "Field")
meta(bac.clean.ss.field)

number_OTU_oligo_copio<- function(phylo, number.field, tab_oligo_copio){
  
phylo.field <- merge_samples(phylo, "Field")

phylo.field <- subset_samples(phylo.field, Field == number.field)
phylo.field <- phyloseq::filter_taxa(phylo.field, function(x) sum(x) != 0, TRUE)

otu.phylo.field <- otu_table(phylo.field)
otu.phylo.field <- data.frame(t(otu.phylo.field))
otu.phylo.field.oligo.copio <- subset(otu.phylo.field , rownames(otu.phylo.field) %in% rownames(tab_oligo_copio))

return(otu.phylo.field.oligo.copio)}

#Bacteria
## CC1
number.field <- "1"
CC1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(CC1.number.otu.oligo)) #301

CC1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(CC1.number.otu.copio)) #183

## CC1.18
number.field <- "2"
CC1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(CC1.18.number.otu.oligo)) #294

CC1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(CC1.18.number.otu.copio)) #141


## CC2
number.field <- "3"
CC2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(CC2.number.otu.oligo)) #266
CC2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(CC2.number.otu.copio)) #207

## CC2.18
number.field <- "4"
CC2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(CC2.18.number.otu.oligo)) #235

CC2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(CC2.18.number.otu.copio)) #155

## CJ1
number.field <- "5"
CJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(CJ1.number.otu.oligo)) #619
CJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(CJ1.number.otu.copio)) #184

## CJ1.18
number.field <- "6"
CJ1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(CJ1.18.number.otu.oligo)) #522
CJ1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(CJ1.18.number.otu.copio)) #150


## CJ2
number.field <- "7"
CJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(CJ2.number.otu.oligo)) #402
CJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(CJ2.number.otu.copio)) #166

## CJ2.18
number.field <- "8"
CJ2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(CJ2.18.number.otu.oligo)) #350
CJ2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(CJ2.18.number.otu.copio)) #135

## DG1
number.field <- "9"
DG1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(DG1.number.otu.oligo)) #195
DG1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(DG1.number.otu.copio)) #522

## DG1.18
number.field <- "10"
DG1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(DG1.18.number.otu.oligo)) #157
DG1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(DG1.18.number.otu.copio)) #374


## DG2
number.field <- "11"
DG2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(DG2.number.otu.oligo)) #159
DG2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(DG2.number.otu.copio)) #538

## DG2.18
number.field <- "12"
DG2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(DG2.18.number.otu.oligo)) #196
DG2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(DG2.18.number.otu.copio)) #425

## IS1
number.field <- "13"
IS1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(IS1.number.otu.oligo)) #162
IS1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(IS1.number.otu.copio)) #530


## IS1.18
number.field <- "14"
IS1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(IS1.18.number.otu.oligo)) #150
IS1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(IS1.18.number.otu.copio)) #348


## IS2
number.field <- "15"
IS2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(IS2.number.otu.oligo)) #199
IS2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(IS2.number.otu.copio)) #461

## IS2.18
number.field <- "16"
IS2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(IS2.18.number.otu.oligo)) #181
IS2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(IS2.18.number.otu.copio)) #336

## JJ1
number.field <- "17"
JJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(JJ1.number.otu.oligo)) #524
JJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(JJ1.number.otu.copio)) #164

## JJ2
number.field <- "18"
JJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(JJ2.number.otu.oligo)) #335
JJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(JJ2.number.otu.copio)) #232


## MY1
number.field <- "19"
MY1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(MY1.number.otu.oligo)) #239
MY1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(MY1.number.otu.copio)) #374

## MY2
number.field <- "20"
MY2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(MY2.number.otu.oligo)) #241
MY2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(MY2.number.otu.copio)) #311


## NJ1
number.field <- "21"
NJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(NJ1.number.otu.oligo)) #266
NJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(NJ1.number.otu.copio)) #261
## NJ2
number.field <- "22"
NJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(NJ2.number.otu.oligo)) #276
NJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(NJ2.number.otu.copio)) #320

## UF1
number.field <- "23"
UF1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(UF1.number.otu.oligo)) #271
UF1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(UF1.number.otu.copio)) #339

## UF1.18
number.field <- "24"
UF1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(UF1.18.number.otu.oligo)) #352
UF1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(UF1.18.number.otu.copio)) #187

## UF2
number.field <- "25"
UF2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(UF2.number.otu.oligo)) #578
UF2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(UF2.number.otu.copio)) #181

## UF2.18
number.field <- "26"
UF2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(UF2.18.number.otu.oligo)) #531
UF2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(UF2.18.number.otu.copio)) #147

## YS1
number.field <- "27"
YS1.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(YS1.number.otu.oligo)) #536
YS1.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(YS1.number.otu.copio)) #182

## YS2
number.field <- "28"
YS2.number.otu.oligo<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.oligo)
length(rownames(YS2.number.otu.oligo)) #262
YS2.number.otu.copio<-number_OTU_oligo_copio(phylo = bac.clean.ss, tab_oligo_copio = tab_ra.b.SOM.TN.copi)
length(rownames(YS2.number.otu.copio)) #206


##Archaea
##CJ1
number.field <- "1"
CC1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(CC1.number.otu.oligo)) #18
CC1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(CC1.number.otu.copio)) #10

## CC1.18
number.field <- "2"
CC1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(CC1.18.number.otu.oligo)) #24

CC1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(CC1.18.number.otu.copio)) #9


## CC2
number.field <- "3"
CC2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(CC2.number.otu.oligo)) #20
CC2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(CC2.number.otu.copio)) #13

## CC2.18
number.field <- "4"
CC2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(CC2.18.number.otu.oligo)) #18

CC2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(CC2.18.number.otu.copio)) #9

## CJ1
number.field <- "5"
CJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(CJ1.number.otu.oligo)) #32
CJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(CJ1.number.otu.copio)) #13

## CJ1.18
number.field <- "6"
CJ1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(CJ1.18.number.otu.oligo)) #31
CJ1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(CJ1.18.number.otu.copio)) #12


## CJ2
number.field <- "7"
CJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(CJ2.number.otu.oligo)) #24
CJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(CJ2.number.otu.copio)) #15

## CJ2.18
number.field <- "8"
CJ2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(CJ2.18.number.otu.oligo)) #22
CJ2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(CJ2.18.number.otu.copio)) #7

## DG1
number.field <- "9"
DG1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(DG1.number.otu.oligo)) #11
DG1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(DG1.number.otu.copio)) #50

## DG1.18
number.field <- "10"
DG1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(DG1.18.number.otu.oligo)) #11
DG1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(DG1.18.number.otu.copio)) #26


## DG2
number.field <- "11"
DG2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(DG2.number.otu.oligo)) #11
DG2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(DG2.number.otu.copio)) #51

## DG2.18
number.field <- "12"
DG2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(DG2.18.number.otu.oligo)) #11
DG2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(DG2.18.number.otu.copio)) #29

## IS1
number.field <- "13"
IS1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(IS1.number.otu.oligo)) #19
IS1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(IS1.number.otu.copio)) #42


## IS1.18
number.field <- "14"
IS1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(IS1.18.number.otu.oligo)) #18
IS1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(IS1.18.number.otu.copio)) #19


## IS2
number.field <- "15"
IS2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(IS2.number.otu.oligo)) #17
IS2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(IS2.number.otu.copio)) #36

## IS2.18
number.field <- "16"
IS2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(IS2.18.number.otu.oligo)) #13
IS2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(IS2.18.number.otu.copio)) #22

## JJ1
number.field <- "17"
JJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(JJ1.number.otu.oligo)) #25
JJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(JJ1.number.otu.copio)) #13

## JJ2
number.field <- "18"
JJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(JJ2.number.otu.oligo)) #19
JJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(JJ2.number.otu.copio)) #12


## MY1
number.field <- "19"
MY1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(MY1.number.otu.oligo)) #11
MY1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(MY1.number.otu.copio)) #22

## MY2
number.field <- "20"
MY2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(MY2.number.otu.oligo)) #7
MY2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(MY2.number.otu.copio)) #20


## NJ1
number.field <- "21"
NJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(NJ1.number.otu.oligo)) #17
NJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(NJ1.number.otu.copio)) #22

## NJ2
number.field <- "22"
NJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(NJ2.number.otu.oligo)) #17
NJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(NJ2.number.otu.copio)) #24

## UF1
number.field <- "23"
UF1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(UF1.number.otu.oligo)) #16
UF1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(UF1.number.otu.copio)) #20

## UF1.18
number.field <- "24"
UF1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(UF1.18.number.otu.oligo)) #19
UF1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(UF1.18.number.otu.copio)) #10

## UF2
number.field <- "25"
UF2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(UF2.number.otu.oligo)) #27
UF2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(UF2.number.otu.copio)) #14

## UF2.18
number.field <- "26"
UF2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(UF2.18.number.otu.oligo)) #24
UF2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(UF2.18.number.otu.copio)) #9

## YS1
number.field <- "27"
YS1.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(YS1.number.otu.oligo)) #31
YS1.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(YS1.number.otu.copio)) #19

## YS2
number.field <- "28"
YS2.number.otu.oligo<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.oligo)
length(rownames(YS2.number.otu.oligo)) #16
YS2.number.otu.copio<-number_OTU_oligo_copio(phylo = arch.clean.ss, tab_oligo_copio = tab_ra.a.SOM.TN.copi)
length(rownames(YS2.number.otu.copio)) #19



##Fungi
##CJ1
number.field <- "1"
CC1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(CC1.number.otu.oligo)) #83
CC1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(CC1.number.otu.copio)) #140

## CC1.18
number.field <- "2"
CC1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(CC1.18.number.otu.oligo)) #110

CC1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(CC1.18.number.otu.copio)) #65


## CC2
number.field <- "3"
CC2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(CC2.number.otu.oligo)) #76
CC2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(CC2.number.otu.copio)) #131

## CC2.18
number.field <- "4"
CC2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(CC2.18.number.otu.oligo)) #112

CC2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(CC2.18.number.otu.copio)) #70

## CJ1
number.field <- "5"
CJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(CJ1.number.otu.oligo)) #178
CJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(CJ1.number.otu.copio)) #109

## CJ1.18
number.field <- "6"
CJ1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(CJ1.18.number.otu.oligo)) #241
CJ1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(CJ1.18.number.otu.copio)) #67


## CJ2
number.field <- "7"
CJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(CJ2.number.otu.oligo)) #105
CJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(CJ2.number.otu.copio)) #119

## CJ2.18
number.field <- "8"
CJ2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(CJ2.18.number.otu.oligo)) #136
CJ2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(CJ2.18.number.otu.copio)) #59

## DG1
number.field <- "9"
DG1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(DG1.number.otu.oligo)) #39
DG1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(DG1.number.otu.copio)) #287

## DG1.18
number.field <- "10"
DG1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(DG1.18.number.otu.oligo)) #88
DG1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(DG1.18.number.otu.copio)) #153


## DG2
number.field <- "11"
DG2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(DG2.number.otu.oligo)) #46
DG2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(DG2.number.otu.copio)) #310

## DG2.18
number.field <- "12"
DG2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(DG2.18.number.otu.oligo)) #101
DG2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(DG2.18.number.otu.copio)) #228

## IS1
number.field <- "13"
IS1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(IS1.number.otu.oligo)) #62
IS1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(IS1.number.otu.copio)) #291


## IS1.18
number.field <- "14"
IS1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(IS1.18.number.otu.oligo)) #96
IS1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(IS1.18.number.otu.copio)) #141


## IS2
number.field <- "15"
IS2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(IS2.number.otu.oligo)) #76
IS2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(IS2.number.otu.copio)) #272

## IS2.18
number.field <- "16"
IS2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(IS2.18.number.otu.oligo)) #109
IS2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(IS2.18.number.otu.copio)) #146

## JJ1
number.field <- "17"
JJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(JJ1.number.otu.oligo)) #171
JJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(JJ1.number.otu.copio)) #130

## JJ2
number.field <- "18"
JJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(JJ2.number.otu.oligo)) #111
JJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(JJ2.number.otu.copio)) #189


## MY1
number.field <- "19"
MY1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(MY1.number.otu.oligo)) #53
MY1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(MY1.number.otu.copio)) #235

## MY2
number.field <- "20"
MY2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(MY2.number.otu.oligo)) #67
MY2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(MY2.number.otu.copio)) #183


## NJ1
number.field <- "21"
NJ1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(NJ1.number.otu.oligo)) #89
NJ1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(NJ1.number.otu.copio)) #207

## NJ2
number.field <- "22"
NJ2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(NJ2.number.otu.oligo)) #78
NJ2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(NJ2.number.otu.copio)) #205

## UF1
number.field <- "23"
UF1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(UF1.number.otu.oligo)) #66
UF1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(UF1.number.otu.copio)) #195

## UF1.18
number.field <- "24"
UF1.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(UF1.18.number.otu.oligo)) #168
UF1.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(UF1.18.number.otu.copio)) #69

## UF2
number.field <- "25"
UF2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(UF2.number.otu.oligo)) #133
UF2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(UF2.number.otu.copio)) #147

## UF2.18
number.field <- "26"
UF2.18.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(UF2.18.number.otu.oligo)) #24
UF2.18.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(UF2.18.number.otu.copio)) #9

## YS1
number.field <- "27"
YS1.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(YS1.number.otu.oligo)) #234
YS1.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(YS1.number.otu.copio)) #149

## YS2
number.field <- "28"
YS2.number.otu.oligo<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.oligo)
length(rownames(YS2.number.otu.oligo)) #108
YS2.number.otu.copio<-number_OTU_oligo_copio(phylo = fun.clean.ss, tab_oligo_copio = tab_ra.f.SOM.TN.copi)
length(rownames(YS2.number.otu.copio)) #181




## Differences of numbers of putative oligo and copiotrophs between the fields

head(ps1.meta)

otu.bac.clean.ss <- otu_table(bac.clean.ss)
otu.bac.clean.ss <- data.frame(otu.bac.clean.ss)
otu.bac.clean.ss.oligo <- subset(otu.bac.clean.ss , rownames(otu.bac.clean.ss) %in% rownames(tab_ra.b.SOM.TN.oligo))
otu.bac.clean.ss.copio <- subset(otu.bac.clean.ss , rownames(otu.bac.clean.ss) %in% rownames(tab_ra.b.SOM.TN.copi))

library(tidyverse)

bac.oligo.each.sample<- otu.bac.clean.ss.oligo%>%
  gather(x, value, CC1.A.1.1:YS2.C.3.3)%>%
  group_by(x)%>%
  tally(value > 0)

names(bac.oligo.each.sample)[1] <- "SampleID"
names(bac.oligo.each.sample)[2] <- "Count_oligo"

bac.copio.each.sample<- otu.bac.clean.ss.copio%>%
  gather(x, value, CC1.A.1.1:YS2.C.3.3)%>%
  group_by(x)%>%
  tally(value > 0)

names(bac.copio.each.sample)[1] <- "SampleID"
names(bac.copio.each.sample)[2] <- "Count_copio"

ps1.oli<-merge(ps1.meta, bac.oligo.each.sample, by = "SampleID")
ps1.oli.copi<-merge(ps1.oli, bac.copio.each.sample, by = "SampleID")
head(ps1.oli.copi)


##Archaea

head(ps2.meta)
otu.arch.clean.ss <- otu_table(arch.clean.ss)
otu.arch.clean.ss <- data.frame(otu.arch.clean.ss)
otu.arch.clean.ss.oligo <- subset(otu.arch.clean.ss , rownames(otu.arch.clean.ss) %in% rownames(tab_ra.a.SOM.TN.oligo))
otu.arch.clean.ss.copio <- subset(otu.arch.clean.ss , rownames(otu.arch.clean.ss) %in% rownames(tab_ra.a.SOM.TN.copi))

library(tidyverse)

arch.oligo.each.sample<- otu.arch.clean.ss.oligo%>%
  gather(x, value, CC1.A.1.1:YS2.C.3.3)%>%
  group_by(x)%>%
  tally(value > 0)

names(arch.oligo.each.sample)[1] <- "SampleID"
names(arch.oligo.each.sample)[2] <- "Count_oligo"

arch.copio.each.sample<- otu.arch.clean.ss.copio%>%
  gather(x, value, CC1.A.1.1:YS2.C.3.3)%>%
  group_by(x)%>%
  tally(value > 0)

names(arch.copio.each.sample)[1] <- "SampleID"
names(arch.copio.each.sample)[2] <- "Count_copio"

ps2.oli<-merge(ps2.meta, arch.oligo.each.sample, by = "SampleID")
ps2.oli.copi<-merge(ps2.oli, arch.copio.each.sample, by = "SampleID")
head(ps2.oli.copi)



##Fungi

head(ps3.meta)
otu.fun.clean.ss <- otu_table(fun.clean.ss)
otu.fun.clean.ss <- data.frame(otu.fun.clean.ss)
otu.fun.clean.ss.oligo <- subset(otu.fun.clean.ss , rownames(otu.fun.clean.ss) %in% rownames(tab_ra.f.SOM.TN.oligo))
otu.fun.clean.ss.copio <- subset(otu.fun.clean.ss , rownames(otu.fun.clean.ss) %in% rownames(tab_ra.f.SOM.TN.copi))

library(tidyverse)

fun.oligo.each.sample<- otu.fun.clean.ss.oligo%>%
  gather(x, value, CC1.A.1.1:YS2.C.3.3)%>%
  group_by(x)%>%
  tally(value > 0)

names(fun.oligo.each.sample)[1] <- "SampleID"
names(fun.oligo.each.sample)[2] <- "Count_oligo"

fun.copio.each.sample<- otu.fun.clean.ss.copio%>%
  gather(x, value, CC1.A.1.1:YS2.C.3.3)%>%
  group_by(x)%>%
  tally(value > 0)

names(fun.copio.each.sample)[1] <- "SampleID"
names(fun.copio.each.sample)[2] <- "Count_copio"

ps3.oli<-merge(ps3.meta, fun.oligo.each.sample, by = "SampleID")
ps3.oli.copi<-merge(ps3.oli, fun.copio.each.sample, by = "SampleID")
head(ps3.oli.copi)

library(ggsignif)
x <- subset(ps1.oli.copi, Field=='CJ1')$Count_oligo
y <- subset(ps1.oli.copi, Field=='CJ2')$Count_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value
p <- ggplot(data = subset(ps1.oli.copi, Field %in% c("CJ1", "CJ2")), aes(x=Field, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Count of putative oligotrophs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("CJ1", "CJ2")), step_increase = 0.1, map_signif_level=TRUE)

p


x <- subset(ps1.oli.copi, Field=='CJ1')$Count_copio
y <- subset(ps1.oli.copi, Field=='CJ2')$Count_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value
p <- ggplot(data = subset(ps1.oli.copi, Field %in% c("CJ1", "CJ2")), aes(x=Field, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Count of putative copiotrophs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("CJ1", "CJ2")), step_increase = 0.1, map_signif_level=TRUE)

p



x <- subset(ps1.oli.copi, Field=='YS1')$Count_oligo
y <- subset(ps1.oli.copi, Field=='YS2')$Count_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value
p <- ggplot(data = subset(ps1.oli.copi, Field %in% c("YS1", "YS2")), aes(x=Field, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Count of putative oligotrophs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("YS1", "YS2")), step_increase = 0.1, map_signif_level=TRUE)

p


x <- subset(ps1.oli.copi, Field=='YS1')$Count_copio
y <- subset(ps1.oli.copi, Field=='YS2')$Count_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value
p <- ggplot(data = subset(ps1.oli.copi, Field %in% c("YS1", "YS2")), aes(x=Field, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Count of putative copiotrophs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("YS1", "YS2")), step_increase = 0.1, map_signif_level=TRUE)

p





x <- subset(ps1.oli.copi, Field=='MY1')$Count_oligo
y <- subset(ps1.oli.copi, Field=='MY2')$Count_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value
p <- ggplot(data = subset(ps1.oli.copi, Field %in% c("MY1", "MY2")), aes(x=Field, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Count of putative oligotrophs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("MY1", "MY2")), step_increase = 0.1, map_signif_level=TRUE)

p


x <- subset(ps1.oli.copi, Field=='MY1')$Count_copio
y <- subset(ps1.oli.copi, Field=='MY2')$Count_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value
p <- ggplot(data = subset(ps1.oli.copi, Field %in% c("MY1", "MY2")), aes(x=Field, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Count of putative copiotrophs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("MY1", "MY2")), step_increase = 0.1, map_signif_level=TRUE)

p


x <- subset(ps1.oli.copi, Field=='NJ1')$Count_oligo
y <- subset(ps1.oli.copi, Field=='NJ2')$Count_oligo
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value
p <- ggplot(data = subset(ps1.oli.copi, Field %in% c("NJ1", "NJ2")), aes(x=Field, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Count of putative oligotrophs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("NJ1", "NJ2")), step_increase = 0.1, map_signif_level=TRUE)

p


x <- subset(ps1.oli.copi, Field=='NJ1')$Count_copio
y <- subset(ps1.oli.copi, Field=='NJ2')$Count_copio
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value
p <- ggplot(data = subset(ps1.oli.copi, Field %in% c("NJ1", "NJ2")), aes(x=Field, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Count of putative copiotrophs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("NJ1", "NJ2")), step_increase = 0.1, map_signif_level=TRUE)

p


##Bacteria

##Oligo
ps1.meta.trophic$Field <- factor(ps1.meta.trophic$Field, levels = lb)
test <- aggregate(ps1.meta.trophic$bac_RA_oligo, by = list(ps1.meta.trophic$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(bac_RA_oligo ~ Field, data = ps1.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(bac_RA_oligo ~ Field,
              data=ps1.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.meta.trophic, aes(x=Field, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps1.meta.trophic$bac_RA_copio, by = list(ps1.meta.trophic$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(bac_RA_copio ~ Field, data = ps1.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(bac_RA_copio ~ Field,
              data=ps1.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.meta.trophic, aes(x=Field, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


##Archaea
##Oligo
ps2.meta.trophic$Field <- factor(ps2.meta.trophic$Field, levels = lb)
test <- aggregate(ps2.meta.trophic$arch_RA_oligo, by = list(ps2.meta.trophic$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(arch_RA_oligo ~ Field, data = ps2.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(arch_RA_oligo ~ Field,
              data=ps2.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.meta.trophic, aes(x=Field, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps2.meta.trophic$arch_RA_copio, by = list(ps2.meta.trophic$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(arch_RA_copio ~ Field, data = ps2.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(arch_RA_copio ~ Field,
              data=ps2.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.meta.trophic, aes(x=Field, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Fungi
##Oligo
ps3.meta.trophic$Field <- factor(ps3.meta.trophic$Field, levels = lb)
test <- aggregate(ps3.meta.trophic$fun_RA_oligo, by = list(ps3.meta.trophic$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(fun_RA_oligo ~ Field, data = ps3.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(fun_RA_oligo ~ Field,
              data=ps3.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.meta.trophic, aes(x=Field, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps3.meta.trophic$fun_RA_copio, by = list(ps3.meta.trophic$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(fun_RA_copio ~ Field, data = ps3.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(fun_RA_copio ~ Field,
              data=ps3.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.meta.trophic, aes(x=Field, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p




### OTU count
ps1.oli.copi$Field <- factor(ps1.oli.copi$Field, levels = lb)
test <- aggregate(ps1.oli.copi$Count_oligo, by = list(ps1.oli.copi$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_oligo ~ Field, data = ps1.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_oligo ~ Field,
              data=ps1.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.oli.copi, aes(x=Field, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps1.oli.copi$Count_copio, by = list(ps1.oli.copi$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_copio ~ Field, data = ps1.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_copio ~ Field,
              data=ps1.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.oli.copi, aes(x=Field, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


##Archaea
##Oligo
ps2.oli.copi$Field <- factor(ps2.oli.copi$Field, levels = lb)
test <- aggregate(ps2.oli.copi$Count_oligo, by = list(ps2.oli.copi$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_oligo ~ Field, data = ps2.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_oligo ~ Field,
              data=ps2.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.oli.copi, aes(x=Field, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps2.oli.copi$Count_copio, by = list(ps2.oli.copi$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_copio ~ Field, data = ps2.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_copio ~ Field,
              data=ps2.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.oli.copi, aes(x=Field, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Fungi
##Oligo
ps3.oli.copi$Field <- factor(ps3.oli.copi$Field, levels = lb)
test <- aggregate(ps3.oli.copi$Count_oligo, by = list(ps3.oli.copi$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_oligo ~ Field, data = ps3.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_oligo ~ Field,
              data=ps3.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.oli.copi, aes(x=Field, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps3.oli.copi$Count_copio, by = list(ps3.oli.copi$Field), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_copio ~ Field, data = ps3.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_copio ~ Field,
              data=ps3.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.oli.copi, aes(x=Field, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p



### Effect of cultural practices
##Bacteria

##Oligo
lb <- c("Conventional", "No_fertilizer", "No_pesticide")
ps1.meta.trophic$Cultural_practice <- factor(ps1.meta.trophic$Cultural_practice, levels = lb)
test <- aggregate(ps1.meta.trophic$bac_RA_oligo, by = list(ps1.meta.trophic$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(bac_RA_oligo ~ Cultural_practice, data = ps1.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(bac_RA_oligo ~ Cultural_practice,
              data=ps1.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.meta.trophic, aes(x=Cultural_practice, y=bac_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps1.meta.trophic$bac_RA_copio, by = list(ps1.meta.trophic$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(bac_RA_copio ~ Cultural_practice, data = ps1.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(bac_RA_copio ~ Cultural_practice,
              data=ps1.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.meta.trophic, aes(x=Cultural_practice, y=bac_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


##Archaea
##Oligo
ps2.meta.trophic$Cultural_practice <- factor(ps2.meta.trophic$Cultural_practice, levels = lb)
test <- aggregate(ps2.meta.trophic$arch_RA_oligo, by = list(ps2.meta.trophic$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(arch_RA_oligo ~ Cultural_practice, data = ps2.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(arch_RA_oligo ~ Cultural_practice,
              data=ps2.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.meta.trophic, aes(x=Cultural_practice, y=arch_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps2.meta.trophic$arch_RA_copio, by = list(ps2.meta.trophic$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(arch_RA_copio ~ Cultural_practice, data = ps2.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(arch_RA_copio ~ Cultural_practice,
              data=ps2.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.meta.trophic, aes(x=Cultural_practice, y=arch_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Fungi
##Oligo
ps3.meta.trophic$Cultural_practice <- factor(ps3.meta.trophic$Cultural_practice, levels = lb)
test <- aggregate(ps3.meta.trophic$fun_RA_oligo, by = list(ps3.meta.trophic$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(fun_RA_oligo ~ Cultural_practice, data = ps3.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(fun_RA_oligo ~ Cultural_practice,
              data=ps3.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.meta.trophic, aes(x=Cultural_practice, y=fun_RA_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps3.meta.trophic$fun_RA_copio, by = list(ps3.meta.trophic$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(fun_RA_copio ~ Cultural_practice, data = ps3.meta.trophic)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(fun_RA_copio ~ Cultural_practice,
              data=ps3.meta.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.meta.trophic, aes(x=Cultural_practice, y=fun_RA_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Relative abundance (%) \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p




### OTU count
ps1.oli.copi$Cultural_practice <- factor(ps1.oli.copi$Cultural_practice, levels = lb)
test <- aggregate(ps1.oli.copi$Count_oligo, by = list(ps1.oli.copi$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_oligo ~ Cultural_practice, data = ps1.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_oligo ~ Cultural_practice,
              data=ps1.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.oli.copi, aes(x=Cultural_practice, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps1.oli.copi$Count_copio, by = list(ps1.oli.copi$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_copio ~ Cultural_practice, data = ps1.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_copio ~ Cultural_practice,
              data=ps1.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps1.oli.copi, aes(x=Cultural_practice, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


##Archaea
##Oligo
ps2.oli.copi$Cultural_practice <- factor(ps2.oli.copi$Cultural_practice, levels = lb)
test <- aggregate(ps2.oli.copi$Count_oligo, by = list(ps2.oli.copi$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_oligo ~ Cultural_practice, data = ps2.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_oligo ~ Cultural_practice,
              data=ps2.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.oli.copi, aes(x=Cultural_practice, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps2.oli.copi$Count_copio, by = list(ps2.oli.copi$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_copio ~ Cultural_practice, data = ps2.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_copio ~ Cultural_practice,
              data=ps2.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps2.oli.copi, aes(x=Cultural_practice, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Fungi
##Oligo
ps3.oli.copi$Cultural_practice <- factor(ps3.oli.copi$Cultural_practice, levels = lb)
test <- aggregate(ps3.oli.copi$Count_oligo, by = list(ps3.oli.copi$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_oligo ~ Cultural_practice, data = ps3.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_oligo ~ Cultural_practice,
              data=ps3.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.oli.copi, aes(x=Cultural_practice, y=Count_oligo)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p

##Copio
test <- aggregate(ps3.oli.copi$Count_copio, by = list(ps3.oli.copi$Cultural_practice), max)

colnames(test) <- c("Group", "MaxAbund")

##Kruskal-Wallis test
kw<-kruskal.test(Count_copio ~ Cultural_practice, data = ps3.oli.copi)
kw$p.value
kw$p.value<- round(kw$p.value, 4)

#library(FSA)
DT = dunnTest(Count_copio ~ Cultural_practice,
              data=ps3.oli.copi,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
p <- ggplot(data =ps3.oli.copi, aes(x=Cultural_practice, y=Count_copio)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
  ylab("Number of OTUs \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")


p


dev.off()


## Supplementary Table
tab_ra.a.SOM.TN.oligo$Category <- "Oligotroph"
tab_ra.a.SOM.TN.copi$Category <- "Copiotroph"

tab_ra.oligo.copio <- rbind(tab_ra.a.SOM.TN.oligo, tab_ra.a.SOM.TN.copi)

tab_ra.oligo.copio <- tab_ra.oligo.copio%>% select('Category')
tab_ra.oligo.copio$OTU <- rownames(tab_ra.oligo.copio)

arch.clean.field <- merge_samples(arch.clean.ss, "Field")
arch.clean.field.rel <- microbiome::transform(arch.clean.field, transform = "compositional")
arch.copio.oligo<- t(otu_table(arch.clean.field.rel))
arch.copio.oligo <- data.frame(arch.copio.oligo)
arch.copio.oligo$OTU <- rownames(arch.copio.oligo)

arch.copio.oligo <- merge(tab_ra.oligo.copio,arch.copio.oligo, by = "OTU" )

otu.list <- rbind(phy.list, arch.list, fun.list)

arch.copio.oligo.id <- merge(otu.list, arch.copio.oligo, by ="OTU")
head(arch.copio.oligo.id)

write.xlsx(arch.copio.oligo.id , "Archaea_oligo and copio.xlsx")


tab_ra.b.SOM.TN.oligo$Category <- "Oligotroph"
tab_ra.b.SOM.TN.copi$Category <- "Copiotroph"

tab_ra.oligo.copio <- rbind(tab_ra.b.SOM.TN.oligo, tab_ra.b.SOM.TN.copi)

tab_ra.oligo.copio <- tab_ra.oligo.copio%>% select('Category')
tab_ra.oligo.copio$OTU <- rownames(tab_ra.oligo.copio)

bac.clean.field <- merge_samples(bac.clean.ss, "Field")
bac.clean.field.rel <- microbiome::transform(bac.clean.field, transform = "compositional")
bac.copio.oligo<- t(otu_table(bac.clean.field.rel))
bac.copio.oligo <- data.frame(bac.copio.oligo)
bac.copio.oligo$OTU <- rownames(bac.copio.oligo)

bac.copio.oligo <- merge(tab_ra.oligo.copio,bac.copio.oligo, by = "OTU" )

otu.list <- rbind(phy.list, bac.list, fun.list)

bac.copio.oligo.id <- merge(otu.list, bac.copio.oligo, by ="OTU")
head(bac.copio.oligo.id)

write.xlsx(bac.copio.oligo.id , "Bacteria_oligo and copio.xlsx")


tab_ra.f.SOM.TN.oligo$Category <- "Oligotroph"
tab_ra.f.SOM.TN.copi$Category <- "Copiotroph"

tab_ra.oligo.copio <- rbind(tab_ra.f.SOM.TN.oligo, tab_ra.f.SOM.TN.copi)

tab_ra.oligo.copio <- tab_ra.oligo.copio%>% select('Category')
tab_ra.oligo.copio$OTU <- rownames(tab_ra.oligo.copio)

fun.clean.field <- merge_samples(fun.clean.ss, "Field")
fun.clean.field.rel <- microbiome::transform(fun.clean.field, transform = "compositional")
fun.copio.oligo<- t(otu_table(fun.clean.field.rel))
fun.copio.oligo <- data.frame(fun.copio.oligo)
fun.copio.oligo$OTU <- rownames(fun.copio.oligo)

fun.copio.oligo <- merge(tab_ra.oligo.copio,fun.copio.oligo, by = "OTU" )

otu.list <- rbind(phy.list, fun.list, fun.list)

fun.copio.oligo.id <- merge(otu.list, fun.copio.oligo, by ="OTU")
head(fun.copio.oligo.id)

write.xlsx(fun.copio.oligo.id , "Fungi_oligo and copio.xlsx")

