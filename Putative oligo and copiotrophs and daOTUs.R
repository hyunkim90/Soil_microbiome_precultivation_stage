### Differentially abundant OTUs between fertilized and non-fertilized conditions
## In oligotrophic condition (CJ1,)
## Defining daOTUs
#Phyloseq files
bac.clean.ss
arch.clean.ss


fun.clean.ss



bac.clean.ss.CJ.YS <- subset_samples(bac.clean.ss, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
bac.clean.ss.CJ.YS <- phyloseq::filter_taxa(bac.clean.ss.CJ.YS, function(x) sum(x) != 0, TRUE)

arch.clean.ss.CJ.YS <- subset_samples(arch.clean.ss, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
arch.clean.ss.CJ.YS <- phyloseq::filter_taxa(arch.clean.ss.CJ.YS, function(x) sum(x) != 0, TRUE)

fun.clean.ss.CJ.YS <- subset_samples(fun.clean.ss, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
fun.clean.ss.CJ.YS <- phyloseq::filter_taxa(fun.clean.ss.CJ.YS, function(x) sum(x) != 0, TRUE)


zero.clean.filt <- phyloseq::filter_taxa(fun.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(zero.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(zero.clean.filt))

## make binary samples
zero.clean.filt.1 <- subset_samples(zero.clean.filt, Field %in% c('CJ1','CJ2'))
zero.clean.filt.2 <- subset_samples(zero.clean.filt, Field %in% c('CJ1','CJ2'))
zero.clean.filt.3 <- subset_samples(zero.clean.filt, Field %in% c('CJ1','CJ2'))

zero.clean.filt.1 <- phyloseq::filter_taxa(zero.clean.filt.1, function(x) sum(x) != 0, TRUE)
zero.clean.filt.2 <- phyloseq::filter_taxa(zero.clean.filt.2, function(x) sum(x) != 0, TRUE)
zero.clean.filt.3 <- phyloseq::filter_taxa(zero.clean.filt.3, function(x) sum(x) != 0, TRUE)

sample_names(zero.clean.filt.1)
# only keep samples over 1000 read
(filt.sample <- sample_sums(zero.clean.filt.1) > 100)
sum(sample_sums(zero.clean.filt.1) <= 100)  ## 0 sample discarded
zero.clean.filt.1 <- prune_samples(filt.sample, zero.clean.filt.1)
## 40x difference in sequencing depth between samples. 
#CSS normalization is not going to ??fix?? this much difference. 
#In this case I would recommend rarefaction, but we will move forward since this is the example dataset in phyloseq
obj <- phyloseq_to_metagenomeSeq(zero.clean.filt.3)

# obj <- filterData(obj, depth = 100)
obj   ## removed 0 samples below 100 reads
##
## fitZig sample
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 30, verbose = TRUE)
Type  <-  pData(obj)$Field
mod  <-  model.matrix(~Type)
colnames(mod)  <-  levels(Type)
colnames(mod)

res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)


zigFit = res$fit
finalMod= res$fit$design


contrast.matrix.1 =makeContrasts(CJ1 - CJ2, levels = finalMod)


fit2_nonfer_fer = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 

fit3.1 = eBayes(fit2_nonfer_fer)

fit3.2 = eBayes(fit2_nonfer_fer)

fit3.3 = eBayes(fit2_nonfer_fer)

topTable(fit3.1, coef="CJ1 - CJ2")
topTable(fit3.2, coef="CJ1 - CJ2")
topTable(fit3.3, coef="CJ1 - CJ2")
res.1 <- topTable(fit3.1,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.1)
dim(res.1)

res.2 <- topTable(fit3.2,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.2)
dim(res.2)

res.3 <- topTable(fit3.3,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.3)
dim(res.3)

write.xlsx(res.1, 'Bac_daOTU_CJ1 and CJ2.xlsx')

write.xlsx(res.2, 'Arch_daOTU_CJ1 and CJ2.xlsx')

write.xlsx(res.3, 'Fun_daOTU_CJ1 and CJ2.xlsx')

res.1$OTU <- rownames(res.1)

res.2$OTU <- rownames(res.2)

res.3$OTU <- rownames(res.3)


res.CJ1.CJ2.bac <- merge(res.1, OTU_id.list, by= "OTU")
res.CJ1.CJ2.arch <- merge(res.2, OTU_id.list, by= "OTU")
res.CJ1.CJ2.fun<- merge(res.3, OTU_id.list, by= "OTU")

write.xlsx(res.CJ1.CJ2.bac, 'Bac_daOTU_CJ1 and CJ2.xlsx')

write.xlsx(res.CJ1.CJ2.arch, 'Arch_daOTU_CJ1 and CJ2.xlsx')

write.xlsx(res.CJ1.CJ2.fun, 'Fun_daOTU_CJ1 and CJ2.xlsx')

res.CJ1_enriched.bac <- res.CJ1.CJ2.bac$OTU_id[which(res.CJ1.CJ2.bac$logFC>0)]
res.CJ2_enriched.bac <- res.CJ1.CJ2.bac$OTU_id[which(res.CJ1.CJ2.bac$logFC<0)]

res.CJ1_enriched.arch <- res.CJ1.CJ2.arch$OTU_id[which(res.CJ1.CJ2.arch$logFC>0)]
res.CJ2_enriched.arch <- res.CJ1.CJ2.arch$OTU_id[which(res.CJ1.CJ2.arch$logFC<0)]

res.CJ1_enriched.fun <- res.CJ1.CJ2.fun$OTU_id[which(res.CJ1.CJ2.fun$logFC>0)]
res.CJ2_enriched.fun <- res.CJ1.CJ2.fun$OTU_id[which(res.CJ1.CJ2.fun$logFC<0)]




## In Copiotrophic condition
zero.clean.filt <- phyloseq::filter_taxa(fun.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(zero.clean.filt) == 0)

zero.clean.filt.1 <- subset_samples(zero.clean.filt, Field %in% c('MY1','MY2'))
zero.clean.filt.2 <- subset_samples(zero.clean.filt, Field %in% c('MY1','MY2'))
zero.clean.filt.3 <- subset_samples(zero.clean.filt, Field %in% c('MY1','MY2'))

zero.clean.filt.1 <- phyloseq::filter_taxa(zero.clean.filt.1, function(x) sum(x) != 0, TRUE)
zero.clean.filt.2 <- phyloseq::filter_taxa(zero.clean.filt.2, function(x) sum(x) != 0, TRUE)
zero.clean.filt.3 <- phyloseq::filter_taxa(zero.clean.filt.3, function(x) sum(x) != 0, TRUE)

sample_names(zero.clean.filt.1)
# only keep samples over 1000 read
(filt.sample <- sample_sums(zero.clean.filt.1) > 100)
sum(sample_sums(zero.clean.filt.1) <= 100)  ## 0 sample discarded
zero.clean.filt.1 <- prune_samples(filt.sample, zero.clean.filt.1)
## 40x difference in sequencing depth between samples. 
#CSS normalization is not going to ??fix?? this much difference. 
#In this case I would recommend rarefaction, but we will move forward since this is the example dataset in phyloseq
obj <- phyloseq_to_metagenomeSeq(zero.clean.filt.3)

# obj <- filterData(obj, depth = 100)
obj   ## removed 0 samples below 100 reads
##
## fitZig sample
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 30, verbose = TRUE)
Type  <-  pData(obj)$Field
mod  <-  model.matrix(~Type)
colnames(mod)  <-  levels(Type)
colnames(mod)

res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)


zigFit = res$fit
finalMod= res$fit$design


contrast.matrix.1 =makeContrasts(MY1 - MY2, levels = finalMod)


fit2_nonfer_fer = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 

fit3.1 = eBayes(fit2_nonfer_fer)

fit3.2 = eBayes(fit2_nonfer_fer)

fit3.3 = eBayes(fit2_nonfer_fer)

topTable(fit3.1, coef="MY1 - MY2")
topTable(fit3.2, coef="MY1 - MY2")
topTable(fit3.3, coef="MY1 - MY2")

res.1 <- topTable(fit3.1,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.1)
dim(res.1)

res.2 <- topTable(fit3.2,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.2)
dim(res.2)

res.3 <- topTable(fit3.3,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.3)
dim(res.3)

res.1$OTU <- rownames(res.1)

res.2$OTU <- rownames(res.2)

res.3$OTU <- rownames(res.3)


res.MY1.MY2.bac <- merge(res.1, OTU_id.list, by= "OTU")
res.MY1.MY2.arch <- merge(res.2, OTU_id.list, by= "OTU")
res.MY1.MY2.fun<- merge(res.3, OTU_id.list, by= "OTU")

write.xlsx(res.MY1.MY2.bac, 'Bac_daOTU_MY1 and MY2.xlsx')

write.xlsx(res.MY1.MY2.arch, 'Arch_daOTU_MY1 and MY2.xlsx')

write.xlsx(res.MY1.MY2.fun, 'Fun_daOTU_MY1 and MY2.xlsx')

res.MY1_enriched.bac <- res.MY1.MY2.bac$OTU_id[which(res.MY1.MY2.bac$logFC>0)]
res.MY2_enriched.bac <- res.MY1.MY2.bac$OTU_id[which(res.MY1.MY2.bac$logFC<0)]

res.MY1_enriched.arch <- res.MY1.MY2.arch$OTU_id[which(res.MY1.MY2.arch$logFC>0)]
res.MY2_enriched.arch <- res.MY1.MY2.arch$OTU_id[which(res.MY1.MY2.arch$logFC<0)]

res.MY1_enriched.fun <- res.MY1.MY2.fun$OTU_id[which(res.MY1.MY2.fun$logFC>0)]
res.MY2_enriched.fun <- res.MY1.MY2.fun$OTU_id[which(res.MY1.MY2.fun$logFC<0)]

bac.clean.ss.CJ12 <- subset_samples(bac.clean.ss, Field %in% c("CJ1", "CJ2"))
bac.clean.ss.CJ12 <- phyloseq::filter_taxa(bac.clean.ss.CJ12, function(x) sum(x) != 0, TRUE)

arch.clean.ss.CJ12 <- subset_samples(arch.clean.ss, Field %in% c("CJ1", "CJ2"))
arch.clean.ss.CJ12 <- phyloseq::filter_taxa(arch.clean.ss.CJ12, function(x) sum(x) != 0, TRUE)

fun.clean.ss.CJ12 <- subset_samples(fun.clean.ss, Field %in% c("CJ1", "CJ2"))
fun.clean.ss.CJ12 <- phyloseq::filter_taxa(fun.clean.ss.CJ12, function(x) sum(x) != 0, TRUE)

bac.clean.ss.MY12 <- subset_samples(bac.clean.ss, Field %in% c("MY1", "MY2"))
bac.clean.ss.MY12 <- phyloseq::filter_taxa(bac.clean.ss.MY12, function(x) sum(x) != 0, TRUE)

arch.clean.ss.MY12 <- subset_samples(arch.clean.ss, Field %in% c("MY1", "MY2"))
arch.clean.ss.MY12 <- phyloseq::filter_taxa(arch.clean.ss.MY12, function(x) sum(x) != 0, TRUE)

fun.clean.ss.MY12 <- subset_samples(fun.clean.ss, Field %in% c("MY1", "MY2"))
fun.clean.ss.MY12 <- phyloseq::filter_taxa(fun.clean.ss.MY12, function(x) sum(x) != 0, TRUE)


get_resSig_CJ <- function(phy.clean.ss.5){
  phy.clean.filt <- phyloseq::filter_taxa(phy.clean.ss.5, function(x) sum(x) != 0, TRUE)
  sum(taxa_sums(phy.clean.filt) == 0)
  
  obj <- phyloseq_to_metagenomeSeq(phy.clean.filt)
  obj <-  cumNorm(obj, p = cumNormStatFast(obj))
  normFactor <-  normFactors(obj)
  normFactor <-  log2(normFactor/median(normFactor) + 1)
  settings <-  zigControl(maxit = 10, verbose = TRUE)
  Type  <-  pData(obj)$Field
  mod  <-  model.matrix(~Type)
  colnames(mod)  <-  levels(Type)
  colnames(mod)
  
  res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)
  res
  res$fit
  zigFit = res$fit
  finalMod= res$fit$design
  finalMod
  contrast.matrix =makeContrasts(CJ1 - CJ2, levels = finalMod)
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit2
  fit3 = eBayes(fit2)
  fit3
  topTable(fit3, coef="CJ1 - CJ2")
  
  
  res <- topTable(fit3,coef=1,adjust="fdr",n=Inf)
  
  log2AverageAbundance <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
  log2AverageAbundance
  Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
  Ta <- unique(Ta)
  Ta
  resSig = res[!is.na(res$adj.P.Val), ]
  resSig = data.frame(resSig)
  head(resSig)
  resSig <- tibble::rownames_to_column(resSig, 'OTU')
  resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
  resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))
  return(resSig)
}

get_resSig_MY <- function(phy.clean.ss.5){
  phy.clean.filt <- phyloseq::filter_taxa(phy.clean.ss.5, function(x) sum(x) != 0, TRUE)
  sum(taxa_sums(phy.clean.filt) == 0)
  
  obj <- phyloseq_to_metagenomeSeq(phy.clean.filt)
  obj <-  cumNorm(obj, p = cumNormStatFast(obj))
  normFactor <-  normFactors(obj)
  normFactor <-  log2(normFactor/median(normFactor) + 1)
  settings <-  zigControl(maxit = 10, verbose = TRUE)
  Type  <-  pData(obj)$Field
  mod  <-  model.matrix(~Type)
  colnames(mod)  <-  levels(Type)
  colnames(mod)
  
  res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)
  res
  res$fit
  zigFit = res$fit
  finalMod= res$fit$design
  finalMod
  contrast.matrix =makeContrasts(MY2 - MY1, levels = finalMod)
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit2
  fit3 = eBayes(fit2)
  fit3
  topTable(fit3, coef="MY2 - MY1")
  
  
  res <- topTable(fit3,coef=1,adjust="fdr",n=Inf)
  
  log2AverageAbundance <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
  log2AverageAbundance
  Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
  Ta <- unique(Ta)
  Ta
  resSig = res[!is.na(res$adj.P.Val), ]
  resSig = data.frame(resSig)
  head(resSig)
  resSig <- tibble::rownames_to_column(resSig, 'OTU')
  resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
  resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))
  return(resSig)
}


resSig_CJ.bac<-get_resSig_CJ(bac.clean.ss.CJ12)
resSig_CJ.arch<-get_resSig_CJ(arch.clean.ss.CJ12)
resSig_CJ.fun<-get_resSig_CJ(fun.clean.ss.CJ12)

resSig_MY.bac<-get_resSig_MY(bac.clean.ss.MY12)
resSig_MY.arch<-get_resSig_MY(arch.clean.ss.MY12)
resSig_MY.fun<-get_resSig_MY(fun.clean.ss.MY12)


resSig_CJ.arch$Enriched <- ifelse(resSig_CJ.arch$adj.P.Val < .05 & resSig_CJ.arch$logFC >0, "Non-fer", ifelse(resSig_CJ.arch$adj.P.Val < .05 & resSig_CJ.arch$logFC <0, "Fer", "ns"))
resSig_CJ.bac$Enriched <- ifelse(resSig_CJ.bac$adj.P.Val < .05 & resSig_CJ.bac$logFC >0, "Non-fer", ifelse(resSig_CJ.bac$adj.P.Val < .05 & resSig_CJ.bac$logFC <0, "Fer", "ns"))
resSig_CJ.fun$Enriched <- ifelse(resSig_CJ.fun$adj.P.Val < .05 & resSig_CJ.fun$logFC >0, "Non-fer", ifelse(resSig_CJ.fun$adj.P.Val < .05 & resSig_CJ.fun$logFC <0, "Fer", "ns"))


resSig_MY.arch$Enriched <- ifelse(resSig_MY.arch$adj.P.Val < .05 & resSig_MY.arch$logFC >0, "Non-fer", ifelse(resSig_MY.arch$adj.P.Val < .05 & resSig_MY.arch$logFC <0, "Fer", "ns"))
resSig_MY.fun$Enriched <- ifelse(resSig_MY.fun$adj.P.Val < .05 & resSig_MY.fun$logFC >0, "Non-fer", ifelse(resSig_MY.fun$adj.P.Val < .05 & resSig_MY.fun$logFC <0, "Fer", "ns"))
resSig_MY.bac$Enriched <- ifelse(resSig_MY.bac$adj.P.Val < .05 & resSig_MY.bac$logFC >0, "Non-fer", ifelse(resSig_MY.bac$adj.P.Val < .05 & resSig_MY.bac$logFC <0, "Fer", "ns"))


resSig_CJ.arch$Enriched = as.factor(resSig_CJ.arch$Enriched)
resSig_CJ.bac$Enriched = as.factor(resSig_CJ.bac$Enriched)
resSig_CJ.fun$Enriched = as.factor(resSig_CJ.fun$Enriched)
resSig_MY.arch$Enriched = as.factor(resSig_MY.arch$Enriched)
resSig_MY.bac$Enriched = as.factor(resSig_MY.bac$Enriched)
resSig_MY.fun$Enriched = as.factor(resSig_MY.fun$Enriched)

resSig_CJ.arch$Enriched <- factor(resSig_CJ.arch$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_CJ.bac$Enriched <- factor(resSig_CJ.bac$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_CJ.fun$Enriched <- factor(resSig_CJ.fun$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_MY.arch$Enriched <- factor(resSig_MY.arch$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_MY.bac$Enriched <- factor(resSig_MY.bac$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_MY.fun$Enriched <- factor(resSig_MY.fun$Enriched, levels = c("Non-fer", "Fer", "ns"))

## Saving files
write.xlsx(resSig_CJ.arch, "Significant daOTU_CJ_Arch.xlsx")
write.xlsx(resSig_CJ.bac, "Significant daOTU_CJ_Bac.xlsx")
write.xlsx(resSig_CJ.fun, "Significant daOTU_CJ_Fun.xlsx")

write.xlsx(resSig_MY.arch, "Significant daOTU_MY_Arch.xlsx")
write.xlsx(resSig_MY.bac, "Significant daOTU_MY_Bac.xlsx")
write.xlsx(resSig_MY.fun, "Significant daOTU_MY_Fun.xlsx")

ggplot(resSig_MY.fun, aes(x=log10(AveExpr), y=logFC, color=Enriched)) +
  geom_point(size = 3, alpha = 0.7) +
  #geom_hline(color = "red3", yintercept = 0) +
  # stat_smooth(se = FALSE, method = "loess", color = "red3") +
  scale_color_manual(values=c("Black","Red", "gray50"))+
  xlab('\n Average Abundance (Log10)')+
  ylab("Fold Change (Log2)\n") +
  #ggtitle("MA Plot \n") +
  theme(plot.title = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(aspect.ratio = 1)+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())




###Indicator OTU analysis
bac.clean.nolog.CJ <- subset_samples(bac.clean.nolog, Field %in% c("CJ1", "CJ2"))
bac.clean.nolog.CJ <- phyloseq::filter_taxa(bac.clean.nolog.CJ, function(x) sum(x) != 0, TRUE)

arch.clean.nolog.CJ <- subset_samples(arch.clean.nolog, Field %in% c("CJ1", "CJ2"))
arch.clean.nolog.CJ <- phyloseq::filter_taxa(arch.clean.nolog.CJ, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.CJ <- subset_samples(fun.clean.nolog, Field %in% c("CJ1", "CJ2"))
fun.clean.nolog.CJ <- phyloseq::filter_taxa(fun.clean.nolog.CJ, function(x) sum(x) != 0, TRUE)

bac.clean.nolog.MY <- subset_samples(bac.clean.nolog, Field %in% c("MY1", "MY2"))
bac.clean.nolog.MY <- phyloseq::filter_taxa(bac.clean.nolog.MY, function(x) sum(x) != 0, TRUE)

arch.clean.nolog.MY <- subset_samples(arch.clean.nolog, Field %in% c("MY1", "MY2"))
arch.clean.nolog.MY <- phyloseq::filter_taxa(arch.clean.nolog.MY, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.MY <- subset_samples(fun.clean.nolog, Field %in% c("MY1", "MY2"))
fun.clean.nolog.MY <- phyloseq::filter_taxa(fun.clean.nolog.MY, function(x) sum(x) != 0, TRUE)


##Oligotrophic environment
otu_norm_soil_16s <- otu_table(bac.clean.nolog.CJ)
design_16s_soil <- meta(bac.clean.nolog.CJ)

indic_soil_16s <- as.data.frame(t(otu_norm_soil_16s))
indic_soil_groups_16s <- design_16s_soil$Field
length(unique(indic_soil_groups_16s))

## Define indicator species for soil bacteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_16s <- multipatt(indic_soil_16s,indic_soil_groups_16s,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_16s,alpha=1,indvalcomp=T)
indic_soil_df_16s <- indicatorsp_soil_16s$sign
write.table(indic_soil_df_16s,"indicsp_bacteria_CJ.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_soil_df_16s <- read.table("indicsp_soil_16s_2018_nolog.txt", header=T, sep="\t")

Non_fer_soil_16s <- as.matrix(indic_soil_df_16s[which(indic_soil_df_16s$s.CJ1 == 1 & indic_soil_df_16s$p.value < 0.05),])
fer_soil_16s <- as.matrix(indic_soil_df_16s[which(indic_soil_df_16s$s.CJ2 == 1 & indic_soil_df_16s$p.value < 0.05),])

soil_r_values_16s_CJ <- rbind(Non_fer_soil_16s,fer_soil_16s)
colnames(soil_r_values_16s_CJ)[1:2] <-c("Non-fer","Fer")



##Fungi
otu_norm_soil_its <- otu_table(fun.clean.nolog.CJ)
design_its_soil <- meta(fun.clean.nolog.CJ)

indic_soil_its <- as.data.frame(t(otu_norm_soil_its))
indic_soil_groups_its <- design_its_soil$Field
length(unique(indic_soil_groups_its))

## Define indicator species for soil funteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_its <- multipatt(indic_soil_its,indic_soil_groups_its,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_its,alpha=1,indvalcomp=T)
indic_soil_df_its <- indicatorsp_soil_its$sign
write.table(indic_soil_df_its,"indicsp_funteria_CJ.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time
Non_fer_soil_its <- as.matrix(indic_soil_df_its[which(indic_soil_df_its$s.CJ1 == 1 & indic_soil_df_its$p.value < 0.05),])
fer_soil_its <- as.matrix(indic_soil_df_its[which(indic_soil_df_its$s.CJ2 == 1 & indic_soil_df_its$p.value < 0.05),])

soil_r_values_its_CJ <- rbind(Non_fer_soil_its,fer_soil_its)
colnames(soil_r_values_its_CJ)[1:2] <-c("Non-fer","Fer")


##Archaea
otu_norm_soil_arch <- otu_table(arch.clean.nolog.CJ)
design_arch_soil <- meta(arch.clean.nolog.CJ)

indic_soil_arch <- as.data.frame(t(otu_norm_soil_arch))
indic_soil_groups_arch <- design_arch_soil$Field
length(unique(indic_soil_groups_arch))

## Define indicator species for soil funteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_arch <- multipatt(indic_soil_arch,indic_soil_groups_arch,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_arch,alpha=1,indvalcomp=T)
indic_soil_df_arch <- indicatorsp_soil_arch$sign
write.table(indic_soil_df_arch,"indicsp_archaea_CJ.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time
Non_fer_soil_arch <- as.matrix(indic_soil_df_arch[which(indic_soil_df_arch$s.CJ1 == 1 & indic_soil_df_arch$p.value < 0.05),])
fer_soil_arch <- as.matrix(indic_soil_df_arch[which(indic_soil_df_arch$s.CJ2 == 1 & indic_soil_df_arch$p.value < 0.05),])

soil_r_values_arch_CJ <- rbind(Non_fer_soil_arch,fer_soil_arch)
colnames(soil_r_values_arch_CJ)[1:2] <-c("Non-fer","Fer")








## Range of correlation coefficients
range(soil_r_values_16s_CJ[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(soil_r_values_16s_CJ)))

## Proportion of soil bacteria OTUs responding to cropping system
length(unique(rownames(soil_r_values_16s_CJ)))/nrow(otu_norm_soil_16s_CJ)

## Proportion of soil bacteria sequences responding to cropping system
otu_16s_soil <- otu_table(phy.clean.ss.2018.f)
tax_soil_16s <- tax_table(phy.clean.ss.2018.f)
soil_16s_ra <- t(t(otu_16s_soil)/colSums(otu_16s_soil)) * 100
sum(colSums(soil_16s_ra[unique(rownames(soil_r_values_16s)),]))/sum(colSums(soil_16s_ra))


## Copiotrophic environment
otu_norm_soil_16s.MY <- otu_table(bac.clean.nolog.MY)
design_16s_soil.MY <- meta(bac.clean.nolog.MY)

indic_soil_16s.MY <- as.data.frame(t(otu_norm_soil_16s.MY))
indic_soil_groups_16s.MY <- design_16s_soil.MY$Field
length(unique(indic_soil_groups_16s.MY))

## Define indicator species for soil bacteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_16s.MY <- multipatt(indic_soil_16s.MY,indic_soil_groups_16s.MY,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_16s.MY,alpha=1,indvalcomp=T)
indic_soil_df_16s.MY <- indicatorsp_soil_16s.MY$sign
write.table(indic_soil_df_16s.MY,"indicsp_bacteria_MY.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time

fer_soil_16s.MY <- as.matrix(indic_soil_df_16s.MY[which(indic_soil_df_16s.MY$s.MY1 == 1 & indic_soil_df_16s.MY$p.value < 0.05),])
Non_fer_soil_16s.MY <- as.matrix(indic_soil_df_16s.MY[which(indic_soil_df_16s.MY$s.MY2 == 1 & indic_soil_df_16s.MY$p.value < 0.05),])

soil_r_values_16s_MY <- rbind(Non_fer_soil_16s.MY,fer_soil_16s.MY)
colnames(soil_r_values_16s_MY)[1:2] <-c("fer","Non_fer")



##Fungi
otu_norm_soil_its.MY <- otu_table(fun.clean.nolog.MY)
design_its_soil.MY <- meta(fun.clean.nolog.MY)

indic_soil_its.MY <- as.data.frame(t(otu_norm_soil_its.MY))
indic_soil_groups_its.MY <- design_its_soil.MY$Field
length(unique(indic_soil_groups_its.MY))

## Define indicator species for soil funteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_its.MY <- multipatt(indic_soil_its.MY,indic_soil_groups_its.MY,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_its.MY,alpha=1,indvalcomp=T)
indic_soil_df_its.MY <- indicatorsp_soil_its.MY$sign
write.table(indic_soil_df_its.MY,"indicsp_fungi_MY.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time

fer_soil_its.MY <- as.matrix(indic_soil_df_its.MY[which(indic_soil_df_its.MY$s.MY1 == 1 & indic_soil_df_its.MY$p.value < 0.05),])
Non_fer_soil_its.MY <- as.matrix(indic_soil_df_its.MY[which(indic_soil_df_its.MY$s.MY2 == 1 & indic_soil_df_its.MY$p.value < 0.05),])

soil_r_values_its_MY <- rbind(Non_fer_soil_its.MY,fer_soil_its.MY)
colnames(soil_r_values_its_MY)[1:2] <-c("fer","Non_fer")

##Archaea
otu_norm_soil_arch.MY <- otu_table(arch.clean.nolog.MY)
design_arch_soil.MY <- meta(arch.clean.nolog.MY)

indic_soil_arch.MY <- as.data.frame(t(otu_norm_soil_arch.MY))
indic_soil_groups_arch.MY <- design_arch_soil.MY$Field
length(unique(indic_soil_groups_arch.MY))

## Define indicator species for soil archteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_arch.MY <- multipatt(indic_soil_arch.MY,indic_soil_groups_arch.MY,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_arch.MY,alpha=1,indvalcomp=T)
indic_soil_df_arch.MY <- indicatorsp_soil_arch.MY$sign
write.table(indic_soil_df_arch.MY,"indicsp_archteria_MY.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time

fer_soil_arch.MY <- as.matrix(indic_soil_df_arch.MY[which(indic_soil_df_arch.MY$s.MY1 == 1 & indic_soil_df_arch.MY$p.value < 0.05),])
Non_fer_soil_arch.MY <- as.matrix(indic_soil_df_arch.MY[which(indic_soil_df_arch.MY$s.MY2 == 1 & indic_soil_df_arch.MY$p.value < 0.05),])

soil_r_values_arch_MY <- rbind(Non_fer_soil_arch.MY,fer_soil_arch.MY)
colnames(soil_r_values_arch_MY)[1:2] <-c("fer","Non_fer")


### Fertilization-sensitive OTU
head(resSig_CJ.bac)


resSig_CJ.bac.sub <- subset(resSig_CJ.bac, adj.P.Val < 0.05)
fs_otu.bac_CJ <- intersect(resSig_CJ.bac.sub$OTU, rownames(soil_r_values_16s_CJ))
fs_otu.bac_CJ.otu.id <- subset(OTU_id.list, OTU%in%resSig_CJ.bac.sub$OTU)

resSig_CJ.fun.sub <- subset(resSig_CJ.fun, adj.P.Val < 0.05)
fs_otu.fun_CJ <- intersect(resSig_CJ.fun.sub$OTU, rownames(soil_r_values_its_CJ))
fs_otu.fun_CJ.otu.id <- subset(OTU_id.list, OTU%in%fs_otu.fun_CJ)

resSig_CJ.arch.sub <- subset(resSig_CJ.arch, adj.P.Val < 0.05)
fs_otu.arch_CJ <- intersect(resSig_CJ.arch.sub$OTU, rownames(soil_r_values_arch_CJ))
fs_otu.arch_CJ.otu.id <- subset(OTU_id.list, OTU%in%fs_otu.arch_CJ)

resSig_MY.bac.sub <- subset(resSig_MY.bac, adj.P.Val < 0.05)
fs_otu.bac_MY <- intersect(resSig_MY.bac.sub$OTU, rownames(soil_r_values_16s_MY))
fs_otu.bac_MY.otu.id <- subset(OTU_id.list, OTU%in%fs_otu.bac_MY)

resSig_MY.fun.sub <- subset(resSig_MY.fun, adj.P.Val < 0.05)
fs_otu.fun_MY <- intersect(resSig_MY.fun.sub$OTU, rownames(soil_r_values_its_MY))
fs_otu.fun_MY.otu.id <- subset(OTU_id.list, OTU%in%fs_otu.fun_MY)

resSig_MY.arch.sub <- subset(resSig_MY.arch, adj.P.Val < 0.05)
fs_otu.arch_MY <- intersect(resSig_MY.arch.sub$OTU, rownames(soil_r_values_arch_MY))
fs_otu.arch_MY.otu.id <- subset(OTU_id.list, OTU%in%fs_otu.arch_MY)



## Oligo and copio in fsOTU
intersect(fs_otu.bac_CJ.otu.id$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Oligo"])
intersect(fs_otu.bac_CJ.otu.id$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Copio"])

intersect(fs_otu.bac_CJ.otu.id$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Oligo"])
intersect(fs_otu.bac_CJ.otu.id$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Copio"])


intersect(fs_otu.fun_CJ.otu.id$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Oligo"])
intersect(fs_otu.fun_CJ.otu.id$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Copio"])

intersect(fs_otu.fun_CJ.otu.id$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Oligo"])
intersect(fs_otu.fun_CJ.otu.id$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Copio"])

##Oligo and copio in daOTUs
resSig_CJ.bac.sub.nonfer <- subset(resSig_CJ.bac.sub, Enriched == "Non-fer")
resSig_CJ.bac.sub.fer <- subset(resSig_CJ.bac.sub, Enriched == "Fer")

daotu.bac_CJ.nonfer<- subset(OTU_id.list, OTU%in%resSig_CJ.bac.sub.nonfer$OTU)
daotu.bac_CJ.fer<- subset(OTU_id.list, OTU%in%resSig_CJ.bac.sub.fer$OTU)

intersect(daotu.bac_CJ.nonfer$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Oligo"]) #169
intersect(daotu.bac_CJ.nonfer$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Copio"]) #57


intersect(daotu.bac_CJ.fer$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Oligo"]) #58
intersect(daotu.bac_CJ.fer$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Copio"]) #29


resSig_CJ.fun.sub.nonfer <- subset(resSig_CJ.fun.sub, Enriched == "Non-fer")
resSig_CJ.fun.sub.fer <- subset(resSig_CJ.fun.sub, Enriched == "Fer")

daotu.fun_CJ.nonfer<- subset(OTU_id.list, OTU%in%resSig_CJ.fun.sub.nonfer$OTU)
daotu.fun_CJ.fer<- subset(OTU_id.list, OTU%in%resSig_CJ.fun.sub.fer$OTU)

intersect(daotu.fun_CJ.nonfer$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Oligo"]) #41
intersect(daotu.fun_CJ.nonfer$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Copio"]) #29


intersect(daotu.fun_CJ.fer$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Oligo"]) #1
intersect(daotu.fun_CJ.fer$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Copio"]) #9


resSig_CJ.arch.sub.nonfer <- subset(resSig_CJ.arch.sub, Enriched == "Non-fer")
resSig_CJ.arch.sub.fer <- subset(resSig_CJ.arch.sub, Enriched == "Fer")

daotu.arch_CJ.nonfer<- subset(OTU_id.list, OTU%in%resSig_CJ.arch.sub.nonfer$OTU)
daotu.arch_CJ.fer<- subset(OTU_id.list, OTU%in%resSig_CJ.arch.sub.fer$OTU)

intersect(daotu.arch_CJ.nonfer$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Oligo"]) #20
intersect(daotu.arch_CJ.nonfer$OTU_id, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Copio"]) #6


intersect(daotu.arch_CJ.fer$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Oligo"]) #1
intersect(daotu.arch_CJ.fer$OTU_id, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Copio"]) #0



resSig_MY.bac.sub.nonfer <- subset(resSig_MY.bac.sub, Enriched == "Non-fer")
resSig_MY.bac.sub.fer <- subset(resSig_MY.bac.sub, Enriched == "Fer")

daotu.bac_MY.nonfer<- subset(OTU_id.list, OTU%in%resSig_MY.bac.sub.nonfer$OTU)
daotu.bac_MY.fer<- subset(OTU_id.list, OTU%in%resSig_MY.bac.sub.fer$OTU)

intersect(daotu.bac_MY.nonfer$OTU_id, rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Oligo"]) #2
intersect(daotu.bac_MY.nonfer$OTU_id, rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Copio"]) #14


intersect(daotu.bac_MY.fer$OTU_id, rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Oligo"]) #1
intersect(daotu.bac_MY.fer$OTU_id, rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Copio"]) #17


resSig_MY.fun.sub.nonfer <- subset(resSig_MY.fun.sub, Enriched == "Non-fer")
resSig_MY.fun.sub.fer <- subset(resSig_MY.fun.sub, Enriched == "Fer")

daotu.fun_MY.nonfer<- subset(OTU_id.list, OTU%in%resSig_MY.fun.sub.nonfer$OTU)
daotu.fun_MY.fer<- subset(OTU_id.list, OTU%in%resSig_MY.fun.sub.fer$OTU)

intersect(daotu.fun_MY.nonfer$OTU_id, rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Oligo"]) #0
intersect(daotu.fun_MY.nonfer$OTU_id, rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Copio"]) #5


intersect(daotu.fun_MY.fer$OTU_id, rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Oligo"]) #1
intersect(daotu.fun_MY.fer$OTU_id, rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Copio"]) #19


resSig_MY.arch.sub.nonfer <- subset(resSig_MY.arch.sub, Enriched == "Non-fer")
resSig_MY.arch.sub.fer <- subset(resSig_MY.arch.sub, Enriched == "Fer")

daotu.arch_MY.nonfer<- subset(OTU_id.list, OTU%in%resSig_MY.arch.sub.nonfer$OTU)
daotu.arch_MY.fer<- subset(OTU_id.list, OTU%in%resSig_MY.arch.sub.fer$OTU)

intersect(daotu.arch_MY.nonfer$OTU_id, rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Oligo"]) #0
intersect(daotu.arch_MY.nonfer$OTU_id, rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Copio"]) #3


intersect(daotu.arch_MY.fer$OTU_id, rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Oligo"]) #0
intersect(daotu.arch_MY.fer$OTU_id, rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Copio"]) #0




### CJ+YS
get_resSig <- function(phy.clean.ss.5){
  phy.clean.filt <- phyloseq::filter_taxa(phy.clean.ss.5, function(x) sum(x) != 0, TRUE)
  sum(taxa_sums(phy.clean.filt) == 0)
  
  obj <- phyloseq_to_metagenomeSeq(phy.clean.filt)
  obj <-  cumNorm(obj, p = cumNormStatFast(obj))
  normFactor <-  normFactors(obj)
  normFactor <-  log2(normFactor/median(normFactor) + 1)
  settings <-  zigControl(maxit = 10, verbose = TRUE)
  Type  <-  pData(obj)$Cultural_practice
  mod  <-  model.matrix(~Type)
  colnames(mod)  <-  levels(Type)
  colnames(mod)
  
  res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)
  res
  res$fit
  zigFit = res$fit
  finalMod= res$fit$design
  finalMod
  contrast.matrix =makeContrasts(No_fertilizer - No_pesticide, levels = finalMod)
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit2
  fit3 = eBayes(fit2)
  fit3
  topTable(fit3, coef="No_fertilizer - No_pesticide")
  
  
  res <- topTable(fit3,coef=1,adjust="fdr",n=Inf)
  
  log2AverageAbundance <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
  log2AverageAbundance
  Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
  Ta <- unique(Ta)
  Ta
  resSig = res[!is.na(res$adj.P.Val), ]
  resSig = data.frame(resSig)
  head(resSig)
  resSig <- tibble::rownames_to_column(resSig, 'OTU')
  resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
  resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))
  return(resSig)
}

resSig_bac.CJ.YS<-get_resSig(bac.clean.ss.CJ.YS)
resSig_arch.CJ.YS<-get_resSig(arch.clean.ss.CJ.YS)
resSig_fun.CJ.YS<-get_resSig(fun.clean.ss.CJ.YS)

resSig_CJ.YS.arch <-resSig_arch.CJ.YS
resSig_CJ.YS.fun <-resSig_fun.CJ.YS


resSig_CJ.YS.arch$Enriched <- ifelse(resSig_CJ.YS.arch$adj.P.Val < .05 & resSig_CJ.YS.arch$logFC >0, "Non-fer", ifelse(resSig_CJ.YS.arch$adj.P.Val < .05 & resSig_CJ.YS.arch$logFC <0, "Fer", "ns"))
resSig_bac.CJ.YS$Enriched <- ifelse(resSig_bac.CJ.YS$adj.P.Val < .05 & resSig_bac.CJ.YS$logFC >0, "Non-fer", ifelse(resSig_bac.CJ.YS$adj.P.Val < .05 & resSig_bac.CJ.YS$logFC <0, "Fer", "ns"))
resSig_CJ.YS.fun$Enriched <- ifelse(resSig_CJ.YS.fun$adj.P.Val < .05 & resSig_CJ.YS.fun$logFC >0, "Non-fer", ifelse(resSig_CJ.YS.fun$adj.P.Val < .05 & resSig_CJ.YS.fun$logFC <0, "Fer", "ns"))


resSig_CJ.YS.arch$Enriched = as.factor(resSig_CJ.YS.arch$Enriched)
resSig_bac.CJ.YS$Enriched = as.factor(resSig_bac.CJ.YS$Enriched)
resSig_CJ.YS.fun$Enriched = as.factor(resSig_CJ.YS.fun$Enriched)

resSig_CJ.YS.arch$Enriched <- factor(resSig_CJ.YS.arch$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_bac.CJ.YS$Enriched <- factor(resSig_bac.CJ.YS$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_CJ.YS.fun$Enriched <- factor(resSig_CJ.YS.fun$Enriched, levels = c("Non-fer", "Fer", "ns"))


## MY + NJ
bac.clean.ss.MY.NJ <- subset_samples(bac.clean.ss, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
bac.clean.ss.MY.NJ <- phyloseq::filter_taxa(bac.clean.ss.MY.NJ, function(x) sum(x) != 0, TRUE)

arch.clean.ss.MY.NJ <- subset_samples(arch.clean.ss, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
arch.clean.ss.MY.NJ <- phyloseq::filter_taxa(arch.clean.ss.MY.NJ, function(x) sum(x) != 0, TRUE)

fun.clean.ss.MY.NJ <- subset_samples(fun.clean.ss, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
fun.clean.ss.MY.NJ <- phyloseq::filter_taxa(fun.clean.ss.MY.NJ, function(x) sum(x) != 0, TRUE)




get_resSig_copio <- function(phy.clean.ss.5){
  phy.clean.filt <- phyloseq::filter_taxa(phy.clean.ss.5, function(x) sum(x) != 0, TRUE)
  sum(taxa_sums(phy.clean.filt) == 0)
  
  obj <- phyloseq_to_metagenomeSeq(phy.clean.filt)
  obj <-  cumNorm(obj, p = cumNormStatFast(obj))
  normFactor <-  normFactors(obj)
  normFactor <-  log2(normFactor/median(normFactor) + 1)
  settings <-  zigControl(maxit = 10, verbose = TRUE)
  Type  <-  pData(obj)$Cultural_practice
  mod  <-  model.matrix(~Type)
  colnames(mod)  <-  levels(Type)
  colnames(mod)
  
  res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)
  res
  res$fit
  zigFit = res$fit
  finalMod= res$fit$design
  finalMod
  contrast.matrix =makeContrasts(No_fertilizer -Conventional, levels = finalMod)
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit2
  fit3 = eBayes(fit2)
  fit3
  topTable(fit3, coef="No_fertilizer - Conventional")
  
  
  res <- topTable(fit3,coef=1,adjust="fdr",n=Inf)
  
  log2AverageAbundance <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
  log2AverageAbundance
  Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
  Ta <- unique(Ta)
  Ta
  resSig = res[!is.na(res$adj.P.Val), ]
  resSig = data.frame(resSig)
  head(resSig)
  resSig <- tibble::rownames_to_column(resSig, 'OTU')
  resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
  resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))
  return(resSig)
}

resSig_bac.MY.NJ<-get_resSig_copio(bac.clean.ss.MY.NJ)
resSig_arch.MY.NJ<-get_resSig_copio(arch.clean.ss.MY.NJ)
resSig_fun.MY.NJ<-get_resSig_copio(fun.clean.ss.MY.NJ)

resSig_MY.NJ.arch <-resSig_arch.MY.NJ
resSig_MY.NJ.fun <-resSig_fun.MY.NJ


resSig_MY.NJ.arch$Enriched <- ifelse(resSig_MY.NJ.arch$adj.P.Val < .05 & resSig_MY.NJ.arch$logFC >0, "Non-fer", ifelse(resSig_MY.NJ.arch$adj.P.Val < .05 & resSig_MY.NJ.arch$logFC <0, "Fer", "ns"))
resSig_bac.MY.NJ$Enriched <- ifelse(resSig_bac.MY.NJ$adj.P.Val < .05 & resSig_bac.MY.NJ$logFC >0, "Non-fer", ifelse(resSig_bac.MY.NJ$adj.P.Val < .05 & resSig_bac.MY.NJ$logFC <0, "Fer", "ns"))
resSig_MY.NJ.fun$Enriched <- ifelse(resSig_MY.NJ.fun$adj.P.Val < .05 & resSig_MY.NJ.fun$logFC >0, "Non-fer", ifelse(resSig_MY.NJ.fun$adj.P.Val < .05 & resSig_MY.NJ.fun$logFC <0, "Fer", "ns"))


resSig_MY.NJ.arch$Enriched = as.factor(resSig_MY.NJ.arch$Enriched)
resSig_bac.MY.NJ$Enriched = as.factor(resSig_bac.MY.NJ$Enriched)
resSig_MY.NJ.fun$Enriched = as.factor(resSig_MY.NJ.fun$Enriched)

resSig_MY.NJ.arch$Enriched <- factor(resSig_MY.NJ.arch$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_bac.MY.NJ$Enriched <- factor(resSig_bac.MY.NJ$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_MY.NJ.fun$Enriched <- factor(resSig_MY.NJ.fun$Enriched, levels = c("Non-fer", "Fer", "ns"))



### MA plot for oligo and copiotrophic fields
resSig_CJ.YS.bac<- read.table("Bacteria_daOTUs_Oligo.tsv",sep = '\t', header = T)
resSig_CJ.YS.arch<- read.table("Archaea_daOTUs_Oligo.tsv",sep = '\t', header = T)
resSig_CJ.YS.fun<- read.table("Fungi_daOTUs_Oligo.tsv",sep = '\t', header = T)

resSig_MY.NJ.bac<- read.table("Bacteria_daOTUs_Copio.tsv",sep = '\t', header = T)
resSig_MY.NJ.arch<- read.table("Archaea_daOTUs_Copio.tsv",sep = '\t', header = T)
resSig_MY.NJ.fun<- read.table("Fungi_daOTUs_Copio.tsv",sep = '\t', header = T)

resSig_CJ.YS.arch$Enriched <- factor(resSig_CJ.YS.arch$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_CJ.YS.bac$Enriched <- factor(resSig_CJ.YS.bac$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_CJ.YS.fun$Enriched <- factor(resSig_CJ.YS.fun$Enriched, levels = c("Non-fer", "Fer", "ns"))

resSig_MY.NJ.arch$Enriched <- factor(resSig_MY.NJ.arch$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_MY.NJ.bac$Enriched <- factor(resSig_MY.NJ.bac$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_MY.NJ.fun$Enriched <- factor(resSig_MY.NJ.fun$Enriched, levels = c("Non-fer", "Fer", "ns"))

get_MAplot<-function(resSig_OTU){
  p<-ggplot(resSig_OTU, aes(x=log10(AveExpr), y=logFC, color=Enriched)) +
  geom_point(size = 3, alpha = 0.7) +
  #geom_hline(color = "red3", yintercept = 0) +
  # stat_smooth(se = FALSE, method = "loess", color = "red3") +
  scale_color_manual(values=c("Black","#6699cc", "#cccccc"))+
  xlab('\n Log 10 Average Abundance')+
  ylab("Log2 Fold Change\n") +
  #ggtitle("MA Plot \n") +
  theme(plot.title = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(aspect.ratio = 1)+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
return(p)}


get_MAplot(resSig_CJ.YS.bac)
get_MAplot(resSig_CJ.YS.arch)
get_MAplot(resSig_CJ.YS.fun)

get_MAplot(resSig_MY.NJ.bac)
get_MAplot(resSig_MY.NJ.arch)
get_MAplot(resSig_MY.NJ.fun)


### Assiging putative oligo and copio
tab_ra.a.SOM.TN.oligo$OTU <- rownames(tab_ra.a.SOM.TN.oligo)
tab_ra.b.SOM.TN.oligo$OTU <- rownames(tab_ra.b.SOM.TN.oligo)
tab_ra.f.SOM.TN.oligo$OTU <- rownames(tab_ra.f.SOM.TN.oligo)

tab_ra.a.SOM.TN.copi$OTU <- rownames(tab_ra.a.SOM.TN.copi)
tab_ra.b.SOM.TN.copi$OTU <- rownames(tab_ra.b.SOM.TN.copi)
tab_ra.f.SOM.TN.copi$OTU <- rownames(tab_ra.f.SOM.TN.copi)


arch_oligo<-merge(tab_ra.a.SOM.TN.oligo, OTU.list, by = 'OTU')
bac_oligo<-merge(tab_ra.b.SOM.TN.oligo, OTU.list, by = 'OTU')
fun_oligo<-merge(tab_ra.f.SOM.TN.oligo, OTU.list, by = 'OTU')


arch_copio<-merge(tab_ra.a.SOM.TN.copi, OTU.list, by = 'OTU')
bac_copio<-merge(tab_ra.b.SOM.TN.copi, OTU.list, by = 'OTU')
fun_copio<-merge(tab_ra.f.SOM.TN.copi, OTU.list, by = 'OTU')

arch_oligo$category <- "Oligo"
bac_oligo$category <- "Oligo"
fun_oligo$category <- "Oligo"

arch_copio$category <- "Copio"
bac_copio$category <- "Copio"
fun_copio$category <- "Copio"


bac_oligo.copio <- rbind(bac_oligo, bac_copio)
arch_oligo.copio <- rbind(arch_oligo, arch_copio)
fun_oligo.copio <- rbind(fun_oligo, fun_copio)

write.xlsx(bac_oligo.copio, "Putative oligo and copiotrophs_Bacteria.xlsx")
write.xlsx(arch_oligo.copio, "Putative oligo and copiotrophs_Archaea.xlsx")
write.xlsx(fun_oligo.copio, "Putative oligo and copiotrophs_Fungi.xlsx")

head(resSig_CJ.YS.bac)


tab_ra.a.SOM.TN.oligo$category <- "Oligo"
tab_ra.a.SOM.TN.copi$category <- "Copio"
arch.oligo.copio.list<-rbind(tab_ra.a.SOM.TN.oligo, tab_ra.a.SOM.TN.copi)
head(arch.oligo.copio.list)

tab_ra.b.SOM.TN.oligo$category <- "Oligo"
tab_ra.b.SOM.TN.copi$category <- "Copio"
bac.oligo.copio.list<-rbind(tab_ra.b.SOM.TN.oligo, tab_ra.b.SOM.TN.copi)
head(bac.oligo.copio.list)


tab_ra.f.SOM.TN.oligo$category <- "Oligo"
tab_ra.f.SOM.TN.copi$category <- "Copio"
fun.oligo.copio.list<-rbind(tab_ra.f.SOM.TN.oligo, tab_ra.f.SOM.TN.copi)
head(fun.oligo.copio.list)

bac.daOTUs.oligo <-merge(resSig_CJ.YS.bac, bac.oligo.copio.list, by = 'OTU', all =T)
arch.daOTUs.oligo <-merge(resSig_CJ.YS.arch, arch.oligo.copio.list, by = 'OTU', all =T)
fun.daOTUs.oligo <-merge(resSig_CJ.YS.fun, fun.oligo.copio.list, by = 'OTU', all =T)


bac.daOTUs.copio <-merge(resSig_MY.NJ.bac, bac.oligo.copio.list, by = 'OTU', all =T)
arch.daOTUs.copio <-merge(resSig_MY.NJ.arch, arch.oligo.copio.list, by = 'OTU', all =T)
fun.daOTUs.copio <-merge(resSig_MY.NJ.fun, fun.oligo.copio.list, by = 'OTU', all =T)


bac.daOTUs.oligo$category[is.na(bac.daOTUs.oligo$category)] <- "Not_predicted"
arch.daOTUs.oligo$category[is.na(arch.daOTUs.oligo$category)] <- "Not_predicted"
fun.daOTUs.oligo$category[is.na(fun.daOTUs.oligo$category)] <- "Not_predicted"


bac.daOTUs.copio$category[is.na(bac.daOTUs.copio$category)] <- "Not_predicted"
arch.daOTUs.copio$category[is.na(arch.daOTUs.copio$category)] <- "Not_predicted"
fun.daOTUs.copio$category[is.na(fun.daOTUs.copio$category)] <- "Not_predicted"

bac.daOTUs.oligo$Enriched[is.na(bac.daOTUs.oligo$Enriched)] <- "ns"
arch.daOTUs.oligo$Enriched[is.na(arch.daOTUs.oligo$Enriched)] <- "ns"
fun.daOTUs.oligo$Enriched[is.na(fun.daOTUs.oligo$Enriched)] <- "ns"


bac.daOTUs.copio$Enriched[is.na(bac.daOTUs.copio$Enriched)] <- "ns"
arch.daOTUs.copio$Enriched[is.na(arch.daOTUs.copio$Enriched)] <- "ns"
fun.daOTUs.copio$Enriched[is.na(fun.daOTUs.copio$Enriched)] <- "ns"


bac.daOTUs.oligo <-bac.daOTUs.oligo.copy

## significant different and oligo or copio

get_sigdiff_oligo.copio<-function(df.daotu){for (i in as.character(df.daotu$OTU))
{
  if (i %in% intersect(df.daotu$OTU[(df.daotu$Enriched) == "Non-fer"], df.daotu$OTU[(df.daotu$category) == "Oligo"]) == TRUE)
  {df.daotu[df.daotu$OTU==i,"category2"] <- "Nonfer_oligo"}
  
  else if (i %in% intersect(df.daotu$OTU[(df.daotu$Enriched) == "Non-fer"], df.daotu$OTU[(df.daotu$category) == "Copio"]) == TRUE)
  {df.daotu[df.daotu$OTU==i,"category2"] <- "Nonfer_copio"} 
  
  else if (i %in% intersect(df.daotu$OTU[(df.daotu$Enriched) == "Fer"], df.daotu$OTU[(df.daotu$category) == "Copio"]) == TRUE)
  {df.daotu[df.daotu$OTU==i,"category2"] <- "Fer_copio"} 
  
  else if (i %in% intersect(df.daotu$OTU[(df.daotu$Enriched) == "Fer"], df.daotu$OTU[(df.daotu$category) == "Oligo"]) == TRUE)
  {df.daotu[df.daotu$OTU==i,"category2"] <- "Fer_oligo"} 
  
  else if (i %in% intersect(df.daotu$OTU[(df.daotu$Enriched) == "Non-fer"], df.daotu$OTU[(df.daotu$category) == "Not_predicted"]) == TRUE)
  {df.daotu[df.daotu$OTU==i,"category2"] <- "Nonfer_NP"} 
  
  else if (i %in% intersect(df.daotu$OTU[(df.daotu$Enriched) == "Fer"], df.daotu$OTU[(df.daotu$category) == "Not_predicted"]) == TRUE)
  {df.daotu[df.daotu$OTU==i,"category2"] <- "Fer_NP"} 
  
  else
  {df.daotu[df.daotu$OTU==i,"category2"] <- "ns"}
}
  return(df.daotu)
}

bac.daOTUs.oligo_sigdiff<-get_sigdiff_oligo.copio(bac.daOTUs.oligo)
arch.daOTUs.oligo_sigdiff<-get_sigdiff_oligo.copio(arch.daOTUs.oligo)
fun.daOTUs.oligo_sigdiff<-get_sigdiff_oligo.copio(fun.daOTUs.oligo)

bac.daOTUs.copio_sigdiff<-get_sigdiff_oligo.copio(bac.daOTUs.copio)
arch.daOTUs.copio_sigdiff<-get_sigdiff_oligo.copio(arch.daOTUs.copio)
fun.daOTUs.copio_sigdiff<-get_sigdiff_oligo.copio(fun.daOTUs.copio)



bac.daOTUs.oligo_sigdiff$category2 <- factor(bac.daOTUs.oligo_sigdiff$category2, levels = c('Nonfer_oligo', 'Nonfer_copio', "Nonfer_NP",'Fer_oligo','Fer_copio',  "Fer_NP",'ns'))
arch.daOTUs.oligo_sigdiff$category2 <- factor(arch.daOTUs.oligo_sigdiff$category2, levels = c('Nonfer_oligo', 'Nonfer_copio', "Nonfer_NP",'Fer_oligo','Fer_copio',  "Fer_NP",'ns'))
fun.daOTUs.oligo_sigdiff$category2 <- factor(fun.daOTUs.oligo_sigdiff$category2, levels = c('Nonfer_oligo', 'Nonfer_copio', "Nonfer_NP",'Fer_oligo','Fer_copio',  "Fer_NP",'ns'))

bac.daOTUs.copio_sigdiff$category2 <- factor(bac.daOTUs.copio_sigdiff$category2, levels = c('Nonfer_oligo', 'Nonfer_copio', "Nonfer_NP",'Fer_oligo','Fer_copio',  "Fer_NP",'ns'))
arch.daOTUs.copio_sigdiff$category2 <- factor(arch.daOTUs.copio_sigdiff$category2, levels = c('Nonfer_oligo', 'Nonfer_copio', "Nonfer_NP",'Fer_oligo','Fer_copio',  "Fer_NP",'ns'))
fun.daOTUs.copio_sigdiff$category2 <- factor(fun.daOTUs.copio_sigdiff$category2, levels = c('Nonfer_oligo', 'Nonfer_copio', "Nonfer_NP",'Fer_oligo','Fer_copio',  "Fer_NP",'ns'))

bac.daOTUs.oligo_sigdiff.t<- subset(bac.daOTUs.oligo_sigdiff, OTU %in% resSig_CJ.YS.bac$OTU)
arch.daOTUs.oligo_sigdiff.t<- subset(arch.daOTUs.oligo_sigdiff, OTU %in% resSig_CJ.YS.arch$OTU)
fun.daOTUs.oligo_sigdiff.t<- subset(fun.daOTUs.oligo_sigdiff, OTU %in% resSig_CJ.YS.fun$OTU)

bac.daOTUs.copio_sigdiff.t<- subset(bac.daOTUs.copio_sigdiff, OTU %in% resSig_MY.NJ.bac$OTU)
arch.daOTUs.copio_sigdiff.t<- subset(arch.daOTUs.copio_sigdiff, OTU %in% resSig_MY.NJ.arch$OTU)
fun.daOTUs.copio_sigdiff.t<- subset(fun.daOTUs.copio_sigdiff, OTU %in% resSig_MY.NJ.fun$OTU)


get_MAplot_trophic<-function(resSig_OTU){
  p<-ggplot(resSig_OTU, aes(x=log10(AveExpr), y=logFC, color=category2)) +
    geom_point(size = 2, alpha = 1) +
    #geom_hline(color = "red3", yintercept = 0) +
    # stat_smooth(se = FALSE, method = "loess", color = "red3") +
    scale_color_manual(values=c('Nonfer_oligo'="#FFCC99",'Nonfer_copio' = "#996633", 'Fer_oligo' = "#CCCCFF",'Fer_copio'= "#6666CC", "Nonfer_NP" = "Black","Fer_NP" = "#6699cc","ns"= "#cccccc"))+
    xlab('\n Log 10 Average Abundance')+
    ylab("Log2 Fold Change\n") +
    #ggtitle("MA Plot \n") +
    theme(plot.title = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(aspect.ratio = 1)+
    theme(legend.position="right")+
    theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=12))+
    theme(axis.text.y = element_text(size=12))+
    
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
  return(p)}

get_MAplot_trophic(bac.daOTUs.oligo_sigdiff.t)
get_MAplot_trophic(arch.daOTUs.oligo_sigdiff.t)
get_MAplot_trophic(fun.daOTUs.oligo_sigdiff.t)


get_MAplot_trophic(bac.daOTUs.copio_sigdiff.t)
get_MAplot_trophic(arch.daOTUs.copio_sigdiff.t)
get_MAplot_trophic(fun.daOTUs.copio_sigdiff.t)


## number of OTUs in each group
length(bac.daOTUs.copio_sigdiff.t$category2[bac.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"]) #105
length(bac.daOTUs.copio_sigdiff.t$category2[bac.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"]) #144
length(bac.daOTUs.copio_sigdiff.t$category2[bac.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"]) #119
length(bac.daOTUs.copio_sigdiff.t$category2[bac.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"]) #117

length(bac.daOTUs.copio_sigdiff.t$Enriched[bac.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #1315
length(bac.daOTUs.copio_sigdiff.t$Enriched[bac.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #1230


length(arch.daOTUs.copio_sigdiff.t$category2[arch.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"]) #1
length(arch.daOTUs.copio_sigdiff.t$category2[arch.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"]) #8
length(arch.daOTUs.copio_sigdiff.t$category2[arch.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"]) #13
length(arch.daOTUs.copio_sigdiff.t$category2[arch.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"]) #16

length(arch.daOTUs.copio_sigdiff.t$Enriched[arch.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #85
length(arch.daOTUs.copio_sigdiff.t$Enriched[arch.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #120

length(fun.daOTUs.copio_sigdiff.t$category2[fun.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"]) #45
length(fun.daOTUs.copio_sigdiff.t$category2[fun.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"]) #70
length(fun.daOTUs.copio_sigdiff.t$category2[fun.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"]) #28
length(fun.daOTUs.copio_sigdiff.t$category2[fun.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"]) #123

length(fun.daOTUs.copio_sigdiff.t$Enriched[fun.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #571
length(fun.daOTUs.copio_sigdiff.t$Enriched[fun.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #335



## Oligotrophic

length(bac.daOTUs.oligo_sigdiff.t$category2[bac.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"]) #225
length(bac.daOTUs.oligo_sigdiff.t$category2[bac.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"]) #95
length(bac.daOTUs.oligo_sigdiff.t$category2[bac.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"]) #126
length(bac.daOTUs.oligo_sigdiff.t$category2[bac.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"]) #83

length(bac.daOTUs.oligo_sigdiff.t$Enriched[bac.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"]) #1246
length(bac.daOTUs.oligo_sigdiff.t$Enriched[bac.daOTUs.oligo_sigdiff.t$Enriched == "Fer"]) #1382


length(arch.daOTUs.oligo_sigdiff.t$category2[arch.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"]) #23
length(arch.daOTUs.oligo_sigdiff.t$category2[arch.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"]) #15
length(arch.daOTUs.oligo_sigdiff.t$category2[arch.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"]) #5
length(arch.daOTUs.oligo_sigdiff.t$category2[arch.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"]) #5

length(arch.daOTUs.oligo_sigdiff.t$Enriched[arch.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"]) #127
length(arch.daOTUs.oligo_sigdiff.t$Enriched[arch.daOTUs.oligo_sigdiff.t$Enriched == "Fer"]) #27

length(fun.daOTUs.oligo_sigdiff.t$category2[fun.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"]) #111
length(fun.daOTUs.oligo_sigdiff.t$category2[fun.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"]) #52
length(fun.daOTUs.oligo_sigdiff.t$category2[fun.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"]) #24
length(fun.daOTUs.oligo_sigdiff.t$category2[fun.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"]) #83

length(fun.daOTUs.oligo_sigdiff.t$Enriched[fun.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"]) #359
length(fun.daOTUs.oligo_sigdiff.t$Enriched[fun.daOTUs.oligo_sigdiff.t$Enriched == "Fer"]) #333


write.xlsx(bac.daOTUs.oligo_sigdiff.t, "Oligo_daOTUs with oligo and copiotrophs_Bacteria.xlsx")
write.xlsx(arch.daOTUs.oligo_sigdiff.t, "Oligo_daOTUs with oligo and copiotrophs_Archaea.xlsx")
write.xlsx(fun.daOTUs.oligo_sigdiff.t, "Oligo_daOTUs with oligo and copiotrophs_Fungi.xlsx")

write.xlsx(bac.daOTUs.copio_sigdiff.t, "Copio_daOTUs with oligo and copiotrophs_Bacteria.xlsx")
write.xlsx(arch.daOTUs.copio_sigdiff.t, "Copio_daOTUs with oligo and copiotrophs_Archaea.xlsx")
write.xlsx(fun.daOTUs.copio_sigdiff.t, "Copio_daOTUs with oligo and copiotrophs_Fungi.xlsx")


### Differences
length(bac.daOTUs.copio_sigdiff.t$category2[bac.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"]) #105
length(bac.daOTUs.copio_sigdiff.t$category2[bac.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"]) #144
length(bac.daOTUs.copio_sigdiff.t$category2[bac.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"]) #119
length(bac.daOTUs.copio_sigdiff.t$category2[bac.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"]) #117

length(bac.daOTUs.copio_sigdiff.t$Enriched[bac.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #1315
length(bac.daOTUs.copio_sigdiff.t$Enriched[bac.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #1230


length(arch.daOTUs.copio_sigdiff.t$category2[arch.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"]) #1
length(arch.daOTUs.copio_sigdiff.t$category2[arch.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"]) #8
length(arch.daOTUs.copio_sigdiff.t$category2[arch.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"]) #13
length(arch.daOTUs.copio_sigdiff.t$category2[arch.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"]) #16

length(arch.daOTUs.copio_sigdiff.t$Enriched[arch.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #85
length(arch.daOTUs.copio_sigdiff.t$Enriched[arch.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #120

length(fun.daOTUs.copio_sigdiff.t$category2[fun.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"]) #45
length(fun.daOTUs.copio_sigdiff.t$category2[fun.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"]) #70
length(fun.daOTUs.copio_sigdiff.t$category2[fun.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"]) #28
length(fun.daOTUs.copio_sigdiff.t$category2[fun.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"]) #123

length(fun.daOTUs.copio_sigdiff.t$Enriched[fun.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #571
length(fun.daOTUs.copio_sigdiff.t$Enriched[fun.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #335



## Oligotrophic
#Total number of OTUs
length(bac.daOTUs.oligo_sigdiff.t$OTU_id) #8045
length(arch.daOTUs.oligo_sigdiff.t$OTU_id) #342
length(fun.daOTUs.oligo_sigdiff.t$OTU_id) #2282

CJYS_nonfer_oligo.bac<-bac.daOTUs.oligo_sigdiff.t$OTU_id[bac.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"] #225
CJYS_nonfer_copio.bac<-bac.daOTUs.oligo_sigdiff.t$OTU_id[bac.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"] #95
CJYS_fer_oligo.bac<-bac.daOTUs.oligo_sigdiff.t$OTU_id[bac.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"] #126
CJYS_fer_copio.bac<-bac.daOTUs.oligo_sigdiff.t$OTU_id[bac.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"] #83

bac.daOTUs.oligo_sigdiff.t$Enriched[bac.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"] #1246
bac.daOTUs.oligo_sigdiff.t$Enriched[bac.daOTUs.oligo_sigdiff.t$Enriched == "Fer"] #1382


CJYS_nonfer_oligo.arch<-arch.daOTUs.oligo_sigdiff.t$OTU_id[arch.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"] #225
CJYS_nonfer_copio.arch<-arch.daOTUs.oligo_sigdiff.t$OTU_id[arch.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"] #95
CJYS_fer_oligo.arch<-arch.daOTUs.oligo_sigdiff.t$OTU_id[arch.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"] #126
CJYS_fer_copio.arch<-arch.daOTUs.oligo_sigdiff.t$OTU_id[arch.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"] #83


arch.daOTUs.oligo_sigdiff.t$Enriched[arch.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"] #127
arch.daOTUs.oligo_sigdiff.t$Enriched[arch.daOTUs.oligo_sigdiff.t$Enriched == "Fer"] #27

CJYS_nonfer_oligo.fun<-fun.daOTUs.oligo_sigdiff.t$OTU_id[fun.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"] #225
CJYS_nonfer_copio.fun<-fun.daOTUs.oligo_sigdiff.t$OTU_id[fun.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"] #95
CJYS_fer_oligo.fun<-fun.daOTUs.oligo_sigdiff.t$OTU_id[fun.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"] #126
CJYS_fer_copio.fun<-fun.daOTUs.oligo_sigdiff.t$OTU_id[fun.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"] #83


fun.daOTUs.oligo_sigdiff.t$Enriched[fun.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"] #359
fun.daOTUs.oligo_sigdiff.t$Enriched[fun.daOTUs.oligo_sigdiff.t$Enriched == "Fer"] #333


#Copiotrophic

length(bac.daOTUs.copio_sigdiff.t$OTU_id) #7969
length(arch.daOTUs.copio_sigdiff.t$OTU_id) #259
length(fun.daOTUs.copio_sigdiff.t$OTU_id) #2636

MYNJ_nonfer_oligo.bac<-bac.daOTUs.copio_sigdiff.t$OTU_id[bac.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"] #105
MYNJ_nonfer_copio.bac<-bac.daOTUs.copio_sigdiff.t$OTU_id[bac.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"] #144
MYNJ_fer_oligo.bac<-bac.daOTUs.copio_sigdiff.t$OTU_id[bac.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"] #119
MYNJ_fer_copio.bac<-bac.daOTUs.copio_sigdiff.t$OTU_id[bac.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"] #117

length(bac.daOTUs.copio_sigdiff.t$Enriched[bac.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #1315
length(bac.daOTUs.copio_sigdiff.t$Enriched[bac.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #1230


MYNJ_nonfer_oligo.arch<-arch.daOTUs.copio_sigdiff.t$OTU_id[arch.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"] #105
MYNJ_nonfer_copio.arch<-arch.daOTUs.copio_sigdiff.t$OTU_id[arch.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"] #144
MYNJ_fer_oligo.arch<-arch.daOTUs.copio_sigdiff.t$OTU_id[arch.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"] #119
MYNJ_fer_copio.arch<-arch.daOTUs.copio_sigdiff.t$OTU_id[arch.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"] #117


length(arch.daOTUs.copio_sigdiff.t$Enriched[arch.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #85
length(arch.daOTUs.copio_sigdiff.t$Enriched[arch.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #120

MYNJ_nonfer_oligo.fun<-fun.daOTUs.copio_sigdiff.t$OTU_id[fun.daOTUs.copio_sigdiff.t$category2 == "Nonfer_oligo"] #105
MYNJ_nonfer_copio.fun<-fun.daOTUs.copio_sigdiff.t$OTU_id[fun.daOTUs.copio_sigdiff.t$category2 == "Nonfer_copio"] #144
MYNJ_fer_oligo.fun<-fun.daOTUs.copio_sigdiff.t$OTU_id[fun.daOTUs.copio_sigdiff.t$category2 == "Fer_oligo"] #119
MYNJ_fer_copio.fun<-fun.daOTUs.copio_sigdiff.t$OTU_id[fun.daOTUs.copio_sigdiff.t$category2 == "Fer_copio"] #117


length(fun.daOTUs.copio_sigdiff.t$Enriched[fun.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]) #571
length(fun.daOTUs.copio_sigdiff.t$Enriched[fun.daOTUs.copio_sigdiff.t$Enriched == "Fer"]) #335



intersect(fun.daOTUs.oligo_sigdiff.t$OTU_id[fun.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"],fun.daOTUs.copio_sigdiff.t$OTU_id[fun.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"])
intersect(bac.daOTUs.oligo_sigdiff.t$OTU_id[bac.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"],bac.daOTUs.copio_sigdiff.t$OTU_id[bac.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"])

intersect(arch.daOTUs.oligo_sigdiff.t$OTU_id[arch.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"],arch.daOTUs.copio_sigdiff.t$OTU_id[arch.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"])


#install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(
  x = list(fun.daOTUs.oligo_sigdiff.t$OTU_id[fun.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"], fun.daOTUs.copio_sigdiff.t$OTU_id[fun.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]),
  category.names = c("CJYS_nonfer" , "MYNJ_nonfer"),
  filename = 'venn_diagramm_common nonfer.sgv',
  output=TRUE
)

venn.diagram(
  x = list(bac.daOTUs.oligo_sigdiff.t$OTU_id[bac.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"],bac.daOTUs.copio_sigdiff.t$OTU_id[bac.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]),
  category.names = c("CJYS_nonfer" , "MYNJ_nonfer"),
  filename = 'venn_diagramm_common nonfer_bac.sgv',
  output=TRUE
)

venn.diagram(
  x = list(arch.daOTUs.oligo_sigdiff.t$OTU_id[arch.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"],arch.daOTUs.copio_sigdiff.t$OTU_id[arch.daOTUs.copio_sigdiff.t$Enriched == "Non-fer"]),
  category.names = c("CJYS_nonfer" , "MYNJ_nonfer"),
  filename = 'venn_diagramm_common nonfer_arch.sgv',
  output=TRUE
)

venn.diagram(
  x = list(fun.daOTUs.oligo_sigdiff.t$OTU_id[fun.daOTUs.oligo_sigdiff.t$Enriched == "Fer"], fun.daOTUs.copio_sigdiff.t$OTU_id[fun.daOTUs.copio_sigdiff.t$Enriched == "Fer"]),
  category.names = c("CJYS_fer" , "MYNJ_fer"),
  filename = 'venn_diagramm_common fer_fun.sgv',
  output=TRUE
)

venn.diagram(
  x = list(bac.daOTUs.oligo_sigdiff.t$OTU_id[bac.daOTUs.oligo_sigdiff.t$Enriched == "Fer"],bac.daOTUs.copio_sigdiff.t$OTU_id[bac.daOTUs.copio_sigdiff.t$Enriched == "Fer"]),
  category.names = c("CJYS_fer" , "MYNJ_fer"),
  filename = 'venn_diagramm_common fer_bac.sgv',
  output=TRUE
)

venn.diagram(
  x = list(arch.daOTUs.oligo_sigdiff.t$OTU_id[arch.daOTUs.oligo_sigdiff.t$Enriched == "Fer"],arch.daOTUs.copio_sigdiff.t$OTU_id[arch.daOTUs.copio_sigdiff.t$Enriched == "Fer"]),
  category.names = c("CJYS_fer" , "MYNJ_fer"),
  filename = 'venn_diagramm_common fer_arch.sgv',
  output=TRUE
)


## Exsistence of common OTUs
intersect(MYNJ_nonfer_oligo.bac, CJYS_nonfer_oligo.bac)
intersect(MYNJ_nonfer_oligo.arch, CJYS_nonfer_oligo.arch)
intersect(MYNJ_nonfer_oligo.fun, CJYS_nonfer_oligo.fun)

intersect(MYNJ_fer_oligo.bac, CJYS_fer_oligo.bac)
intersect(MYNJ_fer_oligo.arch, CJYS_fer_oligo.arch)
intersect(MYNJ_fer_oligo.fun, CJYS_fer_oligo.fun)

intersect(MYNJ_nonfer_copio.bac, CJYS_nonfer_copio.bac)
intersect(MYNJ_nonfer_copio.arch, CJYS_nonfer_copio.arch)
intersect(MYNJ_nonfer_copio.fun, CJYS_nonfer_copio.fun)

intersect(MYNJ_fer_copio.bac, CJYS_fer_copio.bac)
intersect(MYNJ_fer_copio.arch, CJYS_fer_copio.arch)
intersect(MYNJ_fer_copio.fun, CJYS_fer_copio.fun)


tax.bac.nofer.oligo <- subset(OTU.list, OTU_id %in% intersect(MYNJ_nonfer_oligo.bac, CJYS_nonfer_oligo.bac))
tax.bac.nofer.copio <- subset(OTU.list, OTU_id %in% intersect(MYNJ_nonfer_copio.bac, CJYS_nonfer_copio.bac))

tax.arch.nofer.oligo <- subset(OTU.list, OTU_id %in% intersect(MYNJ_nonfer_oligo.arch, CJYS_nonfer_oligo.arch))#none
tax.arch.nofer.copio <- subset(OTU.list, OTU_id %in% intersect(MYNJ_nonfer_copio.arch, CJYS_nonfer_copio.arch))

tax.fun.nofer.oligo <- subset(OTU.list, OTU_id %in% intersect(MYNJ_nonfer_oligo.fun, CJYS_nonfer_oligo.fun))
tax.fun.nofer.copio <- subset(OTU.list, OTU_id %in% intersect(MYNJ_nonfer_copio.fun, CJYS_nonfer_copio.fun))



tax.bac.fer.oligo <- subset(OTU.list, OTU_id %in% intersect(MYNJ_fer_oligo.bac, CJYS_fer_oligo.bac))
tax.bac.fer.copio <- subset(OTU.list, OTU_id %in% intersect(MYNJ_fer_copio.bac, CJYS_fer_copio.bac))

tax.arch.fer.oligo <- subset(OTU.list, OTU_id %in% intersect(MYNJ_fer_oligo.arch, CJYS_fer_oligo.arch))#none
tax.arch.fer.copio <- subset(OTU.list, OTU_id %in% intersect(MYNJ_fer_copio.arch, CJYS_fer_copio.arch))

tax.fun.fer.oligo <- subset(OTU.list, OTU_id %in% intersect(MYNJ_fer_oligo.fun, CJYS_fer_oligo.fun))
tax.fun.fer.copio <- subset(OTU.list, OTU_id %in% intersect(MYNJ_fer_copio.fun, CJYS_fer_copio.fun))




tax.fer.copio <- rbind(tax.bac.fer.copio, tax.arch.fer.copio, tax.fun.fer.copio)
tax.fer.copio$Trophic <- "Copio"
tax.fer.oligo <- rbind(tax.bac.fer.oligo, tax.arch.fer.oligo, tax.fun.fer.oligo)
tax.fer.oligo$Trophic <- "Oligo"

tax.fer <- rbind(tax.fer.copio, tax.fer.oligo)

tax.nofer.copio <- rbind(tax.bac.nofer.copio, tax.arch.nofer.copio, tax.fun.nofer.copio)
tax.nofer.copio$Trophic <- "Copio"
tax.nofer.oligo <- rbind(tax.bac.nofer.oligo, tax.arch.nofer.oligo, tax.fun.nofer.oligo)
tax.nofer.oligo$Trophic <- "Oligo"

tax.nofer <- rbind(tax.nofer.copio, tax.nofer.oligo)


write.table(tax.nofer, "Taxonomy of common response OTUs_Non-fertilized.tsv", sep='\t', quote=F)
write.table(tax.fer, "Taxonomy of common response OTUs_fertilized.tsv", sep='\t', quote=F)


dev.off()

## Different responses depending on condition

intersect(MYNJ_nonfer_oligo.bac, CJYS_nonfer_copio.bac)
intersect(MYNJ_nonfer_copio.bac, CJYS_nonfer_oligo.bac)

intersect(MYNJ_fer_copio.bac, CJYS_nonfer_oligo.bac)
intersect(MYNJ_fer_oligo.bac, CJYS_nonfer_copio.bac)

intersect(MYNJ_nonfer_oligo.arch, CJYS_nonfer_copio.arch)
intersect(MYNJ_nonfer_copio.arch, CJYS_nonfer_oligo.arch)

intersect(MYNJ_fer_copio.arch, CJYS_nonfer_oligo.arch)
intersect(MYNJ_fer_oligo.arch, CJYS_nonfer_copio.arch)

intersect(MYNJ_nonfer_oligo.fun, CJYS_nonfer_copio.fun)
intersect(MYNJ_nonfer_copio.fun, CJYS_nonfer_oligo.fun)

intersect(MYNJ_fer_copio.fun, CJYS_nonfer_oligo.fun)
intersect(MYNJ_fer_oligo.fun, CJYS_nonfer_copio.fun)



intersect(MYNJ_nonfer_oligo.bac, CJYS_fer_oligo.bac)
intersect(MYNJ_fer_oligo.bac, CJYS_nonfer_oligo.bac)


intersect(MYNJ_nonfer_oligo.arch, CJYS_fer_oligo.arch)
intersect(MYNJ_fer_oligo.arch, CJYS_nonfer_oligo.arch)


intersect(MYNJ_nonfer_oligo.fun, CJYS_fer_oligo.fun)
intersect(MYNJ_fer_oligo.fun, CJYS_nonfer_oligo.fun)


intersect(MYNJ_nonfer_copio.bac, CJYS_fer_copio.bac)
intersect(MYNJ_fer_copio.bac, CJYS_nonfer_copio.bac)


intersect(MYNJ_nonfer_copio.arch, CJYS_fer_copio.arch)
intersect(MYNJ_fer_copio.arch, CJYS_nonfer_copio.arch)


intersect(MYNJ_nonfer_copio.fun, CJYS_fer_copio.fun)
intersect(MYNJ_fer_copio.fun, CJYS_nonfer_copio.fun)



## Taxonomic composition of daOTUs and daOTUs with oligo or copio

bac.daOTUs.oligo_sigdiff.t<- read.xlsx("Oligo_daOTUs with oligo and copiotrophs_Bacteria.xlsx",1)
arch.daOTUs.oligo_sigdiff.t <-read.xlsx("Oligo_daOTUs with oligo and copiotrophs_Archaea.xlsx",1)
fun.daOTUs.oligo_sigdiff.t<- read.xlsx( "Oligo_daOTUs with oligo and copiotrophs_Fungi.xlsx",1)

bac.daOTUs.copio_sigdiff.t<-read.xlsx( "Copio_daOTUs with oligo and copiotrophs_Bacteria.xlsx",1)
arch.daOTUs.copio_sigdiff.t<-read.xlsx( "Copio_daOTUs with oligo and copiotrophs_Archaea.xlsx",1)
fun.daOTUs.copio_sigdiff.t<-read.xlsx( "Copio_daOTUs with oligo and copiotrophs_Fungi.xlsx",1)


###Oligotrophic condition
### counts of bacteria OTUs
otu.tax.bac.CJYS <- tax_table(bac.clean.ss.CJ.YS)

length(rownames(otu.tax.bac.CJYS)) #8045
otu.tax.bac.CJYS <- data.frame(otu.tax.bac.CJYS, stringsAsFactors = F)
otu.tax.bac.CJYS$Phylum2 <-otu.tax.bac.CJYS$Phylum

rownames(otu.tax.bac.CJYS)[is.na(rownames(otu.tax.bac.CJYS))]
otu.tax.bac.CJYS$Phylum2[which(otu.tax.bac.CJYS$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
otu.tax.bac.CJYS$Phylum2[which(otu.tax.bac.CJYS$Class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
otu.tax.bac.CJYS$Phylum2[which(otu.tax.bac.CJYS$Class == "Deltaproteobacteria")] <- "Deltaproteobacteria"
otu.tax.bac.CJYS$Phylum2[which(otu.tax.bac.CJYS$Class == "Epsilonproteobacteria")] <- "Epsilonproteobacteria"


bacteria_oligo_nonfer <- as.data.frame(table(otu.tax.bac.CJYS[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"], "Phylum2"] ) )
colnames(bacteria_oligo_nonfer) <- c("Class", "Nonfer_oligo")
bacteria_copio_nonfer <- as.data.frame(table(otu.tax.bac.CJYS[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"], "Phylum2"] ) )
colnames(bacteria_copio_nonfer) <- c("Class", "Nonfer_copio")
bacteria_nonfer <- as.data.frame(table(otu.tax.bac.CJYS[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"], "Phylum2"] ) )
colnames(bacteria_nonfer) <- c("Class", "Nonfer")

bacteria_oligo_fer <- as.data.frame(table(otu.tax.bac.CJYS[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"], "Phylum2"] ) )
colnames(bacteria_oligo_fer) <- c("Class", "Fer_oligo")
bacteria_copio_fer <- as.data.frame(table(otu.tax.bac.CJYS[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"], "Phylum2"] ) )
colnames(bacteria_copio_fer) <- c("Class", "Fer_copio")
bacteria_fer <- as.data.frame(table(otu.tax.bac.CJYS[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$Enriched == "Fer"], "Phylum2"] ) )
colnames(bacteria_fer) <- c("Class", "Fer")


bacteria_copi_oli_nonfer <- merge(bacteria_oligo_nonfer, bacteria_copio_nonfer, all=T, by="Class") 
bacteria_copi_oli_nonfer <- merge(bacteria_copi_oli_nonfer, bacteria_nonfer, all=T, by="Class") 

bacteria_copi_oli_fer <- merge(bacteria_oligo_fer, bacteria_copio_fer, all=T, by="Class") 
bacteria_copi_oli_fer <- merge(bacteria_copi_oli_fer, bacteria_fer, all=T, by="Class") 

bacteria_copi_oli_nonfer_fer <- merge(bacteria_copi_oli_nonfer, bacteria_copi_oli_fer, all=T, by="Class") 

bacteria_all_OTUs <- as.data.frame(table(otu.tax.bac.CJYS[, "Phylum2"] ) )
colnames(bacteria_all_OTUs) <- c("Class", "all bOTUs")
bacteria_modules_daOTUs <- merge(bacteria_copi_oli_nonfer_fer, bacteria_all_OTUs, all=T, by="Class") 
bacteria_modules_daOTUs

bacteria_modules_daOTUs_mat <- bacteria_modules_daOTUs[2:8]
rownames(bacteria_modules_daOTUs_mat) <- bacteria_modules_daOTUs$Class
bacteria_modules_daOTUs_mat[is.na(bacteria_modules_daOTUs_mat)] <- 0
colSums(bacteria_modules_daOTUs_mat)

bacteria_modules_daOTUs_prop <- t(t(bacteria_modules_daOTUs_mat)/colSums(bacteria_modules_daOTUs_mat) ) * 1
bacteria_modules_daOTUs_prop
colSums(bacteria_modules_daOTUs_prop)


bp <- barplot(cbind(bacteria_modules_daOTUs_prop[,1:6], NA, bacteria_modules_daOTUs_prop[,7]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_16s[rownames(bacteria_modules_daOTUs_prop),]$cols )
text(bp, 1.1, labels=c(colSums(bacteria_modules_daOTUs_mat)[1:6], NA,
                       colSums(bacteria_modules_daOTUs_mat)[5]), xpd=T, cex=.6, las=2)
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )

write.xlsx(bacteria_modules_daOTUs_mat, "Taxonomic composition of daOTUs.xlsx")

### counts of archaea OTUs

otu.tax.arch.CJYS <- tax_table(arch.clean.ss.CJ.YS)

length(rownames(otu.tax.arch.CJYS)) #342
otu.tax.arch.CJYS <- data.frame(otu.tax.arch.CJYS, stringsAsFactors = F)
otu.tax.arch.CJYS$Phylum2 <-otu.tax.arch.CJYS$Phylum



archaea_oligo_nonfer <- as.data.frame(table(otu.tax.arch.CJYS[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"], "Phylum2"] ) )
colnames(archaea_oligo_nonfer) <- c("Class", "Nonfer_oligo")
archaea_copio_nonfer <- as.data.frame(table(otu.tax.arch.CJYS[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"], "Phylum2"] ) )
colnames(archaea_copio_nonfer) <- c("Class", "Nonfer_copio")
archaea_nonfer <- as.data.frame(table(otu.tax.arch.CJYS[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"], "Phylum2"] ) )
colnames(archaea_nonfer) <- c("Class", "Nonfer")

archaea_oligo_fer <- as.data.frame(table(otu.tax.arch.CJYS[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"], "Phylum2"] ) )
colnames(archaea_oligo_fer) <- c("Class", "Fer_oligo")
archaea_copio_fer <- as.data.frame(table(otu.tax.arch.CJYS[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"], "Phylum2"] ) )
colnames(archaea_copio_fer) <- c("Class", "Fer_copio")
archaea_fer <- as.data.frame(table(otu.tax.arch.CJYS[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$Enriched == "Fer"], "Phylum2"] ) )
colnames(archaea_fer) <- c("Class", "Fer")


archaea_copi_oli_nonfer <- merge(archaea_oligo_nonfer, archaea_copio_nonfer, all=T, by="Class") 
archaea_copi_oli_nonfer <- merge(archaea_copi_oli_nonfer, archaea_nonfer, all=T, by="Class") 

archaea_copi_oli_fer <- merge(archaea_oligo_fer, archaea_copio_fer, all=T, by="Class") 
archaea_copi_oli_fer <- merge(archaea_copi_oli_fer, archaea_fer, all=T, by="Class") 

archaea_copi_oli_nonfer_fer <- merge(archaea_copi_oli_nonfer, archaea_copi_oli_fer, all=T, by="Class") 

archaea_all_OTUs <- as.data.frame(table(otu.tax.arch.CJYS[, "Phylum2"] ) )
colnames(archaea_all_OTUs) <- c("Class", "all bOTUs")
archaea_modules_daOTUs <- merge(archaea_copi_oli_nonfer_fer, archaea_all_OTUs, all=T, by="Class") 
archaea_modules_daOTUs

archaea_modules_daOTUs_mat <- archaea_modules_daOTUs[2:8]
rownames(archaea_modules_daOTUs_mat) <- archaea_modules_daOTUs$Class
archaea_modules_daOTUs_mat[is.na(archaea_modules_daOTUs_mat)] <- 0
colSums(archaea_modules_daOTUs_mat)

archaea_modules_daOTUs_prop <- t(t(archaea_modules_daOTUs_mat)/colSums(archaea_modules_daOTUs_mat) ) * 1
archaea_modules_daOTUs_prop
colSums(archaea_modules_daOTUs_prop)


ap <- barplot(cbind(archaea_modules_daOTUs_prop[,1:6], NA, archaea_modules_daOTUs_prop[,7]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_arch[rownames(archaea_modules_daOTUs_prop),]$cols )
text(bp, 1.1, labels=c(colSums(archaea_modules_daOTUs_mat)[1:6], NA,
                       colSums(archaea_modules_daOTUs_mat)[5]), xpd=T, cex=.6, las=2)
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_arch_legend$labels), 
       fill=rev(PHYLA_label_cols_arch_legend$cols), 
       border=rev(PHYLA_label_cols_arch_legend$cols) )

write.xlsx(archaea_modules_daOTUs_mat, "Taxonomic composition of daOTUs_archaea.xlsx")


### counts of fungi OTUs
otu.tax.fun.CJYS <- tax_table(fun.clean.ss.CJ.YS)

length(rownames(otu.tax.fun.CJYS)) #2282
otu.tax.fun.CJYS <- data.frame(otu.tax.fun.CJYS, stringsAsFactors = F)
otu.tax.fun.CJYS$Phylum2 <-otu.tax.fun.CJYS$Phylum



fungi_oligo_nonfer <- as.data.frame(table(otu.tax.fun.CJYS[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"], "Phylum2"] ) )
colnames(fungi_oligo_nonfer) <- c("Class", "Nonfer_oligo")
fungi_copio_nonfer <- as.data.frame(table(otu.tax.fun.CJYS[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"], "Phylum2"] ) )
colnames(fungi_copio_nonfer) <- c("Class", "Nonfer_copio")
fungi_nonfer <- as.data.frame(table(otu.tax.fun.CJYS[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"], "Phylum2"] ) )
colnames(fungi_nonfer) <- c("Class", "Nonfer")

fungi_oligo_fer <- as.data.frame(table(otu.tax.fun.CJYS[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"], "Phylum2"] ) )
colnames(fungi_oligo_fer) <- c("Class", "Fer_oligo")
fungi_copio_fer <- as.data.frame(table(otu.tax.fun.CJYS[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"], "Phylum2"] ) )
colnames(fungi_copio_fer) <- c("Class", "Fer_copio")
fungi_fer <- as.data.frame(table(otu.tax.fun.CJYS[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$Enriched == "Fer"], "Phylum2"] ) )
colnames(fungi_fer) <- c("Class", "Fer")


fungi_copi_oli_nonfer <- merge(fungi_oligo_nonfer, fungi_copio_nonfer, all=T, by="Class") 
fungi_copi_oli_nonfer <- merge(fungi_copi_oli_nonfer, fungi_nonfer, all=T, by="Class") 

fungi_copi_oli_fer <- merge(fungi_oligo_fer, fungi_copio_fer, all=T, by="Class") 
fungi_copi_oli_fer <- merge(fungi_copi_oli_fer, fungi_fer, all=T, by="Class") 

fungi_copi_oli_nonfer_fer <- merge(fungi_copi_oli_nonfer, fungi_copi_oli_fer, all=T, by="Class") 

fungi_all_OTUs <- as.data.frame(table(otu.tax.fun.CJYS[, "Phylum2"] ) )
colnames(fungi_all_OTUs) <- c("Class", "all fOTUs")
fungi_modules_daOTUs <- merge(fungi_copi_oli_nonfer_fer, fungi_all_OTUs, all=T, by="Class") 
fungi_modules_daOTUs

fungi_modules_daOTUs_mat <- fungi_modules_daOTUs[2:8]
rownames(fungi_modules_daOTUs_mat) <- fungi_modules_daOTUs$Class
fungi_modules_daOTUs_mat[is.na(fungi_modules_daOTUs_mat)] <- 0
colSums(fungi_modules_daOTUs_mat)

fungi_modules_daOTUs_prop <- t(t(fungi_modules_daOTUs_mat)/colSums(fungi_modules_daOTUs_mat) ) * 1
fungi_modules_daOTUs_prop
colSums(fungi_modules_daOTUs_prop)


fp <- barplot(cbind(fungi_modules_daOTUs_prop[,1:6], NA, fungi_modules_daOTUs_prop[,7]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_its[rownames(fungi_modules_daOTUs_prop),]$cols )
text(bp, 1.1, labels=c(colSums(fungi_modules_daOTUs_mat)[1:6], NA,
                       colSums(fungi_modules_daOTUs_mat)[5]), xpd=T, cex=.6, las=2)
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_its_legend$labels), 
       fill=rev(PHYLA_label_cols_its_legend$cols), 
       border=rev(PHYLA_label_cols_its_legend$cols) )

write.xlsx(fungi_modules_daOTUs_mat, "Taxonomic composition of daOTUs_fungi.xlsx")






###Eutrophic condition
### counts of bacteria OTUs
otu.tax.bac.MYNJ <- tax_table(bac.clean.ss.MY.NJ)

length(rownames(otu.tax.bac.MYNJ)) #7969
otu.tax.bac.MYNJ <- data.frame(otu.tax.bac.MYNJ, stringsAsFactors = F)
otu.tax.bac.MYNJ$Phylum2 <-otu.tax.bac.MYNJ$Phylum

rownames(otu.tax.bac.MYNJ)[is.na(rownames(otu.tax.bac.MYNJ))]
otu.tax.bac.MYNJ$Phylum2[which(otu.tax.bac.MYNJ$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
otu.tax.bac.MYNJ$Phylum2[which(otu.tax.bac.MYNJ$Class == "Gammaproteobacteria")] <- "Gammaproteobacteria"
otu.tax.bac.MYNJ$Phylum2[which(otu.tax.bac.MYNJ$Class == "Deltaproteobacteria")] <- "Deltaproteobacteria"
otu.tax.bac.MYNJ$Phylum2[which(otu.tax.bac.MYNJ$Class == "Epsilonproteobacteria")] <- "Epsilonproteobacteria"


bacteria_oligo_nonfer <- as.data.frame(table(otu.tax.bac.MYNJ[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"], "Phylum2"] ) )
colnames(bacteria_oligo_nonfer) <- c("Class", "Nonfer_oligo")
bacteria_copio_nonfer <- as.data.frame(table(otu.tax.bac.MYNJ[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"], "Phylum2"] ) )
colnames(bacteria_copio_nonfer) <- c("Class", "Nonfer_copio")
bacteria_nonfer <- as.data.frame(table(otu.tax.bac.MYNJ[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"], "Phylum2"] ) )
colnames(bacteria_nonfer) <- c("Class", "Nonfer")

bacteria_oligo_fer <- as.data.frame(table(otu.tax.bac.MYNJ[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"], "Phylum2"] ) )
colnames(bacteria_oligo_fer) <- c("Class", "Fer_oligo")
bacteria_copio_fer <- as.data.frame(table(otu.tax.bac.MYNJ[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"], "Phylum2"] ) )
colnames(bacteria_copio_fer) <- c("Class", "Fer_copio")
bacteria_fer <- as.data.frame(table(otu.tax.bac.MYNJ[bac.daOTUs.oligo_sigdiff.t$OTU[bac.daOTUs.oligo_sigdiff.t$Enriched == "Fer"], "Phylum2"] ) )
colnames(bacteria_fer) <- c("Class", "Fer")


bacteria_copi_oli_nonfer <- merge(bacteria_oligo_nonfer, bacteria_copio_nonfer, all=T, by="Class") 
bacteria_copi_oli_nonfer <- merge(bacteria_copi_oli_nonfer, bacteria_nonfer, all=T, by="Class") 

bacteria_copi_oli_fer <- merge(bacteria_oligo_fer, bacteria_copio_fer, all=T, by="Class") 
bacteria_copi_oli_fer <- merge(bacteria_copi_oli_fer, bacteria_fer, all=T, by="Class") 

bacteria_copi_oli_nonfer_fer <- merge(bacteria_copi_oli_nonfer, bacteria_copi_oli_fer, all=T, by="Class") 

bacteria_all_OTUs <- as.data.frame(table(otu.tax.bac.MYNJ[, "Phylum2"] ) )
colnames(bacteria_all_OTUs) <- c("Class", "all bOTUs")
bacteria_modules_daOTUs <- merge(bacteria_copi_oli_nonfer_fer, bacteria_all_OTUs, all=T, by="Class") 
bacteria_modules_daOTUs

bacteria_modules_daOTUs_mat <- bacteria_modules_daOTUs[2:8]
rownames(bacteria_modules_daOTUs_mat) <- bacteria_modules_daOTUs$Class
bacteria_modules_daOTUs_mat[is.na(bacteria_modules_daOTUs_mat)] <- 0
colSums(bacteria_modules_daOTUs_mat)

bacteria_modules_daOTUs_prop <- t(t(bacteria_modules_daOTUs_mat)/colSums(bacteria_modules_daOTUs_mat) ) * 1
bacteria_modules_daOTUs_prop
colSums(bacteria_modules_daOTUs_prop)


bp <- barplot(cbind(bacteria_modules_daOTUs_prop[,1:6], NA, bacteria_modules_daOTUs_prop[,7]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_16s[rownames(bacteria_modules_daOTUs_prop),]$cols )
text(bp, 1.1, labels=c(colSums(bacteria_modules_daOTUs_mat)[1:6], NA,
                       colSums(bacteria_modules_daOTUs_mat)[5]), xpd=T, cex=.6, las=2)
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )

write.xlsx(bacteria_modules_daOTUs_mat, "Taxonomic composition of daOTUs_bacteria_eutrophic.xlsx")

### counts of archaea OTUs

otu.tax.arch.MYNJ <- tax_table(arch.clean.ss.MY.NJ)

length(rownames(otu.tax.arch.MYNJ)) #342
otu.tax.arch.MYNJ <- data.frame(otu.tax.arch.MYNJ, stringsAsFactors = F)
otu.tax.arch.MYNJ$Phylum2 <-otu.tax.arch.MYNJ$Phylum



archaea_oligo_nonfer <- as.data.frame(table(otu.tax.arch.MYNJ[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"], "Phylum2"] ) )
colnames(archaea_oligo_nonfer) <- c("Class", "Nonfer_oligo")
archaea_copio_nonfer <- as.data.frame(table(otu.tax.arch.MYNJ[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"], "Phylum2"] ) )
colnames(archaea_copio_nonfer) <- c("Class", "Nonfer_copio")
archaea_nonfer <- as.data.frame(table(otu.tax.arch.MYNJ[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"], "Phylum2"] ) )
colnames(archaea_nonfer) <- c("Class", "Nonfer")

archaea_oligo_fer <- as.data.frame(table(otu.tax.arch.MYNJ[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"], "Phylum2"] ) )
colnames(archaea_oligo_fer) <- c("Class", "Fer_oligo")
archaea_copio_fer <- as.data.frame(table(otu.tax.arch.MYNJ[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"], "Phylum2"] ) )
colnames(archaea_copio_fer) <- c("Class", "Fer_copio")
archaea_fer <- as.data.frame(table(otu.tax.arch.MYNJ[arch.daOTUs.oligo_sigdiff.t$OTU[arch.daOTUs.oligo_sigdiff.t$Enriched == "Fer"], "Phylum2"] ) )
colnames(archaea_fer) <- c("Class", "Fer")


archaea_copi_oli_nonfer <- merge(archaea_oligo_nonfer, archaea_copio_nonfer, all=T, by="Class") 
archaea_copi_oli_nonfer <- merge(archaea_copi_oli_nonfer, archaea_nonfer, all=T, by="Class") 

archaea_copi_oli_fer <- merge(archaea_oligo_fer, archaea_copio_fer, all=T, by="Class") 
archaea_copi_oli_fer <- merge(archaea_copi_oli_fer, archaea_fer, all=T, by="Class") 

archaea_copi_oli_nonfer_fer <- merge(archaea_copi_oli_nonfer, archaea_copi_oli_fer, all=T, by="Class") 

archaea_all_OTUs <- as.data.frame(table(otu.tax.arch.MYNJ[, "Phylum2"] ) )
colnames(archaea_all_OTUs) <- c("Class", "all aOTUs")
archaea_modules_daOTUs <- merge(archaea_copi_oli_nonfer_fer, archaea_all_OTUs, all=T, by="Class") 
archaea_modules_daOTUs

archaea_modules_daOTUs_mat <- archaea_modules_daOTUs[2:8]
rownames(archaea_modules_daOTUs_mat) <- archaea_modules_daOTUs$Class
archaea_modules_daOTUs_mat[is.na(archaea_modules_daOTUs_mat)] <- 0
colSums(archaea_modules_daOTUs_mat)

archaea_modules_daOTUs_prop <- t(t(archaea_modules_daOTUs_mat)/colSums(archaea_modules_daOTUs_mat) ) * 1
archaea_modules_daOTUs_prop
colSums(archaea_modules_daOTUs_prop)


ap <- barplot(cbind(archaea_modules_daOTUs_prop[,1:6], NA, archaea_modules_daOTUs_prop[,7]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_arch[rownames(archaea_modules_daOTUs_prop),]$cols )
text(bp, 1.1, labels=c(colSums(archaea_modules_daOTUs_mat)[1:6], NA,
                       colSums(archaea_modules_daOTUs_mat)[5]), xpd=T, cex=.6, las=2)
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_arch_legend$labels), 
       fill=rev(PHYLA_label_cols_arch_legend$cols), 
       border=rev(PHYLA_label_cols_arch_legend$cols) )

write.xlsx(archaea_modules_daOTUs_mat, "Taxonomic composition of daOTUs_archaea_eutrophic.xlsx")


### counts of fungi OTUs
otu.tax.fun.MYNJ <- tax_table(fun.clean.ss.MY.NJ)

length(rownames(otu.tax.fun.MYNJ)) #2282
otu.tax.fun.MYNJ <- data.frame(otu.tax.fun.MYNJ, stringsAsFactors = F)
otu.tax.fun.MYNJ$Phylum2 <-otu.tax.fun.MYNJ$Phylum



fungi_oligo_nonfer <- as.data.frame(table(otu.tax.fun.MYNJ[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_oligo"], "Phylum2"] ) )
colnames(fungi_oligo_nonfer) <- c("Class", "Nonfer_oligo")
fungi_copio_nonfer <- as.data.frame(table(otu.tax.fun.MYNJ[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$category2 == "Nonfer_copio"], "Phylum2"] ) )
colnames(fungi_copio_nonfer) <- c("Class", "Nonfer_copio")
fungi_nonfer <- as.data.frame(table(otu.tax.fun.MYNJ[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$Enriched == "Non-fer"], "Phylum2"] ) )
colnames(fungi_nonfer) <- c("Class", "Nonfer")

fungi_oligo_fer <- as.data.frame(table(otu.tax.fun.MYNJ[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$category2 == "Fer_oligo"], "Phylum2"] ) )
colnames(fungi_oligo_fer) <- c("Class", "Fer_oligo")
fungi_copio_fer <- as.data.frame(table(otu.tax.fun.MYNJ[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$category2 == "Fer_copio"], "Phylum2"] ) )
colnames(fungi_copio_fer) <- c("Class", "Fer_copio")
fungi_fer <- as.data.frame(table(otu.tax.fun.MYNJ[fun.daOTUs.oligo_sigdiff.t$OTU[fun.daOTUs.oligo_sigdiff.t$Enriched == "Fer"], "Phylum2"] ) )
colnames(fungi_fer) <- c("Class", "Fer")


fungi_copi_oli_nonfer <- merge(fungi_oligo_nonfer, fungi_copio_nonfer, all=T, by="Class") 
fungi_copi_oli_nonfer <- merge(fungi_copi_oli_nonfer, fungi_nonfer, all=T, by="Class") 

fungi_copi_oli_fer <- merge(fungi_oligo_fer, fungi_copio_fer, all=T, by="Class") 
fungi_copi_oli_fer <- merge(fungi_copi_oli_fer, fungi_fer, all=T, by="Class") 

fungi_copi_oli_nonfer_fer <- merge(fungi_copi_oli_nonfer, fungi_copi_oli_fer, all=T, by="Class") 

fungi_all_OTUs <- as.data.frame(table(otu.tax.fun.MYNJ[, "Phylum2"] ) )
colnames(fungi_all_OTUs) <- c("Class", "all fOTUs")
fungi_modules_daOTUs <- merge(fungi_copi_oli_nonfer_fer, fungi_all_OTUs, all=T, by="Class") 
fungi_modules_daOTUs

fungi_modules_daOTUs_mat <- fungi_modules_daOTUs[2:8]
rownames(fungi_modules_daOTUs_mat) <- fungi_modules_daOTUs$Class
fungi_modules_daOTUs_mat[is.na(fungi_modules_daOTUs_mat)] <- 0
colSums(fungi_modules_daOTUs_mat)

fungi_modules_daOTUs_prop <- t(t(fungi_modules_daOTUs_mat)/colSums(fungi_modules_daOTUs_mat) ) * 1
fungi_modules_daOTUs_prop
colSums(fungi_modules_daOTUs_prop)


fp <- barplot(cbind(fungi_modules_daOTUs_prop[,1:6], NA, fungi_modules_daOTUs_prop[,7]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_its[rownames(fungi_modules_daOTUs_prop),]$cols )
text(bp, 1.1, labels=c(colSums(fungi_modules_daOTUs_mat)[1:6], NA,
                       colSums(fungi_modules_daOTUs_mat)[5]), xpd=T, cex=.6, las=2)
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_its_legend$labels), 
       fill=rev(PHYLA_label_cols_its_legend$cols), 
       border=rev(PHYLA_label_cols_its_legend$cols) )

write.xlsx(fungi_modules_daOTUs_mat, "Taxonomic composition of daOTUs_fungi_eutrophic.xlsx")
