## Defining daOTUs
#Phyloseq files
bac.clean.ss
arch.clean.ss
fun.clean.ss


zero.clean.filt <- filter_taxa(fun.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(zero.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(zero.clean.filt))

## make binary samples
zero.clean.filt.1 <- subset_samples(zero.clean.filt, Cultural_practice %in% c('Conventional','No_fertilizer'))
zero.clean.filt.2 <- subset_samples(zero.clean.filt, Cultural_practice %in% c('Conventional', 'No_pesticide'))
zero.clean.filt.3 <- subset_samples(zero.clean.filt, Cultural_practice %in% c('No_fertilizer', 'No_pesticide'))

zero.clean.filt.1 <- filter_taxa(zero.clean.filt.1, function(x) sum(x) != 0, TRUE)
zero.clean.filt.2 <- filter_taxa(zero.clean.filt.2, function(x) sum(x) != 0, TRUE)
zero.clean.filt.3 <- filter_taxa(zero.clean.filt.3, function(x) sum(x) != 0, TRUE)

sample_names(zero.clean.filt.1)
# only keep samples over 1000 read
(filt.sample <- sample_sums(zero.clean.filt.1) > 100)
sum(sample_sums(zero.clean.filt.1) <= 100)  ## 0 sample discarded
zero.clean.filt.1 <- prune_samples(filt.sample, zero.clean.filt.1)
## 40x difference in sequencing depth between samples. 
#CSS normalization is not going to ??fix?? this much difference. 
#In this case I would recommend rarefaction, but we will move forward since this is the example dataset in phyloseq
obj <- phyloseq_to_metagenomeSeq(zero.clean.filt.1)

# obj <- filterData(obj, depth = 100)
obj   ## removed 0 samples below 100 reads
##
## fitZig sample
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 30, verbose = TRUE)
Type  <-  pData(obj)$Cultural_practice
mod  <-  model.matrix(~Type)
colnames(mod)  <-  levels(Type)
colnames(mod)

res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)


zigFit = res$fit
finalMod= res$fit$design


contrast.matrix.1 =makeContrasts(Conventional - No_fertilizer, levels = finalMod)
contrast.matrix.2 =makeContrasts(Conventional - No_pesticide, levels = finalMod)
contrast.matrix.3 =makeContrasts(No_fertilizer - No_pesticide, levels = finalMod)

fit2_con_nofertil = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 
fit2_con_nopesti = contrasts.fit(zigFit, contrasts=contrast.matrix.2)
fit2_nofertil_nopesti = contrasts.fit(zigFit, contrasts=contrast.matrix.3)


fit3.1 = eBayes(fit2_con_nofertil)
fit3.2 = eBayes(fit2_con_nopesti)
fit3.3 = eBayes(fit2_nofertil_nopesti)
fit3
topTable(fit3.1, coef="Conventional - No_fertilizer")
topTable(fit3.2, coef="Conventional - No_pesticide")
topTable(fit3.3, coef="No_fertilizer - No_pesticide")

res.1 <- topTable(fit3.1,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res)
dim(res)

write.xlsx(res.1, 'Fun_daOTU_CF and NF.xlsx')

write.xlsx(res.2, 'Fun_daOTU_CF and NP.xlsx')

write.xlsx(res.3, 'Fun_daOTU_NF and NP.xlsx')


res.CF.NF<-read.xlsx('Fun_daOTU_CF and NF.xlsx', 1)
res.CF.NP<-read.xlsx('Fun_daOTU_CF and NP.xlsx', 1)
res.NF.NP<-read.xlsx('Fun_daOTU_NF and NP.xlsx', 1)

names(res.CF.NF)[1] <- c("OTU")
names(res.CF.NP)[1] <- c("OTU")
names(res.NF.NP)[1] <- c("OTU")

res.CF.NF.id <- merge(res.CF.NF, OTU_id.list, by= "OTU")
res.CF.NP.id <- merge(res.CF.NP, OTU_id.list, by= "OTU")
res.NF.NP.id <- merge(res.NF.NP, OTU_id.list, by= "OTU")

write.xlsx(res.CF.NF.id, 'Fun_daOTU_CF and NF.xlsx')

write.xlsx(res.CF.NP.id, 'Fun_daOTU_CF and NP.xlsx')

write.xlsx(res.NF.NP.id, 'Fun_daOTU_NF and NP.xlsx')

res.CF.NF.CF_enriched <- res.CF.NF.id$OTU_id[which(res.CF.NF.id$logFC>0)]
res.CF.NF.NF_enriched <- res.CF.NF.id$OTU_id[which(res.CF.NF.id$logFC<0)]

res.CF.NP.CF_enriched <- res.CF.NP.id$OTU_id[which(res.CF.NP.id$logFC>0)]
res.CF.NP.NP_enriched <- res.CF.NP.id$OTU_id[which(res.CF.NP.id$logFC<0)]

res.NF.NP.NF_enriched <- res.NF.NP.id$OTU_id[which(res.NF.NP.id$logFC>0)]
res.NF.NP.NP_enriched <- res.NF.NP.id$OTU_id[which(res.NF.NP.id$logFC<0)]

##daOTUs list
daOTU_conven_bac.all <-intersect(res.CF.NF.CF_enriched, res.CF.NP.CF_enriched)
daOTU_nofertil_bac.all <- intersect(res.CF.NF.NF_enriched, res.NF.NP.NF_enriched)
daOTU_nopesti_bac.all <- intersect(res.CF.NP.NP_enriched, res.NF.NP.NP_enriched)

daOTU_conven_fun.all <-intersect(res.CF.NF.CF_enriched, res.CF.NP.CF_enriched)
daOTU_nofertil_fun.all <- intersect(res.CF.NF.NF_enriched, res.NF.NP.NF_enriched)
daOTU_nopesti_fun.all <- intersect(res.CF.NP.NP_enriched, res.NF.NP.NP_enriched)

daOTU_conven_arch.all <-intersect(res.CF.NF.CF_enriched, res.CF.NP.CF_enriched)
daOTU_nofertil_arch.all <- intersect(res.CF.NF.NF_enriched, res.NF.NP.NF_enriched)
daOTU_nopesti_arch.all <- intersect(res.CF.NP.NP_enriched, res.NF.NP.NP_enriched)


# Ternary plot
### ternary plots
##substitute otu names to OTU_id
otu.bac.clean.log <- otu_table(bac.clean.log)
df.otu.bac.clean.log <- data.frame(otu.bac.clean.log)
df.otu.bac.clean.log$OTU <-rownames(df.otu.bac.clean.log)
df.otu.bac.clean.log <- merge(df.otu.bac.clean.log, OTU_id.list, by = "OTU")
rownames(df.otu.bac.clean.log) <- df.otu.bac.clean.log$OTU_id
df.otu.bac.clean.log <- df.otu.bac.clean.log[-c(1,254)]

otu.arch.clean.log <- otu_table(arch.clean.log)
df.otu.arch.clean.log <- data.frame(otu.arch.clean.log)
df.otu.arch.clean.log$OTU <-rownames(df.otu.arch.clean.log)
df.otu.arch.clean.log <- merge(df.otu.arch.clean.log, OTU_id.list, by = "OTU")
rownames(df.otu.arch.clean.log) <- df.otu.arch.clean.log$OTU_id
df.otu.arch.clean.log <- df.otu.arch.clean.log[-c(1,254)]

otu.fun.clean.log <- otu_table(fun.clean.log)
df.otu.fun.clean.log <- data.frame(otu.fun.clean.log)
df.otu.fun.clean.log$OTU <-rownames(df.otu.fun.clean.log)
df.otu.fun.clean.log <- merge(df.otu.fun.clean.log, OTU_id.list, by = "OTU")
rownames(df.otu.fun.clean.log) <- df.otu.fun.clean.log$OTU_id
df.otu.fun.clean.log <- df.otu.fun.clean.log[-c(1,254)]

##calling DA otus
daOTU_conven_bac.all 
daOTU_nofertil_bac.all 
daOTU_nopesti_bac.all 

daOTU_conven_fun.all 
daOTU_nofertil_fun.all 
daOTU_nopesti_fun.all 

daOTU_conven_arch.all 
daOTU_nofertil_arch.all 
daOTU_nopesti_arch.all


##Taxonomy of daOTUs
get_taxonomy <- function(daotu_list, front_number, cropping_practice, kingdoms){
              taxonomy_table  <- subset(otu.list, OTU_id%in%daotu_list)
              taxonomy_table$Enriched_in <- cropping_practice
              tablename <- paste0(as.character(front_number),'_',cropping_practice,'_','daOTU_',kingdoms,'.tsv')
              write.table(taxonomy_table, file=tablename, quote=FALSE, sep='\t', row.names = F)
}

get_taxonomy(daOTU_conven_bac.all, 1, "CF", "Bac")
get_taxonomy(daOTU_nofertil_bac.all, 2, "NF", "Bac")
get_taxonomy(daOTU_nopesti_bac.all, 3, "NP", "Bac")

get_taxonomy(daOTU_conven_fun.all, 4, "CF", "Fun")
get_taxonomy(daOTU_nofertil_fun.all, 5, "NF", "Fun")
get_taxonomy(daOTU_nopesti_fun.all, 6, "NP", "Fun")

get_taxonomy(daOTU_conven_arch.all, 7, "CF", "Arch")
get_taxonomy(daOTU_nofertil_arch.all, 8, "NF", "Arch")
get_taxonomy(daOTU_nopesti_arch.all, 9, "NP", "Arch")



otu_table_norm_log_css <- df.otu.fun.clean.log


# create vectors of mean reltive abundances
library(microbiome)
design <- b.meta.all

idx <- design$Cultural_practice %in% c("Conventional", "No_pesticide", "No_fertilizer") 

design_subset <- design[idx, ]


idx <- design_subset$Cultural_practice=="Conventional"
CF_means <- apply(otu_table_norm_log_css[, idx], 1, mean)

idx <- design_subset$Cultural_practice=="No_fertilizer"
NF_means <- apply(otu_table_norm_log_css[, idx], 1, mean)

idx <- design_subset$Cultural_practice=="No_pesticide"
NP_means <- apply(otu_table_norm_log_css[, idx], 1, mean)

# create matrix of average r.a. per group

df <- data.frame(No_fertilizer=NF_means, Conventional=CF_means, No_pesticide=NP_means)
#df <- df[rowSums(df)!=0, ]
#df <- log2(df * scale + 1)

# sort the rows by decreasing abundance (looks better)

idx <- sort(rowSums(df), decreasing=F, index.return=T)$ix
df <- df[idx, ]
df
# create vector of colors according to enrichment

colors <- rep(c_grey, dim(df)[1])
colors[rownames(df) %in% daOTU_nofertil_fun.all] <- c_very_dark_green
colors[rownames(df) %in% daOTU_nopesti_fun.all] <- c_blue
colors[rownames(df) %in% daOTU_conven_fun.all] <- c_red

idx <- sort(colors==c_grey, decreasing=T, index.return=T)$ix
df <- df[idx, ]
colors <- colors[idx]
df

write.table(df, 'ternary plot_cultural practice_fungi.txt', quote=FALSE, sep='\t', row.names = T, col.names = T)
# plot colored by enrichment

tern_e(df, prop_size=T, col=colors, grid_color="grey",
       labels_color="transparent", pch=19, main="Cultural practice")
