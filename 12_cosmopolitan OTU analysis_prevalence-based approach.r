##Read OTU tables and convert to phyloseq object
##Merge samples based on soil sites
bac.clean.ss.field <-  merge_samples(bac.clean.ss, "Field")
t(otu_table(bac.clean.ss.field))

arch.clean.ss.field <-  merge_samples(arch.clean.ss, "Field")
t(otu_table(arch.clean.ss.field))

fun.clean.ss.field <-  merge_samples(fun.clean.ss, "Field")
t(otu_table(fun.clean.ss.field))


# Calculate compositional version of the data
# (relative abundances)
bac.clean.ss.field.rel <- microbiome::transform(bac.clean.ss.field, "compositional")
t(otu_table(bac.clean.ss.field.rel))

arch.clean.ss.field.rel <- microbiome::transform(arch.clean.ss.field, "compositional")
t(otu_table(arch.clean.ss.field.rel))

fun.clean.ss.field.rel <- microbiome::transform(fun.clean.ss.field, "compositional")
t(otu_table(fun.clean.ss.field.rel))

# Filter the data to include only healthy subjects
#pseq.1 <- subset_samples(phy.clean.l, ibd_subtype == "HC" & timepoint == "1") 
#print(pseq.1)
# keep only taxa with positive sums
#pseq.2 <- prune_taxa(taxa_sums(pseq.1) > 0, pseq.1)

#print(pseq.2)

#Relative population frequencies; at 1% compositional abundance threshold:
  
  head(prevalence(bac.clean.ss.field.rel, detection = 0.01, sort = TRUE))
  head(prevalence(arch.clean.ss.field.rel, detection = 0.01, sort = TRUE))
  head(prevalence(fun.clean.ss.field.rel, detection = 0.01, sort = TRUE))
  
  
#We can see that only OTU ids are listed with no taxonomic information. Absolute population frequencies (sample count):
  
  head(prevalence(bac.clean.ss.field.rel, detection = 0.01, sort = TRUE, count = TRUE))
#Core microbiota analysis
#If you only need the names of the core taxa, do as follows. This returns the taxa that exceed the given prevalence and detection thresholds.

core.bac.99 <- core_members(bac.clean.ss.field.rel, detection = 0, prevalence = 99/100)
core.bac.99

core.bac <- subset(OTU_id.list, OTU_id.list$OTU%in% core.bac.99)

core.arch.99 <- core_members(arch.clean.ss.field.rel, detection = 0, prevalence = 99/100)
core.arch.99

core.arch <- subset(OTU_id.list, OTU_id.list$OTU%in% core.arch.99)

core.fun.99 <- core_members(fun.clean.ss.field.rel, detection = 0, prevalence = 99/100)
core.fun.99

core.fun <- subset(OTU_id.list, OTU_id.list$OTU%in% core.fun.99)



### Correlation between core OTU abundance and chemical properties
#2017 full data
df.otu.bac.clean.ss
df.otu.arch.clean.ss
df.fun.clean.ss

df.otu.core.bac <- subset(df.otu.bac.clean.ss, rownames(df.otu.bac.clean.ss)%in% core.bac$OTU_id)
df.otu.core.arch <- subset(df.otu.arch.clean.ss, rownames(df.otu.arch.clean.ss)%in% core.arch$OTU_id)
df.otu.core.fun <- subset(df.fun.clean.ss, rownames(df.fun.clean.ss)%in% core.fun$OTU_id)


b.meta.all
f.meta.all

cor_table.core.bac<-cbind(b.meta[,7:18], t(df.otu.core.bac))
summary(cor_table.core.bac)

cor_table.core.arch<-cbind(b.meta[,7:18], t(df.otu.core.arch))
summary(cor_table.core.arch)

cor_table.core.fun<-cbind(f.meta[,7:18], t(df.otu.core.fun))
summary(cor_table.core.fun)


##Correlation
mat.cor_table.core.bac <- as.matrix(cor_table.core.bac)
mat.cor_table.core.arch <- as.matrix(cor_table.core.arch)
mat.cor_table.core.fun <- as.matrix(cor_table.core.fun)


library(Hmisc)

cor_5 <- rcorr(mat.cor_table.core.bac, type="spearman")
M <- cor_5$r
p_mat <- cor_5$P

library(xlsx)
write.xlsx(M, "Bac_core_correlation_rho.xlsx")
write.xlsx(p_mat, "Bac_core_correlation_p value.xlsx")

corrplot::corrplot(M, type = "upper", 
                   p.mat = p_mat, sig.level = 0.05, insig = "blank", tl.col="black", col=brewer.pal(n=10, name="RdYlBu"))

corrplot::corrplot(M, type = "upper", 
                   tl.col="black", col=brewer.pal(n=10, name="RdYlBu"))

dev.off()