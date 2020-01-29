## Network sorted by fields
all.clean.nolog.id <-all.clean.nolog
otu.all.clean.nolog.id<- otu_table(all.clean.nolog.id)
df.otu.all.clean.nolog.id <- data.frame(otu.all.clean.nolog.id)
df.otu.all.clean.nolog.id$OTU <- rownames(df.otu.all.clean.nolog.id)
df.otu.all.clean.nolog.id <- merge(df.otu.all.clean.nolog.id, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.id) <- df.otu.all.clean.nolog.id$OTU_id
length(colnames(df.otu.all.clean.nolog.id)) #344
df.otu.all.clean.nolog.id.trim <- df.otu.all.clean.nolog.id[-c(1, 344)]

otu_table(all.clean.nolog.id) <- otu_table(df.otu.all.clean.nolog.id.trim, taxa_are_rows = T)

##CC1
all.clean.nolog.cc1 <- subset_samples(all.clean.nolog, Field == "CC1"|Field == "CC1.18")
all.clean.nolog.cc1 <- phyloseq::filter_taxa(all.clean.nolog.cc1, function(x) sum(x) != 0, TRUE)

##CC2
all.clean.nolog.cc2 <- subset_samples(all.clean.nolog, Field == "CC2"|Field == "CC2.18")
all.clean.nolog.cc2 <- phyloseq::filter_taxa(all.clean.nolog.cc2, function(x) sum(x) != 0, TRUE)

##CJ1
all.clean.nolog.cj1 <- subset_samples(all.clean.nolog, Field == "CJ1"|Field == "CJ1.18")
all.clean.nolog.cj1 <- phyloseq::filter_taxa(all.clean.nolog.cj1, function(x) sum(x) != 0, TRUE)

##CJ2
all.clean.nolog.cj2 <- subset_samples(all.clean.nolog, Field == "CJ2"|Field == "CJ2.18")
all.clean.nolog.cj2 <- phyloseq::filter_taxa(all.clean.nolog.cj2, function(x) sum(x) != 0, TRUE)

##DG1
all.clean.nolog.dg1 <- subset_samples(all.clean.nolog, Field == "DG1"|Field == "DG1.18")
all.clean.nolog.dg1 <- phyloseq::filter_taxa(all.clean.nolog.dg1, function(x) sum(x) != 0, TRUE)

##DG2
all.clean.nolog.dg2 <- subset_samples(all.clean.nolog, Field == "DG2"|Field == "DG2.18")
all.clean.nolog.dg2 <- phyloseq::filter_taxa(all.clean.nolog.dg2, function(x) sum(x) != 0, TRUE)

##IS1
all.clean.nolog.is1 <- subset_samples(all.clean.nolog, Field == "IS1"|Field == "IS1.18")
all.clean.nolog.is1 <- phyloseq::filter_taxa(all.clean.nolog.is1, function(x) sum(x) != 0, TRUE)

##IS2
all.clean.nolog.is2 <- subset_samples(all.clean.nolog, Field == "IS2"|Field == "IS2.18")
all.clean.nolog.is2 <- phyloseq::filter_taxa(all.clean.nolog.is2, function(x) sum(x) != 0, TRUE)


## OTU tables for each kingdom and soil site
bac.clean.nolog
fun.clean.nolog
arch.clean.nolog
##Bacteria
##CC1
bac.clean.nolog.cc1 <- subset_samples(bac.clean.nolog, Field == "CC1"|Field == "CC1.18")
bac.clean.nolog.cc1 <- phyloseq::filter_taxa(bac.clean.nolog.cc1, function(x) sum(x) != 0, TRUE)

##CC2
bac.clean.nolog.cc2 <- subset_samples(bac.clean.nolog, Field == "CC2"|Field == "CC2.18")
bac.clean.nolog.cc2 <- phyloseq::filter_taxa(bac.clean.nolog.cc2, function(x) sum(x) != 0, TRUE)

##CJ1
bac.clean.nolog.cj1 <- subset_samples(bac.clean.nolog, Field == "CJ1"|Field == "CJ1.18")
bac.clean.nolog.cj1 <- phyloseq::filter_taxa(bac.clean.nolog.cj1, function(x) sum(x) != 0, TRUE)

##CJ2
bac.clean.nolog.cj2 <- subset_samples(bac.clean.nolog, Field == "CJ2"|Field == "CJ2.18")
bac.clean.nolog.cj2 <- phyloseq::filter_taxa(bac.clean.nolog.cj2, function(x) sum(x) != 0, TRUE)

##DG1
bac.clean.nolog.dg1 <- subset_samples(bac.clean.nolog, Field == "DG1"|Field == "DG1.18")
bac.clean.nolog.dg1 <- phyloseq::filter_taxa(bac.clean.nolog.dg1, function(x) sum(x) != 0, TRUE)

##DG2
bac.clean.nolog.dg2 <- subset_samples(bac.clean.nolog, Field == "DG2"|Field == "DG2.18")
bac.clean.nolog.dg2 <- phyloseq::filter_taxa(bac.clean.nolog.dg2, function(x) sum(x) != 0, TRUE)

##IS1
bac.clean.nolog.is1 <- subset_samples(bac.clean.nolog, Field == "IS1"|Field == "IS1.18")
bac.clean.nolog.is1 <- phyloseq::filter_taxa(bac.clean.nolog.is1, function(x) sum(x) != 0, TRUE)

##IS2
bac.clean.nolog.is2 <- subset_samples(bac.clean.nolog, Field == "IS2"|Field == "IS2.18")
bac.clean.nolog.is2 <- phyloseq::filter_taxa(bac.clean.nolog.is2, function(x) sum(x) != 0, TRUE)

##Archaea
##CC1
arch.clean.nolog.cc1 <- subset_samples(arch.clean.nolog, Field == "CC1"|Field == "CC1.18")
arch.clean.nolog.cc1 <- phyloseq::filter_taxa(arch.clean.nolog.cc1, function(x) sum(x) != 0, TRUE)

##CC2
arch.clean.nolog.cc2 <- subset_samples(arch.clean.nolog, Field == "CC2"|Field == "CC2.18")
arch.clean.nolog.cc2 <- phyloseq::filter_taxa(arch.clean.nolog.cc2, function(x) sum(x) != 0, TRUE)

##CJ1
arch.clean.nolog.cj1 <- subset_samples(arch.clean.nolog, Field == "CJ1"|Field == "CJ1.18")
arch.clean.nolog.cj1 <- phyloseq::filter_taxa(arch.clean.nolog.cj1, function(x) sum(x) != 0, TRUE)

##CJ2
arch.clean.nolog.cj2 <- subset_samples(arch.clean.nolog, Field == "CJ2"|Field == "CJ2.18")
arch.clean.nolog.cj2 <- phyloseq::filter_taxa(arch.clean.nolog.cj2, function(x) sum(x) != 0, TRUE)

##DG1
arch.clean.nolog.dg1 <- subset_samples(arch.clean.nolog, Field == "DG1"|Field == "DG1.18")
arch.clean.nolog.dg1 <- phyloseq::filter_taxa(arch.clean.nolog.dg1, function(x) sum(x) != 0, TRUE)

##DG2
arch.clean.nolog.dg2 <- subset_samples(arch.clean.nolog, Field == "DG2"|Field == "DG2.18")
arch.clean.nolog.dg2 <- phyloseq::filter_taxa(arch.clean.nolog.dg2, function(x) sum(x) != 0, TRUE)

##IS1
arch.clean.nolog.is1 <- subset_samples(arch.clean.nolog, Field == "IS1"|Field == "IS1.18")
arch.clean.nolog.is1 <- phyloseq::filter_taxa(arch.clean.nolog.is1, function(x) sum(x) != 0, TRUE)

##IS2
arch.clean.nolog.is2 <- subset_samples(arch.clean.nolog, Field == "IS2"|Field == "IS2.18")
arch.clean.nolog.is2 <- phyloseq::filter_taxa(arch.clean.nolog.is2, function(x) sum(x) != 0, TRUE)

##Fungi
##CC1
fun.clean.nolog.cc1 <- subset_samples(fun.clean.nolog, Field == "CC1"|Field == "CC1.18")
fun.clean.nolog.cc1 <- phyloseq::filter_taxa(fun.clean.nolog.cc1, function(x) sum(x) != 0, TRUE)

##CC2
fun.clean.nolog.cc2 <- subset_samples(fun.clean.nolog, Field == "CC2"|Field == "CC2.18")
fun.clean.nolog.cc2 <- phyloseq::filter_taxa(fun.clean.nolog.cc2, function(x) sum(x) != 0, TRUE)

##CJ1
fun.clean.nolog.cj1 <- subset_samples(fun.clean.nolog, Field == "CJ1"|Field == "CJ1.18")
fun.clean.nolog.cj1 <- phyloseq::filter_taxa(fun.clean.nolog.cj1, function(x) sum(x) != 0, TRUE)

##CJ2
fun.clean.nolog.cj2 <- subset_samples(fun.clean.nolog, Field == "CJ2"|Field == "CJ2.18")
fun.clean.nolog.cj2 <- phyloseq::filter_taxa(fun.clean.nolog.cj2, function(x) sum(x) != 0, TRUE)

##DG1
fun.clean.nolog.dg1 <- subset_samples(fun.clean.nolog, Field == "DG1"|Field == "DG1.18")
fun.clean.nolog.dg1 <- phyloseq::filter_taxa(fun.clean.nolog.dg1, function(x) sum(x) != 0, TRUE)

##DG2
fun.clean.nolog.dg2 <- subset_samples(fun.clean.nolog, Field == "DG2"|Field == "DG2.18")
fun.clean.nolog.dg2 <- phyloseq::filter_taxa(fun.clean.nolog.dg2, function(x) sum(x) != 0, TRUE)

##IS1
fun.clean.nolog.is1 <- subset_samples(fun.clean.nolog, Field == "IS1"|Field == "IS1.18")
fun.clean.nolog.is1 <- phyloseq::filter_taxa(fun.clean.nolog.is1, function(x) sum(x) != 0, TRUE)

##IS2
fun.clean.nolog.is2 <- subset_samples(fun.clean.nolog, Field == "IS2"|Field == "IS2.18")
fun.clean.nolog.is2 <- phyloseq::filter_taxa(fun.clean.nolog.is2, function(x) sum(x) != 0, TRUE)

## OTU table
#All 
#CC1
otu.all.clean.nolog.cc1<- otu_table(all.clean.nolog.cc1)
df.otu.all.clean.nolog.cc1 <- data.frame(otu.all.clean.nolog.cc1)
df.otu.all.clean.nolog.cc1$OTU <- rownames(df.otu.all.clean.nolog.cc1)
df.otu.all.clean.nolog.cc1 <- merge(df.otu.all.clean.nolog.cc1, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.cc1) <- df.otu.all.clean.nolog.cc1$OTU_id
length(colnames(df.otu.all.clean.nolog.cc1)) #29
df.otu.all.clean.nolog.cc1 <- df.otu.all.clean.nolog.cc1[-c(1, 29)]

write.table(df.otu.all.clean.nolog.cc1, "CC1_all OTU.tsv", sep = '\t', quote = F)

#CC2
otu.all.clean.nolog.cc2<- otu_table(all.clean.nolog.cc2)
df.otu.all.clean.nolog.cc2 <- data.frame(otu.all.clean.nolog.cc2)
df.otu.all.clean.nolog.cc2$OTU <- rownames(df.otu.all.clean.nolog.cc2)
df.otu.all.clean.nolog.cc2 <- merge(df.otu.all.clean.nolog.cc2, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.cc2) <- df.otu.all.clean.nolog.cc2$OTU_id
length(colnames(df.otu.all.clean.nolog.cc2)) #29
df.otu.all.clean.nolog.cc2 <- df.otu.all.clean.nolog.cc2[-c(1, 29)]

write.table(df.otu.all.clean.nolog.cc2, "CC2_all OTU.tsv", sep = '\t', quote = F)

#CJ1
otu.all.clean.nolog.cj1<- otu_table(all.clean.nolog.cj1)
df.otu.all.clean.nolog.cj1 <- data.frame(otu.all.clean.nolog.cj1)
df.otu.all.clean.nolog.cj1$OTU <- rownames(df.otu.all.clean.nolog.cj1)
df.otu.all.clean.nolog.cj1 <- merge(df.otu.all.clean.nolog.cj1, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.cj1) <- df.otu.all.clean.nolog.cj1$OTU_id
length(colnames(df.otu.all.clean.nolog.cj1)) #29
df.otu.all.clean.nolog.cj1 <- df.otu.all.clean.nolog.cj1[-c(1, 29)]

write.table(df.otu.all.clean.nolog.cj1, "CJ1_all OTU.tsv", sep = '\t', quote = F)

#CJ2
otu.all.clean.nolog.cj2<- otu_table(all.clean.nolog.cj2)
df.otu.all.clean.nolog.cj2 <- data.frame(otu.all.clean.nolog.cj2)
df.otu.all.clean.nolog.cj2$OTU <- rownames(df.otu.all.clean.nolog.cj2)
df.otu.all.clean.nolog.cj2 <- merge(df.otu.all.clean.nolog.cj2, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.cj2) <- df.otu.all.clean.nolog.cj2$OTU_id
length(colnames(df.otu.all.clean.nolog.cj2)) #29
df.otu.all.clean.nolog.cj2 <- df.otu.all.clean.nolog.cj2[-c(1, 29)]

write.table(df.otu.all.clean.nolog.cj2, "CJ2_all OTU.tsv", sep = '\t', quote = F)

#DG1
otu.all.clean.nolog.dg1<- otu_table(all.clean.nolog.dg1)
df.otu.all.clean.nolog.dg1 <- data.frame(otu.all.clean.nolog.dg1)
df.otu.all.clean.nolog.dg1$OTU <- rownames(df.otu.all.clean.nolog.dg1)
df.otu.all.clean.nolog.dg1 <- merge(df.otu.all.clean.nolog.dg1, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.dg1) <- df.otu.all.clean.nolog.dg1$OTU_id
length(colnames(df.otu.all.clean.nolog.dg1)) #29
df.otu.all.clean.nolog.dg1 <- df.otu.all.clean.nolog.dg1[-c(1, 29)]

write.table(df.otu.all.clean.nolog.dg1, "DG1_all OTU.tsv", sep = '\t', quote = F)

#DG2
otu.all.clean.nolog.dg2<- otu_table(all.clean.nolog.dg2)
df.otu.all.clean.nolog.dg2 <- data.frame(otu.all.clean.nolog.dg2)
df.otu.all.clean.nolog.dg2$OTU <- rownames(df.otu.all.clean.nolog.dg2)
df.otu.all.clean.nolog.dg2 <- merge(df.otu.all.clean.nolog.dg2, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.dg2) <- df.otu.all.clean.nolog.dg2$OTU_id
length(colnames(df.otu.all.clean.nolog.dg2)) #29
df.otu.all.clean.nolog.dg2 <- df.otu.all.clean.nolog.dg2[-c(1, 29)]

write.table(df.otu.all.clean.nolog.dg2, "DG2_all OTU.tsv", sep = '\t', quote = F)

#IS1
otu.all.clean.nolog.is1<- otu_table(all.clean.nolog.is1)
df.otu.all.clean.nolog.is1 <- data.frame(otu.all.clean.nolog.is1)
df.otu.all.clean.nolog.is1$OTU <- rownames(df.otu.all.clean.nolog.is1)
df.otu.all.clean.nolog.is1 <- merge(df.otu.all.clean.nolog.is1, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.is1) <- df.otu.all.clean.nolog.is1$OTU_id
length(colnames(df.otu.all.clean.nolog.is1)) #29
df.otu.all.clean.nolog.is1 <- df.otu.all.clean.nolog.is1[-c(1, 29)]

write.table(df.otu.all.clean.nolog.is1, "IS1_all OTU.tsv", sep = '\t', quote = F)

#IS2
otu.all.clean.nolog.is2<- otu_table(all.clean.nolog.is2)
df.otu.all.clean.nolog.is2 <- data.frame(otu.all.clean.nolog.is2)
df.otu.all.clean.nolog.is2$OTU <- rownames(df.otu.all.clean.nolog.is2)
df.otu.all.clean.nolog.is2 <- merge(df.otu.all.clean.nolog.is2, OTU_id.list, by = 'OTU')
rownames(df.otu.all.clean.nolog.is2) <- df.otu.all.clean.nolog.is2$OTU_id
length(colnames(df.otu.all.clean.nolog.is2)) #29
df.otu.all.clean.nolog.is2 <- df.otu.all.clean.nolog.is2[-c(1, 29)]

write.table(df.otu.all.clean.nolog.is2, "IS2_all OTU.tsv", sep = '\t', quote = F)

#Bacteria only 
#CC1
otu.bac.clean.nolog.cc1<- otu_table(bac.clean.nolog.cc1)
df.otu.bac.clean.nolog.cc1 <- data.frame(otu.bac.clean.nolog.cc1)
df.otu.bac.clean.nolog.cc1$OTU <- rownames(df.otu.bac.clean.nolog.cc1)
df.otu.bac.clean.nolog.cc1 <- merge(df.otu.bac.clean.nolog.cc1, OTU_id.list, by = 'OTU')
rownames(df.otu.bac.clean.nolog.cc1) <- df.otu.bac.clean.nolog.cc1$OTU_id
length(colnames(df.otu.bac.clean.nolog.cc1)) #20
df.otu.bac.clean.nolog.cc1 <- df.otu.bac.clean.nolog.cc1[-c(1, 20)]

#CC2
otu.bac.clean.nolog.cc2<- otu_table(bac.clean.nolog.cc2)
df.otu.bac.clean.nolog.cc2 <- data.frame(otu.bac.clean.nolog.cc2)
df.otu.bac.clean.nolog.cc2$OTU <- rownames(df.otu.bac.clean.nolog.cc2)
df.otu.bac.clean.nolog.cc2 <- merge(df.otu.bac.clean.nolog.cc2, OTU_id.list, by = 'OTU')
rownames(df.otu.bac.clean.nolog.cc2) <- df.otu.bac.clean.nolog.cc2$OTU_id
length(colnames(df.otu.bac.clean.nolog.cc2)) #20
df.otu.bac.clean.nolog.cc2 <- df.otu.bac.clean.nolog.cc2[-c(1, 20)]

#CJ1
otu.bac.clean.nolog.cj1<- otu_table(bac.clean.nolog.cj1)
df.otu.bac.clean.nolog.cj1 <- data.frame(otu.bac.clean.nolog.cj1)
df.otu.bac.clean.nolog.cj1$OTU <- rownames(df.otu.bac.clean.nolog.cj1)
df.otu.bac.clean.nolog.cj1 <- merge(df.otu.bac.clean.nolog.cj1, OTU_id.list, by = 'OTU')
rownames(df.otu.bac.clean.nolog.cj1) <- df.otu.bac.clean.nolog.cj1$OTU_id
length(colnames(df.otu.bac.clean.nolog.cj1)) #20
df.otu.bac.clean.nolog.cj1 <- df.otu.bac.clean.nolog.cj1[-c(1, 20)]

#CJ2
otu.bac.clean.nolog.cj2<- otu_table(bac.clean.nolog.cj2)
df.otu.bac.clean.nolog.cj2 <- data.frame(otu.bac.clean.nolog.cj2)
df.otu.bac.clean.nolog.cj2$OTU <- rownames(df.otu.bac.clean.nolog.cj2)
df.otu.bac.clean.nolog.cj2 <- merge(df.otu.bac.clean.nolog.cj2, OTU_id.list, by = 'OTU')
rownames(df.otu.bac.clean.nolog.cj2) <- df.otu.bac.clean.nolog.cj2$OTU_id
length(colnames(df.otu.bac.clean.nolog.cj2)) #20
df.otu.bac.clean.nolog.cj2 <- df.otu.bac.clean.nolog.cj2[-c(1, 20)]

#DG1
otu.bac.clean.nolog.dg1<- otu_table(bac.clean.nolog.dg1)
df.otu.bac.clean.nolog.dg1 <- data.frame(otu.bac.clean.nolog.dg1)
df.otu.bac.clean.nolog.dg1$OTU <- rownames(df.otu.bac.clean.nolog.dg1)
df.otu.bac.clean.nolog.dg1 <- merge(df.otu.bac.clean.nolog.dg1, OTU_id.list, by = 'OTU')
rownames(df.otu.bac.clean.nolog.dg1) <- df.otu.bac.clean.nolog.dg1$OTU_id
length(colnames(df.otu.bac.clean.nolog.dg1)) #20
df.otu.bac.clean.nolog.dg1 <- df.otu.bac.clean.nolog.dg1[-c(1, 20)]

#DG2
otu.bac.clean.nolog.dg2<- otu_table(bac.clean.nolog.dg2)
df.otu.bac.clean.nolog.dg2 <- data.frame(otu.bac.clean.nolog.dg2)
df.otu.bac.clean.nolog.dg2$OTU <- rownames(df.otu.bac.clean.nolog.dg2)
df.otu.bac.clean.nolog.dg2 <- merge(df.otu.bac.clean.nolog.dg2, OTU_id.list, by = 'OTU')
rownames(df.otu.bac.clean.nolog.dg2) <- df.otu.bac.clean.nolog.dg2$OTU_id
length(colnames(df.otu.bac.clean.nolog.dg2)) #20
df.otu.bac.clean.nolog.dg2 <- df.otu.bac.clean.nolog.dg2[-c(1, 20)]

#IS1
otu.bac.clean.nolog.is1<- otu_table(bac.clean.nolog.is1)
df.otu.bac.clean.nolog.is1 <- data.frame(otu.bac.clean.nolog.is1)
df.otu.bac.clean.nolog.is1$OTU <- rownames(df.otu.bac.clean.nolog.is1)
df.otu.bac.clean.nolog.is1 <- merge(df.otu.bac.clean.nolog.is1, OTU_id.list, by = 'OTU')
rownames(df.otu.bac.clean.nolog.is1) <- df.otu.bac.clean.nolog.is1$OTU_id
length(colnames(df.otu.bac.clean.nolog.is1)) #20
df.otu.bac.clean.nolog.is1 <- df.otu.bac.clean.nolog.is1[-c(1, 20)]

#IS2
otu.bac.clean.nolog.is2<- otu_table(bac.clean.nolog.is2)
df.otu.bac.clean.nolog.is2 <- data.frame(otu.bac.clean.nolog.is2)
df.otu.bac.clean.nolog.is2$OTU <- rownames(df.otu.bac.clean.nolog.is2)
df.otu.bac.clean.nolog.is2 <- merge(df.otu.bac.clean.nolog.is2, OTU_id.list, by = 'OTU')
rownames(df.otu.bac.clean.nolog.is2) <- df.otu.bac.clean.nolog.is2$OTU_id
length(colnames(df.otu.bac.clean.nolog.is2)) #20
df.otu.bac.clean.nolog.is2 <- df.otu.bac.clean.nolog.is2[-c(1, 20)]


#Archaea only 
#CC1
otu.arch.clean.nolog.cc1<- otu_table(arch.clean.nolog.cc1)
df.otu.arch.clean.nolog.cc1 <- data.frame(otu.arch.clean.nolog.cc1)
df.otu.arch.clean.nolog.cc1$OTU <- rownames(df.otu.arch.clean.nolog.cc1)
df.otu.arch.clean.nolog.cc1 <- merge(df.otu.arch.clean.nolog.cc1, OTU_id.list, by = 'OTU')
rownames(df.otu.arch.clean.nolog.cc1) <- df.otu.arch.clean.nolog.cc1$OTU_id
length(colnames(df.otu.arch.clean.nolog.cc1)) #20
df.otu.arch.clean.nolog.cc1 <- df.otu.arch.clean.nolog.cc1[-c(1, 20)]

#CC2
otu.arch.clean.nolog.cc2<- otu_table(arch.clean.nolog.cc2)
df.otu.arch.clean.nolog.cc2 <- data.frame(otu.arch.clean.nolog.cc2)
df.otu.arch.clean.nolog.cc2$OTU <- rownames(df.otu.arch.clean.nolog.cc2)
df.otu.arch.clean.nolog.cc2 <- merge(df.otu.arch.clean.nolog.cc2, OTU_id.list, by = 'OTU')
rownames(df.otu.arch.clean.nolog.cc2) <- df.otu.arch.clean.nolog.cc2$OTU_id
length(colnames(df.otu.arch.clean.nolog.cc2)) #20
df.otu.arch.clean.nolog.cc2 <- df.otu.arch.clean.nolog.cc2[-c(1, 20)]

#CJ1
otu.arch.clean.nolog.cj1<- otu_table(arch.clean.nolog.cj1)
df.otu.arch.clean.nolog.cj1 <- data.frame(otu.arch.clean.nolog.cj1)
df.otu.arch.clean.nolog.cj1$OTU <- rownames(df.otu.arch.clean.nolog.cj1)
df.otu.arch.clean.nolog.cj1 <- merge(df.otu.arch.clean.nolog.cj1, OTU_id.list, by = 'OTU')
rownames(df.otu.arch.clean.nolog.cj1) <- df.otu.arch.clean.nolog.cj1$OTU_id
length(colnames(df.otu.arch.clean.nolog.cj1)) #20
df.otu.arch.clean.nolog.cj1 <- df.otu.arch.clean.nolog.cj1[-c(1, 20)]

#CJ2
otu.arch.clean.nolog.cj2<- otu_table(arch.clean.nolog.cj2)
df.otu.arch.clean.nolog.cj2 <- data.frame(otu.arch.clean.nolog.cj2)
df.otu.arch.clean.nolog.cj2$OTU <- rownames(df.otu.arch.clean.nolog.cj2)
df.otu.arch.clean.nolog.cj2 <- merge(df.otu.arch.clean.nolog.cj2, OTU_id.list, by = 'OTU')
rownames(df.otu.arch.clean.nolog.cj2) <- df.otu.arch.clean.nolog.cj2$OTU_id
length(colnames(df.otu.arch.clean.nolog.cj2)) #20
df.otu.arch.clean.nolog.cj2 <- df.otu.arch.clean.nolog.cj2[-c(1, 20)]

#DG1
otu.arch.clean.nolog.dg1<- otu_table(arch.clean.nolog.dg1)
df.otu.arch.clean.nolog.dg1 <- data.frame(otu.arch.clean.nolog.dg1)
df.otu.arch.clean.nolog.dg1$OTU <- rownames(df.otu.arch.clean.nolog.dg1)
df.otu.arch.clean.nolog.dg1 <- merge(df.otu.arch.clean.nolog.dg1, OTU_id.list, by = 'OTU')
rownames(df.otu.arch.clean.nolog.dg1) <- df.otu.arch.clean.nolog.dg1$OTU_id
length(colnames(df.otu.arch.clean.nolog.dg1)) #20
df.otu.arch.clean.nolog.dg1 <- df.otu.arch.clean.nolog.dg1[-c(1, 20)]

#DG2
otu.arch.clean.nolog.dg2<- otu_table(arch.clean.nolog.dg2)
df.otu.arch.clean.nolog.dg2 <- data.frame(otu.arch.clean.nolog.dg2)
df.otu.arch.clean.nolog.dg2$OTU <- rownames(df.otu.arch.clean.nolog.dg2)
df.otu.arch.clean.nolog.dg2 <- merge(df.otu.arch.clean.nolog.dg2, OTU_id.list, by = 'OTU')
rownames(df.otu.arch.clean.nolog.dg2) <- df.otu.arch.clean.nolog.dg2$OTU_id
length(colnames(df.otu.arch.clean.nolog.dg2)) #20
df.otu.arch.clean.nolog.dg2 <- df.otu.arch.clean.nolog.dg2[-c(1, 20)]

#IS1
otu.arch.clean.nolog.is1<- otu_table(arch.clean.nolog.is1)
df.otu.arch.clean.nolog.is1 <- data.frame(otu.arch.clean.nolog.is1)
df.otu.arch.clean.nolog.is1$OTU <- rownames(df.otu.arch.clean.nolog.is1)
df.otu.arch.clean.nolog.is1 <- merge(df.otu.arch.clean.nolog.is1, OTU_id.list, by = 'OTU')
rownames(df.otu.arch.clean.nolog.is1) <- df.otu.arch.clean.nolog.is1$OTU_id
length(colnames(df.otu.arch.clean.nolog.is1)) #20
df.otu.arch.clean.nolog.is1 <- df.otu.arch.clean.nolog.is1[-c(1, 20)]

#IS2
otu.arch.clean.nolog.is2<- otu_table(arch.clean.nolog.is2)
df.otu.arch.clean.nolog.is2 <- data.frame(otu.arch.clean.nolog.is2)
df.otu.arch.clean.nolog.is2$OTU <- rownames(df.otu.arch.clean.nolog.is2)
df.otu.arch.clean.nolog.is2 <- merge(df.otu.arch.clean.nolog.is2, OTU_id.list, by = 'OTU')
rownames(df.otu.arch.clean.nolog.is2) <- df.otu.arch.clean.nolog.is2$OTU_id
length(colnames(df.otu.arch.clean.nolog.is2)) #20
df.otu.arch.clean.nolog.is2 <- df.otu.arch.clean.nolog.is2[-c(1, 20)]

#Fungi
#CC1
otu.fun.clean.nolog.cc1<- otu_table(fun.clean.nolog.cc1)
df.otu.fun.clean.nolog.cc1 <- data.frame(otu.fun.clean.nolog.cc1)
df.otu.fun.clean.nolog.cc1$OTU <- rownames(df.otu.fun.clean.nolog.cc1)
df.otu.fun.clean.nolog.cc1 <- merge(df.otu.fun.clean.nolog.cc1, OTU_id.list, by = 'OTU')
rownames(df.otu.fun.clean.nolog.cc1) <- df.otu.fun.clean.nolog.cc1$OTU_id
length(colnames(df.otu.fun.clean.nolog.cc1)) #20
df.otu.fun.clean.nolog.cc1 <- df.otu.fun.clean.nolog.cc1[-c(1, 20)]

#CC2
otu.fun.clean.nolog.cc2<- otu_table(fun.clean.nolog.cc2)
df.otu.fun.clean.nolog.cc2 <- data.frame(otu.fun.clean.nolog.cc2)
df.otu.fun.clean.nolog.cc2$OTU <- rownames(df.otu.fun.clean.nolog.cc2)
df.otu.fun.clean.nolog.cc2 <- merge(df.otu.fun.clean.nolog.cc2, OTU_id.list, by = 'OTU')
rownames(df.otu.fun.clean.nolog.cc2) <- df.otu.fun.clean.nolog.cc2$OTU_id
length(colnames(df.otu.fun.clean.nolog.cc2)) #20
df.otu.fun.clean.nolog.cc2 <- df.otu.fun.clean.nolog.cc2[-c(1, 20)]

#CJ1
otu.fun.clean.nolog.cj1<- otu_table(fun.clean.nolog.cj1)
df.otu.fun.clean.nolog.cj1 <- data.frame(otu.fun.clean.nolog.cj1)
df.otu.fun.clean.nolog.cj1$OTU <- rownames(df.otu.fun.clean.nolog.cj1)
df.otu.fun.clean.nolog.cj1 <- merge(df.otu.fun.clean.nolog.cj1, OTU_id.list, by = 'OTU')
rownames(df.otu.fun.clean.nolog.cj1) <- df.otu.fun.clean.nolog.cj1$OTU_id
length(colnames(df.otu.fun.clean.nolog.cj1)) #20
df.otu.fun.clean.nolog.cj1 <- df.otu.fun.clean.nolog.cj1[-c(1, 20)]

#CJ2
otu.fun.clean.nolog.cj2<- otu_table(fun.clean.nolog.cj2)
df.otu.fun.clean.nolog.cj2 <- data.frame(otu.fun.clean.nolog.cj2)
df.otu.fun.clean.nolog.cj2$OTU <- rownames(df.otu.fun.clean.nolog.cj2)
df.otu.fun.clean.nolog.cj2 <- merge(df.otu.fun.clean.nolog.cj2, OTU_id.list, by = 'OTU')
rownames(df.otu.fun.clean.nolog.cj2) <- df.otu.fun.clean.nolog.cj2$OTU_id
length(colnames(df.otu.fun.clean.nolog.cj2)) #20
df.otu.fun.clean.nolog.cj2 <- df.otu.fun.clean.nolog.cj2[-c(1, 20)]

#DG1
otu.fun.clean.nolog.dg1<- otu_table(fun.clean.nolog.dg1)
df.otu.fun.clean.nolog.dg1 <- data.frame(otu.fun.clean.nolog.dg1)
df.otu.fun.clean.nolog.dg1$OTU <- rownames(df.otu.fun.clean.nolog.dg1)
df.otu.fun.clean.nolog.dg1 <- merge(df.otu.fun.clean.nolog.dg1, OTU_id.list, by = 'OTU')
rownames(df.otu.fun.clean.nolog.dg1) <- df.otu.fun.clean.nolog.dg1$OTU_id
length(colnames(df.otu.fun.clean.nolog.dg1)) #20
df.otu.fun.clean.nolog.dg1 <- df.otu.fun.clean.nolog.dg1[-c(1, 20)]

#DG2
otu.fun.clean.nolog.dg2<- otu_table(fun.clean.nolog.dg2)
df.otu.fun.clean.nolog.dg2 <- data.frame(otu.fun.clean.nolog.dg2)
df.otu.fun.clean.nolog.dg2$OTU <- rownames(df.otu.fun.clean.nolog.dg2)
df.otu.fun.clean.nolog.dg2 <- merge(df.otu.fun.clean.nolog.dg2, OTU_id.list, by = 'OTU')
rownames(df.otu.fun.clean.nolog.dg2) <- df.otu.fun.clean.nolog.dg2$OTU_id
length(colnames(df.otu.fun.clean.nolog.dg2)) #20
df.otu.fun.clean.nolog.dg2 <- df.otu.fun.clean.nolog.dg2[-c(1, 20)]

#IS1
otu.fun.clean.nolog.is1<- otu_table(fun.clean.nolog.is1)
df.otu.fun.clean.nolog.is1 <- data.frame(otu.fun.clean.nolog.is1)
df.otu.fun.clean.nolog.is1$OTU <- rownames(df.otu.fun.clean.nolog.is1)
df.otu.fun.clean.nolog.is1 <- merge(df.otu.fun.clean.nolog.is1, OTU_id.list, by = 'OTU')
rownames(df.otu.fun.clean.nolog.is1) <- df.otu.fun.clean.nolog.is1$OTU_id
length(colnames(df.otu.fun.clean.nolog.is1)) #20
df.otu.fun.clean.nolog.is1 <- df.otu.fun.clean.nolog.is1[-c(1, 20)]

#IS2
otu.fun.clean.nolog.is2<- otu_table(fun.clean.nolog.is2)
df.otu.fun.clean.nolog.is2 <- data.frame(otu.fun.clean.nolog.is2)
df.otu.fun.clean.nolog.is2$OTU <- rownames(df.otu.fun.clean.nolog.is2)
df.otu.fun.clean.nolog.is2 <- merge(df.otu.fun.clean.nolog.is2, OTU_id.list, by = 'OTU')
rownames(df.otu.fun.clean.nolog.is2) <- df.otu.fun.clean.nolog.is2$OTU_id
length(colnames(df.otu.fun.clean.nolog.is2)) #20
df.otu.fun.clean.nolog.is2 <- df.otu.fun.clean.nolog.is2[-c(1, 20)]


## Correlation
otu_norm_soil_combine.cc1 <- read.table('CC1_all OTU_edit.tsv', sep = '\t', header =T)
rownames(otu_norm_soil_combine.cc1) <- otu_norm_soil_combine.cc1$X
otu_norm_soil_combine.cc1 <- otu_norm_soil_combine.cc1[-c(1)]

otu_norm_soil_combine.cc2 <- read.table('CC2_all OTU_edit.tsv', sep = '\t', header =T)
rownames(otu_norm_soil_combine.cc2) <- otu_norm_soil_combine.cc2$X
otu_norm_soil_combine.cc2 <- otu_norm_soil_combine.cc2[-c(1)]

otu_norm_soil_combine.cj1 <- read.table('CJ1_all OTU_edit.tsv', sep = '\t', header =T)
rownames(otu_norm_soil_combine.cj1) <- otu_norm_soil_combine.cj1$X
otu_norm_soil_combine.cj1 <- otu_norm_soil_combine.cj1[-c(1)]

otu_norm_soil_combine.cj2 <- read.table('CJ2_all OTU_edit.tsv', sep = '\t', header =T)
rownames(otu_norm_soil_combine.cj2) <- otu_norm_soil_combine.cj2$X
otu_norm_soil_combine.cj2 <- otu_norm_soil_combine.cj2[-c(1)]

otu_norm_soil_combine.dg1 <- read.table('DG1_all OTU_edit.tsv', sep = '\t', header =T)
rownames(otu_norm_soil_combine.dg1) <- otu_norm_soil_combine.dg1$X
otu_norm_soil_combine.dg1 <- otu_norm_soil_combine.dg1[-c(1)]

otu_norm_soil_combine.dg2 <- read.table('DG2_all OTU_edit.tsv', sep = '\t', header =T)
rownames(otu_norm_soil_combine.dg2) <- otu_norm_soil_combine.dg2$X
otu_norm_soil_combine.dg2 <- otu_norm_soil_combine.dg2[-c(1)]

otu_norm_soil_combine.is1 <- read.table('IS1_all OTU_edit.tsv', sep = '\t', header =T)
rownames(otu_norm_soil_combine.is1) <- otu_norm_soil_combine.is1$X
otu_norm_soil_combine.is1 <- otu_norm_soil_combine.is1[-c(1)]

otu_norm_soil_combine.is2 <- read.table('IS2_all OTU_edit.tsv', sep = '\t', header =T)
rownames(otu_norm_soil_combine.is2) <- otu_norm_soil_combine.is2$X
otu_norm_soil_combine.is2 <- otu_norm_soil_combine.is2[-c(1)]


all_soil_cor.cc1 <- rcorr(t(otu_norm_soil_combine.cc1), type=c("spearman"))
all_soil_cor.cc2 <- rcorr(t(otu_norm_soil_combine.cc2), type=c("spearman"))
all_soil_cor.cj1 <- rcorr(t(otu_norm_soil_combine.cj1), type=c("spearman"))
all_soil_cor.cj2 <- rcorr(t(otu_norm_soil_combine.cj2), type=c("spearman"))
all_soil_cor.dg1 <- rcorr(t(otu_norm_soil_combine.dg1), type=c("spearman"))
all_soil_cor.dg2 <- rcorr(t(otu_norm_soil_combine.dg2), type=c("spearman"))
all_soil_cor.is1 <- rcorr(t(otu_norm_soil_combine.is1), type=c("spearman"))
all_soil_cor.is2 <- rcorr(t(otu_norm_soil_combine.is2), type=c("spearman"))

##CC
all_cor_soil_df.cc1 <- CorrDF(all_soil_cor.cc1$r, all_soil_cor.cc1$P)
all_cor_soil_df.cc1$padj <- p.adjust(all_cor_soil_df.cc1$p, method="none")
all_cor_soil_df_padj.cc1 <- all_cor_soil_df.cc1[which(abs(all_cor_soil_df.cc1$cor) > 0.7),]
all_cor_soil_df_padj.cc1 <- all_cor_soil_df_padj.cc1[which(all_cor_soil_df_padj.cc1$padj < 0.001),]

write.table(all_cor_soil_df_padj.cc1, "all_cor_soil_df_padj.cc1.txt", sep='\t', quote = F)


all_cor_soil_df.cc2 <- CorrDF(all_soil_cor.cc2$r, all_soil_cor.cc2$P)
all_cor_soil_df.cc2$padj <- p.adjust(all_cor_soil_df.cc2$p, method="none")
all_cor_soil_df_padj.cc2 <- all_cor_soil_df.cc2[which(abs(all_cor_soil_df.cc2$cor) > 0.7),]
all_cor_soil_df_padj.cc2 <- all_cor_soil_df_padj.cc2[which(all_cor_soil_df_padj.cc2$padj < 0.001),]

write.table(all_cor_soil_df_padj.cc2, "all_cor_soil_df_padj.cc2.txt", sep='\t', quote = F)

#CJ
all_cor_soil_df.cj1 <- CorrDF(all_soil_cor.cj1$r, all_soil_cor.cj1$P)
all_cor_soil_df.cj1$padj <- p.adjust(all_cor_soil_df.cj1$p, method="none")
all_cor_soil_df_padj.cj1 <- all_cor_soil_df.cj1[which(abs(all_cor_soil_df.cj1$cor) > 0.7),]
all_cor_soil_df_padj.cj1 <- all_cor_soil_df_padj.cj1[which(all_cor_soil_df_padj.cj1$padj < 0.001),]

write.table(all_cor_soil_df_padj.cj1, "all_cor_soil_df_padj.cj1.txt", sep='\t', quote = F)


all_cor_soil_df.cj2 <- CorrDF(all_soil_cor.cj2$r, all_soil_cor.cj2$P)
all_cor_soil_df.cj2$padj <- p.adjust(all_cor_soil_df.cj2$p, method="none")
all_cor_soil_df_padj.cj2 <- all_cor_soil_df.cj2[which(abs(all_cor_soil_df.cj2$cor) > 0.7),]
all_cor_soil_df_padj.cj2 <- all_cor_soil_df_padj.cj2[which(all_cor_soil_df_padj.cj2$padj < 0.001),]

write.table(all_cor_soil_df_padj.cj2, "all_cor_soil_df_padj.cj2.txt", sep='\t', quote = F)


#DG
all_cor_soil_df.dg1 <- CorrDF(all_soil_cor.dg1$r, all_soil_cor.dg1$P)
all_cor_soil_df.dg1$padj <- p.adjust(all_cor_soil_df.dg1$p, method="none")
all_cor_soil_df_padj.dg1 <- all_cor_soil_df.dg1[which(abs(all_cor_soil_df.dg1$cor) > 0.7),]
all_cor_soil_df_padj.dg1 <- all_cor_soil_df_padj.dg1[which(all_cor_soil_df_padj.dg1$padj < 0.001),]

write.table(all_cor_soil_df_padj.dg1, "all_cor_soil_df_padj.dg1.txt", sep='\t', quote = F)


all_cor_soil_df.dg2 <- CorrDF(all_soil_cor.dg2$r, all_soil_cor.dg2$P)
all_cor_soil_df.dg2$padj <- p.adjust(all_cor_soil_df.dg2$p, method="none")
all_cor_soil_df_padj.dg2 <- all_cor_soil_df.dg2[which(abs(all_cor_soil_df.dg2$cor) > 0.7),]
all_cor_soil_df_padj.dg2 <- all_cor_soil_df_padj.dg2[which(all_cor_soil_df_padj.dg2$padj < 0.001),]

write.table(all_cor_soil_df_padj.dg2, "all_cor_soil_df_padj.dg2.txt", sep='\t', quote = F)

#IS
all_cor_soil_df.is1 <- CorrDF(all_soil_cor.is1$r, all_soil_cor.is1$P)
all_cor_soil_df.is1$padj <- p.adjust(all_cor_soil_df.is1$p, method="none")
all_cor_soil_df_padj.is1 <- all_cor_soil_df.is1[which(abs(all_cor_soil_df.is1$cor) > 0.7),]
all_cor_soil_df_padj.is1 <- all_cor_soil_df_padj.is1[which(all_cor_soil_df_padj.is1$padj < 0.001),]

write.table(all_cor_soil_df_padj.is1, "all_cor_soil_df_padj.is1.txt", sep='\t', quote = F)


all_cor_soil_df.is2 <- CorrDF(all_soil_cor.is2$r, all_soil_cor.is2$P)
all_cor_soil_df.is2$padj <- p.adjust(all_cor_soil_df.is2$p, method="none")
all_cor_soil_df_padj.is2 <- all_cor_soil_df.is2[which(abs(all_cor_soil_df.is2$cor) > 0.7),]
all_cor_soil_df_padj.is2 <- all_cor_soil_df_padj.is2[which(all_cor_soil_df_padj.is2$padj < 0.001),]

write.table(all_cor_soil_df_padj.is2, "all_cor_soil_df_padj.is2.txt", sep='\t', quote = F)


#igraph
## CC1
all_cor_soil_df_padj.cc1 <- read.table("all_cor_soil_df_padj.cc1.txt", sep='\t', header = T)


nodeattrib_soil_combine.cc1 <- data.frame(node=union(all_cor_soil_df_padj.cc1$from,all_cor_soil_df_padj.cc1$to))
nodeattrib_soil_combine.cc1$kingdom <- 0

for (i in as.character(nodeattrib_soil_combine.cc1$node))
{
  if (i %in% rownames(df.otu.bac.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.cc1[nodeattrib_soil_combine.cc1$node==i,"kingdom"] <- "Bacteria"}
  
  else if (i %in% rownames(df.otu.arch.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.cc1[nodeattrib_soil_combine.cc1$node==i,"kingdom"] <- "Archaea"} 
  
  else
  {nodeattrib_soil_combine.cc1[nodeattrib_soil_combine.cc1$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_soil_combine.cc1) <- as.character(nodeattrib_soil_combine.cc1$node)
nodeattrib_soil_combine.cc1$kingdom



all_soil_net.cc1 <- graph_from_data_frame(all_cor_soil_df_padj.cc1,direct=F, vertices=nodeattrib_soil_combine.cc1)

## Number of nodes
length(V(all_soil_net.cc1)) #6030

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.cc1)))) #4374
length(grep("^F",names(V(all_soil_net.cc1)))) #1475
length(grep("^A",names(V(all_soil_net.cc1)))) #181

## Connections 
bb_occur_soil.cc1 <- droplevels(all_cor_soil_df_padj.cc1[with(all_cor_soil_df_padj.cc1, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.cc1) #380726

ff_occur_soil.cc1 <- droplevels(all_cor_soil_df_padj.cc1[with(all_cor_soil_df_padj.cc1, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.cc1) #45956

fb_occur_soil.cc1 <- droplevels(all_cor_soil_df_padj.cc1[with(all_cor_soil_df_padj.cc1, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.cc1) #209020

aa_occur_soil.cc1 <- droplevels(all_cor_soil_df_padj.cc1[with(all_cor_soil_df_padj.cc1, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.cc1) #568

ba_occur_soil.cc1 <- droplevels(all_cor_soil_df_padj.cc1[with(all_cor_soil_df_padj.cc1, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.cc1) #28844

fa_occur_soil.cc1 <- droplevels(all_cor_soil_df_padj.cc1[with(all_cor_soil_df_padj.cc1, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.cc1) #7781



## Network properties
meta_degree.cc1 <- sort(igraph::degree(all_soil_net.cc1,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.cc1))) #"max degree = " "399"     
print(c('mean degree = ', mean(meta_degree.cc1))) #"mean degree = "   "223.182421227197"
meta_degree.cc1 <- as.data.frame(meta_degree.cc1)
ggplot(meta_degree.cc1, aes(x=meta_degree.cc1)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.cc1, directed = FALSE))) #"average shortest path length =  3.59523728700828"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.cc1, "global"))) #"mean clustering coefficient =  0.949977913375762"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.cc1, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.cc1, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.cc1, V(all_soil_net.cc1)))))) #"mean number of neighbors =  223.182421227197"

##
net <- all_soil_net.cc1
CC1_all_deg <- igraph::degree(net,mode="all")
CC1_all_betweenness <- betweenness(net, normalized = TRUE)
CC1_all_closeness <- closeness(net, normalized = TRUE)
CC1_all_transitivity <- transitivity(net, "local", vids = V(net))
names(CC1_all_transitivity)<- V(net)$name
CC1_all_transitivity[is.na(CC1_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
CC1_all_deg.1percent <- CC1_all_deg[CC1_all_deg >= quantile(CC1_all_deg,prob=1-n/100)]
length(CC1_all_deg.1percent) #326

CC1_all_betweenness.1percent <- CC1_all_betweenness[CC1_all_betweenness >= quantile(CC1_all_betweenness,prob=1-n/100)]
length(CC1_all_betweenness.1percent) #61

CC1_all_closeness.1percent <- CC1_all_closeness[CC1_all_closeness >= quantile(CC1_all_closeness,prob=1-n/100)]
length(CC1_all_closeness.1percent) #292

intersect(names(CC1_all_deg.1percent), names(CC1_all_betweenness.1percent))
intersect(names(CC1_all_deg.1percent), names(CC1_all_closeness.1percent))
intersect(names(CC1_all_betweenness.1percent), names(CC1_all_closeness.1percent))


## Set node shape soil_all_keystone.2017
V(all_soil_net.cc1)$shape <- V(all_soil_net.cc1)$name
V(all_soil_net.cc1)$shape[V(all_soil_net.cc1)$shape %in% names(CC1_all_betweenness.1percent)] <- "star"
#V(all_soil_net.cc1)$shape[V(all_soil_net.cc1)$shape %in% strictkeystone.cc1] <- "star"
V(all_soil_net.cc1)$shape[V(all_soil_net.cc1)$shape %in% rownames(df.otu.bac.clean.nolog.cc1)] <- "circle"
V(all_soil_net.cc1)$shape[V(all_soil_net.cc1)$shape %in% rownames(df.otu.fun.clean.nolog.cc1)] <- "triangle"
V(all_soil_net.cc1)$shape[V(all_soil_net.cc1)$shape %in% rownames(df.otu.arch.clean.nolog.cc1)] <- "square"

## Set node colors based kingdom
cs <- c("Bacteria", "Archaea", "Fungi")
unique(V(all_soil_net.cc1)$kingdom)
V(all_soil_net.cc1)$color <- V(all_soil_net.cc1)$kingdom
V(all_soil_net.cc1)$color[V(all_soil_net.cc1)$color == "Bacteria"] <- "#EE7600"
V(all_soil_net.cc1)$color[V(all_soil_net.cc1)$color == "Archaea"] <- "#458B74"
V(all_soil_net.cc1)$color[V(all_soil_net.cc1)$color == "Fungi"] <- "#9A32CD"
V(all_soil_net.cc1)$frame.color <- V(all_soil_net.cc1)$color

#Node size
V(all_soil_net.cc1)$size <- V(all_soil_net.cc1)$name
V(all_soil_net.cc1)$size[V(all_soil_net.cc1)$size %in% names(CC1_all_betweenness.1percent)] <- 3
#V(all_soil_net.cc1)$size[V(all_soil_net.cc1)$size %in% strictkeystone.2017] <- 5
V(all_soil_net.cc1)$size[V(all_soil_net.cc1)$size %in% rownames(df.otu.bac.clean.nolog.cc1)] <- 2
V(all_soil_net.cc1)$size[V(all_soil_net.cc1)$size %in% rownames(df.otu.fun.clean.nolog.cc1)] <- 3
V(all_soil_net.cc1)$size[V(all_soil_net.cc1)$size %in% rownames(df.otu.arch.clean.nolog.cc1)] <- 2
soilcombine_nodesizes.cc1 <- as.numeric(V(all_soil_net.cc1)$size)

## Set edge color
E(all_soil_net.cc1)$color <- ifelse(E(all_soil_net.cc1)$cor >0, "#99CCFF","#FF6666")


#Layout
source("star_shape.R")
source("triangle_shape.R")
library(qgraph)
all_soil_net.cc1.integer <- set.vertex.attribute(all_soil_net.cc1, "name", value=paste(1:6030))
edge.cc1 <- get.edgelist(all_soil_net.cc1.integer)
edge.cc1.num<- apply(edge.cc1, 2, as.numeric)

layout.cc1 <- qgraph.layout.fruchtermanreingold(edge.cc1.num,vcount=vcount(all_soil_net.cc1.integer))

#Default
plot(all_soil_net.cc1.integer,layout=layout.cc1,vertex.label=NA, vertex.size=soilcombine_nodesizes.cc1)

mtext("qgraph.layout.fruchtermanreingold default", side=1)

#Modified
layout.cc1.modified <- qgraph.layout.fruchtermanreingold(edge.cc1.num,vcount=vcount(all_soil_net.cc1.integer),
                                                        area=10*(vcount(all_soil_net.cc1.integer)^2),repulse.rad=(vcount(all_soil_net.cc1.integer)^3.1))

#without module
plot(all_soil_net.cc1.integer,layout=layout.cc1.modified,vertex.size=soilcombine_nodesizes.cc1,vertex.label=NA)

memory.limit(200000)

dev.off()

## CC2
all_cor_soil_df_padj.cc2 <- read.table("all_cor_soil_df_padj.cc2.txt", sep='\t', header = T)


nodeattrib_soil_combine.cc2 <- data.frame(node=union(all_cor_soil_df_padj.cc2$from,all_cor_soil_df_padj.cc2$to))
nodeattrib_soil_combine.cc2$kingdom <- 0

for (i in as.character(nodeattrib_soil_combine.cc2$node))
{
  if (i %in% rownames(df.otu.bac.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.cc2[nodeattrib_soil_combine.cc2$node==i,"kingdom"] <- "Bacteria"}
  
  else if (i %in% rownames(df.otu.arch.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.cc2[nodeattrib_soil_combine.cc2$node==i,"kingdom"] <- "Archaea"} 
  
  else
  {nodeattrib_soil_combine.cc2[nodeattrib_soil_combine.cc2$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_soil_combine.cc2) <- as.character(nodeattrib_soil_combine.cc2$node)
nodeattrib_soil_combine.cc2$kingdom



all_soil_net.cc2 <- graph_from_data_frame(all_cor_soil_df_padj.cc2,direct=F, vertices=nodeattrib_soil_combine.cc2)

## Number of nodes
length(V(all_soil_net.cc2)) #5712

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.cc2)))) #4244
length(grep("^F",names(V(all_soil_net.cc2)))) #1275
length(grep("^A",names(V(all_soil_net.cc2)))) #193

## Connections 
bb_occur_soil.cc2 <- droplevels(all_cor_soil_df_padj.cc2[with(all_cor_soil_df_padj.cc2, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.cc2) #366469

ff_occur_soil.cc2 <- droplevels(all_cor_soil_df_padj.cc2[with(all_cor_soil_df_padj.cc2, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.cc2) #32589

fb_occur_soil.cc2 <- droplevels(all_cor_soil_df_padj.cc2[with(all_cor_soil_df_padj.cc2, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.cc2) #168190

aa_occur_soil.cc2 <- droplevels(all_cor_soil_df_padj.cc2[with(all_cor_soil_df_padj.cc2, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.cc2) #729

ba_occur_soil.cc2 <- droplevels(all_cor_soil_df_padj.cc2[with(all_cor_soil_df_padj.cc2, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.cc2) #31907

fa_occur_soil.cc2 <- droplevels(all_cor_soil_df_padj.cc2[with(all_cor_soil_df_padj.cc2, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.cc2) #7099

## Network properties
meta_degree.cc2 <- sort(igraph::degree(all_soil_net.cc2,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.cc2))) #"max degree = " "399"     
print(c('mean degree = ', mean(meta_degree.cc2))) #"mean degree = "   "223.182421227197"
meta_degree.cc2 <- as.data.frame(meta_degree.cc2)
ggplot(meta_degree.cc2, aes(x=meta_degree.cc2)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.cc2, directed = FALSE))) #"average shortest path length =  3.59523728700828"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.cc2, "global"))) #"mean clustering coefficient =  0.949977913375762"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.cc2, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.cc2, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.cc2, V(all_soil_net.cc2)))))) #"mean number of neighbors =  223.182421227197"

##
net <- all_soil_net.cc2
CC2_all_deg <- igraph::degree(net,mode="all")
CC2_all_betweenness <- betweenness(net, normalized = TRUE)
CC2_all_closeness <- closeness(net, normalized = TRUE)
CC2_all_transitivity <- transitivity(net, "local", vids = V(net))
names(CC2_all_transitivity)<- V(net)$name
CC2_all_transitivity[is.na(CC2_all_transitivity)] <- 0


## Defining hub OTUs
n <-1
CC2_all_deg.1percent <- CC2_all_deg[CC2_all_deg >= quantile(CC2_all_deg,prob=1-n/100)]
length(CC2_all_deg.1percent) #305

CC2_all_betweenness.1percent <- CC2_all_betweenness[CC2_all_betweenness >= quantile(CC2_all_betweenness,prob=1-n/100)]
length(CC2_all_betweenness.1percent) #59

CC2_all_closeness.1percent <- CC2_all_closeness[CC2_all_closeness >= quantile(CC2_all_closeness,prob=1-n/100)]
length(CC2_all_closeness.1percent) #305

intersect(names(CC2_all_deg.1percent), names(CC2_all_betweenness.1percent))
intersect(names(CC2_all_deg.1percent), names(CC2_all_closeness.1percent))
intersect(names(CC2_all_betweenness.1percent), names(CC2_all_closeness.1percent))


## Set node shape soil_all_keystone.2017
V(all_soil_net.cc2)$shape <- V(all_soil_net.cc2)$name
V(all_soil_net.cc2)$shape[V(all_soil_net.cc2)$shape %in% names(CC2_all_betweenness.1percent)] <- "star"
#V(all_soil_net.cc2)$shape[V(all_soil_net.cc2)$shape %in% strictkeystone.cc2] <- "star"
V(all_soil_net.cc2)$shape[V(all_soil_net.cc2)$shape %in% rownames(df.otu.bac.clean.nolog.cc2)] <- "circle"
V(all_soil_net.cc2)$shape[V(all_soil_net.cc2)$shape %in% rownames(df.otu.fun.clean.nolog.cc2)] <- "triangle"
V(all_soil_net.cc2)$shape[V(all_soil_net.cc2)$shape %in% rownames(df.otu.arch.clean.nolog.cc2)] <- "square"

## Set node colors based kingdom
cs <- c("Bacteria", "Archaea", "Fungi")
unique(V(all_soil_net.cc2)$kingdom)
V(all_soil_net.cc2)$color <- V(all_soil_net.cc2)$kingdom
V(all_soil_net.cc2)$color[V(all_soil_net.cc2)$color == "Bacteria"] <- "#EE7600"
V(all_soil_net.cc2)$color[V(all_soil_net.cc2)$color == "Archaea"] <- "#458B74"
V(all_soil_net.cc2)$color[V(all_soil_net.cc2)$color == "Fungi"] <- "#9A32CD"
V(all_soil_net.cc2)$frame.color <- V(all_soil_net.cc2)$color

#Node size
V(all_soil_net.cc2)$size <- V(all_soil_net.cc2)$name
V(all_soil_net.cc2)$size[V(all_soil_net.cc2)$size %in% names(CC2_all_betweenness.1percent)] <- 3
#V(all_soil_net.cc2)$size[V(all_soil_net.cc2)$size %in% strictkeystone.2017] <- 5
V(all_soil_net.cc2)$size[V(all_soil_net.cc2)$size %in% rownames(df.otu.bac.clean.nolog.cc2)] <- 2
V(all_soil_net.cc2)$size[V(all_soil_net.cc2)$size %in% rownames(df.otu.fun.clean.nolog.cc2)] <- 3
V(all_soil_net.cc2)$size[V(all_soil_net.cc2)$size %in% rownames(df.otu.arch.clean.nolog.cc2)] <- 2
soilcombine_nodesizes.cc2 <- as.numeric(V(all_soil_net.cc2)$size)

## Set edge color
E(all_soil_net.cc2)$color <- ifelse(E(all_soil_net.cc2)$cor >0, "#99CCFF","#FF6666")


#Layout
source("star_shape.R")
source("triangle_shape.R")
library(qgraph)
all_soil_net.cc2.integer <- set.vertex.attribute(all_soil_net.cc2, "name", value=paste(1:5712))
edge.cc2 <- get.edgelist(all_soil_net.cc2.integer)
edge.cc2.num<- apply(edge.cc2, 2, as.numeric)

layout.cc2 <- qgraph.layout.fruchtermanreingold(edge.cc2.num,vcount=vcount(all_soil_net.cc2.integer))

#Default
plot(all_soil_net.cc2.integer,layout=layout.cc2,vertex.label=NA, vertex.size=soilcombine_nodesizes.cc2)

mtext("qgraph.layout.fruchtermanreingold default", side=1)

#Modified
layout.cc2.modified <- qgraph.layout.fruchtermanreingold(edge.cc2.num,vcount=vcount(all_soil_net.cc2.integer),
                                                         area=10*(vcount(all_soil_net.cc2.integer)^2),repulse.rad=(vcount(all_soil_net.cc2.integer)^3.1))

#without module
plot(all_soil_net.cc2.integer,layout=layout.cc2.modified,vertex.size=soilcombine_nodesizes.cc2,vertex.label=NA)

memory.limit(200000)

dev.off()



## Negative correlations
##Positive and negative connections
bb_occur_soil.cc1.neg<-subset(bb_occur_soil.cc1, bb_occur_soil.cc1$cor < 0)
nrow(bb_occur_soil.cc1.neg) #93

ff_occur_soil.cc1.neg<-subset(ff_occur_soil.cc1, ff_occur_soil.cc1$cor < 0)
nrow(ff_occur_soil.cc1.neg) #4059

aa_occur_soil.cc1.neg<-subset(aa_occur_soil.cc1, aa_occur_soil.cc1$cor < 0)
nrow(aa_occur_soil.cc1.neg) #3

ba_occur_soil.cc1.neg<-subset(ba_occur_soil.cc1, ba_occur_soil.cc1$cor < 0)
nrow(ba_occur_soil.cc1.neg) #37

fa_occur_soil.cc1.neg<-subset(fa_occur_soil.cc1, fa_occur_soil.cc1$cor < 0)
nrow(fa_occur_soil.cc1.neg) #166

fb_occur_soil.cc1.neg<-subset(fb_occur_soil.cc1, fb_occur_soil.cc1$cor < 0)
nrow(fb_occur_soil.cc1.neg) #812

neg.occur.cc1<-rbind(bb_occur_soil.cc1.neg, ff_occur_soil.cc1.neg, aa_occur_soil.cc1.neg,ba_occur_soil.cc1.neg,fa_occur_soil.cc1.neg,fb_occur_soil.cc1.neg)
write.xlsx(neg.occur.cc1, "Negative correlation in cc1 network_with archaea.xlsx")


bb_occur_soil.cc2.neg<-subset(bb_occur_soil.cc2, bb_occur_soil.cc2$cor < 0)
nrow(bb_occur_soil.cc2.neg) #107

ff_occur_soil.cc2.neg<-subset(ff_occur_soil.cc2, ff_occur_soil.cc2$cor < 0)
nrow(ff_occur_soil.cc2.neg) #3172

aa_occur_soil.cc2.neg<-subset(aa_occur_soil.cc2, aa_occur_soil.cc2$cor < 0)
nrow(aa_occur_soil.cc2.neg) #3

ba_occur_soil.cc2.neg<-subset(ba_occur_soil.cc2, ba_occur_soil.cc2$cor < 0)
nrow(ba_occur_soil.cc2.neg) #35

fa_occur_soil.cc2.neg<-subset(fa_occur_soil.cc2, fa_occur_soil.cc2$cor < 0)
nrow(fa_occur_soil.cc2.neg) #187

fb_occur_soil.cc2.neg<-subset(fb_occur_soil.cc2, fb_occur_soil.cc2$cor < 0)
nrow(fb_occur_soil.cc2.neg) #826

neg.occur.cc2<-rbind(bb_occur_soil.cc2.neg, ff_occur_soil.cc2.neg, aa_occur_soil.cc2.neg,ba_occur_soil.cc2.neg,fa_occur_soil.cc2.neg,fb_occur_soil.cc2.neg)
write.xlsx(neg.occur.cc2, "Negative correlation in cc2 network_with archaea.xlsx")


## CJ
#igraph
## CJ1
all_cor_soil_df_padj.cj1 <- read.table("all_cor_soil_df_padj.cj1.txt", sep='\t', header = T)


nodeattrib_soil_combine.cj1 <- data.frame(node=union(all_cor_soil_df_padj.cj1$from,all_cor_soil_df_padj.cj1$to))
nodeattrib_soil_combine.cj1$kingdom <- 0

for (i in as.character(nodeattrib_soil_combine.cj1$node))
{
  if (i %in% rownames(df.otu.bac.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.cj1[nodeattrib_soil_combine.cj1$node==i,"kingdom"] <- "Bacteria"}
  
  else if (i %in% rownames(df.otu.arch.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.cj1[nodeattrib_soil_combine.cj1$node==i,"kingdom"] <- "Archaea"} 
  
  else
  {nodeattrib_soil_combine.cj1[nodeattrib_soil_combine.cj1$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_soil_combine.cj1) <- as.character(nodeattrib_soil_combine.cj1$node)
nodeattrib_soil_combine.cj1$kingdom



all_soil_net.cj1 <- graph_from_data_frame(all_cor_soil_df_padj.cj1,direct=F, vertices=nodeattrib_soil_combine.cj1)

## Number of nodes
length(V(all_soil_net.cj1)) #6159

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.cj1)))) #4420
length(grep("^F",names(V(all_soil_net.cj1)))) #1555
length(grep("^A",names(V(all_soil_net.cj1)))) #184

## Connections 
bb_occur_soil.cj1 <- droplevels(all_cor_soil_df_padj.cj1[with(all_cor_soil_df_padj.cj1, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.cj1) #392591

ff_occur_soil.cj1 <- droplevels(all_cor_soil_df_padj.cj1[with(all_cor_soil_df_padj.cj1, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.cj1) #48816

fb_occur_soil.cj1 <- droplevels(all_cor_soil_df_padj.cj1[with(all_cor_soil_df_padj.cj1, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.cj1) #213090

aa_occur_soil.cj1 <- droplevels(all_cor_soil_df_padj.cj1[with(all_cor_soil_df_padj.cj1, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.cj1) #729

ba_occur_soil.cj1 <- droplevels(all_cor_soil_df_padj.cj1[with(all_cor_soil_df_padj.cj1, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.cj1) #34130

fa_occur_soil.cj1 <- droplevels(all_cor_soil_df_padj.cj1[with(all_cor_soil_df_padj.cj1, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.cj1) #9157

## Network properties
meta_degree.cj1 <- sort(igraph::degree(all_soil_net.cj1,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.cj1))) #"max degree = " "399"     
print(c('mean degree = ', mean(meta_degree.cj1))) #"mean degree = "   "223.182421227197"
meta_degree.cj1 <- as.data.frame(meta_degree.cj1)
ggplot(meta_degree.cj1, aes(x=meta_degree.cj1)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.cj1, directed = FALSE))) #"average shortest path length =  3.59523728700828"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.cj1, "global"))) #"mean clustering coefficient =  0.949977913375762"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.cj1, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.cj1, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.cj1, V(all_soil_net.cj1)))))) #"mean number of neighbors =  223.182421227197"

##
net <- all_soil_net.cj1
CJ1_all_deg <- igraph::degree(net,mode="all")
CJ1_all_betweenness <- betweenness(net, normalized = TRUE)
CJ1_all_closeness <- closeness(net, normalized = TRUE)
CJ1_all_transitivity <- transitivity(net, "local", vids = V(net))
names(CJ1_all_transitivity)<- V(net)$name
CJ1_all_transitivity[is.na(CJ1_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
CJ1_all_deg.1percent <- CJ1_all_deg[CJ1_all_deg >= quantile(CJ1_all_deg,prob=1-n/100)]
length(CJ1_all_deg.1percent) #317

CJ1_all_betweenness.1percent <- CJ1_all_betweenness[CJ1_all_betweenness >= quantile(CJ1_all_betweenness,prob=1-n/100)]
length(CJ1_all_betweenness.1percent) #62

CJ1_all_closeness.1percent <- CJ1_all_closeness[CJ1_all_closeness >= quantile(CJ1_all_closeness,prob=1-n/100)]
length(CJ1_all_closeness.1percent) #311

intersect(names(CJ1_all_deg.1percent), names(CJ1_all_betweenness.1percent))
intersect(names(CJ1_all_deg.1percent), names(CJ1_all_closeness.1percent))
intersect(names(CJ1_all_betweenness.1percent), names(CJ1_all_closeness.1percent))


## Set node shape soil_all_keystone.2017
V(all_soil_net.cj1)$shape <- V(all_soil_net.cj1)$name
V(all_soil_net.cj1)$shape[V(all_soil_net.cj1)$shape %in% names(CJ1_all_betweenness.1percent)] <- "star"
#V(all_soil_net.cj1)$shape[V(all_soil_net.cj1)$shape %in% strictkeystone.cj1] <- "star"
V(all_soil_net.cj1)$shape[V(all_soil_net.cj1)$shape %in% rownames(df.otu.bac.clean.nolog.cj1)] <- "circle"
V(all_soil_net.cj1)$shape[V(all_soil_net.cj1)$shape %in% rownames(df.otu.fun.clean.nolog.cj1)] <- "triangle"
V(all_soil_net.cj1)$shape[V(all_soil_net.cj1)$shape %in% rownames(df.otu.arch.clean.nolog.cj1)] <- "square"

## Set node colors based kingdom
cs <- c("Bacteria", "Archaea", "Fungi")
unique(V(all_soil_net.cj1)$kingdom)
V(all_soil_net.cj1)$color <- V(all_soil_net.cj1)$kingdom
V(all_soil_net.cj1)$color[V(all_soil_net.cj1)$color == "Bacteria"] <- "#EE7600"
V(all_soil_net.cj1)$color[V(all_soil_net.cj1)$color == "Archaea"] <- "#458B74"
V(all_soil_net.cj1)$color[V(all_soil_net.cj1)$color == "Fungi"] <- "#9A32CD"
V(all_soil_net.cj1)$frame.color <- V(all_soil_net.cj1)$color

#Node size
V(all_soil_net.cj1)$size <- V(all_soil_net.cj1)$name
V(all_soil_net.cj1)$size[V(all_soil_net.cj1)$size %in% names(CJ1_all_betweenness.1percent)] <- 3
#V(all_soil_net.cj1)$size[V(all_soil_net.cj1)$size %in% strictkeystone.2017] <- 5
V(all_soil_net.cj1)$size[V(all_soil_net.cj1)$size %in% rownames(df.otu.bac.clean.nolog.cj1)] <- 2
V(all_soil_net.cj1)$size[V(all_soil_net.cj1)$size %in% rownames(df.otu.fun.clean.nolog.cj1)] <- 3
V(all_soil_net.cj1)$size[V(all_soil_net.cj1)$size %in% rownames(df.otu.arch.clean.nolog.cj1)] <- 2
soilcombine_nodesizes.cj1 <- as.numeric(V(all_soil_net.cj1)$size)

## Set edge color
E(all_soil_net.cj1)$color <- ifelse(E(all_soil_net.cj1)$cor >0, "#99CCFF","#FF6666")


#Layout
source("star_shape.R")
source("triangle_shape.R")
library(qgraph)
all_soil_net.cj1.integer <- set.vertex.attribute(all_soil_net.cj1, "name", value=paste(1:6159))
edge.cj1 <- get.edgelist(all_soil_net.cj1.integer)
edge.cj1.num<- apply(edge.cj1, 2, as.numeric)

layout.cj1 <- qgraph.layout.fruchtermanreingold(edge.cj1.num,vcount=vcount(all_soil_net.cj1.integer))

#Default
plot(all_soil_net.cj1.integer,layout=layout.cj1,vertex.label=NA, vertex.size=soilcombine_nodesizes.cj1)

mtext("qgraph.layout.fruchtermanreingold default", side=1)

#Modified
layout.cj1.modified <- qgraph.layout.fruchtermanreingold(edge.cj1.num,vcount=vcount(all_soil_net.cj1.integer),
                                                         area=10*(vcount(all_soil_net.cj1.integer)^2),repulse.rad=(vcount(all_soil_net.cj1.integer)^3.1))

#without module
plot(all_soil_net.cj1.integer,layout=layout.cj1.modified,vertex.size=soilcombine_nodesizes.cj1,vertex.label=NA)

memory.limit(200000)

dev.off()

## CJ2
all_cor_soil_df_padj.cj2 <- read.table("all_cor_soil_df_padj.cj2.txt", sep='\t', header = T)


nodeattrib_soil_combine.cj2 <- data.frame(node=union(all_cor_soil_df_padj.cj2$from,all_cor_soil_df_padj.cj2$to))
nodeattrib_soil_combine.cj2$kingdom <- 0

for (i in as.character(nodeattrib_soil_combine.cj2$node))
{
  if (i %in% rownames(df.otu.bac.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.cj2[nodeattrib_soil_combine.cj2$node==i,"kingdom"] <- "Bacteria"}
  
  else if (i %in% rownames(df.otu.arch.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.cj2[nodeattrib_soil_combine.cj2$node==i,"kingdom"] <- "Archaea"} 
  
  else
  {nodeattrib_soil_combine.cj2[nodeattrib_soil_combine.cj2$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_soil_combine.cj2) <- as.character(nodeattrib_soil_combine.cj2$node)
nodeattrib_soil_combine.cj2$kingdom



all_soil_net.cj2 <- graph_from_data_frame(all_cor_soil_df_padj.cj2,direct=F, vertices=nodeattrib_soil_combine.cj2)

## Number of nodes
length(V(all_soil_net.cj2)) #4830

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.cj2)))) #3405
length(grep("^F",names(V(all_soil_net.cj2)))) #1298
length(grep("^A",names(V(all_soil_net.cj2)))) #127

## Connections 
bb_occur_soil.cj2 <- droplevels(all_cor_soil_df_padj.cj2[with(all_cor_soil_df_padj.cj2, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.cj2) #232581

ff_occur_soil.cj2 <- droplevels(all_cor_soil_df_padj.cj2[with(all_cor_soil_df_padj.cj2, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.cj2) #33007

fb_occur_soil.cj2 <- droplevels(all_cor_soil_df_padj.cj2[with(all_cor_soil_df_padj.cj2, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.cj2) #148087

aa_occur_soil.cj2 <- droplevels(all_cor_soil_df_padj.cj2[with(all_cor_soil_df_padj.cj2, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.cj2) #388

ba_occur_soil.cj2 <- droplevels(all_cor_soil_df_padj.cj2[with(all_cor_soil_df_padj.cj2, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.cj2) #17940

fa_occur_soil.cj2 <- droplevels(all_cor_soil_df_padj.cj2[with(all_cor_soil_df_padj.cj2, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.cj2) #5512

## Network properties
meta_degree.cj2 <- sort(igraph::degree(all_soil_net.cj2,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.cj2))) #"max degree = " "399"     
print(c('mean degree = ', mean(meta_degree.cj2))) #"mean degree = "   "223.182421227197"
meta_degree.cj2 <- as.data.frame(meta_degree.cj2)
ggplot(meta_degree.cj2, aes(x=meta_degree.cj2)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.cj2, directed = FALSE))) #"average shortest path length =  3.59523728700828"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.cj2, "global"))) #"mean clustering coefficient =  0.949977913375762"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.cj2, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.cj2, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.cj2, V(all_soil_net.cj2)))))) #"mean number of neighbors =  223.182421227197"

##
net <- all_soil_net.cj2
CJ2_all_deg <- igraph::degree(net,mode="all")
CJ2_all_betweenness <- betweenness(net, normalized = TRUE)
CJ2_all_closeness <- closeness(net, normalized = TRUE)
CJ2_all_transitivity <- transitivity(net, "local", vids = V(net))
names(CJ2_all_transitivity)<- V(net)$name
CJ2_all_transitivity[is.na(CJ2_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
CJ2_all_deg.1percent <- CJ2_all_deg[CJ2_all_deg >= quantile(CJ2_all_deg,prob=1-n/100)]
length(CJ2_all_deg.1percent) #305

CJ2_all_betweenness.1percent <- CJ2_all_betweenness[CJ2_all_betweenness >= quantile(CJ2_all_betweenness,prob=1-n/100)]
length(CJ2_all_betweenness.1percent) #59

CJ2_all_closeness.1percent <- CJ2_all_closeness[CJ2_all_closeness >= quantile(CJ2_all_closeness,prob=1-n/100)]
length(CJ2_all_closeness.1percent) #305

intersect(names(CJ2_all_deg.1percent), names(CJ2_all_betweenness.1percent))
intersect(names(CJ2_all_deg.1percent), names(CJ2_all_closeness.1percent))
intersect(names(CJ2_all_betweenness.1percent), names(CJ2_all_closeness.1percent))


## Set node shape soil_all_keystone.2017
V(all_soil_net.cj2)$shape <- V(all_soil_net.cj2)$name
V(all_soil_net.cj2)$shape[V(all_soil_net.cj2)$shape %in% names(CJ2_all_betweenness.1percent)] <- "star"
#V(all_soil_net.cj2)$shape[V(all_soil_net.cj2)$shape %in% strictkeystone.cj2] <- "star"
V(all_soil_net.cj2)$shape[V(all_soil_net.cj2)$shape %in% rownames(df.otu.bac.clean.nolog.cj2)] <- "circle"
V(all_soil_net.cj2)$shape[V(all_soil_net.cj2)$shape %in% rownames(df.otu.fun.clean.nolog.cj2)] <- "triangle"
V(all_soil_net.cj2)$shape[V(all_soil_net.cj2)$shape %in% rownames(df.otu.arch.clean.nolog.cj2)] <- "square"

## Set node colors based kingdom
cs <- c("Bacteria", "Archaea", "Fungi")
unique(V(all_soil_net.cj2)$kingdom)
V(all_soil_net.cj2)$color <- V(all_soil_net.cj2)$kingdom
V(all_soil_net.cj2)$color[V(all_soil_net.cj2)$color == "Bacteria"] <- "#EE7600"
V(all_soil_net.cj2)$color[V(all_soil_net.cj2)$color == "Archaea"] <- "#458B74"
V(all_soil_net.cj2)$color[V(all_soil_net.cj2)$color == "Fungi"] <- "#9A32CD"
V(all_soil_net.cj2)$frame.color <- V(all_soil_net.cj2)$color

#Node size
V(all_soil_net.cj2)$size <- V(all_soil_net.cj2)$name
V(all_soil_net.cj2)$size[V(all_soil_net.cj2)$size %in% names(CJ2_all_betweenness.1percent)] <- 3
#V(all_soil_net.cj2)$size[V(all_soil_net.cj2)$size %in% strictkeystone.2017] <- 5
V(all_soil_net.cj2)$size[V(all_soil_net.cj2)$size %in% rownames(df.otu.bac.clean.nolog.cj2)] <- 2
V(all_soil_net.cj2)$size[V(all_soil_net.cj2)$size %in% rownames(df.otu.fun.clean.nolog.cj2)] <- 3
V(all_soil_net.cj2)$size[V(all_soil_net.cj2)$size %in% rownames(df.otu.arch.clean.nolog.cj2)] <- 2
soilcombine_nodesizes.cj2 <- as.numeric(V(all_soil_net.cj2)$size)

## Set edge color
E(all_soil_net.cj2)$color <- ifelse(E(all_soil_net.cj2)$cor >0, "#99CCFF","#FF6666")


#Layout
source("star_shape.R")
source("triangle_shape.R")
library(qgraph)
all_soil_net.cj2.integer <- set.vertex.attribute(all_soil_net.cj2, "name", value=paste(1:4830))
edge.cj2 <- get.edgelist(all_soil_net.cj2.integer)
edge.cj2.num<- apply(edge.cj2, 2, as.numeric)

layout.cj2 <- qgraph.layout.fruchtermanreingold(edge.cj2.num,vcount=vcount(all_soil_net.cj2.integer))

#Default
plot(all_soil_net.cj2.integer,layout=layout.cj2,vertex.label=NA, vertex.size=soilcombine_nodesizes.cj2)

mtext("qgraph.layout.fruchtermanreingold default", side=1)

#Modified
layout.cj2.modified <- qgraph.layout.fruchtermanreingold(edge.cj2.num,vcount=vcount(all_soil_net.cj2.integer),
                                                         area=10*(vcount(all_soil_net.cj2.integer)^2),repulse.rad=(vcount(all_soil_net.cj2.integer)^3.1))

#without module
plot(all_soil_net.cj2.integer,layout=layout.cj2.modified,vertex.size=soilcombine_nodesizes.cj2,vertex.label=NA)

memory.limit(200000)

dev.off()




## Negative correlations
##Positive and negative connections
bb_occur_soil.cj1.neg<-subset(bb_occur_soil.cj1, bb_occur_soil.cj1$cor < 0)
nrow(bb_occur_soil.cj1.neg) #116

ff_occur_soil.cj1.neg<-subset(ff_occur_soil.cj1, ff_occur_soil.cj1$cor < 0)
nrow(ff_occur_soil.cj1.neg) #4918

aa_occur_soil.cj1.neg<-subset(aa_occur_soil.cj1, aa_occur_soil.cj1$cor < 0)
nrow(aa_occur_soil.cj1.neg) #1

ba_occur_soil.cj1.neg<-subset(ba_occur_soil.cj1, ba_occur_soil.cj1$cor < 0)
nrow(ba_occur_soil.cj1.neg) #20

fa_occur_soil.cj1.neg<-subset(fa_occur_soil.cj1, fa_occur_soil.cj1$cor < 0)
nrow(fa_occur_soil.cj1.neg) #94

fb_occur_soil.cj1.neg<-subset(fb_occur_soil.cj1, fb_occur_soil.cj1$cor < 0)
nrow(fb_occur_soil.cj1.neg) #1335

neg.occur.cj1<-rbind(bb_occur_soil.cj1.neg, ff_occur_soil.cj1.neg, aa_occur_soil.cj1.neg,ba_occur_soil.cj1.neg,fa_occur_soil.cj1.neg,fb_occur_soil.cj1.neg)
write.xlsx(neg.occur.cj1, "Negative correlation in cj1 network_with archaea.xlsx")


bb_occur_soil.cj2.neg<-subset(bb_occur_soil.cj2, bb_occur_soil.cj2$cor < 0)
nrow(bb_occur_soil.cj2.neg) #125

ff_occur_soil.cj2.neg<-subset(ff_occur_soil.cj2, ff_occur_soil.cj2$cor < 0)
nrow(ff_occur_soil.cj2.neg) #2454

aa_occur_soil.cj2.neg<-subset(aa_occur_soil.cj2, aa_occur_soil.cj2$cor < 0)
nrow(aa_occur_soil.cj2.neg) #6

ba_occur_soil.cj2.neg<-subset(ba_occur_soil.cj2, ba_occur_soil.cj2$cor < 0)
nrow(ba_occur_soil.cj2.neg) #34

fa_occur_soil.cj2.neg<-subset(fa_occur_soil.cj2, fa_occur_soil.cj2$cor < 0)
nrow(fa_occur_soil.cj2.neg) #157

fb_occur_soil.cj2.neg<-subset(fb_occur_soil.cj2, fb_occur_soil.cj2$cor < 0)
nrow(fb_occur_soil.cj2.neg) #852

neg.occur.cj2<-rbind(bb_occur_soil.cj2.neg, ff_occur_soil.cj2.neg, aa_occur_soil.cj2.neg,ba_occur_soil.cj2.neg,fa_occur_soil.cj2.neg,fb_occur_soil.cj2.neg)
write.xlsx(neg.occur.cj2, "Negative correlation in cj2 network_with archaea.xlsx")


## DG
#igraph
## DG1
all_cor_soil_df_padj.dg1 <- read.table("all_cor_soil_df_padj.dg1.txt", sep='\t', header = T)


nodeattrib_soil_combine.dg1 <- data.frame(node=union(all_cor_soil_df_padj.dg1$from,all_cor_soil_df_padj.dg1$to))
nodeattrib_soil_combine.dg1$kingdom <- 0

for (i in as.character(nodeattrib_soil_combine.dg1$node))
{
  if (i %in% rownames(df.otu.bac.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.dg1[nodeattrib_soil_combine.dg1$node==i,"kingdom"] <- "Bacteria"}
  
  else if (i %in% rownames(df.otu.arch.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.dg1[nodeattrib_soil_combine.dg1$node==i,"kingdom"] <- "Archaea"} 
  
  else
  {nodeattrib_soil_combine.dg1[nodeattrib_soil_combine.dg1$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_soil_combine.dg1) <- as.character(nodeattrib_soil_combine.dg1$node)
nodeattrib_soil_combine.dg1$kingdom



all_soil_net.dg1 <- graph_from_data_frame(all_cor_soil_df_padj.dg1,direct=F, vertices=nodeattrib_soil_combine.dg1)

## Number of nodes
length(V(all_soil_net.dg1)) #5520

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.dg1)))) #3828
length(grep("^F",names(V(all_soil_net.dg1)))) #1492
length(grep("^A",names(V(all_soil_net.dg1)))) #200

## Connections 
bb_occur_soil.dg1 <- droplevels(all_cor_soil_df_padj.dg1[with(all_cor_soil_df_padj.dg1, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.dg1) #280073

ff_occur_soil.dg1 <- droplevels(all_cor_soil_df_padj.dg1[with(all_cor_soil_df_padj.dg1, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.dg1) # 45102

fb_occur_soil.dg1 <- droplevels(all_cor_soil_df_padj.dg1[with(all_cor_soil_df_padj.dg1, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.dg1) #177632

aa_occur_soil.dg1 <- droplevels(all_cor_soil_df_padj.dg1[with(all_cor_soil_df_padj.dg1, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.dg1) #903

ba_occur_soil.dg1 <- droplevels(all_cor_soil_df_padj.dg1[with(all_cor_soil_df_padj.dg1, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.dg1) #30463

fa_occur_soil.dg1 <- droplevels(all_cor_soil_df_padj.dg1[with(all_cor_soil_df_padj.dg1, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.dg1) #9550

## Network properties
meta_degree.dg1 <- sort(igraph::degree(all_soil_net.dg1,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.dg1))) #"max degree = " "399"     
print(c('mean degree = ', mean(meta_degree.dg1))) #"mean degree = "   "223.182421227197"
meta_degree.dg1 <- as.data.frame(meta_degree.dg1)
ggplot(meta_degree.dg1, aes(x=meta_degree.dg1)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.dg1, directed = FALSE))) #"average shortest path length =  3.59523728700828"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.dg1, "global"))) #"mean clustering coefficient =  0.949977913375762"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.dg1, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.dg1, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.dg1, V(all_soil_net.dg1)))))) #"mean number of neighbors =  223.182421227197"

##
net <- all_soil_net.dg1
DG1_all_deg <- igraph::degree(net,mode="all")
DG1_all_betweenness <- betweenness(net, normalized = TRUE)
DG1_all_closeness <- closeness(net, normalized = TRUE)
DG1_all_transitivity <- transitivity(net, "local", vids = V(net))
names(DG1_all_transitivity)<- V(net)$name
DG1_all_transitivity[is.na(DG1_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
DG1_all_deg.1percent <- DG1_all_deg[DG1_all_deg >= quantile(DG1_all_deg,prob=1-n/100)]
length(DG1_all_deg.1percent) #317

DG1_all_betweenness.1percent <- DG1_all_betweenness[DG1_all_betweenness >= quantile(DG1_all_betweenness,prob=1-n/100)]
length(DG1_all_betweenness.1percent) #62

DG1_all_closeness.1percent <- DG1_all_closeness[DG1_all_closeness >= quantile(DG1_all_closeness,prob=1-n/100)]
length(DG1_all_closeness.1percent) #311

intersect(names(DG1_all_deg.1percent), names(DG1_all_betweenness.1percent))
intersect(names(DG1_all_deg.1percent), names(DG1_all_closeness.1percent))
intersect(names(DG1_all_betweenness.1percent), names(DG1_all_closeness.1percent))


## Set node shape soil_all_keystone.2017
V(all_soil_net.dg1)$shape <- V(all_soil_net.dg1)$name
V(all_soil_net.dg1)$shape[V(all_soil_net.dg1)$shape %in% names(DG1_all_betweenness.1percent)] <- "star"
#V(all_soil_net.dg1)$shape[V(all_soil_net.dg1)$shape %in% strictkeystone.dg1] <- "star"
V(all_soil_net.dg1)$shape[V(all_soil_net.dg1)$shape %in% rownames(df.otu.bac.clean.nolog.dg1)] <- "circle"
V(all_soil_net.dg1)$shape[V(all_soil_net.dg1)$shape %in% rownames(df.otu.fun.clean.nolog.dg1)] <- "triangle"
V(all_soil_net.dg1)$shape[V(all_soil_net.dg1)$shape %in% rownames(df.otu.arch.clean.nolog.dg1)] <- "square"

## Set node colors based kingdom
cs <- c("Bacteria", "Archaea", "Fungi")
unique(V(all_soil_net.dg1)$kingdom)
V(all_soil_net.dg1)$color <- V(all_soil_net.dg1)$kingdom
V(all_soil_net.dg1)$color[V(all_soil_net.dg1)$color == "Bacteria"] <- "#EE7600"
V(all_soil_net.dg1)$color[V(all_soil_net.dg1)$color == "Archaea"] <- "#458B74"
V(all_soil_net.dg1)$color[V(all_soil_net.dg1)$color == "Fungi"] <- "#9A32CD"
V(all_soil_net.dg1)$frame.color <- V(all_soil_net.dg1)$color

#Node size
V(all_soil_net.dg1)$size <- V(all_soil_net.dg1)$name
V(all_soil_net.dg1)$size[V(all_soil_net.dg1)$size %in% names(DG1_all_betweenness.1percent)] <- 3
#V(all_soil_net.dg1)$size[V(all_soil_net.dg1)$size %in% strictkeystone.2017] <- 5
V(all_soil_net.dg1)$size[V(all_soil_net.dg1)$size %in% rownames(df.otu.bac.clean.nolog.dg1)] <- 2
V(all_soil_net.dg1)$size[V(all_soil_net.dg1)$size %in% rownames(df.otu.fun.clean.nolog.dg1)] <- 3
V(all_soil_net.dg1)$size[V(all_soil_net.dg1)$size %in% rownames(df.otu.arch.clean.nolog.dg1)] <- 2
soilcombine_nodesizes.dg1 <- as.numeric(V(all_soil_net.dg1)$size)

## Set edge color
E(all_soil_net.dg1)$color <- ifelse(E(all_soil_net.dg1)$cor >0, "#99CCFF","#FF6666")


#Layout
source("star_shape.R")
source("triangle_shape.R")
library(qgraph)
all_soil_net.dg1.integer <- set.vertex.attribute(all_soil_net.dg1, "name", value=paste(1:5520))
edge.dg1 <- get.edgelist(all_soil_net.dg1.integer)
edge.dg1.num<- apply(edge.dg1, 2, as.numeric)

layout.dg1 <- qgraph.layout.fruchtermanreingold(edge.dg1.num,vcount=vcount(all_soil_net.dg1.integer))

#Default
plot(all_soil_net.dg1.integer,layout=layout.dg1,vertex.label=NA, vertex.size=soilcombine_nodesizes.dg1)

mtext("qgraph.layout.fruchtermanreingold default", side=1)

#Modified
layout.dg1.modified <- qgraph.layout.fruchtermanreingold(edge.dg1.num,vcount=vcount(all_soil_net.dg1.integer),
                                                         area=10*(vcount(all_soil_net.dg1.integer)^2),repulse.rad=(vcount(all_soil_net.dg1.integer)^3.1))

#without module
plot(all_soil_net.dg1.integer,layout=layout.dg1.modified,vertex.size=soilcombine_nodesizes.dg1,vertex.label=NA)

memory.limit(200000)

dev.off()

## DG2
all_cor_soil_df_padj.dg2 <- read.table("all_cor_soil_df_padj.dg2.txt", sep='\t', header = T)


nodeattrib_soil_combine.dg2 <- data.frame(node=union(all_cor_soil_df_padj.dg2$from,all_cor_soil_df_padj.dg2$to))
nodeattrib_soil_combine.dg2$kingdom <- 0

for (i in as.character(nodeattrib_soil_combine.dg2$node))
{
  if (i %in% rownames(df.otu.bac.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.dg2[nodeattrib_soil_combine.dg2$node==i,"kingdom"] <- "Bacteria"}
  
  else if (i %in% rownames(df.otu.arch.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.dg2[nodeattrib_soil_combine.dg2$node==i,"kingdom"] <- "Archaea"} 
  
  else
  {nodeattrib_soil_combine.dg2[nodeattrib_soil_combine.dg2$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_soil_combine.dg2) <- as.character(nodeattrib_soil_combine.dg2$node)
nodeattrib_soil_combine.dg2$kingdom



all_soil_net.dg2 <- graph_from_data_frame(all_cor_soil_df_padj.dg2,direct=F, vertices=nodeattrib_soil_combine.dg2)

## Number of nodes
length(V(all_soil_net.dg2)) #5733

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.dg2)))) #3700
length(grep("^F",names(V(all_soil_net.dg2)))) #1838
length(grep("^A",names(V(all_soil_net.dg2)))) #195

## Connections 
bb_occur_soil.dg2 <- droplevels(all_cor_soil_df_padj.dg2[with(all_cor_soil_df_padj.dg2, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.dg2) #247525

ff_occur_soil.dg2 <- droplevels(all_cor_soil_df_padj.dg2[with(all_cor_soil_df_padj.dg2, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.dg2) #62027

fb_occur_soil.dg2 <- droplevels(all_cor_soil_df_padj.dg2[with(all_cor_soil_df_padj.dg2, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.dg2) #206490

aa_occur_soil.dg2 <- droplevels(all_cor_soil_df_padj.dg2[with(all_cor_soil_df_padj.dg2, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.dg2) #903

ba_occur_soil.dg2 <- droplevels(all_cor_soil_df_padj.dg2[with(all_cor_soil_df_padj.dg2, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.dg2) #28312

fa_occur_soil.dg2 <- droplevels(all_cor_soil_df_padj.dg2[with(all_cor_soil_df_padj.dg2, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.dg2) #10716

## Network properties
meta_degree.dg2 <- sort(igraph::degree(all_soil_net.dg2,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.dg2))) #"max degree = " "399"     
print(c('mean degree = ', mean(meta_degree.dg2))) #"mean degree = "   "223.182421227197"
meta_degree.dg2 <- as.data.frame(meta_degree.dg2)
ggplot(meta_degree.dg2, aes(x=meta_degree.dg2)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.dg2, directed = FALSE))) #"average shortest path length =  3.59523728700828"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.dg2, "global"))) #"mean clustering coefficient =  0.949977913375762"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.dg2, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.dg2, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.dg2, V(all_soil_net.dg2)))))) #"mean number of neighbors =  223.182421227197"

##
net <- all_soil_net.dg2
DG2_all_deg <- igraph::degree(net,mode="all")
DG2_all_betweenness <- betweenness(net, normalized = TRUE)
DG2_all_closeness <- closeness(net, normalized = TRUE)
DG2_all_transitivity <- transitivity(net, "local", vids = V(net))
names(DG2_all_transitivity)<- V(net)$name
DG2_all_transitivity[is.na(DG2_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
DG2_all_deg.1percent <- DG2_all_deg[DG2_all_deg >= quantile(DG2_all_deg,prob=1-n/100)]
length(DG2_all_deg.1percent) #305

DG2_all_betweenness.1percent <- DG2_all_betweenness[DG2_all_betweenness >= quantile(DG2_all_betweenness,prob=1-n/100)]
length(DG2_all_betweenness.1percent) #58

DG2_all_closeness.1percent <- DG2_all_closeness[DG2_all_closeness >= quantile(DG2_all_closeness,prob=1-n/100)]
length(DG2_all_closeness.1percent) #305

intersect(names(DG2_all_deg.1percent), names(DG2_all_betweenness.1percent))
intersect(names(DG2_all_deg.1percent), names(DG2_all_closeness.1percent))
intersect(names(DG2_all_betweenness.1percent), names(DG2_all_closeness.1percent))


## Set node shape soil_all_keystone.2017
V(all_soil_net.dg2)$shape <- V(all_soil_net.dg2)$name
V(all_soil_net.dg2)$shape[V(all_soil_net.dg2)$shape %in% names(DG2_all_betweenness.1percent)] <- "star"
#V(all_soil_net.dg2)$shape[V(all_soil_net.dg2)$shape %in% strictkeystone.dg2] <- "star"
V(all_soil_net.dg2)$shape[V(all_soil_net.dg2)$shape %in% rownames(df.otu.bac.clean.nolog.dg2)] <- "circle"
V(all_soil_net.dg2)$shape[V(all_soil_net.dg2)$shape %in% rownames(df.otu.fun.clean.nolog.dg2)] <- "triangle"
V(all_soil_net.dg2)$shape[V(all_soil_net.dg2)$shape %in% rownames(df.otu.arch.clean.nolog.dg2)] <- "square"

## Set node colors based kingdom
cs <- c("Bacteria", "Archaea", "Fungi")
unique(V(all_soil_net.dg2)$kingdom)
V(all_soil_net.dg2)$color <- V(all_soil_net.dg2)$kingdom
V(all_soil_net.dg2)$color[V(all_soil_net.dg2)$color == "Bacteria"] <- "#EE7600"
V(all_soil_net.dg2)$color[V(all_soil_net.dg2)$color == "Archaea"] <- "#458B74"
V(all_soil_net.dg2)$color[V(all_soil_net.dg2)$color == "Fungi"] <- "#9A32CD"
V(all_soil_net.dg2)$frame.color <- V(all_soil_net.dg2)$color

#Node size
V(all_soil_net.dg2)$size <- V(all_soil_net.dg2)$name
V(all_soil_net.dg2)$size[V(all_soil_net.dg2)$size %in% names(DG2_all_betweenness.1percent)] <- 3
#V(all_soil_net.dg2)$size[V(all_soil_net.dg2)$size %in% strictkeystone.2017] <- 5
V(all_soil_net.dg2)$size[V(all_soil_net.dg2)$size %in% rownames(df.otu.bac.clean.nolog.dg2)] <- 2
V(all_soil_net.dg2)$size[V(all_soil_net.dg2)$size %in% rownames(df.otu.fun.clean.nolog.dg2)] <- 3
V(all_soil_net.dg2)$size[V(all_soil_net.dg2)$size %in% rownames(df.otu.arch.clean.nolog.dg2)] <- 2
soilcombine_nodesizes.dg2 <- as.numeric(V(all_soil_net.dg2)$size)

## Set edge color
E(all_soil_net.dg2)$color <- ifelse(E(all_soil_net.dg2)$cor >0, "#99CCFF","#FF6666")


#Layout
source("star_shape.R")
source("triangle_shape.R")
library(qgraph)
all_soil_net.dg2.integer <- set.vertex.attribute(all_soil_net.dg2, "name", value=paste(1:5733))
edge.dg2 <- get.edgelist(all_soil_net.dg2.integer)
edge.dg2.num<- apply(edge.dg2, 2, as.numeric)

layout.dg2 <- qgraph.layout.fruchtermanreingold(edge.dg2.num,vcount=vcount(all_soil_net.dg2.integer))

#Default
plot(all_soil_net.dg2.integer,layout=layout.dg2,vertex.label=NA, vertex.size=soilcombine_nodesizes.dg2)

mtext("qgraph.layout.fruchtermanreingold default", side=1)

#Modified
layout.dg2.modified <- qgraph.layout.fruchtermanreingold(edge.dg2.num,vcount=vcount(all_soil_net.dg2.integer),
                                                         area=10*(vcount(all_soil_net.dg2.integer)^2),repulse.rad=(vcount(all_soil_net.dg2.integer)^3.1))

#without module
plot(all_soil_net.dg2.integer,layout=layout.dg2.modified,vertex.size=soilcombine_nodesizes.dg2,vertex.label=NA)

memory.limit(200000)

dev.off()




## Negative correlations
##Positive and negative connections
bb_occur_soil.dg1.neg<-subset(bb_occur_soil.dg1, bb_occur_soil.dg1$cor < 0)
nrow(bb_occur_soil.dg1.neg) #193

ff_occur_soil.dg1.neg<-subset(ff_occur_soil.dg1, ff_occur_soil.dg1$cor < 0)
nrow(ff_occur_soil.dg1.neg) #5327

aa_occur_soil.dg1.neg<-subset(aa_occur_soil.dg1, aa_occur_soil.dg1$cor < 0)
nrow(aa_occur_soil.dg1.neg) #4

ba_occur_soil.dg1.neg<-subset(ba_occur_soil.dg1, ba_occur_soil.dg1$cor < 0)
nrow(ba_occur_soil.dg1.neg) #49

fa_occur_soil.dg1.neg<-subset(fa_occur_soil.dg1, fa_occur_soil.dg1$cor < 0)
nrow(fa_occur_soil.dg1.neg) #255

fb_occur_soil.dg1.neg<-subset(fb_occur_soil.dg1, fb_occur_soil.dg1$cor < 0)
nrow(fb_occur_soil.dg1.neg) #1714

neg.occur.dg1<-rbind(bb_occur_soil.dg1.neg, ff_occur_soil.dg1.neg, aa_occur_soil.dg1.neg,ba_occur_soil.dg1.neg,fa_occur_soil.dg1.neg,fb_occur_soil.dg1.neg)
write.xlsx(neg.occur.dg1, "Negative correlation in dg1 network_with archaea.xlsx")


bb_occur_soil.dg2.neg<-subset(bb_occur_soil.dg2, bb_occur_soil.dg2$cor < 0)
nrow(bb_occur_soil.dg2.neg) #115

ff_occur_soil.dg2.neg<-subset(ff_occur_soil.dg2, ff_occur_soil.dg2$cor < 0)
nrow(ff_occur_soil.dg2.neg) #4889

aa_occur_soil.dg2.neg<-subset(aa_occur_soil.dg2, aa_occur_soil.dg2$cor < 0)
nrow(aa_occur_soil.dg2.neg) #1

ba_occur_soil.dg2.neg<-subset(ba_occur_soil.dg2, ba_occur_soil.dg2$cor < 0)
nrow(ba_occur_soil.dg2.neg) #22

fa_occur_soil.dg2.neg<-subset(fa_occur_soil.dg2, fa_occur_soil.dg2$cor < 0)
nrow(fa_occur_soil.dg2.neg) #116

fb_occur_soil.dg2.neg<-subset(fb_occur_soil.dg2, fb_occur_soil.dg2$cor < 0)
nrow(fb_occur_soil.dg2.neg) #1230

neg.occur.dg2<-rbind(bb_occur_soil.dg2.neg, ff_occur_soil.dg2.neg, aa_occur_soil.dg2.neg,ba_occur_soil.dg2.neg,fa_occur_soil.dg2.neg,fb_occur_soil.dg2.neg)
write.xlsx(neg.occur.dg2, "Negative correlation in dg2 network_with archaea.xlsx")



##IS
#igraph
## IS1
all_cor_soil_df_padj.is1 <- read.table("all_cor_soil_df_padj.is1.txt", sep='\t', header = T)


nodeattrib_soil_combine.is1 <- data.frame(node=union(all_cor_soil_df_padj.is1$from,all_cor_soil_df_padj.is1$to))
nodeattrib_soil_combine.is1$kingdom <- 0

for (i in as.character(nodeattrib_soil_combine.is1$node))
{
  if (i %in% rownames(df.otu.bac.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.is1[nodeattrib_soil_combine.is1$node==i,"kingdom"] <- "Bacteria"}
  
  else if (i %in% rownames(df.otu.arch.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.is1[nodeattrib_soil_combine.is1$node==i,"kingdom"] <- "Archaea"} 
  
  else
  {nodeattrib_soil_combine.is1[nodeattrib_soil_combine.is1$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_soil_combine.is1) <- as.character(nodeattrib_soil_combine.is1$node)
nodeattrib_soil_combine.is1$kingdom



all_soil_net.is1 <- graph_from_data_frame(all_cor_soil_df_padj.is1,direct=F, vertices=nodeattrib_soil_combine.is1)

## Number of nodes
length(V(all_soil_net.is1)) #5322

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.is1)))) #3488
length(grep("^F",names(V(all_soil_net.is1)))) #1639
length(grep("^A",names(V(all_soil_net.is1)))) #195

## Connections 
bb_occur_soil.is1 <- droplevels(all_cor_soil_df_padj.is1[with(all_cor_soil_df_padj.is1, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.is1) #239534

ff_occur_soil.is1 <- droplevels(all_cor_soil_df_padj.is1[with(all_cor_soil_df_padj.is1, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.is1) #49105

fb_occur_soil.is1 <- droplevels(all_cor_soil_df_padj.is1[with(all_cor_soil_df_padj.is1, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.is1) #162049

aa_occur_soil.is1 <- droplevels(all_cor_soil_df_padj.is1[with(all_cor_soil_df_padj.is1, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.is1) #992

ba_occur_soil.is1 <- droplevels(all_cor_soil_df_padj.is1[with(all_cor_soil_df_padj.is1, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.is1) #29433

fa_occur_soil.is1 <- droplevels(all_cor_soil_df_padj.is1[with(all_cor_soil_df_padj.is1, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.is1) #8524

## Network properties
meta_degree.is1 <- sort(igraph::degree(all_soil_net.is1,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.is1))) #"max degree = " "399"     
print(c('mean degree = ', mean(meta_degree.is1))) #"mean degree = "   "223.182421227197"
meta_degree.is1 <- as.data.frame(meta_degree.is1)
ggplot(meta_degree.is1, aes(x=meta_degree.is1)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.is1, directed = FALSE))) #"average shortest path length =  3.59523728700828"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.is1, "global"))) #"mean clustering coefficient =  0.949977913375762"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.is1, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.is1, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.is1, V(all_soil_net.is1)))))) #"mean number of neighbors =  223.182421227197"

##
net <- all_soil_net.is1
IS1_all_deg <- igraph::degree(net,mode="all")
IS1_all_betweenness <- betweenness(net, normalized = TRUE)
IS1_all_closeness <- closeness(net, normalized = TRUE)
IS1_all_transitivity <- transitivity(net, "local", vids = V(net))
names(IS1_all_transitivity)<- V(net)$name
IS1_all_transitivity[is.na(IS1_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
IS1_all_deg.1percent <- IS1_all_deg[IS1_all_deg >= quantile(IS1_all_deg,prob=1-n/100)]
length(IS1_all_deg.1percent) #317

IS1_all_betweenness.1percent <- IS1_all_betweenness[IS1_all_betweenness >= quantile(IS1_all_betweenness,prob=1-n/100)]
length(IS1_all_betweenness.1percent) #54

IS1_all_closeness.1percent <- IS1_all_closeness[IS1_all_closeness >= quantile(IS1_all_closeness,prob=1-n/100)]
length(IS1_all_closeness.1percent) #311

intersect(names(IS1_all_deg.1percent), names(IS1_all_betweenness.1percent))
intersect(names(IS1_all_deg.1percent), names(IS1_all_closeness.1percent))
intersect(names(IS1_all_betweenness.1percent), names(IS1_all_closeness.1percent))


## Set node shape soil_all_keystone.2017
V(all_soil_net.is1)$shape <- V(all_soil_net.is1)$name
V(all_soil_net.is1)$shape[V(all_soil_net.is1)$shape %in% names(IS1_all_betweenness.1percent)] <- "star"
#V(all_soil_net.is1)$shape[V(all_soil_net.is1)$shape %in% strictkeystone.is1] <- "star"
V(all_soil_net.is1)$shape[V(all_soil_net.is1)$shape %in% rownames(df.otu.bac.clean.nolog.is1)] <- "circle"
V(all_soil_net.is1)$shape[V(all_soil_net.is1)$shape %in% rownames(df.otu.fun.clean.nolog.is1)] <- "triangle"
V(all_soil_net.is1)$shape[V(all_soil_net.is1)$shape %in% rownames(df.otu.arch.clean.nolog.is1)] <- "square"

## Set node colors based kingdom
cs <- c("Bacteria", "Archaea", "Fungi")
unique(V(all_soil_net.is1)$kingdom)
V(all_soil_net.is1)$color <- V(all_soil_net.is1)$kingdom
V(all_soil_net.is1)$color[V(all_soil_net.is1)$color == "Bacteria"] <- "#EE7600"
V(all_soil_net.is1)$color[V(all_soil_net.is1)$color == "Archaea"] <- "#458B74"
V(all_soil_net.is1)$color[V(all_soil_net.is1)$color == "Fungi"] <- "#9A32CD"
V(all_soil_net.is1)$frame.color <- V(all_soil_net.is1)$color

#Node size
V(all_soil_net.is1)$size <- V(all_soil_net.is1)$name
V(all_soil_net.is1)$size[V(all_soil_net.is1)$size %in% names(IS1_all_betweenness.1percent)] <- 3
#V(all_soil_net.is1)$size[V(all_soil_net.is1)$size %in% strictkeystone.2017] <- 5
V(all_soil_net.is1)$size[V(all_soil_net.is1)$size %in% rownames(df.otu.bac.clean.nolog.is1)] <- 2
V(all_soil_net.is1)$size[V(all_soil_net.is1)$size %in% rownames(df.otu.fun.clean.nolog.is1)] <- 3
V(all_soil_net.is1)$size[V(all_soil_net.is1)$size %in% rownames(df.otu.arch.clean.nolog.is1)] <- 2
soilcombine_nodesizes.is1 <- as.numeric(V(all_soil_net.is1)$size)

## Set edge color
E(all_soil_net.is1)$color <- ifelse(E(all_soil_net.is1)$cor >0, "#99CCFF","#FF6666")


#Layout
source("star_shape.R")
source("triangle_shape.R")
library(qgraph)
all_soil_net.is1.integer <- set.vertex.attribute(all_soil_net.is1, "name", value=paste(1:5322))
edge.is1 <- get.edgelist(all_soil_net.is1.integer)
edge.is1.num<- apply(edge.is1, 2, as.numeric)

layout.is1 <- qgraph.layout.fruchtermanreingold(edge.is1.num,vcount=vcount(all_soil_net.is1.integer))

#Default
plot(all_soil_net.is1.integer,layout=layout.is1,vertex.label=NA, vertex.size=soilcombine_nodesizes.is1)

mtext("qgraph.layout.fruchtermanreingold default", side=1)

#Modified
layout.is1.modified <- qgraph.layout.fruchtermanreingold(edge.is1.num,vcount=vcount(all_soil_net.is1.integer),
                                                         area=10*(vcount(all_soil_net.is1.integer)^2),repulse.rad=(vcount(all_soil_net.is1.integer)^3.1))

#without module
plot(all_soil_net.is1.integer,layout=layout.is1.modified,vertex.size=soilcombine_nodesizes.is1,vertex.label=NA)

memory.limit(200000)

dev.off()

## IS2
all_cor_soil_df_padj.is2 <- read.table("all_cor_soil_df_padj.is2.txt", sep='\t', header = T)


nodeattrib_soil_combine.is2 <- data.frame(node=union(all_cor_soil_df_padj.is2$from,all_cor_soil_df_padj.is2$to))
nodeattrib_soil_combine.is2$kingdom <- 0

for (i in as.character(nodeattrib_soil_combine.is2$node))
{
  if (i %in% rownames(df.otu.bac.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.is2[nodeattrib_soil_combine.is2$node==i,"kingdom"] <- "Bacteria"}
  
  else if (i %in% rownames(df.otu.arch.clean.nolog) == TRUE)
  {nodeattrib_soil_combine.is2[nodeattrib_soil_combine.is2$node==i,"kingdom"] <- "Archaea"} 
  
  else
  {nodeattrib_soil_combine.is2[nodeattrib_soil_combine.is2$node==i,"kingdom"]<- "Fungi"}
}

rownames(nodeattrib_soil_combine.is2) <- as.character(nodeattrib_soil_combine.is2$node)
nodeattrib_soil_combine.is2$kingdom



all_soil_net.is2 <- graph_from_data_frame(all_cor_soil_df_padj.is2,direct=F, vertices=nodeattrib_soil_combine.is2)

## Number of nodes
length(V(all_soil_net.is2)) #5601

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.is2)))) #3482
length(grep("^F",names(V(all_soil_net.is2)))) #1913
length(grep("^A",names(V(all_soil_net.is2)))) #206

## Connections 
bb_occur_soil.is2 <- droplevels(all_cor_soil_df_padj.is2[with(all_cor_soil_df_padj.is2, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.is2) #238952

ff_occur_soil.is2 <- droplevels(all_cor_soil_df_padj.is2[with(all_cor_soil_df_padj.is2, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.is2) #73826

fb_occur_soil.is2 <- droplevels(all_cor_soil_df_padj.is2[with(all_cor_soil_df_padj.is2, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.is2) #214487

aa_occur_soil.is2 <- droplevels(all_cor_soil_df_padj.is2[with(all_cor_soil_df_padj.is2, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.is2) #973

ba_occur_soil.is2 <- droplevels(all_cor_soil_df_padj.is2[with(all_cor_soil_df_padj.is2, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.is2) #30169

fa_occur_soil.is2 <- droplevels(all_cor_soil_df_padj.is2[with(all_cor_soil_df_padj.is2, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.is2) #12799

## Network properties
meta_degree.is2 <- sort(igraph::degree(all_soil_net.is2,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.is2))) #"max degree = " "399"     
print(c('mean degree = ', mean(meta_degree.is2))) #"mean degree = "   "223.182421227197"
meta_degree.is2 <- as.data.frame(meta_degree.is2)
ggplot(meta_degree.is2, aes(x=meta_degree.is2)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.is2, directed = FALSE))) #"average shortest path length =  3.59523728700828"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.is2, "global"))) #"mean clustering coefficient =  0.949977913375762"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.is2, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000429103628639027"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.is2, normalized = TRUE)))) #"mean closeness centrality =  0.0735447050148934"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.is2, V(all_soil_net.is2)))))) #"mean number of neighbors =  223.182421227197"

##
net <- all_soil_net.is2
IS2_all_deg <- igraph::degree(net,mode="all")
IS2_all_betweenness <- betweenness(net, normalized = TRUE)
IS2_all_closeness <- closeness(net, normalized = TRUE)
IS2_all_transitivity <- transitivity(net, "local", vids = V(net))
names(IS2_all_transitivity)<- V(net)$name
IS2_all_transitivity[is.na(IS2_all_transitivity)] <- 0


## Defining hub OTUs
n <- 1
IS2_all_deg.1percent <- IS2_all_deg[IS2_all_deg >= quantile(IS2_all_deg,prob=1-n/100)]
length(IS2_all_deg.1percent) #305

IS2_all_betweenness.1percent <- IS2_all_betweenness[IS2_all_betweenness >= quantile(IS2_all_betweenness,prob=1-n/100)]
length(IS2_all_betweenness.1percent) #58

IS2_all_closeness.1percent <- IS2_all_closeness[IS2_all_closeness >= quantile(IS2_all_closeness,prob=1-n/100)]
length(IS2_all_closeness.1percent) #305

intersect(names(IS2_all_deg.1percent), names(IS2_all_betweenness.1percent))
intersect(names(IS2_all_deg.1percent), names(IS2_all_closeness.1percent))
intersect(names(IS2_all_betweenness.1percent), names(IS2_all_closeness.1percent))


## Set node shape soil_all_keystone.2017
V(all_soil_net.is2)$shape <- V(all_soil_net.is2)$name
V(all_soil_net.is2)$shape[V(all_soil_net.is2)$shape %in% names(IS2_all_betweenness.1percent)] <- "star"
#V(all_soil_net.is2)$shape[V(all_soil_net.is2)$shape %in% strictkeystone.is2] <- "star"
V(all_soil_net.is2)$shape[V(all_soil_net.is2)$shape %in% rownames(df.otu.bac.clean.nolog.is2)] <- "circle"
V(all_soil_net.is2)$shape[V(all_soil_net.is2)$shape %in% rownames(df.otu.fun.clean.nolog.is2)] <- "triangle"
V(all_soil_net.is2)$shape[V(all_soil_net.is2)$shape %in% rownames(df.otu.arch.clean.nolog.is2)] <- "square"

## Set node colors based kingdom
cs <- c("Bacteria", "Archaea", "Fungi")
unique(V(all_soil_net.is2)$kingdom)
V(all_soil_net.is2)$color <- V(all_soil_net.is2)$kingdom
V(all_soil_net.is2)$color[V(all_soil_net.is2)$color == "Bacteria"] <- "#EE7600"
V(all_soil_net.is2)$color[V(all_soil_net.is2)$color == "Archaea"] <- "#458B74"
V(all_soil_net.is2)$color[V(all_soil_net.is2)$color == "Fungi"] <- "#9A32CD"
V(all_soil_net.is2)$frame.color <- V(all_soil_net.is2)$color

#Node size
V(all_soil_net.is2)$size <- V(all_soil_net.is2)$name
V(all_soil_net.is2)$size[V(all_soil_net.is2)$size %in% names(IS2_all_betweenness.1percent)] <- 3
#V(all_soil_net.is2)$size[V(all_soil_net.is2)$size %in% strictkeystone.2017] <- 5
V(all_soil_net.is2)$size[V(all_soil_net.is2)$size %in% rownames(df.otu.bac.clean.nolog.is2)] <- 2
V(all_soil_net.is2)$size[V(all_soil_net.is2)$size %in% rownames(df.otu.fun.clean.nolog.is2)] <- 3
V(all_soil_net.is2)$size[V(all_soil_net.is2)$size %in% rownames(df.otu.arch.clean.nolog.is2)] <- 2
soilcombine_nodesizes.is2 <- as.numeric(V(all_soil_net.is2)$size)

## Set edge color
E(all_soil_net.is2)$color <- ifelse(E(all_soil_net.is2)$cor >0, "#99CCFF","#FF6666")


#Layout
source("star_shape.R")
source("triangle_shape.R")
library(qgraph)
all_soil_net.is2.integer <- set.vertex.attribute(all_soil_net.is2, "name", value=paste(1:5601))
edge.is2 <- get.edgelist(all_soil_net.is2.integer)
edge.is2.num<- apply(edge.is2, 2, as.numeric)

layout.is2 <- qgraph.layout.fruchtermanreingold(edge.is2.num,vcount=vcount(all_soil_net.is2.integer))

#Default
plot(all_soil_net.is2.integer,layout=layout.is2,vertex.label=NA, vertex.size=soilcombine_nodesizes.is2)

mtext("qgraph.layout.fruchtermanreingold default", side=1)

#Modified
layout.is2.modified <- qgraph.layout.fruchtermanreingold(edge.is2.num,vcount=vcount(all_soil_net.is2.integer),
                                                         area=10*(vcount(all_soil_net.is2.integer)^2),repulse.rad=(vcount(all_soil_net.is2.integer)^3.1))

#without module
plot(all_soil_net.is2.integer,layout=layout.is2.modified,vertex.size=soilcombine_nodesizes.is2,vertex.label=NA)

memory.limit(200000)

dev.off()



bb_occur_soil.is1.neg<-subset(bb_occur_soil.is1, bb_occur_soil.is1$cor < 0)
nrow(bb_occur_soil.is1.neg) #107

ff_occur_soil.is1.neg<-subset(ff_occur_soil.is1, ff_occur_soil.is1$cor < 0)
nrow(ff_occur_soil.is1.neg) #3629

aa_occur_soil.is1.neg<-subset(aa_occur_soil.is1, aa_occur_soil.is1$cor < 0)
nrow(aa_occur_soil.is1.neg) #0

ba_occur_soil.is1.neg<-subset(ba_occur_soil.is1, ba_occur_soil.is1$cor < 0)
nrow(ba_occur_soil.is1.neg) #20

fa_occur_soil.is1.neg<-subset(fa_occur_soil.is1, fa_occur_soil.is1$cor < 0)
nrow(fa_occur_soil.is1.neg) #9

fb_occur_soil.is1.neg<-subset(fb_occur_soil.is1, fb_occur_soil.is1$cor < 0)
nrow(fb_occur_soil.is1.neg) #836

neg.occur.is1<-rbind(bb_occur_soil.is1.neg, ff_occur_soil.is1.neg, aa_occur_soil.is1.neg,ba_occur_soil.is1.neg,fa_occur_soil.is1.neg,fb_occur_soil.is1.neg)
write.xlsx(neg.occur.is1, "Negative correlation in is1 network_with archaea.xlsx")


bb_occur_soil.is2.neg<-subset(bb_occur_soil.is2, bb_occur_soil.is2$cor < 0)
nrow(bb_occur_soil.is2.neg) #102

ff_occur_soil.is2.neg<-subset(ff_occur_soil.is2, ff_occur_soil.is2$cor < 0)
nrow(ff_occur_soil.is2.neg) #4477

aa_occur_soil.is2.neg<-subset(aa_occur_soil.is2, aa_occur_soil.is2$cor < 0)
nrow(aa_occur_soil.is2.neg) #3

ba_occur_soil.is2.neg<-subset(ba_occur_soil.is2, ba_occur_soil.is2$cor < 0)
nrow(ba_occur_soil.is2.neg) #39

fa_occur_soil.is2.neg<-subset(fa_occur_soil.is2, fa_occur_soil.is2$cor < 0)
nrow(fa_occur_soil.is2.neg) #185

fb_occur_soil.is2.neg<-subset(fb_occur_soil.is2, fb_occur_soil.is2$cor < 0)
nrow(fb_occur_soil.is2.neg) #693

neg.occur.is2<-rbind(bb_occur_soil.is2.neg, ff_occur_soil.is2.neg, aa_occur_soil.is2.neg,ba_occur_soil.is2.neg,fa_occur_soil.is2.neg,fb_occur_soil.is2.neg)
write.xlsx(neg.occur.is2, "Negative correlation in is2 network_with archaea.xlsx")
