##### Identifiying cropping system responsive OTUs with indicator species analysis #####

## Identify indicator species in bacterial community
otu.bac.clean.nolog <- otu_table(bac.clean.nolog)
df.otu.bac.clean.nolog <- data.frame(otu.bac.clean.nolog)
df.otu.bac.clean.nolog$OTU <-rownames(df.otu.bac.clean.nolog)
df.otu.bac.clean.nolog <- merge(df.otu.bac.clean.nolog, OTU_id.list, by = "OTU")
rownames(df.otu.bac.clean.nolog) <- df.otu.bac.clean.nolog$OTU_id
df.otu.bac.clean.nolog <- df.otu.bac.clean.nolog[-c(1,254)]

otu.arch.clean.nolog <- otu_table(arch.clean.nolog)
df.otu.arch.clean.nolog <- data.frame(otu.arch.clean.nolog)
df.otu.arch.clean.nolog$OTU <-rownames(df.otu.arch.clean.nolog)
df.otu.arch.clean.nolog <- merge(df.otu.arch.clean.nolog, OTU_id.list, by = "OTU")
rownames(df.otu.arch.clean.nolog) <- df.otu.arch.clean.nolog$OTU_id
df.otu.arch.clean.nolog <- df.otu.arch.clean.nolog[-c(1,254)]

otu.fun.clean.nolog <- otu_table(fun.clean.nolog)
df.otu.fun.clean.nolog <- data.frame(otu.fun.clean.nolog)
df.otu.fun.clean.nolog$OTU <-rownames(df.otu.fun.clean.nolog)
df.otu.fun.clean.nolog <- merge(df.otu.fun.clean.nolog, OTU_id.list, by = "OTU")
rownames(df.otu.fun.clean.nolog) <- df.otu.fun.clean.nolog$OTU_id
df.otu.fun.clean.nolog <- df.otu.fun.clean.nolog[-c(1,254)]

#bacteria
otu_norm_soil_16s <- df.otu.bac.clean.nolog
design_16s_soil <- meta(bac.clean.nolog)

indic_soil_16s <- as.data.frame(t(otu_norm_soil_16s))
indic_soil_groups_16s <- design_16s_soil$Cultural_practice
length(unique(indic_soil_groups_16s))

#archaea
otu_norm_soil_arch <- df.otu.bac.clean.nolog
design_arch_soil <- meta(bac.clean.nolog)

indic_soil_arch <- as.data.frame(t(otu_norm_soil_arch))
indic_soil_groups_arch <- design_arch_soil$Cultural_practice
length(unique(indic_soil_groups_arch))


#Fungi
otu_norm_soil_its <- df.otu.bac.clean.nolog
design_its_soil <- meta(bac.clean.nolog)

indic_soil_its <- as.data.frame(t(otu_norm_soil_its))
indic_soil_groups_its <- design_its_soil$Cultural_practice
length(unique(indic_soil_groups_its))


## Define indicator species for soil bacteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_16s <- multipatt(indic_soil_16s,indic_soil_groups_16s,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_16s,alpha=1,indvalcomp=T)

indicatorsp_soil_arch <- multipatt(indic_soil_arch,indic_soil_groups_arch,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_arch,alpha=1,indvalcomp=T)

indicatorsp_soil_its <- multipatt(indic_soil_its,indic_soil_groups_its,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_its,alpha=1,indvalcomp=T)


indic_soil_df_16s <- indicatorsp_soil_16s$sign
indic_soil_df_arch<- indicatorsp_soil_arch$sign
indic_soil_df_its<- indicatorsp_soil_its$sign

write.table(indic_soil_df_16s,"indicsp_soil_16s.txt", sep="\t",quote=F)
write.table(indic_soil_df_arch,"indicsp_soil_arch.txt", sep="\t",quote=F)
write.table(indic_soil_df_its,"indicsp_soil_its.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_soil_df_16s <- read.table("indicsp_soil_16s.txt", header=T, sep="\t")
indic_soil_df_arch<- read.table("indicsp_soil_arch.txt", header=T, sep="\t")
indic_soil_df_its<- read.table("indicsp_soil_its.txt", header=T, sep="\t")

#Bacteria
Conventional_soil_16s <- as.matrix(indic_soil_df_16s[which(indic_soil_df_16s$s.Conventional == 1 & indic_soil_df_16s$p.value < 0.05),])
No_fertilizer_soil_16s <- as.matrix(indic_soil_df_16s[which(indic_soil_df_16s$s.No_fertilizer == 1 & indic_soil_df_16s$p.value < 0.05),])
No_pesticide_soil_16s <- as.matrix(indic_soil_df_16s[which(indic_soil_df_16s$s.No_pesticide == 1 & indic_soil_df_16s$p.value < 0.05),])

soil_r_values_16s <- rbind(Conventional_soil_16s,No_fertilizer_soil_16s,No_pesticide_soil_16s)
colnames(soil_r_values_16s)[1:3] <-c("Conventional","No_fertilizer","No_pesticide")
write.xlsx(soil_r_values_16s, "Bac_indicator OTU.xlsx")

#Archaea
Conventional_soil_arch <- as.matrix(indic_soil_df_arch[which(indic_soil_df_arch$s.Conventional == 1 & indic_soil_df_arch$p.value < 0.05),])
No_fertilizer_soil_arch <- as.matrix(indic_soil_df_arch[which(indic_soil_df_arch$s.No_fertilizer == 1 & indic_soil_df_arch$p.value < 0.05),])
No_pesticide_soil_arch <- as.matrix(indic_soil_df_arch[which(indic_soil_df_arch$s.No_pesticide == 1 & indic_soil_df_arch$p.value < 0.05),])

soil_r_values_arch <- rbind(Conventional_soil_arch,No_fertilizer_soil_arch,No_pesticide_soil_arch)
colnames(soil_r_values_arch)[1:3] <-c("Conventional","No_fertilizer","No_pesticide")
write.xlsx(soil_r_values_arch, "Arch_indicator OTU.xlsx")

#Fungi
Conventional_soil_its <- as.matrix(indic_soil_df_its[which(indic_soil_df_its$s.Conventional == 1 & indic_soil_df_its$p.value < 0.05),])
No_fertilizer_soil_its <- as.matrix(indic_soil_df_its[which(indic_soil_df_its$s.No_fertilizer == 1 & indic_soil_df_its$p.value < 0.05),])
No_pesticide_soil_its <- as.matrix(indic_soil_df_its[which(indic_soil_df_its$s.No_pesticide == 1 & indic_soil_df_its$p.value < 0.05),])

soil_r_values_its <- rbind(Conventional_soil_its,No_fertilizer_soil_its,No_pesticide_soil_its)
colnames(soil_r_values_its)[1:3] <-c("Conventional","No_fertilizer","No_pesticide")
write.xlsx(soil_r_values_its, "Fun_indicator OTU.xlsx")

## Range of correlation coefficients
range(soil_r_values_16s[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(soil_r_values_16s)))

## Proportion of soil bacteria OTUs responding to cultural practices
length(unique(rownames(soil_r_values_16s)))/nrow(otu_norm_soil_16s)



## Bipartite network of indicator OTUs
##Extract taxonomy table
phylo_to_taxidtable <- function(mer.clean.ss.5,file_name){
  otutable_mer.clean.ss.5 <- as.data.frame(tax_table(mer.clean.ss.5))
  rownames(otutable_mer.clean.ss.5)
  
  otutable_mer.clean.ss.5 <- rownames_to_column(otutable_mer.clean.ss.5, var="otu_id")
  (otutable_mer.clean.ss.5$otu_id)
  otutable_mer.clean.ss.5 <- left_join(otutable_mer.clean.ss.5,OTU_id.list, by=c('otu_id'='OTU')) %>% select(-otu_id) %>% select(OTU_id,everything()) %>% arrange(OTU_id)
  write.table(otutable_mer.clean.ss.5, file=file_name, quote=FALSE, sep='\t', row.names = F)
}

#bacteria
phylo_to_taxidtable(bac.clean.ss, "tax_soil_16s.txt")
df.tax_soil_16s <- read.table("tax_soil_16s.txt", sep='\t', header=T, stringsAsFactors = F)
rownames(df.tax_soil_16s)<- df.tax_soil_16s$OTU_id
df.tax_soil_16s <- df.tax_soil_16s[-c(1)]


#archaea
phylo_to_taxidtable(arch.clean.ss, "tax_soil_arch.txt")
df.tax_soil_arch <- read.table("tax_soil_arch.txt", sep='\t', header=T, stringsAsFactors = F)
rownames(df.tax_soil_arch)<- df.tax_soil_arch$OTU_id
df.tax_soil_arch <- df.tax_soil_arch[-c(1)]

#Fungi
phylo_to_taxidtable(fun.clean.ss, "tax_soil_its.txt")
df.tax_soil_its <- read.table("tax_soil_its.txt", sep='\t', header=T, stringsAsFactors = F)
rownames(df.tax_soil_its)<- df.tax_soil_its$OTU_id
df.tax_soil_its <- df.tax_soil_its[-c(1)]

### Bipartite network of Bacteria
## Construct node table for bulk soil bacteria communities from indicator OTU data
soil_bipartite_16s <- data.frame(from= c(rep("Conventional",length(which(soil_r_values_16s[,"Conventional"]==1))),
                                         rep("No_fertilizer",length(which(soil_r_values_16s[,"No_fertilizer"]==1))),
                                         rep("No_pesticide",length(which(soil_r_values_16s[,"No_pesticide"]==1)))),
                                 to= c(rownames(soil_r_values_16s)[which(soil_r_values_16s[,"Conventional"]==1)],
                                       rownames(soil_r_values_16s)[which(soil_r_values_16s[,"No_fertilizer"]==1)],
                                       rownames(soil_r_values_16s)[which(soil_r_values_16s[,"No_pesticide"]==1)]),
                                 r= c(soil_r_values_16s[which(soil_r_values_16s[,"Conventional"]==1),"stat"],
                                      soil_r_values_16s[which(soil_r_values_16s[,"No_fertilizer"]==1),"stat"],
                                      soil_r_values_16s[which(soil_r_values_16s[,"No_pesticide"]==1),"stat"]))

## make node attribute table for each OTU
soil_bipartite_attrib_16s <- data.frame(node=unique(rownames(soil_r_values_16s)),indicgroup=0)

for (i in as.character(soil_bipartite_attrib_16s$node))
{
  soil_bipartite_attrib_16s[soil_bipartite_attrib_16s$node==i,"indicgroup"] <- paste(colnames(soil_r_values_16s)[which(soil_r_values_16s[i,1:4]==1)],collapse = "_")
}

soil_bipartite_attrib_16s <- cbind(soil_bipartite_attrib_16s,df.tax_soil_16s[as.character(soil_bipartite_attrib_16s$node),])

## Create bipartite network with igraph
soil_bi_16s <- graph.data.frame(soil_bipartite_16s,directed=F)
V(soil_bi_16s)$type <- V(soil_bi_16s)$name %in% soil_bipartite_16s[,1]
soil_bi_16s <- simplify(soil_bi_16s, remove.multiple=T, remove.loops=T)

## Set labels for nodes
V(soil_bi_16s)$label <- V(soil_bi_16s)$name
V(soil_bi_16s)$label <- gsub("B^",NA,V(soil_bi_16s)$label)

## Set node sizes
V(soil_bi_16s)$size <- c(rep(8,4),rep(3,1430))

## Set node shapes
V(soil_bi_16s)$shape <- c(rep("circle",3),rep("circle",1431))

## Define node colors based upon phylum/class taxonomy assignment
tax_16s <- df.tax_soil_16s
head(tax_16s)

tax_16s[is.na(tax_16s)] <- "Unassigned"
tax_16s[tax_16s=="uncultured bacterium"] <- "Unassigned"
tax_16s[tax_16s=="uncultured soil bacterium"] <- "Unassigned"
tax_16s[tax_16s=="uncultured prokaryote"] <- "Unassigned"
#colnames(tax_16s) <- c("Kingdom","Phylum","Class","Order", "Family", "Genus", "Species","C")
is.na(tax_16s)
# create separate taxonomy label specifying classes of Proteobacteria
tax_16s$labels <- tax_16s$Phylum
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Alphaproteobacteria" ], ]$labels <- "Alphaproteobacteria"
#tax_16s[ rownames(tax_16s)[tax_16s$Class=="Betaproteobacteria" ], ]$labels <- "Betaproteobacteria"
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Gammaproteobacteria" ], ]$labels <- "Gammaproteobacteria"
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Deltaproteobacteria" ], ]$labels <- "Deltaproteobacteria"
table(tax_16s$labels)


#Defining Legend
tax_16s$cols <- tax_16s$labels
table(tax_16s$cols)


## Express OTU counts as relative abunance percent
otu.bac.clean.ss <- otu_table(bac.clean.ss)
df.otu.bac.clean.ss <- data.frame(otu.bac.clean.ss)
df.otu.bac.clean.ss$OTU <- rownames(df.otu.bac.clean.ss)
df.otu.bac.clean.ss <- merge(df.otu.bac.clean.ss, OTU_id.list, by ="OTU")
rownames(df.otu.bac.clean.ss)<- df.otu.bac.clean.ss$OTU_id
colnames(df.otu.bac.clean.ss)
df.otu.bac.clean.ss <- df.otu.bac.clean.ss[-c(1, 254)]

otu_16s <- df.otu.bac.clean.ss

otu_16s_RA <- t(t(otu_16s)/colSums(otu_16s)) * 100

colSums(otu_16s_RA)

nrow(otu_16s_RA)

## Get names of bacteria phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_16s <- names(sort(table(tax_16s[,"labels"]), decr=T))
length(PHYLAnames_16s)
sort(table(tax_16s[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_16s_RA)
for (i in PHYLAnames_16s){
  x <- array(colSums(otu_16s_RA[rownames(tax_16s)[which(tax_16s$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_16s)
colnames(y) <- paste(colnames(otu_16s_RA))
PHYLUM_mat_16s <- y
PHYLUM_mat_16s[,1:5]
colSums(PHYLUM_mat_16s)
PHYLUM_mat_16s_mean <- sort(apply(PHYLUM_mat_16s,1,mean),decr=T)
PHYLUM_mat_16s <- PHYLUM_mat_16s[names(PHYLUM_mat_16s_mean),]


# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_16s, 1, mean) < 1)
low_count_phyla_16s <- rownames(PHYLUM_mat_16s)[sort(apply(PHYLUM_mat_16s, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_16s){
  tax_16s[ rownames(tax_16s)[tax_16s$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_16s$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_16s <- rownames(PHYLUM_mat_16s)[sort(apply(PHYLUM_mat_16s, 1, mean), decr=T) > 1]
abundant_phyla_16s

# Color scheme
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Alphaproteobacteria" ], ]$cols <- "palegreen1"
#tax_16s[ rownames(tax_16s)[tax_16s$labels=="Betaproteobacteria" ], ]$cols <- "palegreen3"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Gammaproteobacteria" ], ]$cols <- "palegreen2"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Deltaproteobacteria" ], ]$cols <- "palegreen3"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Actinobacteria" ], ]$cols <- "indianred2"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Bacteroidetes" ], ]$cols <- "steelblue1"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Firmicutes" ], ]$cols <- "tan1"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Acidobacteria" ], ]$cols <- "lightsalmon4"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Chloroflexi" ], ]$cols <- "gold1"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Verrucomicrobia" ], ]$cols <- "orchid3"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Nitrospirae" ], ]$cols <- "palevioletred2"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Gemmatimonadetes" ], ]$cols <- "peachpuff3"
#tax_16s[ rownames(tax_16s)[tax_16s$labels=="Euryarchaeota" ], ]$cols <- "#7B55E1"
#tax_16s[ rownames(tax_16s)[tax_16s$labels=="Crenarchaeota" ], ]$cols <- "#FF2F2F"
tax_16s[ rownames(tax_16s)[tax_16s$labels=="Planctomycetes" ], ]$cols <- "khaki3"


## collaps OTU colors to prepare Phylum level colors
label_cols_16s <- tax_16s[, c("labels", "cols") ]
#library(plyr)
PHYLA_label_cols_16s <- plyr::ddply(label_cols_16s, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_16s) <- PHYLA_label_cols_16s[,1]
PHYLA_label_cols_16s <- PHYLA_label_cols_16s[c(abundant_phyla_16s, low_count_phyla_16s),]
PHYLA_label_cols_16s

## Legend for Phylum colors
PHYLA_label_cols_16s_legend <- PHYLA_label_cols_16s[1:13,]
PHYLA_label_cols_16s_legend[13,1] <- "Other"
rownames(PHYLA_label_cols_16s_legend)[13] <- "Other"


V(soil_bi_16s)$color <- V(soil_bi_16s)$name
V(soil_bi_16s)$color[1:3] <- "white"
V(soil_bi_16s)$color <- tax_16s[ V(soil_bi_16s)$name, ]$cols
V(soil_bi_16s)$frame.color <- V(soil_bi_16s)$color



### plotting
set.seed(8046)
soil_layout_16s <- layout_with_fr(soil_bi_16s, niter=9999)


pdf(paste0(output,"Figure3.pdf"), encoding="MacRoman", height=6, width=6)
layout(matrix(c(1,2,3,4,5,6), nrow=3, byrow=T), c(2,2), c(4,4,3))

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(soil_bi_16s, vertex.label.cex=F, layout=soil_layout_16s, asp=0)


plot(1, type="n", ann=F, axes=F)
legend("center", pch=19, bty="n", ncol=2, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_16s_legend$labels, 
       col=PHYLA_label_cols_16s_legend$cols)

plot(1,type="n",ann=F,axes=F)
legend("center", pch=19, bty="n", ncol=1, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_its_legend$Phylum, 
       col=PHYLA_label_cols_its_legend$cols)

dev.off()


### Archaea
## Construct node table for bulk soil bacteria communities from indicator species data
soil_bipartite_arch <- data.frame(from= c(rep("Conventional",length(which(soil_r_values_arch[,"Conventional"]==1))),
                                         rep("No_fertilizer",length(which(soil_r_values_arch[,"No_fertilizer"]==1))),
                                         rep("No_pesticide",length(which(soil_r_values_arch[,"No_pesticide"]==1)))),
                                 to= c(rownames(soil_r_values_arch)[which(soil_r_values_arch[,"Conventional"]==1)],
                                       rownames(soil_r_values_arch)[which(soil_r_values_arch[,"No_fertilizer"]==1)],
                                       rownames(soil_r_values_arch)[which(soil_r_values_arch[,"No_pesticide"]==1)]),
                                 r= c(soil_r_values_arch[which(soil_r_values_arch[,"Conventional"]==1),"stat"],
                                      soil_r_values_arch[which(soil_r_values_arch[,"No_fertilizer"]==1),"stat"],
                                      soil_r_values_arch[which(soil_r_values_arch[,"No_pesticide"]==1),"stat"]))

## make node attribute table for each OTU
soil_bipartite_attrib_arch <- data.frame(node=unique(rownames(soil_r_values_arch)),indicgroup=0)

for (i in as.character(soil_bipartite_attrib_arch$node))
{
  soil_bipartite_attrib_arch[soil_bipartite_attrib_arch$node==i,"indicgroup"] <- paste(colnames(soil_r_values_arch)[which(soil_r_values_arch[i,1:4]==1)],collapse = "_")
}

soil_bipartite_attrib_arch <- cbind(soil_bipartite_attrib_arch,df.tax_soil_arch[as.character(soil_bipartite_attrib_arch$node),])

## Create bipartite network with igraph
soil_bi_arch <- graph.data.frame(soil_bipartite_arch,directed=F)
V(soil_bi_arch)$type <- V(soil_bi_arch)$name %in% soil_bipartite_arch[,1]
soil_bi_arch <- simplify(soil_bi_arch, remove.multiple=T, remove.loops=T)

## Set labels for nodes
V(soil_bi_arch)$label <- V(soil_bi_arch)$name
V(soil_bi_arch)$label <- gsub("A^",NA,V(soil_bi_arch)$label)

## Set node sizes
V(soil_bi_arch)$size <- c(rep(8,4),rep(3,84))

## Set node shapes
V(soil_bi_arch)$shape <- c(rep("square",3),rep("square",85))

## Define node colors based upon phylum/class taxonomy assignment
tax_arch <- df.tax_soil_arch
head(tax_arch)

tax_arch[is.na(tax_arch)] <- "Unassigned"
#colnames(tax_arch) <- c("Kingdom","Phylum","Class","Order", "Family", "Genus", "Species","C")
is.na(tax_arch)
# create separate taxonomy label specifying classes of Proteobacteria
tax_arch$labels <- tax_arch$Phylum
table(tax_arch$labels)


#Defining Legend
tax_arch$cols <- tax_arch$labels
table(tax_arch$cols)


## Express arch OTU counts as relative abunance percent
otu.arch.clean.ss <- otu_table(arch.clean.ss)
df.otu.arch.clean.ss <- data.frame(otu.arch.clean.ss)
df.otu.arch.clean.ss$OTU <- rownames(df.otu.arch.clean.ss)
df.otu.arch.clean.ss <- merge(df.otu.arch.clean.ss, OTU_id.list, by ="OTU")
rownames(df.otu.arch.clean.ss)<- df.otu.arch.clean.ss$OTU_id
colnames(df.otu.arch.clean.ss)
df.otu.arch.clean.ss <- df.otu.arch.clean.ss[-c(1, 254)]

otu_arch <- df.otu.arch.clean.ss

otu_arch_RA <- t(t(otu_arch)/colSums(otu_arch)) * 100

colSums(otu_arch_RA)

nrow(otu_arch_RA)

## Get names of bacteria phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_arch <- names(sort(table(tax_arch[,"labels"]), decr=T))
length(PHYLAnames_arch)
sort(table(tax_arch[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_arch_RA)
for (i in PHYLAnames_arch){
  x <- array(colSums(otu_arch_RA[rownames(tax_arch)[which(tax_arch$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_arch)
colnames(y) <- paste(colnames(otu_arch_RA))
PHYLUM_mat_arch <- y
PHYLUM_mat_arch[,1:5]
colSums(PHYLUM_mat_arch)
PHYLUM_mat_arch_mean <- sort(apply(PHYLUM_mat_arch,1,mean),decr=T)
PHYLUM_mat_arch <- PHYLUM_mat_arch[names(PHYLUM_mat_arch_mean),]


# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_arch, 1, mean) < 1)
low_count_phyla_arch <- rownames(PHYLUM_mat_arch)[sort(apply(PHYLUM_mat_arch, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_arch){
  tax_arch[ rownames(tax_arch)[tax_arch$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_arch$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_arch <- rownames(PHYLUM_mat_arch)[sort(apply(PHYLUM_mat_arch, 1, mean), decr=T) > 1]
abundant_phyla_arch


# Color scheme for archaea
tax_arch[ rownames(tax_arch)[tax_arch$labels=="Euryarchaeota" ], ]$cols <- "#7B55E1"
tax_arch[ rownames(tax_arch)[tax_arch$labels=="Crenarchaeota" ], ]$cols <- "#FF9999"
tax_arch[ rownames(tax_arch)[tax_arch$labels=="Thaumarchaeota" ], ]$cols <- "#CC6600"
tax_arch[ rownames(tax_arch)[tax_arch$labels=="Nanoarchaeaeota" ], ]$cols <- "#999933"

## collaps OTU colors to prepare Phylum level colors
label_cols_arch <- tax_arch[, c("labels", "cols") ]
#library(plyr)
PHYLA_label_cols_arch <- plyr::ddply(label_cols_arch, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_arch) <- PHYLA_label_cols_arch[,1]
PHYLA_label_cols_arch <- PHYLA_label_cols_arch[c(abundant_phyla_arch, low_count_phyla_arch),]
PHYLA_label_cols_arch

## Legend for Phylum colors
PHYLA_label_cols_arch_legend <- PHYLA_label_cols_arch[1:5,]
PHYLA_label_cols_arch_legend[5,1] <- "Other"
rownames(PHYLA_label_cols_arch_legend)[5] <- "Other"


V(soil_bi_arch)$color <- V(soil_bi_arch)$name
V(soil_bi_arch)$color[1:3] <- "white"
V(soil_bi_arch)$color <- tax_arch[ V(soil_bi_arch)$name, ]$cols
V(soil_bi_arch)$frame.color <- V(soil_bi_arch)$color

### plotting
set.seed(8046)
soil_layout_arch <- layout_with_fr(soil_bi_arch, niter=9999)


pdf(paste0(output,"Figure3.pdf"), encoding="MacRoman", height=6, width=6)
layout(matrix(c(1,2,3,4,5,6), nrow=3, byrow=T), c(2,2), c(4,4,3))

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(soil_bi_arch, vertex.label.cex=F, layout=soil_layout_arch, asp=0)


plot(1, type="n", ann=F, axes=F)
legend("center", pch=19, bty="n", ncol=2, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_arch_legend$labels, 
       col=PHYLA_label_cols_arch_legend$cols)

plot(1,type="n",ann=F,axes=F)
legend("center", pch=19, bty="n", ncol=1, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_its_legend$Phylum, 
       col=PHYLA_label_cols_its_legend$cols)

dev.off()


### Fungi
## Construct node table for bulk soil bacteria communities from indicator species data
soil_bipartite_its <- data.frame(from= c(rep("Conventional",length(which(soil_r_values_its[,"Conventional"]==1))),
                                          rep("No_fertilizer",length(which(soil_r_values_its[,"No_fertilizer"]==1))),
                                          rep("No_pesticide",length(which(soil_r_values_its[,"No_pesticide"]==1)))),
                                  to= c(rownames(soil_r_values_its)[which(soil_r_values_its[,"Conventional"]==1)],
                                        rownames(soil_r_values_its)[which(soil_r_values_its[,"No_fertilizer"]==1)],
                                        rownames(soil_r_values_its)[which(soil_r_values_its[,"No_pesticide"]==1)]),
                                  r= c(soil_r_values_its[which(soil_r_values_its[,"Conventional"]==1),"stat"],
                                       soil_r_values_its[which(soil_r_values_its[,"No_fertilizer"]==1),"stat"],
                                       soil_r_values_its[which(soil_r_values_its[,"No_pesticide"]==1),"stat"]))

## make node attribute table for each OTU
soil_bipartite_attrib_its <- data.frame(node=unique(rownames(soil_r_values_its)),indicgroup=0)

for (i in as.character(soil_bipartite_attrib_its$node))
{
  soil_bipartite_attrib_its[soil_bipartite_attrib_its$node==i,"indicgroup"] <- paste(colnames(soil_r_values_its)[which(soil_r_values_its[i,1:4]==1)],collapse = "_")
}

soil_bipartite_attrib_its <- cbind(soil_bipartite_attrib_its,df.tax_soil_its[as.character(soil_bipartite_attrib_its$node),])

## Create bipartite network with igraph
soil_bi_its <- graph.data.frame(soil_bipartite_its,directed=F)
V(soil_bi_its)$type <- V(soil_bi_its)$name %in% soil_bipartite_its[,1]
soil_bi_its <- simplify(soil_bi_its, remove.multiple=T, remove.loops=T)

## Set labels for nodes
V(soil_bi_its)$label <- V(soil_bi_its)$name
V(soil_bi_its)$label <- gsub("F^",NA,V(soil_bi_its)$label)

## Set node sizes
V(soil_bi_its)$size <- c(rep(8,4),rep(3,1061))

## Set node shapes
V(soil_bi_its)$shape <- c(rep("triangle",3),rep("triangle",1062))

## Define node colors based upon phylum/class taxonomy assignment
tax_its <- df.tax_soil_its
head(tax_its)

tax_its[is.na(tax_its)] <- "Unassigned"
tax_its[tax_its=="unidentified"] <- "Unassigned"
#colnames(tax_its) <- c("Kingdom","Phylum","Class","Order", "Family", "Genus", "Species","C")
is.na(tax_its)
# create separate taxonomy label specifying classes of Proteobacteria
tax_its$labels <- tax_its$Phylum
table(tax_its$labels)


#Defining Legend
tax_its$cols <- tax_its$labels
table(tax_its$cols)


## Express its OTU counts as relative abunance percent
otu_its <- df.fun.clean.ss

otu_its_RA <- t(t(otu_its)/colSums(otu_its)) * 100

colSums(otu_its_RA)

nrow(otu_its_RA)

## Get names of bacteria phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_its <- names(sort(table(tax_its[,"labels"]), decr=T))
length(PHYLAnames_its)
sort(table(tax_its[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_its_RA)
for (i in PHYLAnames_its){
  x <- array(colSums(otu_its_RA[rownames(tax_its)[which(tax_its$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_its)
colnames(y) <- paste(colnames(otu_its_RA))
PHYLUM_mat_its <- y
PHYLUM_mat_its[,1:5]
colSums(PHYLUM_mat_its)
PHYLUM_mat_its_mean <- sort(apply(PHYLUM_mat_its,1,mean),decr=T)
PHYLUM_mat_its <- PHYLUM_mat_its[names(PHYLUM_mat_its_mean),]


# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_its, 1, mean) < 1)
low_count_phyla_its <- rownames(PHYLUM_mat_its)[sort(apply(PHYLUM_mat_its, 1, mean), decr=T) < 1]
# attribute grey color
for(i in low_count_phyla_its){
  tax_its[ rownames(tax_its)[tax_its$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_its$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_its <- rownames(PHYLUM_mat_its)[sort(apply(PHYLUM_mat_its, 1, mean), decr=T) > 1]
abundant_phyla_its


# Color scheme for fungi
tax_its[ rownames(tax_its)[tax_its$labels=="Ascomycota" ], ]$cols <- "#003366"
tax_its[ rownames(tax_its)[tax_its$labels=="Basidiomycota" ], ]$cols <- "#C84248"
tax_its[ rownames(tax_its)[tax_its$labels=="Unassigned" ], ]$cols <- "#696969"
tax_its[ rownames(tax_its)[tax_its$labels=="Mortierellomycota" ], ]$cols <- "#EEB422"

## collaps OTU colors to prepare Phylum level colors
label_cols_its <- tax_its[, c("labels", "cols") ]
#library(plyr)
PHYLA_label_cols_its <- plyr::ddply(label_cols_its, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_its) <- PHYLA_label_cols_its[,1]
PHYLA_label_cols_its <- PHYLA_label_cols_its[c(abundant_phyla_its, low_count_phyla_its),]
PHYLA_label_cols_its

## Legend for Phylum colors
PHYLA_label_cols_its_legend <- PHYLA_label_cols_its[1:5,]
PHYLA_label_cols_its_legend[5,1] <- "Other"
rownames(PHYLA_label_cols_its_legend)[5] <- "Other"


V(soil_bi_its)$color <- V(soil_bi_its)$name
V(soil_bi_its)$color[1:3] <- "white"
V(soil_bi_its)$color <- tax_its[ V(soil_bi_its)$name, ]$cols
V(soil_bi_its)$frame.color <- V(soil_bi_its)$color

### plotting
set.seed(8046)
soil_layout_its <- layout_with_fr(soil_bi_its, niter=9999)


pdf(paste0(output,"Figure3.pdf"), encoding="MacRoman", height=6, width=6)
layout(matrix(c(1,2,3,4,5,6), nrow=3, byrow=T), c(2,2), c(4,4,3))

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(soil_bi_its, vertex.label.cex=F, layout=soil_layout_its, asp=0)


plot(1, type="n", ann=F, axes=F)
legend("center", pch=19, bty="n", ncol=2, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_its_legend$labels, 
       col=PHYLA_label_cols_its_legend$cols)

plot(1,type="n",ann=F,axes=F)
legend("center", pch=19, bty="n", ncol=1, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_its_legend$Phylum, 
       col=PHYLA_label_cols_its_legend$cols)

dev.off()




##rOTU
daOTU_conven_bac.all #2402
daOTU_nofertil_bac.all #6500
daOTU_nopesti_bac.all #569

daOTU_conven_fun.all #331
daOTU_nofertil_fun.all #1964
daOTU_nopesti_fun.all #1152

daOTU_conven_arch.all #109
daOTU_nofertil_arch.all #328
daOTU_nopesti_arch.all #9


soil_r_values_16s.CF <- subset(soil_r_values_16s, rownames(soil_r_values_16s)%in%daOTU_conven_bac.all) 
length(rownames(soil_r_values_16s.CF)) #420
soil_r_values_16s.NF <- subset(soil_r_values_16s, rownames(soil_r_values_16s)%in%daOTU_nofertil_bac.all) 
length(rownames(soil_r_values_16s.NF)) #732
soil_r_values_16s.NP <- subset(soil_r_values_16s, rownames(soil_r_values_16s)%in%daOTU_nopesti_bac.all) 
length(rownames(soil_r_values_16s.NP)) #81

soil_r_values_16s.rOTU <- rbind(soil_r_values_16s.CF,soil_r_values_16s.NF,soil_r_values_16s.NP)
write.xlsx(soil_r_values_16s.rOTU, "Bac_rOTU.xlsx")

soil_r_values_arch.CF <- subset(soil_r_values_arch, rownames(soil_r_values_arch)%in%daOTU_conven_arch.all) 
length(rownames(soil_r_values_arch.CF)) #30
soil_r_values_arch.NF <- subset(soil_r_values_arch, rownames(soil_r_values_arch)%in%daOTU_nofertil_arch.all) 
length(rownames(soil_r_values_arch.NF)) #54
soil_r_values_arch.NP <- subset(soil_r_values_arch, rownames(soil_r_values_arch)%in%daOTU_nopesti_arch.all) 
length(rownames(soil_r_values_arch.NP)) #2

soil_r_values_arch.rOTU <- rbind(soil_r_values_arch.CF,soil_r_values_arch.NF,soil_r_values_arch.NP)
write.xlsx(soil_r_values_arch.rOTU, "Arch_rOTU.xlsx")


soil_r_values_its.CF <- subset(soil_r_values_its, rownames(soil_r_values_its)%in%daOTU_conven_fun.all) 
length(rownames(soil_r_values_its.CF)) #199
soil_r_values_its.NF <- subset(soil_r_values_its, rownames(soil_r_values_its)%in%daOTU_nofertil_fun.all) 
length(rownames(soil_r_values_its.NF)) #430
soil_r_values_its.NP <- subset(soil_r_values_its, rownames(soil_r_values_its)%in%daOTU_nopesti_fun.all) 
length(rownames(soil_r_values_its.NP)) #205

soil_r_values_its.rOTU <- rbind(soil_r_values_its.CF,soil_r_values_its.NF,soil_r_values_its.NP)
write.xlsx(soil_r_values_its.rOTU, "Fun_rOTU.xlsx")

soil_r_values_its.rOTU.tax<-merge(soil_r_values_its.rOTU, df.tax_soil_its, by ="row.names") 
write.xlsx(soil_r_values_its.rOTU.tax, "Fun_rOTU.xlsx")
soil_r_values_16s.rOTU.tax<-merge(soil_r_values_16s.rOTU, df.tax_soil_16s, by ="row.names") 
write.xlsx(soil_r_values_16s.rOTU.tax, "Bac_rOTU.xlsx")
soil_r_values_arch.rOTU.tax<-merge(soil_r_values_arch.rOTU, df.tax_soil_arch, by ="row.names") 
write.xlsx(soil_r_values_arch.rOTU.tax, "Arch_rOTU.xlsx")



