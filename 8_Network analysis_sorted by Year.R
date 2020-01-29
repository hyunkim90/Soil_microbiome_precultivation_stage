library(igraph)
# The R notebooks ("CorrDF.R", "star_shape.R", and "triangle_shape.R") for constructing co-occurrence networks of soil microbial communities are available 
# at Additional file 4 of Hartman et al. 2018. Microbiome (https://doi.org/10.1186/s40168-017-0389-9).
source("CorrDF.R")
source("star_shape.R")
source("triangle_shape.R")



##Network analysis sorted by year
bac.clean.nolog
arch.clean.nolog
fun.clean.nolog

#All kingdom
all.clean.nolog <- merge_phyloseq(bac.clean.nolog, arch.clean.nolog)
all.clean.nolog <- merge_phyloseq(all.clean.nolog, fun.clean.nolog)

all.clean.nolog.17 <- subset_samples(all.clean.nolog, Year == "2017")
all.clean.nolog.17 <- filter_taxa(all.clean.nolog.17 , function(x) sum(x) != 0, TRUE)
all.clean.nolog.17.5site <- subset_samples(all.clean.nolog.17, Location %in% c("CC", "CJ", "DG", "IS", "UF"))
all.clean.nolog.17.5site <- filter_taxa(all.clean.nolog.17.5site, function(x) sum(x) != 0, TRUE)

all.clean.nolog.18 <- subset_samples(all.clean.nolog, Year == "2018")
all.clean.nolog.18 <- filter_taxa(all.clean.nolog.18 , function(x) sum(x) != 0, TRUE)


#Bac
bac.clean.nolog.17 <- subset_samples(bac.clean.nolog, Year == "2017")
bac.clean.nolog.17 <- filter_taxa(bac.clean.nolog.17 , function(x) sum(x) != 0, TRUE)
bac.clean.nolog.17.5site <- subset_samples(bac.clean.nolog.17, Location %in% c("CC", "CJ", "DG", "IS", "UF"))
bac.clean.nolog.17.5site <- filter_taxa(bac.clean.nolog.17.5site, function(x) sum(x) != 0, TRUE)

bac.clean.nolog.18 <- subset_samples(bac.clean.nolog, Year == "2018")
bac.clean.nolog.18 <- filter_taxa(bac.clean.nolog.18 , function(x) sum(x) != 0, TRUE)

#Fun
fun.clean.nolog.17 <- subset_samples(fun.clean.nolog, Year == "2017")
fun.clean.nolog.17 <- filter_taxa(fun.clean.nolog.17 , function(x) sum(x) != 0, TRUE)
fun.clean.nolog.17.5site <- subset_samples(fun.clean.nolog.17, Location %in% c("CC", "CJ", "DG", "IS", "UF"))
fun.clean.nolog.17.5site <- filter_taxa(fun.clean.nolog.17.5site, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.18 <- subset_samples(fun.clean.nolog, Year == "2018")
fun.clean.nolog.18 <- filter_taxa(fun.clean.nolog.18 , function(x) sum(x) != 0, TRUE)

#Archaea
arch.clean.nolog.17 <- subset_samples(arch.clean.nolog, Year == "2017")
arch.clean.nolog.17 <- filter_taxa(arch.clean.nolog.17 , function(x) sum(x) != 0, TRUE)
arch.clean.nolog.17.5site <- subset_samples(arch.clean.nolog.17, Location %in% c("CC", "CJ", "DG", "IS", "UF"))
arch.clean.nolog.17.5site <- filter_taxa(arch.clean.nolog.17.5site, function(x) sum(x) != 0, TRUE)

arch.clean.nolog.18 <- subset_samples(arch.clean.nolog, Year == "2018")
arch.clean.nolog.18 <- filter_taxa(arch.clean.nolog.18 , function(x) sum(x) != 0, TRUE)



### Multi-Kingdom metanetwork
#Phyloseq sorted by cultural practice
#Bacteria
otu_norm_soil_16s.17 <- otu_table(bac.clean.nolog.17.5site)
rownames(otu_norm_soil_16s.17)
otu_norm_soil_16s.17 <- data.frame(otu_norm_soil_16s.17)
otu_norm_soil_16s.17$OTU <-rownames(otu_norm_soil_16s.17)
head(otu_norm_soil_16s.17)
otu_norm_soil_16s.17 <- merge(otu_norm_soil_16s.17, OTU_id.list, by = "OTU")
rownames(otu_norm_soil_16s.17) <- otu_norm_soil_16s.17$OTU_id
length(colnames(otu_norm_soil_16s.17))
otu_norm_soil_16s.17 <- otu_norm_soil_16s.17[-c(1, 92)]
head(otu_norm_soil_16s.17)


otu_norm_soil_16s.18 <- otu_table(bac.clean.nolog.18)
rownames(otu_norm_soil_16s.18)
otu_norm_soil_16s.18 <- data.frame(otu_norm_soil_16s.18)
otu_norm_soil_16s.18$OTU <-rownames(otu_norm_soil_16s.18)
head(otu_norm_soil_16s.18)
otu_norm_soil_16s.18 <- merge(otu_norm_soil_16s.18, OTU_id.list, by = "OTU")
rownames(otu_norm_soil_16s.18) <- otu_norm_soil_16s.18$OTU_id
length(colnames(otu_norm_soil_16s.18))
otu_norm_soil_16s.18 <- otu_norm_soil_16s.18[-c(1, 92)]
head(otu_norm_soil_16s.18)

#Fungi
otu_norm_soil_its.17 <- otu_table(fun.clean.nolog.17.5site)
rownames(otu_norm_soil_its.17)
otu_norm_soil_its.17 <- data.frame(otu_norm_soil_its.17)
otu_norm_soil_its.17$OTU <-rownames(otu_norm_soil_its.17)
head(otu_norm_soil_its.17)
otu_norm_soil_its.17 <- merge(otu_norm_soil_its.17, OTU_id.list, by = "OTU")
rownames(otu_norm_soil_its.17) <- otu_norm_soil_its.17$OTU_id
length(colnames(otu_norm_soil_its.17))
otu_norm_soil_its.17 <- otu_norm_soil_its.17[-c(1, 92)]
head(otu_norm_soil_its.17)

otu_norm_soil_its.18 <- otu_table(fun.clean.nolog.18)
rownames(otu_norm_soil_its.18)
otu_norm_soil_its.18 <- data.frame(otu_norm_soil_its.18)
otu_norm_soil_its.18$OTU <-rownames(otu_norm_soil_its.18)
head(otu_norm_soil_its.18)
otu_norm_soil_its.18 <- merge(otu_norm_soil_its.18, OTU_id.list, by = "OTU")
rownames(otu_norm_soil_its.18) <- otu_norm_soil_its.18$OTU_id
length(colnames(otu_norm_soil_its.18))
otu_norm_soil_its.18 <- otu_norm_soil_its.18[-c(1, 92)]


##Archaea
otu_norm_soil_arch.17 <- otu_table(arch.clean.nolog.17.5site)
rownames(otu_norm_soil_arch.17)
otu_norm_soil_arch.17 <- data.frame(otu_norm_soil_arch.17)
otu_norm_soil_arch.17$OTU <-rownames(otu_norm_soil_arch.17)
head(otu_norm_soil_arch.17)
otu_norm_soil_arch.17 <- merge(otu_norm_soil_arch.17, OTU_id.list, by = "OTU")
rownames(otu_norm_soil_arch.17) <- otu_norm_soil_arch.17$OTU_id
length(colnames(otu_norm_soil_arch.17))
otu_norm_soil_arch.17 <- otu_norm_soil_arch.17[-c(1, 92)]

otu_norm_soil_arch.18 <- otu_table(arch.clean.nolog.18)
rownames(otu_norm_soil_arch.18)
otu_norm_soil_arch.18 <- data.frame(otu_norm_soil_arch.18)
otu_norm_soil_arch.18$OTU <-rownames(otu_norm_soil_arch.18)
head(otu_norm_soil_arch.18)
otu_norm_soil_arch.18 <- merge(otu_norm_soil_arch.18, OTU_id.list, by = "OTU")
rownames(otu_norm_soil_arch.18) <- otu_norm_soil_arch.18$OTU_id
length(colnames(otu_norm_soil_arch.18))
otu_norm_soil_arch.18 <- otu_norm_soil_arch.18[-c(1, 92)]


##Extract OTU tables (beacuse of sample name) in 2018 samples
otu.all.clean.nolog.18 <- otu_table(all.clean.nolog.18)
write.table(otu.all.clean.nolog.18, "otu.all.clean.nolog.18.txt", sep = '\t', quote = F)

##Import edited OTU_table and make it phyloseq class
otu.all.clean.nolog.18<-read.table("otu.all.clean.nolog.18.txt", sep='\t', header =T)
rownames(otu.all.clean.nolog.18) <- otu.all.clean.nolog.18$X
otu.all.clean.nolog.18 <- otu.all.clean.nolog.18[-c(1)]
head(otu.all.clean.nolog.18)
length(colnames(otu.all.clean.nolog.18))

##### Soil meta co-occurrence network creation and analysis and defining keystone OTUs #####

## Combine OTU counts of both kingdoms together

## Perform Spearman correlation of all OTU pairs
otu_norm_soil_combine.17 <- otu_table(all.clean.nolog.17.5site)
write.table(otu_norm_soil_combine.17, "otu_norm_soil_combine.17.txt", sep='\t', quote = F)
length(rownames(otu_norm_soil_combine.17))
length(rownames(otu_norm_soil_16s.17))
length(rownames(otu_norm_soil_arch.17))
length(rownames(otu_norm_soil_its.17))


otu_norm_soil_combine.18 <- otu.all.clean.nolog.18
write.table(otu_norm_soil_combine.18, "otu_norm_soil_combine.18.txt", sep='\t', quote = F)
length(rownames(otu_norm_soil_combine.18))
length(rownames(otu_norm_soil_16s.18))
length(rownames(otu_norm_soil_arch.18))
length(rownames(otu_norm_soil_its.18))

##Estimating correlation
otu_norm_soil_combine.17 <-read.table("otu_norm_soil_combine.17.txt", sep='\t', header = T)

otu_norm_soil_combine.17$OTU <-rownames(otu_norm_soil_combine.17)
otu_norm_soil_combine.17 <- merge(otu_norm_soil_combine.17, OTU_id.list, by = "OTU")
rownames(otu_norm_soil_combine.17) <- otu_norm_soil_combine.17$OTU_id
otu_norm_soil_combine.17 <- otu_norm_soil_combine.17[-c(1, 92)]


all_soil_cor.17 <- rcorr(t(otu_norm_soil_combine.17), type=c("spearman"))


otu_norm_soil_combine.18 <-read.table("otu_norm_soil_combine.18.txt", sep='\t', header = T)

otu_norm_soil_combine.18$OTU <-rownames(otu_norm_soil_combine.18)
otu_norm_soil_combine.18 <- merge(otu_norm_soil_combine.18, OTU_id.list, by = "OTU")
rownames(otu_norm_soil_combine.18) <- otu_norm_soil_combine.18$OTU_id
otu_norm_soil_combine.18 <- otu_norm_soil_combine.18[-c(1, 92)]


all_soil_cor.18 <- rcorr(t(otu_norm_soil_combine.18), type=c("spearman"))



## Create data frame of co-occurring OTUs
all_cor_soil_df.17 <- CorrDF(all_soil_cor.17$r, all_soil_cor.17$P)
all_cor_soil_df.17$padj <- p.adjust(all_cor_soil_df.17$p, method="none")
all_cor_soil_df_padj.17 <- all_cor_soil_df.17[which(abs(all_cor_soil_df.17$cor) > 0.7),]
all_cor_soil_df_padj.17 <- all_cor_soil_df_padj.17[which(all_cor_soil_df_padj.17$padj < 0.001),]

write.table(all_cor_soil_df_padj.17, "all_cor_soil_df_padj.17.txt", sep='\t', quote = F)


all_cor_soil_df.18 <- CorrDF(all_soil_cor.18$r, all_soil_cor.18$P)
all_cor_soil_df.18$padj <- p.adjust(all_cor_soil_df.18$p, method="none")
all_cor_soil_df_padj.18 <- all_cor_soil_df.18[which(abs(all_cor_soil_df.18$cor) > 0.7),]
all_cor_soil_df_padj.18 <- all_cor_soil_df_padj.18[which(all_cor_soil_df_padj.18$padj < 0.001),]

write.table(all_cor_soil_df_padj.18, "all_cor_soil_df_padj.18.txt", sep='\t', quote = F)

## Network construction
all_cor_soil_df_padj.17 <- read.table("all_cor_soil_df_padj.17.txt", sep='\t', header = T)

#rOTU
soil_r_values_combine.all.rOTU <- rbind(soil_r_values_16s.rOTU,soil_r_values_arch.rOTU,soil_r_values_its.rOTU)
df.soil_r_values_combine.all.rOTU <- data.frame(soil_r_values_combine.all.rOTU)
df.soil_r_values_combine.all.rOTU <- subset(df.soil_r_values_combine.all.rOTU, index == 1|index == 2|index == 3)


soil_r_values_combine.all.rOTU <- df.soil_r_values_combine.all.rOTU
soil_r_values_combine.all.rOTU[,1:3]


indic_edge_soil_combine.all <- rownames(soil_r_values_combine.all.rOTU)

#indicator OTU list
soil_r_values_16s
soil_r_values_arch
soil_r_values_its

soil_r_values_combine.all <- rbind(soil_r_values_16s,soil_r_values_arch,soil_r_values_its)
df.soil_r_values_combine.all <- data.frame(soil_r_values_combine.all)
df.soil_r_values_combine.all <- subset(df.soil_r_values_combine.all, index == 1|index == 2|index == 3)


soil_r_values_combine.all <- df.soil_r_values_combine.all
soil_r_values_combine.all[,1:3]
rownames(soil_r_values_combine.all)

nodeattrib_soil_combine.17 <- data.frame(node=union(all_cor_soil_df_padj.17$from,all_cor_soil_df_padj.17$to))
nodeattrib_soil_combine.17$indicgroup <- 0

for (i in as.character(nodeattrib_soil_combine.17$node))
{
  if (i %in% indic_edge_soil_combine.all == TRUE)
  {nodeattrib_soil_combine.17[nodeattrib_soil_combine.17$node==i,"indicgroup"] <- paste(colnames(soil_r_values_combine.all)[which(soil_r_values_combine.all[i,1:3]==1)],collapse = "_")}
  else
  {nodeattrib_soil_combine.17[nodeattrib_soil_combine.17$node==i,"indicgroup"]<- "NA"}
}

rownames(nodeattrib_soil_combine.17) <- as.character(nodeattrib_soil_combine.17$node)
nodeattrib_soil_combine.17$indicgroup

library(igraph)

all_soil_net.17 <- graph_from_data_frame(all_cor_soil_df_padj.17,direct=F,vertices=nodeattrib_soil_combine.17)


## Number of nodes
length(V(all_soil_net.17)) #15645

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.17)))) #11679
length(grep("^F",names(V(all_soil_net.17)))) #3312
length(grep("^A",names(V(all_soil_net.17)))) #654

## Connections 
bb_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.17) #555775

ff_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.17) #41029

bf_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^B",from) & grepl("^F",to)),])
nrow(bf_occur_soil.17) #160756

fb_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.17) #111041

aa_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.17) #1990

ab_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^A",from) & grepl("^B",to)),])
nrow(ab_occur_soil.17) #27454

af_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^A",from) & grepl("^F",to)),])
nrow(af_occur_soil.17) #8030

ba_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.17) #33999

fa_occur_soil.17 <- droplevels(all_cor_soil_df_padj.17[with(all_cor_soil_df_padj.17, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.17) #7413



##Positive and negative connections
bb_occur_soil.17.neg<-subset(bb_occur_soil.17, bb_occur_soil.17$cor < 0)
nrow(bb_occur_soil.17.neg) #2

ff_occur_soil.17.neg<-subset(ff_occur_soil.17, ff_occur_soil.17$cor < 0)
nrow(ff_occur_soil.17.neg) #5

bf_occur_soil.17.neg<-subset(bf_occur_soil.17, bf_occur_soil.17$cor < 0)
nrow(bf_occur_soil.17.neg) #1

aa_occur_soil.17.neg<-subset(aa_occur_soil.17, aa_occur_soil.17$cor < 0)
nrow(aa_occur_soil.17.neg) #0

ab_occur_soil.17.neg<-subset(ab_occur_soil.17, ab_occur_soil.17$cor < 0)
nrow(ab_occur_soil.17.neg) #0

af_occur_soil.17.neg<-subset(af_occur_soil.17, af_occur_soil.17$cor < 0)
nrow(af_occur_soil.17.neg) #1

ba_occur_soil.17.neg<-subset(ba_occur_soil.17, ba_occur_soil.17$cor < 0)
nrow(ba_occur_soil.17.neg) #1

fa_occur_soil.17.neg<-subset(fa_occur_soil.17, fa_occur_soil.17$cor < 0)
nrow(fa_occur_soil.17.neg) #1

fb_occur_soil.17.neg<-subset(fb_occur_soil.17, fb_occur_soil.17$cor < 0)
nrow(fb_occur_soil.17.neg) #1

neg.occur.17<-rbind(bb_occur_soil.17.neg, ff_occur_soil.17.neg, bf_occur_soil.17.neg,fa_occur_soil.17.neg,fb_occur_soil.17.neg)
write.xlsx(neg.occur.17, "Negative correlation in 17 network_with archaea.xlsx")

## Node degree values
soil_all_deg.17 <- sort(degree(all_soil_net.17,mode="all"),decr=T)
max(soil_all_deg.17) #241
mean(soil_all_deg.17) #121.1233

write.table(soil_all_deg.17,"soil_all_deg_conventional.txt",sep="\t",quote=F)

df.soil_all_deg.17 <- data.frame(soil_all_deg.17)
colnames(df.soil_all_deg.17) <- c("Degree")


## Degree
n <- 1
soil_all_degree.17 <- soil_all_deg.17[soil_all_deg.17 >= quantile(soil_all_deg.17,prob=1-n/100)]
length(soil_all_degree.17) #322
degree.17.1pro<-data.frame(soil_all_degree.17)
colnames(degree.17.1pro) <- c("Degree")
rownames(degree.17.1pro)

##Defining Keystone OTU by betweenness
soil_all_between.17<-betweenness(all_soil_net.17, normalized = TRUE)

between.17<-data.frame(soil_all_between.17)
colnames(between.17) <- c("Betweenness")
head(between.17)

n <- 1 #this threshold was used for predicting hub OTUs
soil_all_keystone.between.17 <- soil_all_between.17[soil_all_between.17 >= quantile(soil_all_between.17,prob=1-n/100)]
length(soil_all_keystone.between.17)#157

between.2017.1pro.cc<-data.frame(soil_all_keystone.between.17)
colnames(between.2017.1pro.cc) <- c("Betweenness")

##Eigenvector
soil_all_eigen.CF<-eigen_centrality(all_soil_net.CF, directed = FALSE, scale = TRUE,
                                    weights = NULL, options = arpack_defaults)
df.soil_all_eigen.CF<-data.frame(soil_all_eigen.CF)
df.soil_all_eigen.CF$OTU <- rownames(df.soil_all_eigen.CF)
df.soil_all_eigen.CF <- df.soil_all_eigen.CF[,-c(2:22)]
write.table(df.soil_all_eigen.CF,"df.soil_all_eigen.CF.txt",sep="\t",quote=F)

n <- 1
soil_all_keystone.eigen.CF <- df.soil_all_eigen.CF$OTU[df.soil_all_eigen.CF$vector> quantile(df.soil_all_eigen.CF$vector,prob=1-n/100)]
length(soil_all_keystone.eigen.CF) #112

###keystone defined by degree,betweenness, and closeness
## Set node shape soil_all_keystone.2017
V(all_soil_net.17)$shape <- V(all_soil_net.17)$name
V(all_soil_net.17)$shape[V(all_soil_net.17)$shape %in% soil_all_keystone.between.17] <- "star"
V(all_soil_net.17)$shape[V(all_soil_net.17)$shape %in% rownames(otu_norm_soil_16s.17)] <- "circle"
V(all_soil_net.17)$shape[V(all_soil_net.17)$shape %in% rownames(otu_norm_soil_its.17)] <- "triangle"
V(all_soil_net.17)$shape[V(all_soil_net.17)$shape %in% rownames(otu_norm_soil_arch.17)] <- "square"

## Set node colors based upon sensitivity to management system
cs <- c("Conventional", "No_fertilizer", "No_pesticide")
unique(V(all_soil_net.17)$indicgroup)
V(all_soil_net.17)$color <- V(all_soil_net.17)$indicgroup
V(all_soil_net.17)$color[V(all_soil_net.17)$color == "NA"] <- "gray50"
V(all_soil_net.17)$color[V(all_soil_net.17)$color == "No_fertilizer"] <- "#0066CC"
V(all_soil_net.17)$color[V(all_soil_net.17)$color == "No_pesticide"] <- "#336633"
V(all_soil_net.17)$color[V(all_soil_net.17)$color == "Conventional"] <- "#CC9900"
V(all_soil_net.17)$frame.color <- V(all_soil_net.17)$color

#node attributed to important rOTU
soil_all_net_csnodes.17 <- rownames(nodeattrib_soil_combine.17[nodeattrib_soil_combine.17$indicgroup %in% cs,]) #559
soil_all_net_csnodes_bac.17 <- soil_all_net_csnodes.17[grep("^B",soil_all_net_csnodes.17)] #259
soil_all_net_csnodes_fun.17 <- soil_all_net_csnodes.17[grep("^F",soil_all_net_csnodes.17)] #281
soil_all_net_csnodes_arch.17 <- soil_all_net_csnodes.17[grep("^A",soil_all_net_csnodes.17)] #19

#Node size
V(all_soil_net.17)$size <- V(all_soil_net.17)$name
V(all_soil_net.17)$size[V(all_soil_net.17)$size %in% soil_all_keystone.between.17] <- 3
V(all_soil_net.17)$size[V(all_soil_net.17)$size %in% soil_all_net_csnodes_bac.17] <- 4
V(all_soil_net.17)$size[V(all_soil_net.17)$size %in% soil_all_net_csnodes_arch.17] <- 4
V(all_soil_net.17)$size[V(all_soil_net.17)$size %in% soil_all_net_csnodes_fun.17] <- 5
V(all_soil_net.17)$size[V(all_soil_net.17)$size %in% rownames(otu_norm_soil_16s.17)] <- 2
V(all_soil_net.17)$size[V(all_soil_net.17)$size %in% rownames(otu_norm_soil_its.17)] <- 3
V(all_soil_net.17)$size[V(all_soil_net.17)$size %in% rownames(otu_norm_soil_arch.17)] <- 2
soilcombine_nodesizes.17 <- as.numeric(V(all_soil_net.17)$size)

## Set edge color
E(all_soil_net.17)$color <- ifelse(E(all_soil_net.17)$cor >0, "#99CCFF","#FF6666")

##### defining modules #####

## Make vectors of network nodes responding to different cropping systems
CF_nodes_soil.17 <- rownames(nodeattrib_soil_combine.17[nodeattrib_soil_combine.17$indicgroup=="Conventional",])
NP_nodes_soil.17 <- rownames(nodeattrib_soil_combine.17[nodeattrib_soil_combine.17$indicgroup=="No_pesticide",])
NF_nodes_soil.17 <- rownames(nodeattrib_soil_combine.17[nodeattrib_soil_combine.17$indicgroup=="No_fertilizer",])

cs_nodes_soil_all.17 <- c(CF_nodes_soil.17,NF_nodes_soil.17, NP_nodes_soil.17)
Bcs_nodes_soil_all.17 <- cs_nodes_soil_all.17[grep("^B",cs_nodes_soil_all.17)] #234
Fcs_nodes_soil_all.17 <- cs_nodes_soil_all.17[grep("^F",cs_nodes_soil_all.17)] #199
Acs_nodes_soil_all.17 <- cs_nodes_soil_all.17[grep("^A",cs_nodes_soil_all.17)] #18

## Perform cluster analysis using greedy clustering algorithm 
cfg_soil.17 <- cluster_fast_greedy(as.undirected(all_soil_net.17))

## Subset for top 20 biggest nodes
soil_modules.17 <- sort(table(membership(cfg_soil.17)),decr=T)
soil_modules_20.17 <- soil_modules.17[1:20]

sum(soil_modules_20.17)/sum(soil_modules.17)
sm20_plot.17 <- soil_modules_20.17
names(sm20_plot.17) <- as.factor(1:20)
soil_modules_cs.17 <- table(factor(membership(cfg_soil.17)[cs_nodes_soil_all.17],levels=names(soil_modules.17)))
soil_modules_cs_20.17 <- soil_modules_cs.17[names(soil_modules_20.17)]
smcs20_plot.17 <- soil_modules_cs_20.17
names(smcs20_plot.17) <- as.factor(1:20)

## Make vector of nodes in top 20 modules
soil_modules_points.17 <- membership(cfg_soil.17)[membership(cfg_soil.17) %in% names(soil_modules_20.17)]

soil_points.17 <- NULL
for(i in soil_modules_points.17){
  soilx <- which(names(soil_modules_20.17)==i)
  soil_points.17 <- c(soil_points.17, soilx)
}
names(soil_points.17) <- names(soil_modules_points.17)

## Set node colors by cropping system sensitivity 
soil_all_cols.17 <- sort(soil_points.17)
soil_all_cols.17[!names(soil_all_cols.17) %in% cs] <- "gray50"
soil_all_cols.17[names(soil_all_cols.17) %in% CF_nodes_soil.17] <- "#CC9900"
soil_all_cols.17[names(soil_all_cols.17) %in% NF_nodes_soil.17] <- "#0066CC"
soil_all_cols.17[names(soil_all_cols.17) %in% NP_nodes_soil.17] <- "#336633"

soil_all_pch.17 <- sort(soil_points.17)
soil_all_pch.17[names(soil_all_pch.17) %in% rownames(otu_norm_soil_16s.17)] <- 1
soil_all_pch.17[names(soil_all_pch.17) %in% rownames(otu_norm_soil_its.17)] <- 2
soil_all_pch.17[names(soil_all_pch.17) %in% rownames(otu_norm_soil_arch.17)] <- 3
soil_all_pch.17[names(soil_all_pch.17) %in% intersect(rownames(otu_norm_soil_16s.17),cs_nodes_soil_all.17)] <- 16
soil_all_pch.17[names(soil_all_pch.17) %in% intersect(rownames(otu_norm_soil_its.17),cs_nodes_soil_all.17)] <- 16
soil_all_pch.17[names(soil_all_pch.17) %in% intersect(rownames(otu_norm_soil_arch.17),cs_nodes_soil_all.17)] <- 16
soil_all_pch.17[names(soil_all_pch.17) %in% soil_all_keystone.between.17] <- 8

soil_all_cex.17 <- sort(soil_points.17)
soil_all_cex.17[!names(soil_all_cex.17) %in% cs_nodes_soil_all.17] <- 1
soil_all_cex.17[names(soil_all_cex.17) %in% cs_nodes_soil_all.17] <- 2

soil_mods_list_cs.17 <- list()
for (i in names(soil_modules_cs_20.17)){
  x1 <- names(membership(cfg_soil.17)[membership(cfg_soil.17)==i])
  x2 <- x1[x1 %in% cs_nodes_soil_all.17]
  soil_mods_list_cs.17[[i]] <- as.numeric(V(all_soil_net.17)[x2])
}


##### meta co-occurrence networks #####
## Note: the permuations for the layout of the network can be very time consuming and processor intensive
# set.seed(8051)
coords_soil_all.CF <- layout_(all_soil_net.CF, with_fr(niter=9999, grid="nogrid"))
write.table(coords_soil_all.CF,"coords_soil_all.CF.txt",sep="\t",row.names=F,col.names=F,quote=F)


## Import pre-calculated FR layout coordinates to save time 
coords_soil_all.CF <- as.matrix(read.table("coords_soil_all.CF.txt"))
dimnames(coords_soil_all.CF) <-  NULL

par(mfrow=c(1,1), mar=c(0,0,0,0))


soil_cols.17 <- c("antiquewhite1", "antiquewhite2","antiquewhite3","wheat3")

plot(all_soil_net.CF,vertex.label=NA,vertex.size=soilcombine_nodesizes.CF, layout=coords_soil_all.CF,
     mark.groups=list(soil_mods_list_cs.CF$`1`, soil_mods_list_cs.CF$`8`, soil_mods_list_cs.CF$`12`),
     mark.col=soil_cols.CF, mark.border=soil_cols.CF)

plot(all_soil_net.2017.cc,vertex.label=NA,vertex.size=soilcombine_nodesizes.2017.cc, layout=coords_soil_all.2017.cc)


legend("bottomleft",legend=c("Module 1", "Module 5", "Module 10", "Module 11", "Module 13", "Module 22"),col=soil_cols.2017,
       bty="n",fill=soil_cols.2017,border=soil_cols.2017)

pdf(p_combine_net,"cor0.7_combine_network_2017_nolog.pdf",height=20,width=20)
pdf(p_combine_legend,"cor0.7_combine_network_2017_nolog_combine_legend.pdf",height=10,width=10)


library(plotrix)
##Distribution of rOTUs in modules
plot(sequence(sm20_plot.17)~jitter(sort(soil_points.17),1.75),pch=soil_all_pch.17, col=soil_all_cols.17, cex=soil_all_cex.17,
     ann=F,axes=F,ylim=c(0,800),xlim=c(0,20))
axis(1,at=1:20, labels=F, tcl=-0.15)
staxlab(side=1,at=1:20,labels=paste(rep("M",20),names(soil_modules_20.17),sep=""),srt=45)
axis(1,at=seq(1,20,1), tick=F, labels=paste(round(smcs20_plot.17/sm20_plot.17*100,1),"%",sep=""),line=1.5,cex.axis=0.7)
axis(2,at=seq(0,700,100))
text(-0.25,-50,"rOTUs\n in module:",xpd=NA)
title(ylab="Number of nodes per module",line=2)


dev.off()


##modified FR layout
# To do this, node names have to be changed into number
library(qgraph)
all_soil_net.17.integer <- set.vertex.attribute(all_soil_net.17, "name", value=paste(1:15645))
edge.17 <- get.edgelist(all_soil_net.17.integer)
edge.17.num<- apply(edge.17, 2, as.numeric)

layout.17.modified <- qgraph.layout.fruchtermanreingold(edge.17.num,vcount=vcount(all_soil_net.17.integer),
                                                        area=10*(vcount(all_soil_net.17.integer)^2),repulse.rad=(vcount(all_soil_net.17.integer)^3.1))

#without module
plot(all_soil_net.17.integer,layout=layout.17.modified,vertex.size=soilcombine_nodesizes.17,vertex.label=NA)
V(all_soil_net.17.integer)$color

#with module

plot(all_soil_net.17.integer,layout=layout.17.modified,vertex.size=soilcombine_nodesizes.17,vertex.label=NA,
     mark.groups=list(soil_mods_list_cs.17$`1`, soil_mods_list_cs.17$`3`, soil_mods_list_cs.17$`4`, soil_mods_list_cs.17$`6`),
     mark.col=soil_cols.17, mark.border=soil_cols.17)

legend("bottomleft",legend=c("Module 1", "Module 3", "Module 4", "Module 6"),col=soil_cols.17,
       bty="n",fill=soil_cols.17,border=soil_cols.17)


dev.off()

mtext("qgraph.layout.fruchtermanreingold modified", side=1)

dev.off()

## 2018 metanetwork

all_cor_soil_df_padj.18 <- read.table("all_cor_soil_df_padj.18.txt", sep='\t', header = T)

#rOTU
soil_r_values_combine.all.rOTU <- rbind(soil_r_values_16s.rOTU,soil_r_values_arch.rOTU,soil_r_values_its.rOTU)
df.soil_r_values_combine.all.rOTU <- data.frame(soil_r_values_combine.all.rOTU)
df.soil_r_values_combine.all.rOTU <- subset(df.soil_r_values_combine.all.rOTU, index == 1|index == 2|index == 3)


soil_r_values_combine.all.rOTU <- df.soil_r_values_combine.all.rOTU
soil_r_values_combine.all.rOTU[,1:3]


indic_edge_soil_combine.all <- rownames(soil_r_values_combine.all.rOTU)

#indicator OTU list
soil_r_values_16s
soil_r_values_arch
soil_r_values_its

soil_r_values_combine.all <- rbind(soil_r_values_16s,soil_r_values_arch,soil_r_values_its)
df.soil_r_values_combine.all <- data.frame(soil_r_values_combine.all)
df.soil_r_values_combine.all <- subset(df.soil_r_values_combine.all, index == 1|index == 2|index == 3)


soil_r_values_combine.all <- df.soil_r_values_combine.all
soil_r_values_combine.all[,1:3]
rownames(soil_r_values_combine.all)

nodeattrib_soil_combine.18 <- data.frame(node=union(all_cor_soil_df_padj.18$from,all_cor_soil_df_padj.18$to))
nodeattrib_soil_combine.18$indicgroup <- 0

for (i in as.character(nodeattrib_soil_combine.18$node))
{
  if (i %in% indic_edge_soil_combine.all == TRUE)
  {nodeattrib_soil_combine.18[nodeattrib_soil_combine.18$node==i,"indicgroup"] <- paste(colnames(soil_r_values_combine.all)[which(soil_r_values_combine.all[i,1:3]==1)],collapse = "_")}
  else
  {nodeattrib_soil_combine.18[nodeattrib_soil_combine.18$node==i,"indicgroup"]<- "NA"}
}

rownames(nodeattrib_soil_combine.18) <- as.character(nodeattrib_soil_combine.18$node)
nodeattrib_soil_combine.18$indicgroup

library(igraph)

all_soil_net.18 <- graph_from_data_frame(all_cor_soil_df_padj.18,direct=F,vertices=nodeattrib_soil_combine.18)


## Number of nodes
length(V(all_soil_net.18)) #10765

## Number of bacteria and fungi nodes
length(grep("^B",names(V(all_soil_net.18)))) #7250
length(grep("^F",names(V(all_soil_net.18)))) #3228
length(grep("^A",names(V(all_soil_net.18)))) #287

## Connections 
bb_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^B",from) & grepl("^B",to)),])
nrow(bb_occur_soil.18) #212287

ff_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^F",from) & grepl("^F",to)),])
nrow(ff_occur_soil.18) #38987

bf_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^B",from) & grepl("^F",to)),])
nrow(bf_occur_soil.18) #98697

fb_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^F",from) & grepl("^B",to)),])
nrow(fb_occur_soil.18) #68294

aa_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^A",from) & grepl("^A",to)),])
nrow(aa_occur_soil.18) #321

ab_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^A",from) & grepl("^B",to)),])
nrow(ab_occur_soil.18) #7546

af_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^A",from) & grepl("^F",to)),])
nrow(af_occur_soil.18) #3883

ba_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^B",from) & grepl("^A",to)),])
nrow(ba_occur_soil.18) #8296

fa_occur_soil.18 <- droplevels(all_cor_soil_df_padj.18[with(all_cor_soil_df_padj.18, grepl("^F",from) & grepl("^A",to)),])
nrow(fa_occur_soil.18) #2924

##Positive and negative connections
bb_occur_soil.18.neg<-subset(bb_occur_soil.18, bb_occur_soil.18$cor < 0)
nrow(bb_occur_soil.18.neg) #9

ff_occur_soil.18.neg<-subset(ff_occur_soil.18, ff_occur_soil.18$cor < 0)
nrow(ff_occur_soil.18.neg) #4

bf_occur_soil.18.neg<-subset(bf_occur_soil.18, bf_occur_soil.18$cor < 0)
nrow(bf_occur_soil.18.neg) #7

aa_occur_soil.18.neg<-subset(aa_occur_soil.18, aa_occur_soil.18$cor < 0)
nrow(aa_occur_soil.18.neg) #0

ab_occur_soil.18.neg<-subset(ab_occur_soil.18, ab_occur_soil.18$cor < 0)
nrow(ab_occur_soil.18.neg) #5

af_occur_soil.18.neg<-subset(af_occur_soil.18, af_occur_soil.18$cor < 0)
nrow(af_occur_soil.18.neg) #1

ba_occur_soil.18.neg<-subset(ba_occur_soil.18, ba_occur_soil.18$cor < 0)
nrow(ba_occur_soil.18.neg) #3

fa_occur_soil.18.neg<-subset(fa_occur_soil.18, fa_occur_soil.18$cor < 0)
nrow(fa_occur_soil.18.neg) #4

fb_occur_soil.18.neg<-subset(fb_occur_soil.18, fb_occur_soil.18$cor < 0)
nrow(fb_occur_soil.18.neg) #12

neg.occur.18<-rbind(bb_occur_soil.18.neg, ff_occur_soil.18.neg, bf_occur_soil.18.neg,fa_occur_soil.18.neg,fb_occur_soil.18.neg)
write.xlsx(neg.occur.18, "Negative correlation in 18 network_with archaea.xlsx")

## Node degree values
soil_all_deg.18 <- sort(degree(all_soil_net.18,mode="all"),decr=T)
max(soil_all_deg.18) #241
mean(soil_all_deg.18) #121.1233

write.table(soil_all_deg.18,"soil_all_deg_conventional.txt",sep="\t",quote=F)

df.soil_all_deg.18 <- data.frame(soil_all_deg.18)
colnames(df.soil_all_deg.18) <- c("Degree")


## Degree
n <- 1
soil_all_degree.18 <- soil_all_deg.18[soil_all_deg.18 >= quantile(soil_all_deg.18,prob=1-n/100)]
length(soil_all_degree.18) #322
degree.18.1pro<-data.frame(soil_all_degree.18)
colnames(degree.18.1pro) <- c("Degree")
rownames(degree.18.1pro)

intersect(rownames(degree.18.1pro), list.glomeromycota)
subset(df.soil_all_deg.18,rownames(df.soil_all_deg.18)%in%list.glomeromycota)

intersect(rownames(degree.18.1pro), list.rhizobiaceae$OTU_id) #B4037_Mesorhizobium

## Get bacteria and fungi keystone OTUs
soil_net_keystone3_bac.18 <- names(soil_all_keystone3.18)[grep("^B",names(soil_all_keystone3.18))]
length(soil_net_keystone3_bac.18)

soil_net_keystone3_fun.18 <- names(soil_all_keystone3.18)[grep("^F",names(soil_all_keystone3.18))]
length(soil_net_keystone3_fun.18)


##Defining Keystone OTU by betweenness
soil_all_between.18<-betweenness(all_soil_net.18, normalized = TRUE)

between.18<-data.frame(soil_all_between.18)
colnames(between.18) <- c("Betweenness")
head(between.18)

n <- 1 #this threshold was used for predicting hub OTUs
soil_all_keystone.between.18 <- soil_all_between.18[soil_all_between.18 >= quantile(soil_all_between.18,prob=1-n/100)]
length(soil_all_keystone.between.18)#108

## Get bacteria and fungi keystone OTUs
soil_net_keystone_bac.between.18 <- names(soil_all_keystone.between.18.score)[grep("^B",names(soil_all_keystone.between.18.score))]
length(soil_net_keystone_bac.between.18)#72

soil_net_keystone_fun.between.18 <- names(soil_all_keystone.between.18.score)[grep("^F",names(soil_all_keystone.between.18.score))]
length(soil_net_keystone_fun.between.18)#35

soil_net_keystone_arch.between.18 <- names(soil_all_keystone.between.18.score)[grep("^A",names(soil_all_keystone.between.18.score))]
length(soil_net_keystone_arch.between.18)#1



##Eigenvector
soil_all_eigen.18<-eigen_centrality(all_soil_net.18, directed = FALSE, scale = TRUE,
                                    weights = NULL, options = arpack_defaults)
df.soil_all_eigen.18<-data.frame(soil_all_eigen.18)
df.soil_all_eigen.18$OTU <- rownames(df.soil_all_eigen.18)
df.soil_all_eigen.18 <- df.soil_all_eigen.18[,-c(2:22)]
write.table(df.soil_all_eigen.18,"df.soil_all_eigen.18.txt",sep="\t",quote=F)

n <- 1
soil_all_keystone.eigen.18 <- df.soil_all_eigen.18$OTU[df.soil_all_eigen.18$vector> quantile(df.soil_all_eigen.18$vector,prob=1-n/100)]
length(soil_all_keystone.eigen.18) #112


## Set node shape soil_all_keystone.2018
V(all_soil_net.18)$shape <- V(all_soil_net.18)$name
V(all_soil_net.18)$shape[V(all_soil_net.18)$shape %in% soil_all_keystone.between.18.score] <- "star"
#V(all_soil_net.18)$shape[V(all_soil_net.18)$shape %in% strictkeystone.18] <- "star"
V(all_soil_net.18)$shape[V(all_soil_net.18)$shape %in% rownames(otu_norm_soil_16s.18)] <- "circle"
V(all_soil_net.18)$shape[V(all_soil_net.18)$shape %in% rownames(otu_norm_soil_its.18)] <- "triangle"
V(all_soil_net.18)$shape[V(all_soil_net.18)$shape %in% rownames(otu_norm_soil_arch.18)] <- "square"

## Set node colors based upon sensitivity to management system
cs <- c("Conventional", "No_fertilizer", "No_pesticide")
unique(V(all_soil_net.18)$indicgroup)
V(all_soil_net.18)$color <- V(all_soil_net.18)$indicgroup
V(all_soil_net.18)$color[V(all_soil_net.18)$color == "NA"] <- "gray50"
V(all_soil_net.18)$color[V(all_soil_net.18)$color == "No_fertilizer"] <- "#0066CC"
V(all_soil_net.18)$color[V(all_soil_net.18)$color == "No_pesticide"] <- "#336633"
V(all_soil_net.18)$color[V(all_soil_net.18)$color == "Conventional"] <- "#CC9900"
V(all_soil_net.18)$frame.color <- V(all_soil_net.18)$color

#node attributed to important rOTU
soil_all_net_csnodes.18 <- rownames(nodeattrib_soil_combine.18[nodeattrib_soil_combine.18$indicgroup %in% cs,]) #421
soil_all_net_csnodes_bac.18 <- soil_all_net_csnodes.18[grep("^B",soil_all_net_csnodes.18)] #209
soil_all_net_csnodes_fun.18 <- soil_all_net_csnodes.18[grep("^F",soil_all_net_csnodes.18)] #199
soil_all_net_csnodes_arch.18 <- soil_all_net_csnodes.18[grep("^A",soil_all_net_csnodes.18)] #13

#Node size
V(all_soil_net.18)$size <- V(all_soil_net.18)$name
V(all_soil_net.18)$size[V(all_soil_net.18)$size %in% soil_all_keystone.between.18.score] <- 3
#V(all_soil_net.18)$size[V(all_soil_net.18)$size %in% strictkeystone.2018] <- 5
V(all_soil_net.18)$size[V(all_soil_net.18)$size %in% soil_all_net_csnodes_bac.18] <- 4
V(all_soil_net.18)$size[V(all_soil_net.18)$size %in% soil_all_net_csnodes_arch.18] <- 4
V(all_soil_net.18)$size[V(all_soil_net.18)$size %in% soil_all_net_csnodes_fun.18] <- 5
V(all_soil_net.18)$size[V(all_soil_net.18)$size %in% rownames(otu_norm_soil_16s.18)] <- 2
V(all_soil_net.18)$size[V(all_soil_net.18)$size %in% rownames(otu_norm_soil_its.18)] <- 3
V(all_soil_net.18)$size[V(all_soil_net.18)$size %in% rownames(otu_norm_soil_arch.18)] <- 2
soilcombine_nodesizes.18 <- as.numeric(V(all_soil_net.18)$size)

## Set edge color
E(all_soil_net.18)$color <- ifelse(E(all_soil_net.18)$cor >0, "#99CCFF","#FF6666")

##### Explore community structure of soil meta-network by defining modules #####

## Make vectors of network nodes responding to different cropping systems
CF_nodes_soil.18 <- rownames(nodeattrib_soil_combine.18[nodeattrib_soil_combine.18$indicgroup=="Conventional",])
NP_nodes_soil.18 <- rownames(nodeattrib_soil_combine.18[nodeattrib_soil_combine.18$indicgroup=="No_pesticide",])
NF_nodes_soil.18 <- rownames(nodeattrib_soil_combine.18[nodeattrib_soil_combine.18$indicgroup=="No_fertilizer",])

cs_nodes_soil_all.18 <- c(CF_nodes_soil.18,NF_nodes_soil.18,NP_nodes_soil.18)
Bcs_nodes_soil_all.18 <- cs_nodes_soil_all.18[grep("^B",cs_nodes_soil_all.18)] #181
Fcs_nodes_soil_all.18 <- cs_nodes_soil_all.18[grep("^F",cs_nodes_soil_all.18)] #129
Acs_nodes_soil_all.18 <- cs_nodes_soil_all.18[grep("^A",cs_nodes_soil_all.18)] #12

## Perform cluster analysis using greedy clustering algorithm 
cfg_soil.18 <- cluster_fast_greedy(as.undirected(all_soil_net.18))

## Subset for top 20 biggest nodes
soil_modules.18 <- sort(table(membership(cfg_soil.18)),decr=T)
soil_modules_20.18 <- soil_modules.18[1:20]

sum(soil_modules_20.18)/sum(soil_modules.18)
sm20_plot.18 <- soil_modules_20.18
names(sm20_plot.18) <- as.factor(1:20)
soil_modules_cs.18 <- table(factor(membership(cfg_soil.18)[cs_nodes_soil_all.18],levels=names(soil_modules.18)))
soil_modules_cs_20.18 <- soil_modules_cs.18[names(soil_modules_20.18)]
smcs20_plot.18 <- soil_modules_cs_20.18
names(smcs20_plot.18) <- as.factor(1:20)

## Make vector of nodes in top 20 modules
soil_modules_points.18 <- membership(cfg_soil.18)[membership(cfg_soil.18) %in% names(soil_modules_20.18)]

soil_points.18 <- NULL
for(i in soil_modules_points.18){
  soilx <- which(names(soil_modules_20.18)==i)
  soil_points.18 <- c(soil_points.18, soilx)
}
names(soil_points.18) <- names(soil_modules_points.18)

## Set node colors by cropping system sensitivity 
soil_all_cols.18 <- sort(soil_points.18)
soil_all_cols.18[!names(soil_all_cols.18) %in% cs] <- "gray50"
soil_all_cols.18[names(soil_all_cols.18) %in% CF_nodes_soil.18] <- "#CC9900"
soil_all_cols.18[names(soil_all_cols.18) %in% NF_nodes_soil.18] <- "#0066CC"
soil_all_cols.18[names(soil_all_cols.18) %in% NP_nodes_soil.18] <- "#336633"


soil_all_pch.18 <- sort(soil_points.18)
soil_all_pch.18[names(soil_all_pch.18) %in% rownames(otu_norm_soil_16s.18)] <- 1
soil_all_pch.18[names(soil_all_pch.18) %in% rownames(otu_norm_soil_its.18)] <- 2
soil_all_pch.18[names(soil_all_pch.18) %in% rownames(otu_norm_soil_arch.18)] <- 3
soil_all_pch.18[names(soil_all_pch.18) %in% intersect(rownames(otu_norm_soil_16s.18),cs_nodes_soil_all.18)] <- 16
soil_all_pch.18[names(soil_all_pch.18) %in% intersect(rownames(otu_norm_soil_its.18),cs_nodes_soil_all.18)] <- 16
soil_all_pch.18[names(soil_all_pch.18) %in% intersect(rownames(otu_norm_soil_arch.18),cs_nodes_soil_all.18)] <- 16
soil_all_pch.18[names(soil_all_pch.18) %in% soil_all_keystone.between.18.score] <- 8

soil_all_cex.18 <- sort(soil_points.18)
soil_all_cex.18[!names(soil_all_cex.18) %in% cs_nodes_soil_all.18] <- 1
soil_all_cex.18[names(soil_all_cex.18) %in% cs_nodes_soil_all.18] <- 2

soil_mods_list_cs.18 <- list()
for (i in names(soil_modules_cs_20.18)){
  x1 <- names(membership(cfg_soil.18)[membership(cfg_soil.18)==i])
  x2 <- x1[x1 %in% cs_nodes_soil_all.18]
  soil_mods_list_cs.18[[i]] <- as.numeric(V(all_soil_net.18)[x2])
}


##### meta co-occurrence networks #####
## Note: the permuations for the layout of the network can be very time consuming and processor intensive
# set.seed(8051)
coords_soil_all.18 <- layout_(all_soil_net.18, with_fr(niter=9999, grid="nogrid"))
write.table(coords_soil_all.18,"coords_soil_all.18.txt",sep="\t",row.names=F,col.names=F,quote=F)


## Import pre-calculated FR layout coordinates to save time 
coords_soil_all.18 <- as.matrix(read.table("coords_soil_all.18.txt"))
dimnames(coords_soil_all.18) <-  NULL

par(mfrow=c(1,1), mar=c(0,0,0,0))


soil_cols.18 <- c("lavenderblush1", "lavenderblush2","lavenderblush3","lavenderblush4")

plot(all_soil_net.18,vertex.label=NA,vertex.size=soilcombine_nodesizes.18, layout=coords_soil_all.18,
     mark.groups=list(soil_mods_list_cs.18$`6`, soil_mods_list_cs.18$`25`, soil_mods_list_cs.18$`26`),
     mark.col=soil_cols.18, mark.border=soil_cols.18)

library(plotrix)

##Distribution of rOTUs in modules
plot(sequence(sm20_plot.18)~jitter(sort(soil_points.18),1.75),pch=soil_all_pch.18, col=soil_all_cols.18, cex=soil_all_cex.18,
     ann=F,axes=F,ylim=c(0,800),xlim=c(0,20))
axis(1,at=1:20, labels=F, tcl=-0.15)
staxlab(side=1,at=1:20,labels=paste(rep("M",20),names(soil_modules_20.18),sep=""),srt=45)
axis(1,at=seq(1,20,1), tick=F, labels=paste(round(smcs20_plot.18/sm20_plot.18*100,1),"%",sep=""),line=1.5,cex.axis=0.7)
axis(2,at=seq(0,700,100))
text(-0.25,-50,"rOTUs\n in module:",xpd=NA)
title(ylab="Number of nodes per module",line=2)


dev.off()


##modified FR layout
# To do this, node names have to be changed into number
library(qgraph)

all_soil_net.18.integer <- set.vertex.attribute(all_soil_net.18, "name", value=paste(1:10765))
edge.18 <- get.edgelist(all_soil_net.18.integer)
edge.18.num<- apply(edge.18, 2, as.numeric)

layout.18.modified <- qgraph.layout.fruchtermanreingold(edge.18.num,vcount=vcount(all_soil_net.18.integer),
                                                        area=10*(vcount(all_soil_net.18.integer)^2),repulse.rad=(vcount(all_soil_net.18.integer)^3.1))

#without module
plot(all_soil_net.18.integer,layout=layout.18.modified,vertex.size=soilcombine_nodesizes.18,vertex.label=NA)
V(all_soil_net.18.integer)$color


#layout_sphere <- layout_on_sphere(all_soil_net.18.integer) 
plot(all_soil_net.18.integer,layout=layout_sphere,vertex.size=soilcombine_nodesizes.18,vertex.label=NA)

#with module
plot(all_soil_net.18.integer,layout=layout.18.modified,vertex.size=soilcombine_nodesizes.18,vertex.label=NA,
     mark.groups=list(soil_mods_list_cs.18$`6`, soil_mods_list_cs.18$`25`, soil_mods_list_cs.18$`26`),
     mark.col=soil_cols.18, mark.border=soil_cols.18)

legend("bottomleft",legend=c("Module 6", "Module 25", "Module 26"),col=soil_cols.18,
       bty="n",fill=soil_cols.18,border=soil_cols.18)


dev.off()



###Investigating the composition of modules
### plotting average module response
soil_all_cols.17[names(soil_all_cols.17) %in% CF_nodes_soil.17] <- "#CC9900"
soil_all_cols.17[names(soil_all_cols.17) %in% NF_nodes_soil.17] <- "#0066CC"
soil_all_cols.17[names(soil_all_cols.17) %in% NP_nodes_soil.17] <- "#336633"


soil_all_cols.18[names(soil_all_cols.18) %in% CF_nodes_soil.18] <- "#CC9900"
soil_all_cols.18[names(soil_all_cols.18) %in% NF_nodes_soil.18] <- "#0066CC"
soil_all_cols.18[names(soil_all_cols.18) %in% NP_nodes_soil.18] <- "#336633"


par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("#CC9900","#0066CC","#336633")
names(CS_cols) <- c("Conventional","No_fertilizer", "No_pesticide")

library(sciplot)
## SOIL module 1
meta.2017 <-sample_data(bac.clean.nolog.17.5site)
class(cfg_soil.17)

bargraph.CI(meta.2017$Cultural_practice, colSums(otu_norm_soil_combine.17[cfg_soil.17[[1]],]), 
            las=2, ylab="Cumulative normalized abundance", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 1", col=CS_cols, border=F)
## SOIL module 3
bargraph.CI(meta.2017$Cultural_practice, colSums(otu_norm_soil_combine.17[cfg_soil.17[[3]],]), 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 3", col=CS_cols, border=F)
## SOIL module 4
bargraph.CI(meta.2017$Cultural_practice, colSums(otu_norm_soil_combine.17[cfg_soil.17[[4]],]), 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 4", col=CS_cols, border=F)

## SOIL module 6
bargraph.CI(meta.2017$Cultural_practice, colSums(otu_norm_soil_combine.17[cfg_soil.17[[6]],]), 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 6", col=CS_cols, border=F)
## Legend
plot.new()
par(mar=c(0.5,0,2,0))
legend("left",  bty="n", cex=1, #x.intersp=0.1, y.intersp=1,
       legend=names(CS_cols), 
       fill=CS_cols, 
       border=CS_cols , xpd=T)



par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("#CC9900","#0066CC","#336633")
names(CS_cols) <- c("Conventional","No_fertilizer", "No_pesticide")

meta.2018 <-sample_data(bac.clean.nolog.18)

#SOIL module 6
bargraph.CI(meta.2018$Cultural_practice, colSums(otu_norm_soil_combine.18[cfg_soil.18[[6]],]), 
            las=2, ylab="Cumulative normalized abundance", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 6", col=CS_cols, border=F)
## SOIL module 25
bargraph.CI(meta.2018$Cultural_practice, colSums(otu_norm_soil_combine.18[cfg_soil.18[[25]],]), 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 25", col=CS_cols, border=F)
## SOIL module 26
bargraph.CI(meta.2018$Cultural_practice, colSums(otu_norm_soil_combine.18[cfg_soil.18[[26]],]), 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 26", col=CS_cols, border=F)

## Legend
plot.new()
par(mar=c(0.5,0,2,0))
legend("left",  bty="n", cex=1, #x.intersp=0.1, y.intersp=1,
       legend=names(CS_cols), 
       fill=CS_cols, 
       border=CS_cols , xpd=T)

dev.off()





#### taxonomies of module OTUs 

## defining bacteria and fungi of soil and root modules

# soil module 1
soil_M1 <- cfg_soil.17[[1]]
bacteria_in_soil_M1 <- soil_M1[grep("^B", as.vector(soil_M1))]
fungi_in_soil_M1 <- soil_M1[grep("^F", as.vector(soil_M1))]
archaea_in_soil_M1 <- soil_M1[grep("^A", as.vector(soil_M1))]
# soil module 3
soil_M3 <- cfg_soil.17[[3]]
bacteria_in_soil_M3 <- soil_M3[grep("^B", as.vector(soil_M3))]
fungi_in_soil_M3 <- soil_M3[grep("^F", as.vector(soil_M3))]
archaea_in_soil_M3 <- soil_M3[grep("^A", as.vector(soil_M3))]
# soil module 4
soil_M4 <- cfg_soil.17[[4]]
bacteria_in_soil_M4 <- soil_M4[grep("^B", as.vector(soil_M4))]
fungi_in_soil_M4 <- soil_M4[grep("^F", as.vector(soil_M4))]
archaea_in_soil_M4 <- soil_M4[grep("^A", as.vector(soil_M4))]

# soil module 6
soil_M6 <- cfg_soil.17[[6]]
bacteria_in_soil_M6 <- soil_M6[grep("^B", as.vector(soil_M6))]
fungi_in_soil_M6 <- soil_M6[grep("^F", as.vector(soil_M6))]
archaea_in_soil_M6 <- soil_M6[grep("^A", as.vector(soil_M6))]

#2018
# soil module 6
soil_M6 <- cfg_soil.18[[1]]
bacteria_in_soil_M6 <- soil_M6[grep("^B", as.vector(soil_M6))]
fungi_in_soil_M6 <- soil_M6[grep("^F", as.vector(soil_M6))]
archaea_in_soil_M6 <- soil_M6[grep("^A", as.vector(soil_M6))]
# soil module 25
soil_M25 <- cfg_soil.18[[3]]
bacteria_in_soil_M25 <- soil_M25[grep("^B", as.vector(soil_M25))]
fungi_in_soil_M25 <- soil_M25[grep("^F", as.vector(soil_M25))]
archaea_in_soil_M25 <- soil_M25[grep("^A", as.vector(soil_M25))]
# soil module 26
soil_M26 <- cfg_soil.18[[4]]
bacteria_in_soil_M26 <- soil_M26[grep("^B", as.vector(soil_M26))]
fungi_in_soil_M26 <- soil_M26[grep("^F", as.vector(soil_M26))]
archaea_in_soil_M26 <- soil_M26[grep("^A", as.vector(soil_M26))]

### counts of bacteria OTUs
bacteria_soil_M1 <- as.data.frame(table(tax_16s[bacteria_in_soil_M1, "labels"] ) )
colnames(bacteria_soil_M1) <- c("Class", "soil_M1")
bacteria_soil_M3 <- as.data.frame(table(tax_16s[bacteria_in_soil_M3, "labels"] ) )
colnames(bacteria_soil_M3) <- c("Class", "soil_M3")
bacteria_soil_M4 <- as.data.frame(table(tax_16s[bacteria_in_soil_M4, "labels"] ) )
colnames(bacteria_soil_M4) <- c("Class", "soil_M4")
bacteria_soil_M6 <- as.data.frame(table(tax_16s[bacteria_in_soil_M6, "labels"] ) )
colnames(bacteria_soil_M6) <- c("Class", "soil_M6")
bacteria_soil_modules <- merge(bacteria_soil_M1, bacteria_soil_M3, all=T, by="Class") 
bacteria_soil_modules <- merge(bacteria_soil_modules, bacteria_soil_M4, all=T, by="Class")
bacteria_soil_modules <- merge(bacteria_soil_modules, bacteria_soil_M6, all=T, by="Class") 
bacteria_soil_modules

bacteria_modules <- bacteria_soil_modules
bacteria_all_OTUs <- as.data.frame(table(tax_16s[, "labels"] ) )
colnames(bacteria_all_OTUs) <- c("Class", "all bOTUs")
bacteria_modules <- merge(bacteria_modules, bacteria_all_OTUs, all=T, by="Class") 
bacteria_modules

bacteria_modules_mat <- bacteria_modules[2:6]
rownames(bacteria_modules_mat) <- bacteria_modules$Class
bacteria_modules_mat[is.na(bacteria_modules_mat)] <- 0
colSums(bacteria_modules_mat)

bacteria_modules_prop <- t(t(bacteria_modules_mat)/colSums(bacteria_modules_mat) ) * 1
bacteria_modules_prop
colSums(bacteria_modules_prop)


### counts of fungi OTUs
fungi_soil_M1 <- as.data.frame(table(tax_its[fungi_in_soil_M1, "labels"] ) )
colnames(fungi_soil_M1) <- c("Class", "soil_M1")
fungi_soil_M3 <- as.data.frame(table(tax_its[fungi_in_soil_M3, "labels"] ) )
colnames(fungi_soil_M3) <- c("Class", "soil_M3")
fungi_soil_M4 <- as.data.frame(table(tax_its[fungi_in_soil_M4, "labels"] ) )
colnames(fungi_soil_M4) <- c("Class", "soil_M4")
fungi_soil_M6 <- as.data.frame(table(tax_its[fungi_in_soil_M6, "labels"] ) )
colnames(fungi_soil_M6) <- c("Class", "soil_M6")
fungi_soil_modules <- merge(fungi_soil_M1, fungi_soil_M3, all=T, by="Class") 
fungi_soil_modules <- merge(fungi_soil_modules, fungi_soil_M4, all=T, by="Class")
fungi_soil_modules <- merge(fungi_soil_modules, fungi_soil_M6, all=T, by="Class") 
fungi_soil_modules

fungi_modules <- fungi_soil_modules
fungi_all_OTUs <- as.data.frame(table(tax_its[, "labels"] ) )
colnames(fungi_all_OTUs) <- c("Class", "all fOTUs")
fungi_modules <- merge(fungi_modules, fungi_all_OTUs, all=T, by="Class") 
fungi_modules

fungi_modules_mat <- fungi_modules[2:5]
rownames(fungi_modules_mat) <- fungi_modules$Class
fungi_modules_mat[is.na(fungi_modules_mat)] <- 0
colSums(fungi_modules_mat)

fungi_modules_prop <- t(t(fungi_modules_mat)/colSums(fungi_modules_mat) ) * 1
fungi_modules_prop
colSums(fungi_modules_prop)



#Archaea
archaea_soil_M1 <- as.data.frame(table(tax_arch[archaea_in_soil_M1, "labels"] ) )
colnames(archaea_soil_M1) <- c("Class", "soil_M1")
archaea_soil_M3 <- as.data.frame(table(tax_arch[archaea_in_soil_M3, "labels"] ) )
colnames(archaea_soil_M3) <- c("Class", "soil_M3")
archaea_soil_M4 <- as.data.frame(table(tax_arch[archaea_in_soil_M4, "labels"] ) )
colnames(archaea_soil_M4) <- c("Class", "soil_M4")
archaea_soil_M6 <- as.data.frame(table(tax_arch[archaea_in_soil_M6, "labels"] ) )
colnames(archaea_soil_M6) <- c("Class", "soil_M6")
archaea_soil_modules <- merge(archaea_soil_M1, archaea_soil_M3, all=T, by="Class") 
archaea_soil_modules <- merge(archaea_soil_modules, archaea_soil_M4, all=T, by="Class")
archaea_soil_modules <- merge(archaea_soil_modules, archaea_soil_M6, all=T, by="Class") 
archaea_soil_modules

archaea_modules <- archaea_soil_modules
archaea_all_OTUs <- as.data.frame(table(tax_arch[, "labels"] ) )
colnames(archaea_all_OTUs) <- c("Class", "all aOTUs")
archaea_modules <- merge(archaea_modules, archaea_all_OTUs, all=T, by="Class") 
archaea_modules

archaea_modules_mat <- archaea_modules[2:6]
rownames(archaea_modules_mat) <- archaea_modules$Class
archaea_modules_mat[is.na(archaea_modules_mat)] <- 0
colSums(archaea_modules_mat)

archaea_modules_prop <- t(t(archaea_modules_mat)/colSums(archaea_modules_mat) ) * 1
archaea_modules_prop
colSums(archaea_modules_prop)



pdf(paste0(output,"Figure4c.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))

### bacteria
# PHYLA_label_cols_16s_legend
table(rownames(bacteria_modules_prop) %in% PHYLA_label_cols_16s$labels) 
bp <- barplot(cbind(bacteria_modules_prop[,1:4], NA, bacteria_modules_prop[,5]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_16s[rownames(bacteria_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(bacteria_modules_mat)[1:4], NA,
                       colSums(bacteria_modules_mat)[5]), xpd=T, cex=.6, las=2)

## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )

### Archaea
# PHYLA_label_cols_its_legend
par(mfrow=c(1,4), mar=c(2,4,1,0))
table(rownames(archaea_modules_prop) %in% PHYLA_label_cols_arch$labels) 

fp <- barplot(cbind(archaea_modules_prop[,1:4], NA, archaea_modules_prop[,5]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_arch[rownames(archaea_modules_prop),]$cols )
text(fp, 1.1, labels=c(colSums(archaea_modules_mat)[1:4], NA,
                       colSums(archaea_modules_mat)[5]), xpd=T, cex=.6, las=2)

plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_arch_legend$labels), 
       fill=rev(PHYLA_label_cols_arch_legend$cols), 
       border=rev(PHYLA_label_cols_arch_legend$cols) )

### Fungi
# PHYLA_label_cols_its_legend

par(mfrow=c(1,4), mar=c(2,4,1,0))
table(rownames(fungi_modules_prop) %in% PHYLA_label_cols_its$labels) 

fp <- barplot(cbind(fungi_modules_prop[,1:3], NA, fungi_modules_prop[,4]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_its[rownames(fungi_modules_prop),]$cols )
text(fp, 1.1, labels=c(colSums(fungi_modules_mat)[1:3], NA,
                       colSums(fungi_modules_mat)[4]), xpd=T, cex=.6, las=2)



## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=1, 
       legend=rev(PHYLA_label_cols_its_legend$labels),
       fill=rev(PHYLA_label_cols_its_legend$cols),
       border=rev(PHYLA_label_cols_its_legend$cols) )

dev.off()

##2018
#2018
# soil module 6
soil_M6 <- cfg_soil.18[[1]]
bacteria_in_soil_M6 <- soil_M6[grep("^B", as.vector(soil_M6))]
fungi_in_soil_M6 <- soil_M6[grep("^F", as.vector(soil_M6))]
archaea_in_soil_M6 <- soil_M6[grep("^A", as.vector(soil_M6))]
# soil module 25
soil_M25 <- cfg_soil.18[[3]]
bacteria_in_soil_M25 <- soil_M25[grep("^B", as.vector(soil_M25))]
fungi_in_soil_M25 <- soil_M25[grep("^F", as.vector(soil_M25))]
archaea_in_soil_M25 <- soil_M25[grep("^A", as.vector(soil_M25))]
# soil module 26
soil_M26 <- cfg_soil.18[[4]]
bacteria_in_soil_M26 <- soil_M26[grep("^B", as.vector(soil_M26))]
fungi_in_soil_M26 <- soil_M26[grep("^F", as.vector(soil_M26))]
archaea_in_soil_M26 <- soil_M26[grep("^A", as.vector(soil_M26))]

### counts of bacteria OTUs

bacteria_soil_M6 <- as.data.frame(table(tax_16s[bacteria_in_soil_M6, "labels"] ) )
colnames(bacteria_soil_M6) <- c("Class", "soil_M6")
bacteria_soil_M25 <- as.data.frame(table(tax_16s[bacteria_in_soil_M25, "labels"] ) )
colnames(bacteria_soil_M25) <- c("Class", "soil_M25")
bacteria_soil_M26 <- as.data.frame(table(tax_16s[bacteria_in_soil_M26, "labels"] ) )
colnames(bacteria_soil_M26) <- c("Class", "soil_M26")

bacteria_soil_modules <- merge(bacteria_soil_M6, bacteria_soil_M25, all=T, by="Class") 
bacteria_soil_modules <- merge(bacteria_soil_modules, bacteria_soil_M26, all=T, by="Class")
bacteria_soil_modules

bacteria_modules <- bacteria_soil_modules
bacteria_all_OTUs <- as.data.frame(table(tax_16s[, "labels"] ) )
colnames(bacteria_all_OTUs) <- c("Class", "all bOTUs")
bacteria_modules <- merge(bacteria_modules, bacteria_all_OTUs, all=T, by="Class") 
bacteria_modules

bacteria_modules_mat <- bacteria_modules[2:5]
rownames(bacteria_modules_mat) <- bacteria_modules$Class
bacteria_modules_mat[is.na(bacteria_modules_mat)] <- 0
colSums(bacteria_modules_mat)

bacteria_modules_prop <- t(t(bacteria_modules_mat)/colSums(bacteria_modules_mat) ) * 1
bacteria_modules_prop
colSums(bacteria_modules_prop)


### counts of fungi OTUs
fungi_soil_M6 <- as.data.frame(table(tax_its[fungi_in_soil_M6, "labels"] ) )
colnames(fungi_soil_M6) <- c("Class", "soil_M6")
fungi_soil_M25 <- as.data.frame(table(tax_its[fungi_in_soil_M25, "labels"] ) )
colnames(fungi_soil_M25) <- c("Class", "soil_M25")
fungi_soil_M26 <- as.data.frame(table(tax_its[fungi_in_soil_M26, "labels"] ) )
colnames(fungi_soil_M26) <- c("Class", "soil_M26")

fungi_soil_modules <- merge(fungi_soil_M6, fungi_soil_M25, all=T, by="Class") 
fungi_soil_modules <- merge(fungi_soil_modules, fungi_soil_M26, all=T, by="Class")
fungi_soil_modules

fungi_modules <- fungi_soil_modules
fungi_all_OTUs <- as.data.frame(table(tax_its[, "labels"] ) )
colnames(fungi_all_OTUs) <- c("Class", "all fOTUs")
fungi_modules <- merge(fungi_modules, fungi_all_OTUs, all=T, by="Class") 
fungi_modules

fungi_modules_mat <- fungi_modules[2:5]
rownames(fungi_modules_mat) <- fungi_modules$Class
fungi_modules_mat[is.na(fungi_modules_mat)] <- 0
colSums(fungi_modules_mat)

fungi_modules_prop <- t(t(fungi_modules_mat)/colSums(fungi_modules_mat) ) * 1
fungi_modules_prop
colSums(fungi_modules_prop)



#Archaea
archaea_soil_M6 <- as.data.frame(table(tax_arch[archaea_in_soil_M6, "labels"] ) )
colnames(archaea_soil_M6) <- c("Class", "soil_M6")
archaea_soil_M25 <- as.data.frame(table(tax_arch[archaea_in_soil_M25, "labels"] ) )
colnames(archaea_soil_M25) <- c("Class", "soil_M25")
archaea_soil_M26 <- as.data.frame(table(tax_arch[archaea_in_soil_M26, "labels"] ) )
colnames(archaea_soil_M26) <- c("Class", "soil_M26")

archaea_soil_modules <- merge(archaea_soil_M6, archaea_soil_M25, all=T, by="Class") 
archaea_soil_modules <- merge(archaea_soil_modules, archaea_soil_M26, all=T, by="Class")
archaea_soil_modules

archaea_modules <- archaea_soil_modules
archaea_all_OTUs <- as.data.frame(table(tax_arch[, "labels"] ) )
colnames(archaea_all_OTUs) <- c("Class", "all aOTUs")
archaea_modules <- merge(archaea_modules, archaea_all_OTUs, all=T, by="Class") 
archaea_modules

archaea_modules_mat <- archaea_modules[2:5]
rownames(archaea_modules_mat) <- archaea_modules$Class
archaea_modules_mat[is.na(archaea_modules_mat)] <- 0
colSums(archaea_modules_mat)

archaea_modules_prop <- t(t(archaea_modules_mat)/colSums(archaea_modules_mat) ) * 1
archaea_modules_prop
colSums(archaea_modules_prop)



pdf(paste0(output,"Figure4c.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))

### bacteria
# PHYLA_label_cols_16s_legend
table(rownames(bacteria_modules_prop) %in% PHYLA_label_cols_16s$labels) 
bp <- barplot(cbind(bacteria_modules_prop[,1:3], NA, bacteria_modules_prop[,4]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_16s[rownames(bacteria_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(bacteria_modules_mat)[1:3], NA,
                       colSums(bacteria_modules_mat)[4]), xpd=T, cex=.6, las=2)

## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )

### Archaea
# PHYLA_label_cols_its_legend
par(mfrow=c(1,4), mar=c(2,4,1,0))
table(rownames(archaea_modules_prop) %in% PHYLA_label_cols_arch$labels) 

fp <- barplot(cbind(archaea_modules_prop[,1:3], NA, archaea_modules_prop[,4]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_arch[rownames(archaea_modules_prop),]$cols )
text(fp, 1.1, labels=c(colSums(archaea_modules_mat)[1:3], NA,
                       colSums(archaea_modules_mat)[4]), xpd=T, cex=.6, las=2)

plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=0.75,
       legend=rev(PHYLA_label_cols_arch_legend$labels), 
       fill=rev(PHYLA_label_cols_arch_legend$cols), 
       border=rev(PHYLA_label_cols_arch_legend$cols) )

### Fungi
# PHYLA_label_cols_its_legend

par(mfrow=c(1,4), mar=c(2,4,1,0))
table(rownames(fungi_modules_prop) %in% PHYLA_label_cols_its$labels) 

fp <- barplot(cbind(fungi_modules_prop[,1:3], NA, fungi_modules_prop[,4]),
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_its[rownames(fungi_modules_prop),]$cols )
text(fp, 1.1, labels=c(colSums(fungi_modules_mat)[1:3], NA,
                       colSums(fungi_modules_mat)[4]), xpd=T, cex=.6, las=2)



## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.6, x.intersp=0.1, y.intersp=1, 
       legend=rev(PHYLA_label_cols_its_legend$labels),
       fill=rev(PHYLA_label_cols_its_legend$cols),
       border=rev(PHYLA_label_cols_its_legend$cols) )

dev.off()

