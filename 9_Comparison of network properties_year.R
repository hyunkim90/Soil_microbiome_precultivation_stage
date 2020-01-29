##2017 and 2018 metacommunity networks
meta_degree.17 <- sort(igraph::degree(all_soil_net.17,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.17))) #"max degree = " "241"     
print(c('mean degree = ', mean(meta_degree.17))) #"mean degree = "   "121.123298178332"
meta_degree.17 <- as.data.frame(meta_degree.17)
ggplot(meta_degree.17, aes(x=meta_degree.17)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.17, directed = FALSE))) #"average shortest path length =  6.2960242811269"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.17, "global"))) #"mean clustering coefficient = 0.951655932776188"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.17, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000334886796822388"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.17, normalized = TRUE)))) #"mean closeness centrality =  0.0108984265897049"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.17, V(all_soil_net.17)))))) #"mean number of neighbors =  121.123298178332"


#18
meta_degree.18 <- sort(igraph::degree(all_soil_net.18,mode="all"),decr=T)
print(c('max degree = ', max(meta_degree.18))) #"max degree = " "160"    
print(c('mean degree = ', mean(meta_degree.18))) #"mean degree = "   "81.9758476544357"
meta_degree.18 <- as.data.frame(meta_degree.18)
ggplot(meta_degree.18, aes(x=meta_degree.18)) + geom_density()

print(paste("average shortest path length = ", mean_distance(all_soil_net.18, directed = FALSE))) #"average shortest path length =  7.3653557114995"
print(paste("mean clustering coefficient = ", transitivity(all_soil_net.18, "global"))) #"mean clustering coefficient =  0.951848391184271"
print(paste("mean betweenness centrality = ", mean(betweenness(all_soil_net.18, directed = FALSE, normalized = TRUE)))) #"mean betweenness centrality =  0.000583744737167767"
print(paste("mean closeness centrality = ", mean(closeness(all_soil_net.18, normalized = TRUE)))) #"mean closeness centrality =  0.0128508095913441"
print(paste("mean number of neighbors = ", mean(lengths(adjacent_vertices(all_soil_net.18, V(all_soil_net.18)))))) #"mean number of neighbors =  81.9758476544357"



#2017
net <- all_soil_net.17
y2017_all_deg <- igraph::degree(net,mode="all")
y2017_all_betweenness <- betweenness(net, normalized = TRUE)
y2017_all_closeness <- closeness(net, normalized = TRUE)
y2017_all_transitivity <- transitivity(net, "local", vids = V(net))
names(y2017_all_transitivity)<- V(net)$name
y2017_all_transitivity[is.na(y2017_all_transitivity)] <- 0
## Bootstrapping degree
print(mean(y2017_all_deg)) #121.1233
set.seed(8046)
y2017_boot_degree = replicate(10000, sample(y2017_all_deg, 1, replace=TRUE))
print(mean(y2017_boot_degree)) #121.4705
## Bootstrapping betweenness
print(mean(y2017_all_betweenness)) #0.0003348868
set.seed(8046)
y2017_boot_betweenness = replicate(10000, sample(y2017_all_betweenness, 1, replace=TRUE))
print(mean(y2017_boot_betweenness)) #0.0003284825
## Bootstrapping closeness
print(mean(y2017_all_closeness)) #0.01089843
set.seed(8046)
y2017_boot_closeness = replicate(10000, sample(y2017_all_closeness, 1, replace=TRUE))
print(mean(y2017_boot_closeness)) #0.01089172
## Bootstrapping transitivity
print(mean(y2017_all_transitivity)) # 0.8426127
set.seed(8046)
y2017_boot_transitivity = replicate(10000, sample(y2017_all_transitivity, 1, replace=TRUE))
print(mean(y2017_boot_transitivity)) #0.8430568
y2017entional_node_characterstics<-cbind(y2017_boot_degree,y2017_boot_betweenness,y2017_boot_closeness,y2017_boot_transitivity)
write.table(y2017entional_node_characterstics, file="y2017_nodes_chars.csv", sep=",")


#2018
net <- all_soil_net.18
y2018_all_deg <- igraph::degree(net,mode="all")
y2018_all_betweenness <- betweenness(net, normalized = TRUE)
y2018_all_closeness <- closeness(net, normalized = TRUE)
y2018_all_transitivity <- transitivity(net, "local", vids = V(net))
names(y2018_all_transitivity)<- V(net)$name
y2018_all_transitivity[is.na(y2018_all_transitivity)] <- 0
## Bootstrapping degree
print(mean(y2018_all_deg)) #81.97585
set.seed(8046)
y2018_boot_degree = replicate(10000, sample(y2018_all_deg, 1, replace=TRUE))
print(mean(y2018_boot_degree)) #81.6968
## Bootstrapping betweenness
print(mean(y2018_all_betweenness)) #0.0005837447
set.seed(8046)
y2018_boot_betweenness = replicate(10000, sample(y2018_all_betweenness, 1, replace=TRUE))
print(mean(y2018_boot_betweenness)) #0.0006091443
## Bootstrapping closeness
print(mean(y2018_all_closeness)) # 0.01285081
set.seed(8046)
y2018_boot_closeness = replicate(10000, sample(y2018_all_closeness, 1, replace=TRUE))
print(mean(y2018_boot_closeness)) #0.01284943
## Bootstrapping transitivity
print(mean(y2018_all_transitivity)) #0.8350327
set.seed(8046)
y2018_boot_transitivity = replicate(10000, sample(y2018_all_transitivity, 1, replace=TRUE))
print(mean(y2018_boot_transitivity)) #0.8335679
y2018ization_node_characterstics<-cbind(y2018_boot_degree,y2018_boot_betweenness,y2018_boot_closeness,y2018_boot_transitivity)
write.table(y2018ization_node_characterstics, file="y2018_nodes_chars.csv", sep=",")



#Betweenness centrality
y2017_all_betweenness
y2018_all_betweenness


y2017.betw<-data.frame(y2017_all_betweenness)
y2017.betw$OTU_id <- rownames(y2017.betw)

y2018.betw<-data.frame(y2018_all_betweenness)
y2018.betw$OTU_id <- rownames(y2018.betw)


norm.betweenness <- data.frame(y2017.betw) %>% full_join(y2018.betw, by="OTU_id")
norm.betweenness[is.na(norm.betweenness)] <- 0
head(norm.betweenness)
class(norm.betweenness)
df.norm.betweenness <- data.frame(norm.betweenness)


##Degree
y2017_all_deg
y2018_all_deg

y2017degree<-data.frame(y2017_all_deg)
y2017degree$OTU_id <- rownames(y2017degree)

y2018degree<-data.frame(y2018_all_deg)
y2018degree$OTU_id <- rownames(y2018degree)

norm.degree <- data.frame(y2017degree) %>% full_join(y2018degree, by="OTU_id")
norm.degree[is.na(norm.degree)] <- 0

df.norm.degree <- data.frame(norm.degree)

df.network.property.year <- merge(df.norm.degree,df.norm.betweenness, by ="OTU_id")
head(df.network.property.year)


#Closeness
y2017_all_closeness
y2018_all_closeness

y2017.close<-data.frame(y2017_all_closeness)
y2017.close$OTU_id <- rownames(y2017.close)

y2018.close<-data.frame(y2018_all_closeness)
y2018.close$OTU_id <- rownames(y2018.close)


norm.closeness <- data.frame(y2017.close) %>% full_join(y2018.close, by="OTU_id")
norm.closeness[is.na(norm.closeness)] <- 0
head(norm.closeness)
class(norm.closeness)
df.norm.closeness <- data.frame(norm.closeness)

df.network.property.year <- merge(df.network.property.year,df.norm.closeness, by ="OTU_id")
head(df.network.property.year)

#transitivity
y2017_all_transitivity
y2018_all_transitivity

y2017.trans<-data.frame(y2017_all_transitivity)
y2017.trans$OTU_id <- rownames(y2017.trans)

y2018.trans<-data.frame(y2018_all_transitivity)
y2018.trans$OTU_id <- rownames(y2018.trans)


norm.transitivity <- data.frame(y2017.trans) %>% full_join(y2018.trans, by="OTU_id")
norm.transitivity[is.na(norm.transitivity)] <- 0
head(norm.transitivity)
class(norm.transitivity)
df.norm.transitivity <- data.frame(norm.transitivity)

df.network.property.year <- merge(df.network.property.year,df.norm.transitivity, by ="OTU_id")
head(df.network.property.year)
length(df.network.property.year$OTU_id)



##Supplementary Table_ rOTUs and network properties
df.network.property
df.network.property.year
head(df.network.property.year)
#Rename column names
names(df.network.property)[1]<- "OTU"
names(df.network.property.year)[1]<- "OTU"

list.total.rOTU<-read.xlsx("list_of_total_rOTU.xlsx",1)
head(list.total.rOTU)
length(rownames(list.total.rOTU))

#Merging rOTU table and network properties
list.total.rOTU.network <- merge(list.total.rOTU, df.network.property.year, by = "OTU", all = T)
list.total.rOTU.network <- merge(list.total.rOTU.network, df.network.property, by = "OTU", all = T)

head(list.total.rOTU.network)


write.table(list.total.rOTU.network, "List of rOTU with network properties_final_revised hub threshold.txt", sep = '\t', quote=F)


## Non rOTU list
list.non.rOTU.network<-read.table("list of non rOTU.txt",sep='\t', header = T)

list.non.rOTU.network <- list.non.rOTU.network[c(1)]
names(list.non.rOTU.network)[1]

tax_all_kingdom <- rbind(tax_16s, tax_its, tax_arch)

list.non.rOTU.network.w.taxonomy <- subset(tax_all_kingdom, rownames(tax_all_kingdom)%in%list.non.rOTU.network$OTU)
head(list.non.rOTU.network.w.taxonomy)

write.table(list.non.rOTU.network.w.taxonomy, "List of non rOTU with taxonomy.txt", sep = '\t', quote=F)


dev.off()



### Abundance of hub OTU

## Soil
soilall_cols <- sort(soil_all_deg)
soilall_cols[!names(soilall_cols) %in% indic_edge_soil_combine] <- "burlywood4"
soilall_cols[names(soilall_cols) %in% indic_edge_soil_combine] <- "orange"

soilall_cex <- sort(soil_all_deg)
soilall_cex[!names(soilall_cex) %in% indic_edge_soil_combine] <- 1
soilall_cex[names(soilall_cex) %in% indic_edge_soil_combine] <- 2

soilall_pch <- sort(soil_all_deg)
soilall_pch[names(soilall_pch) %in% rownames(otu_norm_soil_16s)] <- 16
soilall_pch[names(soilall_pch) %in% rownames(otu_norm_soil_its)] <- 17

# plotting function
otu_norm_soil_16s.17
otu_norm_soil_its.17
otu_norm_soil_arch.17

otu_norm_soil_16s.18
otu_norm_soil_its.18
otu_norm_soil_arch.18

soil_ra_16s <- apply(otu_norm_soil_16s.17,1,mean)
soil_ra_16s <- subset(soil_ra_16s, names(soil_ra_16s)%in%V(all_soil_net.17)$name)


soil_ra_its <- apply(otu_norm_soil_its.17,1,mean)
soil_ra_its <- subset(soil_ra_its, names(soil_ra_its)%in%V(all_soil_net.17)$name)

soil_ra_arch <- apply(otu_norm_soil_arch.17,1,mean)
soil_ra_arch <- subset(soil_ra_arch, names(soil_ra_arch)%in%V(all_soil_net.17)$name)

## 16S
scatterplot_mat <- cbind(abundance=soil_ra_16s,
                         node_degree=names(soil_ra_16s),
                         qualifier=names(soil_ra_16s),
                         taxon=rep("bacteria", length(soil_ra_16s)),
                         pch=16, cex=1.2)

## ITS
scatterplot_mat <- rbind(scatterplot_mat, cbind(abundance=soil_ra_its,
                                                node_degree=names(soil_ra_its),
                                                qualifier=names(soil_ra_its),
                                                taxon=rep("fungi",length(soil_ra_its)),
                                                pch=17, cex=1.2))

#Arch
scatterplot_mat <- rbind(scatterplot_mat, cbind(abundance=soil_ra_arch,
                                                node_degree=names(soil_ra_arch),
                                                qualifier=names(soil_ra_arch),
                                                taxon=rep("archaea",length(soil_ra_arch)),
                                                pch=15, cex=1.2))

length(names(y2017_all_closeness))


## add node degree
for(i in names(y2017_all_closeness)){
  scatterplot_mat[which(rownames(scatterplot_mat)==i), 2] <- y2017_all_closeness[i]
}
scatterplot_mat[which(scatterplot_mat[,2]=="0.011250875227395"),]



## add qualifier as colour
for(i in 1:length(V(all_soil_net.17)$name)){
  scatterplot_mat[which(rownames(scatterplot_mat)==V(all_soil_net.17)$name[i]), 3] <- V(all_soil_net.17)$color[i]
}

## add scaling factor by qualifier
scatterplot_mat[ ! scatterplot_mat[,"qualifier"] =="gray30"  ,"cex"] <- 1.7

## transform to dataframe
all_rownames <- rownames(scatterplot_mat)
scatterplot_mat <- as.data.frame(scatterplot_mat, stringsAsFactors=FALSE)
scatterplot_mat <- data.frame(lapply(scatterplot_mat, type.convert))
rownames(scatterplot_mat) <- all_rownames

##==============
## plot
##==============
plot.new()

## main plot
par(fig=c(0, 0.8, 0, 0.8), mar=c(4,4,0,0), new=T)

xrange <- c(0,max(scatterplot_mat$node_degree))
#xrange <-range(scatterplot_mat$node_degree)
#range(scatterplot_mat$node_degree)
yrange <- range(log(scatterplot_mat$abundance))


plot(c(), ylim=yrange, xlim=xrange, axes=F, bty="no", xlab="", ylab="")
rect(quantile(scatterplot_mat$node_degree, prob=1-n/100),
     min(log(scatterplot_mat$abundance)),
     max(scatterplot_mat$node_degree)+0.5,
     max(log(scatterplot_mat$abundance)), 
     col="lightgoldenrodyellow", border=NA)

# rect(min(scatterplot_mat$node_degree)-0.5,
#      quantile(log(scatterplot_mat$abundance),prob=1-n/100),
#      max(scatterplot_mat$node_degree),
#      max(log(scatterplot_mat$abundance)), 
#      col="lightgoldenrodyellow", border=NA)


par(fig=c(0, 0.8, 0, 0.8), mar=c(4,4,0,0), new=T)
plot(x=jitter(scatterplot_mat$node_degree, 1.1), 
     xlab="Closeness centrality",
     y=log(scatterplot_mat$abundance),
     ylab="Normalized abundance",
     col=as.character(scatterplot_mat$qualifier),
     pch=scatterplot_mat$pch,
     cex=scatterplot_mat$cex,
     ylim=yrange, xlim=xrange, bty="no", new=T)



## right plot
par(fig=c(0.8,1, 0,0.8), mar=c(4,1,0,2), new=T)
dens <- density(log(scatterplot_mat$abundance))
plot(x=dens$y, y=dens$x, ylim=yrange, type="l", col="transparent", axes=F, ann=F) ; par(new=T)
polygon(x=dens$y, y=dens$x, col="grey70", border=NA) ; par(new=T)

par(fig=c(0.81,.99, 0,0.8), mar=c(4,1,0,2), xpd=T, new=T)
ccsOTUs <- rownames(scatterplot_mat)[! scatterplot_mat$qualifier=="gray30"]
plot(y=log(scatterplot_mat[ccsOTUs,]$abundance),
     x=jitter(rep(1,length(ccsOTUs)),1.1),
     pch=scatterplot_mat[ccsOTUs,]$pch, 
     cex=scatterplot_mat[ccsOTUs,]$cex,
     col=as.character(scatterplot_mat[ccsOTUs,]$qualifier), #     col="grey30",
     ylim=yrange, axes=F, ann=F)


## upper plot
par(fig=c(0,0.8, 0.8,1), mar=c(0,4,3,0), new=T)
dens2 <- density(scatterplot_mat$node_degree)
plot(x=dens2$x, y=dens2$y, xlim=xrange, type="l", col="transparent", axes=F, ann=F)
polygon(x=dens2$x, y=dens2$y, col="grey70", border=NA) ; par(new=T)

par(fig=c(0,0.8, 0.81,0.99), mar=c(0,4,3,0), xpd=T, new=T)
plot(x=scatterplot_mat[ccsOTUs,]$node_degree,
     y=jitter(rep(1,length(ccsOTUs)),1.1),
     pch=scatterplot_mat[ccsOTUs,]$pch, 
     cex=scatterplot_mat[ccsOTUs,]$cex,
     col=as.character(scatterplot_mat[ccsOTUs,]$qualifier), #col="grey30",
     xlim=xrange, axes=F, ann=F)





source("node_degree_abundance_plot.R")

pdf(paste0(output,"Figure5b.pdf"), encoding="MacRoman", width=7, height=7)
node_degree_abundance_plot(bOTU_ra=root_ra_16s,
                           fOTU_ra=root_ra_its,
                           allOTU_node_degree=root_all_deg,
                           allOTU_vertices=all_root_net)
dev.off()

pdf(paste0(output,"Figure5a.pdf"), encoding="MacRoman", width=7, height=7)
node_degree_abundance_plot(bOTU_ra=soil_ra_16s,
                           fOTU_ra=soil_ra_its,
                           allOTU_node_degree=soil_all_deg,
                           allOTU_vertices=all_soil_net)
dev.off()


##Comparison of betweenness centrality between two years
df.norm.betweenness$y2017_all_betweenness

quantile(df.norm.betweenness$y2017_all_betweenness,prob=1-1/100) #0.4359687
quantile(df.norm.betweenness$y2018_all_betweenness,prob=1-1/100)

rownames(df.norm.betweenness)
df.norm.betweenness.plot<- df.norm.betweenness
head(df.norm.betweenness.plot)
names(df.norm.betweenness.plot)[2] <- "OTU"
df.norm.betweenness.plot$OTU

df.norm.betweenness.plot2 <- merge(df.norm.betweenness.plot, list.total.rOTU, by="OTU", all = T)
df.norm.betweenness.plot2$Bacteria <- ifelse(grepl("^B",df.norm.betweenness.plot2$OTU),'Bacteria',ifelse(grepl("^A",df.norm.betweenness.plot2$OTU),'Archaea','Fungi'))
df.norm.betweenness.plot2$index[is.na(df.norm.betweenness.plot2$index)] <- "Non"

df.norm.betweenness.plot2$rOTU <- ifelse(df.norm.betweenness.plot2$OTU%in%list.total.rOTU$OTU & df.norm.betweenness.plot2$index == 1, "CF", ifelse(df.norm.betweenness.plot2$OTU%in%list.total.rOTU$OTU & df.norm.betweenness.plot2$index == 2, "NF", ifelse(df.norm.betweenness.plot2$OTU%in%list.total.rOTU$OTU & df.norm.betweenness.plot2$index == 3, "NP", "Non-dif")))
head(df.norm.betweenness.plot2)

df.norm.betweenness.plot2$OTU[which(is.na(df.norm.betweenness.plot2$rOTU) == TRUE)]
df.norm.betweenness.plot2$OTU[which(is.na(df.norm.betweenness.plot2$Bacteria) == TRUE)]

df.norm.betweenness.plot2$Bacteria <-factor(df.norm.betweenness.plot2$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))
df.norm.betweenness.plot2$rOTU <-factor(df.norm.betweenness.plot2$rOTU, levels=c("CF", "NF","NP","Non-dif"))

write.table(df.norm.betweenness.plot2,"Hub properties plot.txt", sep ='\t', quote =F)


ggplot(df.norm.betweenness.plot2, aes(x=y2017_all_betweenness, y=y2018_all_betweenness)) +
  xlab('\n Betweenness centrality\n(2017)')+
  ylab("Betweenness centrality\n(2018)\n") +
  geom_point(aes(colour = rOTU, shape = Bacteria) ,size=5, alpha=0.9) +
  scale_colour_manual(labels = c('CF','NF','NP','Non-dif'), values = c("#CC9900", "#0066CC", '#336633','gray50'))+
  scale_shape_manual(values=c(16,15,17))+
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  geom_hline(yintercept=quantile(df.norm.betweenness$y2018_all_betweenness,prob=1-1/100) , color="maroon4", linetype='dotted')+
  geom_vline(xintercept=quantile(df.norm.betweenness$y2017_all_betweenness,prob=1-1/100), color="maroon4", linetype='dotted')+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  scale_x_continuous(breaks=seq(0,0.008,0.002))+
  scale_y_continuous(breaks=seq(0,0.02,0.005))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


dev.off()


##Comparison of closeness centrality between two years
max(df.norm.closeness$y2017_all_closeness) #0.01108051
max(df.norm.closeness$y2018_all_closeness) #0.01319871

quantile(df.norm.closeness$y2017_all_closeness,prob=1-1/100) #0.4359687
quantile(df.norm.closeness$y2018_all_closeness,prob=1-1/100)

rownames(df.norm.closeness)
df.norm.closeness.plot<- df.norm.closeness
head(df.norm.closeness.plot)
names(df.norm.closeness.plot)[2] <- "OTU"
df.norm.closeness.plot$OTU

df.norm.closeness.plot2 <- merge(df.norm.closeness.plot, list.total.rOTU, by="OTU", all = T)
df.norm.closeness.plot2$Bacteria <- ifelse(grepl("^B",df.norm.closeness.plot2$OTU),'Bacteria',ifelse(grepl("^A",df.norm.closeness.plot2$OTU),'Archaea','Fungi'))
df.norm.closeness.plot2$index[is.na(df.norm.closeness.plot2$index)] <- "Non"

df.norm.closeness.plot2$rOTU <- ifelse(df.norm.closeness.plot2$OTU%in%list.total.rOTU$OTU & df.norm.closeness.plot2$index == 1, "CF", ifelse(df.norm.closeness.plot2$OTU%in%list.total.rOTU$OTU & df.norm.closeness.plot2$index == 2, "NF", ifelse(df.norm.closeness.plot2$OTU%in%list.total.rOTU$OTU & df.norm.closeness.plot2$index == 3, "NP", "Non-dif")))
head(df.norm.closeness.plot2)

df.norm.closeness.plot2$OTU[which(is.na(df.norm.closeness.plot2$rOTU) == TRUE)]
df.norm.closeness.plot2$OTU[which(is.na(df.norm.closeness.plot2$Bacteria) == TRUE)]

df.norm.closeness.plot2$Bacteria <-factor(df.norm.closeness.plot2$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))
df.norm.closeness.plot2$rOTU <-factor(df.norm.closeness.plot2$rOTU, levels=c("CF", "NF","NP","Non-dif"))

write.table(df.norm.closeness.plot2,"Hub properties plot.txt", sep ='\t', quote =F)


ggplot(df.norm.closeness.plot2, aes(x=y2017_all_closeness, y=y2018_all_closeness)) +
  xlab('\n closeness centrality\n(2017)')+
  ylab("closeness centrality\n(2018)\n") +
  geom_point(aes(colour = rOTU, shape = Bacteria) ,size=5, alpha=0.9) +
  scale_colour_manual(labels = c('CF','NF','NP','Non-dif'), values = c("#CC9900", "#0066CC", '#336633','gray50'))+
  scale_shape_manual(values=c(16,15,17))+
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  geom_hline(yintercept=quantile(df.norm.closeness$y2018_all_closeness,prob=1-1/100) , color="maroon4", linetype='dotted')+
  geom_vline(xintercept=quantile(df.norm.closeness$y2017_all_closeness,prob=1-1/100), color="maroon4", linetype='dotted')+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  scale_x_continuous(breaks=seq(0,0.0115,0.002))+
  scale_y_continuous(breaks=seq(0,0.0132,0.002))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


dev.off()



### Degree
max(df.norm.degree$y2017_all_deg)
max(df.norm.degree$y2018_all_deg)

quantile(df.norm.degree$y2017_all_deg,prob=1-1/100) #0.4359687
quantile(df.norm.degree$y2018_all_deg,prob=1-1/100)

rownames(df.norm.degree)
df.norm.degree.plot<- df.norm.degree
head(df.norm.degree.plot)
names(df.norm.degree.plot)[2] <- "OTU"
df.norm.degree.plot$OTU

df.norm.degree.plot2 <- merge(df.norm.degree.plot, list.total.rOTU, by="OTU", all = T)
df.norm.degree.plot2$Bacteria <- ifelse(grepl("^B",df.norm.degree.plot2$OTU),'Bacteria',ifelse(grepl("^A",df.norm.degree.plot2$OTU),'Archaea','Fungi'))
df.norm.degree.plot2$index[is.na(df.norm.degree.plot2$index)] <- "Non"

df.norm.degree.plot2$rOTU <- ifelse(df.norm.degree.plot2$OTU%in%list.total.rOTU$OTU & df.norm.degree.plot2$index == 1, "CF", ifelse(df.norm.degree.plot2$OTU%in%list.total.rOTU$OTU & df.norm.degree.plot2$index == 2, "NF", ifelse(df.norm.degree.plot2$OTU%in%list.total.rOTU$OTU & df.norm.degree.plot2$index == 3, "NP", "Non-dif")))
head(df.norm.degree.plot2)

df.norm.degree.plot2$OTU[which(is.na(df.norm.degree.plot2$rOTU) == TRUE)]
df.norm.degree.plot2$OTU[which(is.na(df.norm.degree.plot2$Bacteria) == TRUE)]

df.norm.degree.plot2$Bacteria <-factor(df.norm.degree.plot2$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))
df.norm.degree.plot2$rOTU <-factor(df.norm.degree.plot2$rOTU, levels=c("CF", "NF","NP","Non-dif"))

write.table(df.norm.degree.plot2,"Hub properties plot.txt", sep ='\t', quote =F)


ggplot(df.norm.degree.plot2, aes(x=y2017_all_deg, y=y2018_all_deg)) +
  xlab('\n Degree\n(2017)')+
  ylab("Degree\n(2018)\n") +
  geom_point(aes(colour = rOTU, shape = Bacteria) ,size=5, alpha=0.9) +
  scale_colour_manual(labels = c('CF','NF','NP','Non-dif'), values = c("#CC9900", "#0066CC", '#336633','gray50'))+
  scale_shape_manual(values=c(16,15,17))+
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  geom_hline(yintercept=quantile(df.norm.degree$y2018_all_deg,prob=1-1/100) , color="maroon4", linetype='dotted')+
  geom_vline(xintercept=quantile(df.norm.degree$y2017_all_deg,prob=1-1/100), color="maroon4", linetype='dotted')+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  scale_x_continuous(breaks=seq(0,250,50))+
  scale_y_continuous(breaks=seq(0,160,40))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


dev.off()
#Edge list
#2017
E(all_soil_net.17)

#2018
E(all_soil_net.18)


edge.17.pos <- droplevels(all_cor_soil_df_padj.17[all_cor_soil_df_padj.17$cor>0,])
edge.17.neg <- droplevels(all_cor_soil_df_padj.17[all_cor_soil_df_padj.17$cor<0,])

edge.17.pos.combine<- edge.17.pos %>% tidyr::unite("from_to", from:to, remove = FALSE)
head(edge.17.pos.combine$from_to)
length(edge.17.pos.combine$from_to)#947475

edge.17.neg.combine<- edge.17.neg %>% tidyr::unite("from_to", from:to, remove = FALSE)
head(edge.17.neg.combine$from_to)
length(edge.17.neg.combine$from_to)#12




edge.18.pos <- droplevels(all_cor_soil_df_padj.18[all_cor_soil_df_padj.18$cor>0,])
edge.18.neg <- droplevels(all_cor_soil_df_padj.18[all_cor_soil_df_padj.18$cor<0,])

edge.18.pos.combine<- edge.18.pos %>% tidyr::unite("from_to", from:to, remove = FALSE)
head(edge.18.pos.combine$from_to)
length(edge.18.pos.combine$from_to)#441190

edge.18.neg.combine<- edge.18.neg %>% tidyr::unite("from_to", from:to, remove = FALSE)
head(edge.18.neg.combine$from_to)
length(edge.18.neg.combine$from_to)#45


intersect(edge.17.pos.combine$from_to, edge.18.pos.combine$from_to) #385
intersect(edge.17.neg.combine$from_to, edge.18.neg.combine$from_to) #5



intersect(names(V(all_soil_net.17)), names(V(all_soil_net.18))) #3671
