##Making order right
order <- c(lb <- c('Conventional', 'No_fertilizer', 'No_pesticide'))

##for full data
sample_data(bac.clean.ss)$Cultural_practice <- factor(sample_data(bac.clean.ss)$Cultural_practice, levels = order)
sample_data(bac.clean.ss)$Cultural_practice

library(plyr)
###Phylum level

##Phylum level
box_phylum_cultural <- function(phy.firstpart.5, phyl){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Phylum")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Phylum","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Phylum == phyl)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    #xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

box_phylum_cultural(bac.clean.ss, "Acidobacteria")
box_phylum_cultural(bac.clean.ss, "Proteobacteria")
box_phylum_cultural(bac.clean.ss, "Actinobacteria")
box_phylum_cultural(bac.clean.ss, "Chloroflexi")
box_phylum_cultural(bac.clean.ss, "Nitrospirae")
box_phylum_cultural(bac.clean.ss, "Verrucomicrobia")
box_phylum_cultural(bac.clean.ss, "Planctomycetes")

box_phylum_cultural(arch.clean.ss, "Euryarchaeota")
box_phylum_cultural(arch.clean.ss, "Crenarchaeota")


##Class level
box_class_cultural <- function(phy.firstpart.5, cla){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Class")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Class","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Class == cla)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}


box_class_cultural(bac.clean.ss,"Gammaproteobacteria")
box_class_cultural(fun.clean.ss,"Agaricomycetes")
box_class_cultural(fun.clean.ss,"Tremellomycetes")
box_class_cultural(fun.clean.ss,"Sordariomycetes")
box_class_cultural(fun.clean.ss,"Mortierellomycetes")

### Order level
box_order_cultural <- function(phy.firstpart.5, ord){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Order")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Order","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Order == ord)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

box_order_cultural(bac.clean.ss, "Rhizobiales")


### Family level
box_family_cultural <- function(phy.firstpart.5, fam){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Family")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Family","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Family == fam)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  
  return(p)
}

box_family_cultural(phy.clean.ss, "Rhizobiales")
box_family_cultural(arch.clean.ss, "Rice Cluster I")

library(plyr)
### Genus level
box_genus_cultural <- function(phy.firstpart.5, gen){
  
  good_phylum <-tax_glom(phy.firstpart.5, taxrank = "Genus")
  good_phylum
  
  ##  Melt it into a dataframe 
  goodsamps_phy_melt <- psmelt(good_phylum)  
  goodsamps_phy_melt
  sub_goodsamps_phy_melt <- subset(goodsamps_phy_melt, select = c("Sample", "Cultural_practice","Genus","Abundance"))
  TOTALS <- plyr::ddply(sub_goodsamps_phy_melt, c("Sample"), summarise, total = sum(Abundance))   
  sub_goodsamps_phy_melt
  ### Merge phylum sums and the total numbers -> so we can calculate the Relative Abundance
  sub_phy_melt_totals <- merge(sub_goodsamps_phy_melt, TOTALS, by = "Sample")  
  sub_phy_melt_totals
  ## Calculate the relative abundance
  sub_phy_melt_totals$RelAbundance <- sub_phy_melt_totals$Abundance/sub_phy_melt_totals$total  
  ##  Calculate the Percent Abundance
  sub_phy_melt_totals$PercentAbund <- sub_phy_melt_totals$RelAbundance * 1000 
  sub_phy_melt_totals
  
  sub_phy_melt_totals.proteo<- subset.data.frame(sub_phy_melt_totals, Genus == gen)
  sub_phy_melt_totals.proteo
  
  test <- aggregate(sub_phy_melt_totals.proteo$PercentAbund, by = list(sub_phy_melt_totals.proteo$Cultural_practice), max)
  
  colnames(test) <- c("Group", "MaxAbund")
  
  ##Kruskal-Wallis test
  kw<-kruskal.test(PercentAbund ~ Cultural_practice, data = sub_phy_melt_totals.proteo)
  kw$p.value
  kw$p.value<- round(kw$p.value, 4)
  
  #library(FSA)
  DT = dunnTest(PercentAbund ~ Cultural_practice,
                data=sub_phy_melt_totals.proteo,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  p <- ggplot(data = sub_phy_melt_totals.proteo, aes(x=Cultural_practice, y=PercentAbund)) + geom_boxplot(fill="white", width = 0.8) +
    theme_bw() +
    geom_point(position='jitter',shape=1, alpha=.5)+
    geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
    xlab(paste('Kruskal-Wallis, P-value = ', kw$p.value))+
    ylab("Relative Abundance () \n") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
  
  return(p)
}
oligotrphic <- c('Novosphingobium', "Candidatus Solibacter","Geobacter","Chryseobacterium","Segetibacter",'Paenibacillus', "Mucilaginibacter","Haliangium",'Pedobacter', 'Chitinophaga', 'Dyella', "Collimonas", "Labrys", "Bradyrhizobium", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Mesorhizobium", "Mycobacterium", "Rhodococcus","Terrabacter", "Bosea", "Delftia", "Afipia", "Nocardia", "Streptomyces", "Staphylococcus", "Kitasatospora", "Micromonospora", 'Kribella', 'Caulobacter', 'Variovorax')

box_genus_cultural(bac.clean.ss, "Haliangium")
box_genus_cultural(bac.clean.ss, "Nostoc")
box_genus_cultural(arch.clean.ss, "Rice Cluster I")
box_genus_cultural(fun.clean.ss, "Schizothecium")
