
## Ridge plot
bac.clean.ss.CJ.YS <- subset_samples(bac.clean.ss, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
bac.clean.ss.CJ.YS <- phyloseq::filter_taxa(bac.clean.ss.CJ.YS, function(x) sum(x) != 0, TRUE)

fun.clean.ss.CJ.YS <- subset_samples(fun.clean.ss, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
fun.clean.ss.CJ.YS <- phyloseq::filter_taxa(fun.clean.ss.CJ.YS, function(x) sum(x) != 0, TRUE)

arch.clean.ss.CJ.YS <- subset_samples(arch.clean.ss, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
arch.clean.ss.CJ.YS <- phyloseq::filter_taxa(arch.clean.ss.CJ.YS, function(x) sum(x) != 0, TRUE)


bac.clean.ss.MY.NJ <- subset_samples(bac.clean.ss, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
bac.clean.ss.MY.NJ <- phyloseq::filter_taxa(bac.clean.ss.MY.NJ, function(x) sum(x) != 0, TRUE)

fun.clean.ss.MY.NJ <- subset_samples(fun.clean.ss, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
fun.clean.ss.MY.NJ <- phyloseq::filter_taxa(fun.clean.ss.MY.NJ, function(x) sum(x) != 0, TRUE)

arch.clean.ss.MY.NJ <- subset_samples(arch.clean.ss, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
arch.clean.ss.MY.NJ <- phyloseq::filter_taxa(arch.clean.ss.MY.NJ, function(x) sum(x) != 0, TRUE)


## testing


## get imp.tax
get_imp_tax <- function(phy.clean.ss.5, classifier){
  imp <- data.frame(importance(classifier))
  imp <- tibble::rownames_to_column(imp, 'OTU')
  imp.desc <- imp %>% arrange(desc(MeanDecreaseGini))
  tax <- tax_table(phy.clean.ss.5)
  tax <- as.data.frame(tax)
  tax
  tax <- tibble::rownames_to_column(tax,'OTU')
  
  df.clean.ss.5 <- psmelt(phy.clean.ss.5)
  head(df.clean.ss.5)
  abun <- df.clean.ss.5 %>% group_by(OTU) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
  imp.tax <- left_join(imp.desc, tax, by=c('OTU', 'OTU')) %>% left_join(abun, by=c('OTU','OTU'))
  return(imp.tax)
}

#Subset of Conventional and no fertilization 
imp_tax_arch <- get_imp_tax(arch.clean.ss.CJ.YS, classifier_all.a)
imp_tax_fun <- get_imp_tax(fun.clean.ss.CJ.YS, classifier_all.f)
imp_tax_bac <- get_imp_tax(bac.clean.ss.CJ.YS, classifier_all)


#head(mydata)
imp_tax_bac.edit <- subset(imp_tax_bac, OTU != "result.all")
head(imp_tax_bac.edit)
imp.arrange.b <- imp_tax_bac.edit %>% arrange(desc(MeanDecreaseGini))
imp.arrange.b$rank <- paste0('RF',1:dim(imp_tax_bac.edit)[1])
another.b <- imp.arrange.b[c('OTU','rank')]
head(another.b)

imp_tax_fun.edit <- subset(imp_tax_fun, OTU != "result.all")

imp.arrange.f <- imp_tax_fun.edit %>% arrange(desc(MeanDecreaseGini))
imp.arrange.f$rank <- paste0('RF',1:dim(imp_tax_fun.edit)[1])
another.f <- imp.arrange.f[c('OTU','rank')]
head(another.f)


imp_tax_arch.edit <- subset(imp_tax_arch, OTU != "result.all")

imp.arrange.a <- imp_tax_arch.edit %>% arrange(desc(MeanDecreaseGini))
imp.arrange.a$rank <- paste0('RF',1:dim(imp_tax_arch.edit)[1])
another.a <- imp.arrange.a[c('OTU','rank')]


write.table(imp_tax_bac.edit, "imp_tax_bac.edit.oligo.tsv", sep='\t', quote =F)
write.table(imp_tax_fun.edit, "imp_tax_fun.edit.oligo.tsv", sep='\t', quote =F)
write.table(imp_tax_arch.edit, "imp_tax_arch.edit.oligo.tsv", sep='\t', quote =F)

##Construct the plot object

imp.total.arrange.b <- imp_tax_bac.edit %>% arrange(desc(total))
imp.total.arrange.b$number <- 1:dim(imp_tax_bac.edit)[1]
imp.total.arrange.b
num.b <- imp.total.arrange.b[c('OTU','number')]
num.b 
head(num.b)

imp.total.arrange.f <- imp_tax_fun.edit %>% arrange(desc(total))
imp.total.arrange.f$number <- 1:dim(imp_tax_fun.edit)[1]
imp.total.arrange.f
num.f <- imp.total.arrange.f[c('OTU','number')]
num.f 


imp.total.arrange.a <- imp_tax_arch.edit %>% arrange(desc(total))
imp.total.arrange.a$number <- 1:dim(imp_tax_arch.edit)[1]
imp.total.arrange.a
num.a <- imp.total.arrange.a[c('OTU','number')]
num.a

##Construct the plot object
#bac
Ab.b <- psmelt(bac.clean.ss.CJ.YS) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
Ta.b <- psmelt(bac.clean.ss.CJ.YS) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
Ta.b <- unique(Ta.b)

#archaea
Ab.a <- psmelt(arch.clean.ss.CJ.YS) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
Ta.a <- psmelt(arch.clean.ss.CJ.YS) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
Ta.a <- unique(Ta.a)

#fun
Ab.f <- psmelt(fun.clean.ss.CJ.YS) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
Ta.f <- psmelt(fun.clean.ss.CJ.YS) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
Ta.f <- unique(Ta.f)



## indexing the enrichment classification (rOTU)

imp.tax <- imp_tax_bac.edit
imp.tax <- imp_tax_arch.edit
imp.tax <- imp_tax_fun.edit

#Bacteria
resSig <- dplyr::left_join(Ab.b,Ta.b,by=c('OTU','OTU')) %>% dplyr::left_join(num.b,by='OTU')
imp.3 <- imp_tax_bac.edit[c('OTU','MeanDecreaseGini','total')]
resSig <- merge(resSig,imp.3,by='OTU')
resSig <- dplyr::left_join(resSig, another.b, by=('OTU'='OTU'))
head(resSig)
unique(resSig)
resSig.arrange <- resSig %>% arrange(desc(total))

OTU.list <- rbind(phy.list, arch.list, fun.list)
OTU_id.list <- OTU.list%>%group_by(OTU)%>%select('OTU', "OTU_id")

resSig <-merge(resSig, OTU_id.list, by = "OTU")
head(resSig)


resSig_bac.CJ.YS <-merge(resSig_bac.CJ.YS , OTU_id.list, by = "OTU")
length(colnames(resSig_bac.CJ.YS))


rotu.table <- resSig_bac.CJ.YS %>% group_by(OTU)%>%select('OTU_id', "Enriched")

rotu.table<-rotu.table[-c(1)]


mydata.b <- dplyr::left_join(resSig, rotu.table, by=('OTU_id'='OTU_id'))

sub.mydata.b <- mydata.b %>% arrange(desc(MeanDecreaseGini)) %>% head(35) #b:head(35), f:head(35), a:head(25)
unique(sub.mydata.b)



#Fungi
resSig <- dplyr::left_join(Ab.f,Ta.f,by=c('OTU','OTU')) %>% dplyr::left_join(num.f,by='OTU')
imp.3 <- imp_tax_fun.edit[c('OTU','MeanDecreaseGini','total')]
resSig <- merge(resSig,imp.3,by='OTU')
resSig <- dplyr::left_join(resSig, another.f, by=('OTU'='OTU'))
head(resSig)
unique(resSig)
resSig.arrange <- resSig %>% arrange(desc(total))

#OTU.list <- rbind(phy.list, arch.list, fun.list)
#OTU_id.list <- OTU.list%>%group_by(OTU)%>%select('OTU', "OTU_id")
resSig <-merge(resSig, OTU_id.list, by = "OTU")
head(resSig)

resSig_fun.CJ.YS <-merge(resSig_CJ.YS.fun , OTU_id.list, by = "OTU")

rotu.table <- resSig_fun.CJ.YS %>% group_by(OTU)%>%select('OTU_id', "Enriched")
rotu.table<-rotu.table[-c(1)]

#rotu.table <- rotu.final.list %>% group_by(OTU)%>%select('OTU', "rOTU")
#names(rotu.table)[1] <- "OTU_id"

mydata.f <- dplyr::left_join(resSig, rotu.table, by=('OTU_id'='OTU_id'))

sub.mydata.f <- mydata.f %>% arrange(desc(MeanDecreaseGini)) %>% head(35) #b:head(1000), f:head(80), a:head(64)
unique(sub.mydata.f$Enriched)


#Archaea
resSig <- dplyr::left_join(Ab.a,Ta.a,by=c('OTU','OTU')) %>% dplyr::left_join(num.a,by='OTU')
imp.3 <- imp_tax_arch.edit[c('OTU','MeanDecreaseGini','total')]
resSig <- merge(resSig,imp.3,by='OTU')
resSig <- dplyr::left_join(resSig, another.a, by=('OTU'='OTU'))
head(resSig)
unique(resSig)
resSig.arrange <- resSig %>% arrange(desc(total))

#OTU.list <- rbind(phy.list, arch.list, fun.list)
#OTU_id.list <- OTU.list%>%group_by(OTU)%>%select('OTU', "OTU_id")
resSig <-merge(resSig, OTU_id.list, by = "OTU")
head(resSig)
length(resSig$OTU)

resSig_arch.CJ.YS <-merge(resSig_CJ.YS.arch , OTU_id.list, by = "OTU")

rotu.table <- resSig_arch.CJ.YS %>% group_by(OTU)%>%select('OTU_id', "Enriched")
rotu.table<-rotu.table[-c(1)]


#rotu.table <- rotu.final.list %>% group_by(OTU)%>%select('OTU', "rOTU")
#names(rotu.table)[1] <- "OTU_id"

mydata.a <- dplyr::left_join(resSig, rotu.table, by=('OTU_id'='OTU_id'))

sub.mydata.a <- mydata.a %>% arrange(desc(MeanDecreaseGini)) %>% head(25) #b:head(1000), f:head(80), a:head(64)
unique(sub.mydata.a$Enriched)


# }

## get df.ridge
get_df_ridge <- function(phy.clean.ss.5, imp.tax.1, sub.mydata.c){
    df.otu <- phy.clean.ss.5 %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(Sample) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
  
  df.otu.rel
  imp.top100 <- imp.tax.1 %>% arrange(desc(MeanDecreaseGini)) %>% head(25) #b: 100, a: 64, f: 80
  
  imp100 <- imp.top100$OTU
  imp100
  
  imp.top100
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% imp100)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  str(df.ridge)
  b.thresh <- sub.mydata.c[c('OTU','Enriched','OTU_id')]
  b.thresh
  df.ridge_2 <- left_join(df.ridge, b.thresh, by=c('OTU','OTU'))
  
  unique(df.ridge$OTU)
  unique(b.thresh$OTU)
  
  df.ridge_2
  str(df.ridge_2)
  dim(df.ridge_2)
  43*20
  
  return(df.ridge_2)
}

df.ridge_2.b<- get_df_ridge(bac.clean.ss.CJ.YS, imp_tax_bac.edit, sub.mydata.b)

df.ridge_2.a<- get_df_ridge(arch.clean.ss.CJ.YS, imp_tax_arch.edit, sub.mydata.a)

df.ridge_2.f<- get_df_ridge(fun.clean.ss.CJ.YS, imp_tax_fun.edit, sub.mydata.f)

order.sample <- c("CJ1.A.1.1","CJ1.A.1.2","CJ1.A.1.3","CJ1.B.2.1","CJ1.B.2.2","CJ1.B.2.3","CJ1.C.3.1","CJ1.C.3.2","CJ1.C.3.3", 
                  "YS1.A.1.1","YS1.A.1.2","YS1.A.1.3","YS1.B.2.1","YS1.B.2.2","YS1.B.2.3","YS1.C.3.1","YS1.C.3.2","YS1.C.3.3",
                  "CJ2.A.1.1","CJ2.A.1.2","CJ2.A.1.3","CJ2.B.2.1","CJ2.B.2.2","CJ2.B.2.3","CJ2.C.3.1","CJ2.C.3.2","CJ2.C.3.3",
                  "YS2.A.1.1","YS2.A.1.2","YS2.A.1.3","YS2.B.2.1","YS2.B.2.2","YS2.B.2.3","YS2.C.3.1","YS2.C.3.2","YS2.C.3.3")
length(order.field )
df.ridge_2.b$Sample <- factor(df.ridge_2.b$Sample, levels=order.sample)
df.ridge_2.f$Sample <- factor(df.ridge_2.f$Sample, levels=order.sample)
df.ridge_2.a$Sample <- factor(df.ridge_2.a$Sample, levels=order.sample)

b.manipulate <- sub.mydata.b[c('rank','Enriched','OTU_id')]
write.xlsx(b.manipulate, "b.manipulate.xlsx")
b.manipulate <-read.xlsx("b.manipulate.xlsx",1)

b.id <- rev(b.manipulate$OTU_id)

a.manipulate <- sub.mydata.a[c('rank','Enriched','OTU_id')]
write.xlsx(a.manipulate, "a.manipulate.xlsx")
a.manipulate <-read.xlsx("a.manipulate.xlsx",1)
a.id <- rev(a.manipulate$OTU_id)

f.manipulate <- sub.mydata.f[c('rank','Enriched','OTU_id')]
write.xlsx(f.manipulate, "f.manipulate.xlsx")
f.manipulate <-read.xlsx("f.manipulate.xlsx",1)
f.id <- rev(f.manipulate$OTU_id)

df.ridge_2.b$OTU_id <- factor(df.ridge_2.b$OTU_id, levels=b.id)
df.ridge_2.b$Enriched <- factor(df.ridge_2.b$Enriched, levels=c('Non-fer', 'Fer', 'ns')) 

df.ridge_2.a$OTU_id <- factor(df.ridge_2.a$OTU_id, levels=a.id)
df.ridge_2.a$Enriched <- factor(df.ridge_2.a$Enriched, levels=c('Non-fer', 'Fer', 'ns')) 

df.ridge_2.f$OTU_id <- factor(df.ridge_2.f$OTU_id, levels=f.id)
df.ridge_2.f$Enriched <- factor(df.ridge_2.f$Enriched, levels=c('Non-fer', 'Fer', 'ns')) 

library(ggplot2)
library(ggridges)

plot_ridge <- function(df_ridge){
  ggplot(df_ridge, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=Enriched)) + 
    # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
    geom_density_ridges2(stat = "identity", scale=5, color='white',size=0.5)+
    xlab('\n Field')+
    ylab("Relative abundance \n") +
    theme(legend.text=element_text(size=12)) + 
    # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    scale_fill_manual(labels = c('Non-fer','Fer','ns'), values = c("Non-fer"= "Black", 'Fer'='Red','ns'='gray50'))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    
    geom_vline(xintercept=18.5, color="slategray3", linetype='dashed',size=1)
  
}

plot_ridge(df.ridge_2.b)
plot_ridge(df.ridge_2.f)
plot_ridge(df.ridge_2.a)

df.ridge_2$threshold


  
  
  write.xlsx(mydata.a, "Oligo_Archaea_RF.xlsx")
  write.xlsx(mydata.b, "Oligo_Bacteria_RF.xlsx")
  write.xlsx(mydata.f, "Oligo_Fungi_RF.xlsx")

  
  

  

  ###Copio
  ## get imp.tax
  get_imp_tax <- function(phy.clean.ss.5, classifier){
    imp <- data.frame(importance(classifier))
    imp <- tibble::rownames_to_column(imp, 'OTU')
    imp.desc <- imp %>% arrange(desc(MeanDecreaseGini))
    tax <- tax_table(phy.clean.ss.5)
    tax <- as.data.frame(tax)
    tax
    tax <- tibble::rownames_to_column(tax,'OTU')
    
    df.clean.ss.5 <- psmelt(phy.clean.ss.5)
    head(df.clean.ss.5)
    abun <- df.clean.ss.5 %>% group_by(OTU) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
    imp.tax <- left_join(imp.desc, tax, by=c('OTU', 'OTU')) %>% left_join(abun, by=c('OTU','OTU'))
    return(imp.tax)
  }
  
  #Subset of Conventional and no fertilization 
  imp_tax_arch <- get_imp_tax(arch.clean.ss.MY.NJ, classifier_all.a)
  imp_tax_fun <- get_imp_tax(fun.clean.ss.MY.NJ, classifier_all.f)
  imp_tax_bac <- get_imp_tax(bac.clean.ss.MY.NJ, classifier_all)
  
  
  #head(mydata)
  imp_tax_bac.edit.copio <- subset(imp_tax_bac, OTU != "result.all")
  head(imp_tax_bac.edit.copio)
  imp.arrange.b.copio <- imp_tax_bac.edit.copio %>% arrange(desc(MeanDecreaseGini))
  imp.arrange.b.copio$rank <- paste0('RF',1:dim(imp_tax_bac.edit.copio)[1])
  another.b.copio <- imp.arrange.b.copio[c('OTU','rank')]
  head(another.b.copio)
  
  write.table(imp_tax_bac.edit.copio, "imp_tax_bac.edit.copio_re.tsv", sep='\t', quote =F)
  
  imp_tax_fun.edit.copio <- subset(imp_tax_fun, OTU != "result.all")
  
  imp.arrange.f.copio <- imp_tax_fun.edit.copio %>% arrange(desc(MeanDecreaseGini))
  imp.arrange.f.copio$rank <- paste0('RF',1:dim(imp_tax_fun.edit.copio)[1])
  another.f.copio <- imp.arrange.f.copio[c('OTU','rank')]
  head(another.f.copio)
  write.table(imp_tax_fun.edit.copio, "imp_tax_fun.edit.copio_re.tsv", sep='\t', quote =F)
  
  imp_tax_arch.edit.copio <- subset(imp_tax_arch, OTU != "result.all")
  
  write.table(imp_tax_arch.edit.copio, "imp_tax_arch.edit.copio_re.tsv", sep='\t', quote =F)
  
  imp.arrange.a.copio <- imp_tax_arch.edit.copio %>% arrange(desc(MeanDecreaseGini))
  imp.arrange.a.copio$rank <- paste0('RF',1:dim(imp_tax_arch.edit.copio)[1])
  another.a.copio <- imp.arrange.a.copio[c('OTU','rank')]
  
  
  ##Construct the plot object
  
  imp.total.arrange.b.copio <- imp_tax_bac.edit.copio %>% arrange(desc(total))
  imp.total.arrange.b.copio$number <- 1:dim(imp_tax_bac.edit.copio)[1]
  imp.total.arrange.b.copio
  num.b.copio <- imp.total.arrange.b.copio[c('OTU','number')]
  num.b.copio 
  head(num.b.copio)
  
  imp.total.arrange.f.copio <- imp_tax_fun.edit.copio %>% arrange(desc(total))
  imp.total.arrange.f.copio$number <- 1:dim(imp_tax_fun.edit.copio)[1]
  imp.total.arrange.f.copio
  num.f.copio <- imp.total.arrange.f.copio[c('OTU','number')]
  num.f.copio 
  
  
  imp.total.arrange.a.copio <- imp_tax_arch.edit.copio %>% arrange(desc(total))
  imp.total.arrange.a.copio$number <- 1:dim(imp_tax_arch.edit.copio)[1]
  imp.total.arrange.a.copio
  num.a.copio <- imp.total.arrange.a.copio[c('OTU','number')]
  num.a.copio
  
  ##Construct the plot object
  #bac
  Ab.b.copio <- psmelt(bac.clean.ss.MY.NJ) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
  Ta.b.copio <- psmelt(bac.clean.ss.MY.NJ) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
  Ta.b.copio <- unique(Ta.b.copio)
  
  #archaea
  Ab.a.copio <- psmelt(arch.clean.ss.MY.NJ) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
  Ta.a.copio <- psmelt(arch.clean.ss.MY.NJ) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
  Ta.a.copio <- unique(Ta.a.copio)
  
  #fun
  Ab.f.copio <- psmelt(fun.clean.ss.MY.NJ) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
  Ta.f.copio <- psmelt(fun.clean.ss.MY.NJ) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
  Ta.f.copio <- unique(Ta.f.copio)
  
  
  
  ## indexing the enrichment classification (rOTU)
  
  imp.tax <- imp_tax_bac.edit.copio
  imp.tax <- imp_tax_arch.edit.copio
  imp.tax <- imp_tax_fun.edit.copio
  
  #Bacteria
  resSig <- dplyr::left_join(Ab.b.copio,Ta.b.copio,by=c('OTU','OTU')) %>% dplyr::left_join(num.b.copio,by='OTU')
  imp.3 <- imp_tax_bac.edit.copio[c('OTU','MeanDecreaseGini','total')]
  resSig <- merge(resSig,imp.3,by='OTU')
  resSig <- dplyr::left_join(resSig, another.b.copio, by=('OTU'='OTU'))
  head(resSig)
  unique(resSig)
  resSig.arrange <- resSig %>% arrange(desc(total))
  
  OTU.list <- rbind(phy.list, arch.list, fun.list)
  OTU_id.list <- OTU.list%>%group_by(OTU)%>%select('OTU', "OTU_id")
  
  resSig <-merge(resSig, OTU_id.list, by = "OTU")
  head(resSig)
  
  resSig_bac.MY.NJ <- read.table("Bacteria_daOTUs_Copio.tsv", sep = '\t', header = T)
  
  resSig_bac.MY.NJ <-merge(resSig_bac.MY.NJ , OTU_id.list, by = "OTU")
  length(colnames(resSig_bac.MY.NJ))
  names(resSig_bac.MY.NJ[15]) <- "OTU_id"
  head(resSig_bac.MY.NJ)

  rotu.table <- resSig_bac.MY.NJ %>% group_by(OTU)%>%select('OTU_id', "Enriched")
  
  rotu.table<-rotu.table[-c(1)]
 
  mydata.b.copio <- dplyr::left_join(resSig, rotu.table, by=('OTU_id'='OTU_id'))
  
  sub.mydata.b.copio <- mydata.b.copio %>% arrange(desc(MeanDecreaseGini)) %>% head(30) #b:head(35), f:head(35), a:head(25)
  unique(sub.mydata.b.copio)
  
  
  
  #Fungi
  resSig <- dplyr::left_join(Ab.f.copio,Ta.f.copio,by=c('OTU','OTU')) %>% dplyr::left_join(num.f.copio,by='OTU')
  imp.3 <- imp_tax_fun.edit.copio[c('OTU','MeanDecreaseGini','total')]
  resSig <- merge(resSig,imp.3,by='OTU')
  resSig <- dplyr::left_join(resSig, another.f.copio, by=('OTU'='OTU'))
  head(resSig)
  unique(resSig)
  resSig.arrange <- resSig %>% arrange(desc(total))
  
  #OTU.list <- rbind(phy.list, arch.list, fun.list)
  #OTU_id.list <- OTU.list%>%group_by(OTU)%>%select('OTU', "OTU_id")
  resSig <-merge(resSig, OTU_id.list, by = "OTU")
  head(resSig)
  
  resSig_fun.MY.NJ <- read.table("Fungi_daOTUs_Copio.tsv", sep = '\t', header = T)
  
  rotu.table <- resSig_fun.MY.NJ %>% group_by(OTU)%>%select('OTU_id', "Enriched")
  rotu.table<-rotu.table[-c(1)]
  
  #rotu.table <- rotu.f.copioinal.list %>% group_by(OTU)%>%select('OTU', "rOTU")
  #names(rotu.table)[1] <- "OTU_id"
  
  mydata.f.copio <- dplyr::left_join(resSig, rotu.table, by=('OTU_id'))
  
  sub.mydata.f.copio <- mydata.f.copio %>% arrange(desc(MeanDecreaseGini)) %>% head(30) #b:head(1000), f:head(80), a:head(64)
  unique(sub.mydata.f.copio$Enriched)
  
  
  #Archaea
  resSig <- dplyr::left_join(Ab.a.copio,Ta.a.copio,by=c('OTU','OTU')) %>% dplyr::left_join(num.a.copio,by='OTU')
  imp.3 <- imp_tax_arch.edit.copio[c('OTU','MeanDecreaseGini','total')]
  resSig <- merge(resSig,imp.3,by='OTU')
  resSig <- dplyr::left_join(resSig, another.a.copio, by=('OTU'='OTU'))
  head(resSig)
  unique(resSig)
  resSig.arrange <- resSig %>% arrange(desc(total))
  
  #OTU.list <- rbind(phy.list, arch.list, fun.list)
  #OTU_id.list <- OTU.list%>%group_by(OTU)%>%select('OTU', "OTU_id")
  resSig <-merge(resSig, OTU_id.list, by = "OTU")
  head(resSig)
  length(resSig$OTU)
  
  resSig_arch.MY.NJ <- read.table("Archaea_daOTUs_Copio.tsv", sep = '\t', header = T)
  
  rotu.table <- resSig_arch.MY.NJ %>% group_by(OTU)%>%select('OTU_id', "Enriched")
  rotu.table<-rotu.table[-c(1)]

  
  #rotu.table <- rotu.f.copioinal.list %>% group_by(OTU)%>%select('OTU', "rOTU")
  #names(rotu.table)[1] <- "OTU_id"
  
  mydata.a.copio <- dplyr::left_join(resSig, rotu.table, by=('OTU_id'='OTU_id'))
  mydata.a.copio$Enriched[is.na(mydata.a.copio$Enriched)] <- "ns"
  
  sub.mydata.a.copio <- mydata.a.copio %>% arrange(desc(MeanDecreaseGini)) %>% head(30) #b:head(1000), f:head(80), a:head(64)
  unique(sub.mydata.a.copio$Enriched)
  
  
  
  get_df_ridge <- function(phy.clean.ss.5, imp.tax.1, sub.mydata.c){
    df.otu <- phy.clean.ss.5 %>% psmelt()
    head(df.otu)
    # we need to group by samples
    df.otu.rel <- df.otu %>%  
      group_by(Sample) %>%                         # Filter out at absolute read of 20       
      mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
    
    df.otu.rel
    imp.top100 <- imp.tax.1 %>% arrange(desc(MeanDecreaseGini)) %>% head(30) #b: 90, a: 70, f: 30
    
    imp100 <- imp.top100$OTU
    imp100
    
    imp.top100
    
    df.selected.rel <- df.otu.rel %>% filter(OTU %in% imp100)
    df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
    str(df.ridge)
    b.thresh <- sub.mydata.c[c('OTU','Enriched','OTU_id')]
    b.thresh
    df.ridge_2 <- left_join(df.ridge, b.thresh, by=c('OTU','OTU'))
    
    unique(df.ridge$OTU)
    unique(b.thresh$OTU)
    
    df.ridge_2
    str(df.ridge_2)
    dim(df.ridge_2)
    43*20
    
    return(df.ridge_2)
  }
  
  df.ridge_2.b.copio<- get_df_ridge(bac.clean.ss.MY.NJ, imp_tax_bac.edit.copio, sub.mydata.b.copio)
  
  df.ridge_2.a.copio<- get_df_ridge(arch.clean.ss.MY.NJ, imp_tax_arch.edit.copio, sub.mydata.a.copio)
  
  df.ridge_2.f.copio<- get_df_ridge(fun.clean.ss.MY.NJ, imp_tax_fun.edit.copio, sub.mydata.f.copio)
  
  order.sample.copio <- c("MY2.A.1.1","MY2.A.1.2","MY2.A.1.3","MY2.B.2.1","MY2.B.2.2","MY2.B.2.3","MY2.C.3.1","MY2.C.3.2","MY2.C.3.3",
                     "NJ2.A.1.1","NJ2.A.1.2","NJ2.A.1.3","NJ2.B.2.1","NJ2.B.2.2","NJ2.B.2.3","NJ2.C.3.1","NJ2.C.3.2","NJ2.C.3.3",
                    "MY1.A.1.1","MY1.A.1.2","MY1.A.1.3","MY1.B.2.1","MY1.B.2.2","MY1.B.2.3","MY1.C.3.1","MY1.C.3.2","MY1.C.3.3", 
                    "NJ1.A.1.1","NJ1.A.1.2","NJ1.A.1.3","NJ1.B.2.1","NJ1.B.2.2","NJ1.B.2.3","NJ1.C.3.1","NJ1.C.3.2","NJ1.C.3.3")
 
  df.ridge_2.b.copio$Sample <- factor(df.ridge_2.b.copio$Sample, levels=order.sample.copio)
  df.ridge_2.f.copio$Sample <- factor(df.ridge_2.f.copio$Sample, levels=order.sample.copio)
  df.ridge_2.a.copio$Sample <- factor(df.ridge_2.a.copio$Sample, levels=order.sample.copio)
  
  b.manipulate_copio <- sub.mydata.b.copio[c('rank','Enriched','OTU_id')]
  write.xlsx(b.manipulate_copio, "b.manipulate_copio_re.xlsx")
  b.manipulate_copio <-read.xlsx("b.manipulate_copio.xlsx",1)
  
  b.id_copio <- rev(b.manipulate_copio$OTU_id)
  
  a.manipulate_copio <- sub.mydata.a.copio[c('rank','Enriched','OTU_id')]
  write.xlsx(a.manipulate_copio, "a.manipulate_copio_re.xlsx")
  a.manipulate_copio <-read.xlsx("a.manipulate_copio.xlsx",1)
  a.id_copio <- rev(a.manipulate_copio$OTU_id)
  
  f.manipulate_copio <- sub.mydata.f.copio[c('rank','Enriched','OTU_id')]
  write.xlsx(f.manipulate_copio, "f.manipulate_copio_re.xlsx")
  f.manipulate_copio <-read.xlsx("f.manipulate_copio.xlsx",1)
  f.id_copio <- rev(f.manipulate_copio$OTU_id)
  
  df.ridge_2.b.copio$OTU_id <- factor(df.ridge_2.b.copio$OTU_id, levels=b.id_copio)
  df.ridge_2.b.copio$Enriched <- factor(df.ridge_2.b.copio$Enriched, levels=c('Non-fer', 'Fer', 'ns')) 
  
  df.ridge_2.a.copio$OTU_id <- factor(df.ridge_2.a.copio$OTU_id, levels=a.id_copio)
  df.ridge_2.a.copio$Enriched <- factor(df.ridge_2.a.copio$Enriched, levels=c('Non-fer', 'Fer', 'ns')) 
  
  df.ridge_2.f.copio$OTU_id <- factor(df.ridge_2.f.copio$OTU_id, levels=f.id_copio)
  df.ridge_2.f.copio$Enriched <- factor(df.ridge_2.f.copio$Enriched, levels=c('Non-fer', 'Fer', 'ns')) 
  
  library(ggplot2)
  library(ggridges)
  
  plot_ridge <- function(df_ridge){
    ggplot(df_ridge, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=Enriched)) + 
      # geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)
      geom_density_ridges2(stat = "identity", scale=5, color='white',size=0.5)+
      xlab('\n Field')+
      ylab("Relative abundance \n") +
      theme(legend.text=element_text(size=12)) + 
      # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
      theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
      theme(legend.title=element_blank()) +
      guides(colour = guide_legend(override.a.copioes = list(size=8)))+
      guides(size=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
      theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
      theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
      theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
      scale_fill_manual(labels = c('Non-fer','Fer','ns'), values = c("Non-fer"= "Black", 'Fer'='Red','ns'='gray50'))+
      theme(panel.grid.major = element_blank())+
      theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
      
      geom_vline(xintercept=18.5, color="slategray3", linetype='dashed',size=1)
    
  }
  
  plot_ridge(df.ridge_2.b.copio)
  plot_ridge(df.ridge_2.f.copio)
  plot_ridge(df.ridge_2.a.copio)
  
  
  write.xlsx(mydata.a.copio, "Copio_Archaea_RF_re.xlsx")
  write.xlsx(mydata.b.copio, "Copio_Bacteria_RF_re.xlsx")
  write.xlsx(mydata.f.copio, "Copio_Fungi_RF_re.xlsx")

  write.xlsx(mydata.a, "Oligo_Archaea_RF.xlsx")
  write.xlsx(mydata.b, "Oligo_Bacteria_RF.xlsx")
  write.xlsx(mydata.f, "Oligo_Fungi_RF.xlsx")
  
  
  ## top 30
  sub.mydata.b.copio$OTU_id[which(sub.mydata.b.copio$OTU %in% rownames(tab_ra.b.SOM.TN.oligo))]  #6
  sub.mydata.b.copio$OTU_id[which(sub.mydata.b.copio$OTU %in% rownames(tab_ra.b.SOM.TN.copi))]
  

  sub.mydata.f.copio$OTU_id[which(sub.mydata.f.copio$OTU %in% rownames(tab_ra.f.SOM.TN.oligo))]  #6
  sub.mydata.f.copio$OTU_id[which(sub.mydata.f.copio$OTU %in% rownames(tab_ra.f.SOM.TN.copi))]
  
  sub.mydata.a.copio$OTU_id[which(sub.mydata.a.copio$OTU %in% rownames(tab_ra.a.SOM.TN.oligo))]  #6
  sub.mydata.a.copio$OTU_id[which(sub.mydata.a.copio$OTU %in% rownames(tab_ra.a.SOM.TN.copi))]
  
  sub.mydata.b$OTU_id[which(sub.mydata.b$OTU %in% rownames(tab_ra.b.SOM.TN.oligo))]  #6
  sub.mydata.b$OTU_id[which(sub.mydata.b$OTU %in% rownames(tab_ra.b.SOM.TN.copi))]
  
  
  sub.mydata.f$OTU_id[which(sub.mydata.f$OTU %in% rownames(tab_ra.f.SOM.TN.oligo))]  #6
  sub.mydata.f$OTU_id[which(sub.mydata.f$OTU %in% rownames(tab_ra.f.SOM.TN.copi))]
  
  sub.mydata.a$OTU_id[which(sub.mydata.a$OTU %in% rownames(tab_ra.a.SOM.TN.oligo))]  #6
  sub.mydata.a$OTU_id[which(sub.mydata.a$OTU %in% rownames(tab_ra.a.SOM.TN.copi))]
  
  
  
  mydata.a<- read.xlsx("Oligo_Archaea_RF.xlsx",1)
  mydata.b<-read.xlsx("Oligo_Bacteria_RF.xlsx",1)
  mydata.f<-read.xlsx("Oligo_Fungi_RF.xlsx",1)
  
  sub.mydata.b <- mydata.b %>% arrange(desc(MeanDecreaseGini)) %>% head(30) #b:head(35), f:head(35), a:head(25)
  unique(sub.mydata.b)
  sub.mydata.a <- mydata.a %>% arrange(desc(MeanDecreaseGini)) %>% head(30) #b:head(35), f:head(35), a:head(25)
  unique(sub.mydata.a)
  sub.mydata.f <- mydata.f %>% arrange(desc(MeanDecreaseGini)) %>% head(30) #b:head(35), f:head(35), a:head(25)
  unique(sub.mydata.f)
  

  #in top 30
  intersect(sub.mydata.b$OTU, rownames(tab_ra.b.SOM.TN.oligo))  #18 #16
  intersect(sub.mydata.b$OTU, rownames(tab_ra.b.SOM.TN.copi)) #8 #7
  
  
  intersect(sub.mydata.f$OTU, rownames(tab_ra.f.SOM.TN.oligo))  #6 #6
  intersect(sub.mydata.f$OTU, rownames(tab_ra.f.SOM.TN.copi)) #18 #16
  
  
  intersect(sub.mydata.a$OTU, rownames(tab_ra.a.SOM.TN.oligo))  #16 #18
  intersect(sub.mydata.a$OTU, rownames(tab_ra.a.SOM.TN.copi)) #3 #5
  
  
  
  write.table(resSig_bac.CJ.YS, "Bacteria_daOTUs_Oligo.tsv", sep = '\t', quote = F)
  write.table(resSig_fun.CJ.YS, "Fungi_daOTUs_Oligo.tsv", sep = '\t', quote = F)
  write.table(resSig_arch.CJ.YS, "Archaea_daOTUs_Oligo.tsv", sep = '\t', quote = F)
  
  
  write.table(resSig_bac.MY.NJ, "Bacteria_daOTUs_Copio.tsv", sep = '\t', quote = F)
  write.table(resSig_fun.MY.NJ, "Fungi_daOTUs_Copio.tsv", sep = '\t', quote = F)
  write.table(resSig_arch.MY.NJ, "Archaea_daOTUs_Copio.tsv", sep = '\t', quote = F)
  
  
  
## Combined ridge plot and phylogenetic tree of bacteria, archaea, and fungi
  sub.mydata.b
  library(seqinr)
  
  pop_taxa_wanted = function(physeq, badTaxa){
    allTaxa = taxa_names(physeq)
    myTaxa <- allTaxa[(allTaxa %in% badTaxa)]
    return(prune_taxa(myTaxa, physeq))
  }
  
  
  phy.clean.ss <- merge_phyloseq(phy.clean.ss, phy_tree(phy_tree))
  fun.clean.ss <- merge_phyloseq(fun.clean.ss, phy_tree(fun_tree))
  
  bac.clean.ss.tree <- subset_taxa(phy.clean.ss, Kingdom == "Bacteria")
  arch.clean.ss.tree <- subset_taxa(phy.clean.ss, Kingdom == "Archaea")
  
  bac.clean.ss.tree.sig <- pop_taxa_wanted(bac.clean.ss.tree, sub.mydata.b$OTU)
plot_tree(bac.clean.ss.tree.sig)


bac_oligo_phy_tree <- phy_tree(bac.clean.ss.tree.sig)

bac_oligo_tax <- tax_table(bac.clean.ss.tree.sig)
bac_oligo_tax <- data.frame(bac_oligo_tax, stringsAsFactors = F)
bac_oligo_tax$Phylum2 <- bac_oligo_tax$Phylum

#bac_oligo_tax[ rownames(bac_oligo_tax)[bac_oligo_tax$Class=="Alphaproteobacteria" ], ]$labels <- "Alphaproteobacteria"
#bac_oligo_tax[ rownames(bac_oligo_tax)[bac_oligo_tax$Class=="Betaproteobacteria" ], ]$labels <- "Betaproteobacteria"
bac_oligo_tax[ rownames(bac_oligo_tax)[bac_oligo_tax$Class=="Gammaproteobacteria" ], ]$Phylum2 <- "Gammaproteobacteria"
bac_oligo_tax[ rownames(bac_oligo_tax)[bac_oligo_tax$Class=="Deltaproteobacteria" ], ]$Phylum2 <- "Deltaproteobacteria"
table(bac_oligo_tax$Phylum2)


phy.colors <- c("Gammaproteobacteria" = "palegreen2","Deltaproteobacteria"= "palegreen3", "Actinobacteria"= "indianred2",
            "Bacteroidetes"="steelblue1", "Firmicutes"="tan1", "Acidobacteria"="lightsalmon4", "Chloroflexi" = "gold1",
            "Verrucomicrobia" = "orchid3", "Nitrospirae"= "palevioletred2", "Gemmatimonadetes"="peachpuff3")

bac_oligo_tax.2 <- bac_oligo_tax

tree.plot <- ggtree(bac_oligo_phy_tree, ladderize = F,branch.length="none") %<+% bac_oligo_tax +
  geom_tippoint(aes(color = Phylum2)) +  scale_color_manual(values = phy.colors)+
  geom_tiplab(aes(color = Phylum2), align = T, size = 0, alpha = 0.2)+
  theme(legend.position = "left")

  plot_data_tree<-tree.plot$data
  plot_data_tree.t<- plot_data_tree[-c(10:17)]
  bac_oligo_tax.2$label <- rownames(bac_oligo_tax.2)
  plot_data_tree.t <- merge(plot_data_tree.t, bac_oligo_tax.2, by = "label", all = T)
  plot_data_tree.t<- plot_data_tree.t %>% arrange(node)

  tree.plot$data<- plot_data_tree.t
  
  tree.plot
  
tax_order_bac.oligo <- data.frame(variable = as.character(bac_oligo_phy_tree$tip.label), 
                        order = 1:length(bac_oligo_phy_tree$tip.label)) # reverse order

## oligo archaea tree
arch.clean.ss.tree <- subset_taxa(phy.clean.ss, Kingdom == "Archaea")

arch.clean.ss.tree.sig <- pop_taxa_wanted(arch.clean.ss.tree, sub.mydata.a$OTU)
plot_tree(arch.clean.ss.tree.sig)


arch_oligo_phy_tree <- phy_tree(arch.clean.ss.tree.sig)

arch_oligo_tax <- tax_table(arch.clean.ss.tree.sig)
arch_oligo_tax <- data.frame(arch_oligo_tax, stringsAsFactors = F)
arch_oligo_tax$Phylum2 <- arch_oligo_tax$Phylum
unique(arch_oligo_tax$Phylum2)

arch.colors<- c("Euryarchaeota" = "#7B55E1","Crenarchaeota"= "#FF9999", "Thaumarchaeota"= "#CC6600",
                "Nanoarchaeaeota"="#999933", "Altiarchaeota"="gray50")

arch_oligo_tax.2 <- arch_oligo_tax

tree.plot <- ggtree(arch_oligo_phy_tree, ladderize = F,branch.length="none") %<+% arch_oligo_tax +
  geom_tippoint(aes(color = Phylum2)) +  scale_color_manual(values = arch.colors)+
  geom_tiplab(aes(color = Phylum2), align = T, size = 0, alpha = 0.2)+
  theme(legend.position = "left")

plot_data_tree<-tree.plot$data
plot_data_tree.t<- plot_data_tree[-c(10:17)]
arch_oligo_tax.2$label <- rownames(arch_oligo_tax.2)
plot_data_tree.t <- merge(plot_data_tree.t, arch_oligo_tax.2, by = "label", all = T)
plot_data_tree.t<- plot_data_tree.t %>% arrange(node)

tree.plot$data<- plot_data_tree.t

tree.plot

tax_order_arch.oligo <- data.frame(variable = as.character(arch_oligo_phy_tree$tip.label), 
                                  order = 1:length(arch_oligo_phy_tree$tip.label)) # reverse order


## oligo fungi tree
fun.clean.ss.tree.sig <- pop_taxa_wanted(fun.clean.ss, sub.mydata.f$OTU)
plot_tree(fun.clean.ss.tree.sig)


fun_oligo_phy_tree <- phy_tree(fun.clean.ss.tree.sig)

fun_oligo_tax <- tax_table(fun.clean.ss.tree.sig)
fun_oligo_tax <- data.frame(fun_oligo_tax, stringsAsFactors = F)
fun_oligo_tax$Phylum2 <- fun_oligo_tax$Phylum
fun_oligo_tax$Phylum2[is.na(fun_oligo_tax$Phylum2)] <- "unidentified"


unique(fun_oligo_tax$Phylum2)

fun.colors<- c("Ascomycota" = "#003366","Basidiomycota"= "#C84248", "Mortierellomycota"= "#EEB422",
                "unidentified"="#696969")

fun_oligo_tax.2 <- fun_oligo_tax

tree.plot <- ggtree(fun_oligo_phy_tree, ladderize = F,branch.length="none") %<+% fun_oligo_tax +
  geom_tippoint(aes(color = Phylum2)) +  scale_color_manual(values = fun.colors)+
  geom_tiplab(aes(color = Phylum2), align = T, size = 0, alpha = 0.2)+
  theme(legend.position = "left")

plot_data_tree<-tree.plot$data
plot_data_tree.t<- plot_data_tree[-c(10:17)]
fun_oligo_tax.2$label <- rownames(fun_oligo_tax.2)
plot_data_tree.t <- merge(plot_data_tree.t, fun_oligo_tax.2, by = "label", all = T)
plot_data_tree.t<- plot_data_tree.t %>% arrange(node)

tree.plot$data<- plot_data_tree.t

tree.plot

tax_order_fun.oligo <- data.frame(variable = as.character(fun_oligo_phy_tree$tip.label), 
                                   order = 1:length(fun_oligo_phy_tree$tip.label)) # reverse order


### Copio

library(ggtree)
bac.clean.ss.tree.sig.copio <- pop_taxa_wanted(bac.clean.ss.tree, sub.mydata.b.copio$OTU)
plot_tree(bac.clean.ss.tree.sig.copio)


bac_copio_phy_tree <- phy_tree(bac.clean.ss.tree.sig.copio)

bac_copio_tax <- tax_table(bac.clean.ss.tree.sig.copio)
bac_copio_tax <- data.frame(bac_copio_tax, stringsAsFactors = F)
bac_copio_tax$Phylum2 <- bac_copio_tax$Phylum
bac_copio_tax$Phylum2[is.na(bac_copio_tax$Phylum2)] <- "Unidentified"
bac_copio_tax$Class[is.na(bac_copio_tax$Class)] <- "Unidentified"
unique(bac_copio_tax$Phylum2)

bac_copio_tax[ rownames(bac_copio_tax)[bac_copio_tax$Class=="Alphaproteobacteria" ], ]$Phylum2 <- "Alphaproteobacteria"
#bac_copio_tax[ rownames(bac_copio_tax)[bac_copio_tax$Class=="Betaproteobacteria" ], ]$labels <- "Betaproteobacteria"
bac_copio_tax[ rownames(bac_copio_tax)[bac_copio_tax$Class=="Gammaproteobacteria" ], ]$Phylum2 <- "Gammaproteobacteria"
bac_copio_tax[ rownames(bac_copio_tax)[bac_copio_tax$Class=="Deltaproteobacteria" ], ]$Phylum2 <- "Deltaproteobacteria"
table(bac_copio_tax$Phylum2)


bac.colors.copio <- c("Alphaproteobacteria"="palegreen1","Gammaproteobacteria" = "palegreen2","Deltaproteobacteria"= "palegreen3", "Actinobacteria"= "indianred2",
                "Bacteroidetes"="steelblue1", "Acidobacteria"="lightsalmon4", "Chloroflexi" = "gold1",
                "Verrucomicrobia" = "orchid3", "Nitrospirae"= "palevioletred2", "Gemmatimonadetes"="peachpuff3",
                "Planctomycetes"= "khaki3","Actinobacteria" = "indianred2", "Unidentified" = "black", "Zixibacteria" = "gray50")

bac_copio_tax.2 <- bac_copio_tax

tree.plot <- ggtree(bac_copio_phy_tree, ladderize = F,branch.length="none") %<+% bac_copio_tax +
  geom_tippoint(aes(color = Phylum2)) +  scale_color_manual(values = bac.colors.copio)+
  geom_tiplab(aes(color = Phylum2), align = T, size = 0, alpha = 0.2)+
  theme(legend.position = "left")

plot_data_tree<-tree.plot$data
plot_data_tree.t<- plot_data_tree[-c(10:17)]
bac_copio_tax.2$label <- rownames(bac_copio_tax.2)
plot_data_tree.t <- merge(plot_data_tree.t, bac_copio_tax.2, by = "label", all = T)
plot_data_tree.t<- plot_data_tree.t %>% arrange(node)

tree.plot$data<- plot_data_tree.t

tree.plot

tax_order_bac.copio <- data.frame(variable = as.character(bac_copio_phy_tree$tip.label), 
                                  order = 1:length(bac_copio_phy_tree$tip.label)) # reverse order

## copio archaea tree
arch.clean.ss.tree.sig.copio <- pop_taxa_wanted(arch.clean.ss.tree, sub.mydata.a.copio$OTU)
plot_tree(arch.clean.ss.tree.sig.copio)


arch_copio_phy_tree <- phy_tree(arch.clean.ss.tree.sig.copio)

arch_copio_tax <- tax_table(arch.clean.ss.tree.sig.copio)
arch_copio_tax <- data.frame(arch_copio_tax, stringsAsFactors = F)
arch_copio_tax$Phylum2 <- arch_copio_tax$Phylum
unique(arch_copio_tax$Phylum2)

arch.colors.copio<- c("Euryarchaeota" = "#7B55E1","Crenarchaeota"= "#FF9999", "Thaumarchaeota"= "#CC6600",
                "Nanoarchaeaeota"="#999933", "Altiarchaeota"="gray50", "Diapherotrites" = "#3366FF")

arch_copio_tax.2 <- arch_copio_tax

tree.plot <- ggtree(arch_copio_phy_tree, ladderize = F,branch.length="none") %<+% arch_copio_tax +
  geom_tippoint(aes(color = Phylum2)) +  scale_color_manual(values = arch.colors.copio)+
  geom_tiplab(aes(color = Phylum2), align = T, size = 0, alpha = 0.2)+
  theme(legend.position = "left")

plot_data_tree<-tree.plot$data
plot_data_tree.t<- plot_data_tree[-c(10:17)]
arch_copio_tax.2$label <- rownames(arch_copio_tax.2)
plot_data_tree.t <- merge(plot_data_tree.t, arch_copio_tax.2, by = "label", all = T)
plot_data_tree.t<- plot_data_tree.t %>% arrange(node)

tree.plot$data<- plot_data_tree.t

tree.plot

tax_order_arch.copio <- data.frame(variable = as.character(arch_copio_phy_tree$tip.label), 
                                   order = 1:length(arch_copio_phy_tree$tip.label)) # reverse order


## copio fungi tree
fun.clean.ss.tree.sig.copio <- pop_taxa_wanted(fun.clean.ss, sub.mydata.f.copio$OTU)
plot_tree(fun.clean.ss.tree.sig.copio)


fun_copio_phy_tree <- phy_tree(fun.clean.ss.tree.sig.copio)

fun_copio_tax <- tax_table(fun.clean.ss.tree.sig.copio)
fun_copio_tax <- data.frame(fun_copio_tax, stringsAsFactors = F)
fun_copio_tax$Phylum2 <- fun_copio_tax$Phylum
fun_copio_tax$Phylum2[is.na(fun_copio_tax$Phylum2)] <- "unidentified"


unique(fun_copio_tax$Phylum2)

tax_its[ rownames(tax_its)[tax_its$labels=="Ascomycota" ], ]$cols <- "#003366"
tax_its[ rownames(tax_its)[tax_its$labels=="Basidiomycota" ], ]$cols <- "#C84248"
tax_its[ rownames(tax_its)[tax_its$labels=="Unassigned" ], ]$cols <- "#696969"
tax_its[ rownames(tax_its)[tax_its$labels=="Mortierellomycota" ], ]$cols <- "#EEB422"

fun.colors.copio<- c("Ascomycota" = "#003366","Chytridiomycota"= "#336600", "Mortierellomycota"= "#EEB422",
               "unidentified"="#696969", "Rozellomycota"= "#663399"  )

fun_copio_tax.2 <- fun_copio_tax

tree.plot <- ggtree(fun_copio_phy_tree, ladderize = F,branch.length="none") %<+% fun_copio_tax +
  geom_tippoint(aes(color = Phylum2)) +  scale_color_manual(values = fun.colors.copio)+
  geom_tiplab(aes(color = Phylum2), align = T, size = 0, alpha = 0.2)+
  theme(legend.position = "left")

plot_data_tree<-tree.plot$data
plot_data_tree.t<- plot_data_tree[-c(10:17)]
fun_copio_tax.2$label <- rownames(fun_copio_tax.2)
plot_data_tree.t <- merge(plot_data_tree.t, fun_copio_tax.2, by = "label", all = T)
plot_data_tree.t<- plot_data_tree.t %>% arrange(node)

tree.plot$data<- plot_data_tree.t

tree.plot

tax_order_fun.copio <- data.frame(variable = as.character(fun_copio_phy_tree$tip.label), 
                                  order = 1:length(fun_copio_phy_tree$tip.label)) # reverse order

