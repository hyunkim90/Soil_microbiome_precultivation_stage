## Distribution of functional taxa
## oligotrophic vs copiotrophic
otu.list <- rbind(phy.list, arch.list, fun.list)

bac.clean.ss
arch.clean.ss
fun.clean.ss

methanogen <- c('Methanobacterium', 'Methanosarcina', 'Methanocella','Methanosaeta')
syntroph <- c("Syntrophomonas", 'Syntrophobacter','Pelotomaculum','Clostridium','Syntrophus')
methanecycling <- c("Vogesella",
                    "Paenibacillus",
                    "Treponema",
                    "Cellulomonas",
                    "Candidatus_Solibacter",
                    "Pedomicrobium",
                    "Anaerolinea",
                    "Opitutus",
                    "Phenylobacterium",
                    "Dechloromonas")
methanotroph <- c("Crenothrix",
                  "Methylosinus","Methylobacter","Methylomonas", "Methylosarcina", "Methylomicrobium", "Methylococcus", "Methylocaldum", "Methylocella", "Methylocapsa", "Methylocystis")
iron_reducer <- c("Geobacter", "Anaeromyxobacter")
otu.list.methano <-subset(otu.list, Genus %in% methanogen)
otu.list.methanotroph <-subset(otu.list, Genus %in% methanotroph)
otu.list.methanecycling <-subset(otu.list, Genus %in% methanecycling)
otu.list.syntroph <-subset(otu.list, Genus %in% syntroph)
otu.list.iron_reducer <-subset(otu.list, Genus %in% iron_reducer)


## profiles of oligo and copio OTUs
list.bac.oligo.copio <- read.table("bac_RA.tsv", sep= '\t', header = T)
list.arch.oligo.copio <- read.table("arch_RA.tsv", sep= '\t', header = T)
list.fun.oligo.copio <- read.table("Fun_RA.tsv", sep= '\t', header = T)
head(list.arch.oligo.copio)

#Methanogen
list.arch.methano <- subset(list.arch.oligo.copio, Genus %in% methanogen)
head(list.arch.methano)

poor.arch.methano <- list.arch.methano %>% select(OTU_id, Category, CJ1, CJ2, YS1, YS2, CJ1.18, CJ2.18)

rich.arch.methano <- list.arch.methano %>% select(OTU_id, Category, MY1, MY2, NJ1, NJ2)

#Methanotroph
list.bac.methanotroph <- subset(list.bac.oligo.copio, Genus %in% methanotroph)
head(list.bac.methanotroph)

poor.bac.methanotroph <- list.bac.methanotroph %>% select(OTU_id, Category, CJ1, CJ2, YS1, YS2, CJ1.18, CJ2.18)

rich.bac.methanotroph <- list.bac.methanotroph %>% select(OTU_id, Category, MY1, MY2, NJ1, NJ2)


#Syntroph
list.bac.syntroph <- subset(list.bac.oligo.copio, Genus %in% syntroph)
head(list.bac.syntroph)

poor.bac.syntroph <- list.bac.syntroph %>% select(OTU_id, Category, CJ1, CJ2, YS1, YS2, CJ1.18, CJ2.18)

rich.bac.syntroph <- list.bac.syntroph %>% select(OTU_id, Category, MY1, MY2, NJ1, NJ2)

#Methane cycling
list.bac.methanecycling <- subset(list.bac.oligo.copio, Genus %in% methanecycling)
head(list.bac.methanecycling)

poor.bac.methanecycling <- list.bac.methanecycling %>% select(OTU_id, Category, CJ1, CJ2, YS1, YS2, CJ1.18, CJ2.18)

rich.bac.methanecycling <- list.bac.methanecycling %>% select(OTU_id, Category, MY1, MY2, NJ1, NJ2)


#Methane cycling
list.bac.iron_reducer <- subset(list.bac.oligo.copio, Genus %in% iron_reducer)
head(list.bac.iron_reducer)

poor.bac.iron_reducer <- list.bac.iron_reducer %>% select(OTU_id, Category, CJ1, CJ2, YS1, YS2, CJ1.18, CJ2.18)

rich.bac.iron_reducer <- list.bac.iron_reducer %>% select(OTU_id, Category, MY1, MY2, NJ1, NJ2)



## Methanogen
##Mean RA of each OTU in nutrient poor and rich
## Nutrient-poor
poor.arch.methano
rownames(poor.arch.methano) <- poor.arch.methano$OTU_id

##Non-fer
poor.arch.methano.nonfer <- poor.arch.methano%>% select(CJ1, YS1, CJ1.18)
poor.arch.methano.nonfer$Mean_RA <- rowMeans(poor.arch.methano.nonfer)

poor.arch.methano.nonfer <- poor.arch.methano.nonfer %>% select(Mean_RA)
poor.arch.methano.oligo.copio <- poor.arch.methano %>% select(Category)

poor.arch.methano.nonfer<-merge(poor.arch.methano.nonfer, poor.arch.methano.oligo.copio, by = 'row.names')
poor.arch.methano.nonfer$Category2 <- "Non-fer"
poor.arch.methano.nonfer$Category3 <- "Poor"

#Fer
poor.arch.methano.fer <- poor.arch.methano%>% select(CJ2, YS2, CJ2.18)
poor.arch.methano.fer$Mean_RA <- rowMeans(poor.arch.methano.fer)

poor.arch.methano.fer <- poor.arch.methano.fer %>% select(Mean_RA)
poor.arch.methano.oligo.copio <- poor.arch.methano %>% select(Category)

poor.arch.methano.fer<-merge(poor.arch.methano.fer, poor.arch.methano.oligo.copio, by = 'row.names')
poor.arch.methano.fer$Category2 <- "Fer"
poor.arch.methano.fer$Category3 <- "Poor"


## Nutrient-rich
rich.arch.methano
rownames(rich.arch.methano) <- rich.arch.methano$OTU_id

##Non-fer
rich.arch.methano.nonfer <- rich.arch.methano%>% select(MY2, NJ2)
rich.arch.methano.nonfer$Mean_RA <- rowMeans(rich.arch.methano.nonfer)

rich.arch.methano.nonfer <- rich.arch.methano.nonfer %>% select(Mean_RA)
rich.arch.methano.oligo.copio <- rich.arch.methano %>% select(Category)

rich.arch.methano.nonfer<-merge(rich.arch.methano.nonfer, rich.arch.methano.oligo.copio, by = 'row.names')
rich.arch.methano.nonfer$Category2 <- "Non-fer"
rich.arch.methano.nonfer$Category3 <- "Rich"


#Fer
rich.arch.methano.fer <- rich.arch.methano%>% select(MY1, NJ1)
rich.arch.methano.fer$Mean_RA <- rowMeans(rich.arch.methano.fer)

rich.arch.methano.fer <- rich.arch.methano.fer %>% select(Mean_RA)
rich.arch.methano.oligo.copio <- rich.arch.methano %>% select(Category)

rich.arch.methano.fer<-merge(rich.arch.methano.fer, rich.arch.methano.oligo.copio, by = 'row.names')
rich.arch.methano.fer$Category2 <- "Fer"
rich.arch.methano.fer$Category3 <- "Rich"

arch.methanogen.trophic <- rbind(poor.arch.methano.nonfer, poor.arch.methano.fer, rich.arch.methano.nonfer, rich.arch.methano.fer)

x <- subset(arch.methanogen.trophic, Category3=='Poor')$Mean_RA
y <- subset(arch.methanogen.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)#0.01504138
sum(x) #0.5349507
length(x) #42
length(x[x>0]) #28
mean(y)#0.005541192
sum(y) #0.2327301
length(y) #42
length(y[y>0]) #15

p <- ggplot(data =arch.methanogen.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


##Syntroph
poor.bac.syntroph
rownames(poor.bac.syntroph) <- poor.bac.syntroph$OTU_id

##Non-fer
poor.bac.syntroph.nonfer <- poor.bac.syntroph%>% select(CJ1, YS1, CJ1.18)
poor.bac.syntroph.nonfer$Mean_RA <- rowMeans(poor.bac.syntroph.nonfer)

poor.bac.syntroph.nonfer <- poor.bac.syntroph.nonfer %>% select(Mean_RA)
poor.bac.syntroph.oligo.copio <- poor.bac.syntroph %>% select(Category)

poor.bac.syntroph.nonfer<-merge(poor.bac.syntroph.nonfer, poor.bac.syntroph.oligo.copio, by = 'row.names')
poor.bac.syntroph.nonfer$Category2 <- "Non-fer"
poor.bac.syntroph.nonfer$Category3 <- "Poor"

#Fer
poor.bac.syntroph.fer <- poor.bac.syntroph%>% select(CJ2, YS2, CJ2.18)
poor.bac.syntroph.fer$Mean_RA <- rowMeans(poor.bac.syntroph.fer)

poor.bac.syntroph.fer <- poor.bac.syntroph.fer %>% select(Mean_RA)
poor.bac.syntroph.oligo.copio <- poor.bac.syntroph %>% select(Category)

poor.bac.syntroph.fer<-merge(poor.bac.syntroph.fer, poor.bac.syntroph.oligo.copio, by = 'row.names')
poor.bac.syntroph.fer$Category2 <- "Fer"
poor.bac.syntroph.fer$Category3 <- "Poor"


## Nutrient-rich
rich.bac.syntroph
rownames(rich.bac.syntroph) <- rich.bac.syntroph$OTU_id

##Non-fer
rich.bac.syntroph.nonfer <- rich.bac.syntroph%>% select(MY2, NJ2)
rich.bac.syntroph.nonfer$Mean_RA <- rowMeans(rich.bac.syntroph.nonfer)

rich.bac.syntroph.nonfer <- rich.bac.syntroph.nonfer %>% select(Mean_RA)
rich.bac.syntroph.oligo.copio <- rich.bac.syntroph %>% select(Category)

rich.bac.syntroph.nonfer<-merge(rich.bac.syntroph.nonfer, rich.bac.syntroph.oligo.copio, by = 'row.names')
rich.bac.syntroph.nonfer$Category2 <- "Non-fer"
rich.bac.syntroph.nonfer$Category3 <- "Rich"


#Fer
rich.bac.syntroph.fer <- rich.bac.syntroph%>% select(MY1, NJ1)
rich.bac.syntroph.fer$Mean_RA <- rowMeans(rich.bac.syntroph.fer)

rich.bac.syntroph.fer <- rich.bac.syntroph.fer %>% select(Mean_RA)
rich.bac.syntroph.oligo.copio <- rich.bac.syntroph %>% select(Category)

rich.bac.syntroph.fer<-merge(rich.bac.syntroph.fer, rich.bac.syntroph.oligo.copio, by = 'row.names')
rich.bac.syntroph.fer$Category2 <- "Fer"
rich.bac.syntroph.fer$Category3 <- "Rich"

bac.syntrophgen.trophic <- rbind(poor.bac.syntroph.nonfer, poor.bac.syntroph.fer, rich.bac.syntroph.nonfer, rich.bac.syntroph.fer)

x <- subset(bac.syntrophgen.trophic, Category3=='Poor')$Mean_RA
y <- subset(bac.syntrophgen.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #12
mean(y)
length(y) #12
sum(x)#0.002639022
sum(y)#0.006893105
length(x[x>0]) #6
length(y[y>0]) #6
p <- ggplot(data =bac.syntrophgen.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


##Methane_cycling
poor.bac.methanecycling
rownames(poor.bac.methanecycling) <- poor.bac.methanecycling$OTU_id

##Non-fer
poor.bac.methanecycling.nonfer <- poor.bac.methanecycling%>% select(CJ1, YS1, CJ1.18)
poor.bac.methanecycling.nonfer$Mean_RA <- rowMeans(poor.bac.methanecycling.nonfer)

poor.bac.methanecycling.nonfer <- poor.bac.methanecycling.nonfer %>% select(Mean_RA)
poor.bac.methanecycling.oligo.copio <- poor.bac.methanecycling %>% select(Category)

poor.bac.methanecycling.nonfer<-merge(poor.bac.methanecycling.nonfer, poor.bac.methanecycling.oligo.copio, by = 'row.names')
poor.bac.methanecycling.nonfer$Category2 <- "Non-fer"
poor.bac.methanecycling.nonfer$Category3 <- "Poor"

#Fer
poor.bac.methanecycling.fer <- poor.bac.methanecycling%>% select(CJ2, YS2, CJ2.18)
poor.bac.methanecycling.fer$Mean_RA <- rowMeans(poor.bac.methanecycling.fer)

poor.bac.methanecycling.fer <- poor.bac.methanecycling.fer %>% select(Mean_RA)
poor.bac.methanecycling.oligo.copio <- poor.bac.methanecycling %>% select(Category)

poor.bac.methanecycling.fer<-merge(poor.bac.methanecycling.fer, poor.bac.methanecycling.oligo.copio, by = 'row.names')
poor.bac.methanecycling.fer$Category2 <- "Fer"
poor.bac.methanecycling.fer$Category3 <- "Poor"


## Nutrient-rich
rich.bac.methanecycling
rownames(rich.bac.methanecycling) <- rich.bac.methanecycling$OTU_id

##Non-fer
rich.bac.methanecycling.nonfer <- rich.bac.methanecycling%>% select(MY2, NJ2)
rich.bac.methanecycling.nonfer$Mean_RA <- rowMeans(rich.bac.methanecycling.nonfer)

rich.bac.methanecycling.nonfer <- rich.bac.methanecycling.nonfer %>% select(Mean_RA)
rich.bac.methanecycling.oligo.copio <- rich.bac.methanecycling %>% select(Category)

rich.bac.methanecycling.nonfer<-merge(rich.bac.methanecycling.nonfer, rich.bac.methanecycling.oligo.copio, by = 'row.names')
rich.bac.methanecycling.nonfer$Category2 <- "Non-fer"
rich.bac.methanecycling.nonfer$Category3 <- "Rich"


#Fer
rich.bac.methanecycling.fer <- rich.bac.methanecycling%>% select(MY1, NJ1)
rich.bac.methanecycling.fer$Mean_RA <- rowMeans(rich.bac.methanecycling.fer)

rich.bac.methanecycling.fer <- rich.bac.methanecycling.fer %>% select(Mean_RA)
rich.bac.methanecycling.oligo.copio <- rich.bac.methanecycling %>% select(Category)

rich.bac.methanecycling.fer<-merge(rich.bac.methanecycling.fer, rich.bac.methanecycling.oligo.copio, by = 'row.names')
rich.bac.methanecycling.fer$Category2 <- "Fer"
rich.bac.methanecycling.fer$Category3 <- "Rich"

bac.methanecycling.trophic <- rbind(poor.bac.methanecycling.nonfer, poor.bac.methanecycling.fer, rich.bac.methanecycling.nonfer, rich.bac.methanecycling.fer)

x <- subset(bac.methanecycling.trophic, Category3=='Poor')$Mean_RA
y <- subset(bac.methanecycling.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #32
mean(y)
length(y) #32

p <- ggplot(data =bac.methanecycling.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter',aes(shape=Category), alpha=.5)+
    #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")

p


##Methanotroph
poor.bac.methanotroph
rownames(poor.bac.methanotroph) <- poor.bac.methanotroph$OTU_id

##Non-fer
poor.bac.methanotroph.nonfer <- poor.bac.methanotroph%>% select(CJ1, YS1, CJ1.18)
poor.bac.methanotroph.nonfer$Mean_RA <- rowMeans(poor.bac.methanotroph.nonfer)

poor.bac.methanotroph.nonfer <- poor.bac.methanotroph.nonfer %>% select(Mean_RA)
poor.bac.methanotroph.oligo.copio <- poor.bac.methanotroph %>% select(Category)

poor.bac.methanotroph.nonfer<-merge(poor.bac.methanotroph.nonfer, poor.bac.methanotroph.oligo.copio, by = 'row.names')
poor.bac.methanotroph.nonfer$Category2 <- "Non-fer"
poor.bac.methanotroph.nonfer$Category3 <- "Poor"

#Fer
poor.bac.methanotroph.fer <- poor.bac.methanotroph%>% select(CJ2, YS2, CJ2.18)
poor.bac.methanotroph.fer$Mean_RA <- rowMeans(poor.bac.methanotroph.fer)

poor.bac.methanotroph.fer <- poor.bac.methanotroph.fer %>% select(Mean_RA)
poor.bac.methanotroph.oligo.copio <- poor.bac.methanotroph %>% select(Category)

poor.bac.methanotroph.fer<-merge(poor.bac.methanotroph.fer, poor.bac.methanotroph.oligo.copio, by = 'row.names')
poor.bac.methanotroph.fer$Category2 <- "Fer"
poor.bac.methanotroph.fer$Category3 <- "Poor"


## Nutrient-rich
rich.bac.methanotroph
rownames(rich.bac.methanotroph) <- rich.bac.methanotroph$OTU_id

##Non-fer
rich.bac.methanotroph.nonfer <- rich.bac.methanotroph%>% select(MY2, NJ2)
rich.bac.methanotroph.nonfer$Mean_RA <- rowMeans(rich.bac.methanotroph.nonfer)

rich.bac.methanotroph.nonfer <- rich.bac.methanotroph.nonfer %>% select(Mean_RA)
rich.bac.methanotroph.oligo.copio <- rich.bac.methanotroph %>% select(Category)

rich.bac.methanotroph.nonfer<-merge(rich.bac.methanotroph.nonfer, rich.bac.methanotroph.oligo.copio, by = 'row.names')
rich.bac.methanotroph.nonfer$Category2 <- "Non-fer"
rich.bac.methanotroph.nonfer$Category3 <- "Rich"


#Fer
rich.bac.methanotroph.fer <- rich.bac.methanotroph%>% select(MY1, NJ1)
rich.bac.methanotroph.fer$Mean_RA <- rowMeans(rich.bac.methanotroph.fer)

rich.bac.methanotroph.fer <- rich.bac.methanotroph.fer %>% select(Mean_RA)
rich.bac.methanotroph.oligo.copio <- rich.bac.methanotroph %>% select(Category)

rich.bac.methanotroph.fer<-merge(rich.bac.methanotroph.fer, rich.bac.methanotroph.oligo.copio, by = 'row.names')
rich.bac.methanotroph.fer$Category2 <- "Fer"
rich.bac.methanotroph.fer$Category3 <- "Rich"

bac.methanotroph.trophic <- rbind(poor.bac.methanotroph.nonfer, poor.bac.methanotroph.fer, rich.bac.methanotroph.nonfer, rich.bac.methanotroph.fer)

x <- subset(bac.methanotroph.trophic, Category3=='Poor')$Mean_RA
y <- subset(bac.methanotroph.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #8
mean(y)
length(y) #8

p <- ggplot(data =bac.methanotroph.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


##Non-acetogenic


poor.bac.iron_reducer
rownames(poor.bac.iron_reducer) <- poor.bac.iron_reducer$OTU_id

##Non-fer
poor.bac.iron_reducer.nonfer <- poor.bac.iron_reducer%>% select(CJ1, YS1, CJ1.18)
poor.bac.iron_reducer.nonfer$Mean_RA <- rowMeans(poor.bac.iron_reducer.nonfer)

poor.bac.iron_reducer.nonfer <- poor.bac.iron_reducer.nonfer %>% select(Mean_RA)
poor.bac.iron_reducer.oligo.copio <- poor.bac.iron_reducer %>% select(Category)

poor.bac.iron_reducer.nonfer<-merge(poor.bac.iron_reducer.nonfer, poor.bac.iron_reducer.oligo.copio, by = 'row.names')
poor.bac.iron_reducer.nonfer$Category2 <- "Non-fer"
poor.bac.iron_reducer.nonfer$Category3 <- "Poor"

#Fer
poor.bac.iron_reducer.fer <- poor.bac.iron_reducer%>% select(CJ2, YS2, CJ2.18)
poor.bac.iron_reducer.fer$Mean_RA <- rowMeans(poor.bac.iron_reducer.fer)

poor.bac.iron_reducer.fer <- poor.bac.iron_reducer.fer %>% select(Mean_RA)
poor.bac.iron_reducer.oligo.copio <- poor.bac.iron_reducer %>% select(Category)

poor.bac.iron_reducer.fer<-merge(poor.bac.iron_reducer.fer, poor.bac.iron_reducer.oligo.copio, by = 'row.names')
poor.bac.iron_reducer.fer$Category2 <- "Fer"
poor.bac.iron_reducer.fer$Category3 <- "Poor"


## Nutrient-rich
rich.bac.iron_reducer
rownames(rich.bac.iron_reducer) <- rich.bac.iron_reducer$OTU_id

##Non-fer
rich.bac.iron_reducer.nonfer <- rich.bac.iron_reducer%>% select(MY2, NJ2)
rich.bac.iron_reducer.nonfer$Mean_RA <- rowMeans(rich.bac.iron_reducer.nonfer)

rich.bac.iron_reducer.nonfer <- rich.bac.iron_reducer.nonfer %>% select(Mean_RA)
rich.bac.iron_reducer.oligo.copio <- rich.bac.iron_reducer %>% select(Category)

rich.bac.iron_reducer.nonfer<-merge(rich.bac.iron_reducer.nonfer, rich.bac.iron_reducer.oligo.copio, by = 'row.names')
rich.bac.iron_reducer.nonfer$Category2 <- "Non-fer"
rich.bac.iron_reducer.nonfer$Category3 <- "Rich"


#Fer
rich.bac.iron_reducer.fer <- rich.bac.iron_reducer%>% select(MY1, NJ1)
rich.bac.iron_reducer.fer$Mean_RA <- rowMeans(rich.bac.iron_reducer.fer)

rich.bac.iron_reducer.fer <- rich.bac.iron_reducer.fer %>% select(Mean_RA)
rich.bac.iron_reducer.oligo.copio <- rich.bac.iron_reducer %>% select(Category)

rich.bac.iron_reducer.fer<-merge(rich.bac.iron_reducer.fer, rich.bac.iron_reducer.oligo.copio, by = 'row.names')
rich.bac.iron_reducer.fer$Category2 <- "Fer"
rich.bac.iron_reducer.fer$Category3 <- "Rich"

bac.iron_reducer.trophic <- rbind(poor.bac.iron_reducer.nonfer, poor.bac.iron_reducer.fer, rich.bac.iron_reducer.nonfer, rich.bac.iron_reducer.fer)

x <- subset(bac.iron_reducer.trophic, Category3=='Poor')$Mean_RA
y <- subset(bac.iron_reducer.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #52
mean(y)
length(y) #52
sum(x) #0.01877899
sum(y) #0.009767985
length(x[x>0]) #34
length(y[y>0]) #19
p <- ggplot(data =bac.iron_reducer.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


### each nutrient condition
###In nutrient-poor
arch.methanogen.trophic$Category2 <- factor(arch.methanogen.trophic$Category2, levels = c("Non-fer", "Fer"))
arch.methanogen.trophic.poor <- subset(arch.methanogen.trophic, Category3 == "Poor")
arch.methanogen.trophic.poor.copio <- subset(arch.methanogen.trophic.poor, Category == "Copiotroph" & Mean_RA>0)
x <- subset(arch.methanogen.trophic.poor.copio, Category2=='Non-fer')$Mean_RA
y <- subset(arch.methanogen.trophic.poor.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value 
mean(x)
length(x) #3
mean(y)
length(y) #3
sum(x) #0.04619469
sum(y) #0.017628
length(x[x>0]) #34
length(y[y>0]) #19

p <- ggplot(data =arch.methanogen.trophic.poor.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


arch.methanogen.trophic.poor <- subset(arch.methanogen.trophic, Category3 == "Poor")
arch.methanogen.trophic.poor.oligo <- subset(arch.methanogen.trophic.poor, Category == "Oligotroph" & Mean_RA >0)
x <- subset(arch.methanogen.trophic.poor.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(arch.methanogen.trophic.poor.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #12
mean(y)
length(y) #10
sum(x) #0.1977114
sum(y) #0.2734166
length(x[x>0]) #12
length(y[y>0]) #10


p <- ggplot(data =arch.methanogen.trophic.poor.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


###In nutrient-rich
arch.methanogen.trophic.rich <- subset(arch.methanogen.trophic, Category3 == "Rich")
arch.methanogen.trophic.rich.copio <- subset(arch.methanogen.trophic.rich, Category == "Copiotroph" & Mean_RA >0)
x <- subset(arch.methanogen.trophic.rich.copio, Category2=='Non-fer')$Mean_RA
y <- subset(arch.methanogen.trophic.rich.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #3
mean(y)
length(y) #3
sum(x) #0.04258724
sum(y) #0.06538185
length(x[x>0]) #3
length(y[y>0]) #3


p <- ggplot(data =arch.methanogen.trophic.rich.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


arch.methanogen.trophic.rich <- subset(arch.methanogen.trophic, Category3 == "Rich")
arch.methanogen.trophic.rich.oligo <- subset(arch.methanogen.trophic.rich, Category == "Oligotroph" & Mean_RA>0)
x <- subset(arch.methanogen.trophic.rich.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(arch.methanogen.trophic.rich.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #4
mean(y)
length(y) #5
sum(x) #0.0604136
sum(y) #0.06434737
length(x[x>0]) #4
length(y[y>0]) #5

p <- ggplot(data =arch.methanogen.trophic.rich.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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

##Methanotroph
bac.methanotroph.trophic$Category2 <- factor(bac.methanotroph.trophic$Category2, levels = c("Non-fer", "Fer"))
bac.methanotroph.trophic.poor <- subset(bac.methanotroph.trophic, Category3 == "Poor")
bac.methanotroph.trophic.poor.copio <- subset(bac.methanotroph.trophic.poor, Category == "Copiotroph")
x <- subset(bac.methanotroph.trophic.poor.copio, Category2=='Non-fer')$Mean_RA
y <- subset(bac.methanotroph.trophic.poor.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value 
mean(x)
length(x) #2
mean(y)
length(y) #2

p <- ggplot(data =bac.methanotroph.trophic.poor.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


bac.methanotroph.trophic.poor <- subset(bac.methanotroph.trophic, Category3 == "Poor")
bac.methanotroph.trophic.poor.oligo <- subset(bac.methanotroph.trophic.poor, Category == "Oligotroph")
x <- subset(bac.methanotroph.trophic.poor.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(bac.methanotroph.trophic.poor.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #13
mean(y)
length(y) #13

p <- ggplot(data =bac.methanotroph.trophic.poor.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


###In nutrient-rich
bac.methanotroph.trophic.rich <- subset(bac.methanotroph.trophic, Category3 == "Rich")
bac.methanotroph.trophic.rich.copio <- subset(bac.methanotroph.trophic.rich, Category == "Copiotroph")
x <- subset(bac.methanotroph.trophic.rich.copio, Category2=='Non-fer')$Mean_RA
y <- subset(bac.methanotroph.trophic.rich.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x)
mean(y)
length(y) #8

p <- ggplot(data =bac.methanotroph.trophic.rich.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


bac.methanotroph.trophic.rich <- subset(bac.methanotroph.trophic, Category3 == "Rich")
bac.methanotroph.trophic.rich.oligo <- subset(bac.methanotroph.trophic.rich, Category == "Oligotroph")
x <- subset(bac.methanotroph.trophic.rich.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(bac.methanotroph.trophic.rich.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #13
mean(y)
length(y)

p <- ggplot(data =bac.methanotroph.trophic.rich.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


##Syntroph
bac.syntroph.trophic$Category2 <- factor(bac.syntroph.trophic$Category2, levels = c("Non-fer", "Fer"))
bac.syntroph.trophic.poor <- subset(bac.syntroph.trophic, Category3 == "Poor")
bac.syntroph.trophic.poor.copio <- subset(bac.syntroph.trophic.poor, Category == "Copiotroph")
x <- subset(bac.syntroph.trophic.poor.copio, Category2=='Non-fer')$Mean_RA
y <- subset(bac.syntroph.trophic.poor.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value 
mean(x)
length(x) #2
mean(y)
length(y) #2

p <- ggplot(data =bac.syntroph.trophic.poor.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


bac.syntroph.trophic.poor <- subset(bac.syntroph.trophic, Category3 == "Poor")
bac.syntroph.trophic.poor.oligo <- subset(bac.syntroph.trophic.poor, Category == "Oligotroph")
x <- subset(bac.syntroph.trophic.poor.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(bac.syntroph.trophic.poor.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #18
mean(y)
length(y) #12

p <- ggplot(data =bac.syntroph.trophic.poor.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


###In nutrient-rich
bac.syntroph.trophic.rich <- subset(bac.syntroph.trophic, Category3 == "Rich")
bac.syntroph.trophic.rich.copio <- subset(bac.syntroph.trophic.rich, Category == "Copiotroph" & Mean_RA >0)
x <- subset(bac.syntroph.trophic.rich.copio, Category2=='Non-fer')$Mean_RA
y <- subset(bac.syntroph.trophic.rich.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #4
mean(y)
length(y) #4

p <- ggplot(data =bac.syntroph.trophic.rich.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


bac.syntroph.trophic.rich <- subset(bac.syntroph.trophic, Category3 == "Rich")
bac.syntroph.trophic.rich.oligo <- subset(bac.syntroph.trophic.rich, Category == "Oligotroph" & Mean_RA >0)
x <- subset(bac.syntroph.trophic.rich.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(bac.syntroph.trophic.rich.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6

p <- ggplot(data =bac.syntroph.trophic.rich.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


##Syntroph
bac.syntroph.trophic$Category2 <- factor(bac.syntroph.trophic$Category2, levels = c("Non-fer", "Fer"))
bac.syntroph.trophic.poor <- subset(bac.syntroph.trophic, Category3 == "Poor")
bac.syntroph.trophic.poor.copio <- subset(bac.syntroph.trophic.poor, Category == "Copiotroph" & Mean_RA >0)
x <- subset(bac.syntroph.trophic.poor.copio, Category2=='Non-fer')$Mean_RA
y <- subset(bac.syntroph.trophic.poor.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value 
mean(x)
length(x) #2
mean(y)
length(y) #2

p <- ggplot(data =bac.syntroph.trophic.poor.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


bac.syntroph.trophic.poor <- subset(bac.syntroph.trophic, Category3 == "Poor")
bac.syntroph.trophic.poor.oligo <- subset(bac.syntroph.trophic.poor, Category == "Oligotroph" & Mean_RA >0)
x <- subset(bac.syntroph.trophic.poor.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(bac.syntroph.trophic.poor.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #18
mean(y)
length(y) #12

p <- ggplot(data =bac.syntroph.trophic.poor.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


###In nutrient-rich
bac.syntroph.trophic.rich <- subset(bac.syntroph.trophic, Category3 == "Rich")
bac.syntroph.trophic.rich.copio <- subset(bac.syntroph.trophic.rich, Category == "Copiotroph" & Mean_RA >0)
x <- subset(bac.syntroph.trophic.rich.copio, Category2=='Non-fer')$Mean_RA
y <- subset(bac.syntroph.trophic.rich.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #4
mean(y)
length(y) #4

p <- ggplot(data =bac.syntroph.trophic.rich.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


bac.syntroph.trophic.rich <- subset(bac.syntroph.trophic, Category3 == "Rich")
bac.syntroph.trophic.rich.oligo <- subset(bac.syntroph.trophic.rich, Category == "Oligotroph" & Mean_RA >0)
x <- subset(bac.syntroph.trophic.rich.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(bac.syntroph.trophic.rich.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6

p <- ggplot(data =bac.syntroph.trophic.rich.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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

##Natural variation
arch.methanogen.trophic.nonfer <- subset(arch.methanogen.trophic, Category2 == "Non-fer")
arch.methanogen.trophic.nonfer.oligo <- subset(arch.methanogen.trophic.nonfer, Category == "Oligotroph")
x <- subset(arch.methanogen.trophic.nonfer.oligo, Category3=='Poor')$Mean_RA
y <- subset(arch.methanogen.trophic.nonfer.oligo, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
mean(y)

p <- ggplot(data =arch.methanogen.trophic.nonfer.oligo, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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

arch.methanogen.trophic.nonfer <- subset(arch.methanogen.trophic, Category2 == "Non-fer")
arch.methanogen.trophic.nonfer.copio <- subset(arch.methanogen.trophic.nonfer, Category == "Copiotroph")
x <- subset(arch.methanogen.trophic.nonfer.copio, Category3=='Poor')$Mean_RA
y <- subset(arch.methanogen.trophic.nonfer.copio, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
mean(y)

p <- ggplot(data =arch.methanogen.trophic.nonfer.copio, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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



arch.methanogen.trophic.fer <- subset(arch.methanogen.trophic, Category2 == "Fer")
arch.methanogen.trophic.fer.oligo <- subset(arch.methanogen.trophic.fer, Category == "Oligotroph")
x <- subset(arch.methanogen.trophic.fer.oligo, Category3=='Poor')$Mean_RA
y <- subset(arch.methanogen.trophic.fer.oligo, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
mean(y)

p <- ggplot(data =arch.methanogen.trophic.fer.oligo, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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

arch.methanogen.trophic.fer <- subset(arch.methanogen.trophic, Category2 == "Fer")
arch.methanogen.trophic.fer.copio <- subset(arch.methanogen.trophic.fer, Category == "Copiotroph")
x <- subset(arch.methanogen.trophic.fer.copio, Category3=='Poor')$Mean_RA
y <- subset(arch.methanogen.trophic.fer.copio, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
mean(y)

p <- ggplot(data =arch.methanogen.trophic.fer.copio, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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

### Fungal saprotroph
saprotroph <- c("Guehomyces", "Solicoccozyma", "Schizothecium", "Papiliotrema")
otu.list.saprotroph <-subset(otu.list, Genus %in% saprotroph)


list.fun.saprotroph <- subset(list.fun.oligo.copio, Genus %in% saprotroph)
head(list.fun.saprotroph)

poor.fun.saprotroph <- list.fun.saprotroph %>% select(OTU_id, Category, CJ1, CJ2, YS1, YS2, CJ1.18, CJ2.18)

rich.fun.saprotroph <- list.fun.saprotroph %>% select(OTU_id, Category, MY1, MY2, NJ1, NJ2)

##Syntroph
poor.fun.saprotroph
rownames(poor.fun.saprotroph) <- poor.fun.saprotroph$OTU_id

##Non-fer
poor.fun.saprotroph.nonfer <- poor.fun.saprotroph%>% select(CJ1, YS1, CJ1.18)
poor.fun.saprotroph.nonfer$Mean_RA <- rowMeans(poor.fun.saprotroph.nonfer)

poor.fun.saprotroph.nonfer <- poor.fun.saprotroph.nonfer %>% select(Mean_RA)
poor.fun.saprotroph.oligo.copio <- poor.fun.saprotroph %>% select(Category)

poor.fun.saprotroph.nonfer<-merge(poor.fun.saprotroph.nonfer, poor.fun.saprotroph.oligo.copio, by = 'row.names')
poor.fun.saprotroph.nonfer$Category2 <- "Non-fer"
poor.fun.saprotroph.nonfer$Category3 <- "Poor"

#Fer
poor.fun.saprotroph.fer <- poor.fun.saprotroph%>% select(CJ2, YS2, CJ2.18)
poor.fun.saprotroph.fer$Mean_RA <- rowMeans(poor.fun.saprotroph.fer)

poor.fun.saprotroph.fer <- poor.fun.saprotroph.fer %>% select(Mean_RA)
poor.fun.saprotroph.oligo.copio <- poor.fun.saprotroph %>% select(Category)

poor.fun.saprotroph.fer<-merge(poor.fun.saprotroph.fer, poor.fun.saprotroph.oligo.copio, by = 'row.names')
poor.fun.saprotroph.fer$Category2 <- "Fer"
poor.fun.saprotroph.fer$Category3 <- "Poor"


## Nutrient-rich
rich.fun.saprotroph
rownames(rich.fun.saprotroph) <- rich.fun.saprotroph$OTU_id

##Non-fer
rich.fun.saprotroph.nonfer <- rich.fun.saprotroph%>% select(MY2, NJ2)
rich.fun.saprotroph.nonfer$Mean_RA <- rowMeans(rich.fun.saprotroph.nonfer)

rich.fun.saprotroph.nonfer <- rich.fun.saprotroph.nonfer %>% select(Mean_RA)
rich.fun.saprotroph.oligo.copio <- rich.fun.saprotroph %>% select(Category)

rich.fun.saprotroph.nonfer<-merge(rich.fun.saprotroph.nonfer, rich.fun.saprotroph.oligo.copio, by = 'row.names')
rich.fun.saprotroph.nonfer$Category2 <- "Non-fer"
rich.fun.saprotroph.nonfer$Category3 <- "Rich"


#Fer
rich.fun.saprotroph.fer <- rich.fun.saprotroph%>% select(MY1, NJ1)
rich.fun.saprotroph.fer$Mean_RA <- rowMeans(rich.fun.saprotroph.fer)

rich.fun.saprotroph.fer <- rich.fun.saprotroph.fer %>% select(Mean_RA)
rich.fun.saprotroph.oligo.copio <- rich.fun.saprotroph %>% select(Category)

rich.fun.saprotroph.fer<-merge(rich.fun.saprotroph.fer, rich.fun.saprotroph.oligo.copio, by = 'row.names')
rich.fun.saprotroph.fer$Category2 <- "Fer"
rich.fun.saprotroph.fer$Category3 <- "Rich"

fun.saprotroph.trophic <- rbind(poor.fun.saprotroph.nonfer, poor.fun.saprotroph.fer, rich.fun.saprotroph.nonfer, rich.fun.saprotroph.fer)

x <- subset(fun.saprotroph.trophic, Category3=='Poor')$Mean_RA
y <- subset(fun.saprotroph.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #26
mean(y)
length(y) #26
sum(x) #0.3748284
sum(y) #0.09987724
length(x[x>0]) #23
length(y[y>0]) #17

p <- ggplot(data =fun.saprotroph.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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



##Saprotroph
fun.saprotroph.trophic$Category2 <- factor(fun.saprotroph.trophic$Category2, levels = c("Non-fer", "Fer"))
fun.saprotroph.trophic.poor <- subset(fun.saprotroph.trophic, Category3 == "Poor")
fun.saprotroph.trophic.poor.copio <- subset(fun.saprotroph.trophic.poor, Category == "Copiotroph")
x <- subset(fun.saprotroph.trophic.poor.copio, Category2=='Non-fer')$Mean_RA
y <- subset(fun.saprotroph.trophic.poor.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value 
mean(x)
length(x) #6
mean(y)
length(y) #6
sum(x) #0.01286591
sum(y) #0.0107054
length(x[x>0]) #6
length(y[y>0]) #6

p <- ggplot(data =fun.saprotroph.trophic.poor.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


fun.saprotroph.trophic.poor <- subset(fun.saprotroph.trophic, Category3 == "Poor")
fun.saprotroph.trophic.poor.oligo <- subset(fun.saprotroph.trophic.poor, Category == "Oligotroph" & Mean_RA>0)
x <- subset(fun.saprotroph.trophic.poor.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(fun.saprotroph.trophic.poor.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #6
mean(y)
length(y) #5
sum(x) #0.1403144
sum(y) #0.2109427
length(x[x>0]) #6
length(y[y>0]) #5


p <- ggplot(data =fun.saprotroph.trophic.poor.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


###In nutrient-rich
fun.saprotroph.trophic.rich <- subset(fun.saprotroph.trophic, Category3 == "Rich")
fun.saprotroph.trophic.rich.copio <- subset(fun.saprotroph.trophic.rich, Category == "Copiotroph" & Mean_RA >0)
x <- subset(fun.saprotroph.trophic.rich.copio, Category2=='Non-fer')$Mean_RA
y <- subset(fun.saprotroph.trophic.rich.copio, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #6
mean(y)
length(y) #6
sum(x) #0.002509408
sum(y) #0.01113845
length(x[x>0]) #6
length(y[y>0]) #6

p <- ggplot(data =fun.saprotroph.trophic.rich.copio, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


fun.saprotroph.trophic.rich <- subset(fun.saprotroph.trophic, Category3 == "Rich")
fun.saprotroph.trophic.rich.oligo <- subset(fun.saprotroph.trophic.rich, Category == "Oligotroph")
x <- subset(fun.saprotroph.trophic.rich.oligo, Category2=='Non-fer')$Mean_RA
y <- subset(fun.saprotroph.trophic.rich.oligo, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.04907371
sum(y) #0.03715566
length(x[x>0]) #2
length(y[y>0]) #3

p <- ggplot(data =fun.saprotroph.trophic.rich.oligo, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


write.xlsx(fun.saprotroph.trophic, "Fungal saprotroph.xlsx")
write.xlsx(bac.iron_reducer.trophic, "Bacterial iron_reducer.xlsx")
write.xlsx(bac.syntroph.trophic, "Bacterial syntroph.xlsx")
write.xlsx(bac.methanotroph.trophic, "Bacterial methanotroph.xlsx")
write.xlsx(bac.methanecycling.trophic, "Bacterial methane cycling.xlsx")
write.xlsx(arch.methanogen.trophic, "Archaeal methanogen.xlsx")

##Ammonia-oxidizing archaea
arch.oxidizer <- c("Nitrosopumilus", "Nitrososphaera")
bac.oxidizer <- c("Nitrosospira", "Nitrosomonas")

otu.list.arch.oxidizer <-subset(otu.list, Genus %in% arch.oxidizer)


list.arch.arch.oxidizer <- subset(list.arch.oligo.copio, Genus %in% arch.oxidizer)
head(list.arch.arch.oxidizer)

poor.arch.arch.oxidizer <- list.arch.arch.oxidizer %>% select(OTU_id, Category, CJ1, CJ2, YS1, YS2, CJ1.18, CJ2.18)

rich.arch.arch.oxidizer <- list.arch.arch.oxidizer %>% select(OTU_id, Category, MY1, MY2, NJ1, NJ2)





fun.saprotroph.trophic.nonfer <- subset(fun.saprotroph.trophic, Category2 == "Non-fer")
fun.saprotroph.trophic.nonfer.oligo <- subset(fun.saprotroph.trophic.nonfer, Category == "Oligotroph" & Mean_RA>0)
x <- subset(fun.saprotroph.trophic.nonfer.oligo, Category3=='Poor')$Mean_RA
y <- subset(fun.saprotroph.trophic.nonfer.oligo, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #6
mean(y)
length(y) #2

p <- ggplot(data =fun.saprotroph.trophic.nonfer.oligo, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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



bac.non_acetogenic_syntroph.trophic.nonfer <- subset(bac.non_acetogenic_syntroph.trophic, Category2 == "Non-fer")
bac.non_acetogenic_syntroph.trophic.nonfer.oligo <- subset(bac.non_acetogenic_syntroph.trophic.nonfer, Category == "Oligotroph" & Mean_RA>0)
x <- subset(bac.non_acetogenic_syntroph.trophic.nonfer.oligo, Category3=='Poor')$Mean_RA
y <- subset(bac.non_acetogenic_syntroph.trophic.nonfer.oligo, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #18
mean(y)
length(y) #5

p <- ggplot(data =bac.non_acetogenic_syntroph.trophic.nonfer.oligo, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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


arch.methanogen.trophic.nonfer <- subset(arch.methanogen.trophic, Category2 == "Non-fer")
arch.methanogen.trophic.nonfer.oligo <- subset(arch.methanogen.trophic.nonfer, Category == "Oligotroph" & Mean_RA>0)
x <- subset(arch.methanogen.trophic.nonfer.oligo, Category3=='Poor')$Mean_RA
y <- subset(arch.methanogen.trophic.nonfer.oligo, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #12
mean(y)
length(y) #4

p <- ggplot(data =arch.methanogen.trophic.nonfer.oligo, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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











fun.saprotroph.trophic.rich <- subset(fun.saprotroph.trophic, Category3 == "Rich")
x <- subset(fun.saprotroph.trophic.rich, Category2=='Non-fer')$Mean_RA
y <- subset(fun.saprotroph.trophic.rich, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.04907371
sum(y) #0.03715566
length(x[x>0]) #2
length(y[y>0]) #3

p <- ggplot(data =fun.saprotroph.trophic.rich, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
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




### new figures
fun.saprotroph.trophic$Category <- factor(fun.saprotroph.trophic$Category, levels = c("Oligotroph", "Copiotroph"))
fun.saprotroph.trophic.poor <- subset(fun.saprotroph.trophic, Category3 == "Poor")
x <- subset(fun.saprotroph.trophic.poor, Category2=='Non-fer')$Mean_RA
y <- subset(fun.saprotroph.trophic.poor, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.1531803
sum(y) #0.2216481
length(x[x>0]) #12
length(y[y>0]) #11

p <- ggplot(data =fun.saprotroph.trophic.poor, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


fun.saprotroph.trophic.rich <- subset(fun.saprotroph.trophic, Category3 == "Rich")
x <- subset(fun.saprotroph.trophic.rich, Category2=='Non-fer')$Mean_RA
y <- subset(fun.saprotroph.trophic.rich, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.05158312
sum(y) #0.04829411
length(x[x>0]) #8
length(y[y>0]) #9

p <- ggplot(data =fun.saprotroph.trophic.rich, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


x <- subset(fun.saprotroph.trophic, Category3=='Poor')$Mean_RA
y <- subset(fun.saprotroph.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.3748284
sum(y) #0.09987724
length(x[x>0]) #23
length(y[y>0]) #17

p <- ggplot(data =fun.saprotroph.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


### methanogen

arch.methanogen.trophic$Category <- factor(arch.methanogen.trophic$Category, levels = c("Oligotroph", "Copiotroph"))
arch.methanogen.trophic$Category2 <- factor(arch.methanogen.trophic$Category2, levels = c("Non-fer", "Fer"))

arch.methanogen.trophic.poor <- subset(arch.methanogen.trophic, Category3 == "Poor")
x <- subset(arch.methanogen.trophic.poor, Category2=='Non-fer')$Mean_RA
y <- subset(arch.methanogen.trophic.poor, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.2439061
sum(y) #0.2910446
length(x[x>0]) #15
length(y[y>0]) #13

p <- ggplot(data =arch.methanogen.trophic.poor, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


arch.methanogen.trophic.rich <- subset(arch.methanogen.trophic, Category3 == "Rich")
x <- subset(arch.methanogen.trophic.rich, Category2=='Non-fer')$Mean_RA
y <- subset(arch.methanogen.trophic.rich, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.1030008
sum(y) #0.1297292
length(x[x>0]) #7
length(y[y>0]) #8

p <- ggplot(data =arch.methanogen.trophic.rich, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


x <- subset(arch.methanogen.trophic, Category3=='Poor')$Mean_RA
y <- subset(arch.methanogen.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.5349507
sum(y) #0.2327301
length(x[x>0]) #28
length(y[y>0]) #15

p <- ggplot(data =arch.methanogen.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


### iron_reducer

bac.iron_reducer.trophic$Category <- factor(bac.iron_reducer.trophic$Category, levels = c("Oligotroph", "Copiotroph"))
bac.iron_reducer.trophic$Category2 <- factor(bac.iron_reducer.trophic$Category2, levels = c("Non-fer", "Fer"))

bac.iron_reducer.trophic.poor <- subset(bac.iron_reducer.trophic, Category3 == "Poor")
x <- subset(bac.iron_reducer.trophic.poor, Category2=='Non-fer')$Mean_RA
y <- subset(bac.iron_reducer.trophic.poor, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) # 0.009544603
sum(y) #0.009234389
length(x[x>0]) #20
length(y[y>0]) #14

p <- ggplot(data =bac.iron_reducer.trophic.poor, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


bac.iron_reducer.trophic.rich <- subset(bac.iron_reducer.trophic, Category3 == "Rich")
x <- subset(bac.iron_reducer.trophic.rich, Category2=='Non-fer')$Mean_RA
y <- subset(bac.iron_reducer.trophic.rich, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.005993997
sum(y) #0.003773988
length(x[x>0]) #9
length(y[y>0]) #10

p <- ggplot(data =bac.iron_reducer.trophic.rich, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


x <- subset(bac.iron_reducer.trophic, Category3=='Poor')$Mean_RA
y <- subset(bac.iron_reducer.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.01877899
sum(y) #0.009767985
length(x[x>0]) #34
length(y[y>0]) #19

p <- ggplot(data =bac.iron_reducer.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


### methanotroph

bac.methanotroph.trophic$Category <- factor(bac.methanotroph.trophic$Category, levels = c("Oligotroph", "Copiotroph"))
bac.methanotroph.trophic$Category2 <- factor(bac.methanotroph.trophic$Category2, levels = c("Non-fer", "Fer"))

bac.methanotroph.trophic.poor <- subset(bac.methanotroph.trophic, Category3 == "Poor")
x <- subset(bac.methanotroph.trophic.poor, Category2=='Non-fer')$Mean_RA
y <- subset(bac.methanotroph.trophic.poor, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) # 0.009544603
sum(y) #0.009234389
length(x[x>0]) #34
length(y[y>0]) #19

p <- ggplot(data =bac.methanotroph.trophic.poor, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


bac.methanotroph.trophic.rich <- subset(bac.methanotroph.trophic, Category3 == "Rich")
x <- subset(bac.methanotroph.trophic.rich, Category2=='Non-fer')$Mean_RA
y <- subset(bac.methanotroph.trophic.rich, Category2=='Fer')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.005993997
sum(y) #0.003773988
length(x[x>0]) #9
length(y[y>0]) #10

p <- ggplot(data =bac.methanotroph.trophic.rich, aes(x=Category2, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p


x <- subset(bac.methanotroph.trophic, Category3=='Poor')$Mean_RA
y <- subset(bac.methanotroph.trophic, Category3=='Rich')$Mean_RA
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.00123406
mean(x)
length(x) #5
mean(y)
length(y) #6
sum(x) #0.01877899
sum(y) #0.009767985
length(x[x>0]) #34
length(y[y>0]) #19

p <- ggplot(data =bac.methanotroph.trophic, aes(x=Category3, y=Mean_RA)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  geom_point(position='jitter', alpha=1, aes(colour = Category, shape = Category))+
  scale_colour_manual(labels = c('Oligotroph','Copiotroph'), values = c("#cccc99", "#669933"))+
  scale_shape_manual(labels = c('Oligotroph','Copiotroph'), values = c(18, 17))+
  #geom_text(data=hsd1,aes(x=Group,y=MaxAbund, label=Letter), vjust=-1) +
  xlab(paste('Wilcoxon rank sum test, P-value = ', wil$p.value))+
  ylab("RA \n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")
p



### Supplementary Table S16
arch.methano <- list.arch.methano %>% select(OTU_id, Category, CJ1, YS1, CJ1.18, CJ2, YS2, CJ2.18, MY1, NJ1, MY2, NJ2)
bac.iron_reducer <- list.bac.iron_reducer %>% select(OTU_id, Category, CJ1, YS1, CJ1.18, CJ2, YS2, CJ2.18, MY1, NJ1, MY2, NJ2)
bac.methanotroph <- list.bac.methanotroph %>% select(OTU_id, Category, CJ1, YS1, CJ1.18, CJ2, YS2, CJ2.18, MY1, NJ1, MY2, NJ2)
fun.saprotroph <- list.fun.saprotroph %>% select(OTU_id, Category, CJ1, YS1, CJ1.18, CJ2, YS2, CJ2.18, MY1, NJ1, MY2, NJ2)

write.xlsx(arch.methano, "Archaeal methanogen.xlsx")
write.xlsx(bac.iron_reducer, "bac.iron_reducer.xlsx")
write.xlsx(bac.methanotroph, "bac.methanotroph.xlsx")
write.xlsx(fun.saprotroph, "fun.saprotroph.xlsx")