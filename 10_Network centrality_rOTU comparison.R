### Comparison of centrality among rOTUs sorted by cultural practice
df.norm.degree.plot2
df.norm.betweenness.plot2
df.norm.closeness.plot2


df.norm.degree.rOTU <- subset(df.norm.degree.plot2, rOTU != "Non-dif")
df.norm.degree.rOTU.2017 <- subset(df.norm.degree.rOTU, y2017_all_deg != "NA")
df.norm.degree.rOTU.2018 <- subset(df.norm.degree.rOTU, y2018_all_deg != "NA")

df.norm.betweenness.rOTU <- subset(df.norm.betweenness.plot2, rOTU != "Non-dif")
df.norm.betweenness.rOTU.2017 <- subset(df.norm.betweenness.rOTU, y2017_all_betweenness != "NA")
df.norm.betweenness.rOTU.2018 <- subset(df.norm.betweenness.rOTU, y2018_all_betweenness != "NA")

df.norm.closeness.rOTU <- subset(df.norm.closeness.plot2, rOTU != "Non-dif")
df.norm.closeness.rOTU.2017 <- subset(df.norm.closeness.rOTU, y2017_all_closeness != "NA")
df.norm.closeness.rOTU.2018 <- subset(df.norm.closeness.rOTU, y2018_all_closeness != "NA")


##Kingdom ordering
df.norm.degree.rOTU.2017$Bacteria <-factor(df.norm.degree.rOTU.2017$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))
df.norm.betweenness.rOTU.2017$Bacteria <-factor(df.norm.betweenness.rOTU.2017$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))
df.norm.closeness.rOTU.2017$Bacteria <-factor(df.norm.closeness.rOTU.2017$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))

df.norm.degree.rOTU.2018$Bacteria <-factor(df.norm.degree.rOTU.2018$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))
df.norm.betweenness.rOTU.2018$Bacteria <-factor(df.norm.betweenness.rOTU.2018$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))
df.norm.closeness.rOTU.2018$Bacteria <-factor(df.norm.closeness.rOTU.2018$Bacteria, levels=c("Bacteria", "Archaea","Fungi"))


#2017 degree
#Kruskal-Wallis
test <- aggregate(df.norm.degree.rOTU.2017$y2017_all_deg, by = list(df.norm.degree.rOTU.2017$rOTU), max)

colnames(test) <- c("Group", "maxdegree")

##Kruskal-Wallis
kw<-kruskal.test(y2017_all_deg ~ rOTU, data = df.norm.degree.rOTU.2017)

#library(FSA)
DT = dunnTest(y2017_all_deg ~ rOTU,
              data=df.norm.degree.rOTU.2017,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "rOTU"

p<-ggplot(data=df.norm.degree.rOTU.2017, aes(x=rOTU, y=y2017_all_deg, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2, aes(shape = Bacteria))+
  scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p

#2017 betweenness

#Kruskal-Wallis
test <- aggregate(df.norm.betweenness.rOTU.2017$y2017_all_betweenness, by = list(df.norm.betweenness.rOTU.2017$rOTU), max)

colnames(test) <- c("Group", "maxbetweenness")

##Kruskal-Wallis
kw<-kruskal.test(y2017_all_betweenness ~ rOTU, data = df.norm.betweenness.rOTU.2017)

#library(FSA)
DT = dunnTest(y2017_all_betweenness ~ rOTU,
              data=df.norm.betweenness.rOTU.2017,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "rOTU"

p<-ggplot(data=df.norm.betweenness.rOTU.2017, aes(x=rOTU, y=y2017_all_betweenness, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2, aes(shape = Bacteria))+
  scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=rOTU,y=maxbetweenness, label=Letter), vjust=-1)+
  labs(x="", y ="betweenness") +
  scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p

#2017 closeness
#Kruskal-Wallis
test <- aggregate(df.norm.closeness.rOTU.2017$y2017_all_closeness, by = list(df.norm.closeness.rOTU.2017$rOTU), max)

colnames(test) <- c("Group", "maxcloseness")

##Kruskal-Wallis
df.norm.closeness.rOTU.2017$rOTU <- as.factor(df.norm.closeness.rOTU.2017$rOTU)

kw<-kruskal.test(y2017_all_closeness ~ rOTU, data = df.norm.closeness.rOTU.2017)

#library(FSA)
DT = dunnTest(y2017_all_closeness ~ rOTU,
              data=df.norm.closeness.rOTU.2017,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "rOTU"

p<-ggplot(data=df.norm.closeness.rOTU.2017, aes(x=rOTU, y=y2017_all_closeness, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2, aes(shape = Bacteria))+
  scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=rOTU,y=maxcloseness, label=Letter), vjust=-1)+
  labs(x="", y ="closeness") +
  scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p




#2018 degree
#Kruskal-Wallis
test <- aggregate(df.norm.degree.rOTU.2018$y2018_all_deg, by = list(df.norm.degree.rOTU.2018$rOTU), max)

colnames(test) <- c("Group", "maxdegree")

##Kruskal-Wallis
kw<-kruskal.test(y2018_all_deg ~ rOTU, data = df.norm.degree.rOTU.2018)
kw<-aov(y2018_all_deg ~ rOTU, data = df.norm.degree.rOTU.2018)
#library(FSA)
DT = dunnTest(y2018_all_deg ~ rOTU,
              data=df.norm.degree.rOTU.2018,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "rOTU"

p<-ggplot(data=df.norm.degree.rOTU.2018, aes(x=rOTU, y=y2018_all_deg, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2, aes(shape = Bacteria))+
  scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p

#2018 betweenness
#Kruskal-Wallis
test <- aggregate(df.norm.betweenness.rOTU.2018$y2018_all_betweenness, by = list(df.norm.betweenness.rOTU.2018$rOTU), max)

colnames(test) <- c("Group", "maxbetweenness")

##Kruskal-Wallis
kw<-kruskal.test(y2018_all_betweenness ~ rOTU, data = df.norm.betweenness.rOTU.2018)

#library(FSA)
DT = dunnTest(y2018_all_betweenness ~ rOTU,
              data=df.norm.betweenness.rOTU.2018,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "rOTU"

p<-ggplot(data=df.norm.betweenness.rOTU.2018, aes(x=rOTU, y=y2018_all_betweenness, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2, aes(shape = Bacteria))+
  scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=rOTU,y=maxbetweenness, label=Letter), vjust=-1)+
  labs(x="", y ="betweenness") +
  scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p

#2018 closeness
#Kruskal-Wallis
test <- aggregate(df.norm.closeness.rOTU.2018$y2018_all_closeness, by = list(df.norm.closeness.rOTU.2018$rOTU), max)

colnames(test) <- c("Group", "maxcloseness")

##Kruskal-Wallis
df.norm.closeness.rOTU.2018$rOTU <- as.factor(df.norm.closeness.rOTU.2018$rOTU)

kw<-kruskal.test(y2018_all_closeness ~ rOTU, data = df.norm.closeness.rOTU.2018)

#library(FSA)
DT = dunnTest(y2018_all_closeness ~ rOTU,
              data=df.norm.closeness.rOTU.2018,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "rOTU"

p<-ggplot(data=df.norm.closeness.rOTU.2018, aes(x=rOTU, y=y2018_all_closeness, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2, aes(shape = Bacteria))+
  scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=rOTU,y=maxcloseness, label=Letter), vjust=-1)+
  labs(x="", y ="closeness") +
  scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p

