## Comparison degree, betweenness, closeness
# Degree
df.CJ1_all_deg<-data.frame(CJ1_all_deg)
df.CJ1_all_deg$category <- "non-fer"
df.CJ1_all_deg$field <- "CJ1"
df.CJ1_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.CJ1_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ1_all_deg)),'Archaea', 'Fungi'))
head(df.CJ1_all_deg)
names(df.CJ1_all_deg)[1] <- "Degree"
mean(df.CJ1_all_deg$Degree)

df.CJ2_all_deg<-data.frame(CJ2_all_deg)
df.CJ2_all_deg$category <- "fer"
df.CJ2_all_deg$field <- "CJ2"
df.CJ2_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.CJ2_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ2_all_deg)),'Archaea', 'Fungi'))
names(df.CJ2_all_deg)[1] <- "Degree"
mean(df.CJ2_all_deg$Degree)

df.CJ_all_deg <- rbind(df.CJ1_all_deg, df.CJ2_all_deg)

x <- subset(df.CJ_all_deg, category=='non-fer')$Degree
y <- subset(df.CJ_all_deg, category=='fer')$Degree
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.CJ1_18_all_deg<-data.frame(CJ1_18_all_deg)
df.CJ1_18_all_deg$category <- "non-fer"
df.CJ1_18_all_deg$field <- "CJ1_18"
df.CJ1_18_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.CJ1_18_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ1_18_all_deg)),'Archaea', 'Fungi'))

names(df.CJ1_18_all_deg)[1] <- "Degree"
mean(df.CJ1_18_all_deg$Degree)

df.CJ2_18_all_deg<-data.frame(CJ2_18_all_deg)
df.CJ2_18_all_deg$category <- "fer"
df.CJ2_18_all_deg$field <- "CJ2_18"
df.CJ2_18_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.CJ2_18_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ2_18_all_deg)),'Archaea', 'Fungi'))

names(df.CJ2_18_all_deg)[1] <- "Degree"
mean(df.CJ2_18_all_deg$Degree)

df.CJ_18_all_deg <- rbind(df.CJ1_18_all_deg, df.CJ2_18_all_deg)

x <- subset(df.CJ_18_all_deg, category=='non-fer')$Degree
y <- subset(df.CJ_18_all_deg, category=='fer')$Degree
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value



df.MY2_all_deg<-data.frame(MY2_all_deg)
df.MY2_all_deg$category <- "non-fer"
df.MY2_all_deg$field <- "MY2"
df.MY2_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.MY2_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.MY2_all_deg)),'Archaea', 'Fungi'))

names(df.MY2_all_deg)[1] <- "Degree"
mean(df.MY2_all_deg$Degree)

df.MY1_all_deg<-data.frame(MY1_all_deg)
df.MY1_all_deg$category <- "fer"
df.MY1_all_deg$field <- "MY1"
df.MY1_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.MY1_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.MY1_all_deg)),'Archaea', 'Fungi'))

names(df.MY1_all_deg)[1] <- "Degree"
mean(df.MY1_all_deg$Degree)

df.MY_all_deg <- rbind(df.MY2_all_deg, df.MY1_all_deg)

x <- subset(df.MY_all_deg, category=='non-fer')$Degree
y <- subset(df.MY_all_deg, category=='fer')$Degree
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.NJ2_all_deg<-data.frame(NJ2_all_deg)
df.NJ2_all_deg$category <- "non-fer"
df.NJ2_all_deg$field <- "NJ2"
df.NJ2_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.NJ2_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.NJ2_all_deg)),'Archaea', 'Fungi'))

names(df.NJ2_all_deg)[1] <- "Degree"
mean(df.NJ2_all_deg$Degree)

df.NJ1_all_deg<-data.frame(NJ1_all_deg)
df.NJ1_all_deg$category <- "fer"
df.NJ1_all_deg$field <- "NJ1"
df.NJ1_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.NJ1_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.NJ1_all_deg)),'Archaea', 'Fungi'))

names(df.NJ1_all_deg)[1] <- "Degree"
mean(df.NJ1_all_deg$Degree)

df.NJ_all_deg <- rbind(df.NJ2_all_deg, df.NJ1_all_deg)

x <- subset(df.NJ_all_deg, category=='non-fer')$Degree
y <- subset(df.NJ_all_deg, category=='fer')$Degree
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.YS1_all_deg<-data.frame(YS1_all_deg)
df.YS1_all_deg$category <- "non-fer"
df.YS1_all_deg$field <- "YS1"
df.YS1_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.YS1_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.YS1_all_deg)),'Archaea', 'Fungi'))

names(df.YS1_all_deg)[1] <- "Degree"
mean(df.YS1_all_deg$Degree)

df.YS2_all_deg<-data.frame(YS2_all_deg)
df.YS2_all_deg$category <- "fer"
df.YS2_all_deg$field <- "YS2"
df.YS2_all_deg$kingdom <- ifelse(grepl("^B",rownames(df.YS2_all_deg)),'Bacteria',ifelse(grepl("^A",rownames(df.YS2_all_deg)),'Archaea', 'Fungi'))

names(df.YS2_all_deg)[1] <- "Degree"
mean(df.YS2_all_deg$Degree)

df.YS_all_deg <- rbind(df.YS1_all_deg, df.YS2_all_deg)

x <- subset(df.YS_all_deg, category=='non-fer')$Degree
y <- subset(df.YS_all_deg, category=='fer')$Degree
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


## 
df.all_deg <- rbind(df.CJ_all_deg, df.CJ_18_all_deg, df.MY_all_deg, df.NJ_all_deg, df.YS_all_deg)
x <- subset(df.all_deg, category=='non-fer')$Degree
mean(subset(df.all_deg, category=='non-fer')$Degree)
y <- subset(df.all_deg, category=='fer')$Degree
mean(subset(df.all_deg, category=='fer')$Degree)
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value

write.table(df.all_deg, "df.all_deg.tsv", sep= '\t', quote=F)


##Betweenness
df.CJ1_all_betweenness<-data.frame(CJ1_all_betweenness)
df.CJ1_all_betweenness$category <- "non-fer"
df.CJ1_all_betweenness$field <- "CJ1"
df.CJ1_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.CJ1_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ1_all_betweenness)),'Archaea', 'Fungi'))
head(df.CJ1_all_betweenness)
names(df.CJ1_all_betweenness)[1] <- "Betweenness"
mean(df.CJ1_all_betweenness$Betweenness)

df.CJ2_all_betweenness<-data.frame(CJ2_all_betweenness)
df.CJ2_all_betweenness$category <- "fer"
df.CJ2_all_betweenness$field <- "CJ2"
df.CJ2_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.CJ2_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ2_all_betweenness)),'Archaea', 'Fungi'))
names(df.CJ2_all_betweenness)[1] <- "Betweenness"
mean(df.CJ2_all_betweenness$Betweenness)

df.CJ_all_betweenness <- rbind(df.CJ1_all_betweenness, df.CJ2_all_betweenness)

x <- subset(df.CJ_all_betweenness, category=='non-fer')$Betweenness
y <- subset(df.CJ_all_betweenness, category=='fer')$Betweenness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.CJ1_18_all_betweenness<-data.frame(CJ1_18_all_betweenness)
df.CJ1_18_all_betweenness$category <- "non-fer"
df.CJ1_18_all_betweenness$field <- "CJ1_18"
df.CJ1_18_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.CJ1_18_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ1_18_all_betweenness)),'Archaea', 'Fungi'))

names(df.CJ1_18_all_betweenness)[1] <- "Betweenness"
mean(df.CJ1_18_all_betweenness$Betweenness)

df.CJ2_18_all_betweenness<-data.frame(CJ2_18_all_betweenness)
df.CJ2_18_all_betweenness$category <- "fer"
df.CJ2_18_all_betweenness$field <- "CJ2_18"
df.CJ2_18_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.CJ2_18_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ2_18_all_betweenness)),'Archaea', 'Fungi'))

names(df.CJ2_18_all_betweenness)[1] <- "Betweenness"
mean(df.CJ2_18_all_betweenness$Betweenness)

df.CJ_18_all_betweenness <- rbind(df.CJ1_18_all_betweenness, df.CJ2_18_all_betweenness)

x <- subset(df.CJ_18_all_betweenness, category=='non-fer')$Betweenness
y <- subset(df.CJ_18_all_betweenness, category=='fer')$Betweenness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value



df.MY2_all_betweenness<-data.frame(MY2_all_betweenness)
df.MY2_all_betweenness$category <- "non-fer"
df.MY2_all_betweenness$field <- "MY2"
df.MY2_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.MY2_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.MY2_all_betweenness)),'Archaea', 'Fungi'))

names(df.MY2_all_betweenness)[1] <- "Betweenness"
mean(df.MY2_all_betweenness$Betweenness)

df.MY1_all_betweenness<-data.frame(MY1_all_betweenness)
df.MY1_all_betweenness$category <- "fer"
df.MY1_all_betweenness$field <- "MY1"
df.MY1_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.MY1_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.MY1_all_betweenness)),'Archaea', 'Fungi'))

names(df.MY1_all_betweenness)[1] <- "Betweenness"
mean(df.MY1_all_betweenness$Betweenness)

df.MY_all_betweenness <- rbind(df.MY2_all_betweenness, df.MY1_all_betweenness)

x <- subset(df.MY_all_betweenness, category=='non-fer')$Betweenness
y <- subset(df.MY_all_betweenness, category=='fer')$Betweenness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.NJ2_all_betweenness<-data.frame(NJ2_all_betweenness)
df.NJ2_all_betweenness$category <- "non-fer"
df.NJ2_all_betweenness$field <- "NJ2"
df.NJ2_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.NJ2_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.NJ2_all_betweenness)),'Archaea', 'Fungi'))

names(df.NJ2_all_betweenness)[1] <- "Betweenness"
mean(df.NJ2_all_betweenness$Betweenness)

df.NJ1_all_betweenness<-data.frame(NJ1_all_betweenness)
df.NJ1_all_betweenness$category <- "fer"
df.NJ1_all_betweenness$field <- "NJ1"
df.NJ1_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.NJ1_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.NJ1_all_betweenness)),'Archaea', 'Fungi'))

names(df.NJ1_all_betweenness)[1] <- "Betweenness"
mean(df.NJ1_all_betweenness$Betweenness)

df.NJ_all_betweenness <- rbind(df.NJ2_all_betweenness, df.NJ1_all_betweenness)

x <- subset(df.NJ_all_betweenness, category=='non-fer')$Betweenness
y <- subset(df.NJ_all_betweenness, category=='fer')$Betweenness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.YS1_all_betweenness<-data.frame(YS1_all_betweenness)
df.YS1_all_betweenness$category <- "non-fer"
df.YS1_all_betweenness$field <- "YS1"
df.YS1_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.YS1_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.YS1_all_betweenness)),'Archaea', 'Fungi'))

names(df.YS1_all_betweenness)[1] <- "Betweenness"
mean(df.YS1_all_betweenness$Betweenness)

df.YS2_all_betweenness<-data.frame(YS2_all_betweenness)
df.YS2_all_betweenness$category <- "fer"
df.YS2_all_betweenness$field <- "YS2"
df.YS2_all_betweenness$kingdom <- ifelse(grepl("^B",rownames(df.YS2_all_betweenness)),'Bacteria',ifelse(grepl("^A",rownames(df.YS2_all_betweenness)),'Archaea', 'Fungi'))

names(df.YS2_all_betweenness)[1] <- "Betweenness"
mean(df.YS2_all_betweenness$Betweenness)

df.YS_all_betweenness <- rbind(df.YS1_all_betweenness, df.YS2_all_betweenness)

x <- subset(df.YS_all_betweenness, category=='non-fer')$Betweenness
y <- subset(df.YS_all_betweenness, category=='fer')$Betweenness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


## 
df.all_betweenness <- rbind(df.CJ_all_betweenness, df.CJ_18_all_betweenness, df.MY_all_betweenness, df.NJ_all_betweenness, df.YS_all_betweenness)
x <- subset(df.all_betweenness, category=='non-fer')$Betweenness
mean(subset(df.all_betweenness, category=='non-fer')$Betweenness)
y <- subset(df.all_betweenness, category=='fer')$Betweenness
mean(subset(df.all_betweenness, category=='fer')$Betweenness)
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


write.table(df.all_betweenness, "df.all_betweenness.tsv", sep= '\t', quote=F)
##Closeness
df.CJ1_all_closeness<-data.frame(CJ1_all_closeness)
df.CJ1_all_closeness$category <- "non-fer"
df.CJ1_all_closeness$field <- "CJ1"
df.CJ1_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.CJ1_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ1_all_closeness)),'Archaea', 'Fungi'))
head(df.CJ1_all_closeness)
names(df.CJ1_all_closeness)[1] <- "Closeness"
mean(df.CJ1_all_closeness$Closeness)

df.CJ2_all_closeness<-data.frame(CJ2_all_closeness)
df.CJ2_all_closeness$category <- "fer"
df.CJ2_all_closeness$field <- "CJ2"
df.CJ2_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.CJ2_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ2_all_closeness)),'Archaea', 'Fungi'))
names(df.CJ2_all_closeness)[1] <- "Closeness"
mean(df.CJ2_all_closeness$Closeness)

df.CJ_all_closeness <- rbind(df.CJ1_all_closeness, df.CJ2_all_closeness)

x <- subset(df.CJ_all_closeness, category=='non-fer')$Closeness
y <- subset(df.CJ_all_closeness, category=='fer')$Closeness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.CJ1_18_all_closeness<-data.frame(CJ1_18_all_closeness)
df.CJ1_18_all_closeness$category <- "non-fer"
df.CJ1_18_all_closeness$field <- "CJ1_18"
df.CJ1_18_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.CJ1_18_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ1_18_all_closeness)),'Archaea', 'Fungi'))

names(df.CJ1_18_all_closeness)[1] <- "Closeness"
mean(df.CJ1_18_all_closeness$Closeness)

df.CJ2_18_all_closeness<-data.frame(CJ2_18_all_closeness)
df.CJ2_18_all_closeness$category <- "fer"
df.CJ2_18_all_closeness$field <- "CJ2_18"
df.CJ2_18_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.CJ2_18_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.CJ2_18_all_closeness)),'Archaea', 'Fungi'))

names(df.CJ2_18_all_closeness)[1] <- "Closeness"
mean(df.CJ2_18_all_closeness$Closeness)

df.CJ_18_all_closeness <- rbind(df.CJ1_18_all_closeness, df.CJ2_18_all_closeness)

x <- subset(df.CJ_18_all_closeness, category=='non-fer')$Closeness
y <- subset(df.CJ_18_all_closeness, category=='fer')$Closeness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value



df.MY2_all_closeness<-data.frame(MY2_all_closeness)
df.MY2_all_closeness$category <- "non-fer"
df.MY2_all_closeness$field <- "MY2"
df.MY2_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.MY2_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.MY2_all_closeness)),'Archaea', 'Fungi'))

names(df.MY2_all_closeness)[1] <- "Closeness"
mean(df.MY2_all_closeness$Closeness)

df.MY1_all_closeness<-data.frame(MY1_all_closeness)
df.MY1_all_closeness$category <- "fer"
df.MY1_all_closeness$field <- "MY1"
df.MY1_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.MY1_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.MY1_all_closeness)),'Archaea', 'Fungi'))

names(df.MY1_all_closeness)[1] <- "Closeness"
mean(df.MY1_all_closeness$Closeness)

df.MY_all_closeness <- rbind(df.MY2_all_closeness, df.MY1_all_closeness)

x <- subset(df.MY_all_closeness, category=='non-fer')$Closeness
y <- subset(df.MY_all_closeness, category=='fer')$Closeness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.NJ2_all_closeness<-data.frame(NJ2_all_closeness)
df.NJ2_all_closeness$category <- "non-fer"
df.NJ2_all_closeness$field <- "NJ2"
df.NJ2_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.NJ2_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.NJ2_all_closeness)),'Archaea', 'Fungi'))

names(df.NJ2_all_closeness)[1] <- "Closeness"
mean(df.NJ2_all_closeness$Closeness)

df.NJ1_all_closeness<-data.frame(NJ1_all_closeness)
df.NJ1_all_closeness$category <- "fer"
df.NJ1_all_closeness$field <- "NJ1"
df.NJ1_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.NJ1_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.NJ1_all_closeness)),'Archaea', 'Fungi'))

names(df.NJ1_all_closeness)[1] <- "Closeness"
mean(df.NJ1_all_closeness$Closeness)

df.NJ_all_closeness <- rbind(df.NJ2_all_closeness, df.NJ1_all_closeness)

x <- subset(df.NJ_all_closeness, category=='non-fer')$Closeness
y <- subset(df.NJ_all_closeness, category=='fer')$Closeness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


df.YS1_all_closeness<-data.frame(YS1_all_closeness)
df.YS1_all_closeness$category <- "non-fer"
df.YS1_all_closeness$field <- "YS1"
df.YS1_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.YS1_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.YS1_all_closeness)),'Archaea', 'Fungi'))

names(df.YS1_all_closeness)[1] <- "Closeness"
mean(df.YS1_all_closeness$Closeness)

df.YS2_all_closeness<-data.frame(YS2_all_closeness)
df.YS2_all_closeness$category <- "fer"
df.YS2_all_closeness$field <- "YS2"
df.YS2_all_closeness$kingdom <- ifelse(grepl("^B",rownames(df.YS2_all_closeness)),'Bacteria',ifelse(grepl("^A",rownames(df.YS2_all_closeness)),'Archaea', 'Fungi'))

names(df.YS2_all_closeness)[1] <- "Closeness"
mean(df.YS2_all_closeness$Closeness)

df.YS_all_closeness <- rbind(df.YS1_all_closeness, df.YS2_all_closeness)

x <- subset(df.YS_all_closeness, category=='non-fer')$Closeness
y <- subset(df.YS_all_closeness, category=='fer')$Closeness
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


## 
df.all_closeness <- rbind(df.CJ_all_closeness, df.CJ_18_all_closeness, df.MY_all_closeness, df.NJ_all_closeness, df.YS_all_closeness)
x <- subset(df.all_closeness, category=='non-fer')$Closeness
mean(subset(df.all_closeness, category=='non-fer')$Closeness)
y <- subset(df.all_closeness, category=='fer')$Closeness
mean(subset(df.all_closeness, category=='fer')$Closeness)
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


write.table(df.all_closeness, "df.all_closeness.tsv", sep= '\t', quote=F)

library(ggsignif)
##Plotting
p <- ggplot(data = df.CJ_all_deg, aes(x=field, y=Degree)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Degree (CJ)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("CJ1", "CJ2")), step_increase = 0.1, map_signif_level=TRUE)

p

p <- ggplot(data = df.CJ_18_all_deg, aes(x=field, y=Degree)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Degree (CJ_18)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("CJ1_18", "CJ2_18")), step_increase = 0.1, map_signif_level=TRUE)

p

p <- ggplot(data = df.MY_all_deg, aes(x=field, y=Degree)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Degree (MY)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("MY1", "MY2")), step_increase = 0.1, map_signif_level=TRUE)

p

p <- ggplot(data = df.NJ_all_deg, aes(x=field, y=Degree)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Degree (NJ)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("NJ1", "NJ2")), step_increase = 0.1, map_signif_level=TRUE)

p

p <- ggplot(data = df.YS_all_deg, aes(x=field, y=Degree)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Degree (YS)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("YS1", "YS2")), step_increase = 0.1, map_signif_level=TRUE)

p



p <- ggplot(data = df.CJ_all_betweenness, aes(x=field, y=Betweenness)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Betweenness (CJ)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("non-fer", "fer")), step_increase = 0.1, map_signif_level=TRUE)

p

p <- ggplot(data = df.CJ_18_all_betweenness, aes(x=field, y=Betweenness)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Betweenness (CJ_18)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("non-fer", "fer")), step_increase = 0.1, map_signif_level=TRUE)

p

p <- ggplot(data = df.MY_all_betweenness, aes(x=field, y=Betweenness)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Betweenness (MY)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("non-fer", "fer")), step_increase = 0.1, map_signif_level=TRUE)

p

p <- ggplot(data = df.NJ_all_betweenness, aes(x=field, y=Betweenness)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Betweenness (NJ)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("non-fer", "fer")), step_increase = 0.1, map_signif_level=TRUE)

p

p <- ggplot(data = df.YS_all_betweenness, aes(x=field, y=Betweenness)) + geom_boxplot(fill="white", width = 0.8) +
  theme_bw() +
  theme(aspect.ratio = 2)+
  #geom_point(position='jitter',shape=1, alpha=.5)+
  #geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  #xlab(paste('Wilcox rank sum test, P-value = ', wil$p.value))+
  ylab("Betweenness (YS)\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  geom_signif(comparisons = list(c("non-fer", "fer")), step_increase = 0.1, map_signif_level=TRUE)

p

dev.off()


## Proportion of positive and negative associations
draw_edge_proportion <- function(filename){
  edge.mer <- read.table(file = filename, sep = '\t', header = TRUE)
  head(edge.mer)
  
  edge.mer$Association <- ifelse(grepl("^B",edge.mer$from) & grepl("^F",edge.mer$to) | grepl("^F",edge.mer$from) & grepl("^B",edge.mer$to),'BF',
                                 ifelse(grepl("^A",edge.mer$from) & grepl("^F",edge.mer$to) | grepl("^F",edge.mer$from) & grepl("^A",edge.mer$to),'AF',
                                        ifelse(grepl("^B",edge.mer$from) & grepl("^A",edge.mer$to) | grepl("^A",edge.mer$from) & grepl("^B",edge.mer$to),'AB',
                                               ifelse(grepl("^B",edge.mer$from) & grepl("^B",edge.mer$to),'B', 
                                                      ifelse(grepl("^F",edge.mer$from) & grepl("^F",edge.mer$to),'F','A')))))
  edge.mer
  
  edge.mer$PN <- ifelse(edge.mer$cor > 0,'positive','negative')
  edge.mer$count <- 1
  edge.mer %>% arrange(desc(cor))
  edge.mer %>% arrange(desc(Association))
  
  df.edge.mer <- edge.mer %>% group_by(Association,PN) %>% summarise(count=sum(count))
  
  df.edge.mer$Proportion <- df.edge.mer$count / sum(df.edge.mer$count)
  df.edge.mer
  print(df.edge.mer %>% filter(Association == 'B') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'F') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'A') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'BF') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'AB') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'AF') %>% mutate(ratio = Proportion / sum(Proportion)))
  
  df.edge.mer$Association <- factor(df.edge.mer$Association, levels = c('B', 'A', 'F', 'AB', 'AF', 'BF'))
  
  p.edge.mer <- ggplot(df.edge.mer, aes(x=Association, y = Proportion, fill = PN)) + 
    geom_bar(stat="identity", width = 0.8, position = 'stack', colour="black") +
    #scale_fill_discrete() +
    scale_fill_manual(values = c('#333333','#99CCFF')) +
    
    xlab('')+
    ylab("Proportion of edges \n") +
    #ggtitle("Phylum Community Composition by Sample \n") +
    ## adjust positions
    guides(fill = guide_legend(ncol = 2,reverse = T))+
    theme(legend.position="bottom",legend.title=element_blank()) +
    theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=15, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(),plot.background=element_blank())
  
  return(p.edge.mer)
}

draw_edge_proportion('sig_correlation_CJ1_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_CJ2_0.6_p0.05.tsv')

draw_edge_proportion('sig_correlation_CJ1_18_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_CJ2_18_0.6_p0.05.tsv')

draw_edge_proportion('sig_correlation_MY1_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_MY2_0.6_p0.05.tsv')

draw_edge_proportion('sig_correlation_NJ1_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_NJ2_0.6_p0.05.tsv')

draw_edge_proportion('sig_correlation_YS1_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_YS2_0.6_p0.05.tsv')


## Network properties of kingdoms
df.CJ1_all_deg$kingdom <- factor(df.CJ1_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ2_all_deg$kingdom <- factor(df.CJ2_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ1_18_all_deg$kingdom <- factor(df.CJ1_18_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ2_18_all_deg$kingdom <- factor(df.CJ2_18_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))

df.MY1_all_deg$kingdom <- factor(df.MY1_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.MY2_all_deg$kingdom <- factor(df.MY2_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))

df.NJ1_all_deg$kingdom <- factor(df.NJ1_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.NJ2_all_deg$kingdom <- factor(df.NJ2_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))

df.YS1_all_deg$kingdom <- factor(df.YS1_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.YS2_all_deg$kingdom <- factor(df.YS2_all_deg$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))

deg_comparison <- function(df.deg){
test <- aggregate(df.deg$Degree, by = list(df.deg$kingdom), max)

colnames(test) <- c("Group", "maxdegree")

kw<-kruskal.test(Degree ~ kingdom, data = df.deg)

#library(FSA)
DT = dunnTest(Degree ~ kingdom,
              data=df.deg,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "kingdom"

p<-ggplot(data=df.deg, aes(x=kingdom, y=Degree, fill=kingdom)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  geom_text(data=hsd1,aes(x=kingdom,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")

return(p)
}


deg_comparison(df.CJ1_all_deg)
deg_comparison(df.CJ1_18_all_deg)
deg_comparison(df.CJ2_all_deg)
deg_comparison(df.CJ2_18_all_deg)

deg_comparison(df.MY1_all_deg)
deg_comparison(df.MY2_all_deg)

deg_comparison(df.NJ1_all_deg)
deg_comparison(df.NJ2_all_deg)

deg_comparison(df.YS1_all_deg)
deg_comparison(df.YS2_all_deg)

## Betweenness
df.CJ1_all_betweenness$kingdom <- factor(df.CJ1_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ2_all_betweenness$kingdom <- factor(df.CJ2_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ1_18_all_betweenness$kingdom <- factor(df.CJ1_18_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ2_18_all_betweenness$kingdom <- factor(df.CJ2_18_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.MY1_all_betweenness$kingdom <- factor(df.MY1_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.MY2_all_betweenness$kingdom <- factor(df.MY2_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))

df.NJ1_all_betweenness$kingdom <- factor(df.NJ1_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.NJ2_all_betweenness$kingdom <- factor(df.NJ2_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))

df.YS1_all_betweenness$kingdom <- factor(df.YS1_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.YS2_all_betweenness$kingdom <- factor(df.YS2_all_betweenness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))



betw_comparison <- function(df.betw){
test <- aggregate(df.betw$Betweenness, by = list(df.betw$kingdom), max)

colnames(test) <- c("Group", "maxbetweenness")

kw<-kruskal.test(Betweenness ~ kingdom, data = df.betw)

#library(FSA)
DT = dunnTest(Betweenness ~ kingdom,
              data=df.betw,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "kingdom"

p<-ggplot(data=df.betw, aes(x=kingdom, y=Betweenness, fill=kingdom)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  geom_text(data=hsd1,aes(x=kingdom,y=maxbetweenness, label=Letter), vjust=-1)+
  labs(x="", y ="Betweenness") +
  scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")

return(p) }

betw_comparison(df.CJ1_all_betweenness)
betw_comparison(df.CJ2_all_betweenness)

betw_comparison(df.CJ1_18_all_betweenness)
betw_comparison(df.CJ2_18_all_betweenness)

betw_comparison(df.MY1_all_betweenness)
betw_comparison(df.MY2_all_betweenness)

betw_comparison(df.NJ1_all_betweenness)
betw_comparison(df.NJ2_all_betweenness)

betw_comparison(df.YS1_all_betweenness)
betw_comparison(df.YS2_all_betweenness)
##Closeness
df.CJ1_all_closeness$kingdom <- factor(df.CJ1_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ2_all_closeness$kingdom <- factor(df.CJ2_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ1_18_all_closeness$kingdom <- factor(df.CJ1_18_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.CJ2_18_all_closeness$kingdom <- factor(df.CJ2_18_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.MY1_all_closeness$kingdom <- factor(df.MY1_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.MY2_all_closeness$kingdom <- factor(df.MY2_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.NJ1_all_closeness$kingdom <- factor(df.NJ1_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.NJ2_all_closeness$kingdom <- factor(df.NJ2_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.YS1_all_closeness$kingdom <- factor(df.YS1_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))
df.YS2_all_closeness$kingdom <- factor(df.YS2_all_closeness$kingdom, levels = c("Bacteria", "Archaea", "Fungi"))


close_comparison<-function(df.close){test <- aggregate(df.close$Closeness, by = list(df.close$kingdom), max)

colnames(test) <- c("Group", "maxcloseness")

kw<-kruskal.test(Closeness ~ kingdom, data = df.close)

#library(FSA)
DT = dunnTest(Closeness ~ kingdom,
              data=df.close,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "kingdom"

p<-ggplot(data=df.close, aes(x=kingdom, y=Closeness, fill=kingdom)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  geom_text(data=hsd1,aes(x=kingdom,y=maxcloseness, label=Letter), vjust=-1)+
  labs(x="", y ="Closeness") +
  scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")

return(p)
}

close_comparison(df.CJ1_all_closeness)
close_comparison(df.CJ2_all_closeness)
close_comparison(df.CJ1_18_all_closeness)
close_comparison(df.CJ2_18_all_closeness)
close_comparison(df.MY1_all_closeness)
close_comparison(df.MY2_all_closeness)
close_comparison(df.NJ1_all_closeness)
close_comparison(df.NJ2_all_closeness)
close_comparison(df.YS1_all_closeness)
close_comparison(df.YS2_all_closeness)


### Comparison network properties each kingdom and field
library(ggsignif)

comparison_field_kingdom_deg<- function(df.com, kin, f1, f2){
  df.com.t <- subset(df.com, kingdom == kin)
  
  x <- subset(df.com.t, field==f1)$Degree
  y <- subset(df.com.t, field==f2)$Degree
  wil<-wilcox.test(x, y, conf.int = TRUE)
  p_value<-wil$p.value
  
  p<-ggplot(data=df.com.t, aes(x=field, y=Degree, fill=kingdom)) + geom_boxplot() +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    geom_signif(comparisons = list(c(f1, f2)), step_increase = 0.1, map_signif_level=TRUE)+
    xlab(paste('Mann Whitney U test, P-value = ', p_value))+
    ylab("Degree\n") +
    scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
    theme(aspect.ratio = 2)+
    theme(plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
  return(p)
}

comparison_field_kingdom_between<- function(df.com, kin, f1, f2){
  df.com.t <- subset(df.com, kingdom == kin)
  
  x <- subset(df.com.t, field==f1)$Betweenness
  y <- subset(df.com.t, field==f2)$Betweenness
  wil<-wilcox.test(x, y, conf.int = TRUE)
  p_value<-wil$p.value
  
  p<-ggplot(data=df.com.t, aes(x=field, y=Betweenness, fill=kingdom)) + geom_boxplot() +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    geom_signif(comparisons = list(c(f1, f2)), step_increase = 0.1, map_signif_level=TRUE)+
    xlab(paste('Mann Whitney U test, P-value = ', p_value))+
    ylab("Betweenness\n") +
    scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
    theme(aspect.ratio = 2)+
    theme(plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
  return(p)
}

comparison_field_kingdom_close<- function(df.com, kin, f1, f2){
df.com.t <- subset(df.com, kingdom == kin)

x <- subset(df.com.t, field==f1)$Closeness
y <- subset(df.com.t, field==f2)$Closeness
wil<-wilcox.test(x, y, conf.int = TRUE)
p_value<-wil$p.value

p<-ggplot(data=df.com.t, aes(x=field, y=Closeness, fill=kingdom)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  geom_signif(comparisons = list(c(f1, f2)), step_increase = 0.1, map_signif_level=TRUE)+
  xlab(paste('Mann Whitney U test, P-value = ', p_value))+
  ylab("Closeness\n") +
  scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
return(p)
}


comparison_field_kingdom_deg(df.CJ_all_deg, "Bacteria", "CJ1", "CJ2")
comparison_field_kingdom_deg(df.CJ_all_deg, "Archaea", "CJ1", "CJ2")
comparison_field_kingdom_deg(df.CJ_all_deg, "Fungi", "CJ1", "CJ2")

comparison_field_kingdom_deg(df.MY_all_deg, "Bacteria", "MY1", "MY2")
comparison_field_kingdom_deg(df.MY_all_deg, "Archaea", "MY1", "MY2")
comparison_field_kingdom_deg(df.MY_all_deg, "Fungi", "MY1", "MY2")

comparison_field_kingdom_deg(df.YS_all_deg, "Bacteria", "YS1", "YS2")
comparison_field_kingdom_deg(df.YS_all_deg, "Archaea", "YS1", "YS2")
comparison_field_kingdom_deg(df.YS_all_deg, "Fungi", "YS1", "YS2")

comparison_field_kingdom_deg(df.NJ_all_deg, "Bacteria", "NJ1", "NJ2")
comparison_field_kingdom_deg(df.NJ_all_deg, "Archaea", "NJ1", "NJ2")
comparison_field_kingdom_deg(df.NJ_all_deg, "Fungi", "NJ1", "NJ2")

comparison_field_kingdom_deg(df.CJ_18_all_deg, "Bacteria", "CJ1_18", "CJ2_18")
comparison_field_kingdom_deg(df.CJ_18_all_deg, "Archaea", "CJ1_18", "CJ2_18")
comparison_field_kingdom_deg(df.CJ_18_all_deg, "Fungi", "CJ1_18", "CJ2_18")


##Betweenness
comparison_field_kingdom_between(df.CJ_all_betweenness, "Bacteria", "CJ1", "CJ2")
comparison_field_kingdom_between(df.CJ_all_betweenness, "Archaea", "CJ1", "CJ2")
comparison_field_kingdom_between(df.CJ_all_betweenness, "Fungi", "CJ1", "CJ2")

comparison_field_kingdom_between(df.MY_all_betweenness, "Bacteria", "MY1", "MY2")
comparison_field_kingdom_between(df.MY_all_betweenness, "Archaea", "MY1", "MY2")
comparison_field_kingdom_between(df.MY_all_betweenness, "Fungi", "MY1", "MY2")

comparison_field_kingdom_between(df.YS_all_betweenness, "Bacteria", "YS1", "YS2")
comparison_field_kingdom_between(df.YS_all_betweenness, "Archaea", "YS1", "YS2")
comparison_field_kingdom_between(df.YS_all_betweenness, "Fungi", "YS1", "YS2")

comparison_field_kingdom_between(df.NJ_all_betweenness, "Bacteria", "NJ1", "NJ2")
comparison_field_kingdom_between(df.NJ_all_betweenness, "Archaea", "NJ1", "NJ2")
comparison_field_kingdom_between(df.NJ_all_betweenness, "Fungi", "NJ1", "NJ2")

comparison_field_kingdom_between(df.CJ_18_all_betweenness, "Bacteria", "CJ1_18", "CJ2_18")
comparison_field_kingdom_between(df.CJ_18_all_betweenness, "Archaea", "CJ1_18", "CJ2_18")
comparison_field_kingdom_between(df.CJ_18_all_betweenness, "Fungi", "CJ1_18", "CJ2_18")


##Closeness
comparison_field_kingdom_close(df.CJ_all_closeness, "Bacteria", "CJ1", "CJ2")
comparison_field_kingdom_close(df.CJ_all_closeness, "Archaea", "CJ1", "CJ2")
comparison_field_kingdom_close(df.CJ_all_closeness, "Fungi", "CJ1", "CJ2")

comparison_field_kingdom_close(df.MY_all_closeness, "Bacteria", "MY1", "MY2")
comparison_field_kingdom_close(df.MY_all_closeness, "Archaea", "MY1", "MY2")
comparison_field_kingdom_close(df.MY_all_closeness, "Fungi", "MY1", "MY2")

comparison_field_kingdom_close(df.YS_all_closeness, "Bacteria", "YS1", "YS2")
comparison_field_kingdom_close(df.YS_all_closeness, "Archaea", "YS1", "YS2")
comparison_field_kingdom_close(df.YS_all_closeness, "Fungi", "YS1", "YS2")

comparison_field_kingdom_close(df.NJ_all_closeness, "Bacteria", "NJ1", "NJ2")
comparison_field_kingdom_close(df.NJ_all_closeness, "Archaea", "NJ1", "NJ2")
comparison_field_kingdom_close(df.NJ_all_closeness, "Fungi", "NJ1", "NJ2")

comparison_field_kingdom_close(df.CJ_18_all_closeness, "Bacteria", "CJ1_18", "CJ2_18")
comparison_field_kingdom_close(df.CJ_18_all_closeness, "Archaea", "CJ1_18", "CJ2_18")
comparison_field_kingdom_close(df.CJ_18_all_closeness, "Fungi", "CJ1_18", "CJ2_18")


###Hub OTUs defined by degree and betweenness centrality
library(ggrepel)

hub_plot <- function(df.deg, df.betw){
  df.deg$OTU_id <- rownames(df.deg)
  df.betw$OTU_id <- rownames(df.betw)

df.deg.betw<-merge(df.deg, df.betw, by = intersect(names(df.deg), names(df.betw)))

n<-1
x<-df.deg.betw$OTU_id[df.deg.betw$Degree >= quantile(df.deg.betw$Degree,prob=1-n/100)]
y<-df.deg.betw$OTU_id[df.deg.betw$Betweenness >= quantile(df.deg.betw$Betweenness,prob=1-n/100)]

df.deg.betw.hub <-subset(df.deg.betw, OTU_id %in% intersect(x,y))

df.deg.betw$kingdom <- factor(df.deg.betw$kingdom, levels = c('Bacteria', 'Archaea', 'Fungi'))

p<-ggplot(df.deg.betw, aes(x=Degree, y=Betweenness)) +
  xlab('\n Degree \n')+
  ylab("Betweenness centrality \n") +
  geom_point(aes(colour = kingdom) ,size=5, alpha=0.7) +
  scale_colour_manual(labels = c('Bacteria','Archaea','Fungi'), values = c("#EE7600","#458B74","#9A32CD"))+
  theme(aspect.ratio = 1)+
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  geom_hline(yintercept=quantile(df.deg.betw$Betweenness,prob=1-n/100)  , color="maroon4", linetype='dotted')+
  geom_vline(xintercept=quantile(df.deg.betw$Degree,prob=1-n/100), color="maroon4", linetype='dotted')+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  geom_text_repel(data=df.deg.betw.hub,aes(label=OTU_id), size=4)

return(p)}

hub_plot(df.CJ1_all_deg, df.CJ1_all_betweenness)
hub_plot(df.CJ2_all_deg, df.CJ2_all_betweenness)

hub_plot(df.CJ1_18_all_deg, df.CJ1_18_all_betweenness)
hub_plot(df.CJ2_18_all_deg, df.CJ2_18_all_betweenness)

hub_plot(df.MY1_all_deg, df.MY1_all_betweenness)
hub_plot(df.MY2_all_deg, df.MY2_all_betweenness)

hub_plot(df.NJ1_all_deg, df.NJ1_all_betweenness)
hub_plot(df.NJ2_all_deg, df.NJ2_all_betweenness)

hub_plot(df.YS1_all_deg, df.YS1_all_betweenness)
hub_plot(df.YS2_all_deg, df.YS2_all_betweenness)


hub_plot_close <- function(df.deg, df.close){
  df.deg$OTU_id <- rownames(df.deg)
  df.close$OTU_id <- rownames(df.close)
  
  df.deg.close<-merge(df.deg, df.close, by = intersect(names(df.deg), names(df.close)))
  
  n<-1
  x<-df.deg.close$OTU_id[df.deg.close$Degree >= quantile(df.deg.close$Degree,prob=1-n/100)]
  y<-df.deg.close$OTU_id[df.deg.close$Closeness >= quantile(df.deg.close$Closeness,prob=1-n/100)]
  
  df.deg.close.hub <-subset(df.deg.close, OTU_id %in% intersect(x,y))
  
  df.deg.close$kingdom <- factor(df.deg.close$kingdom, levels = c('Bacteria', 'Archaea', 'Fungi'))
  
  p<-ggplot(df.deg.close, aes(x=Degree, y=Closeness)) +
    xlab('\n Degree \n')+
    ylab("Closeness centrality \n") +
    geom_point(aes(colour = kingdom) ,size=5, alpha=0.7) +
    scale_colour_manual(labels = c('Bacteria','Archaea','Fungi'), values = c("#EE7600","#458B74","#9A32CD"))+
    theme(aspect.ratio = 1)+
    theme(legend.text=element_text(size=13)) + 
    theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    geom_hline(yintercept=quantile(df.deg.close$Closeness,prob=1-n/100)  , color="maroon4", linetype='dotted')+
    geom_vline(xintercept=quantile(df.deg.close$Degree,prob=1-n/100), color="maroon4", linetype='dotted')+
    theme(legend.position="top") +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
    guides(size=FALSE) +
    theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
    
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+ geom_text_repel(data=df.deg.close.hub,aes(label=OTU_id), size=4)
  
  return(p)}

hub_plot_close(df.CJ1_all_deg, df.CJ1_all_closeness)
hub_plot_close(df.CJ2_all_deg, df.CJ2_all_closeness)

hub_plot_close(df.CJ1_18_all_deg, df.CJ1_18_all_closeness)
hub_plot_close(df.CJ2_18_all_deg, df.CJ2_18_all_closeness)

hub_plot_close(df.MY1_all_deg, df.MY1_all_closeness)
hub_plot_close(df.MY2_all_deg, df.MY2_all_closeness)

hub_plot_close(df.NJ1_all_deg, df.NJ1_all_closeness)
hub_plot_close(df.NJ2_all_deg, df.NJ2_all_closeness)

hub_plot_close(df.YS1_all_deg, df.YS1_all_closeness)
hub_plot_close(df.YS2_all_deg, df.YS2_all_closeness)


# comparison of network properties
deg_comparison_plot <- function(df.deg1, df.deg2){
  df.deg1$OTU_id <- rownames(df.deg1)
  df.deg2$OTU_id <- rownames(df.deg2)
  
  names(df.deg1)[1] <- "Degree1"
  names(df.deg2)[1] <- "Degree2"
  
  df.deg.com<-data.frame(df.deg1) %>% full_join(df.deg2, by="OTU_id")
  
  df.deg.com$Degree1[is.na(df.deg.com$Degree1)] <- 0
  df.deg.com$Degree2[is.na(df.deg.com$Degree2)] <- 0
  
  df.deg.com$kingdom <- ifelse(grepl("^B",df.deg.com$OTU_id),'Bacteria',ifelse(grepl("^A",df.deg.com$OTU_id),'Archaea', 'Fungi'))
  
  
  n<-1
  x<-df.deg.com$OTU_id[df.deg.com$Degree1 >= quantile(df.deg.com$Degree1,prob=1-n/100)]
  y<-df.deg.com$OTU_id[df.deg.com$Degree2 >= quantile(df.deg.com$Degree2,prob=1-n/100)]
  
  df.deg.com.inter <-subset(df.deg.com, OTU_id %in% intersect(x,y))
  
  df.deg.com$kingdom <- factor(df.deg.com$kingdom, levels = c('Bacteria', 'Archaea', 'Fungi'))
  
  p<-ggplot(df.deg.com, aes(x=Degree1, y=Degree2)) +
    xlab('\n Degree \n')+
    ylab("Degree \n") +
    geom_point(aes(colour = kingdom) ,size=5, alpha=0.7) +
    scale_colour_manual(labels = c('Bacteria','Archaea','Fungi'), values = c("#EE7600","#458B74","#9A32CD"))+
    theme(aspect.ratio = 1)+
    theme(legend.text=element_text(size=13)) + 
    theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    geom_hline(yintercept=quantile(df.deg.com$Degree2,prob=1-n/100)  , color="maroon4", linetype='dotted')+
    geom_vline(xintercept=quantile(df.deg.com$Degree1,prob=1-n/100), color="maroon4", linetype='dotted')+
    theme(legend.position="top") +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
    guides(size=FALSE) +
    theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
    
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    geom_text_repel(data=df.deg.com.inter,aes(label=OTU_id), size=4)
  
  return(p)}

deg_comparison_plot(df.CJ1_all_deg, df.CJ2_all_deg)
deg_comparison_plot(df.CJ1_18_all_deg, df.CJ2_18_all_deg)
deg_comparison_plot(df.MY1_all_deg, df.MY2_all_deg)
deg_comparison_plot(df.NJ1_all_deg, df.NJ2_all_deg)
deg_comparison_plot(df.YS1_all_deg, df.YS2_all_deg)


betw_comparison_plot <- function(df.betw1, df.betw2){
  df.betw1$OTU_id <- rownames(df.betw1)
  df.betw2$OTU_id <- rownames(df.betw2)
  
  names(df.betw1)[1] <- "Betweenness1"
  names(df.betw2)[1] <- "Betweenness2"
  
  df.betw.com<-data.frame(df.betw1) %>% full_join(df.betw2, by="OTU_id")
  
  df.betw.com$Betweenness1[is.na(df.betw.com$Betweenness1)] <- 0
  df.betw.com$Betweenness2[is.na(df.betw.com$Betweenness2)] <- 0
  
  df.betw.com$kingdom <- ifelse(grepl("^B",df.betw.com$OTU_id),'Bacteria',ifelse(grepl("^A",df.betw.com$OTU_id),'Archaea', 'Fungi'))
  
  
  n<-1
  x<-df.betw.com$OTU_id[df.betw.com$Betweenness1 >= quantile(df.betw.com$Betweenness1,prob=1-n/100)]
  y<-df.betw.com$OTU_id[df.betw.com$Betweenness2 >= quantile(df.betw.com$Betweenness2,prob=1-n/100)]
  
  df.betw.com.inter <-subset(df.betw.com, OTU_id %in% intersect(x,y))
  
  df.betw.com$kingdom <- factor(df.betw.com$kingdom, levels = c('Bacteria', 'Archaea', 'Fungi'))
  
  p<-ggplot(df.betw.com, aes(x=Betweenness1, y=Betweenness2)) +
    xlab('\n Betweenness \n')+
    ylab("Betweenness \n") +
    geom_point(aes(colour = kingdom) ,size=5, alpha=0.7) +
    scale_colour_manual(labels = c('Bacteria','Archaea','Fungi'), values = c("#EE7600","#458B74","#9A32CD"))+
    theme(aspect.ratio = 1)+
    theme(legend.text=element_text(size=13)) + 
    theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    geom_hline(yintercept=quantile(df.betw.com$Betweenness2,prob=1-n/100)  , color="maroon4", linetype='dotted')+
    geom_vline(xintercept=quantile(df.betw.com$Betweenness1,prob=1-n/100), color="maroon4", linetype='dotted')+
    theme(legend.position="top") +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
    guides(size=FALSE) +
    theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
    
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    geom_text_repel(data=df.betw.com.inter,aes(label=OTU_id), size=4)
  
  return(p)}

betw_comparison_plot(df.CJ1_all_betweenness, df.CJ2_all_betweenness)
betw_comparison_plot(df.CJ1_18_all_betweenness, df.CJ2_18_all_betweenness)
betw_comparison_plot(df.MY1_all_betweenness, df.MY2_all_betweenness)
betw_comparison_plot(df.NJ1_all_betweenness, df.NJ2_all_betweenness)
betw_comparison_plot(df.YS1_all_betweenness, df.YS2_all_betweenness)

dev.off()




## Alpha diversity of subsamples
bac.clean.ss
arch.clean.ss
fun.clean.ss


bac.rarefy.t
arch.rarefy.t


tab_all.bac <- microbiome::alpha(bac.rarefy.t, index = "all")
tab_all.fun <- microbiome::alpha(fun.rarefy.t, index = "all")
tab_all.arch <- microbiome::alpha(arch.rarefy.t, index = "all")

#write.table(tab_all.2018, "Alpha diversity_all.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")

library(microbiome)
ps1.meta <- b.meta.all
ps2.meta <- f.meta.all
ps3.meta <- b.meta.all

## using Observed OTUs
ps1.meta$Observed <- tab_all.bac$observed 
ps1.meta

## using Shannon
ps1.meta$Shannon <- tab_all.bac$diversity_shannon 
ps1.meta$Shannon


## using evenness Simpson
ps1.meta$simp <- tab_all.bac$evenness_simpson 
ps1.meta$simp

ps1.meta.CJ <- subset(ps1.meta, Field %in% c('CJ1', "CJ2"))
ps1.meta.CJ_18<- subset(ps1.meta, Field %in% c('CJ1.18', "CJ2.18"))
ps1.meta.MY<- subset(ps1.meta, Field %in% c('MY1', "MY2"))
ps1.meta.NJ<- subset(ps1.meta, Field %in% c('NJ1', "NJ2"))
ps1.meta.YS<- subset(ps1.meta, Field %in% c('YS1', "YS2"))

#Fungi
## using Observed OTUs
ps2.meta$Observed <- tab_all.fun$observed 
ps2.meta

## using Shannon
ps2.meta$Shannon <- tab_all.fun$diversity_shannon 
ps2.meta$Shannon


## using evenness Simpson
ps2.meta$simp <- tab_all.fun$evenness_simpson 
ps2.meta$simp

ps2.meta.CJ <- subset(ps2.meta, Field %in% c('CJ1', "CJ2"))
ps2.meta.CJ_18<- subset(ps2.meta, Field %in% c('CJ1.18', "CJ2.18"))
ps2.meta.MY<- subset(ps2.meta, Field %in% c('MY1', "MY2"))
ps2.meta.NJ<- subset(ps2.meta, Field %in% c('NJ1', "NJ2"))
ps2.meta.YS<- subset(ps2.meta, Field %in% c('YS1', "YS2"))


#Archaea
#Fungi
## using Observed OTUs
ps3.meta$Observed <- tab_all.arch$observed 
ps3.meta

## using Shannon
ps3.meta$Shannon <- tab_all.arch$diversity_shannon 
ps3.meta$Shannon


## using evenness Simpson
ps3.meta$simp <- tab_all.arch$evenness_simpson 
ps3.meta$simp

ps3.meta.CJ <- subset(ps3.meta, Field %in% c('CJ1', "CJ2"))
ps3.meta.CJ_18<- subset(ps3.meta, Field %in% c('CJ1.18', "CJ2.18"))
ps3.meta.MY<- subset(ps3.meta, Field %in% c('MY1', "MY2"))
ps3.meta.NJ<- subset(ps3.meta, Field %in% c('NJ1', "NJ2"))
ps3.meta.YS<- subset(ps3.meta, Field %in% c('YS1', "YS2"))



diversity_comparison_observed<- function(df.com.t, f1, f2){
x <- subset(df.com.t, Field==f1)$Observed
y <- subset(df.com.t, Field==f2)$Observed
wil<-wilcox.test(x, y, conf.int = TRUE)
p_value<-wil$p.value

p<-ggplot(data=df.com.t, aes(x=Field, y=Observed, fill=Field)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 2)+
  geom_signif(comparisons = list(c(f1, f2)), step_increase = 0.1, map_signif_level=TRUE)+
  xlab(paste('Mann Whitney U test, P-value = ', p_value))+
  ylab("Observed OTUs\n") +
  #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")


return(p)}

diversity_comparison_observed(ps1.meta.CJ, "CJ1", "CJ2")
diversity_comparison_observed(ps1.meta.CJ_18, "CJ1.18", "CJ2.18")
diversity_comparison_observed(ps1.meta.MY, "MY1", "MY2")
diversity_comparison_observed(ps1.meta.NJ, "NJ1", "NJ2")
diversity_comparison_observed(ps1.meta.YS, "YS1", "YS2")

##Shannon
diversity_comparison_shannon<- function(df.com.t, f1, f2){
  x <- subset(df.com.t, Field==f1)$Shannon
  y <- subset(df.com.t, Field==f2)$Shannon
  wil<-wilcox.test(x, y, conf.int = TRUE)
  p_value<-wil$p.value
  
  p<-ggplot(data=df.com.t, aes(x=Field, y=Shannon, fill=Field)) + geom_boxplot() +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 2)+
    geom_signif(comparisons = list(c(f1, f2)), step_increase = 0.1, map_signif_level=TRUE)+
    xlab(paste('Mann Whitney U test, P-value = ', p_value))+
    ylab("Shannon index\n") +
    #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
    theme(aspect.ratio = 2)+
    theme(plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
  
  
  return(p)}

diversity_comparison_shannon(ps1.meta.CJ, "CJ1", "CJ2")
diversity_comparison_shannon(ps1.meta.CJ_18, "CJ1.18", "CJ2.18")
diversity_comparison_shannon(ps1.meta.MY, "MY1", "MY2")
diversity_comparison_shannon(ps1.meta.NJ, "NJ1", "NJ2")
diversity_comparison_shannon(ps1.meta.YS, "YS1", "YS2")

diversity_comparison_shannon(ps2.meta.CJ, "CJ1", "CJ2")
diversity_comparison_shannon(ps2.meta.CJ_18, "CJ1.18", "CJ2.18")
diversity_comparison_shannon(ps2.meta.MY, "MY1", "MY2")
diversity_comparison_shannon(ps2.meta.NJ, "NJ1", "NJ2")
diversity_comparison_shannon(ps2.meta.YS, "YS1", "YS2")

diversity_comparison_shannon(ps3.meta.CJ, "CJ1", "CJ2")
diversity_comparison_shannon(ps3.meta.CJ_18, "CJ1.18", "CJ2.18")
diversity_comparison_shannon(ps3.meta.MY, "MY1", "MY2")
diversity_comparison_shannon(ps3.meta.NJ, "NJ1", "NJ2")
diversity_comparison_shannon(ps3.meta.YS, "YS1", "YS2")



## Taxa contributing to higher degree in non-fertilized condition
df.CJ_all_deg



draw_edge_proportion <- function(filename){
  edge.mer <- read.table(file = filename, sep = '\t', header = TRUE)
  head(edge.mer)
  
  edge.mer$Association <- ifelse(grepl("^B",edge.mer$from) & grepl("^F",edge.mer$to) | grepl("^F",edge.mer$from) & grepl("^B",edge.mer$to),'BF',
                                 ifelse(grepl("^A",edge.mer$from) & grepl("^F",edge.mer$to) | grepl("^F",edge.mer$from) & grepl("^A",edge.mer$to),'AF',
                                        ifelse(grepl("^B",edge.mer$from) & grepl("^A",edge.mer$to) | grepl("^A",edge.mer$from) & grepl("^B",edge.mer$to),'AB',
                                               ifelse(grepl("^B",edge.mer$from) & grepl("^B",edge.mer$to),'B', 
                                                      ifelse(grepl("^F",edge.mer$from) & grepl("^F",edge.mer$to),'F','A')))))
  edge.mer
  
  edge.mer$PN <- ifelse(edge.mer$cor > 0,'positive','negative')
  edge.mer$count <- 1
  edge.mer %>% arrange(desc(cor))
  edge.mer %>% arrange(desc(Association))
  
  df.edge.mer <- edge.mer %>% group_by(Association,PN) %>% summarise(count=sum(count))
  
  df.edge.mer$Proportion <- df.edge.mer$count / sum(df.edge.mer$count)
  df.edge.mer
  print(df.edge.mer %>% filter(Association == 'B') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'F') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'A') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'BF') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'AB') %>% mutate(ratio = Proportion / sum(Proportion)))
  print(df.edge.mer %>% filter(Association == 'AF') %>% mutate(ratio = Proportion / sum(Proportion)))
  
  df.edge.mer$Association <- factor(df.edge.mer$Association, levels = c('B', 'A', 'F', 'AB', 'AF', 'BF'))
  
  p.edge.mer <- ggplot(df.edge.mer, aes(x=Association, y = Proportion, fill = PN)) + 
    geom_bar(stat="identity", width = 0.8, position = 'stack', colour="black") +
    #scale_fill_discrete() +
    scale_fill_manual(values = c('#333333','#99CCFF')) +
    
    xlab('')+
    ylab("Proportion of edges \n") +
    #ggtitle("Phylum Community Composition by Sample \n") +
    ## adjust positions
    guides(fill = guide_legend(ncol = 2,reverse = T))+
    theme(legend.position="bottom",legend.title=element_blank()) +
    theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
    theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
    theme(axis.text.x = element_text(size=15, face='bold',color='black'))+
    theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
    scale_y_continuous(breaks=seq(0,1,0.1))+
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(),plot.background=element_blank())
  
  return(p.edge.mer)
}

draw_edge_proportion('sig_correlation_CJ1_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_CJ2_0.6_p0.05.tsv')

draw_edge_proportion('sig_correlation_CJ1_18_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_CJ2_18_0.6_p0.05.tsv')

draw_edge_proportion('sig_correlation_MY1_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_MY2_0.6_p0.05.tsv')

draw_edge_proportion('sig_correlation_NJ1_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_NJ2_0.6_p0.05.tsv')

draw_edge_proportion('sig_correlation_YS1_0.6_p0.05.tsv')
draw_edge_proportion('sig_correlation_YS2_0.6_p0.05.tsv')



CJ1_edge <- read.table('sig_correlation_CJ1_0.6_p0.05.tsv', sep = '\t', header = TRUE)
CJ2_edge <- read.table('sig_correlation_CJ2_0.6_p0.05.tsv', sep = '\t', header = TRUE)
CJ1_18_edge <- read.table('sig_correlation_CJ1_18_0.6_p0.05.tsv', sep = '\t', header = TRUE)
CJ2_18_edge <- read.table('sig_correlation_CJ2_18_0.6_p0.05.tsv', sep = '\t', header = TRUE)
MY1_edge <- read.table('sig_correlation_MY1_0.6_p0.05.tsv', sep = '\t', header = TRUE)
MY2_edge <- read.table('sig_correlation_MY2_0.6_p0.05.tsv', sep = '\t', header = TRUE)
NJ1_edge <- read.table('sig_correlation_NJ1_0.6_p0.05.tsv', sep = '\t', header = TRUE)
NJ2_edge <- read.table('sig_correlation_NJ2_0.6_p0.05.tsv', sep = '\t', header = TRUE)
YS1_edge <- read.table('sig_correlation_YS1_0.6_p0.05.tsv', sep = '\t', header = TRUE)
YS2_edge <- read.table('sig_correlation_YS2_0.6_p0.05.tsv', sep = '\t', header = TRUE)

fun.list$Kingdom <- "Fungi"

otu.list <- rbind(phy.list, arch.list, fun.list)

otu_id.list.2 <- otu.list %>% select("OTU", "OTU_id")

all.list<-otu_id.list.2
all.list <- data.frame(all.list, stringsAsFactors = FALSE)
all.list[is.na(all.list)] <- "Unassigned"

for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ1_edge$from |i %in% CJ1_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_CJ1"] <- as.numeric(length(CJ1_edge$from[which(CJ1_edge$from == i)])+length(CJ1_edge$to[which(CJ1_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_CJ1"]<- 0}
}

head(all.list)

##CJ2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ2_edge$from |i %in% CJ2_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_CJ2"] <- as.numeric(length(CJ2_edge$from[which(CJ2_edge$from == i)])+length(CJ2_edge$to[which(CJ2_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_CJ2"]<- 0}
}

##CJ1_18
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ1_18_edge$from |i %in% CJ1_18_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_CJ1_18"] <- as.numeric(length(CJ1_18_edge$from[which(CJ1_18_edge$from == i)])+length(CJ1_18_edge$to[which(CJ1_18_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_CJ1_18"]<- 0}
}

##CJ2_18
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ2_18_edge$from |i %in% CJ2_18_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_CJ2_18"] <- as.numeric(length(CJ2_18_edge$from[which(CJ2_18_edge$from == i)])+length(CJ2_18_edge$to[which(CJ2_18_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_CJ2_18"]<- 0}
}

##MY1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% MY1_edge$from |i %in% MY1_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_MY1"] <- as.numeric(length(MY1_edge$from[which(MY1_edge$from == i)])+length(MY1_edge$to[which(MY1_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_MY1"]<- 0}
}

##MY2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% MY2_edge$from |i %in% MY2_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_MY2"] <- as.numeric(length(MY2_edge$from[which(MY2_edge$from == i)])+length(MY2_edge$to[which(MY2_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_MY2"]<- 0}
}

##NJ1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% NJ1_edge$from |i %in% NJ1_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_NJ1"] <- as.numeric(length(NJ1_edge$from[which(NJ1_edge$from == i)])+length(NJ1_edge$to[which(NJ1_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_NJ1"]<- 0}
}

##NJ2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% NJ2_edge$from |i %in% NJ2_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_NJ2"] <- as.numeric(length(NJ2_edge$from[which(NJ2_edge$from == i)])+length(NJ2_edge$to[which(NJ2_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_NJ2"]<- 0}
}

##YS1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% YS1_edge$from |i %in% YS1_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_YS1"] <- as.numeric(length(YS1_edge$from[which(YS1_edge$from == i)])+length(YS1_edge$to[which(YS1_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_YS1"]<- 0}
}

##YS2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% YS2_edge$from |i %in% YS2_edge$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_YS2"] <- as.numeric(length(YS2_edge$from[which(YS2_edge$from == i)])+length(YS2_edge$to[which(YS2_edge$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_YS2"]<- 0}
}



head(all.list)
all.list.copy <- all.list
all.list.copy.bac <- subset(all.list.copy, Kingdom == "Bacteria")
all.list.copy.arch <- subset(all.list.copy, Kingdom == "Archaea")
all.list.copy.fun <- subset(all.list.copy, Kingdom == "Fungi")

all.list.CJ1 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_CJ1 = sum(Edge_count_CJ1))
all.list.CJ2 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_CJ2 = sum(Edge_count_CJ2))

all.list.bac.order <- merge(all.list.CJ1, all.list.CJ2, by = "Order")

write.xlsx(all.list.bac.order, "CJ_edge count_order.xlsx")



all.list.CJ1_18 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_CJ1_18 = sum(Edge_count_CJ1_18))
all.list.CJ2_18 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_CJ2_18 = sum(Edge_count_CJ2_18))

all.list.bac.order <- merge(all.list.CJ1_18, all.list.CJ2_18, by = "Order")

write.xlsx(all.list.bac.order, "CJ_18_edge count_order.xlsx")



all.list.MY1 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_MY1 = sum(Edge_count_MY1))
all.list.MY2 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_MY2 = sum(Edge_count_MY2))

all.list.bac.order <- merge(all.list.MY1, all.list.MY2, by = "Order")

write.xlsx(all.list.bac.order, "MY_18_edge count_order.xlsx")


all.list.NJ1 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_NJ1 = sum(Edge_count_NJ1))
all.list.NJ2 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_NJ2 = sum(Edge_count_NJ2))

all.list.bac.order <- merge(all.list.NJ1, all.list.NJ2, by = "Order")

write.xlsx(all.list.bac.order, "NJ_18_edge count_order.xlsx")


all.list.YS1 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_YS1 = sum(Edge_count_YS1))
all.list.YS2 <- all.list.copy.bac %>% group_by(Order) %>% summarise(Edge_count_YS2 = sum(Edge_count_YS2))

all.list.bac.order <- merge(all.list.YS1, all.list.YS2, by = "Order")

write.xlsx(all.list.bac.order, "YS_18_edge count_order.xlsx")



##Archaea
all.list.CJ1 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_CJ1 = sum(Edge_count_CJ1))
all.list.CJ2 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_CJ2 = sum(Edge_count_CJ2))

all.list.arch.order <- merge(all.list.CJ1, all.list.CJ2, by = "Order")

write.xlsx(all.list.arch.order, "CJ_edge count_order_arch.xlsx")



all.list.CJ1_18 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_CJ1_18 = sum(Edge_count_CJ1_18))
all.list.CJ2_18 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_CJ2_18 = sum(Edge_count_CJ2_18))

all.list.arch.order <- merge(all.list.CJ1_18, all.list.CJ2_18, by = "Order")

write.xlsx(all.list.arch.order, "CJ_18_edge count_order_arch.xlsx")



all.list.MY1 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_MY1 = sum(Edge_count_MY1))
all.list.MY2 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_MY2 = sum(Edge_count_MY2))

all.list.arch.order <- merge(all.list.MY1, all.list.MY2, by = "Order")

write.xlsx(all.list.arch.order, "MY_18_edge count_order_arch.xlsx")


all.list.NJ1 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_NJ1 = sum(Edge_count_NJ1))
all.list.NJ2 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_NJ2 = sum(Edge_count_NJ2))

all.list.arch.order <- merge(all.list.NJ1, all.list.NJ2, by = "Order")

write.xlsx(all.list.arch.order, "NJ_18_edge count_order_arch.xlsx")


all.list.YS1 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_YS1 = sum(Edge_count_YS1))
all.list.YS2 <- all.list.copy.arch %>% group_by(Order) %>% summarise(Edge_count_YS2 = sum(Edge_count_YS2))

all.list.arch.order <- merge(all.list.YS1, all.list.YS2, by = "Order")

write.xlsx(all.list.arch.order, "YS_18_edge count_order_arch.xlsx")



##Fungi
all.list.CJ1 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_CJ1 = sum(Edge_count_CJ1))
all.list.CJ2 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_CJ2 = sum(Edge_count_CJ2))

all.list.fun.order <- merge(all.list.CJ1, all.list.CJ2, by = "Order")

write.xlsx(all.list.fun.order, "CJ_edge count_order_fun.xlsx")



all.list.CJ1_18 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_CJ1_18 = sum(Edge_count_CJ1_18))
all.list.CJ2_18 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_CJ2_18 = sum(Edge_count_CJ2_18))

all.list.fun.order <- merge(all.list.CJ1_18, all.list.CJ2_18, by = "Order")

write.xlsx(all.list.fun.order, "CJ_18_edge count_order_fun.xlsx")



all.list.MY1 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_MY1 = sum(Edge_count_MY1))
all.list.MY2 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_MY2 = sum(Edge_count_MY2))

all.list.fun.order <- merge(all.list.MY1, all.list.MY2, by = "Order")

write.xlsx(all.list.fun.order, "MY_18_edge count_order_fun.xlsx")


all.list.NJ1 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_NJ1 = sum(Edge_count_NJ1))
all.list.NJ2 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_NJ2 = sum(Edge_count_NJ2))

all.list.fun.order <- merge(all.list.NJ1, all.list.NJ2, by = "Order")

write.xlsx(all.list.fun.order, "NJ_18_edge count_order_fun.xlsx")


all.list.YS1 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_YS1 = sum(Edge_count_YS1))
all.list.YS2 <- all.list.copy.fun %>% group_by(Order) %>% summarise(Edge_count_YS2 = sum(Edge_count_YS2))

all.list.fun.order <- merge(all.list.YS1, all.list.YS2, by = "Order")

write.xlsx(all.list.fun.order, "YS_18_edge count_order_fun.xlsx")

all.list.bac.order %>% arrange(desc(Edge_count_CJ1))

p<-ggplot(data=df, aes(x=dose, y=len)) +
  geom_bar(stat="identity")


## Total
all.list.total.fer <- all.list.copy.fun %>% group_by(Order) %>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ1))
all.list.total.non_fer <- all.list.copy.fun %>% group_by(Order) %>% summarise(Total_non_fer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2))
all.list.fun.total <- merge(all.list.total.fer, all.list.total.non_fer, by = "Order")

write.xlsx(all.list.fun.total, "Total_edge count_order_fun.xlsx")




#### Microbial taxa involved in Postive and negative associations
all.list<-otu_id.list.2
all.list <- data.frame(all.list, stringsAsFactors = FALSE)
all.list[is.na(all.list)] <- "Unassigned"




##subset Postive and negative associations
CJ1_edge.pos <- subset(CJ1_edge, cor > 0)
CJ1_edge.neg <- subset(CJ1_edge, cor < 0)

CJ2_edge.pos <- subset(CJ2_edge, cor > 0)
CJ2_edge.neg <- subset(CJ2_edge, cor < 0)

CJ1_18_edge.pos <- subset(CJ1_18_edge, cor > 0)
CJ1_18_edge.neg <- subset(CJ1_18_edge, cor < 0)

CJ2_18_edge.pos <- subset(CJ2_18_edge, cor > 0)
CJ2_18_edge.neg <- subset(CJ2_18_edge, cor < 0)

MY1_edge.pos <- subset(MY1_edge, cor > 0)
MY1_edge.neg <- subset(MY1_edge, cor < 0)

MY2_edge.pos <- subset(MY2_edge, cor > 0)
MY2_edge.neg <- subset(MY2_edge, cor < 0)

NJ1_edge.pos <- subset(NJ1_edge, cor > 0)
NJ1_edge.neg <- subset(NJ1_edge, cor < 0)

NJ2_edge.pos <- subset(NJ2_edge, cor > 0)
NJ2_edge.neg <- subset(NJ2_edge, cor < 0)

YS1_edge.pos <- subset(YS1_edge, cor > 0)
YS1_edge.neg <- subset(YS1_edge, cor < 0)

YS2_edge.pos <- subset(YS2_edge, cor > 0)
YS2_edge.neg <- subset(YS2_edge, cor < 0)

### Estimate numbers of positive associations
#CJ1

for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ1_edge.pos$from |i %in% CJ1_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_CJ1"] <- as.numeric(length(CJ1_edge.pos$from[which(CJ1_edge.pos$from == i)])+length(CJ1_edge.pos$to[which(CJ1_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_CJ1"]<- 0}
}

head(all.list)

##CJ2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ2_edge.pos$from |i %in% CJ2_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_CJ2"] <- as.numeric(length(CJ2_edge.pos$from[which(CJ2_edge.pos$from == i)])+length(CJ2_edge.pos$to[which(CJ2_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_CJ2"]<- 0}
}

##CJ1_18
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ1_18_edge.pos$from |i %in% CJ1_18_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_CJ1_18"] <- as.numeric(length(CJ1_18_edge.pos$from[which(CJ1_18_edge.pos$from == i)])+length(CJ1_18_edge.pos$to[which(CJ1_18_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_CJ1_18"]<- 0}
}

##CJ2_18
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ2_18_edge.pos$from |i %in% CJ2_18_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_CJ2_18"] <- as.numeric(length(CJ2_18_edge.pos$from[which(CJ2_18_edge.pos$from == i)])+length(CJ2_18_edge.pos$to[which(CJ2_18_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_CJ2_18"]<- 0}
}

##MY1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% MY1_edge.pos$from |i %in% MY1_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_MY1"] <- as.numeric(length(MY1_edge.pos$from[which(MY1_edge.pos$from == i)])+length(MY1_edge.pos$to[which(MY1_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_MY1"]<- 0}
}

##MY2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% MY2_edge.pos$from |i %in% MY2_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_MY2"] <- as.numeric(length(MY2_edge.pos$from[which(MY2_edge.pos$from == i)])+length(MY2_edge.pos$to[which(MY2_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_MY2"]<- 0}
}

##NJ1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% NJ1_edge.pos$from |i %in% NJ1_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_NJ1"] <- as.numeric(length(NJ1_edge.pos$from[which(NJ1_edge.pos$from == i)])+length(NJ1_edge.pos$to[which(NJ1_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_NJ1"]<- 0}
}

##NJ2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% NJ2_edge.pos$from |i %in% NJ2_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_NJ2"] <- as.numeric(length(NJ2_edge.pos$from[which(NJ2_edge.pos$from == i)])+length(NJ2_edge.pos$to[which(NJ2_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_NJ2"]<- 0}
}

##YS1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% YS1_edge.pos$from |i %in% YS1_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_YS1"] <- as.numeric(length(YS1_edge.pos$from[which(YS1_edge.pos$from == i)])+length(YS1_edge.pos$to[which(YS1_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_YS1"]<- 0}
}

##YS2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% YS2_edge.pos$from |i %in% YS2_edge.pos$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_count_YS2"] <- as.numeric(length(YS2_edge.pos$from[which(YS2_edge.pos$from == i)])+length(YS2_edge.pos$to[which(YS2_edge.pos$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_count_YS2"]<- 0}
}



## Estimate negative count
### Estimate numbers of positive associations
#CJ1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ1_edge.neg$from |i %in% CJ1_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_CJ1"] <- as.numeric(length(CJ1_edge.neg$from[which(CJ1_edge.neg$from == i)])+length(CJ1_edge.neg$to[which(CJ1_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_CJ1"]<- 0}
}

head(all.list)

##CJ2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ2_edge.neg$from |i %in% CJ2_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_CJ2"] <- as.numeric(length(CJ2_edge.neg$from[which(CJ2_edge.neg$from == i)])+length(CJ2_edge.neg$to[which(CJ2_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_CJ2"]<- 0}
}

##CJ1_18
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ1_18_edge.neg$from |i %in% CJ1_18_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_CJ1_18"] <- as.numeric(length(CJ1_18_edge.neg$from[which(CJ1_18_edge.neg$from == i)])+length(CJ1_18_edge.neg$to[which(CJ1_18_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_CJ1_18"]<- 0}
}

##CJ2_18
for (i in as.character(all.list$OTU_id))
{
  if (i %in% CJ2_18_edge.neg$from |i %in% CJ2_18_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_CJ2_18"] <- as.numeric(length(CJ2_18_edge.neg$from[which(CJ2_18_edge.neg$from == i)])+length(CJ2_18_edge.neg$to[which(CJ2_18_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_CJ2_18"]<- 0}
}

##MY1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% MY1_edge.neg$from |i %in% MY1_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_MY1"] <- as.numeric(length(MY1_edge.neg$from[which(MY1_edge.neg$from == i)])+length(MY1_edge.neg$to[which(MY1_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_MY1"]<- 0}
}

##MY2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% MY2_edge.neg$from |i %in% MY2_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_MY2"] <- as.numeric(length(MY2_edge.neg$from[which(MY2_edge.neg$from == i)])+length(MY2_edge.neg$to[which(MY2_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_MY2"]<- 0}
}

##NJ1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% NJ1_edge.neg$from |i %in% NJ1_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_NJ1"] <- as.numeric(length(NJ1_edge.neg$from[which(NJ1_edge.neg$from == i)])+length(NJ1_edge.neg$to[which(NJ1_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_NJ1"]<- 0}
}

##NJ2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% NJ2_edge.neg$from |i %in% NJ2_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_NJ2"] <- as.numeric(length(NJ2_edge.neg$from[which(NJ2_edge.neg$from == i)])+length(NJ2_edge.neg$to[which(NJ2_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_NJ2"]<- 0}
}

##YS1
for (i in as.character(all.list$OTU_id))
{
  if (i %in% YS1_edge.neg$from |i %in% YS1_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_YS1"] <- as.numeric(length(YS1_edge.neg$from[which(YS1_edge.neg$from == i)])+length(YS1_edge.neg$to[which(YS1_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_YS1"]<- 0}
}

##YS2
for (i in as.character(all.list$OTU_id))
{
  if (i %in% YS2_edge.neg$from |i %in% YS2_edge.neg$to== TRUE)
  {all.list[all.list$OTU_id==i,"Edge_neg_count_YS2"] <- as.numeric(length(YS2_edge.neg$from[which(YS2_edge.neg$from == i)])+length(YS2_edge.neg$to[which(YS2_edge.neg$to == i)]))}
  
  else
  {all.list[all.list$OTU_id==i,"Edge_neg_count_YS2"]<- 0}
}

head(all.list)

all.list.class<- tidyr::unite(all.list, "Class", Phylum, Class, sep = "_")
head(all.list.class)
all.list.order<- tidyr::unite(all.list.class, "Order", Class, Order, sep = "_")
all.list.family<- tidyr::unite(all.list.order, "Family", Order, Family, sep = "_")
all.list.genus<- tidyr::unite(all.list.family, "Genus", Family, Genus, sep = "_")

head(all.list.genus)


### Kingdom
all.list.bac <- subset(all.list, Kingdom == "Bacteria")
all.list.arch <- subset(all.list, Kingdom == "Archaea")
all.list.fun <- subset(all.list, Kingdom == "Fungi")

all.list.class.bac <- subset(all.list.class, Kingdom == "Bacteria")
all.list.class.arch <- subset(all.list.class, Kingdom == "Archaea")
all.list.class.fun <- subset(all.list.class, Kingdom == "Fungi")

all.list.order.bac <- subset(all.list.order, Kingdom == "Bacteria")
all.list.order.arch <- subset(all.list.order, Kingdom == "Archaea")
all.list.order.fun <- subset(all.list.order, Kingdom == "Fungi")

all.list.family.bac <- subset(all.list.family, Kingdom == "Bacteria")
all.list.family.arch <- subset(all.list.family, Kingdom == "Archaea")
all.list.family.fun <- subset(all.list.family, Kingdom == "Fungi")

all.list.genus.bac <- subset(all.list.genus, Kingdom == "Bacteria")
all.list.genus.arch <- subset(all.list.genus, Kingdom == "Archaea")
all.list.genus.fun <- subset(all.list.genus, Kingdom == "Fungi")

head(all.list.class.bac)

all.list.class.bac.merged<- all.list.class.bac %>% group_by(Class) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1.pos= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1.neg= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.class.bac.merged, 'Bacteria_class_edges.xlsx')
all.list.class.arch.merged<- all.list.class.arch %>% group_by(Class) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1.pos= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1.neg= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.class.arch.merged, 'Archaea_class_edges.xlsx')
all.list.class.fun.merged<- all.list.class.fun %>% group_by(Class) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1.pos= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1.neg= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.class.fun.merged, 'Fungi_class_edges.xlsx')

all.list.order.bac.merged<- all.list.order.bac %>% group_by(Order) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1.pos= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1.neg= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.order.bac.merged, 'Bacteria_order_edges.xlsx')
all.list.order.arch.merged<- all.list.order.arch %>% group_by(Order) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1.pos= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1.neg= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.order.arch.merged, 'Archaea_order_edges.xlsx')
all.list.order.fun.merged<- all.list.order.fun %>% group_by(Order) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1.pos= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1.neg= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.order.fun.merged, 'Fungi_order_edges.xlsx')

all.list.family.bac.merged<- all.list.family.bac %>% group_by(Family) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.family.bac.merged, 'Bacteria_family_edges.xlsx')
all.list.family.arch.merged<- all.list.family.arch %>% group_by(Family) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.family.arch.merged, 'Archaea_family_edges.xlsx')
all.list.family.fun.merged<- all.list.family.fun %>% group_by(Family) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1= sum(Edge_count_NJ1), NJ2.pos = sum(Edge_count_NJ2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1= sum(Edge_neg_count_NJ1), NJ2.neg = sum(Edge_neg_count_NJ2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.family.fun.merged, 'Fungi_family_edges.xlsx')

all.list.genus.bac.merged<- all.list.genus.bac %>% group_by(Genus) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1= sum(Edge_count_YS1), NJ2.pos = sum(Edge_count_YS2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1= sum(Edge_neg_count_YS1), NJ2.neg = sum(Edge_neg_count_YS2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.genus.bac.merged, 'Bacteria_genus_edges.xlsx')
all.list.genus.arch.merged<- all.list.genus.arch %>% group_by(Genus) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1= sum(Edge_count_YS1), NJ2.pos = sum(Edge_count_YS2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1= sum(Edge_neg_count_YS1), NJ2.neg = sum(Edge_neg_count_YS2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.genus.arch.merged, 'Archaea_genus_edges.xlsx')
all.list.genus.fun.merged<- all.list.genus.fun %>% group_by(Genus) %>% summarise(CJ1.pos = sum(Edge_count_CJ1), CJ2.pos = sum(Edge_count_CJ2),CJ1_18.pos = sum(Edge_count_CJ1_18), CJ2_18.pos = sum(Edge_count_CJ2_18), MY1.pos= sum(Edge_count_MY1), MY2.pos = sum(Edge_count_MY2),NJ1= sum(Edge_count_YS1), NJ2.pos = sum(Edge_count_YS2),YS1.pos= sum(Edge_count_YS1), YS2.pos = sum(Edge_count_YS2), CJ1.neg = sum(Edge_neg_count_CJ1), CJ2.neg = sum(Edge_neg_count_CJ2),CJ1_18.neg = sum(Edge_neg_count_CJ1_18), CJ2_18.neg = sum(Edge_neg_count_CJ2_18), MY1.neg= sum(Edge_neg_count_MY1), MY2.neg = sum(Edge_neg_count_MY2),NJ1= sum(Edge_neg_count_YS1), NJ2.neg = sum(Edge_neg_count_YS2),YS1.neg= sum(Edge_neg_count_YS1), YS2.neg = sum(Edge_neg_count_YS2))
write.xlsx(all.list.genus.fun.merged, 'Fungi_genus_edges.xlsx')

## Subset tables based on the theshold
## Total
all.list.copy <- all.list
all.list<-all.list.copy

total.nonfer.bac.order<-all.list.order.bac %>% group_by(Order)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.bac.order<-all.list.order.bac %>% group_by(Order)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.bac.order<-all.list.order.bac %>% group_by(Order)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

length(total.bac.order$Order[total.bac.order$Edge_total >20])

bac.edge.order <- merge(total.nonfer.bac.order, total.fer.bac.order, by = "Order")


total.nonfer.fun.order<-all.list.order.fun %>% group_by(Order)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.fun.order<-all.list.order.fun %>% group_by(Order)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.fun.order<-all.list.order.fun %>% group_by(Order)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

fun.edge.order <- merge(total.nonfer.fun.order, total.fer.fun.order, by = "Order")

total.nonfer.arch.order<-all.list.order.arch %>% group_by(Order)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.arch.order<-all.list.order.arch %>% group_by(Order)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.arch.order<-all.list.order.arch %>% group_by(Order)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

arch.edge.order <- merge(total.nonfer.arch.order, total.fer.arch.order, by = "Order")


#Family
total.nonfer.bac.family<-all.list.family.bac %>% group_by(Family)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.bac.family<-all.list.family.bac %>% group_by(Family)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.bac.family<-all.list.family.bac %>% group_by(Family)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

bac.edge.family <- merge(total.nonfer.bac.family, total.fer.bac.family, by = "Family")


total.nonfer.fun.family<-all.list.family.fun %>% group_by(Family)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.fun.family<-all.list.family.fun %>% group_by(Family)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.fun.family<-all.list.family.fun %>% group_by(Family)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

fun.edge.family <- merge(total.nonfer.fun.family, total.fer.fun.family, by = "Family")

total.nonfer.arch.family<-all.list.family.arch %>% group_by(Family)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.arch.family<-all.list.family.arch %>% group_by(Family)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.arch.family<-all.list.family.arch %>% group_by(Family)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

arch.edge.family <- merge(total.nonfer.arch.family, total.fer.arch.family, by = "Family")



#Genus
total.nonfer.bac.genus<-all.list.genus.bac %>% group_by(Genus)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.bac.genus<-all.list.genus.bac %>% group_by(Genus)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.bac.genus<-all.list.genus.bac %>% group_by(Genus)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

bac.edge.genus <- merge(total.nonfer.bac.genus, total.fer.bac.genus, by = "Genus")


total.nonfer.fun.genus<-all.list.genus.fun %>% group_by(Genus)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.fun.genus<-all.list.genus.fun %>% group_by(Genus)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.fun.genus<-all.list.genus.fun %>% group_by(Genus)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

fun.edge.genus <- merge(total.nonfer.fun.genus, total.fer.fun.genus, by = "Genus")

total.nonfer.arch.genus<-all.list.genus.arch %>% group_by(Genus)%>% summarise(Total_nonfer = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2))
total.fer.arch.genus<-all.list.genus.arch %>% group_by(Genus)%>% summarise(Total_fer = sum(Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))
total.arch.genus<-all.list.genus.arch %>% group_by(Genus)%>% summarise(Edge_total = sum(Edge_count_YS1, Edge_count_CJ1, Edge_count_CJ1_18, Edge_count_MY2, Edge_count_NJ2,Edge_neg_count_YS1, Edge_neg_count_CJ1, Edge_neg_count_CJ1_18, Edge_neg_count_MY2, Edge_neg_count_NJ2,Edge_count_YS2, Edge_count_CJ2, Edge_count_CJ2_18, Edge_count_MY1, Edge_count_NJ2,Edge_neg_count_YS2, Edge_neg_count_CJ2, Edge_neg_count_CJ2_18, Edge_neg_count_MY1, Edge_neg_count_NJ1))

arch.edge.genus <- merge(total.nonfer.arch.genus, total.fer.arch.genus, by = "Genus")

## Correlation table
bac.order.corr<-read.xlsx("Supplementary Table S5.xlsx", 7)
fun.order.corr<-read.xlsx("Supplementary Table S5.xlsx", 8)
arch.order.corr<-read.xlsx("Supplementary Table S5.xlsx", 9)
bac.family.corr<-read.xlsx("Supplementary Table S5.xlsx", 10)
fun.family.corr<-read.xlsx("Supplementary Table S5.xlsx", 11)
arch.family.corr<-read.xlsx("Supplementary Table S5.xlsx", 12)

bac.genus.corr<-read.table("Bac_genus_corr.tsv", sep = '\t', header = T)
fun.genus.corr<-read.table("Fun_genus_corr.tsv", sep = '\t', header = T)
arch.genus.corr<-read.table("Arch_genus_corr.tsv", sep = '\t', header = T)


#Bac
bac.order.corr.edge<-merge(bac.edge.order, bac.order.corr, by = "Order")
head(bac.order.corr.edge)

bac.order.corr.edge.target<-subset(bac.order.corr.edge, Total_nonfer > Total_fer)
write.xlsx(bac.order.corr.edge.target, "bac.order.corr.edge.target.xlsx")

bac.order.corr.edge.fer<-subset(bac.order.corr.edge, Total_nonfer < Total_fer)
write.xlsx(bac.order.corr.edge.fer, "bac.order.corr.edge.target_fer.xlsx")

#Fun
fun.order.corr.edge<-merge(fun.edge.order, fun.order.corr, by = "Order")
head(fun.order.corr.edge)

fun.order.corr.edge.target<-subset(fun.order.corr.edge, Total_nonfer > Total_fer)
write.xlsx(fun.order.corr.edge.target, "fun.order.corr.edge.target.xlsx")

fun.order.corr.edge.fer<-subset(fun.order.corr.edge, Total_nonfer < Total_fer)
write.xlsx(fun.order.corr.edge.fer, "fun.order.corr.edge.target_fer.xlsx")

#Arch
arch.order.corr.edge<-merge(arch.edge.order, arch.order.corr, by = "Order")
head(arch.order.corr.edge)

arch.order.corr.edge.target<-subset(arch.order.corr.edge, Total_nonfer > Total_fer)
write.xlsx(arch.order.corr.edge.target, "arch.order.corr.edge.target.xlsx")

arch.order.corr.edge.fer<-subset(arch.order.corr.edge, Total_nonfer < Total_fer)
write.xlsx(arch.order.corr.edge.fer, "arch.order.corr.edge.target_fer.xlsx")



#### At the Family level
#Bac
bac.family.corr.edge<-merge(bac.edge.family, bac.family.corr, by = "Family")
head(bac.family.corr.edge)

bac.family.corr.edge.target<-subset(bac.family.corr.edge, Total_nonfer > Total_fer)
write.xlsx(bac.family.corr.edge.target, "bac.family.corr.edge.target.xlsx")

bac.family.corr.edge.fer<-subset(bac.family.corr.edge, Total_nonfer < Total_fer)
write.xlsx(bac.family.corr.edge.fer, "bac.family.corr.edge.target_fer.xlsx")

#Fun
fun.family.corr.edge<-merge(fun.edge.family, fun.family.corr, by = "Family")
head(fun.family.corr.edge)

fun.family.corr.edge.target<-subset(fun.family.corr.edge, Total_nonfer > Total_fer)
write.xlsx(fun.family.corr.edge.target, "fun.family.corr.edge.target.xlsx")

fun.family.corr.edge.fer<-subset(fun.family.corr.edge, Total_nonfer < Total_fer)
write.xlsx(fun.family.corr.edge.fer, "fun.family.corr.edge.target_fer.xlsx")

#Arch
arch.family.corr.edge<-merge(arch.edge.family, arch.family.corr, by = "Family")
head(arch.family.corr.edge)

arch.family.corr.edge.target<-subset(arch.family.corr.edge, Total_nonfer > Total_fer)
write.xlsx(arch.family.corr.edge.target, "arch.family.corr.edge.target.xlsx")

arch.family.corr.edge.fer<-subset(arch.family.corr.edge, Total_nonfer < Total_fer)
write.xlsx(arch.family.corr.edge.fer, "arch.family.corr.edge.target_fer.xlsx")



#### At the Genus level
#Bac
bac.genus.corr.edge<-merge(bac.edge.genus, bac.genus.corr, by = "Genus")
head(bac.genus.corr.edge)

bac.genus.corr.edge.target<-subset(bac.genus.corr.edge, Total_nonfer > Total_fer)
write.xlsx(bac.genus.corr.edge.target, "bac.genus.corr.edge.target.xlsx")

bac.genus.corr.edge.fer<-subset(bac.genus.corr.edge, Total_nonfer < Total_fer)
write.xlsx(bac.genus.corr.edge.fer, "bac.genus.corr.edge.target_fer.xlsx")

#Fun
fun.genus.corr.edge<-merge(fun.edge.genus, fun.genus.corr, by = "Genus")
head(fun.genus.corr.edge)

fun.genus.corr.edge.target<-subset(fun.genus.corr.edge, Total_nonfer > Total_fer)
write.xlsx(fun.genus.corr.edge.target, "fun.genus.corr.edge.target.xlsx")

fun.genus.corr.edge.fer<-subset(fun.genus.corr.edge, Total_nonfer < Total_fer)
write.xlsx(fun.genus.corr.edge.fer, "fun.genus.corr.edge.target_fer.xlsx")

#Arch
arch.genus.corr.edge<-merge(arch.edge.genus, arch.genus.corr, by = "Genus")
head(arch.genus.corr.edge)

arch.genus.corr.edge.target<-subset(arch.genus.corr.edge, Total_nonfer > Total_fer)
write.xlsx(arch.genus.corr.edge.target, "arch.genus.corr.edge.target.xlsx")

arch.genus.corr.edge.fer<-subset(arch.genus.corr.edge, Total_nonfer < Total_fer)
write.xlsx(arch.genus.corr.edge.fer, "arch.genus.corr.edge.target_fer.xlsx")

##Correlation plot for putative oligotrophic microbes affecting microbial network
library(stringr)

ps1.meta.all.order.copy <- ps1.meta.all.order
ps2.meta.all.genus
ps3.meta.all.genus

names(ps1.meta.all.order)<-str_replace_all(names(ps1.meta.all.order), c(" " = "_" ))


ggplot(ps1.meta.all.order, aes(x=TN, y=`Cyanobacteria_Oxyphotobacteria_Nostocales`)) +
  xlab('TN')+
  ylab("RA") +
  geom_point(size=3, alpha=0.7) +
  #scale_colour_manual(labels = c('high','medium','low'), values = c("#CC9900", "#0066CC", '#336633'))+
  #scale_shape_manual(values=c(16,15,17))+
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #geom_hline(yintercept=quantile(df.norm.degree$y2018_all_deg,prob=1-1/100) , color="maroon4", linetype='dotted')+
  #geom_vline(xintercept=quantile(df.norm.degree$y2017_all_deg,prob=1-1/100), color="maroon4", linetype='dotted')+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred")


ggplot(ps2.meta.all.order, aes(x=SOM, y=Tremellomycetes_Filobasidiales)) +
  xlab('SOM')+
  ylab("RA") +
  geom_point(size=3, alpha=0.7) +
  #scale_colour_manual(labels = c('high','medium','low'), values = c("#CC9900", "#0066CC", '#336633'))+
  #scale_shape_manual(values=c(16,15,17))+
  theme(aspect.ratio = 1)+
  # ggtitle("Volcano Plot \n") +
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  #geom_hline(yintercept=quantile(df.norm.degree$y2018_all_deg,prob=1-1/100) , color="maroon4", linetype='dotted')+
  #geom_vline(xintercept=quantile(df.norm.degree$y2017_all_deg,prob=1-1/100), color="maroon4", linetype='dotted')+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  #scale_x_continuous(breaks=seq(0,1,0.2))+
  #scale_y_continuous(breaks=seq(-20,0,-5))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=12, face='bold',color='black'))+
  
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred")


ps1.meta.all.order$Planctomycetacia_Gemmatales

## correlation with linear regression
lmMod <- lm(Gammaproteobacteria_Methylococcales ~ SOM, data=ps1.meta.all.order)  # build the model
summary (lmMod) #R2: 0.1609, p-value: 2.228e-11

lmMod <- lm(NC10_Rokubacteriales ~ SOM, data=ps1.meta.all.order)  # build the model
summary (lmMod) #R2: 0.07755, p-value: 4.278e-06

lmMod <- lm(Planctomycetacia_Planctomycetales ~ SOM, data=ps1.meta.all.order)  # build the model
summary (lmMod) #R2: 0.1786, p-value: 1.47e-12

lmMod <- lm(`Proteobacteria_Gammaproteobacteria_Methylococcales` ~ SOM, data=ps1.meta.all.order)  # build the model
summary (lmMod) #R2:0.1609, p-value: 2.228e-11

lmMod <- lm(`Proteobacteria_Deltaproteobacteria_Myxococcales` ~ TN, data=ps1.meta.all.order)  # build the model
summary (lmMod) #R2: 0.1941 , p-value: 1.315e-13

lmMod <- lm(`Cyanobacteria_Oxyphotobacteria_Nostocales` ~ TN, data=ps1.meta.all.order)  # build the model
summary (lmMod) #R2:0.06879 , p-value: 1.469e-05

lmMod <- lm(Planctomycetacia_Gemmatales ~ SOM, data=ps1.meta.all.order)  # build the model
summary (lmMod) #R2:0.1151 , p-value: 1.977e-08

lmMod <- lm(Archaeosporomycetes_Archaeosporales ~ TN, data=ps2.meta.all.order)  # build the model
summary (lmMod) #R2:0.03924 , p-value: 0.006601

lmMod <- lm(Tremellomycetes_Filobasidiales ~ SOM, data=ps2.meta.all.order)  # build the model
summary (lmMod) #R2:0.1282  , p-value: 2.918e-09

lmMod <- lm(`Ascomycota_Eurotiomycetes_Eurotiales` ~ SOM, data=ps2.meta.all.order)  # build the model
summary (lmMod) #R2:0.04421 , p-value: 0.004173


lmMod <- lm(`Euryarchaeota_Methanomicrobia_Methanosarcinales` ~ SOM, data=ps3.meta.all.order)  # build the model
summary (lmMod) #R2:0.04082  , p-value: 0.0007366

lmMod <- lm(`Euryarchaeota_Thermoplasmata_Methanomassiliicoccales` ~ SOM, data=ps3.meta.all.order)  # build the model
summary (lmMod) #R2:0.07211  , p-value: 9.207e-06



## previously identified rOTU and network properties
head(rotu.final.list)

rOTU.list <- rotu.final.list%>% select("OTU", "rOTU")

head(rOTU.list)


assign_rOTU<-function(df.deg, rOTU.list.form){
  
  df.deg$rOTU <- 0

for (i in rownames(df.deg))
{
  if (i %in% as.character(rOTU.list.form$OTU) == TRUE)
  {df.deg$rOTU[which(rownames(df.deg)==i)] <- as.character(rOTU.list.form$rOTU[which(rOTU.list.form$OTU == i)])}
  
    else
  {df.deg$rOTU[which(rownames(df.deg)==i)] <- "Non-rOTU"}
}
return(df.deg)
}

df.CJ1_all_deg.rOTU<-assign_rOTU(df.CJ1_all_deg, rOTU.list)
df.CJ2_all_deg.rOTU<-assign_rOTU(df.CJ2_all_deg, rOTU.list)
df.CJ1_18_all_deg.rOTU<-assign_rOTU(df.CJ1_18_all_deg, rOTU.list)
df.CJ2_18_all_deg.rOTU<-assign_rOTU(df.CJ2_18_all_deg, rOTU.list)
df.MY1_all_deg.rOTU<-assign_rOTU(df.MY1_all_deg, rOTU.list)
df.MY2_all_deg.rOTU<-assign_rOTU(df.MY2_all_deg, rOTU.list)
df.NJ1_all_deg.rOTU<-assign_rOTU(df.NJ1_all_deg, rOTU.list)
df.NJ2_all_deg.rOTU<-assign_rOTU(df.NJ2_all_deg, rOTU.list)
df.YS1_all_deg.rOTU<-assign_rOTU(df.YS1_all_deg, rOTU.list)
df.YS2_all_deg.rOTU<-assign_rOTU(df.YS2_all_deg, rOTU.list)

## number
length(rownames(df.CJ1_all_deg.rOTU)[which(df.CJ1_all_deg.rOTU$rOTU == "CF")]) #32
length(rownames(df.CJ1_all_deg.rOTU)[which(df.CJ1_all_deg.rOTU$rOTU == "NF")]) #223
length(rownames(df.CJ1_all_deg.rOTU)[which(df.CJ1_all_deg.rOTU$rOTU == "NP")]) #31

length(rownames(df.CJ2_all_deg.rOTU)[which(df.CJ2_all_deg.rOTU$rOTU == "CF")]) #29
length(rownames(df.CJ2_all_deg.rOTU)[which(df.CJ2_all_deg.rOTU$rOTU == "NF")]) #94
length(rownames(df.CJ2_all_deg.rOTU)[which(df.CJ2_all_deg.rOTU$rOTU == "NP")]) #45


length(rownames(df.MY1_all_deg.rOTU)[which(df.MY1_all_deg.rOTU$rOTU == "CF")]) #58
length(rownames(df.MY1_all_deg.rOTU)[which(df.MY1_all_deg.rOTU$rOTU == "NF")]) #102
length(rownames(df.MY1_all_deg.rOTU)[which(df.MY1_all_deg.rOTU$rOTU == "NP")]) #29

length(rownames(df.MY2_all_deg.rOTU)[which(df.MY2_all_deg.rOTU$rOTU == "CF")]) #34
length(rownames(df.MY2_all_deg.rOTU)[which(df.MY2_all_deg.rOTU$rOTU == "NF")]) #223
length(rownames(df.MY2_all_deg.rOTU)[which(df.MY2_all_deg.rOTU$rOTU == "NP")]) #23


length(rownames(df.NJ1_all_deg.rOTU)[which(df.NJ1_all_deg.rOTU$rOTU == "CF")]) #46
length(rownames(df.NJ1_all_deg.rOTU)[which(df.NJ1_all_deg.rOTU$rOTU == "NF")]) #134
length(rownames(df.NJ1_all_deg.rOTU)[which(df.NJ1_all_deg.rOTU$rOTU == "NP")]) #34

length(rownames(df.NJ2_all_deg.rOTU)[which(df.NJ2_all_deg.rOTU$rOTU == "CF")]) #29
length(rownames(df.NJ2_all_deg.rOTU)[which(df.NJ2_all_deg.rOTU$rOTU == "NF")]) #245
length(rownames(df.NJ2_all_deg.rOTU)[which(df.NJ2_all_deg.rOTU$rOTU == "NP")]) #30

length(rownames(df.YS1_all_deg.rOTU)[which(df.YS1_all_deg.rOTU$rOTU == "CF")]) #24
length(rownames(df.YS1_all_deg.rOTU)[which(df.YS1_all_deg.rOTU$rOTU == "NF")]) #314
length(rownames(df.YS1_all_deg.rOTU)[which(df.YS1_all_deg.rOTU$rOTU == "NP")]) #32

length(rownames(df.YS2_all_deg.rOTU)[which(df.YS2_all_deg.rOTU$rOTU == "CF")]) #42
length(rownames(df.YS2_all_deg.rOTU)[which(df.YS2_all_deg.rOTU$rOTU == "NF")]) #68
length(rownames(df.YS2_all_deg.rOTU)[which(df.YS2_all_deg.rOTU$rOTU == "NP")]) #65

## contribution of rOTUs in network connectivity
deg_comparison_rOTU <- function(df.deg){
  test <- aggregate(df.deg$Degree, by = list(df.deg$rOTU), max)
  
  colnames(test) <- c("Group", "maxdegree")
  
  kw<-kruskal.test(Degree ~ rOTU, data = df.deg)
  
  #library(FSA)
  DT = dunnTest(Degree ~ rOTU,
                data=df.deg,
                method="bh")
  PT = DT$res
  #library(rcompanion)
  dunn<-cldList(P.adj ~ Comparison,
                data = PT,
                threshold = 0.05)
  hsd1 <- merge(dunn,test, by = 'Group')
  names(hsd1)[1] <- "rOTU"
  
  p<-ggplot(data=df.deg, aes(x=rOTU, y=Degree, fill=rOTU)) + geom_boxplot() +
    stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = .75, linetype = "dashed") + 
    geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
    geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
    labs(x="", y ="Degree") +
    #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
    theme(aspect.ratio = 2)+
    theme(plot.background = element_blank()
          ,panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
  
  return(p)
}

deg_comparison_rOTU(df.CJ1_all_deg.rOTU)


p<-ggplot(data=df.CJ2_all_deg.rOTU, aes(x=rOTU, y=Degree, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p


p<-ggplot(data=df.MY1_all_deg.rOTU, aes(x=rOTU, y=Degree, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p

p<-ggplot(data=df.MY2_all_deg.rOTU, aes(x=rOTU, y=Degree, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p

p<-ggplot(data=df.NJ1_all_deg.rOTU, aes(x=rOTU, y=Degree, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p


p<-ggplot(data=df.YS2_all_deg.rOTU, aes(x=rOTU, y=Degree, fill=rOTU)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p


###histogram
df.CJ1_all_deg


ggplot(data = df.MY1_all_deg.rOTU, aes(x = Degree, fill = kingdom)) + geom_area(stat = "bin")

ggplot(data = df.MY2_all_deg.rOTU, aes(x = Degree, fill = kingdom)) + geom_area(stat = "bin")

ggplot(data = df.YS2_all_deg.rOTU, aes(x = Degree, fill = kingdom)) + geom_area(stat = "bin")


ggplot(data = subset(df.YS2_all_deg.rOTU, rOTU != "Non-rOTU"), aes(x = Degree, fill = rOTU)) + 
  geom_histogram(colour = 'white') 


ggplot(data = subset(df.YS1_all_deg.rOTU, rOTU != "Non-rOTU"), aes(x = Degree, fill = rOTU)) + 
  geom_histogram(colour = 'white') 

ggplot(data = subset(df.NJ2_all_deg.rOTU, rOTU != "Non-rOTU"), aes(x = Degree, fill = rOTU)) + 
  geom_histogram(colour = 'white') 

ggplot(data = subset(df.NJ1_all_deg.rOTU, rOTU != "Non-rOTU"), aes(x = Degree, fill = rOTU)) + 
  geom_histogram(colour = 'white') 


ggplot(data = subset(df.MY2_all_deg.rOTU, rOTU != "Non-rOTU"), aes(x = Degree, fill = rOTU)) + 
  geom_histogram(colour = 'white') 

ggplot(data = subset(df.MY1_all_deg.rOTU, rOTU != "Non-rOTU"), aes(x = Degree, fill = rOTU)) + 
  geom_histogram(colour = 'white') 


#Kingdom histogram (area plot)

ggplot(data = subset(df.CJ1_all_deg.rOTU, kingdom == "Archaea"), aes(x = Degree)) + geom_area(stat = "bin")
ggplot(data = subset(df.CJ1_all_deg.rOTU, kingdom == "Fungi"), aes(x = Degree)) + geom_area(stat = "bin")
ggplot(data = subset(df.CJ1_all_deg.rOTU, kingdom == "Bacteria"), aes(x = Degree)) + geom_area(stat = "bin")


ggplot(data = subset(df.CJ2_all_deg.rOTU, kingdom == "Archaea"), aes(x = Degree)) + geom_area(stat = "bin")
ggplot(data = subset(df.CJ2_all_deg.rOTU, kingdom == "Fungi"), aes(x = Degree)) + geom_area(stat = "bin")
ggplot(data = subset(df.CJ2_all_deg.rOTU, kingdom == "Bacteria"), aes(x = Degree)) + geom_area(stat = "bin")


ggplot(data = subset(df.MY2_all_deg.rOTU, kingdom == "Archaea"), aes(x = Degree)) + geom_area(stat = "bin")
ggplot(data = subset(df.MY2_all_deg.rOTU, kingdom == "Fungi"), aes(x = Degree)) + geom_area(stat = "bin")
ggplot(data = subset(df.MY2_all_deg.rOTU, kingdom == "Bacteria"), aes(x = Degree)) + geom_area(stat = "bin")


ggplot(data = subset(df.MY1_all_deg.rOTU, rOTU != "Non-rOTU"), aes(x = Degree, fill = rOTU)) + 
  geom_area(stat = "bin", colour = "white")
ggplot(data = df.MY_all_deg, aes(x = Degree, fill = category)) + 
  geom_area(stat = "bin", colour = "white")


ggplot(data = df.CJ_all_deg, aes(x = Degree, fill = category)) + 
  geom_area(stat = "bin", colour = "white")


ggplot(data = df.YS_all_deg, aes(x = Degree, fill = category)) + 
  geom_area(stat = "bin", colour = "white")

ggplot(data = df.MY_all_deg, aes(x = Degree, fill = category)) + 
  geom_area(stat = "bin", colour = "white")

ggplot(data = df.NJ_all_deg, aes(x = Degree, fill = category)) + 
  geom_area(stat = "bin", colour = "white")



## 
df.all_deg<-read.table("df.all_deg.tsv", sep= '\t', header =T)
df.all_betweenness<-read.table("df.all_betweenness.tsv", sep= '\t', header =T)

df.all_closeness<-read.table("df.all_closeness.tsv", sep= '\t', header =T)

head(df.all_deg)

p<-ggplot(data=df.all_deg, aes(x=category, y=Degree, fill=category)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',shape=1, alpha=.3, size = 0.5)+
  #geom_text(data=hsd1,aes(x=rOTU,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#EE7600","#458B74","#9A32CD"), labels=c("Bacteria", "Archaea", "Fungi")) + 
  theme(aspect.ratio = 2)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p


