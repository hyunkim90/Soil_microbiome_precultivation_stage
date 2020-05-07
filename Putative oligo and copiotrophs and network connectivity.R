## rOTU and oligo and copiotrophs
rotu.final.list.otu_id <- rotu.final.list %>% select("OTU", "rOTU")
names(rotu.final.list.otu_id)[1] <- "OTU_id"

rotu.final.list.otu_id<- merge(OTU_id.list, rotu.final.list.otu_id, by = "OTU_id")
head(rotu.final.list.otu_id)
head(tab_ra.b.SOM.TN.oligo)

tab_ra.oligo <- rbind(tab_ra.b.SOM.TN.oligo,tab_ra.a.SOM.TN.oligo, tab_ra.f.SOM.TN.oligo)


tab_ra.oligo$OTU <- rownames(tab_ra.oligo)

rotu.oligo <- merge(tab_ra.oligo,rotu.final.list.otu_id, by='OTU')
head(rotu.oligo)
length(rotu.oligo$OTU_id[rotu.oligo$rOTU == "CF"]) #49
sum(rotu.oligo$RA[rotu.oligo$rOTU == "CF"]) #0.04193766
length(rotu.oligo$OTU_id[rotu.oligo$rOTU == "NF"]) #364
sum(rotu.oligo$RA[rotu.oligo$rOTU == "NF"]) # 0.06688968
length(rotu.oligo$OTU_id[rotu.oligo$rOTU == "NP"]) #22
sum(rotu.oligo$RA[rotu.oligo$rOTU == "NP"]) #0.1362102
length(rotu.oligo$OTU_id[rotu.oligo$rOTU == "CF-NP"]) #8
length(rotu.oligo$OTU_id[rotu.oligo$rOTU == "CF-NF"]) #16
length(rotu.oligo$OTU_id[rotu.oligo$rOTU == "NF-NP"]) #10
length(rotu.oligo$OTU_id[rotu.oligo$rOTU == "Non-rOTU"]) # 1246
sum(rotu.oligo$RA[rotu.oligo$rOTU == "Non-rOTU"]) #0.3804444

tab_ra.copi <- rbind(tab_ra.b.SOM.TN.copi,tab_ra.a.SOM.TN.copi, tab_ra.f.SOM.TN.copi)


tab_ra.copi$OTU <- rownames(tab_ra.copi)

rotu.copi <- merge(tab_ra.copi,rotu.final.list.otu_id, by='OTU')
head(rotu.copi)
length(rotu.copi$OTU_id[rotu.copi$rOTU == "CF"]) #50
sum(rotu.copi$RA[rotu.copi$rOTU == "CF"]) #0.02440122
length(rotu.copi$OTU_id[rotu.copi$rOTU == "NF"]) #33
sum(rotu.copi$RA[rotu.copi$rOTU == "NF"]) #0.126174
length(rotu.copi$OTU_id[rotu.copi$rOTU == "NP"]) #77
sum(rotu.copi$RA[rotu.copi$rOTU == "NP"]) #0.04098322
length(rotu.copi$OTU_id[rotu.copi$rOTU == "Non-rOTU"]) #1432
sum(rotu.copi$RA[rotu.copi$rOTU == "Non-rOTU"]) #0.565698
length(rotu.copi$OTU_id[rotu.copi$rOTU == "CF-NP"]) #42
length(rotu.copi$OTU_id[rotu.copi$rOTU == "CF-NF"]) #6
length(rotu.copi$OTU_id[rotu.copi$rOTU == "NF-NP"]) #0



### Network structure and oligo/copiotrophs - fertilized and non-fertilized condition
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

df.all_deg <- rbind(df.CJ_all_deg, df.CJ_18_all_deg, df.MY_all_deg, df.NJ_all_deg, df.YS_all_deg)
x <- subset(df.all_deg, category=='non-fer')$Degree
mean(subset(df.all_deg, category=='non-fer')$Degree)
y <- subset(df.all_deg, category=='fer')$Degree
mean(subset(df.all_deg, category=='fer')$Degree)
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value


for (i in as.character(rownames(df.all_deg)))
{
  if (i %in% rotu.oligo$OTU_id == TRUE)
  {df.all_deg[rownames(df.all_deg)==i,"Trophic"] <- "Oligo"}
  
  else if (i %in% rotu.copi$OTU_id == TRUE)
  {df.all_deg[rownames(df.all_deg)==i,"Trophic"] <- "Copio"} 
  
  else
  {df.all_deg[rownames(df.all_deg)==i,"Trophic"] <- "Not_predicted"}
}


head(df.all_deg)



df.deg.trophic.CJ1<- subset(df.all_deg, field == "CJ1")
df.deg.trophic.CJ1.18<- subset(df.all_deg, field == "CJ1_18")
df.deg.trophic.CJ2<- subset(df.all_deg, field == "CJ2")
df.deg.trophic.CJ2.18<- subset(df.all_deg, field == "CJ2_18")
df.deg.trophic.MY1<- subset(df.all_deg, field == "MY1")
df.deg.trophic.MY2<- subset(df.all_deg, field == "MY2")
df.deg.trophic.NJ1<- subset(df.all_deg, field == "NJ1")
df.deg.trophic.NJ2<- subset(df.all_deg, field == "NJ2")
df.deg.trophic.YS1<- subset(df.all_deg, field == "YS1")
df.deg.trophic.YS2<- subset(df.all_deg, field == "YS2")


x <- df.deg.trophic.CJ1$Degree[df.deg.trophic.CJ1$Trophic=="Oligo"]
mean(df.deg.trophic.CJ1$Degree[df.deg.trophic.CJ1$Trophic=="Oligo"])
y <- df.deg.trophic.CJ1$Degree[df.deg.trophic.CJ1$Trophic=="Copio"]
mean(df.deg.trophic.CJ1$Degree[df.deg.trophic.CJ1$Trophic=="Copio"])
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.003401826


x <- df.deg.trophic.CJ2$Degree[df.deg.trophic.CJ2$Trophic=="Oligo"]
y <- df.deg.trophic.CJ2$Degree[df.deg.trophic.CJ2$Trophic=="Copio"]
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.4877211

x <- df.deg.trophic.MY1$Degree[df.deg.trophic.MY1$Trophic=="Oligo"]
y <- df.deg.trophic.MY1$Degree[df.deg.trophic.MY1$Trophic=="Copio"]
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.9782145

x <- df.deg.trophic.MY2$Degree[df.deg.trophic.MY2$Trophic=="Oligo"]
mean(df.deg.trophic.MY2$Degree[df.deg.trophic.MY2$Trophic=="Oligo"])
y <- df.deg.trophic.MY2$Degree[df.deg.trophic.MY2$Trophic=="Copio"]
mean(df.deg.trophic.MY2$Degree[df.deg.trophic.MY2$Trophic=="Copio"])
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.9822516


x <- df.deg.trophic.NJ1$Degree[df.deg.trophic.NJ1$Trophic=="Oligo"]
y <- df.deg.trophic.NJ1$Degree[df.deg.trophic.NJ1$Trophic=="Copio"]
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.3693181

x <- df.deg.trophic.NJ2$Degree[df.deg.trophic.NJ2$Trophic=="Oligo"]
y <- df.deg.trophic.NJ2$Degree[df.deg.trophic.NJ2$Trophic=="Copio"]
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.3693181

x <- df.deg.trophic.YS1$Degree[df.deg.trophic.YS1$Trophic=="Oligo"]
y <- df.deg.trophic.YS1$Degree[df.deg.trophic.YS1$Trophic=="Copio"]
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.3693181

x <- df.deg.trophic.YS2$Degree[df.deg.trophic.YS2$Trophic=="Oligo"]
y <- df.deg.trophic.YS2$Degree[df.deg.trophic.YS2$Trophic=="Copio"]
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.3693181

x <- df.deg.trophic$Degree[df.deg.trophic$Trophic=="Copio"]
y <- df.deg.trophic.fer$Degree[df.deg.trophic.fer$Trophic=="Copio"]
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.2526985


## Fertilized vs non-fertilized

x <- df.deg.trophic.CJ1$Degree[df.deg.trophic.CJ1$Trophic=="Oligo"]
mean(df.deg.trophic.CJ1$Degree[df.deg.trophic.CJ1$Trophic=="Oligo"]) #29.87921
rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Oligo"]#505

y <- df.deg.trophic.CJ2$Degree[df.deg.trophic.CJ2$Trophic=="Oligo"]
mean(df.deg.trophic.CJ2$Degree[df.deg.trophic.CJ2$Trophic=="Oligo"]) #33.76562
rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Oligo"]#128

x <- df.deg.trophic.CJ1$Degree[df.deg.trophic.CJ1$Trophic=="Copio"]
mean(df.deg.trophic.CJ1$Degree[df.deg.trophic.CJ1$Trophic=="Copio"])
rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Copio"]#192
y <- df.deg.trophic.CJ2$Degree[df.deg.trophic.CJ2$Trophic=="Copio"]
mean(df.deg.trophic.CJ2$Degree[df.deg.trophic.CJ2$Trophic=="Copio"])
rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Copio"]#82

wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.2535292

intersect(hub.spar.CJ1, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Oligo"])
intersect(hub.spar.CJ1, rownames(df.deg.trophic.CJ1)[df.deg.trophic.CJ1$Trophic=="Copio"])

intersect(hub.spar.CJ2, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Oligo"])
intersect(hub.spar.CJ2, rownames(df.deg.trophic.CJ2)[df.deg.trophic.CJ2$Trophic=="Copio"])

x <- df.deg.trophic.MY1$Degree[df.deg.trophic.MY1$Trophic=="Oligo"]
mean(df.deg.trophic.MY1$Degree[df.deg.trophic.MY1$Trophic=="Oligo"]) #25.36364
rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Oligo"]#11
y <- df.deg.trophic.MY2$Degree[df.deg.trophic.MY2$Trophic=="Oligo"]
mean(df.deg.trophic.MY2$Degree[df.deg.trophic.MY2$Trophic=="Oligo"]) #28.38889
rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Oligo"]#18
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.7527209


x <- df.deg.trophic.MY1$Degree[df.deg.trophic.MY1$Trophic=="Copio"]
mean(df.deg.trophic.MY1$Degree[df.deg.trophic.MY1$Trophic=="Copio"]) #21.6789
rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Copio"]#109
y <- df.deg.trophic.MY2$Degree[df.deg.trophic.MY2$Trophic=="Copio"]
mean(df.deg.trophic.MY2$Degree[df.deg.trophic.MY2$Trophic=="Copio"]) #29.10256
rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Copio"]#156

intersect(rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Copio"],rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Copio"] )

intersect(hub.spar.MY1, rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Oligo"])
intersect(hub.spar.MY1, rownames(df.deg.trophic.MY1)[df.deg.trophic.MY1$Trophic=="Copio"])

intersect(hub.spar.MY2, rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Oligo"])
intersect(hub.spar.MY2, rownames(df.deg.trophic.MY2)[df.deg.trophic.MY2$Trophic=="Copio"])

wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.04417336



x <- df.deg.trophic.NJ1$Degree[df.deg.trophic.NJ1$Trophic=="Oligo"] #15.61111
mean(df.deg.trophic.NJ1$Degree[df.deg.trophic.NJ1$Trophic=="Oligo"])
rownames(df.deg.trophic.NJ1)[df.deg.trophic.NJ1$Trophic=="Oligo"]#18
y <- df.deg.trophic.NJ2$Degree[df.deg.trophic.NJ2$Trophic=="Oligo"]
mean(df.deg.trophic.NJ2$Degree[df.deg.trophic.NJ2$Trophic=="Oligo"]) #34.56522
rownames(df.deg.trophic.NJ2)[df.deg.trophic.NJ2$Trophic=="Oligo"]#23
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.0781547
intersect(rownames(df.deg.trophic.NJ1)[df.deg.trophic.NJ1$Trophic=="Oligo"],rownames(df.deg.trophic.NJ2)[df.deg.trophic.NJ2$Trophic=="Oligo"] )

x <- df.deg.trophic.NJ1$Degree[df.deg.trophic.NJ1$Trophic=="Copio"]
mean(df.deg.trophic.NJ1$Degree[df.deg.trophic.NJ1$Trophic=="Copio"]) #19.35484
rownames(df.deg.trophic.NJ1)[df.deg.trophic.NJ1$Trophic=="Copio"]#31
y <- df.deg.trophic.NJ2$Degree[df.deg.trophic.NJ2$Trophic=="Copio"]
mean(df.deg.trophic.NJ2$Degree[df.deg.trophic.NJ2$Trophic=="Copio"]) #21.7907
rownames(df.deg.trophic.NJ2)[df.deg.trophic.NJ2$Trophic=="Copio"]#86

intersect(rownames(df.deg.trophic.NJ1)[df.deg.trophic.NJ1$Trophic=="Copio"],rownames(df.deg.trophic.NJ2)[df.deg.trophic.NJ2$Trophic=="Copio"] )

wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.7177141


intersect(hub.spar.NJ1, rownames(df.deg.trophic.NJ1)[df.deg.trophic.NJ1$Trophic=="Oligo"])
intersect(hub.spar.NJ1, rownames(df.deg.trophic.NJ1)[df.deg.trophic.NJ1$Trophic=="Copio"])

intersect(hub.spar.NJ2, rownames(df.deg.trophic.NJ2)[df.deg.trophic.NJ2$Trophic=="Oligo"])
intersect(hub.spar.NJ2, rownames(df.deg.trophic.NJ2)[df.deg.trophic.NJ2$Trophic=="Copio"])



x <- df.deg.trophic.YS1$Degree[df.deg.trophic.YS1$Trophic=="Oligo"] #20.18405
mean(df.deg.trophic.YS1$Degree[df.deg.trophic.YS1$Trophic=="Oligo"])
rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Oligo"]#163
y <- df.deg.trophic.YS2$Degree[df.deg.trophic.YS2$Trophic=="Oligo"]
mean(df.deg.trophic.YS2$Degree[df.deg.trophic.YS2$Trophic=="Oligo"]) #30.5
rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Oligo"]#14
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.0781547
intersect(rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Oligo"],rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Oligo"] )

x <- df.deg.trophic.YS1$Degree[df.deg.trophic.YS1$Trophic=="Copio"]
mean(df.deg.trophic.YS1$Degree[df.deg.trophic.YS1$Trophic=="Copio"]) #17.84615
rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Copio"]#13
y <- df.deg.trophic.YS2$Degree[df.deg.trophic.YS2$Trophic=="Copio"]
mean(df.deg.trophic.YS2$Degree[df.deg.trophic.YS2$Trophic=="Copio"]) #26
rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Copio"]#34

intersect(rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Copio"],rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Copio"] )

wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.3407657

intersect(hub.spar.YS1, rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Oligo"])
intersect(hub.spar.YS1, rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Copio"])

intersect(hub.spar.YS2, rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Oligo"])
intersect(hub.spar.YS2, rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Copio"])



x <- df.deg.trophic.CJ1.18$Degree[df.deg.trophic.CJ1.18$Trophic=="Oligo"]
mean(df.deg.trophic.CJ1.18$Degree[df.deg.trophic.CJ1.18$Trophic=="Oligo"]) # 20.51837
rownames(df.deg.trophic.CJ1.18)[df.deg.trophic.CJ1.18$Trophic=="Oligo"]#245
y <- df.deg.trophic.CJ2.18$Degree[df.deg.trophic.CJ2.18$Trophic=="Oligo"]
mean(df.deg.trophic.CJ2.18$Degree[df.deg.trophic.CJ2.18$Trophic=="Oligo"]) #18.35714
rownames(df.deg.trophic.CJ2.18)[df.deg.trophic.CJ2.18$Trophic=="Oligo"]#28
wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value # 0.0781547
intersect(rownames(df.deg.trophic.CJ1.18)[df.deg.trophic.CJ1.18$Trophic=="Oligo"],rownames(df.deg.trophic.CJ2.18)[df.deg.trophic.CJ2.18$Trophic=="Oligo"] )

x <- df.deg.trophic.CJ1.18$Degree[df.deg.trophic.CJ1.18$Trophic=="Copio"]
mean(df.deg.trophic.CJ1.18$Degree[df.deg.trophic.CJ1.18$Trophic=="Copio"]) #20.42105
rownames(df.deg.trophic.CJ1.18)[df.deg.trophic.CJ1.18$Trophic=="Copio"]#38
y <- df.deg.trophic.CJ2.18$Degree[df.deg.trophic.CJ2.18$Trophic=="Copio"]
mean(df.deg.trophic.CJ2.18$Degree[df.deg.trophic.CJ2.18$Trophic=="Copio"]) #19.65517
rownames(df.deg.trophic.CJ2.18)[df.deg.trophic.CJ2.18$Trophic=="Copio"]#29

intersect(rownames(df.deg.trophic.CJ1.18)[df.deg.trophic.CJ1.18$Trophic=="Copio"],rownames(df.deg.trophic.CJ2.18)[df.deg.trophic.CJ2.18$Trophic=="Copio"] )

wil<-wilcox.test(x, y, conf.int = TRUE)
wil$p.value #0.7177141


intersect(hub.spar.CJ1_18, rownames(df.deg.trophic.CJ1.18)[df.deg.trophic.CJ1.18$Trophic=="Oligo"]) #"F335_unidentified"
intersect(hub.spar.CJ1_18, rownames(df.deg.trophic.CJ1.18)[df.deg.trophic.CJ1.18$Trophic=="Copio"])

intersect(hub.spar.CJ2_18, rownames(df.deg.trophic.CJ2.18)[df.deg.trophic.CJ2.18$Trophic=="Oligo"])
intersect(hub.spar.CJ2_18, rownames(df.deg.trophic.CJ2.18)[df.deg.trophic.CJ2.18$Trophic=="Copio"]) #"F386_Scutellinia"



## Degree histogram
df.deg.trophic.CJ1.trophic <- subset(df.deg.trophic.CJ1, Trophic != "Not_predicted")
df.deg.trophic.CJ2.trophic <- subset(df.deg.trophic.CJ2, Trophic != "Not_predicted")

df.deg.trophic.CJ1.18.trophic <- subset(df.deg.trophic.CJ1.18, Trophic != "Not_predicted")
df.deg.trophic.CJ2.18.trophic <- subset(df.deg.trophic.CJ2.18, Trophic != "Not_predicted")

df.deg.trophic.MY1.trophic <- subset(df.deg.trophic.MY1, Trophic != "Not_predicted")
df.deg.trophic.MY2.trophic <- subset(df.deg.trophic.MY2, Trophic != "Not_predicted")

df.deg.trophic.NJ1.trophic <- subset(df.deg.trophic.NJ1, Trophic != "Not_predicted")
df.deg.trophic.NJ2.trophic <- subset(df.deg.trophic.NJ2, Trophic != "Not_predicted")

df.deg.trophic.YS1.trophic <- subset(df.deg.trophic.YS1, Trophic != "Not_predicted")
df.deg.trophic.YS2.trophic <- subset(df.deg.trophic.YS2, Trophic != "Not_predicted")


ggplot(df.deg.trophic.CJ1.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.CJ2.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.CJ1.18.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.CJ2.18.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.MY1.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.MY2.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.NJ1.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.NJ2.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.YS1.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)

ggplot(df.deg.trophic.YS2.trophic, aes(Degree, fill = Trophic))+geom_histogram(position = "stack", binwidth=2)



test <- aggregate(df.deg.trophic$Degree, by = list(df.deg.trophic$Trophic), max)

colnames(test) <- c("Group", "maxdegree")

##Kruskal-Wallis
kw<-kruskal.test(Degree ~ Trophic, data = df.deg.trophic)

#library(FSA)
DT = dunnTest(Degree ~ Trophic,
              data=df.deg.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "Trophic"

p<-ggplot(data=df.deg.trophic, aes(x=Trophic, y=Degree, fill=Trophic)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2)+
  #scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=Trophic,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p




df.deg.trophic<- subset(df.all_deg, category == "non-fer")

test <- aggregate(df.deg.trophic$Degree, by = list(df.deg.trophic$Trophic), max)

colnames(test) <- c("Group", "maxdegree")

df.deg.trophic$Trophic <- as.factor(df.deg.trophic$Trophic)

##Kruskal-Wallis
kw<-kruskal.test(Trophic ~ Degree, data = df.deg.trophic)

#library(FSA)
DT = dunnTest(Degree ~ Trophic,
              data=df.deg.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "Trophic"

p<-ggplot(data=df.deg.trophic, aes(x=Trophic, y=Degree, fill=Trophic)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2)+
  #scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=Trophic,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p



df.deg.trophic.fer<- subset(df.all_deg, category == "fer")

mean(df.deg.trophic$Degree[df.deg.trophic$Trophic=="Copio"]) #31.05319
mean(df.deg.trophic$Degree[df.deg.trophic$Trophic=="Oligo"]) #36.87931



test <- aggregate(df.deg.trophic$Degree, by = list(df.deg.trophic$Trophic), max)

colnames(test) <- c("Group", "maxdegree")

df.deg.trophic$Trophic <- as.factor(df.deg.trophic$Trophic)

##Kruskal-Wallis
kw<-kruskal.test(Trophic ~ Degree, data = df.deg.trophic)

#library(FSA)
DT = dunnTest(Degree ~ Trophic,
              data=df.deg.trophic,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
hsd1 <- merge(dunn,test, by = 'Group')
names(hsd1)[1] <- "Trophic"

p<-ggplot(data=df.deg.trophic, aes(x=Trophic, y=Degree, fill=Trophic)) + geom_boxplot() +
  stat_summary(fun.y=mean, colour="darkred", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = .75, linetype = "dashed") + 
  geom_point(position='jitter',alpha=.3, size = 2)+
  #scale_shape_manual(values=c(1,0,2),labels=c("Bacteria", "Archaea", "Fungi") ) +
  geom_text(data=hsd1,aes(x=Trophic,y=maxdegree, label=Letter), vjust=-1)+
  labs(x="", y ="Degree") +
  #scale_fill_manual(values=c("#CC9900","#0066CC","#336633"), labels=c("CF", "NF", "NP")) + 
  theme(aspect.ratio = 1/1)+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank())+ theme(legend.position = "none")
p





df.deg.trophic.CJ1.oligo.copio <- subset(df.deg.trophic.CJ1, Trophic != "Not_predicted")
df.deg.trophic.CJ1.oligo.copio$Trophic <- factor(df.deg.trophic.CJ1.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.CJ1.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


df.deg.trophic.CJ2.oligo.copio <- subset(df.deg.trophic.CJ2, Trophic != "Not_predicted")
df.deg.trophic.CJ2.oligo.copio$Trophic <- factor(df.deg.trophic.CJ2.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.CJ2.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


df.deg.trophic.YS1.oligo.copio <- subset(df.deg.trophic.YS1, Trophic != "Not_predicted")
df.deg.trophic.YS1.oligo.copio$Trophic <- factor(df.deg.trophic.YS1.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.YS1.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


df.deg.trophic.YS2.oligo.copio <- subset(df.deg.trophic.YS2, Trophic != "Not_predicted")
df.deg.trophic.YS2.oligo.copio$Trophic <- factor(df.deg.trophic.YS2.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.YS2.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


df.deg.trophic.MY1.oligo.copio <- subset(df.deg.trophic.MY1, Trophic != "Not_predicted")
df.deg.trophic.MY1.oligo.copio$Trophic <- factor(df.deg.trophic.MY1.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.MY1.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


df.deg.trophic.MY2.oligo.copio <- subset(df.deg.trophic.MY2, Trophic != "Not_predicted")
df.deg.trophic.MY2.oligo.copio$Trophic <- factor(df.deg.trophic.MY2.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.MY2.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


df.deg.trophic.NJ1.oligo.copio <- subset(df.deg.trophic.NJ1, Trophic != "Not_predicted")
df.deg.trophic.NJ1.oligo.copio$Trophic <- factor(df.deg.trophic.NJ1.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.NJ1.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


df.deg.trophic.NJ2.oligo.copio <- subset(df.deg.trophic.NJ2, Trophic != "Not_predicted")
df.deg.trophic.NJ2.oligo.copio$Trophic <- factor(df.deg.trophic.NJ2.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.NJ2.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))



df.deg.trophic.CJ1_18.oligo.copio <- subset(df.deg.trophic.CJ1.18, Trophic != "Not_predicted")
df.deg.trophic.CJ1_18.oligo.copio$Trophic <- factor(df.deg.trophic.CJ1_18.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.CJ1_18.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


df.deg.trophic.CJ2_18.oligo.copio <- subset(df.deg.trophic.CJ2.18, Trophic != "Not_predicted")
df.deg.trophic.CJ2_18.oligo.copio$Trophic <- factor(df.deg.trophic.CJ2_18.oligo.copio$Trophic, levels = c("Copio", "Oligo"))
ggplot(data = df.deg.trophic.CJ2_18.oligo.copio, aes(x = Degree, fill = Trophic)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))


head(df.CJ_all_deg)

ggplot(data =df.CJ_all_deg, aes(x = Degree, fill = field)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))
ggplot(data =df.CJ_18_all_deg, aes(x = Degree, fill = field)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))

ggplot(data =df.MY_all_deg, aes(x = Degree, fill = field)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))
ggplot(data =df.NJ_all_deg, aes(x = Degree, fill = field)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))
ggplot(data =df.YS_all_deg, aes(x = Degree, fill = field)) + geom_area(stat = "bin")+ scale_fill_manual(values = c("#669933", "#CCCC99"))
