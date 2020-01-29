### Soil physicochemical properties
#Kruskal-Wallis
b.meta

##Kruskal-Wallis
#pH
kw<-kruskal.test(pH ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(pH ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn

#SOM
kw<-kruskal.test(SOM ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(SOM ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn


#TN
kw<-kruskal.test(TN ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(TN ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn



#CEC
kw<-kruskal.test(CEC ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(CEC ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn

#Na
kw<-kruskal.test(Na ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Na ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn


#K
kw<-kruskal.test(K ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(K ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn


#Mg
kw<-kruskal.test(Mg ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Mg ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn

#Ca
kw<-kruskal.test(Ca ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Ca ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn

#P2O5
kw<-kruskal.test(P2O5 ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(P2O5 ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn


#Sand
kw<-kruskal.test(Sand ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Sand ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn

#Silt
kw<-kruskal.test(Silt ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Silt ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn

#Clay
kw<-kruskal.test(Clay ~ Cultural_practice, data = b.meta)
kw$p.value

#library(FSA)
DT = dunnTest(Clay ~ Cultural_practice,
              data=b.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn