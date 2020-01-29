##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
order.location <- c("CC", "CJ","DG", "IS", "JJ", "MY", "NJ", "UF", "YS", "CC.18", "CJ.18", "DG.18", "IS.18", "UF.18")
order.Soil_texture <- c('Sandy loam', 'Loam', "Silt loam", 'Clay loam', 'Silt clay loam', 'Silty Clay')

b.meta$Location <- factor(b.meta$Location, levels=order.location)
f.meta$Location <- factor(f.meta$Location, levels=order.location)

b.meta$Soil_texture <- factor(b.meta$Soil_texture, levels=order.Soil_texture)
f.meta$Soil_texture <- factor(f.meta$Soil_texture, levels=order.Soil_texture)

bray1.bac <-  ordinate(bac.clean.log, 'PCoA', 'bray')
bray1.arch <-  ordinate(arch.clean.log, 'PCoA', 'bray')
bray1.fun <-  ordinate(fun.clean.log, 'PCoA', 'bray')


library(Hmisc)
corr.all<-cbind(b.meta[,-c(1:3,5:6)], bray1.bac$vector)
class(corr.all)
corr.all<-as.matrix(corr.all)

cor_5 <- rcorr(corr.all, type="spearman")
M <- cor_5$r
p_mat <- cor_5$P

write.xlsx(M, "spearman correlation_PCoA axis and soil properties_bacteria_all.xlsx")
write.xlsx(p_mat, "P value spearman correlation_PCoA axis and soil properties_bacteria_all.xlsx")

#Fungi
corr.all<-cbind(f.meta[,-c(1:3,5:6)], bray1.fun$vector)
class(corr.all)
corr.all<-as.matrix(corr.all)

cor_5 <- rcorr(corr.all, type="spearman")
M <- cor_5$r
p_mat <- cor_5$P

write.xlsx(M, "spearman correlation_PCoA axis and soil properties_fungi_all.xlsx")
write.xlsx(p_mat, "P value spearman correlation_PCoA axis and soil properties_fungi_all.xlsx")

#Archaea
corr.all<-cbind(b.meta[,-c(1:3,5:6)], bray1.arch$vector)
class(corr.all)
corr.all<-as.matrix(corr.all)

cor_5 <- rcorr(corr.all, type="spearman")
M <- cor_5$r
p_mat <- cor_5$P

write.xlsx(M, "Spearman correlation_PCoA axis and soil properties_archaea_all.xlsx")
write.xlsx(p_mat, "P value Spearman correlation_PCoA axis and soil properties_archaea_all.xlsx")



##Unconstrained PCoA
#Location
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='Location',shape = "Location", title = "PCoA (Bray-curtis)", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+scale_shape_manual(values=c(5,6,15,16,17,18,7,13,21,5,6,15,16,13))+scale_colour_manual(values=c("#658B89","#C9C063",  "#FF5534",'#104E8D','#C88843','#9AC94A','#9FB4CE','#CF96CD','#7661EE', "#AEC4C2", "#7E762C", "#FF937B", "#85BBF1", "#662F65"))

plot_ordination(arch.clean.log, bray1.arch, type = "samples", color='Location', shape = "Location",title = "PCoA (Bray-curtis)", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+scale_shape_manual(values=c(5,6,15,16,17,18,7,13,21,5,6,15,16,13))+scale_colour_manual(values=c("#658B89","#C9C063",  "#FF5534",'#104E8D','#C88843','#9AC94A','#9FB4CE','#CF96CD','#7661EE', "#AEC4C2", "#7E762C", "#FF937B", "#85BBF1", "#662F65"))

plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='Location', shape = "Location",title = "PCoA (Bray-curtis)", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+scale_shape_manual(values=c(5,6,15,16,17,18,7,13,21,5,6,15,16,13))+scale_colour_manual(values=c("#658B89","#C9C063",  "#FF5534",'#104E8D','#C88843','#9AC94A','#9FB4CE','#CF96CD','#7661EE', "#AEC4C2", "#7E762C", "#FF937B", "#85BBF1", "#662F65"))

#Soil_texture 
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='Soil_texture',shape = "Soil_texture", title = "PCoA (Bray-curtis)", axes = c(2,3))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+scale_shape_manual(values=c(3,9,10,15,16,17))+scale_colour_manual(values=c('#000000','#00AFBB','#E7B800','#FC4E08','#0072B2','#CC79A7'))

plot_ordination(arch.clean.log, bray1.arch, type = "samples", color='Soil_texture', shape = "Soil_texture",title = "PCoA (Bray-curtis)", axes = c(3,4))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+scale_shape_manual(values=c(3,9,10,15,16,17))+scale_colour_manual(values=c('#000000','#00AFBB','#E7B800','#FC4E08','#0072B2','#CC79A7'))

plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='Soil_texture', shape = "Soil_texture",title = "PCoA (Bray-curtis)", axes = c(2,3))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+scale_shape_manual(values=c(3,9,10,15,16,17))+scale_colour_manual(values=c('#000000','#00AFBB','#E7B800','#FC4E08','#0072B2','#CC79A7'))


#P2O5
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='P2O5', title = "PCoA (Bray-curtis)", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)

plot_ordination(arch.clean.log, bray1.arch, type = "samples", color='P2O5',title = "PCoA (Bray-curtis)", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)

plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='P2O5',title = "PCoA (Bray-curtis)", axes = c(1,3))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)


#Mg
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='Mg', title = "PCoA (Bray-curtis)", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)

plot_ordination(arch.clean.log, bray1.arch, type = "samples", color='Mg',title = "PCoA (Bray-curtis)", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)

plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='Mg',title = "PCoA (Bray-curtis)", axes = c(1,3))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)

dev.off()
