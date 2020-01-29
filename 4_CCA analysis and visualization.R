#####Canonical Correspondence Analysis #######
library(microbiome)

#Bacteria
bac.log <- microbiome::transform(bac.clean.log, "identity")
b.otu <- abundances(bac.log)
b.meta <- meta(bac.log)
sapply(b.meta, class) ## soil properties are class 'factor'
indx.div <- sapply(b.meta[7:18], is.factor)
b.meta[7:18][indx.div] <- lapply(b.meta[7:18][indx.div], function(x) as.numeric(as.character(x)))
head(b.meta)
sapply(b.meta, class)

#Use adonis to find significant environmental variables PERMANOVA
abund_table.adonis <- adonis(formula = t(b.otu) ~ (pH+SOM+TN+CEC+Ca+Mg+K+Na+P2O5+Sand+Silt+Clay), data = b.meta, permutations=9999, method = "bray")
abund_table.adonis

#Extract the best variables

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=0.01]
bestEnvVariables

#Last two are NA entries, so we have to remove them
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(t(b.otu) ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=b.meta)",sep="")))

#You can use the following to use all the environmental variables
#sol<-cca(t(otu) ~ ., data=meta)
sol

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
scrs

#Check the attributes
attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,t(as.data.frame(strsplit(rownames(scrs$sites),"_"))))
df_sites

df_sites$Location <- b.meta$Location
df_sites <- df_sites[,c("CCA1", "CCA2", "Location")]

df_sites

class(df_sites)

colnames(df_sites)
colnames(df_sites)<-c("CCA1", "CCA2", "Location")


#Draw sites
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(CCA1,CCA2, colour=Location), size =4) +scale_color_manual(values = my_color_collection)
p


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
multiplier

# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines 
# which are the coordinates of the heads of unit length vectors. In plot these are 
# scaled by their correlation (square root of the column r2) so that "weak" predictors 
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths 
# using command scores. The plotted (and scaled) arrows are further adjusted to the 
# current graph using a constant multiplier: this will keep the relative r2-scaled 
# lengths of the arrows but tries to fill the current plot. You can see the multiplier 
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul. 

df_arrows<- scrs$biplot*multiplier
df_arrows

colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.3, "cm")),color="black",alpha=1)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="black",alpha=1)

p

dev.off()



## Archaea
arch.log <- microbiome::transform(arch.clean.log, "identity")
a.otu <- abundances(arch.log)
a.meta <- meta(arch.log)
sapply(a.meta, class) ## soil properties are class 'factor'
indx.div <- sapply(a.meta[7:18], is.factor)
a.meta[7:18][indx.div] <- lapply(a.meta[7:18][indx.div], function(x) as.numeric(as.character(x)))
head(a.meta)
sapply(a.meta, class)

#Use adonis to find significant environmental variables PERMANOVA
abund_table.adonis <- adonis(formula = t(a.otu) ~ (pH+SOM+TN+CEC+Ca+Mg+K+Na+P2O5+Sand+Silt+Clay), data = a.meta, permutations=9999, method = "bray")
abund_table.adonis

#Extract the best variables

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=0.01]
bestEnvVariables

#Last two are NA entries, so we have to remove them
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(t(a.otu) ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=a.meta)",sep="")))

#You can use the following to use all the environmental variables
#sol<-cca(t(otu) ~ ., data=meta)
sol

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
scrs

#Check the attributes
attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,t(as.data.frame(strsplit(rownames(scrs$sites),"_"))))
df_sites

df_sites$Location <- a.meta$Location
df_sites <- df_sites[,c("CCA1", "CCA2", "Location")]

df_sites

class(df_sites)

colnames(df_sites)
colnames(df_sites)<-c("CCA1", "CCA2", "Location")


#Draw sites
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(CCA1,CCA2, colour=Location), size =4) +scale_color_manual(values = my_color_collection)
p


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
multiplier

# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines 
# which are the coordinates of the heads of unit length vectors. In plot these are 
# scaled by their correlation (square root of the column r2) so that "weak" predictors 
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths 
# using command scores. The plotted (and scaled) arrows are further adjusted to the 
# current graph using a constant multiplier: this will keep the relative r2-scaled 
# lengths of the arrows but tries to fill the current plot. You can see the multiplier 
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul. 

df_arrows<- scrs$biplot*multiplier
df_arrows

colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.3, "cm")),color="black",alpha=1)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="black",alpha=1)

p

dev.off()


#Fungi
fun.log <- microbiome::transform(fun.clean.log, "identity")
f.otu <- abundances(fun.log)
f.meta <- meta(fun.log)
sapply(f.meta, class) ## soil properties are class 'factor'
indx.div <- sapply(f.meta[7:18], is.factor)
f.meta[7:18][indx.div] <- lapply(f.meta[7:18][indx.div], function(x) as.numeric(as.character(x)))
head(f.meta)
sapply(f.meta, class)

#Use adonis to find significant environmental variables PERMANOVA
abund_table.adonis <- adonis(formula = t(f.otu) ~ (pH+SOM+TN+CEC+Ca+Mg+K+Na+P2O5+Sand+Silt+Clay), data = f.meta, permutations=9999, method = "bray")
abund_table.adonis

#Extract the best variables

bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=0.01]
bestEnvVariables

#Last two are NA entries, so we have to remove them
bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables

#We are now going to use only those environmental variables in cca that were found significant
eval(parse(text=paste("sol <- cca(t(f.otu) ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=f.meta)",sep="")))

#You can use the following to use all the environmental variables
#sol<-cca(t(otu) ~ ., data=meta)
sol

scrs<-scores(sol,display=c("sp","wa","lc","bp","cn"))
scrs

#Check the attributes
attributes(scrs)
# $names
# [1] "species"     "sites"       "constraints" "biplot"     
# [5] "centroids"  

#Extract site data first
df_sites<-data.frame(scrs$sites,t(as.data.frame(strsplit(rownames(scrs$sites),"_"))))
df_sites

df_sites$Location <- f.meta$Location
df_sites <- df_sites[,c("CCA1", "CCA2", "Location")]

df_sites

class(df_sites)

colnames(df_sites)
colnames(df_sites)<-c("CCA1", "CCA2", "Location")


#Draw sites
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(CCA1,CCA2, colour=Location), size =4) +scale_color_manual(values = my_color_collection)
p


#Draw biplots
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
multiplier

# Reference: http://www.inside-r.org/packages/cran/vegan/docs/envfit
# The printed output of continuous variables (vectors) gives the direction cosines 
# which are the coordinates of the heads of unit length vectors. In plot these are 
# scaled by their correlation (square root of the column r2) so that "weak" predictors 
# have shorter arrows than "strong" predictors. You can see the scaled relative lengths 
# using command scores. The plotted (and scaled) arrows are further adjusted to the 
# current graph using a constant multiplier: this will keep the relative r2-scaled 
# lengths of the arrows but tries to fill the current plot. You can see the multiplier 
# using vegan:::ordiArrowMul(result_of_envfit), and set it with the argument arrow.mul. 

df_arrows<- scrs$biplot*multiplier
df_arrows

colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)

p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.3, "cm")),color="black",alpha=1)

p<-p+geom_text(data=as.data.frame(df_arrows*1.1),aes(x, y, label = rownames(df_arrows)),color="black",alpha=1)

p

dev.off()