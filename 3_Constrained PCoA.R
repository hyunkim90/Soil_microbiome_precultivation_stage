# script to reproduce the analysis and figures from Bulgarelli et al., 2015
#
# if you use any of the following code, please cite:
#
# Davide Bulgarelli, Ruben Garrido-Oter, Philipp C. Münch, Aaron Weiman,
# Johannes Dröge, Yao Pan, Alice C. McHardy, Paul Schulze-Lefert
# Structure and Function of the Bacterial Root Microbiota in Wild
# and Domesticated Barley. Cell Host and Microbe, 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

### ananlysis of beta-diversity (CPCoA)

# load functions
source("cpcoa.plot.func.R")


## Bacteria
#OTU count
count_matrix <- otu_table(bac.clean.log)
count_matrix

##Design file
library(microbiome)
d <- meta(bac.clean.log)
d <- data.frame(Soil_texture=d$Soil_texture,
                Field=d$Field,
                Cultural_practice=d$Cultural_practice,
                Location=d$Location,
                row.names=d$SampleID)
d




## Archaea
#OTU count
count_matrix <- otu_table(arch.clean.log)
count_matrix

##Design file
library(microbiome)
d <- meta(arch.clean.log)
d <- data.frame(Soil_texture=d$Soil_texture,
                Field=d$Field,
                Cultural_practice=d$Cultural_practice,
                Location=d$Location,
                row.names=d$SampleID)
d



## Fungi
#OTU count
count_matrix <- otu_table(fun.clean.log)
count_matrix

##Design file
library(microbiome)
d <- meta(fun.clean.log)
d <- data.frame(Soil_texture=d$Soil_texture,
                Field=d$Field,
                Cultural_practice=d$Cultural_practice,
                Location=d$Location,
                row.names=d$SampleID)
d


# apply transformations to OTU count table
matrix_log_transformed <- count_matrix
count_matrix <- count_matrix[, match(colnames(matrix_log_transformed), colnames(count_matrix))]


Norm.count <- matrix_log_transformed

# CPCoA + compartment and genotype using bray-curtis distance

sqrt_transform <- T

# bray-curtis

capscale.soil <- capscale(t(Norm.count) ~ Soil_texture,
                              data=d, add=F, sqrt.dist=sqrt_transform, dist="bray")

capscale.culture <- capscale(t(Norm.count) ~ Cultural_practice,
                             data=d, add=F, sqrt.dist=sqrt_transform, dist="bray")

capscale.location <- capscale(t(Norm.count) ~ Location,
                              data=d, add=F, sqrt.dist=sqrt_transform, dist="bray")

# ANOVA-like permutation analysis
# anova.cca is a vegan wrapper for CCA, which uses the function permutest

perm_anova.soil <- anova.cca(capscale.soil)
print(perm_anova.soil)

perm_anova.culture <- anova.cca(capscale.culture)
print(perm_anova.culture)

perm_anova.location <- anova.cca(capscale.location)
print(perm_anova.location)


# generate variability tables and calculate confidence intervals for the variance

var_tbl.soil <- variability_table(capscale.soil)
cap_var.soil <- cap_var_props(capscale.soil)
ci.soil <-  cca_ci(capscale.soil) 

var_tbl.culture <- variability_table(capscale.culture)
cap_var.culture <- cap_var_props(capscale.culture)
ci.culture <-  cca_ci(capscale.culture) 


var_tbl.location <- variability_table(capscale.location)
cap_var.location <- cap_var_props(capscale.location)
ci.location <-  cca_ci(capscale.location) 

# extract the weighted average (sample) scores

wa.soil <- capscale.soil.bac$CCA$wa
wa.soil

wa.culture <- capscale.culture$CCA$wa
wa.culture

wa.location <- capscale.location$CCA$wa
wa.location
# extract OTU scores from CCA object
capscale.soil
otu.scores.soil <- capscale.soil$CCA$v[, 1:2] 

otu.scores.culture <- capscale.culture$CCA$v[, 1:2] 

otu.scores.location <- capscale.location$CCA$v[, 1:2]

# extract centroids of constrained factor

centroids.soil.bac <- capscale.soil.bac$CCA$centroids[, 1:2] 
centroids.soil.bac

centroids.culture <- capscale.culture$CCA$centroids[, 1:2]

centroids.location <- capscale.location$CCA$centroids[, 1:2]

# plot CPCoA (sample scores)
plot_cap(p=wa.soil, d=d, col_var=c('Sandy loam', 'Loam', "Silt loam", 'Clay loam', 'Silt clay loam', 'Silty Clay'),
         pch_var=c('Sandy loam', 'Loam', "Silt loam", 'Clay loam', 'Silt clay loam', 'Silty Clay'),
         col_comp=c("Soil_texture", "Soil_texture","Soil_texture","Soil_texture","Soil_texture","Soil_texture"),
         pch_comp=c("Soil_texture", "Soil_texture","Soil_texture","Soil_texture","Soil_texture","Soil_texture"), shapes=c(3,9,10,15,16,17),
         colors=c('#000000','#00AFBB','#E7B800','#FC4E08','#0072B2','#CC79A7'), file_name="Soil Texture_CPCoA_all_rooted.pdf",
         svg=F, constraint="Bray-Curtis - constrained by Soil texture",
         ci=ci.soil, var_tbl=var_tbl.soil,
         cap_var=cap_var.soil, perm_anova=perm_anova.soil,
         cex.main=0.8, cex=1.3)


plot_cap(p=wa.culture, d=d, col_var=c('Conventional', 'No_pesticide', "No_fertilizer"),
         pch_var=c('Conventional', 'No_pesticide', "No_fertilizer"),
         col_comp=c("Cultural_practice", "Cultural_practice","Cultural_practice"),
         pch_comp=c("Cultural_practice", "Cultural_practice","Cultural_practice"), shapes=c(15,16,17),
         colors=c("Red","Blue",  "darkgreen"), file_name="Cultural practice_CPCoA_all_rooted.pdf",
         svg=F, constraint="Bray-Curtis - constrained by Cultural Practice",
         ci=ci.culture, var_tbl=var_tbl.culture,
         cap_var=cap_var.culture, perm_anova=perm_anova.culture,
         cex.main=0.8, cex=1.3)

plot_cap(p=wa.location, d=d, col_var=c('CC', 'CJ', "DG", 'IS', 'JJ', 'MY', 'NJ', 'UF', 'YS','CC.18', 'CJ.18', "DG.18", 'IS.18', 'UF.18'),
         pch_var=c('CC', 'CJ', "DG", 'IS', 'JJ', 'MY', 'NJ', 'UF', 'YS'),
         col_comp=c('Location', 'Location', 'Location', 'Location', 'Location', 'Location', 'Location', 'Location', 'Location','Location', 'Location', 'Location', 'Location', 'Location'),
         pch_comp=c('Location', 'Location', 'Location', 'Location', 'Location', 'Location', 'Location', 'Location', 'Location','Location', 'Location', 'Location', 'Location', 'Location', 'Location'), shapes=c(5,6,15,16,17,18,7,13,21,5,6,15,16,7),
         colors=c("#658B89","#C9C063",  "#FF5534",'#104E8D','#C88843','#9AC94A','#9FB4CE','#CF96CD','#7661EE', '#B38184', '#A9DBA8','#E94E76','#403E4B','#C6A49A'), file_name="Location_CPCoA_all_rooted.pdf",
         svg=F, constraint="Bray-Curtis - constrained by location",
         ci=ci.location, var_tbl=var_tbl.location,
         cap_var=cap_var.location, perm_anova=perm_anova.location,
         cex.main=0.8, cex=1.3)
