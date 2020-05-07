### YS1 and YS2
bac.clean.ss.YS12 <- subset_samples(bac.clean.ss, Field %in% c("YS1", "YS2"))
bac.clean.ss.YS12 <- phyloseq::filter_taxa(bac.clean.ss.YS12, function(x) sum(x) != 0, TRUE)

arch.clean.ss.YS12 <- subset_samples(arch.clean.ss, Field %in% c("YS1", "YS2"))
arch.clean.ss.YS12 <- phyloseq::filter_taxa(arch.clean.ss.YS12, function(x) sum(x) != 0, TRUE)

fun.clean.ss.YS12 <- subset_samples(fun.clean.ss, Field %in% c("YS1", "YS2"))
fun.clean.ss.YS12 <- phyloseq::filter_taxa(fun.clean.ss.YS12, function(x) sum(x) != 0, TRUE)

get_resSig_YS <- function(phy.clean.ss.5){
  phy.clean.filt <- phyloseq::filter_taxa(phy.clean.ss.5, function(x) sum(x) != 0, TRUE)
  sum(taxa_sums(phy.clean.filt) == 0)
  
  obj <- phyloseq_to_metagenomeSeq(phy.clean.filt)
  obj <-  cumNorm(obj, p = cumNormStatFast(obj))
  normFactor <-  normFactors(obj)
  normFactor <-  log2(normFactor/median(normFactor) + 1)
  settings <-  zigControl(maxit = 10, verbose = TRUE)
  Type  <-  pData(obj)$Field
  mod  <-  model.matrix(~Type)
  colnames(mod)  <-  levels(Type)
  colnames(mod)
  
  res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)
  res
  res$fit
  zigFit = res$fit
  finalMod= res$fit$design
  finalMod
  contrast.matrix =makeContrasts(YS1 - YS2, levels = finalMod)
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit2
  fit3 = eBayes(fit2)
  fit3
  topTable(fit3, coef="YS1 - YS2")
  
  
  res <- topTable(fit3,coef=1,adjust="fdr",n=Inf)
  
  log2AverageAbundance <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
  log2AverageAbundance
  Ta <- psmelt(phy.clean.ss.5) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
  Ta <- unique(Ta)
  Ta
  resSig = res[!is.na(res$adj.P.Val), ]
  resSig = data.frame(resSig)
  head(resSig)
  resSig <- tibble::rownames_to_column(resSig, 'OTU')
  resSig <- left_join(resSig, Ab,by= c('OTU','OTU'))
  resSig <- left_join(resSig,Ta,by=c('OTU','OTU'))
  return(resSig)
}

resSig_YS.bac<-get_resSig_YS(bac.clean.ss.YS12)
resSig_YS.arch<-get_resSig_YS(arch.clean.ss.YS12)
resSig_YS.fun<-get_resSig_YS(fun.clean.ss.YS12)


resSig_YS.arch$Enriched <- ifelse(resSig_YS.arch$adj.P.Val < .05 & resSig_YS.arch$logFC >0, "Non-fer", ifelse(resSig_YS.arch$adj.P.Val < .05 & resSig_YS.arch$logFC <0, "Fer", "ns"))
resSig_YS.bac$Enriched <- ifelse(resSig_YS.bac$adj.P.Val < .05 & resSig_YS.bac$logFC >0, "Non-fer", ifelse(resSig_YS.bac$adj.P.Val < .05 & resSig_YS.bac$logFC <0, "Fer", "ns"))
resSig_YS.fun$Enriched <- ifelse(resSig_YS.fun$adj.P.Val < .05 & resSig_YS.fun$logFC >0, "Non-fer", ifelse(resSig_YS.fun$adj.P.Val < .05 & resSig_YS.fun$logFC <0, "Fer", "ns"))

resSig_YS.arch$Enriched = as.factor(resSig_YS.arch$Enriched)
resSig_YS.bac$Enriched = as.factor(resSig_YS.bac$Enriched)
resSig_YS.fun$Enriched = as.factor(resSig_YS.fun$Enriched)


resSig_YS.arch$Enriched <- factor(resSig_YS.arch$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_YS.bac$Enriched <- factor(resSig_YS.bac$Enriched, levels = c("Non-fer", "Fer", "ns"))
resSig_YS.fun$Enriched <- factor(resSig_YS.fun$Enriched, levels = c("Non-fer", "Fer", "ns"))


write.xlsx(resSig_YS.arch, "Significant daOTU_YS_Arch.xlsx")
write.xlsx(resSig_YS.bac, "Significant daOTU_YS_Bac.xlsx")
write.xlsx(resSig_YS.fun, "Significant daOTU_YS_Fun.xlsx")

##Oligotrophic environment
otu_norm_soil_16s <- otu_table(bac.clean.nolog.YS)
design_16s_soil <- meta(bac.clean.nolog.YS)

indic_soil_16s <- as.data.frame(t(otu_norm_soil_16s))
indic_soil_groups_16s <- design_16s_soil$Field
length(unique(indic_soil_groups_16s))

## Define indicator species for soil bacteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_16s <- multipatt(indic_soil_16s,indic_soil_groups_16s,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_16s,alpha=1,indvalcomp=T)
indic_soil_df_16s <- indicatorsp_soil_16s$sign
write.table(indic_soil_df_16s,"indicsp_bacteria_YS.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_soil_df_16s <- read.table("indicsp_soil_16s_2018_nolog.txt", header=T, sep="\t")

Non_fer_soil_16s <- as.matrix(indic_soil_df_16s[which(indic_soil_df_16s$s.YS1 == 1 & indic_soil_df_16s$p.value < 0.05),])
fer_soil_16s <- as.matrix(indic_soil_df_16s[which(indic_soil_df_16s$s.YS2 == 1 & indic_soil_df_16s$p.value < 0.05),])

soil_r_values_16s_YS <- rbind(Non_fer_soil_16s,fer_soil_16s)
colnames(soil_r_values_16s_YS)[1:2] <-c("Non-fer","Fer")



##Fungi
otu_norm_soil_its <- otu_table(fun.clean.nolog.YS)
design_its_soil <- meta(fun.clean.nolog.YS)

indic_soil_its <- as.data.frame(t(otu_norm_soil_its))
indic_soil_groups_its <- design_its_soil$Field
length(unique(indic_soil_groups_its))

## Define indicator species for soil funteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_its <- multipatt(indic_soil_its,indic_soil_groups_its,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_its,alpha=1,indvalcomp=T)
indic_soil_df_its <- indicatorsp_soil_its$sign
write.table(indic_soil_df_its,"indicsp_funteria_YS.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time
Non_fer_soil_its <- as.matrix(indic_soil_df_its[which(indic_soil_df_its$s.YS1 == 1 & indic_soil_df_its$p.value < 0.05),])
fer_soil_its <- as.matrix(indic_soil_df_its[which(indic_soil_df_its$s.YS2 == 1 & indic_soil_df_its$p.value < 0.05),])

soil_r_values_its_YS <- rbind(Non_fer_soil_its,fer_soil_its)
colnames(soil_r_values_its_YS)[1:2] <-c("Non-fer","Fer")


##Archaea
otu_norm_soil_arch <- otu_table(arch.clean.nolog.YS)
design_arch_soil <- meta(arch.clean.nolog.YS)

indic_soil_arch <- as.data.frame(t(otu_norm_soil_arch))
indic_soil_groups_arch <- design_arch_soil$Field
length(unique(indic_soil_groups_arch))

## Define indicator species for soil funteria community.
## Note: These calculations can be time and processor intensive
# set.seed(8046)
indicatorsp_soil_arch <- multipatt(indic_soil_arch,indic_soil_groups_arch,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_soil_arch,alpha=1,indvalcomp=T)
indic_soil_df_arch <- indicatorsp_soil_arch$sign
write.table(indic_soil_df_arch,"indicsp_archaea_YS.txt", sep="\t",quote=F)

## Import data frame of indicator species to save time
Non_fer_soil_arch <- as.matrix(indic_soil_df_arch[which(indic_soil_df_arch$s.YS1 == 1 & indic_soil_df_arch$p.value < 0.05),])
fer_soil_arch <- as.matrix(indic_soil_df_arch[which(indic_soil_df_arch$s.YS2 == 1 & indic_soil_df_arch$p.value < 0.05),])

soil_r_values_arch_YS <- rbind(Non_fer_soil_arch,fer_soil_arch)
colnames(soil_r_values_arch_YS)[1:2] <-c("Non-fer","Fer")



resSig_YS.bac.sub <- subset(resSig_YS.bac, adj.P.Val < 0.05)

resSig_YS.fun.sub <- subset(resSig_YS.fun, adj.P.Val < 0.05)


resSig_YS.arch.sub <- subset(resSig_YS.arch, adj.P.Val < 0.05)




##Oligo and copio in daOTUs
resSig_YS.bac.sub.nonfer <- subset(resSig_YS.bac.sub, Enriched == "Non-fer")
resSig_YS.bac.sub.fer <- subset(resSig_YS.bac.sub, Enriched == "Fer")

daotu.bac_YS.nonfer<- subset(OTU_id.list, OTU%in%resSig_YS.bac.sub.nonfer$OTU)
daotu.bac_YS.fer<- subset(OTU_id.list, OTU%in%resSig_YS.bac.sub.fer$OTU)

intersect(daotu.bac_YS.nonfer$OTU_id, rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Oligo"]) #15
intersect(daotu.bac_YS.nonfer$OTU_id, rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Copio"]) #1


intersect(daotu.bac_YS.fer$OTU_id, rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Oligo"]) #1
intersect(daotu.bac_YS.fer$OTU_id, rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Copio"]) #0


resSig_YS.fun.sub.nonfer <- subset(resSig_YS.fun.sub, Enriched == "Non-fer")
resSig_YS.fun.sub.fer <- subset(resSig_YS.fun.sub, Enriched == "Fer")

daotu.fun_YS.nonfer<- subset(OTU_id.list, OTU%in%resSig_YS.fun.sub.nonfer$OTU)
daotu.fun_YS.fer<- subset(OTU_id.list, OTU%in%resSig_YS.fun.sub.fer$OTU)

intersect(daotu.fun_YS.nonfer$OTU_id, rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Oligo"]) #50
intersect(daotu.fun_YS.nonfer$OTU_id, rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Copio"]) #3


intersect(daotu.fun_YS.fer$OTU_id, rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Oligo"]) #2
intersect(daotu.fun_YS.fer$OTU_id, rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Copio"]) #8


resSig_YS.arch.sub.nonfer <- subset(resSig_YS.arch.sub, Enriched == "Non-fer")
resSig_YS.arch.sub.fer <- subset(resSig_YS.arch.sub, Enriched == "Fer")

daotu.arch_YS.nonfer<- subset(OTU_id.list, OTU%in%resSig_YS.arch.sub.nonfer$OTU)
daotu.arch_YS.fer<- subset(OTU_id.list, OTU%in%resSig_YS.arch.sub.fer$OTU)

intersect(daotu.arch_YS.nonfer$OTU_id, rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Oligo"]) #4
intersect(daotu.arch_YS.nonfer$OTU_id, rownames(df.deg.trophic.YS1)[df.deg.trophic.YS1$Trophic=="Copio"]) #1


intersect(daotu.arch_YS.fer$OTU_id, rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Oligo"]) #0
intersect(daotu.arch_YS.fer$OTU_id, rownames(df.deg.trophic.YS2)[df.deg.trophic.YS2$Trophic=="Copio"]) #0




### CJ and YS

intersect(daotu.bac_YS.nonfer$OTU_id, daotu.bac_CJ.nonfer$OTU_id) #71

intersect(daotu.bac_YS.fer$OTU_id, daotu.bac_CJ.fer$OTU_id) #66


intersect(daotu.fun_YS.nonfer$OTU_id, daotu.fun_CJ.nonfer$OTU_id) #44

intersect(daotu.fun_YS.fer$OTU_id, daotu.fun_CJ.fer$OTU_id) #6


intersect(daotu.arch_YS.nonfer$OTU_id, daotu.arch_CJ.nonfer$OTU_id) #31

intersect(daotu.arch_YS.fer$OTU_id, daotu.arch_CJ.fer$OTU_id) #0



## Random forest (?)


bac.clean.nolog.CJ.YS <- subset_samples(bac.clean.nolog, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
bac.clean.nolog.CJ.YS <- phyloseq::filter_taxa(bac.clean.nolog.CJ.YS, function(x) sum(x) != 0, TRUE)

#### All machine learning
# (1) Random forest
library(randomForest)




#All cultural practices and all years
dim(table)
table.all <- otu_table(bac.clean.nolog.CJ.YS)
table.all <- t(table.all)
rownames(table.all)

# 0, nonfer 
# 1, fer 


result.all <- c(rep(0,9),rep(1, 9), rep(0,9), rep(1,9))

#result1 <- c(rep(0,90),rep(1, 9), rep(0, 9), rep(1,9), rep(0,36), rep(1,9)) #conven vs. no fertil 2018

table.all <- cbind(table.all, result.all)
#head(table.all,20)
dim(table.all)
table.all[,c(2)]
head(table.all)
dataset_all <- as.data.frame(table.all)
head(dataset_all$result)
sort(colnames(dataset_all))

## get rid of X in front
#colnames(dataset) <- gsub("^X(.{32})", "\\1",colnames(dataset))
sort(colnames(dataset))

#write.xlsx(dataset,'rf_dataset.xlsx')

# Encoding the target feature as factor
dataset_all$result = factor(dataset_all$result, levels = c(0, 1))
dataset_all$result

colnames(dataset)
# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
library(caTools)
set.seed(123)
split = sample.split(dataset_all$result, SplitRatio = 0.66)
training_set_all = subset(dataset_all, split == TRUE)
test_set_all = subset(dataset_all, split == FALSE)
dim(training_set_all)   # when split ratio 0.75 -> 96
dim(test_set_all)  # when split ratio 0.25 -> 33


set.seed(123)
classifier_all = randomForest(x = training_set_all[-ncol(dataset_all)],
                              y = training_set_all$result, ntree=2000)
rf_pred_all = predict(classifier_all, newdata = test_set_all[-ncol(dataset_all)])
cm = table(test_set_all[, ncol(dataset_all)], rf_pred_all)
print(cm)

##Accuracy calculation based on AUC
predictions.rf=as.vector(as.numeric(rf_pred_all))
pred.rf=prediction(predictions.rf, test_set_all[,ncol(dataset_all)])
AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
(AUC.rf=AUC.rf@y.values[[1]])


##K-fold cross validation
## (1) random forest cross validation
library(caret)
folds = createFolds(training_set_all$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set_all[-x, ]
  test_fold = training_set_all[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset_all)],
                            y = training_fold$result, ntree=2000)
  
  # classifier = svm(formula = result ~ .,
  #                  data = training_fold,
  #                  type = 'C-classification',
  #                  kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset_all)])
  cm = table(test_fold[, ncol(dataset_all)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   ##1
(sd <- sd(as.numeric(cv))) ## 0


library(randomForest)


#Fungi
#All cultural practices and all years

fun.clean.nolog.CJ.YS <- subset_samples(fun.clean.nolog, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
fun.clean.nolog.CJ.YS <- phyloseq::filter_taxa(fun.clean.nolog.CJ.YS, function(x) sum(x) != 0, TRUE)
dim(table)
table.all.f <- otu_table(fun.clean.nolog.CJ.YS)
table.all.f <- t(table.all.f)
rownames(table.all.f)

# 0, CF 
# 1, NF 
# 2, NP 

result.all <- c(rep(0,9),rep(1, 9), rep(0,9), rep(1,9)) #conven vs. no fertil 2017

#result1 <- c(rep(0,90),rep(1, 9), rep(0, 9), rep(1,9), rep(0,36), rep(1,9)) #conven vs. no fertil 2018

table.all.f <- cbind(table.all.f, result.all)
#head(table.all.f,20)
dim(table.all.f)
#table.all.f[,c(365)]
head(table.all.f)
dataset_all.f <- as.data.frame(table.all.f)
dataset_all.f$result
sort(colnames(dataset_all.f))

## get rid of X in front
#colnames(dataset) <- gsub("^X(.{32})", "\\1",colnames(dataset))
#sort(colnames(dataset))

#write.xlsx(dataset,'rf_dataset.xlsx')

# Encoding the target feature as factor
dataset_all.f$result = factor(dataset_all.f$result, levels = c(0, 1))
dataset_all.f$result

colnames(dataset)
# Splitting the dataset into the Training set and Test set
# install.f.packages('caTools')
library(caTools)
set.seed(123)
split = sample.split(dataset_all.f$result, SplitRatio = 0.66)
training_set_all.f = subset(dataset_all.f, split == TRUE)
test_set_all.f = subset(dataset_all.f, split == FALSE)
dim(training_set_all.f)   # when split ratio 0.75 -> 96
dim(test_set_all.f)  # when split ratio 0.25 -> 33


set.seed(123)
classifier_all.f = randomForest(x = training_set_all.f[-ncol(dataset_all.f)],
                                y = training_set_all.f$result, ntree=2000)
rf_pred_all.f = predict(classifier_all.f, newdata = test_set_all.f[-ncol(dataset_all.f)])
cm = table(test_set_all.f[, ncol(dataset_all.f)], rf_pred_all.f)
print(cm)

##Accuracy calculation based on AUC
predictions.rf=as.vector(as.numeric(rf_pred_all.f))
pred.rf=prediction(predictions.rf, test_set_all.f[,ncol(dataset_all.f)])
AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
(AUC.rf=AUC.rf@y.values[[1]])


##K-fold cross validation
## (1) random forest cross validation
library(caret)
folds = createFolds(training_set_all.f$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set_all.f[-x, ]
  test_fold = training_set_all.f[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset_all.f)],
                            y = training_fold$result, ntree=2000)
  
  # classifier = svm(formula = result ~ .,
  #                  data = training_fold,
  #                  type = 'C-classification',
  #                  kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset_all.f)])
  cm = table(test_fold[, ncol(dataset_all.f)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   ## 0.9445455
(sd <- sd(as.numeric(cv))) ## 0.07719247


#Archaea
#All cultural practices and all years
arch.clean.nolog.CJ.YS <- subset_samples(arch.clean.nolog, Field %in% c("CJ1", "CJ2", "YS1", "YS2"))
arch.clean.nolog.CJ.YS <- phyloseq::filter_taxa(arch.clean.nolog.CJ.YS, function(x) sum(x) != 0, TRUE)

dim(table)
table.all.a <- otu_table(arch.clean.nolog.CJ.YS)
table.all.a <- t(table.all.a)
rownames(table.all.a)

# 0, CF 
# 1, NF 
# 2, NP 

result.all <-c(rep(0,9),rep(1, 9), rep(0,9), rep(1,9)) #conven vs. no fertil 2017

#result1 <- c(rep(0,90),rep(1, 9), rep(0, 9), rep(1,9), rep(0,36), rep(1,9)) #conven vs. no fertil 2018

table.all.a <- cbind(table.all.a, result.all)

#head(table.all.a,20)
dim(table.all.a)
#table.all.a[,c(365)]
head(table.all.a)
dataset_all.a <- as.data.frame(table.all.a)
dataset_all.a$result
sort(colnames(dataset_all.a))

## get rid of X in front
#colnames(dataset) <- gsub("^X(.{32})", "\\1",colnames(dataset))
#sort(colnames(dataset))

#write.xlsx(dataset,'rf_dataset.xlsx')

# Encoding the target feature as factor
dataset_all.a$result = factor(dataset_all.a$result, levels = c(0, 1))
dataset_all.a$result

colnames(dataset)
# Splitting the dataset into the Training set and Test set
# install.a.packages('caTools')
library(caTools)
set.seed(123)
split = sample.split(dataset_all.a$result, SplitRatio = 0.66)
training_set_all.a = subset(dataset_all.a, split == TRUE)
test_set_all.a = subset(dataset_all.a, split == FALSE)
dim(training_set_all.a)   # when split ratio 0.75 -> 96
dim(test_set_all.a)  # when split ratio 0.25 -> 33


set.seed(123)
classifier_all.a = randomForest(x = training_set_all.a[-ncol(dataset_all.a)],
                                y = training_set_all.a$result, ntree=2000)
rf_pred_all.a = predict(classifier_all.a, newdata = test_set_all.a[-ncol(dataset_all.a)])
cm = table(test_set_all.a[, ncol(dataset_all.a)], rf_pred_all.a)
print(cm)

##Accuracy calculation based on AUC
predictions.rf=as.vector(as.numeric(rf_pred_all.a))
pred.rf=prediction(predictions.rf, test_set_all.a[,ncol(dataset_all.a)])
AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
(AUC.rf=AUC.rf@y.values[[1]]) #0.8333333


##K-fold cross validation
## (1) random forest cross validation
library(caret)
folds = createFolds(training_set_all.a$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set_all.a[-x, ]
  test_fold = training_set_all.a[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset_all.a)],
                            y = training_fold$result, ntree=2000)
  
  # classifier = svm(formula = result ~ .,
  #                  data = training_fold,
  #                  type = 'C-classification',
  #                  kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset_all.a)])
  cm = table(test_fold[, ncol(dataset_all.a)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   ##  0.95
(sd <- sd(as.numeric(cv))) ## 0.1581139



##Error rate
colnames(classifier_all$err.rate) <- c('OOB','Non-fer','fer')
plot(classifier_all, main = 'Out-of-bag Error estimate')
legend("topright", colnames(classifier_all$err.rate),col=1:4,cex=0.8,fill=1:4)
oob <- classifier_all$err.rate[,1]
oob.oob <- oob[length(oob)]
legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob), digits = 4)*100,'%'))
dev.off()


colnames(classifier_all.f$err.rate) <- c('OOB','CF','NF','NP')
plot(classifier_all.f, main = 'Out-of-bag Error estimate')
legend("topright", colnames(classifier_all.f$err.rate),col=1:4,cex=0.8,fill=1:4)
oob.f <- classifier_all.f$err.rate[,1]
oob.oob.f <- oob.f[length(oob.f)]
legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob.f), digits = 4)*100,'%'))


colnames(classifier_all.a$err.rate) <- c('OOB','CF','NF','NP')
plot(classifier_all.a, main = 'Out-of-bag Error estimate')
legend("topright", colnames(classifier_all.a$err.rate),col=1:4,cex=0.8,fill=1:4)
oob.a <- classifier_all.a$err.rate[,1]
oob.oob.a <- oob.a[length(oob.a)]
legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob.a), digits = 4)*100,'%'))


### Mininum number of OTUs for explaining the differences between conventional farming and no fertilization
### bacteria
cv.1718.bac <- rfcv(training_set_all[-ncol(dataset_all)],training_set_all$result, cv.fold=10, step=0.9)
cv.1718.bac$error.cv 
mean(cv.1718.bac$error.cv) #0.006849315 #35

### fungi
cv2.1718.fun <- rfcv(training_set_all.f[-ncol(dataset_all.f)],training_set_all.f$result, cv.fold=10, step=0.9)
cv2.1718.fun$error.cv
mean(cv2.1718.fun$error.cv) #0.01142473 

##Archaea
cv3.1718.arch <- rfcv(training_set_all.a[-ncol(dataset_all.a)],training_set_all.a$result, cv.fold=10, step=0.9)
cv3.1718.arch$error.cv
mean(cv3.1718.arch$error.cv) # 0.03314394

### plot

par(mfrow=c(1,1), mar=c(5,5,5,2)+0.1)

plot(cv.1718.bac$n.var,cv.1718.bac$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,8100,1000))
lines(cv2.1718.fun$n.var,cv2.1718.fun$error.cv,type='o',lty=2, col=4)
lines(cv3.1718.arch$n.var,cv3.1718.arch$error.cv,type='o',lty=2, col=3)
legend('topright',legend=c("Bacteria","Fungi","Archaea"),lwd=c(1,1),col=c("red","blue", "green"),pch=1,lty=2)
abline(v=300, col="limegreen")

abline(h=mean(cv.1718.bac$error.cv), col="indianred2")
abline(h=mean(cv2.1718.fun$error.cv), col="firebrick2")
abline(h=mean(cv3.1718.arch$error.cv), col="dodgerblue4")

##Split plot
#Bacteria
plot(cv.1718.bac$n.var,cv.1718.bac$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,24000,2000))
abline(h=mean(cv.1718.bac$error.cv), col="firebrick2")
abline(v=35, col="limegreen")

#Fungi
plot(cv2.1718.fun$n.var,cv2.1718.fun$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,9500,500))
abline(h=mean(cv2.1718.fun$error.cv), col="firebrick2")
abline(v=35, col="limegreen")

#Archaea
plot(cv3.1718.arch$n.var,cv3.1718.arch$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,1200,200))
abline(h=mean(cv3.1718.arch$error.cv), col="firebrick2")
abline(v=25, col="limegreen")


df.cv1 <- data.frame(cv.1718.bac$n.var,cv.1718.bac$error.cv)
df.cv2 <- data.frame(cv2.1718.fun$n.var,cv2.1718.fun$error.cv)
df.cv3 <- data.frame(cv3.1718.arch$n.var,cv3.1718.arch$error.cv)

write.xlsx(df.cv1, 'bacteria-otu-cverror_oligo.xlsx')
write.xlsx(df.cv2, 'fungi-otu-cverror_oligo.xlsx')
write.xlsx(df.cv3,'archaea-otu_cv_error_oligo.xlsx')




## Copio
## Random forest (?)

#### All machine learning
# (1) Random forest
library(randomForest)
bac.clean.nolog.MY.NJ <- subset_samples(bac.clean.nolog, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
bac.clean.nolog.MY.NJ <- phyloseq::filter_taxa(bac.clean.nolog.MY.NJ, function(x) sum(x) != 0, TRUE)

arch.clean.nolog.MY.NJ <- subset_samples(arch.clean.nolog, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
arch.clean.nolog.MY.NJ <- phyloseq::filter_taxa(arch.clean.nolog.MY.NJ, function(x) sum(x) != 0, TRUE)

fun.clean.nolog.MY.NJ <- subset_samples(fun.clean.nolog, Field %in% c("MY1", "MY2", "NJ1", "NJ2"))
fun.clean.nolog.MY.NJ <- phyloseq::filter_taxa(fun.clean.nolog.MY.NJ, function(x) sum(x) != 0, TRUE)

#All cultural practices and all years
dim(table)
table.all <- otu_table(bac.clean.nolog.MY.NJ)
table.all <- t(table.all)
rownames(table.all)

# 0, nonfer 
# 1, fer 


result.all <- c(rep(1,9),rep(0, 9), rep(1,9), rep(0,9))

#result1 <- c(rep(0,90),rep(1, 9), rep(0, 9), rep(1,9), rep(0,36), rep(1,9)) #conven vs. no fertil 2018

table.all <- cbind(table.all, result.all)
#head(table.all,20)
dim(table.all)
#table.all[,c(365)]
head(table.all)
dataset_all <- as.data.frame(table.all)
dataset_all$result
sort(colnames(dataset_all))

## get rid of X in front
#colnames(dataset) <- gsub("^X(.{32})", "\\1",colnames(dataset))
sort(colnames(dataset))

#write.xlsx(dataset,'rf_dataset.xlsx')

# Encoding the target feature as factor
dataset_all$result = factor(dataset_all$result, levels = c(1, 0))
dataset_all$result

colnames(dataset)
# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
library(caTools)
set.seed(123)
split = sample.split(dataset_all$result, SplitRatio = 0.66)
training_set_all = subset(dataset_all, split == TRUE)
test_set_all = subset(dataset_all, split == FALSE)
dim(training_set_all)   # when split ratio 0.75 -> 96
dim(test_set_all)  # when split ratio 0.25 -> 33


set.seed(123)
classifier_all = randomForest(x = training_set_all[-ncol(dataset_all)],
                              y = training_set_all$result, ntree=2000)
rf_pred_all = predict(classifier_all, newdata = test_set_all[-ncol(dataset_all)])
cm = table(test_set_all[, ncol(dataset_all)], rf_pred_all)
print(cm)

##Accuracy calculation based on AUC
predictions.rf=as.vector(as.numeric(rf_pred_all))
pred.rf=prediction(predictions.rf, test_set_all[,ncol(dataset_all)])
AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
(AUC.rf=AUC.rf@y.values[[1]]) #0.1666667


##K-fold cross validation
## (1) random forest cross validation
library(caret)
folds = createFolds(training_set_all$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set_all[-x, ]
  test_fold = training_set_all[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset_all)],
                            y = training_fold$result, ntree=2000)
  
  # classifier = svm(formula = result ~ .,
  #                  data = training_fold,
  #                  type = 'C-classification',
  #                  kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset_all)])
  cm = table(test_fold[, ncol(dataset_all)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   ##0.85
(sd <- sd(as.numeric(cv))) ##  0.1995365


#Fungi
#All cultural practices and all years

dim(table)
table.all.f <- otu_table(fun.clean.nolog.MY.NJ)
table.all.f <- t(table.all.f)
rownames(table.all.f)

# 0, CF 
# 1, NF 
# 2, NP 

result.all <- c(rep(1,9),rep(0, 9), rep(1,9), rep(0,9)) #conven vs. no fertil 2017

#result1 <- c(rep(0,90),rep(1, 9), rep(0, 9), rep(1,9), rep(0,36), rep(1,9)) #conven vs. no fertil 2018

table.all.f <- cbind(table.all.f, result.all)
#head(table.all.f,20)
dim(table.all.f)
#table.all.f[,c(365)]
head(table.all.f)
dataset_all.f <- as.data.frame(table.all.f)
dataset_all.f$result
sort(colnames(dataset_all.f))

## get rid of X in front
#colnames(dataset) <- gsub("^X(.{32})", "\\1",colnames(dataset))
#sort(colnames(dataset))

#write.xlsx(dataset,'rf_dataset.xlsx')

# Encoding the target feature as factor
dataset_all.f$result = factor(dataset_all.f$result, levels = c(1, 0))
dataset_all.f$result

colnames(dataset)
# Splitting the dataset into the Training set and Test set
# install.f.packages('caTools')
library(caTools)
set.seed(123)
split = sample.split(dataset_all.f$result, SplitRatio = 0.66)
training_set_all.f = subset(dataset_all.f, split == TRUE)
test_set_all.f = subset(dataset_all.f, split == FALSE)
dim(training_set_all.f)   # when split ratio 0.75 -> 96
dim(test_set_all.f)  # when split ratio 0.25 -> 33


set.seed(123)
classifier_all.f = randomForest(x = training_set_all.f[-ncol(dataset_all.f)],
                                y = training_set_all.f$result, ntree=2000)
rf_pred_all.f = predict(classifier_all.f, newdata = test_set_all.f[-ncol(dataset_all.f)])
cm = table(test_set_all.f[, ncol(dataset_all.f)], rf_pred_all.f)
print(cm)

##Accuracy calculation based on AUC
predictions.rf=as.vector(as.numeric(rf_pred_all.f))
pred.rf=prediction(predictions.rf, test_set_all.f[,ncol(dataset_all.f)])
AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
(AUC.rf=AUC.rf@y.values[[1]])


##K-fold cross validation
## (1) random forest cross validation
library(caret)
folds = createFolds(training_set_all.f$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set_all.f[-x, ]
  test_fold = training_set_all.f[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset_all.f)],
                            y = training_fold$result, ntree=2000)
  
  # classifier = svm(formula = result ~ .,
  #                  data = training_fold,
  #                  type = 'C-classification',
  #                  kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset_all.f)])
  cm = table(test_fold[, ncol(dataset_all.f)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   ## 0.8666667
(sd <- sd(as.numeric(cv))) ## 0.2194269


#Archaea
#All cultural practices and all years
dim(table)
table.all.a <- otu_table(arch.clean.nolog.MY.NJ)
table.all.a <- t(table.all.a)
rownames(table.all.a)

# 0, CF 
# 1, NF 
# 2, NP 

result.all <-c(rep(1,9),rep(0, 9), rep(1,9), rep(0,9)) #conven vs. no fertil 2017

#result1 <- c(rep(0,90),rep(1, 9), rep(0, 9), rep(1,9), rep(0,36), rep(1,9)) #conven vs. no fertil 2018

table.all.a <- cbind(table.all.a, result.all)

#head(table.all.a,20)
dim(table.all.a)
#table.all.a[,c(365)]
head(table.all.a)
dataset_all.a <- as.data.frame(table.all.a)
dataset_all.a$result
sort(colnames(dataset_all.a))

## get rid of X in front
#colnames(dataset) <- gsub("^X(.{32})", "\\1",colnames(dataset))
#sort(colnames(dataset))

#write.xlsx(dataset,'rf_dataset.xlsx')

# Encoding the target feature as factor
dataset_all.a$result = factor(dataset_all.a$result, levels = c(1, 0))
dataset_all.a$result

colnames(dataset)
# Splitting the dataset into the Training set and Test set
# install.a.packages('caTools')
library(caTools)
set.seed(123)
split = sample.split(dataset_all.a$result, SplitRatio = 0.66)
training_set_all.a = subset(dataset_all.a, split == TRUE)
test_set_all.a = subset(dataset_all.a, split == FALSE)
dim(training_set_all.a)   # when split ratio 0.75 -> 96
dim(test_set_all.a)  # when split ratio 0.25 -> 33


set.seed(123)
classifier_all.a = randomForest(x = training_set_all.a[-ncol(dataset_all.a)],
                                y = training_set_all.a$result, ntree=2000)
rf_pred_all.a = predict(classifier_all.a, newdata = test_set_all.a[-ncol(dataset_all.a)])
cm = table(test_set_all.a[, ncol(dataset_all.a)], rf_pred_all.a)
print(cm)

##Accuracy calculation based on AUC
predictions.rf=as.vector(as.numeric(rf_pred_all.a))
pred.rf=prediction(predictions.rf, test_set_all.a[,ncol(dataset_all.a)])
AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
(AUC.rf=AUC.rf@y.values[[1]])


##K-fold cross validation
## (1) random forest cross validation
library(caret)
folds = createFolds(training_set_all.a$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set_all.a[-x, ]
  test_fold = training_set_all.a[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset_all.a)],
                            y = training_fold$result, ntree=2000)
  
  # classifier = svm(formula = result ~ .,
  #                  data = training_fold,
  #                  type = 'C-classification',
  #                  kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset_all.a)])
  cm = table(test_fold[, ncol(dataset_all.a)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   ## 0.9666667
(sd <- sd(as.numeric(cv))) ## 0.1054093



##Error rate
colnames(classifier_all$err.rate) <- c('OOB','Non-fer','fer')
plot(classifier_all, main = 'Out-of-bag Error estimate')
legend("topright", colnames(classifier_all$err.rate),col=1:4,cex=0.8,fill=1:4)
oob <- classifier_all$err.rate[,1]
oob.oob <- oob[length(oob)]
legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob), digits = 4)*100,'%'))
dev.off()


colnames(classifier_all.f$err.rate) <- c('OOB','Non-fer','fer')
plot(classifier_all.f, main = 'Out-of-bag Error estimate')
legend("topright", colnames(classifier_all.f$err.rate),col=1:4,cex=0.8,fill=1:4)
oob.f <- classifier_all.f$err.rate[,1]
oob.oob.f <- oob.f[length(oob.f)]
legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob.f), digits = 4)*100,'%'))


colnames(classifier_all.a$err.rate) <- c('OOB','Non-fer','fer')
plot(classifier_all.a, main = 'Out-of-bag Error estimate')
legend("topright", colnames(classifier_all.a$err.rate),col=1:4,cex=0.8,fill=1:4)
oob.a <- classifier_all.a$err.rate[,1]
oob.oob.a <- oob.a[length(oob.a)]
legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob.a), digits = 4)*100,'%'))


### Mininum number of OTUs for explaining the differences between conventional farming and no fertilization
### bacteria
cv.1718.bac <- rfcv(training_set_all[-ncol(dataset_all)],training_set_all$result, cv.fold=10, step=0.9)
cv.1718.bac$error.cv 
mean(cv.1718.bac$error.cv) #0.09018265

### fungi
cv2.1718.fun <- rfcv(training_set_all.f[-ncol(dataset_all.f)],training_set_all.f$result, cv.fold=10, step=0.9)
cv2.1718.fun$error.cv
mean(cv2.1718.fun$error.cv) #0.07275132

##Archaea
cv3.1718.arch <- rfcv(training_set_all.a[-ncol(dataset_all.a)],training_set_all.a$result, cv.fold=10, step=0.9)
cv3.1718.arch$error.cv
mean(cv3.1718.arch$error.cv) # 0.01018519

### plot

par(mfrow=c(1,1), mar=c(5,5,5,2)+0.1)

plot(cv.1718.bac$n.var,cv.1718.bac$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,23000,1000))
lines(cv2.1718.fun$n.var,cv2.1718.fun$error.cv,type='o',lty=2, col=4)
lines(cv3.1718.arch$n.var,cv3.1718.arch$error.cv,type='o',lty=2, col=3)
legend('topright',legend=c("Bacteria","Fungi","Archaea"),lwd=c(1,1),col=c("red","blue", "green"),pch=1,lty=2)
abline(v=300, col="limegreen")

abline(h=mean(cv.1718.bac$error.cv), col="indianred2")
abline(h=mean(cv2.1718.fun$error.cv), col="firebrick2")
abline(h=mean(cv3.1718.arch$error.cv), col="dodgerblue4")

##Split plot
#Bacteria
plot(cv.1718.bac$n.var,cv.1718.bac$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,24000,2000))
abline(h=mean(cv.1718.bac$error.cv), col="firebrick2")
abline(v=30, col="limegreen")

#Fungi
plot(cv2.1718.fun$n.var,cv2.1718.fun$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,9500,500))
abline(h=mean(cv2.1718.fun$error.cv), col="firebrick2")
abline(v=30, col="limegreen")

#Archaea
plot(cv3.1718.arch$n.var,cv3.1718.arch$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,1200,200))
abline(h=mean(cv3.1718.arch$error.cv), col="firebrick2")
abline(v=30, col="limegreen")


df.cv1 <- data.frame(cv.1718.bac$n.var,cv.1718.bac$error.cv)
df.cv2 <- data.frame(cv2.1718.fun$n.var,cv2.1718.fun$error.cv)
df.cv3 <- data.frame(cv3.1718.arch$n.var,cv3.1718.arch$error.cv)

write.xlsx(df.cv1, 'bacteria-otu-cverror_copio_re.xlsx')
write.xlsx(df.cv2, 'fungi-otu-cverror_copio_re.xlsx')
write.xlsx(df.cv3,'archaea-otu_cv_error_copio_re.xlsx')
