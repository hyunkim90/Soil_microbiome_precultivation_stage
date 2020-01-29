### Normalization

## Let's do normalization with CSS
## phyloseq to metagenomeSeq

#phy.clean? or phy.clean --> let's start from phy.clean

bac.clean.ss

arch.clean.ss

fun.clean.ss


## Remove residual taxa that do not have any sequences
#Bacteria
sum(taxa_sums(bac.clean.ss) == 0)
taxa_sums(bac.clean.ss)

bac.clean.filt <- filter_taxa(bac.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(bac.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(bac.clean.filt))

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
bac.clean.filt    ## use all samples
met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization
#https://github.com/joey711/phyloseq/issues/814
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

p <- cumNormStatFast(met.bac.clean)
p
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
#2017+2018
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
#2017+2018
bac.clean.nolog <- bac.clean.filt
otu_table(bac.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log <- bac.clean.filt
otu_table(bac.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


#Archaea
sum(taxa_sums(arch.clean.ss) == 0)
taxa_sums(arch.clean.ss)

arch.clean.filt <- filter_taxa(arch.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(arch.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(arch.clean.filt))

arch.clean.filt    ## use all samples
met.arch.clean <- phyloseq_to_metagenomeSeq(arch.clean.filt)

# normalization
p <- cumNormStatFast(met.arch.clean)
p
met.arch.norm <- cumNorm(met.arch.clean, p =p)

# returns normalized factors for each sample
normFactors(met.arch.norm)
sort(normFactors(met.arch.norm))

# To export normalized count matrix
#2017+2018
met.CSS.log <- MRcounts(met.arch.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.arch.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
#2017+2018
arch.clean.nolog <- arch.clean.filt
otu_table(arch.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

arch.clean.log <- arch.clean.filt
otu_table(arch.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

#Fungi
sum(taxa_sums(fun.clean.ss) == 0)
taxa_sums(fun.clean.ss)

fun.clean.filt <- filter_taxa(fun.clean.ss, function(x) sum(x) != 0, TRUE)
sum(taxa_sums(fun.clean.filt) == 0)

## CODE for CSS normalization using preloaded data
sort(sample_sums(fun.clean.filt))

fun.clean.filt    ## use all samples
met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization
p <- cumNormStatFast(met.fun.clean)
p
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
#2017+2018
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
#2017+2018
fun.clean.nolog <- fun.clean.filt
otu_table(fun.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log <- fun.clean.filt
otu_table(fun.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


##rarefy
bac.rarefy <- rarefy_even_depth(bac.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
arch.rarefy <- rarefy_even_depth(arch.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
fun.rarefy <- rarefy_even_depth(fun.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)