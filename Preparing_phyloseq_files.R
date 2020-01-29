## import libraries : phyloseq and microbiome
library(dplyr)
library(forcats) 
library(metagenomeSeq)
library(vegan)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(seqtime)
library(agricolae)
library(RColorBrewer)
library(xlsx)
library(magrittr)
library(indicspecies)
library(Hmisc)
library(igraph)
library(qgraph)
library(randomForest)
library(multifunc)
library(OTUtable)
library(FSA)
library(rcompanion)


# Set plotting theme
theme_set(theme_bw())

## Setting working directory
setwd("C:/Users/hyunk/Desktop/Our paper/Soil/Raw data/Bacterial biome/The latest raw data set")

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
bac_phylo=import_biom("OTU_table_final.biom")

### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
map <- read.table(file = 'Soil_metadata_with_environment_bac.tsv', sep = '\t', header = TRUE)
map <- sample_data(map)

head(map)
dim(map)
summary(map)
str(map)

summary(map)
colnames(map)
rownames(map)
nrow(map)

# Assign rownames to be Sample ID's
map$SampleID
rownames(map) <- map$SampleID
rownames(map)
dim(map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
phy_tree = read_tree("rooted_tree.nwk")
phy <- merge_phyloseq(bac_phylo, map, phy_tree)

class(phy)
phy   ## 629 otus

## changing rank names
colnames(tax_table(phy)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy  ## 524 OTUs



## Fungal community

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
fun_phylo=import_biom("otu_table_final.biom")

### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
f.map <- read.table(file = 'Soil_metadata_with_environment_fun.tsv', sep = '\t', header = TRUE)
f.map <- sample_data(f.map)

head(f.map)
dim(f.map)
summary(f.map)
str(f.map)

summary(f.map)
colnames(f.map)
rownames(f.map)
nrow(f.map)

# Assign rownames to be Sample ID's
f.map$SampleID
rownames(f.map) <- f.map$SampleID
rownames(f.map)
dim(f.map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
fun_tree = read_tree("rooted_tree.nwk")
fun <- merge_phyloseq(fun_phylo, f.map, fun_tree)

class(fun)
fun   ## 10981 otus

## changing rank names
colnames(tax_table(fun)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fun



### 

# Set colors for plotting
phylum_colors <- c(
  "gray",'black', "#DA5724", "#5F7FC7","#508578", "#CD9BCD", "orange",
  "#5F7FC7","#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#D1A33D", "#8A7C64", "#599861","#5E738F"
)

my_color_collection <- c(
  "#CBD588", "#5F7FC7", "orange", "#AD6F3B", "#673770", 
  "#D14285", "#652926", "#C84248", "#8569D5", "#5E738F",
  "#D1A33D", "#8A7C64", "#599861","#616163", "#FFCDB2",
  "#6D9F71", "#242F40",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")

my_color_Class <- c(
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "orange", "#5F7FC7", "#CBD588", "#AD6F3B", "#673770",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")



my_color_Order <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B", 
  "#F92A82", "#ED7B84", "#5F7FC7", 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)

my_color_Family <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B", 
  "#F92A82", "#ED7B84", "#5F7FC7", 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)
my_color_Family2 <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B", 
  "#F92A82", "#ED7B84", "#5F7FC7", 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)

# up to 300
my_color_OTU <- c(
  "gray",'black', "#5F7FC7", "orange",  "#AD6F3B", 
  "#673770","#D14285", "#652926", "#C84248",  "#8569D5",
  "#5E738F","#D1A33D", "#8A7C64", "#599861","#616163", 
  "#FFCDB2", "#242F40", "#6D9F71",  "#CCA43B", "#F92A82",
  "#ED7B84", "#7EB77F", "#DEC4A1", "#E5D1D0", '#0E8482', 
  '#C9DAEA', '#337357', '#95C623', '#E55812', '#04471C',
  '#F2D7EE', '#D3BCC0', '#A5668B', '#69306D', 'navy', 
  '#1A535C', '#4ECDC4', 'orange', '#FF6B6B', "orchid1",
  'cyan2', '#FFF275', 'springgreen', '#FF3C38', '#A23E48', 
  '#000000', '#CF5C36', '#EEE5E9', '#7C7C7C', '#EFC88B',
  
  '#2E5266', '#6E8898', '#9FB1BC', '#D3D0CB', '#E2C044', 
  '#5BC0EB', '#FDE74C', '#9BC53D', '#E55934', '#FA7921',
  "#CD9BCD", "#508578", "#CBD588","#CBD588", "#5F7FC7", 
  "orange",   "#AD6F3B", "#673770","#D14285", "#652926", 
  "#C84248",  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", 
  "#599861","#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", "#DEC4A1", 
  "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', '#95C623', 
  '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', '#A5668B', 
  '#69306D', '#0E103D', '#1A535C', '#4ECDC4', '#F7FFF7', 
  '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange")

my_color_otu2 <- c(
  "gray",'black', "#AD6F3B", "#673770","#D14285", 
  "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", 
  "#8A7C64", "#599861", "#616163",  "#FFCDB2", "#242F40", 
  "#6D9F71", "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#9FB1BC', 'springgreen', '#E2C044', '#5BC0EB', 'pink', 
  "orange", "#CBD588", "#5F7FC7",  
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898', 
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C', 
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "orange")


my_color_gen <- c(
  "white",'yellow', "#AD6F3B", "#673770","#D14285", 
  "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", 
  "#8A7C64", "#599861", "#616163",  "#FFCDB2", "#242F40", 
  "#6D9F71", "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', 
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', 
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4', 
  '#F7FFF7', '#FF6B6B', '#CD9BCD', '#6699CC', 'pink', 
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36', 
  '#9FB1BC', 'springgreen', '#E2C044', '#5BC0EB',
  
  "orange", "#CBD588", "#5F7FC7",  
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",  
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71", 
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", 
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357')


### Bacteria

#### get unassigned vectors
#### get CP and MT phyloseq obj and vectors
# (1) CP
phy.cp <- subset_taxa(phy, Class == "D_2__Chloroplast") ## just confirming 
phy.cp <- subset_taxa(phy, Order == "D_3__Chloroplast")
vec.cp <- rownames(otu_table(phy.cp))
length(rownames(otu_table(phy.cp))) ## 244 otus of CP
vec.cp

# (2) MT
#phy.mt <- subset_taxa(phy, Family == "D_4__Mitochondria")
phy.mt <- subset_taxa(phy, Order == "D_3__Rickettsiales")
vec.mt <- rownames(otu_table(phy.mt))
tax_table(phy.mt)
length(rownames(otu_table(phy.mt))) ## 309 otus of CP -> 66 after chimera filtering

# (3) Unassigned
unique(tax_table(phy)[,'Kingdom']) ## only bacteria, then no need to exclude
phy.un <- subset_taxa(phy, Kingdom == "Unassigned")
vec.un <- rownames(otu_table(phy.un))
tax_table(phy.un)
length(rownames(otu_table(phy.un))) ## 18

### exclude those vectors
phy  # 328 28901 taxa

sample_names(phy)
sample_variables(phy)

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### let's do it!!!

phy.clean <- pop_taxa(phy, c(vec.cp, vec.mt, vec.un))
phy  ## 28901
phy.clean ## 28330
sum(otu_table(phy.clean)) # 18456654
# checking procedure of whether the MT and CP otus are cleaned
taxa_names(subset_taxa(phy, Order == "D_3__Rickettsiales"))  # before 66
taxa_names(subset_taxa(phy, Order=="D_3__Chloroplast")) # before 34

taxa_names(subset_taxa(phy.clean , Order == "D_3__Rickettsiales")) # after 0
taxa_names(subset_taxa(phy.clean , Family=="D_4__Mitochondria")) # after 0
taxa_names(subset_taxa(phy.clean , Order == "D_3__Chloroplast")) # after 0
# phy.clean good to go sir.

sample_data(phy)$Soil_texture
#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(phy.two)[,colnames(tax_table(phy.two))] <- gsub(tax_table(phy.two)[,colnames(tax_table(phy.two))],pattern="[A-Z]_[0-9]__",replacement="")
# phy.test4 <- phy.two %>% psmelt() 
# phy.test4

tax_table(phy.clean)[,colnames(tax_table(phy.clean))] <- gsub(tax_table(phy.clean)[,colnames(tax_table(phy.clean))],pattern="[A-Z]_[0-9]__",replacement="")

#' sample_data(phy.clean)$SampleID <- factor(sample_data(phy.clean)$SampleID, levels =target_PAB)


tax_table(phy.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

phy.clean    ## 28330
str(phy.clean)
otu_table(phy.clean)

phy.clean  ## 28330


# ## fix it in phy.clean object!!! pop_taxa does the work
# phy.clean <- pop_taxa(phy.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(phy.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement 
phy.clean.otu <- otu_table(phy.clean)
head(phy.clean.otu)
df.clean.otu <- data.frame(phy.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


sample_names(phy.clean)
# how about we delete major taxa in the negative sequence?
negative <- df.clean.otu[,c('OTU','N1B', 'N2B', 'N3B', 'DWB')]
negative
negative$total <- apply(negative[,-1], 1, sum)
negative_0 <- subset(negative, negative$total > 0) 
neg.otu <- negative_0$OTU  ### 21 otu to be eliminated
neg.otu
length(neg.otu)
tax_table(phy)[neg.otu]
library(xlsx)
write.xlsx(tax_table(phy)[neg.otu], 'neg.otu_taxonomy.xlsx')

## checking 
## exclude all negative taxa
phy.clean <- pop_taxa(phy.clean, negative_0$OTU)
phy.clean  ### 28178 otu <- 28330
sum(otu_table(phy.clean))  ### 18114239


## Get the samples we need
lb <- c('CC1',
        'CC2',
        'CJ1',
        'CJ2',
        'DG1',
        'DG2',
        'JJ1',
        'JJ2',
        'MY1',
        'MY2',
        'NJ1',
        'NJ2',
        'IS1',
        'IS2',
        'YS1',
        'YS2',
        'UF1',
        'UF2',
        'CC1.18',
        'CC2.18',
        'CJ1.18',
        'CJ2.18',
        'DG1.18',
        'DG2.18',
        'IS1.18',
        'IS2.18',
        'UF1.18',
        'UF2.18')
length(lb)
28*9
phy.clean.ss <- subset_samples(phy.clean, Field %in% lb)
(phy.clean.ss)
phy.clean.ss <- filter_taxa(phy.clean.ss, function(x) sum(x) != 0, TRUE)
phy.clean.ss  ## 23839 otu, 252 samples
sum(otu_table(phy.clean.ss)) # 14500875


## Remove reads with over 300 bp
library(seqinr)
bac.seq <- read.fasta(file = "dna-sequences_bac.fasta", as.string = TRUE, seqtype = "DNA")
bac.seq[which(getLength(bac.seq)>300)]
otu_over_300bp <- attr(bac.seq[which(getLength(bac.seq)>300)], "names")
phy.clean.ss
phy.clean.ss <- pop_taxa(phy.clean.ss,otu_over_300bp)
phy.clean.ss  ## 23762
sum(otu_table(phy.clean.ss)) # 14497979

summarize_phyloseq(phy.clean.ss)
phy.clean.ss


## Divide bacteria and archaea
arch.clean.ss <- subset_taxa(phy.clean.ss, Kingdom == "Archaea")
arch.clean.ss <- phyloseq::filter_taxa(arch.clean.ss, function(x) sum(x) != 0, TRUE)

bac.clean.ss <- subset_taxa(phy.clean.ss, Kingdom != "Archaea")
bac.clean.ss <- phyloseq::filter_taxa(bac.clean.ss, function(x) sum(x) != 0, TRUE)

###Fungal community
#### get unassigned vectors

# (3) Unassigned
unique(tax_table(fun)[,'Kingdom']) ## "k__Chromista", "k__Plantae", "Unassigned"
tax_table(fun)[,'Kingdom']
fun.un <- subset_taxa(fun, Kingdom %in% c("Unassigned","k__Chromista","k__Plantae"))
vec.un <- rownames(otu_table(fun.un))
tax_table(fun.un)
length(rownames(otu_table(fun.un))) ##  116

### exclude those vectors
fun  # 328 samples, 10981 taxa

sample_names(fun)
sample_variables(fun)

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### Clean it up!!!
fun.clean <- pop_taxa(fun, c(vec.un))

#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(fun.two)[,colnames(tax_table(fun.two))] <- gsub(tax_table(fun.two)[,colnames(tax_table(fun.two))],pattern="[A-Z]_[0-9]__",replacement="")
# fun.test4 <- fun.two %>% psmelt() 
# fun.test4
tax_table(fun.clean)[,colnames(tax_table(fun.clean))] <- gsub(tax_table(fun.clean)[,colnames(tax_table(fun.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.clean)$SampleID <- factor(sample_data(fun.clean)$SampleID, levels =target_PAB)

tax_table(fun.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.clean    ## 10865 taxa
str(fun.clean)
otu_table(fun.clean)

fun.clean  ## 10865


# ## fix it in fun.clean object!!! pop_taxa does the work
# fun.clean <- pop_taxa(fun.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement 
fun.clean.otu <- otu_table(fun.clean)
head(fun.clean.otu)
df.clean.otu <- data.frame(fun.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


sample_names(fun.clean)
# how about we delete major taxa in the negative sequence?
negative <- df.clean.otu[,c('OTU','N1F','N2F','N3F','DWF')]
negative
negative$total <- apply(negative[,-1], 1, sum)
negative_0 <- subset(negative, negative$total > 0) 
neg.otu <- negative_0$OTU  ### 58 otu to be eliminated
neg.otu
length(neg.otu)
tax_table(fun.clean)[neg.otu]
library(xlsx)
write.xlsx(tax_table(fun.clean)[neg.otu], 'neg.otu_taxonomy.xlsx')

## checking oats
# oats <- df.clean.otu[,c('OTU','F_Oats1','F_Oats2','F_Oats3')]
# oats$total <- apply(negative[,-1], 1, sum)
# oats_0 <- subset(oats, oats$total > 0) 
# oats_0$OTU
# oats.otu <- oats_0$OTU
# 
# intersect(oats.otu, neg.otu)
# setdiff(oats.otu, neg.otu)
# union(oats.otu, neg.otu)

## exclude all negative taxa
fun.clean <- pop_taxa(fun.clean, negative_0$OTU)
fun.clean  ### 10865 -> 10845 otu
class(fun.clean)
sum(otu_table(fun.clean)) ## 15818852

### 

## Phew!!
## only get the samples of 
lab <- c('CC1',
         'CC2',
         'CJ1',
         'CJ2',
         'DG1',
         'DG2',
         'JJ1',
         'JJ2',
         'MY1',
         'MY2',
         'NJ1',
         'NJ2',
         'IS1',
         'IS2',
         'YS1',
         'YS2',
         'UF1',
         'UF2',
         'CC1.18',
         'CC2.18',
         'CJ1.18',
         'CJ2.18',
         'DG1.18',
         'DG2.18',
         'IS1.18',
         'IS2.18',
         'UF1.18',
         'UF2.18')
length(lab)
43*3
fun.clean.ss <- subset_samples(fun.clean, Field %in% lab)

fun.clean.ss <- filter_taxa(fun.clean.ss, function(x) sum(x) != 0, TRUE)
fun.clean.ss  ## 9450
sum(otu_table(fun.clean.ss)) # 12201972

##### get rid of otu of less than 100 reads

fun.seq <- read.fasta(file = "dna-sequences_fun.fasta", as.string = TRUE, seqtype = "DNA")
fun.seq[which(getLength(fun.seq)<100)]

otu_less_than_100bp <- attr(fun.seq[which(getLength(fun.seq)<100)], "names")
fun.clean.ss
fun.clean.ss <- pop_taxa(fun.clean.ss,otu_less_than_100bp)
fun.clean.ss ## 9437 
sum(otu_table(fun.clean.ss)) ## 12191641



## Designating OTU id
bac.list <- bac.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarize(total=sum(Abundance)) %>% arrange(desc(total))
arch.list <- arch.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarize(total=sum(Abundance)) %>% arrange(desc(total))
fun.list <- fun.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarize(total=sum(Abundance)) %>% arrange(desc(total))

bac.list$number <- paste0('B',1:dim(bac.list)[1])
bac.list

bac.list$OTU_id <- ifelse(is.na(bac.list$Genus),ifelse(is.na(bac.list$Family),paste0(bac.list$number,'_o_',bac.list$Order),paste0(bac.list$number,'_f_',bac.list$Family)),paste0(bac.list$number,'_',bac.list$Genus))
bac.list$OTU_id

arch.list$number <- paste0('A',1:dim(arch.list)[1])
arch.list

arch.list$OTU_id <- ifelse(is.na(arch.list$Genus),ifelse(is.na(arch.list$Family),paste0(arch.list$number,'_o_',arch.list$Order),paste0(arch.list$number,'_f_',arch.list$Family)),paste0(arch.list$number,'_',arch.list$Genus))
arch.list$OTU_id

fun.list$number <- paste0('F',1:dim(fun.list)[1])
fun.list

fun.list$OTU_id <- ifelse(is.na(fun.list$Genus),ifelse(is.na(fun.list$Family),paste0(fun.list$number,'_o_',fun.list$Order),paste0(fun.list$number,'_f_',fun.list$Family)),paste0(fun.list$number,'_',fun.list$Genus))
fun.list$OTU_id

bac.list
arch.list
fun.list

otu.list <- rbind(bac.list, arch.list, fun.list)
dim(otu.list)
# write.xlsx(otu.list,'otu_id_book.xlsx')


OTU_id.list <- rbind(bac.list[c('OTU','OTU_id')],arch.list[c('OTU','OTU_id')],fun.list[c('OTU','OTU_id')])
OTU_id.list$OTU_id
