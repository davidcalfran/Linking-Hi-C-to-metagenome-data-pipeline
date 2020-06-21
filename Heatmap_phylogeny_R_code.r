#!/usr/bin/env Rscript
library(ggplot2)
library(taxize)
library(dplyr)

# Documento del script ARG


dataARG<-read.csv("./ARG_R_clean.csv", header = TRUE, sep="\t")
dataARG<-subset(dataARG,Number.of.Contigs>=100) #Modify, depends on case.


LOG.normalizedabun<-log10(dataARG$NormalizedSum)
mine.heatmap<-ggplot(data=dataARG, mapping = aes(x=ARG,
                                              y=SPECIES,
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "Log(Normalized sum of interactions)",
                      low = "#003366",
                      high = "#EFEDC2") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=20, hjust=1),
    axis.text.y = element_text(size=20, hjust=1)
  )


mine.heatmap

ggsave("heatmapARGtest.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')


#Taxize was used for retrieve the species phylum and genus names from NCBI

names<-as.character(dataARG$SPECIES)
class(names)
test<-tax_name(query=c(names), get= c("phylum","genus"), db="ncbi")


#Subsetting genus and phylum to be added to the original output file


ARG_genus<-test$genus
ARG_phylum<-test$phylum
dataARG1<-cbind(dataARG, ARG_genus, ARG_phylum)
dataARG1


LOG.normalizedabun<-log10(dataARG1$NormalizedSum)
mine.heatmap<-ggplot(data=dataARG1, mapping = aes(x=ARG,
                                              y=SPECIES,  #This can be modified by changin SPECIES to Genus or Phylum depending on the level of detail needed.
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "Log(Normalized sum of interactions)",
                      low = "#003366",
                      high = "#EFEDC2") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=20, hjust=1),
    axis.text.y = element_text(size=20, hjust=1)
  )


mine.heatmap

ggsave("heatmapARG1test.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')


#Remove duplicates for tree

library(dplyr)
names<-as.data.frame(dataARG1$SPECIES)
names<-distinct(names)

names<-as.character(names[,1])


##Generate the classification vector

tree_ARG<-classification(names, db="ncbi")

## Generat tree

tree_ARG_real<-class2tree(tree_ARG, check=TRUE, varstep = TRUE)

tiff(file="tree_ARG.tiff", width=15, height=30, units="in", res=300)

plot(tree_ARG_real)

dev.off()


tree_ARG_real$phylo$tip.label


## Obtain a vector with the species order of appearance in the tree


order_ARG<-rev(tree_ARG_real$phylo$tip.label[tree_ARG_real$phylo$edge[tree_ARG_real$phylo$edge[,2] <= 91,2]]) #This 91 has to be adjusted depending on the number of entries. Case-dependent. Manually.
order_ARG

df<-data.frame(dataARG1)
df_corrected<-gsub("\\b[sp]{1,2}\\b", "sp.", df$SPECIES) #This will add a point after every sp
dataARG2<-cbind(df,df_corrected)
dataARG2

# Reverse order to fit the way the tree is generated
level_order<-rev(order_ARG)
level_order


#Representation of the heatmap in order of appearance for linking it to the phylogenetic tree

LOG.normalizedabun<-log10(dataARG2$NormalizedSum)
mine.heatmap<-ggplot(data=dataARG2, mapping = aes(x=ARG,
                                              y=factor(df_corrected, level=level_order),
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "",
                      low = "#003366",
                      high = "#FFD700",
                      ) +

  labs(x="Antibiotic Resistance Genes (ARG)", y="Species") +

  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=20, hjust=1),
    axis.text.y = element_text(size=30, hjust=1),
    legend.key.size = unit(1,"in"),
    legend.text = element_text(size=30),
    axis.title  = element_text(size=30)
  )

mine.heatmap

ggsave("heatmapARG_ordeded.tiff", units="in", width=40, height=30, dpi=300, compression = 'lzw')


# Documento del script plasmid

dataplasmid<-read.csv("./plasmid_R_clean.csv", header = TRUE, sep="\t")


LOG.normalizedabun<-log10(dataplasmid$NormalizedSum)
mine.heatmap<-ggplot(data=dataplasmid, mapping = aes(x=EVENT,
                                              y=SPECIES,
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "Log(Normalized sum of interactions)",
                      low = "#003366",
                      high = "#EFEDC2") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=20, hjust=1),
    axis.text.y = element_text(size=20, hjust=1)
  )


mine.heatmap

ggsave("heatmapplasmid.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')

# Taxize. Retrieve names with phylum and genus from NCBI

names<-as.character(dataplasmid$SPECIES)
class(names)
test<-tax_name(query=c(names), get= c("phylum","genus"), db="ncbi")

# Add genus and phylum to original file

plasmid_genus<-test$genus
plasmid_phylum<-test$phylum
dataplasmid1<-cbind(dataplasmid, plasmid_genus, plasmid_phylum)

# Heatmap test

LOG.normalizedabun<-log10(dataplasmid1$NormalizedSum)
mine.heatmap<-ggplot(data=dataplasmid1, mapping = aes(x=EVENT,
                                              y=SPECIES,
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "Log(Normalized sum of interactions)",
                      low = "#003366",
                      high = "#EFEDC2") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=20, hjust=1),
    axis.text.y = element_text(size=20, hjust=1)
  )


mine.heatmap

ggsave("heatmapplasmid1.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')

#Remove duplicates for tree

names<-as.data.frame(dataplasmid1$SPECIES)
names<-distinct(names)
names<-as.character(names[,1])


##Generate the classification vector

tree_plasmid<-classification(names, db="ncbi")

## Generat tree

tree_plasmid_real<-class2tree(tree_plasmid, check=TRUE, varstep = TRUE)

tiff(file="tree_plasmid.tiff", width=7, height=15, units="in", res=300)

plot(tree_plasmid_real)

dev.off()


## Obtain a vector with the species order of appearance in the tree

order_plasmid<-rev(tree_plasmid_real$phylo$tip.label[tree_plasmid_real$phylo$edge[tree_plasmid_real$phylo$edge[,2] <= 12,2]])

order_plasmid

df<-data.frame(dataplasmid1)
df_corrected<-gsub("\\b[sp]{1,2}\\b", "sp.", df$SPECIES) #This will add a point after every sp
dataplasmid2<-cbind(df,df_corrected)
level_order<-rev(order_plasmid)



LOG.normalizedabun<-log10(dataplasmid2$NormalizedSum)
mine.heatmap<-ggplot(data=dataplasmid2, mapping = aes(x=EVENT,
                                              y=factor(df_corrected, level=level_order),
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "",
                      low = "#003366",
                      high = "#FFD700",
                      ) +

  labs(x="Plasmid", y="Species") +

  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=30, hjust=1),
    axis.text.y = element_text(size=30, hjust=1),
    legend.key.size = unit(1,"in"),
    legend.text = element_text(size=30),
    axis.title  = element_text(size=30)
  )

mine.heatmap

ggsave("heatmapplasmid_ordered.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')


# Documento del script integrones


dataint<-read.csv("./integrases_R_clean.csv", header = TRUE, sep="\t")

# Test heatmap

LOG.normalizedabun<-log10(dataint$NormalizedSum)
mine.heatmap<-ggplot(data=dataint, mapping = aes(x=EVENT,
                                              y=SPECIES,
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "Log(Normalized sum of interactions)",
                      low = "#003366",
                      high = "#EFEDC2") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=20, hjust=1),
    axis.text.y = element_text(size=20, hjust=1)
  )


mine.heatmap

ggsave("heatmapint.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')

# Retrieve species names with phylum and genus from NCBI

names<-as.character(dataint$SPECIES)
class(names)
test<-tax_name(query=c(names), get= c("phylum","genus"), db="ncbi")

# Add phylum and genus info to initial document

int_genus<-test$genus
int_phylum<-test$phylum
dataint1<-cbind(dataint, int_genus, int_phylum)

library(ggplot2)
LOG.normalizedabun<-log10(dataint1$NormalizedSum)
mine.heatmap<-ggplot(data=dataint1, mapping = aes(x=EVENT,
                                              y=SPECIES,
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "Log(Normalized sum of interactions)",
                      low = "#003366",
                      high = "#EFEDC2") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=20, hjust=1),
    axis.text.y = element_text(size=20, hjust=1)
  )


mine.heatmap

ggsave("heatmapint1.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')


#Remove duplicates for tree

names<-as.data.frame(dataint1$SPECIES)
names<-distinct(names)
names<-as.character(names[,1])


##Generate the classification vector

tree_int<-classification(names, db="ncbi")

## Generat tree


tree_int_real<-class2tree(tree_int, check=TRUE, varstep = TRUE)

tiff(file="tree_int.tiff", width=7, height=15, units="in", res=300)

plot(tree_int_real)

dev.off()


## Obtain a vector with the species order of appearance in the tree


order_int<-rev(tree_int_real$phylo$tip.label[tree_int_real$phylo$edge[tree_int_real$phylo$edge[,2] <= 11,2]])


df<-data.frame(dataint1)
df_corrected<-gsub("\\b[sp]{1,2}\\b", "sp.", df$SPECIES) #This will add a point after every sp
dataint2<-cbind(df,df_corrected)
level_order<-rev(order_int)

#Ordered heatmap

LOG.normalizedabun<-log10(dataint2$NormalizedSum)
mine.heatmap<-ggplot(data=dataint2, mapping = aes(x=EVENT,
                                              y=factor(df_corrected, level=level_order),
                                              fill=LOG.normalizedabun)) +
  geom_tile() +
  scale_fill_gradient(name = "",
                      low = "#003366",
                      high = "#FFD700",
                      ) +

  labs(x="Integrons", y="Species") +

  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle=45, size=30, hjust=1),
    axis.text.y = element_text(size=30, hjust=1),
    legend.key.size = unit(1,"in"),
    legend.text = element_text(size=30),
    axis.title  = element_text(size=30)
  )

mine.heatmap

ggsave("heatmapint_ordered.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')
