#!/usr/bin/env Rscript
library(ggplot2)

#ARG_genes

dataARG<-read.csv("/Volumes/TOSHIIBA/TFM BIOINFORMATICA/ARG_R_clean.csv", header = TRUE, sep="\t")

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

ggsave("heatmapARG_gene.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')


#ARG_family

mine.heatmap<-ggplot(data=dataARG, mapping = aes(x=ARG.family,
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

ggsave("heatmapARG_family.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')





#PLASMIDS

dataARG<-read.csv("/Volumes/TOSHIIBA/TFM BIOINFORMATICA/plasmid_R_clean.csv", header = TRUE, sep="\t")

LOG.normalizedabun<-log10(dataARG$NormalizedSum)
mine.heatmap<-ggplot(data=dataARG, mapping = aes(x=EVENT,
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

ggsave("heatmapPLASMID.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')


#INTEGRASES

dataARG<-read.csv("/Volumes/TOSHIIBA/TFM BIOINFORMATICA/integrases_R_clean.csv", header = TRUE, sep="\t")

LOG.normalizedabun<-log10(dataARG$NormalizedSum)
mine.heatmap<-ggplot(data=dataARG, mapping = aes(x=EVENT,
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

ggsave("heatmapIntegrases.tiff", units="in", width=30, height=30, dpi=300, compression = 'lzw')
