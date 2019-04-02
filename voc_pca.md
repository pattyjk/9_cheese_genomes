## Analysis of VOC data vai PCA
```
library(ggplot2)
library(vegan)

#read in data
voc_rel_peak<-read.delim("9_cheese_genomes/VOC/VOC_rel_peak_area.txt", header=T, row.names=1)

#transpose data
voc_t<-t(voc_rel_peak)

#calculate a PCA
voc_pca<-rda(voc_t)

#extract coordinates to plot
voc.scores<-scores(voc_pca)
voc.scores<-as.data.frame(voc.scores$sites)
plot(voc.scores$PC1, voc.scores$PC2)

#make metadata
voc.scores$community<-row.names(voc.scores)
voc.scores$community<-gsub("X", "", voc.scores$community)
voc.scores$community<-gsub("_[[:digit:]]+", "", voc.scores$community)

#plot data
pally<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
ggplot(voc.scores, aes(PC1, PC2, colour=community))+
  geom_point(aes(size=2))+
  theme_bw()+
  scale_colour_manual(values=pally)

#extract VOC (species) scores
voc.sp.scores<-as.data.frame(scores(voc_pca, display = "species"))
voc.sp.scores$VOC<-row.names(voc.sp.scores)

#basic plot
ggplot(voc.sp.scores, aes(VOC, PC1))+
  geom_point()+
  theme_bw()+
  coord_flip()

write.table(voc.sp.scores, "VOC_sp_scores.txt", row.names=F, quote=F, sep="\t")

#get eigen vales
voc.summary<-summary(voc_pca, scaling = 2)
voc.eigen<-voc.summary$cont$importance

write.table(voc.eigen, "VOC_eign_values.txt", row.names=T, quote=F, sep="\t")

```


