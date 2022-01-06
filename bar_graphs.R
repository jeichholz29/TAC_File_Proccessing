#import libraries 
library(readxl)
library(ggsci)
library(ggplot2)
library(gridExtra)

#load files
altb <- read_excel('Downloads/TAC_Files/alt_bar.xlsx')
eltb <- read_excel('Downloads/TAC_Files/elt_bar.xlsx')

#bar graph for group A
agb <- ggplot(data = altb, aes(x = Genes, y = expr, fill = ID)) + 
    geom_bar(stat="identity", position = "dodge")+ 
    labs(x =  "Genes", y = "Expression", title = "A Left Tumor (All Timepoints)") + 
  theme_classic() + 
  scale_fill_igv()

#bar graph for group E
egb <- ggplot(data = eltb, aes(x = genes, y = expr, fill = ID)) + 
  geom_bar(stat="identity", position = "dodge")  + 
  labs(x =  "Genes", y = "Expression", title = "E Left Tumor (All Timepoints)") + 
    theme_classic() + 
  scale_fill_igv()

#format graphs to be together
ga <- grid.arrange(agb,egb)
