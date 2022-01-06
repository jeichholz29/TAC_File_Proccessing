#import libraries
library(dplyr)
library(readr)
library(ComplexHeatmap)
library(circlize)

#Load csv
B <- read_csv('group_b.csv')

#Remove ID var and select genes of interest
BG <- B %>%
  select("KRAS", "EGFR", "STK11", "BRAF")

#apply rownames
colnames(BG) <- c("KRAS", "EGFR", "STK11", "BRAF")
rownames(BG) <- B$ID

#create var to store colors
col_fun = colorRamp2(c(0, 5, 10), c("dodgerblue", "lightyellow", " red"))

#print Heatmap Object
Heatmap(BG,
        name = "Expression",
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = ("left"),
        column_names_side =("bottom"),
        col = col_fun,
        row_order = sort(rownames(BG)))

