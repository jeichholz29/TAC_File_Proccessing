library(dplyr)
library(readr)
library(ComplexHeatmap)
library(reticulate)
library(circlize)
library(ggplot2)
library(ggrepel)

#omit N/As and remove gene symbols
both <- py$df
both <- select(both, -GeneSymbols)
both <- na.omit(both)

#apply p-values from t-test
ttestRat <- function(df,grp1,grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)
  results = t.test(x,y)
  results$p.value
}
rawpvalue = apply(both,1,ttestRat,grp1 = c(1:4), grp2 = c(5:8))

#Normalize dataset & calculate fold change
both <- log2(both)
control <- apply(both[,1:4],1,mean)
test <- apply(both[,5:8],1,mean)
foldchange <- control- test

#bind p val and log change data
results = cbind(foldchange,rawpvalue)
results = as.data.frame(results)
results$log10 <- -log10(rawpvalue)
results$probename <- py$ab1wk_S[,1]

#add columns showing whether genes were differentially expressed
results$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
results$diffexpressed[results$foldchange > 0.6 & results$rawpvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
results$diffexpressed[results$foldchange < -0.6 & results$rawpvalue < 0.05] <- "DOWN"

#set colors to diffexpressed values
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

#create a column so only non-N/A Values are shown in volcano plot
results$delabel <- NA
results$delabel[results$diffexpressed != "NO"] <- results$probename[results$diffexpressed != "NO"]

#plot ggplot object so final figure loads faster
ggplot(data=results, aes(x=foldchange, y=-log10(rawpvalue), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text()

# plot adding up all layers we have seen so far
 ggplot(data=results, aes(x=foldchange, y=-log10(rawpvalue), col=diffexpressed, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.6, 0.6), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")
