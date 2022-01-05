library(dplyr)
library(readr)
library(ComplexHeatmap)
library(reticulate)
library(circlize)
library(ggplot2)
library(ggrepel)

#omit N/As and remove gene symbols
both <- py$ab24_RLN
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

#apply log change to dataset & calculate fold change
both <- log2(both)
control <- apply(both[,1:4],1,mean)
test <- apply(both[,5:8],1,mean)
foldchange <- control- test

#bind p val and log change data
results = cbind(foldchange,rawpvalue)
results = as.data.frame(results)
results$log10 <- -log10(rawpvalue)
results$probename <- py$ab1wk_S[,1]

#Volcano Plot
volcano = ggplot(data = results,(aes( x = foldchange, y = -log10(rawpvalue))))

p <- volcano + geom_point() + theme_classic()

p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")

#add a column of NAs
results$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
results$diffexpressed[results$foldchange > 0.6 & results$rawpvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
results$diffexpressed[results$foldchange < -0.6 & results$rawpvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=results, aes(x=foldchange, y=-log10(rawpvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")

p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: create a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

results$delabel <- NA
results$delabel[results$diffexpressed != "NO"] <- results$probename[results$diffexpressed != "NO"]

# plot adding up all layers we have seen so far
 ggplot(data=results, aes(x=foldchange, y=-log10(rawpvalue), col=diffexpressed, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.6, 0.6), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")
