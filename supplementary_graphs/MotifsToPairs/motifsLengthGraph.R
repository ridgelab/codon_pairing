#!/usr/bin/env Rscript
#Input 1-3: cotrna, dna, and combined files
#Input 4: name of taxonomic group
#Input 5: name of output pdf file

args = commandArgs(trailingOnly=TRUE)

# Load packages (To install, run install.packages("ggplot2"), etc.
library(ggplot2)
library(readr)
library(grid)
library(gridExtra)


remove_outliers <- function(myFile){
  #Identify and remove outliers for Gene Length
  lower_lim = quantile(myFile$Gene_length,.25,na.rm=TRUE) - (1.5 * IQR(myFile$Gene_length, na.rm = TRUE))  
  upper_lim = quantile(myFile$Gene_length,.75,na.rm=TRUE) + (1.5 * IQR(myFile$Gene_length,na.rm = TRUE))
  filtered1 <- subset(myFile,Gene_length > lower_lim & Gene_length < upper_lim)
  
  #Identify and remove outliers for Num Pairs
  lower_lim2 = quantile(myFile$Num_Pairs,.25,na.rm=TRUE) - (1.5 * IQR(myFile$Num_Pairs,na.rm = TRUE))  
  upper_lim2 = quantile(myFile$Num_Pairs,.75,na.rm=TRUE) + (1.5 * IQR(myFile$Num_Pairs,na.rm = TRUE))
  filtered <- subset(filtered1,Num_Pairs > lower_lim2 & Num_Pairs < upper_lim2)
  return(filtered)
}


# DATA: CSV files of pairing to length for co-tRNA, identical codon, and combined pairing

coFile <- read_csv(args[1])
pairingFile <- read_csv(args[2])
combFile <- read_csv(args[3])

#Add the pairing type column
coFile$Type = rep("A.  Co-tRNA",nrow(coFile))
pairingFile$Type = rep("B.  Identical Codon",nrow(pairingFile))
combFile$Type = rep("C.  Combined",nrow(combFile))

# Calculate linear models
co_linearMod <- lm(Num_Pairs ~ Gene_length, data=coFile)
co_result = summary(co_linearMod)
co_r_squared = co_result$r.squared

pairing_linearMod <- lm(Num_Pairs ~ Gene_length, data=pairingFile)
pairing_result = summary(pairing_linearMod)
pairing_r_squared = pairing_result$r.squared

comb_linearMod <- lm(Num_Pairs ~ Gene_length, data=combFile)
comb_result = summary(comb_linearMod)
comb_r_squared = comb_result$r.squared

# Remove outliers
co_filtered = remove_outliers(coFile)
comb_filtered = remove_outliers(combFile)
pairing_filtered = remove_outliers(pairingFile)

# Combine into one dataframe
bigFile = rbind(co_filtered,pairing_filtered,comb_filtered)

# Subset to a reasonable size
rand.rows <- sample(nrow(bigFile),5000)
smallFile <- bigFile[rand.rows,]

# Make the nice plot
l1 <- paste("Co-tRNA R-squared = ",round(co_r_squared,3))
l2 <- paste("Identical Codon R-squared = ",round(pairing_r_squared,3))
l3 <- paste("Combined R-squared = ",round(comb_r_squared,3))
labels.minor <- paste(l1,l2,l3,sep= "\n")      

regression_data = data.frame(Type = c("A.  Co-tRNA","B.  Identical Codon","C.  Combined"),  slopes = c(co_linearMod$coefficients[[2]],pairing_linearMod$coefficients[[2]],comb_linearMod$coefficients[[2]]), intercepts = c(co_linearMod$coefficients[[1]],pairing_linearMod$coefficients[[1]],comb_linearMod$coefficients[[1]]))

g <- ggplot(data = smallFile, aes(x = Gene_length, y = Num_Pairs)) +
  facet_wrap(~Type) + 
  geom_point(size = 1,alpha = 0.5) +
  geom_abline(data = regression_data,aes(slope = slopes, intercept = intercepts)) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  ggtitle(paste(args[4],": Number of Codons in Pairing Motif vs. Gene Length",sep = "")) +
  labs(x = "Gene Length (Number of Codons)", y = "Number of Codon Pairs") +
  theme_bw() + 
  theme(plot.title = element_text(size=14,hjust=0.5)) +
  theme(axis.text.y = element_text(size = 9)) +
  theme(axis.text.x = element_text(size = 9)) +
  theme(axis.title = element_text(size = 12))


g2 <- grid.arrange(g,bottom = textGrob(labels.minor, x = 0.55, y = .6, gp = gpar(fontsize = 9)))

#Save the plot
ggsave(filename = args[5], plot = g2,  dpi = 300)
