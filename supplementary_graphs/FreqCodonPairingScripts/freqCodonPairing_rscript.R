#Create a plot for the frequency of codong pairing by codon

# Uncomment these lines to install needed packages
#install.packages("ggplot2")
#install.packages("gridExtra")
#install.packages("grid")

args <- commandArgs(TRUE)
inputFile <- args[1]
outputFile <- args[2]
graphTitle <- args[3]
if (length(args) != 3) {
print("Please supply 3 arguments: an input file, an output file, and the name of the clade")
}else {

library(ggplot2)
library(gridExtra)
library(grid)

# Load the data
myData <- read.csv(file=inputFile,header=TRUE,sep=",");  
# Replace the * symbol for Stop with "Stop"
myData$amino_acid = gsub("\\*","Stop",myData$amino_acid)

# Order the codons in the order of the amino acids
myData$codon = factor(myData$codon, levels = c("ATT","ATC","ATA","CTT","CTC","CTA","CTG","TTA","TTG","GTT","GTC","GTA","GTG","TTT","TTC","ATG","TGT","TGC","GCT","GCC","GCA","GCG","GGT","GGC","GGA","GGG","CCT","CCC","CCA","CCG","ACT","ACC","ACA","ACG","TCT","TCC","TCA","TCG","AGT","AGC","TAT","TAC","TGG","CAA","CAG","AAT","AAC","CAT","CAC","GAA","GAG","GAT","GAC","AAA","AAG","CGT","CGC","CGA","CGG","AGA","AGG","TAA","TAG","TGA"))
  
# Make the figure                      
plot <- ggplot(data = myData, aes(x = codon, y = frequency)) + 
  ylim(-0.05,1.0) +
  theme(axis.text.x = element_text(angle = 50,vjust=0.6)) +
  theme(axis.text = element_text(size = 10)) +
  theme(axis.title = element_text(size=18)) +
  theme(plot.title = element_text(size=28,hjust=0.5)) +
  geom_boxplot(outlier.size = 0.2) +
  ggtitle(paste(graphTitle, " - Frequency of Codon Pairing by Codon")) + 
  labs(x = "Codon", y = "Frequency")

# Add the amino acid label
agg <- aggregate(frequency ~ amino_acid + codon,data = myData, min)
agg$frequency = rep(-0.05,nrow(agg))
plot2 <- plot + geom_text(data = agg, aes(label = amino_acid,vjust = 0,size = 10))
print(plot2)
plot3 <- grid.arrange(plot2,bottom = textGrob("Amino Acid", x = 0.5, y = 6.5, gp = gpar(fontsize = 18)))

# Save the figure
ggsave(filename = outputFile, plot = plot3, dpi = 300)

}
