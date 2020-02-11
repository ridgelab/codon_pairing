#Plot for number of codon pairs in each  motif

#Uncomment this line to install needed package
#install.packages("ggplot2")

library(ggplot2)

args <- commandArgs(TRUE)
inputFile <- args[1]
outputFile <- args[2]
graphTitle <- args[3]
if (length(args) != 3) {
print("Please supply 3 arguments: an input file, an output file, and the name of the clade")
}else {

myData <- read.csv(file=inputFile,header=TRUE,sep=","); #Insert path to input file here

codonPlot <- ggplot(data = myData) +
 geom_point(mapping = aes(x = Codons_excluded, y = frequency), color = "black") +
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title = element_text(size=24)) +
  theme(plot.title = element_text(size=28,hjust=0.5))

codonPlot_2 <- codonPlot + ggtitle(paste(graphTitle," - Number of Codon Pairs in Each Motif")) + labs(y = "ln(frequency)", x = "Number of Codon Pairs") #Modify title to the graph here
ggsave(filename = outputFile, plot = codonPlot_2, dpi = 300)#Insert path to output graph here
}
