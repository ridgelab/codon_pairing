#Create a plot for the frequency of codong pairing by codong

#install.packages("ggplot2")

<<<<<<< HEAD
args <- commandArgs(TRUE)
inputFile <- args[1]
outputFile <- args[2]
graphTitle <- args[3]
if (length(args) != 3) {
print("Please supply 3 arguments: an input file, an output file, and the name of the clade")
}else {

library(ggplot2)

 myData <- read.csv(file=inputFile,header=TRUE,sep=",");  

 png(outputFile, width = 1400, height = 600) 
=======
library(ggplot2)

 myData <- read.csv(file="./viruses_codonsUsed.txt",header=TRUE,sep=",");  #insert the path here

 png("viruses_codonsUsed.png", width = 1400, height = 600)
>>>>>>> 1cb174ba31a54f513c839296f90651f9dc93da67

 plot <- ggplot(data = myData, aes(x = species_codon, y = frequency)) + ylim(0,1.0) +
   theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
   theme(axis.text = element_text(size = 18)) +
   theme(axis.title = element_text(size=24)) +
   theme(plot.title = element_text(size=28,hjust=0.5)) +
   geom_boxplot(outlier.size = 0.2)

<<<<<<< HEAD
 print(plot + ggtitle(paste(graphTitle, " - Frequency of Codon Pairing by Codon")) + labs(x = "Codon", y = "Frequency")) #Change the title of the graph here

dev.off()
}
=======
 print(plot + ggtitle("Viruses - Frequency of Codon Pairing by Codon") + labs(x = "Codon", y = "Frequency"))

dev.off()

>>>>>>> 1cb174ba31a54f513c839296f90651f9dc93da67
