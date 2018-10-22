#Create a plot for number of times a motif is repeated

<<<<<<< HEAD
=======

>>>>>>> 1cb174ba31a54f513c839296f90651f9dc93da67
#install.packages("ggplot2")

library(ggplot2)

<<<<<<< HEAD
args <- commandArgs(TRUE)
inputFile <- args[1]
outputFile <- args[2]
graphTitle <- args[3]
if (length(args) != 3) {
print("Please supply 3 arguments: an input file, an output file, and the name of the clade")
}else {

 myData <- read.csv(file=inputFile, header = TRUE,sep=","); #Insert path to input file here

  png(outputFile, width = 1400, height = 600)#Insert path to output graph here

 motifPlot <- ggplot(data = myData) + geom_point(mapping = aes(x = motif_frequency, y = count), color = "black", size = 2) + #xlim(0,120) + ##Add this to set an x-limit, removing outliers
   theme(axis.text = element_text(size = 18)) +
   theme(axis.title = element_text(size=24)) +
   theme(plot.title = element_text(size=28,hjust=0.5))
 print(motifPlot + ggtitle(paste(graphTitle, " - Repeated Motifs")) + labs(x = "Number of Times a Motif is Repeated", y = "ln(frequency)")) #Adjust this line to change the title

dev.off()

}
=======

 myData <- read.csv(file="./viruses_repeatedMotifs.txt", header = TRUE,sep=",");

  png("./viruses_repeatedMotifs.png", width = 1400, height = 600)

 motifPlot <- ggplot(data = myData) + geom_point(mapping = aes(x = motif_frequency, y = count), color = "black", size = 2) + #xlim(0,120) +
   theme(axis.text = element_text(size = 18)) +
   theme(axis.title = element_text(size=24)) +
   theme(plot.title = element_text(size=28,hjust=0.5))
 print(motifPlot + ggtitle("Viruses - Repeated Motifs") + labs(x = "Number of Times a Motif is Repeated", y = "ln(frequency)"))

dev.off()


>>>>>>> 1cb174ba31a54f513c839296f90651f9dc93da67

