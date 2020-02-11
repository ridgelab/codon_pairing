library(readr)
library(ggplot2)

args <- commandArgs(TRUE)
inputFile <- args[1]
outputFile <- args[2]
graphTitle <- args[3]
if (length(args) != 3) {
print("Please supply 3 arguments: an input file, an output file, and the name of the clade")
}else {

perms <- read_csv(inputFile)

g <- ggplot(perms,aes(x = meanRI))+
  geom_histogram(bins = 80) +
  theme_bw() +
  theme(text = element_text(size=25),axis.text.x = element_text(angle = 0, hjust = 1),legend.title = element_blank(),plot.title = element_text(size=22, hjust=0.5))  +
  labs(x = "Mean Retention Index", y = "Frequency") +
  ggtitle(paste(c("Retention Index of",graphTitle,"Ramps\nCompared to Random Permutations"))+
  geom_vline(xintercept = 0.278,colour = "black",linetype = "dashed",size = 1.3))

ggsave(filename=outputFile,plot=g,dpi=300)
}

