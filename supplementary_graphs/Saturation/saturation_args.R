#This script makes a graph of the saturation of pairing distances
#Input 1-3: cotrna, dna, and combined matrix files
#Input 4: Name of taxonomic group
#Input 5: Output pdf file

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
print("Please supply 5 arguments: 3 nput files of Co-tRNA, DNA, and combined pairing, the name of the clade, and an output file.")
}else {

library(readr)
library(ggplot2)
library(tidyr)


# DATA for a taxonomic group. A csv file for co-tRNA, dna, and pairing
aa_matrix = read_csv(args[1])
dna_matrix = read_csv(args[2])
comb_matrix = read_csv(args[3])

process_file <- function(df){
  tempNames = df$X1  #[1:5]
  df$X1 = NULL
  
  #Sort the row values in increasing order. This also transforms the matrix for some reason.
  df = apply(df, 1, sort)
  df= as.data.frame(df)
  colnames(df) = tempNames
  
  #The index keeps track of what species comparison it is.
  df$Index = seq_along(df[,1])
  
  df.g = gather(df,key = "Species", value = "value",-Index)
  
  # #Remove 0 values (where index is 1) then decrement the index
  df.n = df.g[df.g$Index != 1,]
  df.n$Index = df.n$Index - 1
  
  matrixMelted <- reshape2::melt(df.n, id.vars = c("Index", "Species"))
  head(matrixMelted)
  matrixMelted$Index = as.numeric(matrixMelted$Index)
  
  return(matrixMelted)
}

aa_melted = process_file(aa_matrix)
aa_melted$Type = rep("Co-tRNA", nrow(aa_melted))
dna_melted = process_file(dna_matrix)
dna_melted$Type = rep("Identical Codon", nrow(dna_melted))
comb_melted = process_file(comb_matrix)
comb_melted$Type = rep("Combined", nrow(comb_melted))

big_melted = rbind(aa_melted,dna_melted,comb_melted)

# Subset to a reasonable size
rand.rows <- sample(nrow(big_melted),min(c(nrow(big_melted),40000)))
smallFile <- big_melted[rand.rows,]

# Make a line plot of the saturation. Each line is a species, and the index is the closest neighbor
g <- ggplot(smallFile, aes(x = Index, y = value, group = Species)) +
  facet_wrap(~Type) +
  geom_line() +
  theme(legend.position = "none") +
    ggtitle(paste("Saturation of",args[4],"Pairwise Distances")) +
    labs(x = "Taxonomic Distance to Next Species", y = "Computed Distance") +
  theme_bw() + 
    theme(plot.title = element_text(size=18,hjust=0.5)) +
    theme(axis.text = element_text(size = 12)) + 
  theme(axis.title = element_text(size = 14)) 

#Save the figure
ggsave(filename = args[5],plot=g, dpi = 300)

}
