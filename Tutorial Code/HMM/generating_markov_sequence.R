#A Markov model of DNA sequence evolution

# The transition matrix for a Markov model

# Define the alphabet of nucleotides
nucleotides <- c("A", "C", "G", "T")
# Set the values of the probabilities, where the previous
afterAprobs <- c(0.2, 0.3, 0.3, 0.2)
# Set the values of the probabilities, where the previous
afterCprobs <- c(0.1, 0.41, 0.39, 0.1)
# Set the values of the probabilities, where the previous
afterGprobs <- c(0.25, 0.25, 0.25, 0.25)
# Set the values of the probabilities, where the previous
afterTprobs <- c(0.5, 0.17, 0.17, 0.17)

mergedprobs = c(afterAprobs,afterCprobs,afterGprobs,afterTprobs)

mytransitionmatrix <- matrix(c(afterAprobs, afterCprobs, afterGprobs, afterTprobs), 4, 4, byrow=TRUE)
rownames(mytransitionmatrix) <- nucleotides
colnames(mytransitionmatrix) <- nucleotides
# Print out the transition matrix
mytransitionmatrix

# Generating a DNA sequence using a Markov model

generatemarkovseq <- function(transitionmatrix, firstnucleotide, seqlength)
{
  nucleotides <- c("A", "C", "G", "T")
  mysequence <- character()
  
  mysequence[1] <- firstnucleotide
  for (i in 2:seqlength)
  {
    prevnucleotide <- mysequence[i-1]
    probabilities <- transitionmatrix[prevnucleotide,]
    nucleotide <- sample(nucleotides, 1, prob=probabilities)
    mysequence[i] <- nucleotide 
  }
  return(mysequence)
}

a =generatemarkovseq(mytransitionmatrix, 'A', 30)
