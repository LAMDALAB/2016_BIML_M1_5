states <- c("AT-rich", "GC-rich") 
ATrichprobs <- c(0.7, 0.3) 
GCrichprobs <- c(0.1, 0.9) 
thetransitionmatrix <- matrix(c(ATrichprobs, GCrichprobs), 2, 2, byrow = TRUE) 
rownames(thetransitionmatrix) <- states
colnames(thetransitionmatrix) <- states
thetransitionmatrix


nucleotides <- c("A", "C", "G", "T") 
ATrichstateprobs <- c(0.39, 0.1, 0.1, 0.41) 
GCrichstateprobs <- c(0.1, 0.41, 0.39, 0.1) 
theemissionmatrix <- matrix(c(ATrichstateprobs, GCrichstateprobs), 2, 4, byrow = TRUE)
rownames(theemissionmatrix) <- states
colnames(theemissionmatrix) <- nucleotides
theemissionmatrix 



viterbi <- function(sequence, transitionmatrix, emissionmatrix)
{
  # Get the names of the states in the HMM:
  states <- rownames(theemissionmatrix)
  
  # Make the Viterbi matrix v:
  v <- makeViterbimat(sequence, transitionmatrix, emissionmatrix)
  
  for( i in 1:(dim(v))[1]){
    print(v[i,])
  }
}

makeViterbimat <- function(sequence, transitionmatrix, emissionmatrix)
{
  v <- matrix(NA, nrow = length(sequence), ncol = dim(transitionmatrix)[1])
  colnames(v) = c('AT-rich','GC-rich')
  
  v[1,'AT-rich'] = 1
  v[1,'GC-rich'] = 0
  
  for (i in 2:length(sequence)) # For each position in the DNA sequence:
  {
    statelprobnucleotidei <- emissionmatrix[1,sequence[i]]
    v[i,1] <-  statelprobnucleotidei * max(v[(i-1),] * transitionmatrix[,1])
      
    statelprobnucleotidei = emissionmatrix[2,sequence[i]]
    v[i,2] <-  statelprobnucleotidei * max(v[(i-1),] * transitionmatrix[,2])
  
  }
  return(v)
}

myseq = c("A", "A", "G", "C", "G", "T", "G", "G", "G", "G", "C", "C", "C", "C", "G", "G", "C", "G", "A", "C", "A")
myseq = c("A","A","G","C","G","T","G","G","T","G","C","T","A","C","G","G","C")
myseq = c("A","T","A","T","T","T","A","T","A","T","A","A","A","T","A","A","T","T")
myseq = c("C","G","G","C","C","G","G","C","G","C","G","C","G","C","G","C","G","C","G")
viterbi(myseq, thetransitionmatrix, theemissionmatrix)
