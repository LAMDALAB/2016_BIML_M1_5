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

generatemarkovseq <- function(transitionmatrix, initial_nucleotide, seqlength)
{
  nucleotides <- c("A", "C", "G", "T")
  mysequence <- character()
  firstnucleotide = initial_nucleotide
  
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

myinitialprobs <- c(0.25, 0.25, 0.25, 0.25)
generatemarkovseq(mytransitionmatrix, myinitialprobs, 30)

# A Hidden Markov Model of DNA sequence evolution

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


nucleotides = c("A","C","G","T")
states = c("AT-rich","GC-rich")
first_state = 'AT-rich'

probabilities = theemissionmatrix[first_state,]
seq = sample(nucleotides,1,prob=probabilities)

prev_state = first_state
next_state = sample(states,1,prob=thetransitionmatrix[prev_state,])
probabilities = theemissionmatrix[next_state,]
new_nucleotide = sample(nucleotides,1,prob=probabilities)
seq = c(seq,new_nucleotide)

prev_state = next_state
next_state = sample(states,1,prob=thetransitionmatrix[prev_state,])
probabilities = theemissionmatrix[next_state,]
new_nucleotide = sample(nucleotides,1,prob=probabilities)
seq = c(seq,new_nucleotide)


# Function to generate a DNA sequence, given a HMM and the length of the sequence to be generated.
generatehmmseq <- function(transitionmatrix, emissionmatrix, initial_state, seqlength)
{
  nucleotides = c("A", "C", "G", "T")   
  states = c("AT-rich", "GC-rich") 
  mysequence = character()             
  mystates = character()             

  firststate = initial_state
  probabilities = emissionmatrix[firststate,]
  firstnucleotide = sample(nucleotides, 1, prob=probabilities)
  mysequence[1] = firstnucleotide         
  mystates[1] = firststate              
  
  for (i in 2:seqlength)
  {
    prevstate    <- mystates[i-1]          
    stateprobs   <- transitionmatrix[prevstate,]
    state        <- sample(states, 1, prob=stateprobs)
    probabilities <- emissionmatrix[state,]
    nucleotide   <- sample(nucleotides, 1, prob=probabilities)
    mysequence[i] <- nucleotide             
    mystates[i]  <- state     
  }
  
  for (i in 1:length(mysequence))
  {
    nucleotide   <- mysequence[i]
    state        <- mystates[i]
    print(paste("Position", i, ", State", state, ", Nucleotide = ", nucleotide))
  }
}

theinitialprobs <- c(0.5, 0.5)
generatehmmseq(thetransitionmatrix, theemissionmatrix, theinitialprobs, 30)

# Inferring the states of a HMM that generated a DNA sequence

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