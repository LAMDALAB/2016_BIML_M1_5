# A Hidden Markov Model of DNA sequence evolution

nucleotides = c("A","C","G","T")

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

generatehmmseq(thetransitionmatrix, theemissionmatrix, 'AT-rich', 30)
