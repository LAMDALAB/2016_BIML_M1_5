# Gene expression data from Khan et al. (2001) - SRBCT Data set

## The small round blue cell tumors (SRBCTs) are 4 different childhood tumors named,
## so because of their similar appearance on routine histology, 
## which makes correct clinical diagnosis extremely challenging. 
## However, accurate diagnosis is essential because the treatment options, 
## responses to therapy and prognoses vary widely depending on the diagnosis. 
## They include Ewing's family of tumors (EWS), neuroblastoma (NB), 
## non-Hodgkin lymphoma (in our case Burkitt's lymphoma, BL) and rhabdomyosarcoma (RMS). 
## Our classification model was built to distinguish between these four tumors based on gene expression values.

## 두 가지 라이브러리를 통해서 데이터 로드 가능
## 1. bioconductor에 존재하는 made4 라이브러리 사용
## 2. R 레파지토리에 존재하는 plsgenomics 라이브러리 사용

## 1. bioconductor에 made4 라이브러리 사용

## In order to reduce the size of the MADE4 package, and produce small example datasets, the top 50
## genes from the ends of 3 axes following bga were selected. This produced a reduced datasets of 306 genes.

### 라이브러리 설치
source("https://bioconductor.org/biocLite.R")
biocLite("made4")

### 라이브러리 로드
library(made4)

### 데이터 불러오기
data(khan)
summary(khan)
khan$train
khan$test

## 2. plsgenomics 라이브러리 사용

## This data set contains 83 samples with 2308 genes: 29 cases of Ewing sarcoma (EWS), coded 1,
## 11 cases of Burkitt lymphoma (BL), coded 2, 18 cases of neuroblastoma (NB), coded 3, 25 cases
## of rhabdomyosarcoma (RMS), coded 4

### 라이브러리 설치
install.packages('plsgenomics')

### 라이브러리 로드
library(plsgenomics)

### 데이터 불러오기
data(SRBCT)
SRBCT$X
SRBCT$Y
dim(SRBCT$X)
sum(SRBCT$Y==1)
sum(SRBCT$Y==2)
sum(SRBCT$Y==3)
sum(SRBCT$Y==4)