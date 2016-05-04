#Test all

#Resulting files are written to the workdir1 directory and 
#have names name_of_tested_dir + ext + '.txt' and contain two column. 
#The first column contains time in seconds 
#and the second column contains error.

#load used libraries
load("pcaL1")
load("pcaPP")

#Define required values
meth = 1 #1 is pcal1, 2 is l1pca, 3 is l2pca, 4 is l1pca*, 5 is pcaPP

switch(meth,
{ nRep = 1000	#pcal1
  ext = "_pcaL1.txt"},
{ nRep = 5	#l1pca
  ext = "_l1pca.txt"},
{ nRep = 5000	#l2pca
  ext = "_l2pca.txt"},
{ nRep = 2	#l1pca*
  ext = "_l1pcaS.txt"},
{ nRep = 25	#pcaPP
  ext = "_pcaPP.txt"})

#Directory with directories with test files
workdir1="C:\\LocalData\\em322\\Articles\\PCA norm\\Matlab\\"

#Work with each directory

workdir = paste(workdir1,"n1000m0010q5p0mu0f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p1mu01f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p1mu05f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p1mu10f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p1mu25f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p2mu01f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p2mu05f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p2mu10f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p2mu25f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p3mu01f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p3mu05f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p3mu10f01",sep="")
source("testDirectory.R")

workdir = paste(workdir1,"n1000m0010q5p3mu25f01",sep="")
source("testDirectory.R")
