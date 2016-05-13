#Work with one directory of test files

#Specify following variables as necessary

#Directory with test files
#workdir="C:\\LocalData\\em322\\Articles\\PCA norm\\Matlab\\n1000m0010q5p3mu01f01"
#Extension of resulting file
#ext = "_pcaL1.txt"
#Number of repetition of tests
#nRep = 1000
#meth = 1

#Pure dimension
q = 5
pc = 5

#Read file list
fNames=as.matrix(read.csv(paste(workdir,"\\fileList",sep=""),header=F))

#Number of files 
N = length(fNames)

#Create array for results
res = matrix(0.0, nrow=N, ncol=2)

#Loop of file processing
for (n in 1:N){
	#Load data from file
	test <- as.matrix(read.table(paste(workdir,"\\",fNames[n],sep=""), header=FALSE, sep = "\t", as.is=TRUE))

	#Estimate number of repetitions
	#Save initial time
	tim1 = proc.time()[3]
	switch(meth,
	{ #pcal1
		centre = colMedians(test)
		test1 = t(t(test)-centre)
		myres <- pcal1(test1, projDim=pc, center=FALSE, scores=TRUE, projPoints=TRUE, dispExp=TRUE, initialize="l2pca")
		myres$projPoints = t(t(myres$projPoints)+centre)},
	{ #l1pca
		centre = colMedians(test)
		test1 = t(t(test)-centre)
		myres <- l1pca(test1, projDim=pc, center=FALSE, projPoints=TRUE, initialize="l2pca", tolerance=0.0001, iterations=10)
		myres$projPoints = t(t(myres$projPoints)+centre)},
	{ #l2pca
		#For L2 we have no one function and use small script
		#Centralize
		centre = colMeans(test)
		test1 = t(t(test)-centre)
		#Compute PCs
		myres =svd(test1, nu = pc, nv = pc)
		#Compute reconstruction
		test1=myres$u %*% diag(myres$d[1:pc],pc,pc) %*% t(myres$v)
		myres$projPoints = t(t(test1)+centre)},
	{ #l1pca*
		centre = colMedians(test)
		test1 = t(t(test)-centre)
		myres = l1pcastar(test1, projDim=pc, center=FALSE, scores=TRUE, projPoints=TRUE, dispExp=TRUE)
		myres$projPoints = t(t(myres$projPoints)+centre)},
	{ #pcaPP
		#Calculate projects because it is nor calculated in PCAgrid
		myres = PCAgrid (test, k = pc, method = c ("mad", "sd", "qn"), maxiter = 10, splitcircle = 25, scores = TRUE, zero.tol = 1e-16, center = l1median, trace = 0, store.call = TRUE)
		#Compute reconstruction
		test1=myres$scores %*% t(myres$loadings)
		myres$projPoints = t(t(test1)+myres$center)}  #
	)
	#Save final time
	tim2 = proc.time()[3]
	tims = tim2 - tim1
	if (tims<0.0001)
		nRep = 10000
	else if (tims>=1)
		nRep = 1
	else
		nRep = round(1/tims)
	end

	print(cat("#of pc",pc,"Singular time",tims,"nRep",nRep,"   "))

nRep=1;

	if (nRep>1){
		#Real estimations
		#Save initial time
		tim1 = proc.time()[3]
		for (rep in 1:nRep){
			switch(meth,
			{ #pcal1
				centre = colMedians(test)
				test1 = t(t(test)-centre)
				myres <- pcal1(test1, projDim=pc, center=FALSE, scores=TRUE, projPoints=TRUE, dispExp=TRUE, initialize="l2pca")
				myres$projPoints = t(t(myres$projPoints)+centre)},
			{ #l1pca
				centre = colMedians(test)
				test1 = t(t(test)-centre)
				myres <- l1pca(test1, projDim=pc, center=FALSE, projPoints=TRUE, initialize="l2pca", tolerance=0.0001, iterations=10)
				myres$projPoints = t(t(myres$projPoints)+centre)},
			{ #l2pca
				#For L2 we have no one function and use small script
				#Centralize
				centre = colMeans(test)
				test1 = t(t(test)-centre)
				#Compute PCs
				myres =svd(test1, nu = pc, nv = pc)
				#Compute reconstruction
				test1=myres$u %*% diag(myres$d[1:pc],pc,pc) %*% t(myres$v)
				myres$projPoints = t(t(test1)+centre)},
			{ #l1pca*
				centre = colMedians(test)
				test1 = t(t(test)-centre)
				myres = l1pcastar(test1, projDim=pc, center=FALSE, scores=TRUE, projPoints=TRUE, dispExp=TRUE)
				myres$projPoints = t(t(myres$projPoints)+centre)},
			{ #pcaPP
				#Calculate projects because it is nor calculated in PCAgrid
				myres = PCAgrid (test, k = pc, method = c ("mad", "sd", "qn"), maxiter = 10, splitcircle = 25, scores = TRUE, zero.tol = 1e-16, center = l1median, trace = 0, store.call = TRUE)
				#Compute reconstruction
				test1=myres$scores %*% t(myres$loadings)
				myres$projPoints = t(t(test1)+myres$center)} ##
			)
		}
		#Save final time
		tim2 = proc.time()[3]
	}
	#Save used time to array
	res[n,1]=(tim2-tim1)/nRep
	#calculate error of reconstruction
	col = ncol(myres$projPoints)
	row = nrow(myres$projPoints)
	res[n,2]=sum(abs(myres$projPoints[,(q+1):col]))/row
}

print(paste0("Current working dir: ", workdir))
#Save final file
write.table(res, file=paste(workdir,ext,sep=""), row.names=FALSE, col.names=FALSE)
