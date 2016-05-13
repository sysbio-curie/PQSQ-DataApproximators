#Test all

#Resulting file are written to the current directory and 
#has name Test_name + ext and contain four column. 
#The first column contains number of PCS, the second column 
#contains time in seconds, the third column contains L1 error
#and the fourth column contains L2 error.

#load used libraries
#load("pcaL1")
#load("pcaPP")

#Define required values
Test_name = "wdbc"
#meth = 1 #1 is pcal1, 2 is l1pca, 3 is l2pca, 4 is l1pca*, 5 is pcaPP
nnorm = ncol(test)*nrow(test)

#Create array for results
res = matrix(0.0, nrow=10, ncol=4)

for (meth in 1:5){
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

	#Number of PC loop
	for (pc in 1:10){
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
			myres$projPoints = t(t(test1)+myres$center)}
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
					myres$projPoints = t(t(test1)+myres$center)}
				)
			}
			#Save final time
			tim2 = proc.time()[3]
		}
		#Save results
		res[pc,1] = pc
		#Save used time to array
		res[pc,2] = (tim2-tim1)/nRep
		#calculate L1 error of reconstruction
		res[pc,3]=sum(abs(myres$projPoints-test))/nnorm
		#calculate L2 error of reconstruction
		res[pc,4]=var(as.vector(myres$projPoints-test))
	}

	#Save final file
	write.table(res, file=paste(Test_name,ext,sep=""), row.names=FALSE, col.names=FALSE)
}
