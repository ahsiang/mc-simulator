#	Mc model simulator
#
# 	Need to add ability to vary rho by cluster
#	Cannot deal with independent characters (i.e., clusters with size == 1)


writeData <- function(sim, file = NULL, numChars)
	{
	if (is.null(file)) file <- "data.csv"
	
	space <- data.frame(rep("\t",length(sim$tipNames)))
	toWrite <- cbind(sim$tipNames,space,sim$data)
	write(t(toWrite),file=file,sep=",",ncolumns=numChars+2)
	}


mcSim <- function(tree, numChars = NULL, allocGen = NULL, numClusters = NULL, charsPerCluster = NULL, rho = NULL, alphaDir = NULL)
	{
	# Set default values as necessary
	if (is.null(numChars)) numChars <- 10
	if (is.null(allocGen)) allocGen <- "prob"
	if (is.null(numClusters)) numClusters <- 3
	if (is.null(charsPerCluster)) charsPerCluster <- c(3,3,4)
	if (is.null(rho)) rho <- 0.05
	if (is.null(alphaDir)) alphaDir <- 0.05
		
	if (allocGen == "det")
		alloc <- genAllocDet(numChars, numClusters, charsPerCluster)
	if (allocGen == "prob")
		{
		alloc <- genAllocProb(numChars, alphaDir)
		numClusters <- max(alloc)
		charsPerCluster <- tabulate(alloc)
		}
		
	rates <- getAlphaBeta(rho)
	stateFreqs <- getStateFreqs(rho)
	root <- drawRootStates(numClusters, stateFreqs)
	
	numTaxa <- length(tree$tip.label)
	numNodes <- tree$Nnode

	latent <- evolveLatentStates(numTaxa, numNodes, root, numClusters, rates, tree)
	
	latentMatrix <- genLatentMatrix(latent$tipStates,numClusters,numTaxa)
	
	data <- emitData(latentMatrix, numClusters, charsPerCluster, numTaxa, numChars, alloc)
	
	reordered <- reorder(data$matrix, data$allocationVector, numTaxa)
	
	return(list("tree" = tree,
				"newick" = write.tree(tree),
				"data" = reordered$data,
				"allocationVector" = reordered$alloc, 
				"numClusters" = data$numClusters, 
				"numChars" = length(data$allocationVector), 
				"latentMatrix" = latentMatrix, 
				"nodeStates" = latent$nodeStates, 
				"tipStates" = latent$tipStates, 
				"nodeNames" = latent$nodeNames, 
				"tipNames" = latent$tipNames))
	}


genLatentMatrix <- function(tipStates, numClusters, numTaxa)
	{
	# Put latent states into matrix format
	latentMatrix <- matrix(nrow = numTaxa, ncol = numClusters)
	for (i in seq(length(tipStates)))
		{
		for (j in seq(numClusters))
			{
			latentMatrix[i,j] <- tipStates[[i]][j]
			}
		}
	return(latentMatrix)
	}


genAllocProb <- function(numChars, alphaDir)
	{
	# Creates probabilistic allocation vector based on Dirichlet concentration parameter alphaDir
	alloc <- rep(NA,numChars)
	alloc[1] <- 1
	highest <- 1
	for (i in 2:numChars)
		{
		newTableProb <- alphaDir / (alphaDir + i - 1)
		
		probs = c()
		for (j in seq(highest))
			{
			numAtTable <- 0
			for (k in seq(i-1))
				{
				if (alloc[k] == j)
					{
					numAtTable <- numAtTable + 1
					}
				}
			probs <- c(probs, numAtTable / (alphaDir + i - 1))
			}
			
		allProbs <- c(probs,newTableProb)
		newTable <- sample(seq(highest+1),size=1,prob=allProbs)		
			
		alloc[i] = newTable
		
		if (newTable == highest + 1)
			{
			highest <- highest + 1
			}
		}
	return (alloc)
	}


genAllocDet <- function(numChars, numClusters, charsPerCluster)
	{
	# Creates deterministic allocation vector based on number of characters, number of clusters, and number of characters per cluster
	alloc <- vector(mode="integer", length=numChars)
	cluster <- 1
	idx <- 1
	for (i in 1:numClusters)
		{
		end <- idx + charsPerCluster[i] - 1
		for (j in idx:end)
			{
			alloc[j] <- cluster
			}
		idx <- j + 1
		cluster <- cluster + 1
		}
	return(alloc)
	}
	
	
getAlphaBeta <- function(rho)
	{
	# Get alpha (fast rate) and beta (slow rate) from rho
	#beta <- rho / (2 * rho + 1)
	#alpha <- beta / rho
	
	alpha <- 1 / (2 + rho)
	beta <- 1 - (2 * alpha)
	
	return(list("alpha" = alpha,"beta" = beta))
	}
	
	
getStateFreqs <- function(rho)
	{
	# Get stationary state frequencies for 3-state model
	freq0 <- 1 / (2 + rho) # End state
	freq1 <- 1 - 2 * freq0 # Intermediate state
	freq2 <- freq0 # Opposite end state
	stateFreqs <- c(freq0,freq1,freq2)
	
	return(stateFreqs)
	}
	
	
drawRootStates <- function(numClusters, stateFreqs)
	{
	# Draws latent states at root
	root <- sample(seq(3), size=numClusters, prob=stateFreqs, replace=TRUE)
	
	return(root)
	}


transProbCalc <- function(t,a,b)
	{
	# Calculates transition probability matrix for a given branch length t, alpha, and beta
	pMatrix <- matrix(nrow=3,ncol=3)
	pMatrix[1,1] <- pMatrix[3,3] <- (exp(-b * t) * ((b * (1 + exp(-2 * a * t))) + ((2 * a) * (1 + exp(b * t))))) / 2
	pMatrix[1,2] <- pMatrix[3,2] <- b - (b * exp(-t))
	pMatrix[1,3] <- pMatrix[3,1] <- (exp(-b * t) * ((b * (-1 + exp(-2 * a * t))) + ((2 * a) * (-1 + exp(b * t))))) / 2
	pMatrix[2,1] <- pMatrix[2,3] <- a - (a * exp(-t))
	pMatrix[2,2] <- b + (2 * a * exp(-t))
	
	return(pMatrix) 
	}
	
	
evolveLatentStates <- function(numTaxa,numNodes,rootState,numClusters,rates,tree)
	{
	# Evolves latent states for all clusters up the tree from given root state
	intNodeStart <- numTaxa + 1
	intNodeEnd <- numTaxa + numNodes
	intNodes <- intNodeStart:intNodeEnd
	
	# Set up data structure to hold states at internal nodes
	nodeStates <- vector("list",numNodes)
	nodeStateNames <- array(dim=numNodes)
	for (i in seq(numNodes))
		{
		nodeStateNames[i] <- paste("node",toString(i + numTaxa),sep="")
		}
	
	# Set up data structure to hold states at tips
	tipStates <- vector("list",numTaxa)
	tipStateNames <- tree$tip.label
	
	edges <- as.data.frame(tree$edge)

	# Loop through interior nodes to get state changes
	for (i in intNodeStart:intNodeEnd)
		{
		currEdges <- subset(edges,edges["V1"] == i)
		
		currEdgeIndices <- as.integer(rownames(currEdges))
		for (j in 1:length(currEdgeIndices))
			{
			t <- tree$edge.length[currEdgeIndices[j]]
			p <- transProbCalc(t,rates$alpha,rates$beta)
			
			# Set starting state for current branch being processed
			if (i == intNodeStart)
				{
				startState <- rootState
				nodeStates[[i - numTaxa]] <- rootState
				} else {
				startState <- nodeStates[[i - numTaxa]]
				}
				
			# Get end state after evolution along current branch
			endState <- array(dim = numClusters)
			for (k in seq(numClusters))
				{
				probs <- p[startState[k],]
				newState <- sample(seq(3),size=1,prob=probs)
				endState[k] <- newState
				}				
			
			# Fill in data structures with appropriate state
			currIdx <- currEdges[j,"V2"]
			if (currIdx <= numTaxa)
				{
				idx <- currEdges[j,"V2"]
				tipStates[[idx]] <- endState
				} else {
				idx <- (currEdges[j,"V2"] - numTaxa)
				nodeStates[[idx]] <- endState
				}
			}
		}
		
	#names(nodeStates) <- nodeStateNames
	#names(tipStates) <- tipStateNames

	return(list("nodeStates" = nodeStates, 
				"tipStates" = tipStates, 
				"nodeNames" = nodeStateNames, 
				"tipNames" = tipStateNames))
	}	


emitData <- function(latentMatrix, numClusters, charsPerCluster, numTaxa, numChars, alloc)
	{
	# Probabilistically emit data at tips from latent states
	dataMatrix <- matrix(nrow = numTaxa, ncol = numChars)
	
	for (i in seq(numClusters))
		{
		numCharToSim <- charsPerCluster[i]
		latentColumn <- latentMatrix[,i]
		
		# Find offset for filling in data matrix
		if (i == 1)
			{
			offset <- 0
			} else
			{
			offset <- offset + charsPerCluster[i-1]
			}
			
		# Randomly draw a state to be the end state
		endState <- sample(c(0,1),size=numCharToSim,replace=TRUE)
		
		# Get opposite end state by flipping all states of chosen end state
		oppEndState <- array(dim=numCharToSim)
		for (j in seq(numCharToSim))
			{
			if (endState[j] == 0)
				{
				oppEndState[j] <- 1
				} else
				{
				oppEndState[j] <- 0
				}
			}
					
		# Loop through latent column and generate states
		for (j in seq(numTaxa))
			{
			if (latentColumn[j] == 1) # 0 end state
				{
				for (k in seq(numCharToSim))
					{
					dataMatrix[j,k+offset] <- endState[k]
					}
				} else if (latentColumn[j] == 3) # 2 end state
				{
				for (k in seq(numCharToSim))
					{
					dataMatrix[j,k+offset] <- oppEndState[k]
					}
				} else # Intermediate state
				{
				# Pick a random intermediate state
				NOTENDSTATE <- FALSE
				while (NOTENDSTATE == FALSE)		
					{
					tempState <- sample(c(0,1),size=numCharToSim,replace=TRUE)
					if (!all(tempState == endState) && !all(tempState == oppEndState))
						{
						NOTENDSTATE <- TRUE
						}
					}
				for (k in seq(numCharToSim))
					{
					dataMatrix[j,k+offset] <- tempState[k]
					}
				}
			}
		}
		
	# Remove invariable characters
	toRemove = c()
	for (j in seq(numChars))
		{
		if (length(unique(dataMatrix[,j])) == 1)
			{
			toRemove <- c(toRemove, j)
			}
		}
	finalMatrix <- dataMatrix[,-toRemove]
	allocationVector <- alloc[-toRemove]
		
	return(list("matrix" = finalMatrix, 
				"allocationVector" = allocationVector, 
				"numClusters" = max(allocationVector)))
	}
	
	
reorder <- function(data, allocationVector, numTaxa)
	{
	# Reorders characters so that all characters in the same cluster are together
	numChars <- length(allocationVector)
	numClusters <- max(allocationVector)
	charsPerCluster <- tabulate(allocationVector)
	
	rData <- matrix(data=NA, nrow=numTaxa, ncol=numChars)
	rAlloc <- rep(NA,numChars)
	
	trackIdx <- 1
	for (i in seq(numClusters))
		{		
		indices <- which(allocationVector %in% i)
		for (idx in indices)
			{
			rAlloc[trackIdx] <- i
			rData[,trackIdx] <- data[,idx]
			trackIdx <- trackIdx + 1
			}
		}
	return(list("data" = rData, 
				"alloc" = rAlloc))
	}

	