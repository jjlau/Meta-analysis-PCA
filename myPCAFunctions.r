# Functions file: store PCA functions
# Last Updated: 2011.08.30
# By Jason Lau
######################################

setupPCA <- function()
{
	myfile <- txt()
	mydata <- perMale(myfile)
	
	mydata <- naToMedian(mydata)
	mydata <- factorize(mydata) 
	mydata.list <- normZ(mydata) 
	print("List contains: $data.unnorm, $data.norm, $data.sd")
	print("To do PCA, type: _new list_ <- runPCA(__setup list__)")
	return (mydata.list)
}

runPCA <- function(datalist)
{
	data.norm <- datalist[[2]] #factorized data
	
	#data is a matrix with columns and rows of interest
	#returns coefficients of comp1 and comp2
	if (!is.matrix(data.norm)) data.norm <- as.matrix(data.norm)
	
	# get the coefficients of the principal components - these are the eigenvectors
	pca.dist <- prcomp(na.omit(data.norm))
	pca.rot <- pca.dist$rotation
	
	# get components 1 and 2 
	coef.12 <- pca.rot[,1:2]

	# project data onto the 2 dimensional surface defined by the first two eigenvectors
	# this is the same results as with classical MDS
	pca.coord <- data.norm %*% coef.12 
	
	dev.new()
	screeplot(pca.dist)
	
	pca.dist <- pca.dist$rotation
	dev.new()
	barplot(abs(pca.dist[,1]), main = "Distribution of Comp1")
	dev.new()
	barplot(abs(pca.dist[,2]), main = "Distribution of Comp2")
	
	dev.new()
	plot(pca.coord)
	#text(pca.coord, rownames(data.norm), cex = .75)
	
	print("List contains: $data.unnorm, $data.norm, $data.sd, $pca.dist, $coef.12, #pca.coord")
	print("Follow the addRefRange instructions to add default or custom reference ranges to the graph")
	return (list(data.unnorm = datalist[[1]],data.norm = datalist[[2]],data.sd = datalist[[3]], pca.dist = pca.dist, coef.12 = coef.12, pca.coord = pca.coord))	
}
addRefRange <- function (myList, userInput = NULL)
{
	
	data.unnorm <- myList[[1]]
	data.norm <- myList[[2]]
	data.sd <- myList[[3]]
	coef.12 <- myList[[5]]
	pca.coord <- myList[[6]]
	ncols <- nrow(coef.12)
	
	choice <- 0
	device <- as.numeric(readline("Which graphics device would you like to add to? (0=new):"))
	
	
	m <- genBinaryMatrix(ncols) #generate binary matrix
		
		maxX <- mat.or.vec(1,ncols)
		minX <- mat.or.vec(1,ncols)
		maxM <- m
		minM <- abs(m-1)
	
	refM <- mat.or.vec(1,ncols)
	
	
	if (is.null(userInput)){
	choice <- readline("Options- 1=add your own ranges, 2=show max and min: ")
	}
	else{
		refM <- userInput
		print("Taking in user input row vector")
	}
	if (any(choice == 1, is.numeric(userInput))){
		
		if (choice ==1){
			print(summary(data.unnorm))
			for (columnNum in 1:ncols){
				refM[1,columnNum]<- as.numeric(readline(paste("Input reference range for ",colnames(data.unnorm)[columnNum], ":")))
				
			}
		}		
		for (column in 1:ncols){
			maxX <- refM[,column]/data.sd[1,column]
			maxM[,column] <- maxM[,column]*maxX
		}
		
		refrange <- maxM%*% coef.12
	}
	
	#This is a sanity check. This takes the max and min of each variable
	#projects the "volume" of 2^n onto the plot
	if (choice == 2){
		refM <- mat.or.vec(2,ncols)
		for (column in 1:ncols){
			maxX <- max(data.norm[,column])
			maxM[,column] <- maxM[,column]*maxX
			refM[1,column] <- maxX
			
			minX <- min(data.norm[,column])
			minM[,column] <- minM[,column]*minX
			refM[2,column] <- minX
		}
		mrange <- maxM + minM
		refrange <- mrange%*% coef.12
		
		refM <- as.data.frame(refM)
		colnames(refM) <- colnames(data.unnorm)
		rownames(refM) <- c("max", "min")
		
	}
	
	refrange <- cbind(refrange, matrix(TRUE, nrow(refrange),byrow= TRUE))
	references <- gramScan(refrange,0)
	refPlot <- rbind(references, references[1,])
	
	if (choice <=1) {
		refPlot[,1] <- refPlot[,1] + abs(min(refPlot[,1]))
		refPlot[,2] <- refPlot[,2] + abs(min(refPlot[,2]))
		
		#sets the mean of the custom reference range to the mean of the pca coordinates
		if (mean(pca.coord[,1]) < 0) refPlot[,1] <- refPlot[,1] - (abs(mean(refPlot[,1])) + abs(mean(pca.coord[,1])))
		else refPlot[,1] <- refPlot[,1] - abs(mean(refPlot[,1])) + abs(mean(pca.coord[,1]))
		
		if (mean(pca.coord[,2]) < 0) refPlot[,2] <- refPlot[,2] - (abs(mean(refPlot[,2])) + abs(mean(pca.coord[,2])))
		else refPlot[,2] <- refPlot[,2] - abs(mean(refPlot[,2])) + abs(mean(pca.coord[,2]))
		
		references <- refPlot
	}
	
	minxlim <- min(c(refPlot[,1], pca.coord[,1]))
	maxxlim <- max(c(refPlot[,1], pca.coord[,1]))
	
	minylim <- min(c(refPlot[,2], pca.coord[,2]))
	maxylim <- max(c(refPlot[,2], pca.coord[,2]))
	
	require(stats)
	
	if (device > 1 ){
		dev.set(which = device)
		color <- readline("What color line would you like? (red, green, black...):")
		points(refPlot[,1],refPlot[,2], type ="l", col =color)
		
	}
	
	if (device == 0){
		dev.new()
		plot(pca.coord, col = "black", xlim = c(minxlim, maxxlim), ylim = c(minylim, maxylim) )
		points(refPlot[,1],refPlot[,2], type ="l", col ="red")
	}
	else if (device == 1){
		
		plot(pca.coord, col = "black", xlim = c(minxlim, maxxlim), ylim = c(minylim, maxylim) )
		points(refPlot[,1],refPlot[,2], type ="l", col ="red")
	}
	refrange.coord  <- references[,1:2]
	
	return (list(refrange.coord = refrange.coord, refrange.input = refM))

	
}

#pca.coord <- myList$pca.coord
#refrange <- refrange3[[1]]
addLabels <- function(myList, refrange = NULL)
{
	if (length(myList) == 6) pca.coord <- myList[[6]]
	else pca.coord <- myList
	
	nrows <- nrow(pca.coord)
	#orig.ref <- refrange
	
	device <- as.numeric(readline("Which graphics device would you like to add to?:"))
	dev.set(which = device)
	
	if (is.null(refrange)){
		color <- readline("What color points would you like? (red, green, black...):")
		points(pca.coord[,1],pca.coord[,2], col =color)
		textsize <- as.numeric(readline("Labeling all text. What size font? (0.5, 0.75, 1...) :"))
		x <- as.numeric(readline("How far off the point in the x (-1, 0, 1...) :"))
		y <- as.numeric(readline("How far off the point in the y (-1, 0, 1...) :"))
		text(pca.coord, rownames(pca.coord), cex = textsize, adj = c(x,y))
	}
	else{
		if (length(refrange) == 2) refrange <- refrange[[1]]
		if (refrange[1,1] == refrange[nrow(refrange),1]) refrange <- refrange[1:(nrow(refrange)-1),]
		tempRef <- refrange
		pca.coord <- cbind(pca.coord, matrix(TRUE, nrow(pca.coord),byrow= TRUE))
		for (i in 1:nrows){
			tempRef <- cbind(tempRef, matrix(TRUE, nrow(tempRef),byrow= TRUE))
			tempRef <- rbind(tempRef,pca.coord[i,])
			
			#Run gramScan on point with refrange
			checkRef <- gramScan(tempRef,2)
			
			#If point is inside new list, then point is outside, set false
			for (j in 1:nrow(checkRef)){
				if (pca.coord[i,1]==checkRef[j,1]){
					pca.coord[i,3] <- FALSE
				}
			
			}
			
			tempRef <- refrange
		}
		
		choice <- readline("Label outside or inside points? (0 = outside, 1 = inside) :")
		if (choice == 0){
			pca.coord <- pca.coord[(pca.coord[,3]==0),]
		}
		else if (choice == 1){
			pca.coord <- pca.coord[(pca.coord[,3]==1),]
		}
		color <- readline("What color points would you like? (red, green, black...):")
		points(pca.coord[,1],pca.coord[,2], col =color)
		textsize <- as.numeric(readline("Labeling all text. What size font? (0.5, 0.75, 1...) :"))
		x <- as.numeric(readline("How far off the point in the x (-1, 0, 1...) :"))
		y <- as.numeric(readline("How far off the point in the y (-1, 0, 1...) :"))
		text(pca.coord, rownames(pca.coord), cex = textsize, adj = c(x,y))
		
		return(pca.coord)
	}
}

#To createPDF, requires the list from runPCA and a list of refranges
#refrange.list needs to be in the format of  createPDF(myList,list(__refrange1__$refrange.coord,__refrange2__$refrange.coord))
createPDF <- function(pcaList, refrange.list)
{
	if (length(pcaList) == 6) pca.coord <- pcaList[[6]]
	else pca.coord <- pcaList
	refPlot <- refrange.list[[1]]
	nranges <- length(refrange.list)[1]
	minxlim <- min(c(refPlot[,1], pca.coord[,1]))
	maxxlim <- max(c(refPlot[,1], pca.coord[,1]))
	
	minylim <- min(c(refPlot[,2], pca.coord[,2]))
	maxylim <- max(c(refPlot[,2], pca.coord[,2]))
	refPlot <- rbind(refPlot, refPlot[1,])
	
	name <- readline("File of PDF:")
	pdf(file = name)
	
	
	plot(pca.coord, col = "black", xlim = c(minxlim, maxxlim), ylim = c(minylim, maxylim) )
	
	
	for (i in 1:nranges){
		refPlot <- refrange.list[[i]]
		refPlot <- rbind(refPlot, refPlot[1,])
		color <- readline(paste("What color line would you like for refrange ", i,"? (red, green, black...):"))
		points(refPlot[,1],refPlot[,2], type ="l", col =color)
		
	}
	dev.off()
}

txt <- function()
{
	myfile <- readline("Text File Name: ")
	mydata <-  read.table(myfile,  header=TRUE )
	
	rownames(mydata)<- paste(mydata[,1], c(1:nrow(mydata)), sep="_")
	mydata[,1] <- NULL
	return (mydata)
}

#Takes "n" and "male" column and replaces with "percent male"
perMale <-function(mydata)
{
	
	#replace male with percent male
	tempdata <- mydata$male/mydata$n
	med <- median(tempdata, na.rm=TRUE)
	tempdata <- as.matrix(tempdata)
	
	
	for (i in 1:nrow(tempdata)){
		if (is.na(tempdata[i])) tempdata[i] <- med
	}
	
	mydata$male <- tempdata
	mydata$n <-NULL #remove n column with total population
	
	return(mydata)

}

#Replaces all NA with medians of the each column
#percent male needs to be calculated before hand - why?
naToMedian <- function(mydata)
{
	ncol <- dim(mydata)[2]
	
	for (i in 1:ncol) {
		a.column <- mydata[, i]
		if (is.numeric(a.column)) {
			the.NAs <- is.na(a.column)
			if (any(the.NAs)) {
				med <- median(a.column, na.rm=TRUE)
				a.column[the.NAs] <- rep(med, sum(the.NAs))
				mydata[,i] <- a.column
			}
		}
	} 
	
	return(mydata)
	
}

#Finds Factor Columns, creates dummy columns, removes original cols
factorize <- function(mydata)
{
	ncols <- dim(mydata)[2]
	nrows <- dim(mydata)[1]
	
	for (k in 1:ncols) {
		if (is.factor(mydata[,k])){
			dumlist <- createDum(mydata[,k], colnames(mydata)[k], nrows)
			mydata<- cbind(mydata, dumlist)
			mydata[,k] <- NULL
		}
			
	}
	return (mydata)
}

#Creates the dummy lists 
createDum <- function(catcol, variableName, nrows)
{
	tempdata <- as.matrix(catcol)
	nlevel <- nlevels(catcol)
	
	dumlist <- mat.or.vec(nrows,nlevel) #blank matrix representing all dummies
	
	for (j in 1:nlevel) { #move through each column
		tempdata <- as.matrix(catcol)	
		myvar <- levels(catcol)[j]
	
		for (i in 1:nrows){ #move down each row
			if (catcol[i] == myvar) tempdata[i] <- 1
			else tempdata[i] <- 0
		}
		dumlist[,j] <- as.numeric(tempdata)
	}
	#dummy columns get a special "..." in their header name
	colnames(dumlist) <- paste(rep(variableName,nlevel), levels(catcol), sep="...")
	return (dumlist)
}

#Normalizes data by subtracting mean and divding by SD
normZ <- function(data.unnorm)
{
	data.norm <- data.unnorm
	nrows <-dim(data.unnorm)[1]
	ncols <-dim(data.unnorm)[2]
	#mymean<-matrix(rep(mean(as.data.frame(data), na.rm=TRUE),rows),ncol=column, byrow=TRUE)
	#dsd<-matrix(rep(sd(data, na.rm=TRUE),rows),ncol=column, byrow=TRUE)
	data.sd <- mat.or.vec(1,ncols)
	
	for (k in 1:ncols) {
		header <- colnames(data.unnorm)[k]
		#refrange[,i] <- 1 #change default for factor columns
		
		if(!grepl("\\.\\.\\.", header)){
			colmean<-mean(as.data.frame(data.norm[,k]), na.rm=TRUE)
			colsd<-sd(data.norm[,k], na.rm=TRUE)
			data.sd[,k] <- colsd
			data.norm[,k] <- (data.norm[,k]-colmean)/colsd
		}	
		else data.sd[,k] <- 1
	}
	
	return(list(data.unnorm = data.unnorm, data.norm = data.norm,data.sd =data.sd))
}

genBinaryMatrix <- function (lvl)
{
	m <- rbind(1,0)
	
	for (i in 1:(lvl-1)){
		m <- rbind(m,m)
		nrows <- nrow(m)
		n<-mat.or.vec(nrows,1)
		
		for (half in 1:(nrows/2)){
			n[half] <- 1
		}
		
		m <- cbind(n,m)
	}
	return (m)
}

gramScan <- function (refrange, repetitions)
{	
	continueLoop <- TRUE
	#outsidePoint <- FALSE
	if (repetitions == 0){
		while (continueLoop == TRUE){
			continueLoop <- FALSE
			
		
			
			colnames(refrange)[3] <- "border"
			refrange[,"border"] <- FALSE
	
			#find coordinates of smallest y value
			ysort <- sort(refrange[,2], index.return = TRUE)
			smallestID = ysort$ix[1]
			ysort <- NULL
	
			shiftM <- refrange[,1:2]
			#shift to smallestID = origin
			#shift X
			if (refrange[smallestID,1] < 0)shiftM[,1] <- shiftM[,1] + abs(refrange[smallestID,1])
			else shiftM[,1] <- shiftM[,1] - abs(refrange[smallestID,1])
			
			#shift Y
			if (refrange[smallestID,2] <0) shiftM[,2] <- shiftM[,2] + abs(refrange[smallestID,2])
			else shiftM[,2] <- shiftM[,2] - abs(refrange[smallestID,2])
				
				
			#create polar coordinate matrix
			polarM <- shiftM
			xCoords <- shiftM[,1]
			yCoords <- shiftM[,2]
			polarM[,1] <- ((shiftM[,1])^2 + (shiftM[,2])^2)^0.5
			polarM[,2] <- acos(xCoords/polarM[,1]) 
			#polarM[,2] <- atan(xCoords/yCoords)
			polarM[smallestID,2] <- 0
			colnames(polarM) <- c("radius", "theta")
			xCoords <- NULL
			yCoords <- NULL
	
			#create index with increasing theta and decreasing radius
			nrows <- nrow(shiftM)
			negPol<- polarM[,"radius"]*rep((-1),nrows)
			thetaIndex <- order(polarM[,2],negPol)
			inverseIndex <- sort(thetaIndex, index.return = TRUE) #this can be used to revert matrix
	
			polarM <- polarM[thetaIndex,] 
			shiftM <- shiftM[thetaIndex,]
			refrange <- refrange[thetaIndex,]
		
			#set the first two rows in theta as true points
			refrange[1,"border"] <- TRUE
			refrange[2, "border"] <- TRUE
	
			thetaN <- 3
			convexChange <- FALSE
			for (thetaN in 3:(nrows)){
		
				pNew <- refrange[thetaN,1:2]
		
				if (!(polarM[thetaN-1,2] == polarM[(thetaN),2])){
			
					pVert <- refrange[thetaN-1,1:2] #set new Vertex
					refrange[thetaN,"border"] <- TRUE #pNew gets set as TRUE
			
			
					if (!convexChange) {
						pOld <- refrange[thetaN-2,1:2] #set new Old only if no convex was found
					}		
		
					if (det(rbind(pNew-pVert,pOld-pVert)) < (1*10^(-7))){
						refrange[thetaN-1,"border"] <- FALSE #turn off ones that cause convex
						convexChange <- TRUE
						continueLoop <- TRUE
						#print("found a convex:")
						#print(refrange[theaN-1,])
					}
					else convexChange <- FALSE
				}
			}
			refrange <- refrange[(refrange[,3]==1),]
	
		}
	}
	else{
		for(i in 1: repetitions){
			
			colnames(refrange)[3] <- "border"
			refrange[,"border"] <- FALSE
	
			#find coordinates of smallest y value
			ysort <- sort(refrange[,2], index.return = TRUE)
			smallestID = ysort$ix[1]
			ysort <- NULL
	
			shiftM <- refrange[,1:2]
			#shift to smallestID = origin
			#shift X
			if (refrange[smallestID,1] < 0)shiftM[,1] <- shiftM[,1] + abs(refrange[smallestID,1])
			else shiftM[,1] <- shiftM[,1] - abs(refrange[smallestID,1])
			
			#shift Y
			if (refrange[smallestID,2] <0) shiftM[,2] <- shiftM[,2] + abs(refrange[smallestID,2])
			else shiftM[,2] <- shiftM[,2] - abs(refrange[smallestID,2])
				
	
			#create polar coordinate matrix
			polarM <- shiftM
			xCoords <- shiftM[,1]
			yCoords <- shiftM[,2]
			polarM[,1] <- ((shiftM[,1])^2 + (shiftM[,2])^2)^0.5
			polarM[,2] <- acos(xCoords/polarM[,1]) 
			#polarM[,2] <- atan(xCoords/yCoords)
			polarM[smallestID,2] <- 0
			colnames(polarM) <- c("radius", "theta")
			xCoords <- NULL
			yCoords <- NULL
	
			#create index with increasing theta and decreasing radius
			nrows <- nrow(shiftM)
			negPol<- polarM[,"radius"]*rep((-1),nrows)
			thetaIndex <- order(polarM[,2],negPol)
			inverseIndex <- sort(thetaIndex, index.return = TRUE) #this can be used to revert matrix
	
			polarM <- polarM[thetaIndex,] 
			shiftM <- shiftM[thetaIndex,]
			refrange <- refrange[thetaIndex,]
		
			#set the first two rows in theta as true points
			refrange[1,"border"] <- TRUE
			refrange[2, "border"] <- TRUE
	
			thetaN <- 3
			convexChange <- FALSE
			for (thetaN in 3:(nrows)){
		
				pNew <- refrange[thetaN,1:2]
		
				if (!(polarM[thetaN-1,2] == polarM[(thetaN),2])){
			
					pVert <- refrange[thetaN-1,1:2] #set new Vertex
					refrange[thetaN,"border"] <- TRUE #pNew gets set as TRUE
			
			
					if (!convexChange) {
						pOld <- refrange[thetaN-2,1:2] #set new Old only if no convex was found
					}		
		
					if (det(rbind(pNew-pVert,pOld-pVert)) < (1*10^(-7))){
						refrange[thetaN-1,"border"] <- FALSE #turn off ones that cause convex
						convexChange <- TRUE
						continueLoop <- TRUE
						#print("found a convex:")
						#print(refrange[theaN-1,])
					}
					else convexChange <- FALSE
				}
			}
			refrange <- refrange[(refrange[,3]==1),]
	
		}
		#return (outsidePoint)
	}
	#dev.new()
	refrange <- refrange[(refrange[,3]==1),]
	#plot(refrange)
	
	return(refrange)
	
}



#continue <- as.numeric(readline("Would you like to add another plot? (1 = yes, 2 = no):"))
	#while(continue ==1){
		
	#	if (continue == 1){
	#		addRefrange <- as.numeric(readline("Additional Refrange:"))
	#		addRefPlot <- addRefrange[[1]]
	#		addRefPlot <- rbind(addRefPlot, addRefPlot[1,])
	#		color <- readline("What color line would you like? (red, green, black...):")
	#		points(addRefPlot[,1],addRefPlot[,2], type ="l", col =color)
	#		continue <- as.numeric(readline("Would you like to add another plot? (1 = yes, 2 = no):"))
	#	}
	#}

#points(refrange2[,1],refrange2[,2], type ="p", col ="red")
#points(refPlot[,1],refPlot[,2], type ="p", col ="red")
	#dev.new()
	#plot(pcaCoord, col = "green", xlim = c(minxlim, maxxlim), ylim = c(minylim, maxylim) )
	#text(pcaCoord, rownames(normdata), cex = .5)
	#points(mean(refrange[,1]), mean(refrange[,2]), col = "green")
#refrange[1:20,]
#plot(refrange)
#temp <- refrange[(refrange[,3]==1),]
#points(temp[,1],temp[,2], col = "red")
# det(rbind(p0-p1,p2-p1))
#c$ix[1]
#c <- sort(refrange[,1], decreasing = TRUE, index.return = TRUE)






dmean <- function(data)
{
	rows <-dim(data)[1]
	column <-dim(data)[2]
	mymean<-matrix(rep(mean(as.data.frame(data), na.rm=TRUE),rows),ncol=column, byrow=TRUE)
	#dsd<-matrix(rep(sd(data, na.rm=TRUE),rows),ncol=column, byrow=TRUE)
	return(data-mymean)
}

myMDS <- function(data,names)
{

	d <- dist(data)
	fit <- cmdscale(d, eig=TRUE, k=2)

	x <-fit$points[,1]
	y <-fit$points[,2]
	plot(x,y)
	text(x,y, labels = names, cex=.7)

}


 #thetaSort <- order(polarM[,2],negPol)
 #polarM[thetaSort,]




#test <- mat.or.vec(1,2)
#for (i in 1:nrow(refrange)){
#	if (refrange[i,3]==1) {
#	add <- c(refrange[i,1],refrange[i,2])
#	test <- rbind(test,add)
#	}
#
#}
