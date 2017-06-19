checkCombosOccurEqually<- function(df,colNames,dropZeros=FALSE) {
	#in data.frame df, check whether the factors in the list colNames reflect full factorial design (all combinations of levels occur equally often)
	#
	#dropZeros is useful if one of the factors nested in the others. E.g. testing different speeds for each level of    
	# something else, then a lot of the combos will occur 0 times because that speed not exist for that level.
	#but it's dangerous to dropZeros because it won't pick up on 0's that occur for the wrong reason- not fully crossed
	#
	#Returns:
	# true/false, and prints informational message
	#
	listOfCols <- as.list( df[colNames] )
	t<- table(listOfCols)
	
	if (dropZeros) {  
		t<- t[t!=0]   
	}           
	colNamesStr <- paste(colNames,collapse=",")
	if ( length(unique(t)) == 1 ) { #if fully crossed, all entries in table should be identical (all combinations occur equally often)
		  print(paste(colNamesStr,"fully crossed- each combination occurred",unique(t)[1],'times'))
		  ans <- TRUE
	  } else {
		  print(paste(colNamesStr,"NOT fully crossed,",length(unique(t)),'distinct repetition numbers.'  ))
		  ans <- FALSE
	  }	
	return(ans)
}


# # checkAllCombosOccurEquallyOften(dat, c("numObjects","numTargets","ringToQuery") )

# checkAllCombosOccurEquallyOften(dat, c("numObjects","numTargets") )

# checkAllCombosOccurEquallyOften(dat, c("numTargets","ringToQuery") )

# checkAllCombosOccurEquallyOften(dat, c("numObjects","numTargets","speed") ) #not counterbalanced (speed nested)

# #If instead of using raw speed, I rank the speed within each numObjects*numTargets, then from that perspective everything should
# #be perfectly counterbalanced, because each numObjects*numTargets combination has the same number of speeds tested
# #But the rank for a speed depends on what numObjects-numTargets condition it's in. Should be easy with ddply
# ordinalSpeedAssign <- function(df) {
	# #df$speedRank <- rank(df$speed)  #Rank won't work, always wants to break all ties. Whereas I want to preserve ties.
	# df$speedRank <- match(xx$speed,unique(xx$speed)) 
	# df
# }

# dat<- ddply(dat,.(numObjects,numTargets),ordinalSpeedAssign)
# #now check whether counterbalanced given this context. HERE, I LEFT OFF HERE
# checkAllCombosOccurEquallyOften(dat, c("numObjects","numTargets","speed") )


# dat<- ddply(dat,.(numObjects,numTargets),function(df) xx<<-df)



# #numTargets and numObjects and speed
# #create a quick unique id for each numTargets/numObjects combo, then make sure within each have same #trials at each speed
# numTargetsObjsUid<- dat$numTargets * length(unique(dat$numObjects)) + dat$numObjects #Give each combination a unique number.
# speedPerTargetsAndObjs <- table(numTargetsObjsUid,dat$speed)
# #Contains a lot of zeros corresponding to speeds not tested for that particular targets*objects combo. But for those that are tested,
# #number of trials should be equal.
# if ( length(unique(speedPerTargetsAndObjs[ speedPerTargetsAndObjs !=0 ]))>1 ) {
	# paste("Speed not counterbalanced with numTargets*numObjects combos for",dataDir)
	
# }

# #Check whether ringToQuery is counterbalanced with all combinations of numTargets,numObjects, and speed
# table(dat$numTargets,dat$numObjects,dat$speed)
# dat$firstThreeFactorsAssumedOk <- do.call(paste, c(dat[c("numTargets","numObjects", "speed")], sep = " "))
# tb<-table(dat$firstThreeFactorsAssumedOk,dat$ringToQuery)
# if (length(unique(tb)) == 1) {  #if counterbalanced, all entries in table should be identical (all combinations occur equally often)
	# paste("ringToQuery counterbalanced with all unique combos of numTargets*numObjects*speed data in",dataDir)
# } else {
	# paste("ringToQuery NOT counterbalanced with all unique combos of numTargets*numObjects*speed data in",dataDir)
# }

# fullCounterBalanceKey <- do.call(paste, c(dat[c("numObjects", "speed")], sep = " "))
# #could also check which rings have targets, by going through whichIsTarget0  !=-999
# dat$targetInRing0 <- (dat$whichIsTarget0  != -999) *1
# dat$targetInRing1 <- (dat$whichIsTarget1  != -999) *1
# dat$targetInRing2 <- (dat$whichIsTarget2  != -999) *1
# timesAtargetIsInEachRing<- c(sum(dat$targetInRing0),sum(dat$targetInRing1),sum(dat$targetInRing2))
# if (length(unique(timesAtargetIsInEachRing)) != 1) {#if counterbalanced, all entries in table should be identical (all combos occur equally often)
	# paste("timesAtargetIsInEachRing not the same for all rings for ",dataDir)
# } else {
	# paste("timesAtargetIsInEachRing is the same for all rings for ",dataDir)
# }
# #check counterbalancing of 
# #Often want to check counterbalancing of some factor(s), conditional on some other factor.
# #For example, ringToQuery conditional on dat$firstThreeFactorsAssumedOk. Also targetInRing0 conditional on numTargets

# #end check of counterbalancing
