#expecting current working directory to be top level of this git-indexed project, and this file to be in top level - dataPreprocess/

#Set to your top-level directory
setwd("/Users/alexh/Documents/attention_tempresltn/EEGRSVP/EEG_RSVP")

destinatnDir<-"dataAnonymized/" #where the anonymized data will be exported to
anonymiseData <- TRUE

thisExpFolder = "dataRaw/"
#In the dataRAw folder, each participant is expected to have their own folder with all their files.
#  In many cases there a participant will have only one .txt file, but could have multiple in case of multiple sessions.

foldersThisExp <- list.dirs(path=thisExpFolder,recursive=FALSE) #each folder should be a subject
print("Loading data from folders:"); print(foldersThisExp)
for (i in 1:length(foldersThisExp)) {
  thisSubjectDir <- foldersThisExp[i]
  files <- dir(path=thisSubjectDir,pattern='.txt')  #find all data files in this directory
  for (j in 1:length(files)) { #read in the sessions of this subject
      file = files[j]
      fname = paste0(thisSubjectDir,"/",file)
	    rawDataLoad=tryCatch( 
      	    		read.table(fname,sep='\t',header=TRUE), 
      	    		error=function(e) { 
      	    	   			stop( paste0("ERROR reading the file ",fname," :",e) )
          		 } )
      rawDataLoad$file <- file
      numTrials<- length(rawDataLoad$trialnum)
      msg=''
      rawDataThis<- rawDataLoad
      cat(paste0("Loaded file ",file,msg))
      #omit first trial if total trials are odd, last probably a repeat. And first trial people often discombobulated      
      msg=""
      removeFirstTrialIfOdd = TRUE
      if (numTrials %% 2 ==1) {
      	msg=paste0(" Odd number of trials (",numTrials,"); was session incomplete, or extra trial at end?")  
        if (removeFirstTrialIfOdd) {
      	  rawDataThis <- subset(rawDataThis, !trialnum %in% c(0))
      	  cat("\tRemoved first trial- assuming it's a repeat")
        }
      }
    cat(paste0(", now contains ",length(rawDataThis$trialnum)," trials ",msg))
    if (expi==1 & i==1 & j==1) { #first file of the first subject
        rawData<- rawDataThis
    } else {  #not the first file of the first subject, so combine it with previously-loaded data
        prevColNames<- colnames(rawData)
        newCols <- setdiff( colnames(rawDataThis),prevColNames )
        oldColsNotInNew <- setdiff( prevColNames,colnames(rawDataThis) )
        if (length(newCols) >0) {
          cat( "newCols are:")
          print( paste(newCols,collapse=','))
          for (n in 1:length(newCols)) {#add newCol to old data.frame with dummy value
            newCol = newCols[n]
            rawData[,newCol] <- NA 
            #if (is.numeric(rawDataThis[,newCol]))   #This seems too risky, might forget have -999 values
            #  rawData[,newCol] <- -999 #dummy value
          }
        }
        if (length(oldColsNotInNew) >0)
          for (n in 1:length(oldColsNotInNew)) { #add old col to new data.frame that doesn't have it
            if (n==1) {
              cat("Adding to new data the old columns:")
              print( paste(oldColsNotInNew,collapse=',') )
            }
            oldCol = oldColsNotInNew[n]
            rawDataThis[,oldCol]<- NA #dummy value
            #if (is.numeric(rawData[,oldCol]))  #seems too risky- might forget it is -999
            #  rawDataThis[,oldCol] <- -999 #dummy value
          }
        #Try to merge new data file with already-loaded
        colnamesNew <- colnames(rawDataThis)
        colnamesOld <- colnames(rawData)
		    #colnamesNewMsg <- paste(colnamesNew,collapse=",")
        #colnamesOldMsg <- paste(colnamesOld,collapse=",")
        #writeLines( paste('colnamesNew=',colnamesNewMsg,'\n colnamesOld=', colnamesOldMsg))
		    if ( length(setdiff(colnamesNew,colnamesOld)) >0 )
          writeLines( paste('New columns not in old are ', setdiff(colnamesNew,colnamesOld)) )
        tryCatch( rawData<-rbind(rawData,rawDataThis), #if fail to bind new with old,
                  error=function(e) { #Give feedback about how the error happened
                    cat(paste0("Tried to merge but error:",e))
                    colnamesNewFile <- colnames(rawDataThis)
                    colnamesOldFiles <- colnames(rawData)
                    #colnamesNewFileMsg <- paste(colnamesNewFile,collapse=",")
                    #colnamesOldFilesMsg <- paste(colnamesOldFiles,collapse=",")
                    #writeLines( paste('colnamesNew=',colnamesNewMsg,'\n colnamesOld=', colnamesOldMsg))
                    #c( 'New cols: ', setdiff(colnamesNewFile,colnamesOldFiles) )
                    newCols <- setdiff(colnamesNewFile,colnamesOld)
                    oldColsNotInNew<- setdiff(colnamesOldFiles,colnamesNew)
                    if (length(newCols)>0) {
                      writeLines( paste('New cols not in old: ', paste(newCols,collapse=",") ) ) 
                    }
                    writeLines( paste('Old cols not in new file: ', paste(oldColsNotInNew,collapse=",") ) )        
                    stop(paste0("ERROR merging, error reported as ",e))
                  } 
        )
      } #if not first subject      
    }		
  }
 rawData = subset(rawData, subject != "auto") #get rid of any autopilot data
 #check data counterbalancing of this exp
 source("analysis/helpers/checkCounterbalancing.R")
 checkCombosOccurEqually(rawData, c("subject","trialnum") )
 checkCombosOccurEqually(subset(rawData,task=="T1"), c("targetLeftRightIfOne") )
 checkCombosOccurEqually(rawData, c("condition","leftOrRight") )
 checkCombosOccurEqually(rawData, c("condition","leftOrRight","offsetXYring0") ) #NO?
 checkCombosOccurEqually(rawData, c("numObjects","numTargets","speed") )  
}

dat <-rawData
#end data importation

#If instead of using raw speed, I rank the speed within each numObjects*numTargets, then from that perspective everything should
#be perfectly counterbalanced, because each numObjects*numTargets combination has the same number of speeds tested
#But the rank for a speed depends on what numObjects-numTargets condition it's in. Should be easy with ddply
ordinalSpeedAssign <- function(df) {
#df$speedRank <- rank(df$speed)  #Rank won't work, always wants to break all ties. Whereas I want to preserve ties.
  df$speedRank <- match(df$speed,unique(df$speed))
  df
}
d<- plyr::ddply(dat,.(numObjects,numTargets),ordinalSpeedAssign)
#grouped<- group_by(dat,numObjects,numTargets) #can't get this to work with dpylr but involves something with .  http://stackoverflow.com/questions/22182442/dplyr-how-to-apply-do-on-result-of-group-by
#d<- dplyr::summarise(grouped, speedRank= match(speed,unique(.$speed)))
#dat %>% group_by(numObjects,numTargets) %>% do(match(speed,unique(.$speed)))
#check whether counterbalanced with for each speed list for a particular condition, did each equally often
#Might not be if ran multiple sessions with different speeds
checkCombosOccurEqually(d, c("numObjects","numTargets","speedRank") )

sanityCheckEyeTracking=TRUE
if (sanityCheckEyeTracking) {
  library(ggplot2)
  h<-ggplot(filter(dat,exp=="circleOrSquare_twoTargets"),
            aes(x=maxXdev,y=maxYdev,color=file)) + geom_point() +facet_grid(~subject)  #Have a look at fixation positions
  quartz("circleOrSquare_twoTargets"); show(h)
  h<-ggplot(filter(dat,exp=="offCenter"),
            aes(x=maxXdev)) + geom_histogram()+ facet_grid(~subject) #Have a look at fixation positions
  quartz("offCenter"); show(h)
}
dat$correct = dat$orderCorrect /3
dat$chanceRate= 1 / dat$numObjects

rotX <- function(ch,x) 
{ #rotate each letter of a string ch by x letters thru the alphabet, as long as x<=13
  old <- paste(letters, LETTERS, collapse="", sep="")
  new <- paste(substr(old, 2*x+1, 26*2), substr(old, 1, 26), sep="")
  chartr(old, new, ch)
}
if (anonymiseData) {
  keyFile = paste0('dataPreprocess/',"anonymisationKey.txt")
  if ( !file.exists(keyFile) ) {
  	stop(paste0('The file ',keyFile, ' does not exist!'))
  }
  linesFromFile= readLines(keyFile,warn=FALSE)
  key = as.numeric(linesFromFile[1]) #the key to encrypt the data with
  subjectNotanonymised<- dat$subject
  dat$subject <- rotX(subjectNotanonymised,key) #anonymise subject initials by rotating them by key characters
  print('Mapping from name to anonymised:')
  print(table(subjectNotanonymised,dat$subject))
}
	
#table(d$speedRank,d$numObjects,d$numTargets,d$subject)

#Save anonymised data for loading by doAllAnalyses.R
fname=paste(destinatnDir,destinationName,sep="")
save(dat, file = paste(fname,".RData",sep=""))
write.csv(dat, file = paste(fname,".csv",sep=""))
print(paste("saved data in ",fname,".RData and ",fname,".csv",sep=""))