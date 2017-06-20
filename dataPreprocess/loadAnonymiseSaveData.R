#expecting current working directory to be top level of this git-indexed project, and this file to be in top level - dataPreprocess/

#Set to your top-level directory
setwd("/Users/alexh/Documents/attention_tempresltn/EEGRSVP/EEG_RSVP")

destinatnDir<-"dataAnonymized/" #where the anonymized data will be exported to
expNameForFile<-"EEG_RSVP"
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
      	    		read.table(fname,sep='\t',header=TRUE, strip.white=TRUE), #strip whitespace because otherwise task has whitespaces before and after it 
      	    		error=function(e) { 
      	    	   			stop( paste0("ERROR reading the file ",fname," :",e) )
          		 } )
      #rawDataLoad$file <- file  Shouldn't do this because first the subject name and data stamp part would need to be anonymized
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
    if (i==1 & j==1) { #first file of the first subject
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
checkCombosOccurEqually(subset(rawData,task=="T1"), c("subject","targetLeftRightIfOne") )
checkCombosOccurEqually(rawData, c("subject","task") )

dat <-rawData
#end data importation

rotX <- function(ch,x) 
{ #rotate each letter of a string ch by x letters thru the alphabet, as long as x<=13
  old <- paste(letters, LETTERS, collapse="", sep="")
  new <- paste(substr(old, 2*x+1, 26*2), substr(old, 1, 26), sep="")
  chartr(old, new, ch)
}
#to improve on this, perhaps generate a random list of unique initials that are as long as the number of Ss

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

#Save anonymised data for loading by doAllAnalyses.R
fname=paste(destinatnDir,expNameForFile,sep="")
save(dat, file = paste(fname,".RData",sep=""))
write.csv(dat, file = paste(fname,".csv",sep=""))
print(paste("Saved ",ifelse(anonymiseData, "anonymised", "")," data in ",fname,".RData and ",fname,".csv",sep=""))
