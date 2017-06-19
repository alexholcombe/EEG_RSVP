library(splines)
library(ggplot2) 
library(boot) 
library(modelfree)
library(brglm)
library(plyr)
#library(Cairo) #for windows
#Alex Holcombe, started November 2010
#previous versions of this file were called things like lapseRateSearchDataSeparate
source('helpers/psychometricRobustify4.R') #load my custom version of binomfit_lims

#global variables this code expects:  
#threshCriterion, 
#linkf
###end global variables

calcGlmCustomLinkDeviance<- function(fitModel,numCorrect,numTrials,xs,
  				                  	 guessing=guessing,lapsing=l) {
	#Note that this means you'll be using penalized deviance within a lapse rate but comparing standard deviance across lapse rates
	#dataFit <- model.frame(fitModel) #Don't use this because it is scaled.
	
	#numCorrect<<-numCorrect; numTrials<<-numTrials; xs<<-xs; guessing<<-guessing; lapsing<<-lapsing
	modelToCalcDeviance<- suppressWarnings( 
    		binomfit_limsAlex(numCorrect, numTrials, xs, guessing=guessing, lapsing=lapsing, 
    						  initial="glmCustomlink", tryOthers=FALSE)   )

	modelToCalcDeviance$coefficients[1] <- fitModel$fit$coefficients[1]
	modelToCalcDeviance$coefficients[2] <- fitModel$fit$coefficients[2]
	modelCheck <- modelToCalcDeviance
	mu<- suppressWarnings( predict(modelToCalcDeviance$fit,data.frame(x=xs),type="response") )
	#If you call the above without suppressWarnings, then it returns a "prediction from a rank-deficient fit may be misleading" error. After some investigation, I find nothing wrong with the design matrix. So I think it has something unimportant to do with the xs I'm calling it with for the prediction
	
	#Calc raw deviance - null model deviance for each point
	deviances<- binomial()$dev.resids(numCorrect/numTrials,mu,numTrials) 
	deviance<- sum(deviances)
	return(deviance)
}

fitBrglmKludge<- function( df, lapseMinMax, returnAsDataframe, initialMethod, verbosity=0 ) {
#for each possible lapse rate, from lapseMinMax[1] to lapseMinMax[2] in steps of .01
#do probit regression. Return estimate of lapseRate and mean and stdev of underlying Gaussian
#requires countWarnings() function that I wrote and put in psychometricHelpersAlex.R
#If returnAsDataframe, that means return the fit params as data frame. This allows returning the text of the warning.
#Boot can't deal with a dataframe or list, so boot-related functions will call with returnAsDataframe=FALSE.
  if (!is.null(df$chanceRate))
  	chanceRate<- df$chanceRate[1] #assume it's same for all these trials
  else stop("df dataframe passed in must have chanceRate")
  #round min up to the nearest .01
  min01 = ceiling(lapseMinMax[1]*100) /100
  #round max down to the nearest .01
  max01 = floor(lapseMinMax[2]*100) /100
  lapseRates= seq(min01,max01,.01)
  if (min01 != lapseMinMax[1])
    lapseRates= c(lapseMinMax[1],lapseRates) #tack on the minimum 
  if (max01 != lapseMinMax[2])
    lapseRates= c(lapseRates,lapseMinMax[2]) #tack on the minimum 
  deviances = c(); predictors=c()
  parms=c()
  warned=c(); errored=c(); warnMsgs=c(); linkFxs=c(); methods=c()
  for (l in lapseRates)	{
  	#print(paste('trying lapseRate=',l)) #debugOFF
	  cntrl = list(epsilon=1e-6,maxit=10000) 
    iv <- df[,1] #assume first column is independent variable, be it speed or be it tf
	  #problem below is that if warning occurs, we don't get the value returned because of error-catching stupidity
  	fitModel<-countWarnings( 
  				binomfit_limsAlex(df$numCorrect,df$numTrials,iv,link="logit",
  				                  guessing=chanceRate,lapsing=l,control=cntrl,initial=initialMethod)
  			  ) 
	  importantWarns <- attributes(fitModel)$warningMsgs
	  importantWarns<- importantWarns[ which(attributes(fitModel)$warningMsgs!="non-integer #successes in a binomial glm!") ]
  	if (length(importantWarns) >0) {
  		if (verbosity>0) {
	  		cat("fitModel WARNED with: "); cat(  paste(importantWarns, " ") ); 
  		}
  		#"glm.fit: fitted probabilities numerically 0 or 1 occurred" suggests separation or quasi-separation
  		# Was ubiquitous with standard logistic regression, so stopped counting it I think
  		#if (attributes(fitModel)$warningMsgs[[1]] != "glm.fit: fitted probabilities numerically 0 or 1 occurred") {
	  	#	warned<-c(warned,TRUE);   
	  	#	firstWarnMsg<- attributes(fitModel)$warningMsgs[1]
	  	#	warnMsgs<-c(warnMsgs,firstWarnMsg)
	  	#} else warned<-c(warned,FALSE)
	  	warned<-c(warned,TRUE)   
	  	firstWarnMsg<- importantWarns[1]
	  	warnMsgs<-c(warnMsgs,firstWarnMsg)	  	
  	} else warned<-c(warned,FALSE)
  	
  	if(is(fitModel,"error")) { 
  		errored=c(errored,TRUE); print(fitModel)   		
  	}  else errored=c(errored,FALSE)

  	#To compare models with different lapse rates, use deviance which is "up to a constant, minus twice the maximized log-likelihood"
  	#Save all this stuff so know details of winning model at end.
  	method<- fitModel$fit$method
  	if (method=="brglm.fit") {
  		#deviances=c(deviances,fitModel$fit$penalized.deviance)  #use penalized deviance
  		#Deviance, penalized deviance calculated by brglm is wrong because is after data has been rescaled and truncated (chance to 1-lapseRate) 
  		#It's ok for arriving at the best model for a particular lapse rate, but no good for comparing models across lapse rates.
  		#cat('about to calcGlmCustomLinkDeviance')
  		deviance<- calcGlmCustomLinkDeviance(fitModel,df$numCorrect,df$numTrials,iv,guessing=chanceRate,lapsing=l)
  		#cat('leaving calcGlmCustomLinkDeviance with deviance=',deviance)
  		deviances<-c(deviances,deviance)
  	}
  	else
  		deviances<-c(deviances,fitModel$fit$deviance)
  	addsigma=c(fitModel$b,fitModel$sigma)
  	#parms= rbind(parms,fitModel$b)
  	parms= rbind(parms,addsigma)
  	predictors=c(predictors,fitModel$fit)
  	# if (!exists(fitModel$fit$method)) #because glmrob doesn't have that field
  	  # method<- fitModel$method #glmrob puts it here
  	methods<- c(methods,method)

  	linkFxs<- c(linkFxs,fitModel$fit$family$link)
  	if (verbosity>1)
  		cat(paste("With lapseRate=",l,", fitModel$fit$deviance=", fitModel$fit$deviance, "\n"))
  } #END lapseRates loop
  if (length(which(warned)) >0)
  	cat( 'WARNED on lapseRates (',lapseRates[which(warned)], ') '  )
  if (length(which(warned))==length(lapseRates))
  	cat("Got WARNING (possibly nonconvergence) for every lapse rate, but assuming nonconvergence is OK ")
  if (length(which(errored)) >0)
  	cat( 'ERRORed on lapseRates (',lapseRates[which(errored)], ') '  )
  if (length(which(errored))==length(lapseRates))
  	warning("Got ERROR for every lapse rate, so we will CRASH SOON ")
  
  notErrored= which(!errored) #list of indexes of lapse rates for which no error
  #of those where did not error, determine which had best fit
  #unfortunately glmrob doesn't provide any deviance measure, so I would have to calculate it myself
  bestIofNotErrored = which.min( deviances[notErrored] )  #index of lowest deviance, among those which converged
  #cat('bestIofNotErrored=', bestIofNotErrored)
  bestI= notErrored[bestIofNotErrored]
  lapseRate = lapseRates[bestI]
  bestParms = parms[bestI,]
  bestPredictor = predictors[bestI]
  #cat('all parms for this run of lapseRates=',parms)
  mean=bestParms[1]; slope=bestParms[2]; 
  sigma=bestParms[3]
  nWarns<- length(which(warned));  nErrs<- length(which(errored))
  if (is.null(warnMsgs))
  	firstWarn<- "NA"
  else firstWarn<- warnMsgs[[1]]
  #cat('method[bestI]=',methods[bestI]," ") #debugOFF
  #cat('linkFxs[bestI]=',linkFxs[bestI]," ") #debugOFF
  if (returnAsDataframe)
  	dg<- data.frame(mean,slope,chanceRate,lapseRate,sigma,method=methods[bestI],linkFx=linkFxs[bestI],nWarns,nErrs,firstWarn)
  else  #boot wants only a vector back. Can't handle a dataframe. So, cant pass text warning message back because all vec vals
  	dg<- cbind(mean,slope,chanceRate,lapseRate,sigma,nWarns,nErrs) #have to be same type
  
  if (verbosity>1)
  	cat('exiting fitBrglmKludge with:\n'); print(dg)
  return( dg )  	#before I had the following which eventually crapped out inside boot return( list(dg,bestPredictor) )  
}

summarizNumTrials<-function(df) {
  if ( !("correct" %in% names(df)) )
   warning("your dataframe must have a column named 'correct'",immediate.=TRUE)

  numCorrect<-sum(df$correct==1)
  numTrials<- sum(complete.cases(df$correct))
  chanceRate<- sum(df$chanceRate)/numTrials
  
  df= data.frame(numCorrect,numTrials,chanceRate) #,correctY,correctN)
  return(df)
}  

#construct a function to use for one-shot (non-bootstrapping) fit
makeParamFit <- function(iv, lapseMinMax, initialMethod, verbosity=0) {
  #verbosity passed to binomFitChanceRateFromDf
  #iv is independent variable
  fn2 <- function(df) {
    #data comes in one row per trial, but binomFit wants total correct, numTrials
    #so now I have to count number of correct, incorrect trials for each speed
    #assuming there's no other factors to worry about
    if (iv=="speed")
      sumry = ddply(df,.(speed),summarizNumTrials) #also calculates chanceRate
    else if (iv=="tf")
      sumry = ddply(df,.(tf),summarizNumTrials) #also calculates chanceRate
  	#curveFit(sumry$speed,sumry$correct,sumry$numTrials,subjectname,lapsePriors,meanPriors,widthPriors,'MAPEstimation')  
	returnAsDataframe=TRUE #this allows keeping the text of the warning messages. (Boot can't do this)
  	fitParms = fitBrglmKludge(sumry,lapseMinMax, returnAsDataframe,initialMethod,verbosity)
  	#print( paste('fitParms=',fitParms) )
  	return( fitParms )
  }
  return (fn2)
}

#construct a function to use for function fitting and bootstrapping. Will be sent one row per trial
makeParamFitForBoot <- function(chanceRate=0.5,lapseMinMax,verbosity=0) { #default chancePerformanceRate=.5
    #so boot function will provide a random list of idxs. The problem is partialling these out among speeds. Old way of doing it is putting the whole experiment in a single hat, so you can end up with
    #fake datasets that don't even test at certain speeds
    fn2 <- function(df,idxs) {
    	thisData <- df[idxs,]  #perhaps I should not be randomly subsampling from entire dataset
    	#instead, randomly subsample equating number of trials per speed
    	#print('thisData, str='); print(str(thisData))
#    	print(c("idxs=",idxs))
    	#data comes in one row per trial, but binomFit wants total correct, numTrials
    	#so now I have to count number of correct, incorrect trials for each speed
    	#assuming there's no other factors to worry about
    	if ( !("speed" %in% names(df)) )
    	  warning("your dataframe must contain speed as an independent variable",immediate.=TRUE)
	    sumry = ddply(thisData,.(speed),summarizNumTrials)
	    if (verbosity>1) {
	    	print('sumry='); print(sumry)
	    }	
	    if ( length(sumry$speed)==1 )
	    	print('boot has unluckily drawn a bootstrapped experiment with only one speed. Not sure what to do') #actually, bootstrapping should have separate hats for each speed. this is called stratified bstrapping in R terms
	    returnAsDataframe=FALSE #Boot can't handle dataframes. Would allow keeping the text of the warning messages.
        fitParms<- fitBrglmKludge(sumry,chanceRate,lapseMinMax, returnAsDataframe, verbosity)
	    if (verbosity>0) 
  			print( paste('fitParms=',fitParms) )
        return( fitParms )
        #return( c(fitParms$mean,fitParms$slope,fitParms$lapseRate) )
    }
    return( fn2 )
}

makeMyPsychoCorr2<- function(iv) { #Very similar to makeMyPlotCurve below, only for just one x
  fnToReturn<-function(df) {
    #expecting to be passed df with fields:
    # mean, slope, lapseRate, chanceRate, method, linkFx
    df = data.frame(df) #in case it wasn't a dataframe yet
    #set up example model with fake data
    #I don't know why the below didn't work with example01 but it doesn't work
    dh=data.frame(speed=c(.7,1.0,1.4,1.7,2.2),tf=c(3.0,4.0,5.0,6.0,7.0),
                  numCorrect=c(46,45,35,26,32),numTrials=c(48,48,48,48,49))
    dh$lapseRate=df$lapseRate
    dh<-dh
    if(iv=="speed") {
      exampleModel<-suppressWarnings( 
        binomfit_limsAlex(dh$numCorrect, dh$numTrials, dh$speed, link=as.character(df$linkFx), 
                          guessing=df$chanceRate, lapsing=df$lapseRate, initial=as.character(df$method))  #, tryAlts=FALSE  ) 
      ) } else if (iv=="tf") {
        exampleModel<-suppressWarnings( 
          binomfit_limsAlex(dh$numCorrect, dh$numTrials, dh$tf, link=as.character(df$linkFx), 
                            guessing=df$chanceRate, lapsing=df$lapseRate, initial=as.character(df$method))  #, tryAlts=FALSE  ) 
        ) } else {
          print(paste("iv must be either speed or tf, but what was passed was",tf))
        }    
    exampleModel=exampleModel$fit
    #modify example fit, use its predictor only plus parameters I've found by fitting
    exampleModel[1]$coefficients[1] = df$mean
    exampleModel[1]$coefficients[2] = df$slope

    if (iv=="speed") {
      pfit= suppressWarnings( predict( exampleModel, data.frame(x=df$speed), type = "response" ) ) #because of bad previous fit, generates warnings
    } else if (iv=="tf")
      pfit= suppressWarnings( predict( exampleModel, data.frame(x=df$tf), type = "response" ) ) #because of bad previous fit, generates warnings
    
    if (df$method=="brglm.fit" | df$method=="glm.fit") {#Doesn't support custom link function, so had to scale from guessing->1-lapsing manually
      pfit<-unscale0to1(pfit,df$chanceRate,df$lapseRate)
    }
    if(df$numTargets=="2P"){ #Parameters were duplicate of numTargets==1, and p's are corresponding prediction averaged with chance
      pfit<-0.5*(df$chanceRate+pfit)
    }  
    return (pfit)
  }
  return (fnToReturn)
}

makeMyPsychoCorr<- function(iv) { #Very similar to makeMyPlotCurve below, only for just one x
  fnToReturn<-function(df,x) {
    #expecting to be passed df with fields:
    # mean, slope, lapseRate, chanceRate, method, linkFx
    df = data.frame(df) #in case it wasn't a dataframe yet
    #set up example model with fake data
    #I don't know why the below didn't work with example01 but it doesn't work
    dh=data.frame(speed=c(.7,1.0,1.4,1.7,2.2),tf=c(3.0,4.0,5.0,6.0,7.0),
                  numCorrect=c(46,45,35,26,32),numTrials=c(48,48,48,48,49))
    dh$lapseRate=df$lapseRate
    if(iv=="speed") {
      exampleModel<-suppressWarnings( 
        binomfit_limsAlex(dh$numCorrect, dh$numTrials, dh$speed, link=as.character(df$linkFx), 
                          guessing=df$chanceRate, lapsing=df$lapseRate, initial=as.character(df$method))  #, tryAlts=FALSE  ) 
      ) } else if (iv=="tf") {
        exampleModel<-suppressWarnings( 
          binomfit_limsAlex(dh$numCorrect, dh$numTrials, dh$tf, link=as.character(df$linkFx), 
                            guessing=df$chanceRate, lapsing=df$lapseRate, initial=as.character(df$method))  #, tryAlts=FALSE  ) 
        ) } else {
          print(paste("iv must be either speed or tf, but what was passed was",tf))
        }    
    exampleModel=exampleModel$fit
    #modify example fit, use its predictor only plus parameters I've found by fitting
    exampleModel[1]$coefficients[1] = df$mean
    exampleModel[1]$coefficients[2] = df$slope
    pfit<- suppressWarnings( predict( exampleModel, data.frame(x=x), type = "response" ) ) #because of bad previous fit, generates warnings
    if (df$method=="brglm.fit" | df$method=="glm.fit") {#Doesn't support custom link function, so had to scale from guessing->1-lapsing manually
      pfit<-unscale0to1(pfit,df$chanceRate,df$lapseRate)
    }
    if(df$numTargets=="2P"){ #Parameters were duplicate of numTargets==1, and p's are corresponding prediction averaged with chance
      pfit<-0.5*(df$chanceRate+pfit)
    }	
    return (pfit)
  }
return (fnToReturn)
}
  
makeMyPlotCurve4<- function(iv,xmin,xmax,numxs) {#create psychometric curve plotting function over specified domain
  fnToReturn<-function(df) {
  	#expecting to be passed df with fields:
  	# mean, slope, lapseRate, chanceRate,
  	# method, linkFx
    df = data.frame(df) #in case it wasn't a dataframe yet
    #set up example model with fake data
    #I don't know why the below didn't work with example01 but it doesn't work
    dh=data.frame(speed=c(.7,1.0,1.4,1.7,2.2),tf=c(3.0,4.0,5.0,6.0,7.0),
                  numCorrect=c(46,45,35,26,32),numTrials=c(48,48,48,48,49))
    dh$lapseRate=df$lapseRate
    #binomfit_limsAlex(df$numCorrect,df$numTrials,df$speed,link=linkf,guessing=chanceRate,lapsing=l,control=cntrl,initial="brglm.fit")
    if(iv=="speed") {
      exampleModel<-suppressWarnings( 
        binomfit_limsAlex(dh$numCorrect, dh$numTrials, dh$speed, link=as.character(df$linkFx), 
                          guessing=df$chanceRate, lapsing=df$lapseRate, initial=as.character(df$method))  #, tryAlts=FALSE  ) 
      ) } else if (iv=="tf") {
      exampleModel<-suppressWarnings( 
        binomfit_limsAlex(dh$numCorrect, dh$numTrials, dh$tf, link=as.character(df$linkFx), 
                          guessing=df$chanceRate, lapsing=df$lapseRate, initial=as.character(df$method))  #, tryAlts=FALSE  ) 
      ) } else {
        print(paste("iv must be either speed or tf, but what was passed was",tf))
      }    
    exampleModel=exampleModel$fit
    #modify example fit, use its predictor only plus parameters I've found by fitting
    exampleModel[1]$coefficients[1] = df$mean
    exampleModel[1]$coefficients[2] = df$slope
    xs = (xmax-xmin) * (0:numxs)/numxs + xmin
    pfit<- suppressWarnings( predict( exampleModel, data.frame(x=xs), type = "response" ) ) #because of bad previous fit, generates warnings
    if (df$method=="brglm.fit" | df$method=="glm.fit") {#Doesn't support custom link function, so had to scale from guessing->1-lapsing manually
		  pfit<-unscale0to1(pfit,df$chanceRate,df$lapseRate)
	  }
    if ("numTargets" %in% names(fitParms))
      if(df$numTargets=="2P"){ #Parameters were duplicate of numTargets==1, and p's are corresponding prediction averaged with chance
        pfit<-0.5*(df$chanceRate+pfit)
      }	
    #returning the dependent variable with two names because some functions expect one
    #Reason is that want to be able to plot it with same ggplot stat_summary as use for raw
    #data that expects "correct"
    if (iv=="tf") {
      data.frame(tf=xs,pCorr=pfit,correct=pfit)
    } else if (iv=="speed")        
      data.frame(speed=xs,pCorr=pfit,correct=pfit)
  }
  return (fnToReturn)
}

makeMyThreshGetNumerically<- function(iv,threshCriterion) {#create function that can use with ddply once have psychometric curves for each condition
  fnToReturn<-function(df) { #after function has been fit, determine x-value needed for criterion performance
    #So if there's an error, return info about what it errored on. And also indicate there was an error
    #in the dataframe.
    ans<- tryCatch( {
      threshSlop<- threshold_slope(df$correct,df[,iv],criterion= threshCriterion)
      return( data.frame(thresh=threshSlop$x_th, slopeThisCrit=threshSlop$slope, error=FALSE) )
    }, 
                    error = function(e) {
                      cat("\nERROR occurred with")  
                      if ("separatnDeg" %in% names(df))
                        cat(paste(' separatnDeg=',df$separatnDeg[1]),' ') #debugON
                      if ("exp" %in% names(df))
                        cat(paste('exp=',df$exp[1]),' ') #debugON
                      if ("subject" %in% names(df))
                        cat(paste('subject=',df$subject[1]),' ') #debugON
                      if ("numObjects" %in% names(df))
                        cat(paste('numObjects=',df$numObjects[1]),' ') #debugON
                      if ("numTargets" %in% names(df))
                        cat(paste('numTargets=',df$numTargets[1])) #debugON
                      print(e)
                      return( data.frame(thresh=NA, slopeThisCrit=NA, error=TRUE) )
                    }#,
                    #       finally = function(e) { #just return the normal answer
                    #         return( data.frame(thresh=threshSlop$x_th, error=FALSE) )
                    #       }
    )
    
    print(ans)
    return (ans)
  }
}

makeMyThreshGetAnalytically<- function(threshCriterion,linkingFunctionType) { #create function that can use with 
#ddply once have function fits for each condition
#based on knowledge of function fit,
#calculate x-value corresponding to threshCriterion (threshold)
  fnToReturn<-function(df) {
  	 dg= calcThresh(df, linkingFunctionType, chanceRate, threshCriterion) #should be in psychometricHelpersAlex file
  	 dg
  }
}

threshLine <- function(df) {   #should be sent a one-row piece of data frame with threshold speed the last column
	#assumes that df has column "thresh"
	threshes = df$thresh
	speeds=c(0,threshes[1],threshes[1])
	#print('speeds=');print(speeds)
	yMin=0
	corrects=c(threshCriterion,threshCriterion,yMin-.2) #draw down to horizontal axis. The -.2 makes sure it extends into margin
	
	grid<-data.frame(speed=speeds,correct=corrects)
	#print('grid='); print (grid)
	return (grid) 
}

options(warn=0) #needs to be 0, otherwise no way to catch warnings while also letting function to continue 
#options(warn=1)#instead of default zero which waits until top-level function returns, this will print warnings as we go
#ALERT THIS OPTIONS STUFF SHOULD BE INSIDE A FUNCTION, IF IT'S NEEDED AT ALL

########################do bootstrapping####################
#get confidence interval on parameters, so can draw confidence region
#bootstrap for each subset of experiment sent by ddply
makeMyBootForDdply<- function(getFitParmsForBoot,iteratns,confInterval,verbosity=0) {  #create a psychometric curve plotting function over specified domain
	#assumes getFitParmsForBoot has already been constructed
#For the parametric bootstrap it is necessary for the user to specify how the resampling is to be #conducted. The best way of accomplishing this is to specify the function ran.gen which will return a #simulated data set from the observed data set and a set of parameter estimates specified in mle	
	fnToReturn<-function(df) {
  		#send to boot the dataframe piece with one row per trial
  		print(paste('boostrapping with',names(df)[1],'=',df[1,1],names(df)[2],'=',df[1,2]))
  		#lastDfForBoot <<-df
  		b<-boot(df,getFitParmsForBoot,R=iteratns,strata=factor(df$speed)) 
  		print('finished boot call, and boot returned:'); print(b)
		
  		ciMethod= 'perc'  #'bca' don't use until investigate why sometimes get error. #with 'bca' boot method using Christina's data, get this error: Error in bca.ci(boot.out, conf, index[1L], L = L, t = t.o, t0 = t0.o,  :   estimated adjustment 'a' is NA
  		
  		ciMeanWithWarnings<- countWarnings(   boot.ci(b,conf=confInterval,index=1,type= ciMethod)    ) 
  	  	
  		if (length(attributes(ciMeanWithWarnings)$warningMsgs) >0) {
  			print("ciMean boot.ci warned with:"); print(attributes(ciMeanWithWarnings)$warningMsgs); 
  			if (verbosity>0) {
  				cat("but gave value of:"); print(ciMeanWithWarnings);
  			}
  	    }
  	    ciMean <- ciMeanWithWarnings
  	    if (is.null(ciMeanWithWarnings)) #calculating statistic on resamplings always yielded the same value
  			ciMean= data.frame(percent=c(-1,-1,-1,b$t[1,1],b$t[2,1]), bca=c(-1,-1,-1,b$t[1,1],b$t[2,1])) 
 #set both ends of CI to that value. Should take notice of warning

  		ciSlopeWithWarnings<- countWarnings(   boot.ci(b,conf=confInterval,index=2,type=ciMethod)    ) 
  		if (length(attributes(ciSlopeWithWarnings)$warningMsgs) >0) {
  			print("ciSlope boot.ci warned with:"); print(attributes(ciSlopeWithWarnings)$warningMsgs); 
  			if (verbosity>0) {
  				cat("but gave value of:"); print(ciSlopeWithWarnings);
  			}
  	    }
  		
  		ciSlope <- ciSlopeWithWarnings
  	    if (is.null(ciSlopeWithWarnings)) #calculating statistic on resamplings always yielded the same value
  			#lapseRates not successfully bootstrapped, so set up dummy bootstrap value return
  			ciSlope= data.frame(percent=c(-1,-1,-1,b$t[1,2],b$t[2,2]), bca=c(-1,-1,-1,b$t[1,2],b$t[2,2])) 

		slopesEachResampling=b$t[,2]
  		failures = which( is.na(slopesEachResampling) )
		numFailures = length( failures )
		#how many times was NA returned. CI calculation will gag if any
		if (numFailures>0) {
		  print(paste('There were',numFailures,' cases where fitting function returned NaN for the slope'))
		  #b$t[failures,2]=-1
		}
		if (lapseMinMax[1] != lapseMinMax[2]) {

  			ciLapseRateWithWarnings<- countWarnings(   boot.ci(b,conf=confInterval,index=3,type= ciMethod)     ) 
  			if (length(attributes(ciLapseRateWithWarnings)$warningMsgs) >0) {
  				print("ciSlope boot.ci warned with:"); print(attributes(ciLapseRateWithWarnings)$warningMsgs); 
  				if (verbosity>0) {
  					cat("but gave value of:"); print(ciLapseRateWithWarnings);
  				}
  	    	}
  			ciLapseRate <- ciLapseRateWithWarnings
  	    	if (is.null(ciLapseRateWithWarnings)) #calculating statistic on resamplings always yielded the same value
  				#lapseRates not successfully bootstrapped, so set up dummy bootstrap value return
  				ciLapseRate = data.frame(percent=c(-1,-1,-1,b$t[1,3],b$t[2,3]), bca=c(-1,-1,-1,b$t[1,3],b$t[2,3])) 
  			
  		} else { 
  			#lapseRates not bootstrapped, so set up dummy bootstrap value return
  			ciLapseRate= data.frame(percent=c(-1,-1,-1,lapseMinMax[1],lapseMinMax[2]), bca=c(-1,-1,-1,lapseMinMax[1],lapseMinMax[2])) 
  		}  

		#pull confidence interval out of objects produced by boot.ci
		if (ciMethod=='perc') { #then confidence interval is in 'percent' field
  			ciMean <- ciMean$percent[4:5]
  			ciSlope <- ciSlope$percent[4:5]
  			ciLapseRate <- ciLapseRate$percent[4:5]
  		} else if (ciMethod=='bca') { #then confidence interval is in 'bca' field
  			ciMean <- ciMean$bca[4:5]
  			ciSlope <- ciSlope$bca[4:5]  	
  			ciLapseRate <- ciLapseRate$bca[4:5]
  		} else print('unexpected ciMethod')
  		print(c('ciMean$percent[4:5]=',ciMean,' ciSlope$percent[4:5]=',ciSlope))
  		
  		data.frame(meanLo=ciMean[1],meanHi=ciMean[2],slopeLo=ciSlope[1],slopeHi=ciSlope[2],lapserateLo=ciLapseRate[1],lapserateHi=ciLapseRate[2])
  	}
  	return (fnToReturn)
}

makeMyMinMaxWorstCaseCurves<- function(myPlotCurve) {
	fn2<-function(df) { 
		#calculate curves for factorial combination of the confidence interval parameters. Then, plot the most extreme of them all by assigning them to response.inf and response.sup
  		allCombos= expand.grid( mean=c(df$meanLo,df$meanHi), slope=c(df$slopeLo,df$slopeHi), lapseRate=c(df$lapserateLo,df$lapserateHi) )
  		#have to add lapse rates to this dataframe
  		
  		#print(c('allCombos=',allCombos))
  		worstCasePsychometrics= adply(allCombos,1,myPlotCurve)  
  		minmaxCIpsychometrics= ddply(worstCasePsychometrics,.(speed), 
		     myMinMax<- function(df) { data.frame(lower=min(df$correct),
		  									       upper=max(df$correct)) } )
		return(minmaxCIpsychometrics)
	}
	return (fn2)
}



threshStatistic= function(df) {	# Wing ADD20101111		
 thresh=mean(df$thresh)
 sethresh=qnorm(0.975)*sd(df$thresh)/sqrt(as.numeric(length(df$thresh)))
 u95thresh=thresh+sethresh
 l95thresh=thresh-sethresh
 meanSlope=mean(df$slopeAtThresh)
 seSlope=qnorm(0.975)*sd(df$slopeAtThresh)/sqrt(as.numeric(length(df$slopeAtThresh)))
 u95Slope=meanSlope+seSlope
 l95Slope=meanSlope-seSlope
 dg=cbind(thresh,meanSlope,sethresh,seSlope,u95thresh,l95thresh,u95Slope,l95Slope)
 return( dg ) 
}

threshStatisticNum= function(df) {	# Wing ADD20101111		
 thresh=mean(df$thresh)
 sethresh=qnorm(0.975)*sd(df$thresh)/sqrt(as.numeric(length(df$thresh)))
 u95thresh=thresh+sethresh
 l95thresh=thresh-sethresh
 dg=cbind(thresh,sethresh,u95thresh,l95thresh)
 return( dg ) 
}


