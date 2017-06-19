#psychometricHelpersAlex.R
library(robustbase)

threshold_slope <- function( pfit, xfit, criterion) {
	#this function will only work for decreasing curves (that decrease as x-variable increases)
    if( missing("pfit") || missing("xfit") ) {
        stop("Check input. First 2 arguments are mandatory");}
    if( !is.double( criterion ) || length( criterion ) > 1 ) {
        stop( "Criterion level must be scalar" );}
    if( ( criterion < 0 ) || (criterion > 1 )) {
        stop( "Criterion level must be between 0 and 1" );} 
	# Check that the input variables match 
	if ( length( pfit ) != length( xfit ) )
    	stop( 'Length of fitted values pfit must be the same as length of xfit' );
	end
	if ( length(which(pfit>=criterion)) ==0 )
		stop( paste('All values are below criterion,',criterion,' so cannot calculate threshold'))
    end
	if ( length(which(pfit<=criterion)) ==0 )
		stop( paste('All values are above criterion,',criterion,' so cannot calculate threshold') )
    end    
    value <- NULL;
	# threshold
	#xs for which difference between curve and threshold are minimized
    value$x_th <- xfit[ which( abs( pfit - criterion ) == min( abs( pfit - criterion ) ) ) ];
    #note that for step-function psychometrics, both those below and above thresh could be equally distant
    if( length( value$x_th ) > 1 ) {  
		# if there are many point for the same threshold value, then function is flat
		# in this point and slope=0
		diffs = diff(pfit)
	    #print(c("diffs are",diffs))
		#find middle one of 0-slope flat area
		idxsFlat= which(diffs==0)
		err=FALSE
		if (length(idxsFlat)==length(diffs)) {
		  err=TRUE
		  warning('Error: psychometric function completely flat, so impossible to estimate threshold')
		} else { #look for last flat area that's higher than threshold
		  #find rightmost pfit that's higher than threshold
		  aboveCriterionIdxs= which(pfit>criterion)
		  lastAboveCriterion = aboveCriterionIdxs[ length(aboveCriterionIdxs) ]
		  idxsFlat=subset(idxsFlat,idxsFlat<lastAboveCriterion)
		  #find highest one of flat areas that's less than aboveCriterionIdxs
		  idxLast = idxsFlat[length(idxsFlat)] +1 #determine the last one of flat area
		  #cat('idxLast=',idxLast,'idxsFlat=',idxsFlat,'xfit=',xfit,'pfit=',round(pfit,digits=2))
		  if (idxLast==length(xfit)) {
		  	err=TRUE
		  	#print ('idxLast==length(xfit)') #debug
		  	warning("Error: never crossed criterion, it seems and is flat at last part")
		  } else {
		      #print(c('idxLast=',idxLast,' aboveCriterionIdxs=',aboveCriterionIdxs))
		      idxChange = idxLast 
		    }
         }
         if (err) {
        	value$x_th = NA; value$slope = NA
        	return (value)
        }
		#Have to add one because came from diff()
		# print(paste('idxsFlat=',idxsFlat))
        slopeAtLaunchOff =  (abs(xfit[idxChange +1]-xfit[idxChange]))/
                            (abs(pfit[idxChange +1]-pfit[idxChange]))
        #print(c('length of xfit=',length(xfit),"idxChange =", idxChange,'xfit[idxChange]=',xfit[idxChange],' slopeAtLaunchOff=', slopeAtLaunchOff)) #AH 
        value$x_th <-xfit[idxChange]+
                     (abs(criterion-pfit[idxChange]))*slopeAtLaunchOff
        print(value$x_th)
        value$slope <-(pfit[idxChange +1]-pfit[idxChange])/
                      (xfit[idxChange +1]-xfit[idxChange]);       
    }
    else {
# slope
        ind <- which( xfit == value$x_th );
        value$slope <- ( pfit[pmin( ind + 1, length( pfit ) )] - 
                       pfit[pmax( ind - 1, 1 )] ) /
                     ( xfit[pmin( ind + 1, length( xfit ) )] -
                       xfit[pmax( ind - 1, 1 )] );
    }
    return( value );
}

getMyLinkingFunctions <- function(linkf,guessing,lapsing) {
 	linkcall= paste(linkf,"_link_private",sep="")
 	eval( call(linkcall,guessing,lapsing) )  #e.g. logit_link_private(guessing,lapsing)
}
#getMyLinkingFunctions("logit",.5,.01)

calcThresh <- function(df,linkingFunctionType,chanceRate,threshCriterion=.75) {
 #use linking function together with regression parameters to calculate threshold
 #expects dataframe with slope,mean,lapseRate of glm, linkingFunctionType like "logit" or "probit", and 
 #threshCriterion (e.g. .75)
 slope=df$slope; intercept=df$mean; lapseRate=df$lapseRate
 myLinks= getMyLinkingFunctions(linkingFunctionType,chanceRate,lapsing=lapseRate) 
 #print(myLinks)
 #y=ax+b  Determine y modified logit when untransformed = threshCriterion (such as .75)
 #myLinks<<-myLinks; threshCriterion<<-threshCriterion
 critTransformd = myLinks$linkfun(threshCriterion)
 #critTransformd = slope*x + intercept      y=mx+b
 # x = (critTransformd - intercept) / slope
 xCrit = (critTransformd - intercept) / slope
 #calculate slope
 x1=xCrit-.005 #estimate slope from small region around xCrit
 x2=xCrit+.005
 y1transformd= x1*slope + intercept
 y2transformd= x2*slope + intercept
 y1= myLinks$linkinv(y1transformd)
 y2= myLinks$linkinv(y2transformd)
 slopeAtThresh= (y2-y1)/(x2-x1)
 #eventually should also calculate slope, and add it to this function
 data.frame(thresh=xCrit, slopeAtThresh=slopeAtThresh)
}

countWarnings <- function(expr) #make a function that will detect how many warnings were thrown, and return them as well as the value returned by the expression. Value is value returned, other aspects are encoded as R 'attributes' of the variable
#if expr yields error, will print "AH error:" followed by the error
{
    .number_of_warnings <- 0L
    .warningMsgs <- list()
    .warningCalls <- list()
    frame_number <- sys.nframe()

	# insanely, it seems you can't define the functions outside of the call, making this very long
    ans <- withCallingHandlers(expr, warning =  function(w) {
      	#message((paste("mywarning:",w)))
      	assign(".number_of_warnings", .number_of_warnings + 1L,  envir = sys.frame(frame_number))
      	assign(".warningMsgs",c(.warningMsgs,w$message),   envir = sys.frame(frame_number) )  #append to list of warnings
      	assign(".warningCalls",c(.warningCalls,w$call),   envir = sys.frame(frame_number) )  #append to list of warnings
      	invokeRestart("muffleWarning")     
       }, error=function(e) {		cat("AH error:"); print(e)     } )
    
    #message(paste("No. of warnings thrown:", .number_of_warnings))
    #attr(ans,"number_of_warnings")<-.number_of_warnings
    #value= c(list(ans),list(.warningMsgs),list(.warningCalls))
    value <- ans
    if (!is.null(value)) { #this is a major flaw in this scheme it seems, if value is null, assigning
    	#attributes yields an error. Should probably give it some dummy value, or change function to return a list
    	#first element could be null
    	attr(value, 'warningMsgs') <- .warningMsgs
    	attr(value, 'warningCalls') <- .warningCalls
    	#value= list(ans,.number_of_warnings, unlist(.warnings))
    }
    return (value)
}

#I guess if expr is NULL, it won't work?
#
#ans=countWarnings(fool())
#ans=countWarnings(foo())
#
#attributes(ans)$warningMsgs  #recover warning messages
#numWarns= length( attributes(ans)$warningMsgs ) #recover number of warnings
#attributes(ans)$warningCalls  #recover list of functions that generated the warnings
#attributes(ans)$warningCalls[1] #recover first one

###################################################################################################
################brglm does not support custom link functions. 
#So I'll have to scale response range -> 0,1. Then pass with logit. Then take predicted probabilities and scale them back
#function testing
scaleTo0to1 <- function(data,chanceRate,lapseRate)
{
	max <- 1-lapseRate
	min <- chanceRate
	ans <- (data-min) / (max-min)
	ans[which(ans>1)] <-1  #Enforce ceiling of 1. Otherwise glm will reject the data
	ans[which(ans<0)] <-0  #Enforce floor of 0. Otherwise glm will reject the data
	return (ans)
}
unscale0to1 <- function(regressnOutput,chanceRate,lapseRate)
{
	max <- 1-lapseRate
	min <- chanceRate
	ans <- min + regressnOutput*(max-min)
	return(ans)
} #unscale0to1( predict(br,type="response"),guessing,lapsing )
###################################################################################################

binomfit_limsAlex <- function(r,m,x,p=1, link="logit", guessing=0, lapsing=0, K=2,
								control=list(),initial="brglm.fit",tryOthers=TRUE) {

#

# The function fits a binomial generalised linear model with fixed guessing and lapsing rates.

# It is based closely on a function of similar name in the modelfree package
# INPUT

#
# r    - number of successes at points x
# m    - number of trials at points x 
# x    - stimulus levels
#

# OPTIONAL INPUT

#

# p    - degree of the polynomial; default is p = 1 

# link - name of the link function; default is "logit"

# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
# K    - power parameter for Weibull and reverse Weibull link; default is 2
# control - control parameters to pass to fitting algorithm

# initial - initial fitting method to try, "brglm.fit", "glmrob", "glmCustomlink", or "glm.fit"

# OUTPUT
# 
# Object with 2 components: 
# b - vector of estiamted coefficients for the linear part 
# fit - glm object to be used in evaluation of fitted values


# MAIN PROGRAM

# First 3 arguments are mandatory

    # First 3 arguments are mandatory
    if( missing("x") || missing("r") || missing("m") ) {
        stop("Check input. First 3 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS

    checkdata<-list();

    checkdata[[1]] <- x;

    checkdata[[2]] <- r

    checkdata[[3]] <- m

    checkinput( "psychometricdata", checkdata );

    rm( checkdata )

    checkinput( "linkfunction", link );

    pn <- list()
	  pn[[1]] <- p
	  pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( "guessingandlapsing", c( guessing, lapsing ) );
    if( link == "weibull"  || link == "revweibull") {
    	checkinput( 'exponentk', K );
    	}
	if( initial!="brglm.fit" & initial!="glmrob" & initial!="glmCustomlink" & initial!="glm.fit" )
		stop(paste("Unexpected value of initial-",initial))

# GLM settings

    glmdata <- data.frame( cbind( r/m , m , x ) );

    names( glmdata ) <- c( "resp", "m", "x" );
    #to pass list of parameters to function accepting only individual named parameters, must use do.call

   	ctrl<- do.call("glm.control",control)
   	cntrl.brglm<- do.call("brglm.control",control) 


# formula

    glmformula <- c( "resp ~ x" );

    if( p > 1 ) {

        for( pp in 2:p ) {

            glmformula <- paste( glmformula, " + I(x^", pp,")", sep = "");

        }

    }

    fit <- NULL;



# GLM fit
	if( link != "logit" && link != "probit" && link != "loglog" && link != "comploglog" &&
        link != "weibull"  && link != "revweibull" ) { #custom link function
       linkfun <- link
     }  else {
               linkfun <- paste( link, "_link_private", sep = "" );
              }


   if( linkfun != "weibull_link_private" && linkfun != "revweibull_link_private" ) {

    	
    	handleErr <- function(e) {
    		eventType = class(e)[2]
    		cat(paste(eventType,"in AH function '", e$call, "': '",e$message,"'\n"))
    		e
    	  }

	method <- initial
	#brglm.fit requires rescaling because brglm doesn't support custom link functions
	#Others do support custom link functions
	#cat("Trying method",method,"\n") 
	if (method=="brglm.fit" | method=="glm.fit") {
		#"glm.fit means use brglm, but have it simulate standard. Only diff is have to scale data because
		#brglm doesn't support custom link function
		glmdataScaled <- glmdata #have to rescale because brglm doesn't support custom link function
		glmdataScaled$resp <- scaleTo0to1(glmdata$resp,guessing,lapsing) 
		fit<- brglm(glmformula, data=glmdataScaled,weights=m,method=method,family=binomial( "logit" ),
					x=T,y=T,control.glm=ctrl, control.brglm=cntrl.brglm) 
		#cat("fit="); print(fit)
	}
	if (method=="glmCustomlink") {
   	    assign("last.warning", NULL, envir = baseenv())  #clear warnings buffer 
   	    #linkfun<<-linkfun; guessing<<-guessing; lapsing<<-lapsing; ctrl<<-ctrl; glmformula<<-glmformula; glmdata<<-glmdata;
  		fit<- glm( glmformula, data = glmdata, weights = m, family = binomial( eval( call( linkfun, guessing, lapsing ) ) ),x=T,y=T,control=ctrl )
		fit$family$link<-"logit" #is something like "logit( c(0.5, 1) )" which will make myPlotCurve crash, so needed to change it
		#Other algorithms return "logit"
	}
	if (method=="glmrob") {
	 	print("Trying glmrob robust")
	 	stop("glmrob not supported for multiple lapse rates, because doesn't provide deviance to compare, 
	 	       you'll have to add code to calculate it manually")
		linkfunRobust="logit_link" #seems better than logit_link_private but dunno why
	   	ctrl<-glmrobMqle.control(maxit=1000) 
	   	assign("last.warning", NULL, envir = baseenv())  #clear warnings buffer  
	 	#Need to fix below to take into account guessing and lapsing rates
	 	fit<- glmrob( glmformula, data = glmdata, weights=m, family=binomial( link="logit" ),x=T,y=T,control=ctrl )  	
	    #unfortunately, glmrob doesn't return deviance, so no way to compare it with other methods
	    #About not returning deviance: https://stat.ethz.ch/pipermail/r-sig-robust/2010/000308.html
	    #that's why I have to test for is.null(fit$call)
	} 	
	
	#test whether need to try another algorithm because above fit attempt crapped out
   	if (is.null(fit) | is.null(fit$call) | class(fit)[2]=="error" & tryOthers) #we'll have to try another algorithm
    	{   
    		cat("Due to previous failure of method",method,", trying conventional glm") #debugON
    		stop("temp error to be removed!") #debugON
    	    ctrl<-glm.control(maxit=1000,trace=FALSE)
    	    assign("last.warning", NULL, envir = baseenv())  #clear warnings buffer 
  			fit = glm( glmformula, data = glmdata, weights = m, 
  					   family = binomial( eval( call( linkfun, guessing, lapsing ) ) ),x=T,y=T,control=ctrl )
        	numWarnings=length(warnings())
			if (numWarnings>0){ 
	    			print("binomfit_limsAlex glm warnings:"); print(warnings())
			}                    
    	}
    	if (is.null(fit) | class(fit)[2]=="error" & tryOthers) #we'll have to try another algorithm
    	{   
    		if (link=='probit')
    			altLink='logit'
    		else if (link=='logit')
    			altLink='probit'
    		else stop("unexpected link function specified in the first instance (link parameter of this function)")
    		
    		print(paste("retrying glm, with",altLink, "instead of",link))
    	    ctrl<-glm.control(maxit=1000,trace=FALSE)
    	    linkfun=paste(altLink, "_link_private", sep = "" )
    	    assign("last.warning", NULL, envir = baseenv())  #clear warnings buffer 
  			fit = glm( glmformula, data = glmdata, weights = m, family = 
  						binomial( eval( call( linkfun, guessing, lapsing ) ) ),x=T,y=T,control=ctrl )    
        	numWarnings=length(warnings())
			if(numWarnings>0){ 
	    		cat(paste("binomfit_limsAlex glm(",altLink," warnings:")); print(warnings())
			}                    
    	} 	    
    	if (is.null(fit) | class(fit)[2]=="error") #we'll have to try another algorithm
    	   print("all fitting algorithms tried yielded errors, YOU IN TROUBLE")
     } #end if not weibull

    else {  #if weibull, have to have K argument
        fit <- glmrob( glmformula, data = glmdata, weights = m, x=T,y=T,
                    family = binomial( eval( call( linkfun, K, guessing, lapsing ) ) ) )
     }
    value <- NULL       
    value$b <- fit$coeff
    value$fit <- fit
    value$sigma <- sqrt(fit$deviance/fit$df.residual)
    return( value );

}