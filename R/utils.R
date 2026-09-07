#' Rescaled Logistic Function
#'
#' A logistic function rescaled between -1 and 1, L(x)=(exp(x)-1)/(exp(x)+1)
#' @param x numeric, values at which the function is evaluated
#' @return a numeric equal to L(x)
#' @examples
#' x=seq(-10,10,0.01)
#' y=rescaledLogistic(x)
#' plot(x,y)
#' @export
rescaledLogistic <- function(x){
  u=exp(x)
  out=(u-1)/(u+1)
  return(out)
}

#' Predictor Transformation
#'
#' Get the transformation function that can be applied to a predictor
#'
#' @param values numeric vector, values used to establish the transformation. Typically,
#'     the predictor values that were used for calibration
#' @param trans character, the requested transformation, one of 'none', 'normal' or 'uniform'
#' @return a function that can be applied to some predictor values
#'     (typically, the ones used for prediction)
#' @examples
#' n=1000
#' values=rlnorm(n)
#' f=getTransform(values,'uniform')
#' plot(values,f(values))
#' hist(f(values))
#' x=rnorm(n,0,10)
#' plot(x);points(values,col='red')
#' plot(x,f(x))
#' hist(f(x))
#' @export
#' @importFrom stats qnorm rnorm rlnorm
getTransform <- function(values,trans=c('none','normal','uniform')){
  tr=match.arg(trans)
  if(tr=='none'){
    transFunk <- function(x){return(x)}
  } else {
    z=sort(values[!is.na(values)])
    n=sum(!is.na(values))
    transFunk <- function(x){
      r=findInterval(x,z)
      r[r<=0]=1
      out=(r-0.5)/n
      if(tr=='normal'){out=qnorm(out)}
      return(out)
    }
  }
  return(transFunk)
}

#' Nash-Sutcliffe Efficiency
#'
#' Compute the Nash-Sutcliffe Efficiency (NSE)
#'
#' @param obs numeric vector, observed values
#' @param sim numeric vector, simulated values
#' @return a numeric value equal to the NSE
#' @examples
#' n=50
#' obs=rnorm(n,0,4)
#' sim=obs+rnorm(n,1,1)
#' NSE(obs,sim)
#' plot(obs);lines(sim)
#' @export
NSE <- function(obs,sim){
  num=sum((obs-sim)^2,na.rm=TRUE)
  m=mean(obs,na.rm=TRUE)
  den=sum((obs-m)^2,na.rm=TRUE)
  out=1-num/den
  return(out)
}
