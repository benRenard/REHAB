#***************************************************************************----
# Uncertainty Estimation ----

#' Residuals Analyzer
#'
#' Analyse the residuals between obs and sim.
#' Residuals are assumed to be realizations from an AR(1) Gaussian process,
#' whose mean mu, standard deviation sigma and lag-1 autocorrelation phi
#' vary according to polynomial regressions with predictors.
#'
#' @param obs numeric vector, observed values.
#' @param sim numeric vector, simulated values.
#' @param x_mu numeric vector, matrix or data frame, predictors for the mean mu.
#' @param x_sigma numeric vector, matrix or data frame, predictors for the sdev sigma.
#' @param x_phi numeric vector, matrix or data frame, predictors for the autocorrelation phi.
#' @param polyFunk function, function used to compute predictor polynomials.
#' @param trans character string, predictor transformation, one of c('uniform','normal','none').
#' @param deg_mu Integer vector, size NCOL(x_mu). Polynomial degree of each predictor
#'    used to compute the mean. Trick: passing a single value as a character will force
#'    mu to remain equal to this value.
#' @param deg_sigma integer vector, size NCOL(x_sigma). Polynomial degree of each predictor
#'    used to compute the standard deviation. Trick: passing a single value as a character will force
#'    sigma to remain equal to this value.
#' @param deg_phi integer vector, size NCOL(x_phi). Polynomial degree of each predictor
#'    used to compute the autocorrelation. Trick: passing a single value as a character will force
#'    phi to remain equal to this value.
#' @param ... additional arguments passed to function [polynomialRegression()]
#' @return An object of class 'anaRes', which is an object of class 'polyRegFit'
#'     augmented with the following components:
#'   \item{obs}{numeric vector, observed values}
#'   \item{sim}{numeric vector, simulated values}
#'   \item{xtrans}{list of size 3, transformed predictors for mu, sigma and phi}
#'   \item{trans}{list of size 3, each element being a list of functions:
#'       transformation function used for each predictor (column) in x_mu, x_sigma and x_phi}
#' @examples
#' w=analyseResiduals(ArdecheRiver$obs,ArdecheRiver$sim)
#' plot(w)
#' @export
analyseResiduals <- function(obs,sim,
                             x_mu=data.frame(sim=sim),
                             x_sigma=data.frame(sim=sim),
                             x_phi=data.frame(dsim=c(0,diff(sim))),
                             polyFunk=LegendrePolynomials,
                             trans=c('uniform','normal','none'),
                             deg_mu=5,deg_sigma=5,deg_phi=5,...){
  transform=match.arg(trans)
  res=obs-sim
  # transform predictors
  tx_mu=NA*x_mu;tx_sigma=NA*x_sigma;tx_phi=NA*x_phi
  trans_mu=lapply(x_mu,getTransform,trans=transform)
  for(i in 1:NCOL(x_mu)){tx_mu[,i]=trans_mu[[i]](x_mu[,i])}
  trans_sigma=lapply(x_sigma,getTransform,trans=transform)
  for(i in 1:NCOL(x_sigma)){tx_sigma[,i]=trans_sigma[[i]](x_sigma[,i])}
  trans_phi=lapply(x_phi,getTransform,trans=transform)
  for(i in 1:NCOL(x_phi)){tx_phi[,i]=trans_phi[[i]](x_phi[,i])}
  # Apply polynomial regression
  out=polynomialRegression(y=res,x_mu=tx_mu,x_sigma=tx_sigma,x_phi=tx_phi,polyFunk=polyFunk,
                           deg_mu=deg_mu,deg_sigma=deg_sigma,deg_phi=deg_phi,...)
  # Add obs, sim and transformation info to polyRegFit object and extend its class
  out$obs=obs;out$sim=sim
  out$xtrans[['mu']]=tx_mu
  out$xtrans[['sigma']]=tx_sigma
  out$xtrans[['phi']]=tx_phi
  out$trans[['mu']]=trans_mu
  out$trans[['sigma']]=trans_sigma
  out$trans[['phi']]=trans_phi
  class(out) <- c('anaRes',class(out))
  return(out)
}

#' Plot Result of Residual Analysis
#'
#' Plot results of a residual analysis.
#' @param x object of class 'anaRes', resulting from a call to function [analyseResiduals()]
#' @param predictorType character string, one of 'transformed' (show transformed predictors),
#'     'raw' (show raw predictors) or 'scaled' (show rescaled predictors).
#' @param alpha numeric in [0,1], transparency of uncertainty interval.
#' @param allPlots, boolean, produce all plots or only the obs vs. sim one?
#' @param ... further plotting arguments passed to or from other methods.
#' @return A list of 4 ggplots (if allPlots=FALSE, only the first one is returned):
#' \enumerate{
#'   \item 'ObsVsSim': predicted values (with uncertainties) and observations
#'   \item 'mu': for each mu predictor, a 3-panel plot with the predictor on the x-axis,
#'       and on the y-axis: estimated mu, predicted values and standardized residuals
#'   \item 'sigma': same as previously but for sigma
#'   \item 'phi': same as previously but for phi
#' }
#' @inherit analyseResiduals examples
#' @export
#' @import ggplot2
plot.anaRes <- function(x,
                        predictorType=c('transformed','raw','scaled'),
                        alpha=0.25,allPlots=FALSE,...){
  xt=match.arg(predictorType)
  if(allPlots){
    out=plot.polyRegFit(x,predictorType=xt,alpha=alpha,...)
  } else {
    out=list()
  }
  # Add a obs vs. sim plot
  DF=data.frame(x=1:NROW(x$sim),sim=x$sim,obs=x$obs,mu=x$pred$mu,sigma=x$pred$sigma)
  g=ggplot(DF)+
    geom_point(aes(.data$x,.data$obs))+
    geom_ribbon(aes(x=.data$x,ymin=.data$sim+.data$mu-1.64*.data$sigma,ymax=.data$sim+.data$mu+1.64*.data$sigma),
                fill='red',alpha=alpha)+
    geom_line(aes(.data$x,.data$sim+.data$mu),col='red')+
    geom_line(aes(.data$x,.data$sim),col='gray')+
    coord_cartesian(ylim=range(c(DF$sim,DF$obs),na.rm=TRUE))+
    labs(x='Index',y='Obs. vs. predicted')+
    theme_bw()
  out[[1]]=g
  names(out)[1] <- 'ObsVsSim'
  return(out)
}

#' Predictive Uncertainty
#'
#' Compute predictive uncertainty around simulated values, after having performed a residual analysis.
#' @param object object of class anaRes, resulting from a call to function [analyseResiduals()]
#' @param sim numeric vector, simulated values
#' @param x  list of size 3 with names 'mu', 'sigma' and 'phi', predictors (matrix or data frame).
#' @param nrep integer, number of replications used to compute parametric and predictive
#'     uncertainties.
#' @param intervalLevel numeric in (0,1), level of the uncertainty intervals.
#' @param ... further arguments passed to or from other methods.
#' @return An object of class 'PredictiveU', which is an object of class 'polyRegPred'
#'     augmented with the following components:
#'   \item{sim}{numeric vector, simulated values}
#'   \item{spag}{numeric matrix, size nrep*length(sim): replications representing predictive uncertainty}
#'   \item{env}{data frame, with columns 'median' (predictive median) and 'low' + 'high' (predictive uncertainty envelope). }
#' @examples
#' w=analyseResiduals(ArdecheRiver$obs[1:366],ArdecheRiver$sim[1:366])
#' u=getUncertainty(w,sim=ArdecheRiver$sim[367:731])
#' plot(u,axisValues=ArdecheRiver$date[367:731])
#' @export
#' @importFrom stats quantile
getUncertainty <-function(object,sim=object$sim,
                          x=list(mu=data.frame(sim=sim),
                                 sigma=data.frame(sim=sim),
                                 phi=data.frame(dsim=c(0,diff(sim)))),
                          nrep=1000,intervalLevel=0.9,...){
  # Transform x and predict
  tx=x
  for(i in 1:length(x)){
    for(j in 1:NCOL(x[[i]])){
      tx[[i]][,j]=object$trans[[i]][[j]](x[[i]][,j])
    }
  }
  out=predict.polyRegFit(object,tx,nrep,...)
  # Complete prediction object
  out$sim=sim
  out$spag=rep(sim,each=NROW(out$predictive))+out$predictive
  probs=c((1-intervalLevel)/2,0.5,1-(1-intervalLevel)/2)
  foo=apply(out$spag,2,quantile,probs=probs)
  env=data.frame(t(foo)[,c(2,1,3)])
  names(env) <- c('median','low','high')
  out$env=env
  class(out) <- c('predictiveUncertainty',class(out))
  return(out)
}

#' Plot Predictive Uncertainty
#'
#' Plot predictive uncertainty estimated after a residual analysis.
#' @param x object of class 'predictiveUncertainty', resulting from a call to function [getUncertainty()]
#' @param axisValues numeric or date vector, what to display on the x-axis.
#' @param labs character vector of size 2, x-axis and y-axis labels.
#' @param col color specification, color used for predictions.
#' @param alpha, numeric in (0,1), transparency of uncertainty interval.
#' @param showRawSim boolean, add raw simulations (before mu-correction) to plot ?
#' @param ... further plotting arguments passed to or from other methods.
#' @return A ggplot showing the simulated values + uncertainties
#' @inherit getUncertainty examples
#' @export
#' @import ggplot2
#' @importFrom stats quantile
plot.predictiveUncertainty <- function(x,
                                       axisValues=data.frame(Index=1:NROW(x$sim)),
                                       labs=c(ifelse(is.null(names(axisValues)),' ',names(axisValues)),'Prediction'),
                                       col='red',alpha=0.25,showRawSim=TRUE,...){
  DF=cbind(x$env,x=as.data.frame(axisValues)[,1],sim=x$sim)
  col='red'
  out=ggplot(DF,aes(x=.data$x))+
    geom_ribbon(aes(ymin=.data$low,ymax=.data$high),fill=col,alpha=alpha)+
    geom_line(aes(y=.data$median),color=col)+
    labs(x=labs[1],y=labs[2])+
    theme_bw()
  if(showRawSim){out=out+geom_line(aes(y=.data$sim))}
  return(out)
}
