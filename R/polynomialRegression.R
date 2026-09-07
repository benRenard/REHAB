#***************************************************************************----
# Polynomials ----

#' Standard Polynomials
#'
#' Compute standard Polynomials (i.e. monomials) from degree 0 up to degree n
#' @param x numeric vector, values at which polynomials are computed.
#' @param n Integer, maximum polynomial degree.
#' @return A data frame with n+1 columns, where column k is the standard polynomial
#'     with degree k-1.
#' @examples
#' x=seq(-1,1,0.01)
#' DF=StandardPolynomials(x,5)
#' plot(NA,xlab='x',ylab='Hn(x)',xlim=range(x),ylim=range(DF))
#' for(i in 1:NCOL(DF)){lines(x,DF[,i],col=i)}
#' @export
StandardPolynomials <- function(x,n=3){
  if(n<0){stop('n can not be negative.',call.=FALSE)}
  out=data.frame(H0=1+0*x)
  if(n==0){return(out)}
  for(i in 1:n){
    out=cbind(out,x^i)
    names(out)[i+1] <- paste0('P',i)
  }
  return(out)
}

#' Hermite Polynomials
#'
#' Compute Hermite Polynomials from degree 0 up to degree n
#' @param x numeric vector, values at which polynomials are computed.
#' @param n Integer, maximum polynomial degree.
#' @return A data frame with n+1 columns, where column k is the Hermite polynomial
#'     with degree k-1.
#' @examples
#' x=seq(-3,3,0.01)
#' DF=HermitePolynomials(x,5)
#' plot(NA,xlab='x',ylab='Hn(x)',xlim=range(x),ylim=range(DF))
#' for(i in 1:NCOL(DF)){lines(x,DF[,i],col=i)}
#' @export
HermitePolynomials <- function(x,n=3){
  if(n<0){stop('n can not be negative.',call.=FALSE)}
  out=data.frame(H0=1+0*x)
  if(n==0){return(out)}
  out=cbind(out,H1=x)
  if(n==1){return(out)}
  for(i in 2:n){
    foo=x*out[,i]-(i-1)*out[,i-1]
    out=cbind(out,foo)
    names(out)[i+1] <- paste0('H',i)
  }
  return(out)
}

#' Legendre Polynomials
#'
#' Compute Legendre Polynomials from degree 0 up to degree n
#' @param x numeric vector, values at which polynomials are computed.
#' @param n Integer, maximum polynomial degree.
#' @return A data frame with n+1 columns, where column k is the Legendre polynomial
#'     with degree k-1.
#' @examples
#' x=seq(-1,1,0.01)
#' DF=LegendrePolynomials(x,5)
#' plot(NA,xlab='x',ylab='Ln(x)',xlim=range(x),ylim=range(DF))
#' for(i in 1:NCOL(DF)){lines(x,DF[,i],col=i)}
#' @export
LegendrePolynomials <- function(x,n=3){
  if(n<0){stop('n can not be negative.',call.=FALSE)}
  out=data.frame(L0=1+0*x)
  if(n==0){return(out)}
  out=cbind(out,L1=x)
  if(n==1){return(out)}
  for(i in 2:n){
    foo=((2*i-1)/i)*x*out[,i]-((i-1)/i)*out[,i-1]
    out=cbind(out,foo)
    names(out)[i+1] <- paste0('L',i)
  }
  return(out)
}

#***************************************************************************----
# Polynomial Regression ----

#' Polynomial Regression
#'
#' Fit a polynomial regression for the mean mu, the standard deviation sigma
#' and the autocorrelation phi an AR(1) a Gaussian distribution
#'
#' @param y numeric vector, predictand.
#' @param x_mu numeric vector, matrix or data frame, predictors for the mean mu.
#' @param x_sigma numeric vector, matrix or data frame, predictors for the sdev sigma.
#' @param x_phi numeric vector, matrix or data frame, predictors for the autocorrelation phi.
#' @param polyFunk function, function used to compute predictor polynomials.
#' @param deg_mu Integer vector, size NCOL(x_mu). Polynomial degree of each predictor
#'    used to compute the mean. Trick: passing a single value as a character will force
#'    mu to remain equal to this value.
#' @param deg_sigma integer vector, size NCOL(x_sigma). Polynomial degree of each predictor
#'    used to compute the standard deviation. Trick: passing a single value as a character will force
#'    sigma to remain equal to this value.
#' @param deg_phi integer vector, size NCOL(x_phi). Polynomial degree of each predictor
#'    used to compute the autocorrelation. Trick: passing a single value as a character will force
#'    phi to remain equal to this value.
#' @param link_mu function, link function used for the mean mu.
#' @param link_sigma function, link function used for the standard deviation sigma.
#' @param link_phi function, link function used for the autocorrelation phi.
#' @param useApproach2 logical, if TRUE "approach 2" is used (based on center-scaled residuals), otherwise approach 1 is used
#' @param optim_control list, list of controls to be passed to the optimizer, see ?optim.
#' @param optim_method character string, optimization method, see ?optim.
#' @return An object of class 'PolyRegFit', which is a list with the following components:
#'   \item{theta}{numeric vector, estimated parameter vector},
#'   \item{coeff}{list of size 3, subvector of theta corresponding to
#'       the polynomial coefficients for mu, sigma and phi}
#'   \item{pred}{list of size 3, estimated mean mu, sdev sigma and autocorrelation phi (same size as y)}
#'   \item{residuals}{numeric vector, standardized (centered-scaled) residuals (same size as y)}
#'   \item{poly}{list of size 3, polynomials used to compute mu and sigma}
#'   \item{xscaled}{list of size 3, scaled predictors for mu, sigma and phi, with scaling depending
#'       on polyFunk (Hermite and Standard: center-scale, Legendre: send between -1 and 1)}
#'   \item{scaler}{list of size 3, each element being a list of functions:
#'       scaling function used for each predictor in x_mu, x_sigma and x_phi}
#'   \item{scalerInv}{list of size 3, each element being a list of functions:
#'       inverse scaling function for each predictor in x_mu, x_sigma and x_phi}
#'   \item{scalerPar}{list of size 2, parameters of each scaling function}
#'   \item{optim}{list, optimisation results returned by [optim()]}
#'   \item{y}{numeric vector, predictand}
#'   \item{x}{list of size 3, predictors for mu, sigma and phi}
#'   \item{deg}{list of size 3, polynomial degree associated to each predictor for computing mu, sigma and phi}
#'   \item{polyFunk}{function, function to compute polynomials}
#' @examples
#' n=1000
#' predictor=rnorm(n)
#' y=rnorm(n,mean=2*sin(predictor),sd=0.1*(abs(predictor)+0.5))
#' plot(predictor,y)
#' reg=polynomialRegression(y,predictor)
#' plot(reg)
#' @export
#' @importFrom stats sd lm na.exclude residuals optim
polynomialRegression <- function(y,x_mu,x_sigma=x_mu,x_phi=x_mu,
                                 polyFunk=LegendrePolynomials,
                                 deg_mu=rep(3,NCOL(x_mu)),deg_sigma=rep(2,NCOL(x_sigma)),deg_phi=rep(0,NCOL(x_phi)),
                                 link_mu=identity,link_sigma=exp,link_phi=rescaledLogistic,
                                 useApproach2=FALSE,
                                 optim_control=list(maxit=100000),
                                 optim_method=c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")){
  deg=list(mu=deg_mu,sigma=deg_sigma,phi=deg_phi)
  coeff=pred=x=xscaled=poly=scaler=scalerInv=scalerPar=list()
  n=length(y)
  for(what in c('mu','sigma','phi')){
    if(what=='mu'){
      xDF=as.data.frame(x_mu)
    } else if(what=='sigma') {
      xDF=as.data.frame(x_sigma)
    } else {
      xDF=as.data.frame(x_phi)
    }
    x[[what]]=xDF
    p=NCOL(xDF)
    n=NROW(xDF)
    # Rescale x by either centering-scaling (Hermite) or sending to [-1,1] (Legendre)
    scaler[[what]]=scalerInv[[what]]=scalerPar[[what]]=vector('list',p)
    xscaled[[what]]=NA*xDF
    for(i in 1:p){
      if(identical(polyFunk,HermitePolynomials) | identical(polyFunk,StandardPolynomials)){
        m=mean(xDF[,i],na.rm=TRUE)
        s=sd(xDF[,i],na.rm=TRUE)
        scalerPar[[what]][[i]]=c(m,s)
        scaler[[what]][[i]] <- function(z){(z-m)/s}
        scalerInv[[what]][[i]] <- function(z){m+z*s}
      } else if(identical(polyFunk,LegendrePolynomials)){
        mini=min(xDF[,i],na.rm=TRUE)
        maxi=max(xDF[,i],na.rm=TRUE)
        scalerPar[[what]][[i]]=c(mini,maxi)
        scaler[[what]][[i]] <- function(z){2*(z-mini)/(maxi-mini)-1}
        scalerInv[[what]][[i]] <- function(z){mini+0.5*(z+1)*(maxi-mini)}
      }
      xscaled[[what]][,i]=scaler[[what]][[i]](xDF[,i])
    }
    # Assemble polynomials
    poly[[what]]=polyReg_getPolynomials(xscaled[[what]],polyFunk,deg[[what]])
    if(what=='mu'){
      if(is.character(deg[[what]])){
        # mu is fixed
        coeff[[what]]=NULL
        if(identical(link_sigma,exp)){
          lres=log(abs(y-as.numeric(deg[[what]])))
        } else {
          lres=abs(y-as.numeric(deg[[what]]))
        }
      } else {
        # Apply standard linear regression to get first estimate of mu coefficients
        w1=lm(y~as.matrix(poly[[what]])-1,na.action=na.exclude)
        coeff[[what]]=w1$coefficients
        if(identical(link_sigma,exp)){
          lres=log(abs(residuals(w1)))
        } else {lres=abs(residuals(w1))}
      }
    } else if(what=='sigma') {
      if(is.character(deg[[what]])){
        # sigma is fixed
        coeff[[what]]=NULL
      } else {
        # Apply standard linear regressions to get first estimate of sigma coefficients
        mask=is.finite(lres) & !(is.na(lres))
        w2=lm(lres[mask]~as.matrix(poly[[what]][mask,])-1)
        coeff[[what]]=w2$coefficients
      }
    } else{
      if(is.character(deg[[what]])){
        # phi is fixed
        coeff[[what]]=NULL
      } else {
        # start at 0 for autocorrelation
        coeff[[what]]=rep(0,NCOL(poly[[what]]))
      }
    }
  }
  # Minimize negative log-likelihood
  theta0=c(coeff[['mu']],coeff[['sigma']],coeff[['phi']])
  method=match.arg(optim_method)
  op=optim(par=theta0,fn=polyReg_nll,y=y,
           poly_mu=poly[['mu']],poly_sigma=poly[['sigma']],poly_phi=poly[['phi']],
           link_mu=link_mu,link_sigma=link_sigma,link_phi=link_phi,
           useApproach2=useApproach2,
           method=method,control=optim_control,hessian=TRUE)
  theta=op$par
  p_mu=NCOL(poly[['mu']])
  p_sigma=NCOL(poly[['sigma']])
  p_phi=NCOL(poly[['phi']])
  if(p_mu>0){
    coeff[['mu']]=theta[1:p_mu]
    pred[['mu']]=link_mu(as.matrix(poly[['mu']]) %*% coeff[['mu']])
  } else {
    coeff[['mu']]=NULL
    pred[['mu']]=rep(as.numeric(row.names(poly[['mu']])),n)
  }
  if(p_sigma>0){
    coeff[['sigma']]=theta[p_mu+(1:p_sigma)]
    pred[['sigma']]=link_sigma(as.matrix(poly[['sigma']]) %*% coeff[['sigma']])
  } else {
    coeff[['sigma']]=NULL
    pred[['sigma']]=rep(as.numeric(row.names(poly[['sigma']])),n)
  }
  if(p_phi>0){
    coeff[['phi']]=theta[p_mu+p_sigma+(1:p_phi)]
    pred[['phi']]=link_phi(as.matrix(poly[['phi']]) %*% coeff[['phi']])
  } else {
    coeff[['phi']]=NULL
    pred[['phi']]=rep(as.numeric(row.names(poly[['phi']])),n)
  }
  res=(y-pred[['mu']])/pred[['sigma']]
  innov=0*res
  if(useApproach2){
    innov[2:n]=res[2:n]-pred[['phi']][2:n]*res[1:(n-1)]
  } else {
    innov[2:n]=(y[2:n]-pred[['phi']][2:n]*y[1:(n-1)]-pred[['mu']][2:n])/pred[['sigma']][2:n]
  }
  out=list(
    # Main results
    theta=theta,coeff=coeff,pred=pred,residuals=res,innovations=innov,
    # Advanced results
    poly=poly,
    xscaled=xscaled,scaler=scaler,scalerInv=scalerInv,scalerPar=scalerPar, #
    optim=op,
    # Copy of input arguments
    y=y,x=x,deg=deg,polyFunk=polyFunk,useApproach2=useApproach2,
    link_mu=link_mu,link_sigma=link_sigma,link_phi=link_phi)
  class(out) <- 'polyRegFit'
  return(out)
}

#' Polynomial Regression Prediction
#'
#' Predict from a polynomial regression for both the mean, the standard deviation
#' and the autocorrelation of a AR(1) Gaussian distribution.
#' @param object object of class polyRegFit, resulting from a call to function [polynomialRegression()]
#' @param newdata  list of size 3 with names 'mu', 'sigma' and 'phi', predictors (matrix or data frame).
#'     If NULL, the predictors used for calibration are used.
#' @param nsim integer, number of replications used to compute parametric and predictive
#'     uncertainties. Uncertainties are skipped if nsim<=0
#' @param ... further arguments passed to or from other methods.
#' @return An object of class 'PolyRegPred', which is a list with the following components:
#'   \item{pred}{list of size 3, estimated mean mu, sdev sigma and autocorrelation phi}
#'   \item{replicates}{list of size 3, nsim replicates of mu, sigma and phi,
#'       representing parametric uncertainty (size nsim*NROW(newdata[['mu']]))}
#'   \item{predictive}{numeric matrix, nsim replicates of predictand values, representing
#'       total predictive uncertainty (size nsim*NROW(newdata[['mu']]))}
#'   \item{x}{list of size 3, predictors for mu, sigma and phi}
#'   \item{xscaled}{list of size 3, scaled predictors for mu, sigma and phi}
#'   \item{poly}{list of size 3, polynomials used to compute mu, sigma and phi}
#' @examples
#' n=1000
#' predictor=rnorm(n)
#' y=rnorm(n,mean=2*sin(predictor),sd=0.1*(abs(predictor)+0.5))
#' reg=polynomialRegression(y,predictor)
#' newpred=rnorm(100,1.5,1)
#' res=predict(reg,newdata=list(mu=newpred,sigma=newpred,phi=newpred))
#' plot(res)
#' @export
#' @importFrom mvtnorm rmvnorm
predict.polyRegFit <-function(object,newdata=NULL,nsim=1000,...){
  poly=pred=xscaled=list()
  if(is.null(newdata)){newdat=object$x} else {newdat=newdata}
  if(!identical(names(newdat),c('mu','sigma','phi'))){stop('newdata: should be a named list with names mu, sigma and phi',call.=FALSE)}
  if(NROW(newdat[['mu']])!=NROW(newdat[['sigma']])){stop('newdata: number of rows should be identical for mu and sigma',call.=FALSE)}
  if(NROW(newdat[['mu']])!=NROW(newdat[['phi']])){stop('newdata: number of rows should be identical for mu and phi',call.=FALSE)}
  for(what in c('mu','sigma','phi')){
    newdat[[what]]=as.data.frame(newdat[[what]])
    xDF=newdat[[what]]
    if(NCOL(xDF)!=NCOL(object$x[[what]])){stop('newdata: number of columns does not match that of calibration data x',call.=FALSE)}
    p=NCOL(xDF)
    n=NROW(xDF)
    # Rescale x
    xscaled[[what]]=NA*xDF
    for(i in 1:p){xscaled[[what]][,i]=object$scaler[[what]][[i]](xDF[,i])}
    # Assemble polynomials
    poly[[what]]=polyReg_getPolynomials(xscaled[[what]],object$polyFunk,object$deg[[what]])
    # Get predicted mu/sigma/phi
    if(is.character(object$deg[[what]])){
      pred[[what]]=rep(as.numeric(row.names(poly[[what]])),n)
    } else {
      pred[[what]]=as.matrix(poly[[what]]) %*% object$coeff[[what]]
    }
    if(what=='mu'){pred[[what]]=object$link_mu(pred[[what]])}
    if(what=='sigma'){pred[[what]]=object$link_sigma(pred[[what]])}
    if(what=='phi'){pred[[what]]=object$link_phi(pred[[what]])}
  }
  # Get uncertainties
  mus=sigmas=predictive=NULL
  if(nsim>0){
    # Parametric uncertainty
    C=tryCatch(solve(object$optim$hessian),error = function(e){NULL})
    if(is.null(C)){
      replicates=matrix(rep(object$theta,nsim),nsim,length(object$theta),byrow=TRUE)
    } else if (!isSymmetric(C)){
      replicates=matrix(rep(object$theta,nsim),nsim,length(object$theta),byrow=TRUE)
    } else if (any(is.na(C))){
      replicates=matrix(rep(object$theta,nsim),nsim,length(object$theta),byrow=TRUE)
    } else if (any(diag(C)<=0)){
      replicates=matrix(rep(object$theta,nsim),nsim,length(object$theta),byrow=TRUE)
    } else if (any(eigen(x=C)$values<=0)){
      replicates=matrix(rep(object$theta,nsim),nsim,length(object$theta),byrow=TRUE)
    } else {
      replicates=mvtnorm::rmvnorm(nsim,object$theta,C)
    }
    p_mu=NCOL(object$poly[['mu']])
    p_sigma=NCOL(object$poly[['sigma']])
    p_phi=NCOL(object$poly[['phi']])
    mus=sigmas=phis=predictive=matrix(NA,nsim,n)
    for(i in 1:nsim){
      if(p_mu>0){
        mus[i,]=object$link_mu(as.matrix(poly[['mu']]) %*% replicates[i,1:p_mu])
      } else {
        mus[i,]=rep(as.numeric(row.names(poly[['mu']])),n)
      }
      if(p_sigma>0){
        sigmas[i,]=object$link_sigma(as.matrix(poly[['sigma']]) %*% replicates[i,p_mu+(1:p_sigma)])
      } else {
        sigmas[i,]=rep(as.numeric(row.names(poly[['sigma']])),n)
      }
      if(p_phi>0){
        phis[i,]=object$link_phi(as.matrix(poly[['phi']]) %*% replicates[i,p_mu+p_sigma+(1:p_phi)])
      } else {
        phis[i,]=rep(as.numeric(row.names(poly[['phi']])),n)
      }
      innov=rnorm(n)
      foo=0*innov
      if(object$useApproach2){
        for(j in 2:n){
          ratio=sigmas[i,j]/sigmas[i,j-1]
          foo[j]=phis[i,j]*ratio*foo[j-1]+mus[i,j]-
            phis[i,j]*ratio*mus[i,j-1]+sigmas[i,j]*innov[j]
        }
      } else {
        for(j in 2:n){
          foo[j]=phis[i,j]*foo[j-1]+mus[i,j]+sigmas[i,j]*innov[j]
        }
      }
      predictive[i,]=rnorm(n,mean=mus[i,],sd=sigmas[i,])
    }
  }
  replicates=list(mu=mus,sigma=sigmas,phi=phis)
  out=list(pred=pred,replicates=replicates,predictive=predictive,
           x=newdat,xscaled=xscaled,poly=poly,C=C)
  class(out) <- 'polyRegPred'
  return(out)
}

#' Plot Polynomial Regression
#'
#' Plot results of a polynomial regression.
#' @param x object of class polyRegFit, resulting from a call to function [polynomialRegression()]
#' @param predictorType character string, one of 'scaled' (show rescaled predictor), 'raw' (show raw predictor)
#'     or 'transformed' (show transformed predictor, only available when x is of class anaRes).
#' @param alpha numeric in [0,1], transparency of uncertainty intervals.
#' @param ... further plotting arguments passed to or from other methods.
#' @return A list of 4 ggplots:
#' \enumerate{
#'   \item 'obs': predicted values (with uncertainties) and standardized residuals
#'   \item 'mu': for each mu predictor, a 3-panel plot with the predictor on the x-axis,
#'       and on the y-axis: estimated mu, predicted values and standardized residuals
#'   \item 'sigma': same as previously but for sigma
#'   \item 'phi': same as previously but for phi
#' }
#' @inherit polynomialRegression examples
#' @export
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom dplyr arrange
plot.polyRegFit <- function(x,predictorType=c('scaled','raw','transformed'),alpha=0.25,...){
  out=list()
  wx=match.arg(predictorType)
  whichx=switch(wx,scaled='xscaled',raw='x',transformed='xtrans')
  mu=x$pred[['mu']]
  sigma=x$pred[['sigma']]
  phi=x$pred[['phi']]
  DF=data.frame(y=x$y,x=1:length(x$y),mu=mu,sigma=sigma,phi=phi,
                res=x$residuals,innov=x$innovations)
  g1=ggplot(DF)+
    geom_point(aes(.data$x,.data$y))+
    geom_ribbon(aes(x=.data$x,ymin=.data$mu-1.64*.data$sigma,ymax=.data$mu+1.64*.data$sigma),
                fill='red',alpha=alpha)+
    geom_line(aes(.data$x,.data$mu),col='red')+
    coord_cartesian(ylim=range(DF$y,na.rm=TRUE))+
    labs(x='Index',y='Predictand y')+
    theme_bw()
  g2=ggplot(DF)+
    geom_point(aes(.data$x,.data$res))+
    geom_ribbon(aes(x=.data$x,ymin=-1.64,ymax=1.64),fill='red',alpha=alpha)+
    geom_line(aes(x,0),col='red')+
    coord_cartesian(ylim=max(abs(DF$res),na.rm=TRUE)*c(-1,1))+
    labs(x='Index',y='Standardized residuals')+
    theme_bw()
  g3=ggplot(DF)+
    geom_point(aes(.data$x,.data$innov))+
    geom_ribbon(aes(x=.data$x,ymin=-1.64,ymax=1.64),fill='red',alpha=alpha)+
    geom_line(aes(x,0),col='red')+
    coord_cartesian(ylim=max(abs(DF$innov),na.rm=TRUE)*c(-1,1))+
    labs(x='Index',y='Innovations')+
    theme_bw()
  out[['obs']]=patchwork::wrap_plots(g1,g2,g3,ncol=1)+plot_layout(axes='collect')
  for(what in c('mu','sigma','phi')){
    gs=vector('list',NCOL(x[[whichx]][[what]]))
    xnames=names(x[[whichx]][[what]])
    for(i in 1:NCOL(x[[whichx]][[what]])){
      DF=data.frame(y=x$y,x=x[[whichx]][[what]][,i],mu=mu,sigma=sigma,phi=phi,
                    pred=x$pred[[what]],res=x$residuals)
      DF=arrange(DF,x)
      g1=ggplot(DF)+
        geom_point(aes(.data$x,.data$y))+
        geom_ribbon(aes(x=.data$x,ymin=.data$mu-1.64*.data$sigma,ymax=.data$mu+1.64*.data$sigma),
                    fill='red',alpha=alpha)+
        geom_line(aes(.data$x,.data$mu),col='red')+
        coord_cartesian(ylim=range(DF$y,na.rm=TRUE))+
        labs(x=xnames[i],y='Predictand y')+
        theme_bw()
      g2=ggplot(DF)+
        geom_line(aes(.data$x,.data$pred),col='red')+
        labs(x=xnames[i],y=what)+
        theme_bw()
      if(what=='phi') {
        g2=g2+coord_cartesian(ylim=c(-1,1))+geom_hline(yintercept=0)
      } else {
        g2=g2+coord_cartesian(ylim=range(DF$pred,na.rm=TRUE))
      }
      g3=ggplot(DF)+
        geom_point(aes(.data$x,.data$res))+
        geom_ribbon(aes(x=.data$x,ymin=-1.64,ymax=1.64),fill='red',alpha=alpha)+
        geom_line(aes(x,0),col='red')+
        coord_cartesian(ylim=max(abs(DF$res),na.rm=TRUE)*c(-1,1))+
        labs(x=xnames[i],y='Standardized residuals')+
        theme_bw()
      gs[[i]]=patchwork::wrap_plots(g2,g1,g3,ncol=1)+plot_layout(axes='collect')
    }
    out[[what]]=gs
  }
  return(out)
}

#' Plot Polynomial Regression Prediction
#'
#' Plot prediction from a fitted polynomial regression.
#' @param x object of class polyRegPred, resulting from a call to function [predict.polyRegFit()]
#' @param predictorType character string, one of 'scaled' (show rescaled predictor), 'raw' (show raw predictor)
#'     or 'transformed' (show transformed predictor, only available when x is of class anaResPred).
#' @param ... further plotting arguments passed to or from other methods.
#' @return A list of 3 ggplots:
#' \enumerate{
#'   \item 'mu': for each mu predictor, a 2-panel plot with the predictor on the x-axis,
#'       and on the y-axis: estimated mu, predicted values
#'   \item 'sigma': same as previously but for sigma
#'   \item 'phi': same as previously but for phi
#' }
#' @inherit predict.polyRegFit examples
#' @export
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom dplyr arrange
#' @importFrom stats quantile
plot.polyRegPred <- function(x,predictorType=c('scaled','raw','transformed'),...){
  out=list()
  # Predictive
  foo=apply(x$predictive,2,quantile,probs=c(0.025,0.5,0.975))
  DF=data.frame(med=foo[2,],low=foo[1,],high=foo[3,],x=1:NCOL(foo))
  g=ggplot(DF)+
    geom_ribbon(aes(x=.data$x,ymin=.data$low,ymax=.data$high),fill='red',alpha=0.5)+
    geom_line(aes(.data$x,.data$med),col='red')+
    labs(x='Index',y='Predictand y')+
    theme_bw()
  out[['pred']]=g
  # Parameters
  wx=match.arg(predictorType)
  whichx=switch(wx,scaled='xscaled',raw='x',transformed='xtrans')
  for(what in c('mu','sigma','phi')){
    gs=vector('list',NCOL(x[[whichx]][[what]]))
    if(is.null(x$pred[[what]])){
      high=low=x$pred[[what]]
      pred_med=x$pred[['mu']]
      pred_low=pred_med-1.64*x$pred[['sigma']]
      pred_high=pred_med+1.64*x$pred[['sigma']]
    } else {
      foo=apply(x$replicates[[what]],2,quantile,probs=c(0.025,0.975))
      low=foo[1,];high=foo[2,]
      foo=apply(x$predictive,2,quantile,probs=c(0.025,0.5,0.975))
      pred_low=foo[1,];pred_med=foo[2,];pred_high=foo[3,]
    }
    DF=data.frame(par=x$pred[[what]],low=low,high=high,
                  pred=pred_med,pred_low=pred_low,pred_high=pred_high,x=NA)
    xnames=names(x[[whichx]][[what]])
    for(i in 1:NCOL(x[[whichx]][[what]])){
      DF$x=x[[whichx]][[what]][,i]
      DF=arrange(DF,x)
      g1=ggplot(DF)+
        geom_ribbon(aes(x=.data$x,ymin=.data$low,ymax=.data$high),fill='red',alpha=0.5)+
        geom_line(aes(.data$x,.data$par),col='red')+
        labs(x=xnames[i],y=what)+
        theme_bw()
      g2=ggplot(DF)+
        geom_ribbon(aes(x=.data$x,ymin=.data$pred_low,ymax=.data$pred_high),fill='red',alpha=0.5)+
        geom_line(aes(.data$x,.data$pred),col='red')+
        labs(x=xnames[i],y='Predictand y')+
        theme_bw()
      gs[[i]]=patchwork::wrap_plots(g1,g2,ncol=1)+plot_layout(axes='collect')
    }
    out[[what]]=gs
  }
  return(out)
}

#***************************************************************************----
# Private functions ----

# Function to assemble polynomial matrix
polyReg_getPolynomials <-function(xscaled,polyFunk,deg){
  if(is.character(deg)){
    # Fixed parameter value, return an empty data frame with the value stored in the row name
    poly=data.frame(row.names=deg)
  } else {
    p=NCOL(xscaled)
    n=NROW(xscaled)
    # Assemble polynomials
    poly=data.frame(P0=rep(1,n))
    for(i in 1:p){
      if(deg[i]>0){
        foo=polyFunk(xscaled[,i],deg[i])[,-1]
        names(foo) <- paste0('P',1:deg[i],'_X',i)
        poly=cbind(poly,foo)
      }
    }
  }
  return(poly)
}

# Negative log-likelihood to be minimized
polyReg_nll <- function(theta,y,poly_mu,poly_sigma,poly_phi,link_mu,link_sigma,link_phi,useApproach2){
  p_mu=NCOL(poly_mu)
  p_sigma=NCOL(poly_sigma)
  p_phi=NCOL(poly_phi)
  n=length(y)
  if(length(theta) != p_mu+p_sigma+p_phi){
    stop('size of theta does not match that of poly_mu and/or poly_sigma and/or poly_phi',call.=FALSE)
  }
  if(p_mu>0){
    mu=link_mu(as.matrix(poly_mu) %*% theta[1:p_mu])
  } else {
    mu=rep(as.numeric(row.names(poly_mu)),n)
  }
  if(p_sigma>0){
    sigma=link_sigma(as.matrix(poly_sigma) %*% theta[p_mu+(1:p_sigma)])
  } else {
    sigma=rep(as.numeric(row.names(poly_sigma)),n)
  }
  if(p_phi>0){
    phi=link_phi(as.matrix(poly_phi) %*% theta[p_mu+p_sigma+(1:p_phi)])
  } else {
    phi=rep(as.numeric(row.names(poly_phi)),n)
  }
  # Remove 1st data since conditioning on the previous one is not possible
  ix=2:n
  yt=y[ix]
  ytm1=y[ix-1]
  mut=mu[ix]
  mutm1=mu[ix-1]
  sigmat=sigma[ix]
  sigmatm1=sigma[ix-1]
  phit=phi[ix]
  # Get mean depending on the requested approach
  if(useApproach2){
    ratio=sigmat/sigmatm1
    moy=mut-phit*ratio*mutm1+phit*ratio*ytm1
  } else {
    moy=mut+phit*ytm1
  }
  if(any(sigmat<=0)){return(Inf)}
  terms=stats::dnorm(yt,mean=moy,sd=sigmat,log=TRUE)
  return(-1*sum(terms,na.rm=TRUE))
}

