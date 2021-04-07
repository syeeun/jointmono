
fit_ramsay = function(tobs, yobs, ini.theta, trng = NULL, norder = 5,
                      lmu, elim, maxiter){
  
  if(is.null(trng)){ trng = range(tobs) }
  rng = trng 
  nobs = length(tobs)
  
  tfine = seq(rng[1], rng[2], l = 500)
  
  #### initialize and set storage to track updates####
  theta = theta.storage = ini.theta

  nbasis = length(ini.theta);
  
  bsp = o.bsp.mat(tfine, nbasis = nbasis, norder = norder)
  Bsp = o.bsp.mat(tfine, nbasis = nbasis, norder = norder, integ = T)
  bspmat = bsp$obspmat; Bspmat = Bsp$obspmat
  penmat = diag(nbasis) 
  
  b0.storage = b1.storage = NULL
  bcoef = c(NA, 1)
  
  re = 0
  repeat{
    #### temporal space ####
    XX_mu = matrix(0, nbasis, nbasis); Xr_mu = matrix(0, nbasis, 1)
    
    #### update theta_mu ####
      wcoef = matrix(theta, nbasis, 1); 
      Hmat = Hfn(tobs, wcoef, Bsp)
      bcoef = lm(yobs ~ Hmat)$coef
      
      res = yobs-bcoef[1]-bcoef[2]*Hmat
      dH = dHfn(tobs, wcoef, Bsp)
      d2H = crossprod(dH)#t(dH) %*% dH; 
      dHy = crossprod(dH, res)
      
      XX_mu = bcoef[2]^2 * d2H / nobs; Xr_mu = bcoef[2] * dHy / nobs
      theta_new = theta + ginv(XX_mu + lmu*penmat) %*%  (Xr_mu - lmu*penmat%*%theta)
  
    eps.t = abs(theta_new-theta)/(1+abs(theta))

    #### store updates ####
    theta = theta_new; 
    
    theta.storage = cbind(theta.storage, theta)
    b0.storage = cbind(b0.storage, bcoef[1]); 
    b1.storage = cbind(b1.storage, bcoef[2])
    
    #### iterate or stop ####
    eps = max(eps.t)
    re = re + 1; #print(re)
    if(eps < elim | re > (maxiter-1)) break
  }
  
  error = sqrt(sum((yobs - bcoef[1] - bcoef[2]*Hmat)^2)/nobs)
  
  #print(paste(re, "iterated"))
  result = list(dat = data.frame(tobs, yobs),
                basis = list(bsp = bsp, Bsp = Bsp),
                est = list(theta = theta, bcoef = bcoef, sd = error),
                track = list(theta.storage = theta.storage, b0.storage = b0.storage, b1.storage = b1.storage))
  
  return(result)
}


std_fitted_ramsay = function(obj){
  
  theta = obj$est$theta; bcoef = obj$est$bcoef
  yobs = obj$dat$yobs; tobs = obj$dat$tobs; nobs = length(tobs)

  nbasis = length(theta);
  Bsp = obj$basis$Bsp;  #bsp = obj$basis$bsp; 
  Bspmat = Bsp$obspmat #bspmat = bsp$obspmat;
  
  wcoef = matrix(theta, nbasis, 1); 
  Hmat = Hfn(tobs, wcoef, Bsp)
  yhat = bcoef[1] + bcoef[2]*Hmat
  
  result = sqrt(sum((yobs - yhat)^2)/nobs)#/(nobs-2))
  return(result)
}



plot_fitted_ramsay = function(obj, curvetype = c("y", "mu")){
  
  theta = obj$est$theta; bcoef = obj$est$bcoef
  tfine = obj$basis$bsp$tgrid
  
  wcoef = matrix(theta, length(theta), 1); 
  Wmat = Wfn(tfine, wcoef, obj$basis$bsp); 
  Hmat = Hfn(tfine, wcoef, obj$basis$Bsp)
  yhat = bcoef[1] + bcoef[2]*Hmat
  
  if(curvetype == "y"){
    plot(tfine, yhat, type="l", ylab = "y(t)")#, ylim = range(obj$dat$yobs))
    with(obj$dat, points(tobs, yobs, pch = 16, cex = .4))
  }
  
  if(curvetype == "mu"){
    plot(tfine, Wmat, type="l", ylab=expression(paste(mu, "(t)")));abline(h=0)
  }
}


track_fitted_ramsay = function(obj, trackpar = c("mu", "beta")){
  
  re = ncol(obj$track$b0.storage); 
  
  if(trackpar == "mu"){ 
    par(mfrow=c(1,1))
    TRACK = obj$track$theta.storage; 
    matplot(0:re, t(TRACK), type="l", ylab = expression(theta[mu])) 
    }

  if(trackpar == "beta"){ 
    par(mfrow=c(1,2))
    TRACK = obj$track$b0.storage; 
    matplot(1:re, t(TRACK), type="l", ylab = expression(beta[0]))
    TRACK = obj$track$b1.storage; 
    matplot(1:re, t(TRACK), type="l", ylab = expression(beta[1]))
    }
  
  par(mfrow=c(1,1))
}



