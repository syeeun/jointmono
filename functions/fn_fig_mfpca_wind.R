

plot_fitted_mfpca = function(obj, curvetype = c("y", "w", "fpc", "mu"), ptit = F, ...){
  
  theta = obj$est$theta; 
  Theta = as.matrix(obj$est$Theta); 
  ascore = as.matrix(obj$est$ascore); 
  bcoef = obj$est$bcoef
  
  ncurve = max(obj$dat$id);
  
  tfine = obj$basis$bsp$tgrid
  nbasis = nrow(theta); npc = ncol(Theta)
  bsp = obj$basis$bsp; Bsp = obj$basis$Bsp
  bspmat = bsp$obspmat; Bspmat = Bsp$obspmat
  
  W0 = What = Hhat = yhat = NULL
  for(i in 1:ncurve){
    wcoef = theta + Theta %*% ascore[i,];
    #wcoef0 = with(true.par, theta + Theta %*% ascore[i,])
    Wmat = Wfn(tfine, wcoef, bsp); 
    #W0mat = Wfn(tfine, wcoef0, bsp); 
    Hmat = Hfn(tfine, wcoef, Bsp)
    
    What = cbind(What, Wmat); 
    #W0 = cbind(W0, W0mat); 
    Hhat = cbind(Hhat, Hmat)
    yhat = cbind(yhat, bcoef[i,1] + bcoef[i,2]*Hmat)
  }
  
  if(curvetype == "w"){
    #### plot estimated w(t) ####
    matplot(tfine, What, type="l", col = 1:ncurve, lty = 1, ...)#, W0)))
    #matlines(tfine, W0, col = 1:ncurve, lty = 2)
    legend("topright", c("estimated", "true"), lty = 1:2, bty = "n")
  }
  
  if(curvetype == "y"){
    #### plot estimated y(t) ####
    matplot(tfine, yhat, type="l", col = 1:ncurve, lty = 1, ...)
    if(ptit ==T ) {with(obj$dat, points(tobs, yobs, col=id, pch = 16, cex = .4))}
    #legend("topleft", c("observed", "estimated"), lty = c(NA, 1), pch = c(16, NA), bty = "n")
  }
  
  if(curvetype == "fpc"){
    plot(range(tfine), range(c(bspmat%*%Theta)), type="n", ...); abline(h=0, col="gray")
    for(k in 1:npc){
      points(tfine, bspmat%*%Theta[,k], type="l", col=k, lty = k)
      #points(tfine, bspmat%*%true.par$Theta[,k], type="l", lty=2, col=k);
    }
    #legend("topright", c(paste("estimated f", 1:npc, sep = ""), 
    #                     paste("true f", 1:npc, sep = "")), 
    #       col=1:npc, lty=c(rep(1,npc),rep(2,npc)), bty="n")
  }
  
  if(curvetype == "mu"){
    plot(tfine, bspmat%*%theta, type="l", ...);
    #lines(tfine, bspmat%*%true.par$theta, lty=2);abline(h=0)
    #legend("topright", c("estimated", "true"), lty = 1:2, bty = "n")
  }
}



#### plot result of fit_shin ####
track_fitted_mfpca = function(obj, trackpar = c("mu", "fpc", "alpha", "beta", "std")){
  
  re = ncol(obj$track$theta.storage)-1; 
  npc = ncol(obj$est$Theta)
  
  par(mfrow=c(1,1))
  if(trackpar == "mu"){ TRACK = obj$track$theta.storage; 
  matplot(0:re, t(TRACK), type="l", ylab=trackpar) }
  if(trackpar == "fpc"){ par(mfrow=c(1,npc))
    for(k in 1:npc){
      TRACK = obj$track$Theta.storage[[k]];
      matplot(0:re, t(TRACK), type="l", ylab=paste(trackpar, k, sep=""))} }
  if(trackpar == "alpha"){ par(mfrow=c(1,npc))
    for(k in 1:npc){
      TRACK = obj$track$a.storage[[k]]
      matplot(0:re, t(TRACK), type="l", ylab=paste(trackpar, k, sep=""), ylim = range(obj$track$a.storage))} }
  if(trackpar == "beta"){ par(mfrow=c(1,2))
    TRACK = obj$track$b0.storage; 
    matplot(0:re, t(TRACK), type="l", ylab="beta0")
    TRACK = obj$track$b1.storage; 
    matplot(0:re, t(TRACK), type="l", ylab="beta1")}
  if(trackpar == "std"){ TRACK = obj$track$s.storage
  plot(1:re, TRACK, type="l", ylab=trackpar) }
  par(mfrow=c(1,1))
}

plot.yfit = function(obj, idx, rnginput = NULL, ...){
  
  wcoef = with(obj$est, theta + Theta %*% ascore[idx,]);
  Hmat = with(obj$basis, Hfn(bsp$tgrid, wcoef, Bsp))
  yhat = with(obj$est, bcoef[idx,1] + bcoef[idx,2]*Hmat)
  
  #### plot estimated y(t) ####
  with(obj$dat, plot(tobs[id==idx], yobs[id==idx], 
                     xlim = range(obj$dat$tobs), #range(obj$basis$bsp$tgrid),
                     ylim = range(c(yobs[id==idx], rnginput)), 
                     pch = 16, cex = .4, xlab = "t", ylab = paste("power curve", idx), ...))
  points(obj$basis$bsp$tgrid, yhat, type="l", lwd = 2, ...)
  abline(h=c(0,1), lty = 3)
  
  #legend("topleft", c("observed", "estimated"), lty = c(NA, 1), pch = c(16, NA), bty = "n")
}


plot.muf = function(obj, K, hlim, wlim, ...){
  
  tlim = range(obj$dat$tobs)
  
  par(mfrow = c(K, 2))
  for(k in 1:K){
    with(obj, plot(basis$Bsp$tgrid, Hfn(basis$Bsp$tgrid, est$theta, basis$Bsp), 
                   type = "l", lty = 1, xlim = tlim, ylim = hlim, xlab = "t", ylab = "H(t)",...));
    with(obj, points(basis$Bsp$tgrid, Hfn(basis$Bsp$tgrid, est$theta + est$Theta[,k], basis$Bsp), type = "l", lty = 2)); 
    with(obj, points(basis$Bsp$tgrid, Hfn(basis$Bsp$tgrid, est$theta - est$Theta[,k], basis$Bsp), type = "l", lty = 2)); 
    
    with(obj, plot(basis$Bsp$tgrid, Wfn(basis$Bsp$tgrid, est$theta, basis$bsp), 
                   type = "l", lty = 1, xlim = tlim, ylim = wlim, xlab = "t", ylab = "w(t)", ...));
    with(obj, points(basis$Bsp$tgrid, Wfn(basis$Bsp$tgrid, est$theta + est$Theta[,k], basis$bsp), type = "l", lty = 2)); 
    with(obj, points(basis$Bsp$tgrid, Wfn(basis$Bsp$tgrid, est$theta - est$Theta[,k], basis$bsp), type = "l", lty = 2)); 
    abline(h = 0, lty = 3)
  }
  par(mfrow = c(1,1))
}


plot.ym = function(obj, idx, trng, yrnginput = NULL, ...){
  
  tfine = seq(trng[1], trng[2], l = 100)
  
  wcoef = with(obj$est, theta + Theta %*% ascore[idx,]);
  Hmat = with(obj$basis, Hfn(tfine, wcoef, Bsp))
  yhat = with(obj$est, bcoef[idx,1] + bcoef[idx,2]*Hmat)
  
  #### plot estimated y(t) ####
  with(obj$dat, plot(tobs[id==idx], yobs[id==idx], 
                     xlim = range(obj$dat$tobs), #range(obj$basis$bsp$tgrid),
                     ylim = range(c(yobs[id==idx], yrnginput)), 
                     pch = 16, cex = .3, col = "grey70", yaxt = "n", ...))
  axis(2, at = 0:5/5, paste(0:5*20, "%"))
  points(tfine, yhat, type="l", lty = 2, lwd = 2, ...)
  abline(h=c(0,1), lty = 3)
  
  points(tfine, obj$est$bcoef[idx,1] + obj$est$bcoef[idx,2] * Hfn(tfine, obj$est$theta, oBsp), 
         type="l", lty = 1, lwd = 1.5, col = "grey40")
  
  #legend("topleft", c("observed", "estimated"), lty = c(NA, 1), pch = c(16, NA), bty = "n")
}
