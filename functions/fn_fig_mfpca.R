
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
      matplot(0:re, t(TRACK), type="l", ylab=paste(trackpar, k, sep=""))} }
  if(trackpar == "beta"){ par(mfrow=c(1,2))
    TRACK = obj$track$b0.storage; 
    matplot(0:re, t(TRACK), type="l", ylab="beta0")
    TRACK = obj$track$b1.storage; 
    matplot(0:re, t(TRACK), type="l", ylab="beta1")}
  if(trackpar == "std"){ TRACK = obj$track$s.storage
  plot(1:re, TRACK, type="l", ylab=trackpar) }
  par(mfrow=c(1,1))
}


plot_fitted_mfpca = function(obj, curvetype = c("y", "w", "fpc", "mu"), true.par = NULL){
  
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
    wcoef0 = with(true.par, theta + Theta %*% ascore[i,])
    Wmat = Wfn(tfine, wcoef, bsp); 
    W0mat = Wfn(tfine, wcoef0, bsp); 
    Hmat = Hfn(tfine, wcoef, Bsp)
    
    What = cbind(What, Wmat); 
    W0 = cbind(W0, W0mat); 
    Hhat = cbind(Hhat, Hmat)
    yhat = cbind(yhat, bcoef[i,1] + bcoef[i,2]*Hmat)
  }

  if(curvetype == "w"){
    #### plot estimated w(t) ####
    matplot(tfine, What, type="l", col = 1:ncurve, lty = 1, ylab = "w(t)",
            ylim = range(c(What, W0)))
    matlines(tfine, W0, col = 1:ncurve, lty = 2)
    legend("topright", c("estimated", "true"), lty = 1:2, bty = "n")
  }
  
  if(curvetype == "y"){
    #### plot estimated y(t) ####
    matplot(tfine, yhat, type="l", col = 1:ncurve, 
            ylim = range(c(obj$dat$yobs, yhat)), lty = 1)
    with(obj$dat, points(tobs, yobs, col=id, pch = 16, cex = .4))
    legend("topleft", c("observed", "estimated"), lty = c(NA, 1), pch = c(16, NA), bty = "n")
  }
  
  if(curvetype == "fpc"){
    plot(range(tfine), range(c(bspmat%*%Theta, bspmat%*%true.par$Theta)), type="n", 
         xlab = "tfine", ylab=curvetype); abline(h=0, col="gray")
    for(k in 1:npc){
      points(tfine, bspmat%*%Theta[,k], type="l", col=k)
      points(tfine, bspmat%*%true.par$Theta[,k], type="l", lty=2, col=k);
    }
    legend("topright", c(paste("estimated f", 1:npc, sep = ""), 
                         paste("true f", 1:npc, sep = "")), 
           col=1:npc, lty=c(rep(1,npc),rep(2,npc)), bty="n")
  }
  
  if(curvetype == "mu"){
    plot(tfine, bspmat%*%theta, type="l", 
         ylim = range(c(bspmat%*%theta, bspmat%*%true.par$theta)),
         ylab=expression(paste(mu, "(t)")));
    lines(tfine, bspmat%*%true.par$theta, lty=2);abline(h=0)
    legend("topright", c("estimated", "true"), lty = 1:2, bty = "n")
  }
}



plot_est_true_ini = function(obj, truepar, inipar, npc = 2){
  
  t = obj$basis$bsp$tgrid
  
  par(mfrow = c(2, 2))
  plot(t, Wfn(t, truepar$theta, obj$basis$bsp), ylab = expression(paste(mu,"(t)",sep="")), type="l", lty=2)
  lines(t, Wfn(t, obj$est$theta, obj$basis$bsp), col="red")
  lines(t, Wfn(t, inipar$theta, obj$basis$bsp), col = "blue")

  for(k in 1:npc){
  plot(t, Wfn(t, truepar$Theta[,k], obj$basis$bsp), ylab = paste("f",k,"(t)",sep=""), type="l", lty=2)
  lines(t, Wfn(t, obj$est$Theta[,k], obj$basis$bsp), col="red")
  lines(t, Wfn(t, inipar$Theta[,k], obj$basis$bsp), col = "blue")
  }

  plot(t, Wfn(t, truepar$theta, obj$basis$Bsp), ylab = "int w(t)", type="l", lty=2)
  lines(t, Wfn(t, obj$est$theta, obj$basis$Bsp), col="red")
  lines(t, Wfn(t, inipar$theta, obj$basis$Bsp), col = "blue")

  plot(t, Hfn(t, truepar$theta, obj$basis$Bsp), ylab = "h(t)", type="l", lty=2)
  lines(t, Hfn(t, obj$est$theta, obj$basis$Bsp), col="red")
  lines(t, Hfn(t, inipar$theta, obj$basis$Bsp), col = "blue")
  par(mfrow = c(1,1))
}


