
fit_mfpca = function(dat.obs, ini.par, trng = c(0,1), norder = 5,
                     lmu = 0, lpc = 0, elim, maxiter, printit = F){ 
  
  ncurve = max(dat.obs$id)
  if(is.null(trng)) {trng = range(dat.obs$tobs)}
  #rng = c(trng[1]-diff(trng)*.2, trng[2]+diff(trng)*.2)
  rng = trng
  tfine = seq(rng[1], rng[2], l = 1000)
  
  #### initialize and set storage to track updates####
  theta = theta.storage = ini.par$theta; 
  Theta = as.matrix(ini.par$Theta); 
  ascore = as.matrix(ini.par$ascore); 
  bcoef = ini.par$bcoef
  
  nbasis = nrow(theta); npc = ncol(Theta)
  
  Theta.storage = a.storage = list()
  for(k in 1:npc){
    Theta.storage[[k]] = as.matrix(Theta[,k])
    a.storage[[k]] = as.matrix(ascore[,k])
  }
  
  b0.storage = bcoef[, 1]; b1.storage = bcoef[, 2]
  s.storage = ini.par$sigma
  
  bsp = o.bsp.mat(tfine, nbasis = nbasis, norder = norder)
  Bsp = o.bsp.mat(tfine, nbasis = nbasis, norder = norder, integ = T)
  bspmat = bsp$obspmat; Bspmat = Bsp$obspmat
  
  re = 0
  repeat{
    #### temporal space ####
    XX_mu = matrix(0, nbasis, nbasis); Xr_mu = matrix(0, nbasis, 1)
    XX_f = array(0, dim = c(nbasis, nbasis, npc)); Xr_f = array(0, dim = c(nbasis, 1, npc))
    
    Theta_new = matrix(0, nbasis, npc)
    ascore_new = matrix(0, ncurve, npc)
    bcoef_new = bcoef = matrix(c(0,1), ncurve, 2, byrow= T)
    
    #### update theta_mu ####
    for(i in 1:ncurve){
      
      wcoef = theta+Theta%*%ascore[i,]; Hmat = with(dat.obs, Hfn(tobs[id==i], wcoef, Bsp))
      
      res = with(dat.obs, yobs[id==i]-bcoef[i,1]-bcoef[i,2]*Hmat) 
      dH = with(dat.obs, dHfn(tobs[id==i], wcoef, Bsp))
      d2H = t(dH) %*% dH; dHy = t(dH) %*% res
      
      XX_mu = XX_mu + bcoef[i,2]^2 * d2H; 
      Xr_mu = Xr_mu + bcoef[i,2] * dHy
      
      for(k in 1:npc){
        XX_f[,,k] = XX_f[,,k] + bcoef[i,2]^2 * ascore[i,k]^2 * d2H
        Xr_f[,,k] = Xr_f[,,k] + bcoef[i,2] * ascore[i,k] * dHy
      }
      
      #### update ascore ####
      ni = length(res)
      ascore_new[i,] = ascore[i,] + solve(bcoef[i,2] *t(Theta) %*% (d2H) %*% Theta, t(Theta) %*% dHy)
      
      #### update bcoef ####
      bcoef_new[i,] = with(dat.obs, lm(yobs[id==i] ~ Hmat))$coef
    }
    
    #### update theta_mu ####
    theta_new = theta + solve(XX_mu + lmu*diag(nbasis), Xr_mu - lmu*theta)
  
    #### update Theta_f ####
    for(k in 1:npc){
      Theta_new[,k] = Theta[,k] + solve(XX_f[,,k] + lpc*diag(nbasis), Xr_f[,,k])# - lpc*Theta[,k])
    }
    Theta_new = qr.Q(qr(Theta_new))
    #### avoid sign flipped ####
    for(k in 1:npc){
      flip_sum = as.matrix(Theta_new[,k] + Theta[,k]); 
      flip_diff = as.matrix(Theta_new[,k] - Theta[,k])
      if(norm(flip_sum,"f") < norm(flip_diff,"f")) { Theta_new[,k] = -Theta_new[,k] }
    }
    
    #### update ascore ####
    tempa = apply(ascore_new, 2, mean)
    ascore_new = ascore_new - matrix(tempa, ncurve, npc, byrow=T)
    #### avoid order of variance flipped ####
    aorder = order(apply(ascore_new, 2, var), decreasing = T);
    ascore_new = matrix(ascore_new[, aorder], ncurve, npc);
    Theta_new = matrix(Theta_new[, aorder], nbasis, npc);

    #### update sd ####
    ressq = 0
    for(i in 1:ncurve){
      wcoef = theta_new + Theta_new %*% ascore_new[i,]
      Hmat = Hfn(tfine, wcoef, Bsp)
      Hval = with(dat.obs, approx(tfine, Hmat, tobs[id==i])$y)
      ressq = ressq + with(dat.obs, sum((yobs[id==i] - bcoef[i,1] - bcoef[i,2]*Hval)^2))
    }
    sigma = sqrt(ressq/nrow(dat.obs)); 
    
    #### measure change ####
    eps.t = norm(theta_new-theta, "f")/(1+norm(theta, "f"))
    eps.T = norm(as.matrix(Theta_new-Theta), "f")/(1+norm(as.matrix(Theta), "f"))
    eps.a = norm(as.matrix(ascore_new-ascore), "f")/(1+norm(as.matrix(ascore), "f"))
    #eps.b = norm(bcoef_new-bcoef, "f")/(1+norm(bcoef, "f"))
    
    #### store updates ####
    theta = theta_new; Theta = Theta_new; ascore = ascore_new; bcoef = bcoef_new
    
    theta.storage = cbind(theta.storage, theta)
    for(k in 1:npc){
      Theta.storage[[k]] = cbind(Theta.storage[[k]], Theta[,k])
      a.storage[[k]] = cbind(a.storage[[k]], ascore[,k])
    }
    b0.storage = cbind(b0.storage, bcoef[,1]); b1.storage = cbind(b1.storage, bcoef[,2])
    s.storage = c(s.storage, sigma)
    
    #### iterate or stop ####
    eps = max(c(eps.t, eps.T, eps.a))#, eps.b))
    re = re + 1; if(printit == T) {print(re)}
    if(eps < elim | re > (maxiter-1)) break
  }
  
  colnames(Theta) = paste("fpc", 1:npc, sep="")
  rownames(Theta) = rownames(theta) #paste("bspl", norder, ".",1:nbasis, sep = "")
  colnames(ascore) = paste("ascore", 1:npc, sep="")
  colnames(bcoef) = paste("beta", 0:1, sep="")
  
  result = list(dat = dat.obs,
                est = list(theta=theta, Theta=Theta, ascore=ascore, bcoef=bcoef, sigma=sigma),
                track = list(theta.storage=theta.storage, Theta.storage=Theta.storage, a.storage=a.storage, 
                             b0.storage=b0.storage, b1.storage=b1.storage, s.storage=s.storage),
                basis = list(bsp = bsp, Bsp = Bsp))
  return(result)
}




