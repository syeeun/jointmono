
ini_james = function(dat, nbasis, norder, npc){
  
  tfine = with(dat, tgrid[id==1]) # same t-grid for all curves;
  rng = range(tfine); ncurve = max(dat$id)
  
  B = o.bsp.mat(tfine, nbasis = nbasis, norder = norder)$obspmat
  
  BY0 = 0
  for(i in 1:ncurve){
    whatitmp = with(dat, what[id==i])
    BY0 = BY0 + t(B[!is.na(whatitmp),]) %*% whatitmp[!is.na(whatitmp)]
  }
  theta0 = solve(t(B) %*% B * ncurve) %*% BY0
  
  wfd.res = with(dat, what - rep(B %*% theta0, ncurve))
  Gam = matrix(0, nbasis, ncurve)
  for(i in 1:ncurve){
    wresitmp = wfd.res[dat$id==i]
    Gam[,i] = lm(wresitmp[!is.na(wresitmp)] ~ B[!is.na(wresitmp),]-1)$coef
  }
  
  Gam = ifelse(is.na(Gam), 0, Gam)
 
  Theta0 = with(svd(t(Gam)), v %*% diag(d)[, 1:npc])
  Theta0 = qr.Q(qr(Theta0))
  Tsign0 = sign(Theta0[1,]); Theta0 = t(Tsign0*t(Theta0)) 
  
  ascore0 = NULL
  for(i in 1:ncurve){
    ascore0 = rbind(ascore0, lm(Gam[,i] ~ Theta0 - 1)$coef)
  }
  
  const = apply(ascore0, 2, mean)
  ascore0 = ascore0 - matrix(const, nrow(ascore0), ncol(ascore0), byrow=T) #centering
  aorder = order(apply(ascore0, 2, var), decreasing = T); 
  ascore0 = ascore0[, aorder] ; Theta0 = Theta0[, aorder] #reordering
  # theta0 = theta0+Theta0*const ###
  
  ini = list(theta = theta0, Theta = Theta0, ascore = ascore0)
  
  return(ini)
}







fit_james = function(dat, ini, nbasis, norder, npc, elim = 0.05, maxiter = 30){
  
  id = dat$id
  w = dat$what
  y = dat$yhat
  tfine = dat$tgrid
  
  ncurve = max(id) #number of curves
  
  theta.up = matrix(0, nbasis, 1); 
  Theta.up = matrix(0, nbasis, npc);  
  ascore.up = matrix(0, ncurve, npc)
  
  theta = matrix(ini$theta, nbasis, 1); 
  Theta = matrix(ini$Theta, nbasis, npc); 
  ascore = matrix(ini$ascore, ncurve, npc)
  theta.storage = theta; Theta.storage = ascore.storage = list()
  for(k in 1:npc){
    Theta.storage[[k]] = Theta[,k]; ascore.storage[[k]] = ascore[,k]
  }
  
  bsp = o.bsp.mat(tfine[id==1], nbasis, norder)
  Bsp = o.bsp.mat(tfine[id==1], nbasis, norder, integ = T)
  B = bsp$obspmat; BB = t(B)%*%B; 
  
  re = 0
  repeat{
    
    # theta
    BY = 0
    for(i in 1:ncurve){
      whatitmp = with(dat, what[id==i])
      BY = BY + t(B[!is.na(whatitmp),]) %*% (w[id==i] - B %*% Theta %*% ascore[i, ])[!is.na(whatitmp)]
    }
    theta.up = solve(ncurve*BB, BY)
    
    for(k in 1:npc){
      aBY = 0; 
      for(i in 1:ncurve){  
        whatitmp = with(dat, what[id==i])
        aBY = aBY + ascore[i,k] * t(B[!is.na(whatitmp),]) %*% (w[id==i] - B%*%theta.up - B%*%Theta[,-k]%*%ascore[i,-k])[!is.na(whatitmp),]
      }
      Theta.up[,k] = solve(sum(ascore[,k]^2)*BB, aBY) #) %*% aBY
    }
    Theta.up = qr.Q(qr(Theta.up))
    #### avoid sign flipped ####
    for(k in 1:npc){
      flip_sum = as.matrix(Theta.up[,k] + Theta[,k]); 
      flip_diff = as.matrix(Theta.up[,k] - Theta[,k])
      if(norm(flip_sum,"f") < norm(flip_diff,"f")) { Theta.up[,k] = -Theta.up[,k] }
    }
    
    for (i in 1:ncurve){
      whatitmp = with(dat, what[id==i])
      ascore.up[i, ] = solve(t(Theta.up)%*%BB%*%Theta.up, t(B[!is.na(whatitmp),]%*%Theta.up)%*%((w[id==i] - B%*%theta.up)[!is.na(whatitmp),]))
    }
    const = apply(ascore.up, 2, mean) ###
    ascore.up = ascore.up - matrix(const, ncurve, npc, byrow=T) #center a
    aorder = order(apply(ascore.up, 2, var), decreasing = T)
    ascore.up = matrix(ascore.up[, aorder], ncurve, npc); 
    Theta.up = matrix(Theta.up[, aorder], nbasis, npc); 
    #theta.up = theta.up+Theta.up*const ###
    
    eps.t = max(abs(theta.up-theta)/(abs(theta)))
    eps.T = max(abs(Theta.up-Theta)/(abs(Theta)))
    eps.a = max(abs(ascore.up-ascore)/(abs(ascore)))
    eps = max(c(eps.t, eps.T, eps.a))
    
    theta.storage = cbind(theta.storage, theta.up)
    for(k in 1:npc){
      Theta.storage[[k]] = cbind(Theta.storage[[k]], Theta.up[,k])
      ascore.storage[[k]] = cbind(ascore.storage[[k]], ascore.up[,k])
    }
    theta = theta.up; Theta = Theta.up; ascore = ascore.up
    
    re = re+1; #print(re)
    if(eps < elim | re > maxiter -1) break # stop when all of parameters changes less than 10%
  }
  
  if(eps > elim) {print("not converge")}
  #else {print(paste(re, "iterated"))}
  # 
  # w_jam = h_jam = y_jam = bcoef = NULL
  # if(eps < .10) {print("converged")
  #   for(i in 1:ncurve){
  #     wcoef_i = theta + Theta%*%ascore[i,]
  #     w_jam_i = Wfn(tfine[id==i], wcoef_i, bsp)
  #     h_jam_i = Hfn(tfine[id==i], wcoef_i, Bsp)
  #     bcoef_i = lm(y[id==i] ~ h_jam_i)$coef
  #     y_jam_i = bcoef_i[1] + bcoef_i[2]*h_jam_i
  #     
  #     y_jam = c(y_jam, y_jam_i);
  #     h_jam = c(h_jam, h_jam_i); 
  #     w_jam = c(w_jam, w_jam_i); 
  #     bcoef = rbind(bcoef, bcoef_i); 
  #   }
  # }
  # 
  
  colnames(Theta) = paste("fpc", 1:npc, sep="")
  rownames(Theta) = rownames(theta) #paste("bspl", norder, ".",1:nbasis, sep = "")
  colnames(ascore) = paste("ascore", 1:npc, sep="")
  
  result = list(dat = dat, 
                basis = list(bsp = bsp, Bsp = Bsp),
                est = list(theta = theta, Theta = Theta, ascore = ascore), 
                track = list(theta.storage = theta.storage, Theta.storage = Theta.storage, ascore.storage = ascore.storage))
  return(result)
}