fn_data_simul = function(ncurve, Sig = c(1, .5)^2, sig = .05, seedn = 1234,
                         nbasis = 10, norder = 5, ntime = 1000){
  
  #### true curves ####
  mu = function(t){ 5*(1-2*t)}
  f1 = function(t){sqrt(2)*sin(2*pi*(t))}
  f2 = function(t){sqrt(2)*cos(2*pi*(t))}
  
  #### basis ####
  tgrid = seq(0, 1, l = ntime)
  obsp = o.bsp.mat(tgrid, nbasis = nbasis, norder = norder)
  oBsp = o.bsp.mat(tgrid, nbasis = nbasis, norder = norder, integ = T)
  bsp_sim = obsp$obspmat
  
  #### true parameters ####
  theta_mu = lm(mu(tgrid) ~ bsp_sim - 1)$coef
  theta_f1 = lm(f1(tgrid) ~ bsp_sim - 1)$coef
  theta_f2 = lm(f2(tgrid) ~ bsp_sim - 1)$coef
  theta_f = cbind(theta_f1, theta_f2)
  
  npc = length(Sig)
  set.seed(1)
  alpha = mvrnorm(ncurve, mu = rep(0, npc), Sigma = diag(Sig^2, npc, npc)); 
  alpha = alpha - matrix(apply(alpha, 2, mean), ncurve, npc, byrow=T)
  
  Tsign = sign(theta_f[1,]); theta_f = t(Tsign*t(theta_f)) ; alpha = t(Tsign*t(alpha))
  
  sigma = sig
  
  true.par = list(theta = matrix(theta_mu, nbasis, 1), 
                  Theta = theta_f[, 1:npc, drop = F], 
                  ascore = alpha[, 1:npc, drop = F], 
                  bcoef = cbind(rep(0, ncurve), rep(1, ncurve)),
                  sigma = sigma)
  
  #### simulate observations ####
  set.seed(seedn)
  n = sample(50:100, ncurve, replace = T); 
  id = rep(1:ncurve, times = n)
  tobs = wobs = hobs = NULL
  for (m in 1:ncurve) { 
    set.seed(seedn + m)
    tsample = sort(sample(0:ntime, n[m])/ntime)
    wcoef = theta_mu + theta_f[, 1:npc, drop = F]%*%t(alpha[m, 1:npc, drop = F])
    tobs = c(tobs, tsample) 
    wobs = c(wobs, Wfn(tsample, wcoef, obsp))
    hobs = c(hobs, Hfn(tsample, wcoef, oBsp)) 
  }
  
  set.seed(seedn)
  yobs = hobs + rnorm(length(hobs), 0, sigma)
  dat.sim = data.frame(id, tobs, wobs, hobs, yobs)
  
  return(list(dat = dat.sim, par = true.par))
}