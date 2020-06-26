source("fn_ortho.R")
source("fn_ftnlz.R")
source("fn_data_simul.R")

source("fn_fit_ramsay.R")
source("fn_fit_james.R")
source("fn_fit_mfpca.R")

source("fn_cvfit_ramjames.R")
source("fn_cvfit_mfpca.R")

library(fda)
library(MASS)
library(cvTools)

obj_ram = obj_jam = obj_mnf = list()

RAND = 100
theta_ram = theta_jam = theta_mnf = vector("list", RAND)
Theta_jam = Theta_mnf = vector("list", RAND)
ascore_jam = ascore_mnf = vector("list", RAND)
beta_ram = beta_mnf = vector("list", RAND)

for(rand in c(1:RAND)){
  tryCatch({
    obj_simul = fn_data_simul(50, seedn = 1234 + rand - 1, sig = 0.01)
    
    dat.sim = obj_simul$dat
    true.par = obj_simul$par
    
    ncurve = max(dat.sim$id)
    c_lmu = 10^{-5:1}; 
    
    print(paste("simulation", rand, "started"))
    
    ##### Approach 1: Ramsay's #####
    for(i in 1:ncurve){
      obj_cv = rep(NA, length(c_lmu))
      for(j in 1:length(c_lmu)){
        obj_cv[j] = cv.ram_opt(c_lmu[j], dat = dat.sim[dat.sim$id == i,], 
                               ini = true.par$theta, elim = .1, maxiter = 50)
      }
      opt_lmu = c_lmu[which(obj_cv == min(obj_cv, na.rm = TRUE))]
      obj_ram[[i]] = tryCatch(with(dat.sim, fit_ramsay(tobs[id == i], yobs[id == i], true.par$theta, norder = 5,
                                                       trng = c(0,1), lmu = opt_lmu, elim = 0.1, maxiter = 50)),
                              error = function(err) NULL)
    }
    
    theta_ram[[rand]] = t(sapply(obj_ram, function(x){x$est$theta})); 
    colnames(theta_ram[[rand]]) = paste("bspl5.", 1:10, sep = "")
    beta_ram[[rand]] = t(sapply(obj_ram, function(x){x$est$bcoef}))
    colnames(beta_ram[[rand]]) = paste("beta", 0:1, sep = "")
    
    print("  approach 1 estimated")
    
    ##### Approach 2: two-step #####
    dat.est = data_james(obj_ram)
    inipar_2 = ini_james(dat.est, nbasis = 10, norder = 5, npc = 2)
    
    obj_jam[[rand]] = fit_james(dat.est, inipar_2, nbasis = 10, norder = 5, npc = 2, elim = 0.1, maxiter = 30)
    
    theta_jam[[rand]] = obj_jam[[rand]]$est$theta
    Theta_jam[[rand]] = obj_jam[[rand]]$est$Theta;
    ascore_jam[[rand]] = obj_jam[[rand]]$est$ascore
    
    print("  approach 2 estimated")
    
    ##### Approach 3: integrated ####
    cv_2pc = cv.search(c_lmu = c_lmu, c_lpc = 10^{0:5}, K = 5,
                       dat = dat.sim, ini = true.par, elim = 0.1, maxiter = 30)
    obj_mnf[[rand]] = cv_2pc$opt_fit
    
    theta_mnf[[rand]] = obj_mnf[[rand]]$est$theta
    Theta_mnf[[rand]] = obj_mnf[[rand]]$est$Theta
    ascore_mnf[[rand]] = obj_mnf[[rand]]$est$ascore
    beta_mnf[[rand]] = obj_mnf[[rand]]$est$bcoef
    
    
    print("  approach 3 estimated")
    
    m1_ramsay = list(theta = theta_ram, beta = beta_ram)
    m2_twostep = list(theta = theta_jam, Theta = Theta_jam, ascore = ascore_jam, beta = beta_ram)
    m3_integrated = list(theta = theta_mnf, Theta = Theta_mnf, ascore = ascore_mnf, beta = beta_mnf)
    
    saveRDS(m1_ramsay, "simresult_compare/m1_01.rds")
    saveRDS(m2_twostep, "simresult_compare/m2_01.rds")
    saveRDS(m3_integrated, "simresult_compare/m3_01.rds")
    
  }, error = function(e){print(paste("skip simulation:", rand))})
}
