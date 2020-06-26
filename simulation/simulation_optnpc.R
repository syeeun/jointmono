source("fn_ortho.R")
source("fn_ftnlz.R")
source("fn_data_simul.R")

source("fn_fit_mfpca.R")

source("fn_cvfit_mfpca.R")

library(fda)
library(MASS)
library(cvTools)

opK.01 = opK.05 = opK.10 = NULL

RAND = 100

c_lmu = 10^{-5:1}; c_lpc = 10^{1:5}

for(rand in c(1:RAND)){
  tryCatch({
    obj_simul = fn_data_simul(50, seedn = 1234 + rand - 1, sig = 0.01)
    
    dat.sim = obj_simul$dat
    true.par = obj_simul$par
    ini1 = ini2 = ini3 = true.par
    
    ini1$Theta = ini1$Theta[,1,drop=F]; ini1$ascore = ini1$ascore[,1,drop=F]
    ini3$Theta = cbind(ini3$Theta, theta_f3 = matrix(.1, 10, 1)); 
    ini3$ascore = cbind(ini3$ascore, matrix(0, 50, 1))
    
    ncurve = max(dat.sim$id)
    
    print(paste("[sigma = 0.01]", "simulation =", rand, "started"))
    ##### PROPOSED ####
    cv_1pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini1, elim = 0.1, maxiter = 30)
    cv_2pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini2, elim = 0.1, maxiter = 30)
    cv_3pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini3, elim = 0.1, maxiter = 30)
    cv_mnf = list(cv_1pc, cv_2pc, cv_3pc)
    
    errmnf = c(min(cv_1pc$opt_cv_res, na.rm=T), min(cv_2pc$opt_cv_res, na.rm=T), min(cv_3pc$opt_cv_res, na.rm=T))
    opK2 = which(errmnf == min(errmnf, na.rm = T))

    opK.01 = c(opK.01, opK2)
    
    saveRDS(opK.01, "simresult_optnpc/optnpc_01.rds")
    
  }, error = function(e){print(paste("skip simulation:", rand))})
}


for(rand in c(1:RAND)){
  tryCatch({
    obj_simul = fn_data_simul(50, seedn = 1234 + rand - 1, sig = 0.05)
    
    dat.sim = obj_simul$dat
    true.par = obj_simul$par
    ini1 = ini2 = ini3 = true.par
    
    ini1$Theta = ini1$Theta[,1,drop=F]; ini1$ascore = ini1$ascore[,1,drop=F]
    ini3$Theta = cbind(ini3$Theta, theta_f3 = matrix(.1, 10, 1)); 
    ini3$ascore = cbind(ini3$ascore, matrix(0, 50, 1))
    
    ncurve = max(dat.sim$id)
    
    print(paste("[sigma = 0.05]", "simulation =", rand, "started"))
    ##### PROPOSED ####
    cv_1pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini1, elim = 0.1, maxiter = 30)
    cv_2pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini2, elim = 0.1, maxiter = 30)
    cv_3pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini3, elim = 0.1, maxiter = 30)
    cv_mnf = list(cv_1pc, cv_2pc, cv_3pc)
    
    errmnf = c(min(cv_1pc$opt_cv_res, na.rm=T), min(cv_2pc$opt_cv_res, na.rm=T), min(cv_3pc$opt_cv_res, na.rm=T))
    opK2 = which(errmnf == min(errmnf, na.rm = T))
    
    opK.05 = c(opK.05, opK2)
    
    saveRDS(opK.05, "simresult_optnpc/optnpc_05.rds")
    
  }, error = function(e){print(paste("skip simulation:", rand))})
}


for(rand in c(1:RAND)){
  tryCatch({
    obj_simul = fn_data_simul(50, seedn = 1234 + rand - 1, sig = 0.10)
    
    dat.sim = obj_simul$dat
    true.par = obj_simul$par
    ini1 = ini2 = ini3 = true.par
    
    ini1$Theta = ini1$Theta[,1,drop=F]; ini1$ascore = ini1$ascore[,1,drop=F]
    ini3$Theta = cbind(ini3$Theta, theta_f3 = matrix(.1, 10, 1)); 
    ini3$ascore = cbind(ini3$ascore, matrix(0, 50, 1))
    
    ncurve = max(dat.sim$id)
    
    print(paste("[sigma = 0.10]", "simulation =", rand, "started"))
    ##### PROPOSED ####
    cv_1pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini1, elim = 0.1, maxiter = 30)
    cv_2pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini2, elim = 0.1, maxiter = 30)
    cv_3pc = cv.search(c_lmu = c_lmu, c_lpc = c_lpc, K = 5, dat = dat.sim, ini = ini3, elim = 0.1, maxiter = 30)
    cv_mnf = list(cv_1pc, cv_2pc, cv_3pc)
    
    errmnf = c(min(cv_1pc$opt_cv_res, na.rm=T), min(cv_2pc$opt_cv_res, na.rm=T), min(cv_3pc$opt_cv_res, na.rm=T))
    opK2 = which(errmnf == min(errmnf, na.rm = T))
    
    opK.10 = c(opK.10, opK2)
    
    saveRDS(opK.10, "simresult_optnpc/optnpc_10.rds")
    
  }, error = function(e){print(paste("skip simulation:", rand))})
}

