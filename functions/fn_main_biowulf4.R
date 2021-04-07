
simul_study = function(nsimul, output1, seed, sigv, alphv){
  
  RAND = nsimul
  
  true_par = M1 = M2 = M3 = list()
  true_par$beta = M1$beta = M2$beta = M3$beta = array(NA, dim = c(50, 2, RAND))
  M1$theta = array(NA, dim = c(50, 10, RAND))
  true_par$theta = M2$theta = M3$theta = array(NA, dim = c(10,1, RAND))
  true_par$Theta = M2$Theta = M3$Theta = array(NA, dim = c(10,2, RAND))
  true_par$ascore = M2$ascore = M3$ascore = array(NA, dim = c(50,2, RAND))
  
  rand = 1; re = 0
  set.seed(seed)
  repeat{
    re = re + 1
    # tryCatch({
      obj_simul = fn_data_simul(50, sig = sigv, alpha = alphv)
      
      dat.sim = obj_simul$dat
      true.par = obj_simul$par
      
      true_par$beta[,,rand] = true.par$bcoef
      true_par$theta[,,rand] = true.par$theta
      true_par$Theta[,,rand] = true.par$Theta
      true_par$ascore[,,rand] = true.par$ascore
      
      ncurve = max(dat.sim$id)
      c_lmu = 10^{-10:1};
      
      test.id = test.idi = list(NULL, NULL, NULL, NULL, NULL); K = 5
      # for(i in 1:ncurve){
      #   #set.seed(seedn)
      #   group.id = cvFolds(nrow(dat.sim[dat.sim$id == i,]), K)
      #   for(k in 1:K){
      #     test.id[[k]] = c(test.id[[k]], 
      #                      sum(dat.sim$id < i) + with(group.id, sort(subsets[which == k])))
      #   }
      # }
      
      cat("===== simulation", rand, "started =====", "\n")
      
      ##### Approach 1: Ramsay's #####
      obj_ram = vector("list", ncurve)
      exseq = rep(F, ncurve)
      for(i in 1:ncurve){
        obj_cv = rep(NA, length(c_lmu))
        group.id = cvFolds(nrow(dat.sim[dat.sim$id == i,]), K)
            for(k in 1:K){
              test.idi[[k]] = with(group.id, sort(subsets[which == k]))
              test.id[[k]] = c(test.id[[k]], 
                               sum(dat.sim$id < i) + with(group.id, sort(subsets[which == k])))
            }
        for(j in 1:length(c_lmu)){
          obj_cv[j] = cv.ram_opt(c_lmu[j], tid = test.idi, dat = dat.sim[dat.sim$id == i,],
                                 ini = true.par$theta, elim = .1, maxiter = 50)
        }
        opt_lmu = c_lmu[which(obj_cv == min(obj_cv, na.rm = TRUE))]
        obj_ram[[i]] = tmp = tryCatch(with(dat.sim, fit_ramsay(tobs[id == i], yobs[id == i], true.par$theta, norder = 5,
                                                               trng = c(0,1), lmu = opt_lmu, elim = 0.1, maxiter = 50)),
                                      error = function(err) NULL)
        
        if(is.null(tmp)) {exseq[i] = T}
        if(!is.null(tmp)) M1$beta[i,,rand] = obj_ram[[i]]$est$bcoef
        if(!is.null(tmp)) M1$theta[i,,rand] = obj_ram[[i]]$est$theta
        cat(i, ifelse(class(tmp)!="list", "x ", "- "))
      }
      cat("\n")
      if(sum(exseq)>0) next
      cat("  approach 1 estimated", "\n")
       
      ##### Approach 2: two-step #####
      # exid = which(exseq)
      dat.est = tmp = tryCatch(data_james(obj_ram), error = function(err) NULL)
      if(is.null(tmp)) next

      inipar_2 = ini_james(dat.est, nbasis = 10, norder = 5, npc = 2)

      # obj_jam = fit_james(dat.est, inipar_2, nbasis = 10, norder = 5, npc = 2, elim = 0.1, maxiter = 30)

      obj_jam = tmp = tryCatch(fit_james(dat.est, inipar_2, nbasis = 10, norder = 5, npc = 2, elim = 0.1, maxiter = 30),
                               error = function(err) NULL)

      if(is.null(tmp)) next

      M2$beta[,,rand] = bcoef.jam(obj_jam)
      M2$theta[,,rand] = obj_jam$est$theta
      M2$Theta[,,rand] = obj_jam$est$Theta;
      M2$ascore[,,rand] = obj_jam$est$ascore

      cat("  approach 2 estimated", "\n")

      ##### Approach 3: integrated ####
      cv_2pc = cv.search(c_lmu = c_lmu, c_lpc = 10^{-1:5}, K = 5,
                         tid = test.id, dat = dat.sim, 
                         ini = true.par, elim = 0.1, maxiter = 30)

      if(is.null(cv_2pc$opt_fit)) next

      obj_mnf = cv_2pc$opt_fit

      M3$theta[,,rand] = obj_mnf$est$theta
      M3$Theta[,,rand] = obj_mnf$est$Theta
      M3$ascore[,,rand] = obj_mnf$est$ascore
      M3$beta[,,rand] = obj_mnf$est$bcoef


      cat("  approach 3 estimated", "\n")


    # }, error = function(e){next})
    if (is.na(obj_mnf$est$theta[1,])) next
    if (rand == RAND | re == 10) break
    rand = rand + 1
  }
  res = list(par0 = true_par, m1 = M1, m2 = M2, m3 = M3)
  if(re == 10) res$ex = T
  saveRDS(res, file=output1)
}

simul_study_t = function(nsimul, output1, seed, sigv, alphv){
  
  RAND = nsimul
  
  true_par = M1 = M2 = M3 = list()
  true_par$beta = M1$beta = M2$beta = M3$beta = array(NA, dim = c(50, 2, RAND))
  M1$theta = array(NA, dim = c(50, 10, RAND))
  true_par$theta = M2$theta = M3$theta = array(NA, dim = c(10,1, RAND))
  true_par$Theta = M2$Theta = M3$Theta = array(NA, dim = c(10,2, RAND))
  true_par$ascore = M2$ascore = M3$ascore = array(NA, dim = c(50,2, RAND))
  
  rand = 1; re = 0
  set.seed(seed)
  repeat{
    re = re + 1
    # tryCatch({
    obj_simul = fn_data_tsimul(50, sig = sigv, alpha = alphv)
    
    dat.sim = obj_simul$dat
    true.par = obj_simul$par
    
    true_par$beta[,,rand] = true.par$bcoef
    true_par$theta[,,rand] = true.par$theta
    true_par$Theta[,,rand] = true.par$Theta
    true_par$ascore[,,rand] = true.par$ascore
    
    ncurve = max(dat.sim$id)
    c_lmu = 10^{-10:1};
    
    test.id = test.idi = list(NULL, NULL, NULL, NULL, NULL); K = 5
    cat("===== simulation", rand, "started =====", "\n")
    
    ##### Approach 1: Ramsay's #####
    obj_ram = vector("list", ncurve)
    exseq = rep(F, ncurve)
    for(i in 1:ncurve){
      obj_cv = rep(NA, length(c_lmu))
      group.id = cvFolds(nrow(dat.sim[dat.sim$id == i,]), K)
      for(k in 1:K){
        test.idi[[k]] = with(group.id, sort(subsets[which == k]))
        test.id[[k]] = c(test.id[[k]], 
                         sum(dat.sim$id < i) + with(group.id, sort(subsets[which == k])))
      }
      for(j in 1:length(c_lmu)){
        obj_cv[j] = cv.ram_opt(c_lmu[j], tid = test.idi, dat = dat.sim[dat.sim$id == i,],
                               ini = true.par$theta, elim = .1, maxiter = 50)
      }
      opt_lmu = c_lmu[which(obj_cv == min(obj_cv, na.rm = TRUE))]
      obj_ram[[i]] = tmp = tryCatch(with(dat.sim, fit_ramsay(tobs[id == i], yobs[id == i], true.par$theta, norder = 5,
                                                             trng = c(0,1), lmu = opt_lmu, elim = 0.1, maxiter = 50)),
                                    error = function(err) NULL)
      
      if(is.null(tmp)) {exseq[i] = T}
      if(!is.null(tmp)) M1$beta[i,,rand] = obj_ram[[i]]$est$bcoef
      if(!is.null(tmp)) M1$theta[i,,rand] = obj_ram[[i]]$est$theta
      cat(i, ifelse(class(tmp)!="list", "x ", "- "))
    }
    cat("\n")
    if(sum(exseq)>0) next
    cat("  approach 1 estimated", "\n")
    
    ##### Approach 2: two-step #####
    # exid = which(exseq)
    dat.est = tmp = tryCatch(data_james(obj_ram), error = function(err) NULL)
    if(is.null(tmp)) next
    
    inipar_2 = ini_james(dat.est, nbasis = 10, norder = 5, npc = 2)
    
    # obj_jam = fit_james(dat.est, inipar_2, nbasis = 10, norder = 5, npc = 2, elim = 0.1, maxiter = 30)
    
    obj_jam = tmp = tryCatch(fit_james(dat.est, inipar_2, nbasis = 10, norder = 5, npc = 2, elim = 0.1, maxiter = 30),
                             error = function(err) NULL)
    
    if(is.null(tmp)) next
    
    M2$beta[,,rand] = bcoef.jam(obj_jam)
    M2$theta[,,rand] = obj_jam$est$theta
    M2$Theta[,,rand] = obj_jam$est$Theta;
    M2$ascore[,,rand] = obj_jam$est$ascore
    
    cat("  approach 2 estimated", "\n")
    
    ##### Approach 3: integrated ####
    cv_2pc = cv.search(c_lmu = c_lmu, c_lpc = 10^{-1:5}, K = 5,
                       tid = test.id, dat = dat.sim, 
                       ini = true.par, elim = 0.1, maxiter = 30)
    
    if(is.null(cv_2pc$opt_fit)) next
    
    obj_mnf = cv_2pc$opt_fit
    
    M3$theta[,,rand] = obj_mnf$est$theta
    M3$Theta[,,rand] = obj_mnf$est$Theta
    M3$ascore[,,rand] = obj_mnf$est$ascore
    M3$beta[,,rand] = obj_mnf$est$bcoef
    
    
    cat("  approach 3 estimated", "\n")
    
    
    # }, error = function(e){next})
    if (rand == RAND | re == 10) break
    if (is.na(obj_mnf$est$theta[1,])) next
    rand = rand + 1
  }
  res = list(par0 = true_par, m1 = M1, m2 = M2, m3 = M3)
  if(re == 10) res$ex = T
  saveRDS(res, file=output1)
}


