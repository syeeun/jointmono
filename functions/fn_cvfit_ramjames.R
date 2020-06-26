cv.ram_opt = function(par0, K = 5, seedn = 1234, dat, ini, ...){
  
  set.seed(seedn)
  group.id = cvFolds(nrow(dat), K)
  
  obj_fold = rep(NA, K)
  for(k in 1:K){
    test.id = with(group.id, sort(subsets[which == k]))
    dat_train = dat[-test.id, ] 
    dat_test = dat[test.id, ]
    
    obj_fold[k] = tryCatch(cv.ram_error_test(obj_fitted = fit_ramsay(dat_train$tobs, dat_train$yobs, ini, norder = 5,
                                                                 lmu = par0, ...),#, ...), 
                                         dat_test = dat_test),
                           error = function(err) NA)
  }
  
  obj_mean = mean(obj_fold, na.rm = T)
  return(obj_mean)
}


cv.ram_error_test = function(obj_fitted, dat_test, trng = NULL){
  # curveid = as.numeric(levels(as.factor(dat_test$id)))
  # ncurve = length(curveid)
  if(is.null(trng)) {trng = range(dat_test$tobs)}
  rng = trng
  tfine = seq(rng[1], rng[2], l = 1000)
  
  theta = obj_fitted$est$theta
  bcoef = obj_fitted$est$bcoef
  Bsp = obj_fitted$basis$Bsp
  
    wcoef = theta# + Theta %*% ascore[i,]
    Hmat = Hfn(tfine, wcoef, Bsp)
    Hval = with(dat_test, approx(tfine, Hmat, tobs)$y)
    ressq = with(dat_test, mean((yobs - bcoef[1] - bcoef[2]*Hval)^2))
    
  return(sqrt(ressq))
}

f.mse.mu = function(list_obj_fitted, true.par){
  
  result = 0
  for(rand in 1:length(list_obj_fitted)){
    muifn = with(list_obj_fitted[[rand]],
                 function(t){(Wfn(t, est$theta-true.par$theta, basis$bsp))^2})
    result = result + integrate(muifn, 0, 1, subdivisions = 500)$value
  }
  mse = result/length(list_obj_fitted)
  
  return(mse)
}

f.mse.w = function(list_obj_fitted, true.par, isramsay = F){
  
  mse = NULL
  result = k = 0
  for(rand in 1:length(list_obj_fitted)){
    if(!(is.null(list_obj_fitted[[rand]]))){
    resultcomp = 0
    ncurve = ifelse(isramsay, length(list_obj_fitted[[rand]]), max(list_obj_fitted[[rand]]$dat$id))
    for(i in 1:ncurve){
      if(isramsay == T){
        wcoef = with(list_obj_fitted[[rand]][[i]]$est, (theta-true.par$theta))
      }
      if(isramsay == F){
        wcoef = with(list_obj_fitted[[rand]]$est, 
                 (theta-true.par$theta) + (Theta%*%ascore[i,] - true.par$Theta %*% true.par$ascore[i,]))
    #muifn = with(list_obj_fitted[[rand]]$basis, function(t){(Wfn(t, wcoef, bsp))^2})
      }
    resultcomp = resultcomp + t(wcoef)%*%wcoef #integrate(muifn, 0, 1, subdivisions = 500)$value
    }
    result = resultcomp/ncurve
    mse = c(mse, as.numeric(result))
    }
    #if(is.null(list_obj_fitted[[rand]])){k = k+1}
  }
  #/(length(list_obj_fitted)-k)
  
  return(mse)
}

f.mse.y = function(list_obj_fitted, true.par, isramsay = F){

  mse = NULL
  result = k = 0
  for(rand in 1:length(list_obj_fitted)){
    if(!(is.null(list_obj_fitted[[rand]]))){
    resultcomp = 0
    ncurve = ifelse(isramsay, length(list_obj_fitted[[rand]]), max(list_obj_fitted[[rand]]$dat$id))
    for(i in 1:ncurve){
      if(isramsay == T){
        delta = diff(list_obj_fitted[[rand]][[i]]$basis$Bsp$tgrid)[1]
        whatcoef = list_obj_fitted[[rand]][[i]]$est$theta
        wcoef = true.par$theta
        yhat = with(list_obj_fitted[[rand]][[i]], est$bcoef[1] + est$bcoef[2]*Hfn(basis$Bsp$tgrid, whatcoef, basis$Bsp))
        ytrue = with(list_obj_fitted[[rand]][[i]], true.par$bcoef[i,1] + true.par$bcoef[i,2]*Hfn(basis$Bsp$tgrid, wcoef, basis$Bsp))
      }
      if(isramsay == F){
        delta = diff(list_obj_fitted[[rand]]$basis$Bsp$tgrid)[1]
        whatcoef = with(list_obj_fitted[[rand]]$est, theta + Theta%*%ascore[i,])
        wcoef = with(true.par, theta + Theta%*%ascore[i,])
        yhat = with(list_obj_fitted[[rand]], est$bcoef[i,1] + est$bcoef[i,2]*Hfn(basis$Bsp$tgrid, whatcoef, basis$Bsp))
        ytrue = with(list_obj_fitted[[rand]], true.par$bcoef[i,1] + true.par$bcoef[i,2]*Hfn(basis$Bsp$tgrid, wcoef, basis$Bsp))
      }
      resultcomp = resultcomp + sum((yhat-ytrue)^2)*delta
    }
    result = resultcomp/ncurve
    mse = c(mse, result)
    }
    #if(is.null(list_obj_fitted[[rand]])){k = k+1}
  }
  #mse = result/(length(list_obj_fitted)-k)
  
  
  return(mse)
}


f_mse_w = function(estw, truew){ as.numeric(crossprod(estw-truew)) }

