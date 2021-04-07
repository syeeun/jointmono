cv.search = function(c_lmu, c_lpc, K = 5, dat, tid, ini, ...){
  
  # ncurve = max(dat$id)
  # test.id = list(NULL, NULL, NULL, NULL, NULL)
  # for(i in 1:ncurve){
  # #set.seed(seedn)
  # group.id = cvFolds(nrow(dat[dat$id == i,]), K)
  # for(k in 1:K){
  # test.id[[k]] = c(test.id[[k]], 
  #                  sum(dat$id < i) + with(group.id, sort(subsets[which == k])))
  # }
  # }
  
  obj_fold = list()
  for(k in 1:K){
    # test.id = with(group.id, sort(subsets[which == k]))
    dat_train = dat[-tid[[k]], ] 
    dat_test = dat[tid[[k]], ]
  
    obj_cv = matrix(NA, length(c_lmu), length(c_lpc))
    for(i in 1:length(c_lmu)){
      for(j in 1:length(c_lpc)){
        obj_cv[i,j] = tryCatch(cv.error_test(obj_fitted = fit_mfpca(dat_train, ini, norder = 5,
                                                          lmu = c_lmu[i], lpc = c_lpc[j], ...),#, ...), 
                                              dat_test = dat_test),
                               error = function(err) NA)
        #print(paste("K = ", k, "; lmu = ", round(c_lmu[i], 4), ", lpc = ", round(c_lpc[j], 4), 
        #            "; error = ", round(obj_cv[i,j], 4), sep = ""))
      }
    }
    if(prod(is.na(obj_cv)) == 1) {
      print("nothing converges; try other grids")
      obj_cv = NA
    }
    obj_fold[[k]] = obj_cv
    #print(paste("proposed; K = ", k))
  }

  obj_fold = simplify2array(obj_fold)
  if(!is.null(dim(obj_fold))){
  obj_mean = apply(obj_fold, 1:2, mean, na.rm = T)
    rownames(obj_mean) = round(c_lmu, 4); colnames(obj_mean) = round(c_lpc, 4)
    id = which(obj_mean == min(obj_mean, na.rm = T), T)
    opt_lambda = c(c_lmu[id[1]], c_lpc[id[2]])
    opt_fit = tryCatch(fit_mfpca(dat, ini, norder = 5,
                       lmu = c_lmu[id[1]], lpc = c_lpc[id[2]], elim = 0.1, maxiter = 30),
                       error = function(err) NULL)
  }
  if(is.null(dim(obj_fold))){
    opt_lambda = obj_mean = opt_fit = NULL
  }
    return(list(opt_lambda = opt_lambda,
                opt_cv_res = obj_mean,
                opt_fit = opt_fit))
  }


cv.error_test = function(obj_fitted, dat_test, trng = NULL){
  curveid = as.numeric(levels(as.factor(dat_test$id)))
  ncurve = length(curveid)
  if(is.null(trng)) {trng = range(dat_test$tobs)}
  rng = trng
  tfine = seq(rng[1], rng[2], l = 1000)
  
  theta = obj_fitted$est$theta
  Theta = obj_fitted$est$Theta
  ascore = obj_fitted$est$ascore
  bcoef = obj_fitted$est$bcoef
  Bsp = obj_fitted$basis$Bsp
  
  ressq = 0
  for(i in curveid){
    wcoef = theta + Theta %*% ascore[i,]
    Hmat = Hfn(tfine, wcoef, Bsp)
    Hval = with(dat_test, approx(tfine, Hmat, tobs[id==i])$y)
    ressq = ressq + with(dat_test, sum((yobs[id==i] - bcoef[i,1] - bcoef[i,2]*Hval)^2))
  }
  cv_error = sqrt(ressq/nrow(dat_test)); 
  
  return(cv_error)
}

plot_cv_res = function(obj){
  image(obj, axes = F, xlab = expression(lambda[mu]), ylab = expression(lambda[f]))
  axis(1, (0:(nrow(obj)-1))/(nrow(obj)-1), rownames(obj))
  axis(2, (0:(ncol(obj)-1))/(ncol(obj)-1), colnames(obj))
  contour(obj, add = TRUE, drawlabels = TRUE)
  
  id = which(obj == min(obj, na.rm = T), T)
  points((id[1]-1)/(nrow(obj)-1), (id[2]-1)/(ncol(obj)-1), pch = 8)
}


fit.error= function(obj_w, obj_r, dat){
  if(ncol(obj_w$track$theta.storage)<30){
    tfine = obj_w$basis$bsp$tgrid
    ressq = 0
    for(i in 1:ncurve){
      wcoef = with(obj_w$est, theta + Theta %*% ascore[i,])
      Hmat = Hfn(tfine, wcoef, obj_w$basis$Bsp)
      Hval = with(dat, approx(tfine, Hmat, tobs[id==i])$y)
      ressq = ressq + with(dat, sum((yobs[id==i] - obj_r[[i]]$est$bcoef[1] - obj_r[[i]]$est$bcoef[2]*Hval)^2))
    }
    sigma = sqrt(ressq/nrow(dat));
  }
  else sigma = Inf
  return(sigma)
}

# 
# 
# f.mse.mu = function(list_obj_fitted, true.par){
#   
#   seqlist = which(sapply(list_obj_fitted, function(x){is.null(x$opt_fit)}))
#   result = 0
#   for(rand in (1:length(list_obj_fitted))[-seqlist]){
#     muifn = with(list_obj_fitted[[rand]],
#                  function(t){(Wfn(t, opt_fit$est$theta-true.par$theta, opt_fit$basis$bsp))^2})
#     result = result + integrate(muifn, 0, 1, subdivisions = 500)$value
#   }
#   mse = result/(length(list_obj_fitted)-length(seqlist))
#   
#   return(mse)
# }