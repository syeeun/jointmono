
Hfn = function(tval, wcoef, Bsp){
  Wobj = (Bsp$obspmat)%*%wcoef
  eWobj = exp(Wobj)
  dlt = diff((Bsp$tgrid))[1]
  
  eWobj_n = eWobj[Bsp$tgrid < 0]; eWobj_p = eWobj[Bsp$tgrid >= 0]
  result_n = dlt*(cumsum(eWobj_n) - .5*(eWobj_n[length(eWobj_n)]+eWobj_n))
  result_p = dlt*(cumsum(eWobj_p) - .5*(eWobj_p[1]+eWobj_p))
  Hval = c(result_n - result_n[length(result_n)], result_p)
  
  result = approx(Bsp$tgrid, Hval, tval)$y
  return(result)
}

Wfn = function(tval, wcoef, bsp){
  wobj = (bsp$obspmat)%*%wcoef
  result = approx(bsp$tgrid, wobj, tval)$y
  return(result)
}


dHfn = function(tval, wcoef, Bsp){
    
  nbasis = ncol(wcoef)
  Wobj = (Bsp$obspmat)%*%wcoef
  eWobj = matrix(exp(Wobj), nrow(Bsp$obspmat), nbasis)
  BeWobj = apply(Bsp$obspmat, 2, function(M){M*eWobj})
  
    dlt = diff(Bsp$tgrid)[1]
    fcumsum = function(fval, tfine){ 
      fval_n = fval[tfine < 0]; fval_p = fval[tfine >= 0]
      result_n = cumsum(fval_n) - .5*(fval_n[length(fval_n)]+fval_n)
      result_p = cumsum(fval_p) - .5*(fval_p[1]+fval_p)
      return(c(result_n - result_n[length(result_n)], result_p))
    }
    dHval = dlt*apply(BeWobj, 2, fcumsum, tfine = Bsp$tgrid) 
  
  fapprox = function(fval){ approx(Bsp$tgrid, fval, tval)$y }
  result = apply(dHval, 2, fapprox)

  return(result)
}