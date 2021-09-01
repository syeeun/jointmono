#### orthogonalizing operator ####

orthobsp = function(rng, nbasis, norder, L = 1000){
  
  #rng = c(rng0[1]-diff(rng0)*.2, rng0[2]+diff(rng0)*.2)
  bsp = create.bspline.basis(rng, nbasis, norder) 
  
  # Orthogonalizing B-spline
  twd = diff(rng)
  tgrid = seq(from = rng[1], to = rng[2], length.out=L)
  Bgrid = eval.basis(tgrid, bsp)
  G = t(Bgrid) %*% Bgrid
  opG = solve(chol(G))*sqrt(L/twd)
  
  return(opG)
}

#### b(t) & B(t): nbasis x tgrid ####
o.bsp.mat = function(t, nbasis, norder, integ = F, nderiv = 0, finelevel = 4){
  
  rng = range(t) 
  
  bsp_obj = create.bspline.basis(rng, nbasis, norder)
  ortho_op = orthobsp(rng, nbasis, norder)
  
  if(integ == F){
    bspmat = getbasismatrix(t, bsp_obj, nderiv, returnMatrix = T)
    result = bspmat%*%ortho_op
  }
  
  if(integ == T){
    tfine = seq(rng[1],rng[2],l=10^finelevel)
    delta = diff(tfine)[1]
    bspmat_fine = getbasismatrix(tfine, bsp_obj, nderiv, returnMatrix = T)
    fcumsum = function(fval, tfine){ 
      fval_n = fval[tfine < 0]; fval_p = fval[tfine > 0 | tfine == 0]
      result_n = cumsum(fval_n) - .5*(fval_n[length(fval_n)]+fval_n)
      result_p = cumsum(fval_p) - .5*(fval_p[1]+fval_p)
      return(c(result_n - result_n[length(result_n)], result_p))
      }
    Bspmat_fine = delta*apply(bspmat_fine%*%ortho_op, 2, fcumsum, tfine = tfine)
    approx_y = function(x, y, xout){ approx(x, y, xout)$y }
    result = apply(Bspmat_fine, 2, approx_y, x=tfine, xout = t)
  }
  
  return(list(tgrid = t, obspmat = result))
}

# 
# #### Penalty matrix of orthogonal basis ####
# o.pen.mat = function(rng0, nbasis, norder, L = 10000){
# 
#   rng = c(rng0[1]-diff(rng0)*.2, rng0[2]+diff(rng)*.2)
#   tfine = seq(rng[1], rng[2], l = L)
#   twidth = diff(rng)
#   
#   o_bsp_mat = o.bsp.mat(tfine, nbasis, norder, nderiv = 2)$obspmat
#   
#   pen_mat = (t(o_bsp_mat)%*%o_bsp_mat)*(twidth/L)
#   
#   return(pen_mat)
# }