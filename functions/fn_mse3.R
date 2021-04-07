f.mse.y_ram = function(m_est, true, tau){
  
  RAND = dim(m_est$beta)[3]; ncurve = dim(true$beta)[1]
  tgrid = seq(0,1,l=500); delta = diff(tgrid)[1]
  oBsp = o.bsp.mat(tgrid, 10, 5, T)
  
  y0 = array(NA, dim = c(ncurve, RAND, length(tgrid)))
  YHAT = array(NA, dim = c(ncurve, RAND, length(tgrid)))
  for(i in 1:ncurve){
    for(j in 1:RAND){
      # if (any(is.na(m_est$beta[i,,j]))) next
      YHAT[i,j,] = with(m_est, beta[i,1,j] + beta[i,2,j]*Hfn(tgrid, theta[i,,j], oBsp))
      y0[i,j,] = with(true, Hfn(tgrid, theta[,,j] + Theta[,,j] %*% ascore[i,,j], oBsp))# + Theta[,,j] %*% ascore[i,,j]
    }
  }
  
  Bias2 = apply((apply(YHAT - y0, c(1,3), mean))^2, 1,
                function(x){sum(x[tgrid>tau[1] & tgrid<tau[2]])*delta})
  Var = apply(apply((sweep(YHAT, c(1,3), apply(YHAT, c(1,3), mean), "-"))^2, 
                    c(1,3), mean), 1, function(x){sum(x[tgrid>tau[1] & tgrid<tau[2]])*delta})
  MISE = apply(apply((YHAT-y0)^2, c(1,3), mean),1,function(x){sum(x[tgrid>tau[1] & tgrid<tau[2]])*delta})
  
  # Bias2 = mean((apply(apply(sweep(YHAT, c(1,3), y0, "-"), c(1,3), mean), 1, function(x){sum(x)*delta}))^2)
  # Var = mean(apply(apply((sweep(YHAT, c(1,3), apply(YHAT, 2, mean), "-"))^2, c(1,3), mean), 1, function(x){sum(x)*delta}))
  
  return(list(Bias2=Bias2, Var = Var, MISE = MISE))
}

f.mse.y = function(m_est, true, tau, ntime){
  
  RAND = dim(m_est$beta)[3]; ncurve = dim(true$beta)[1]
  tgrid = seq(0,1,l=ntime); delta = diff(tgrid)[1]
  oBsp = o.bsp.mat(tgrid, 10, 5, T)
  
  y0 = array(NA, dim = c(ncurve, RAND, length(tgrid)))
  YHAT = array(NA, dim = c(ncurve, RAND, length(tgrid)))
  for(i in 1:ncurve){
    for(j in 1:RAND){
      # if (any(is.na(m_est$beta[i,,j]), is.na(m_est$theta[,,j]))) next
      YHAT[i,j,] = with(m_est, beta[i,1,j] + beta[i,2,j]*Hfn(tgrid, theta[,,j]+ Theta[,,j] %*% ascore[i,,j], oBsp)) # 
      y0[i,j,] = with(true, Hfn(tgrid, theta[,,j]+ Theta[,,j] %*% ascore[i,,j], oBsp)) # + Theta[,,j] %*% ascore[i,,j]
    }
  }
  
  Bias2 = apply((apply(YHAT - y0, c(1,3), mean))^2, 1,
                function(x){sum(x[tgrid>tau[1] & tgrid<tau[2]])*delta})
  Var = apply(apply((sweep(YHAT, c(1,3), apply(YHAT, c(1,3), mean), "-"))^2, 
                    c(1,3), mean), 1, function(x){sum(x[tgrid>tau[1] & tgrid<tau[2]])*delta})
  MISE = apply(apply((YHAT-y0)^2, c(1,3), mean),1,function(x){sum(x[tgrid>tau[1] & tgrid<tau[2]])*delta})

  return(list(Bias2=Bias2, Var = Var, MISE = MISE))
}


f.mse.w_ram = function(m_est, true){
  
  RAND = dim(m_est$beta)[3]; 
  ncurve = dim(true$beta)[1]
  norder = dim(m_est$theta)[2]
  
  w0 = array(NA, dim = c(ncurve, RAND, norder))
  WHAT = array(NA, dim = c(ncurve, RAND, norder))
  for(i in 1:ncurve){
    for(j in 1:RAND){
      WHAT[i,j,] = with(m_est, theta[i,,j])
      w0[i,j,] = with(true, theta[,,j]+ Theta[,,j] %*% ascore[i,,j]) # + Theta[,,j] %*% ascore[i,,j]
    }
  }
  
  Bias2 = apply(apply(WHAT - w0, c(1,3), mean), 1, crossprod)
  Var = rowMeans(apply(sweep(WHAT, c(1,3), apply(WHAT, c(1,3), mean), "-"), c(1,2), crossprod))
  MISE = apply(apply(WHAT-w0,c(1,2),crossprod), 1, mean)
  
  return(list(Bias2 = Bias2, Var = Var, MISE = MISE))
}

f.mse.w = function(m_est, true){
  
  RAND = dim(m_est$beta)[3]; 
  ncurve = dim(m_est$beta)[1]
  norder = dim(m_est$theta)[1]
  
  w0 = array(NA, dim = c(ncurve, RAND, norder))
  WHAT = array(NA, dim = c(ncurve, RAND, norder))
  for(i in 1:ncurve){
    for(j in 1:RAND){
      WHAT[i,j,] = with(m_est, theta[,,j] + Theta[,,j] %*% ascore[i,,j])
      w0[i,j,] = with(true, theta[,,j] + Theta[,,j] %*% ascore[i,,j]) # + Theta %*% ascore[i,]
    }
  }
  
  Bias2 = apply(apply(WHAT - w0, c(1,3), mean), 1, crossprod)
  Var = rowMeans(apply(sweep(WHAT, c(1,3), apply(WHAT, c(1,3), mean), "-"), c(1,2), crossprod))
  MISE = apply(apply(WHAT-w0,c(1,2),crossprod), 1, mean)
  
  return(list(Bias2 = Bias2, Var = Var, MISE = MISE))
}


table.mse = function(m1, m2, m3){
  res = rbind(c(mean(m1$Bias2), mean(m2$Bias2), mean(m3$Bias2)), 
                c(mean(m1$Var), mean(m2$Var), mean(m3$Var)),
                  c(mean(m1$MISE), mean(m2$MISE), mean(m3$MISE)))
  colnames(res) = c("Ramsay", "Two-step", "Integrated")
  rownames(res) = c("Bias squared", "Variance", "MISE")
  return(res)
}


boxplot.mse = function(m1, m2, m3, plotit = c("bias2", "var", "mse"), sigv, ...){
  if ("bias2" %in% plotit){
    boxplot(sqrt(c(m1$Bias2)), sqrt(c(m2$Bias2)), sqrt(c(m3$Bias2)), 
            outline = F, col = "grey", 
            main = sigv, #expression(paste(sigma, " = ", sigv)), 
            # ylim = range(c(log(c(m1$Bias2)), log(c(m2$Bias2)), log(c(m3$Bias2)))), #cex.main = 1.5,
            ylab = expression(paste("BIAS","(",hat(m),")")),
            ...)}
  if ("var" %in% plotit){
    boxplot(sqrt(c(m1$Var)), sqrt(c(m2$Var)), sqrt(c(m3$Var)), 
            outline = F, col = "grey", 
            main = sigv, #expression(paste(sigma, " = ", sigv)), 
            # ylim = range(c(log(c(m1$Var)), log(c(m2$Var)), log(c(m3$Var)))), #cex.main = 1.5,
            ylab = expression(paste("SD","(",hat(m),")")),
            ...)}
  if ("mse" %in% plotit){
    boxplot(sqrt(m1$MISE), sqrt(m2$MISE), sqrt(m3$MISE), 
            outline = F, col = "grey", 
            main = sigv, #expression(paste(sigma, " = ", sigv)), 
            # ylim = range(c(log(m1$Bias2+m2$Var), log(m2$Bias2+m2$Var), log(m3$Bias2+m3$Var))), #cex.main = 1.5,
            ylab = expression(paste("sqrt-MISE","(",hat(m),")")),
            ...)}
  
}
