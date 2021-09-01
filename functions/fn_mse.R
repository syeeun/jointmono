f.mse.y_ram = function(m_est, true){
  
  RAND = length(m_est$beta); ncurve = nrow(true[[1]]$bcoef)
  tgrid = seq(0,1,l=120); delta = diff(tgrid)[1]
  oBsp = o.bsp.mat(tgrid, 10, 5, T)
  
  y0 = matrix(NA, ncurve, length(tgrid))
  YHAT = array(NA, dim = c(ncurve, RAND, length(tgrid)))
  for(i in 1:ncurve){
    for(j in 1:RAND){
      YHAT[i,j,] = with(m_est, beta[[j]][i,1] + beta[[j]][i,2]*Hfn(tgrid, theta[[j]][i,], oBsp))
    }
    y0[i,] = with(true[[1]], bcoef[i,1] + bcoef[i,2]*Hfn(tgrid, theta, oBsp))# + Theta %*% ascore[i,]
  }

  Bias2 = apply((apply(YHAT, c(1,3), mean) - y0)^2, 1, function(x){sum(x)*delta})
  Var = apply(apply((sweep(YHAT, c(1,3), apply(YHAT, c(1,3), mean), "-"))^2, c(1,3), mean), 1, function(x){sum(x)*delta})
  
  # Bias2 = mean((apply(apply(sweep(YHAT, c(1,3), y0, "-"), c(1,3), mean), 1, function(x){sum(x)*delta}))^2)
  # Var = mean(apply(apply((sweep(YHAT, c(1,3), apply(YHAT, 2, mean), "-"))^2, c(1,3), mean), 1, function(x){sum(x)*delta}))
  
  return(list(Bias2=Bias2, Var = Var))
}

f.mse.y = function(m_est, true){
  
  RAND = length(m_est$beta); ncurve = nrow(true[[1]]$bcoef)
  tgrid = seq(0,1,l=120); delta = diff(tgrid)[1]
  oBsp = o.bsp.mat(tgrid, 10, 5, T)
  
  y0 = matrix(NA, ncurve, length(tgrid))
  YHAT = array(NA, dim = c(ncurve, RAND, length(tgrid)))
  for(i in 1:ncurve){
    for(j in 1:RAND){
      YHAT[i,j,] = with(m_est, beta[[j]][i,1] + beta[[j]][i,2]*Hfn(tgrid, theta[[j]] + Theta[[j]] %*% ascore[[j]][i,], oBsp))
    }
    y0[i,] = with(true[[1]], bcoef[i,1] + bcoef[i,2]*Hfn(tgrid, theta + Theta %*% ascore[i,], oBsp))
  }
  
  Bias2 = apply((apply(YHAT, c(1,3), mean) - y0)^2, 1, function(x){sum(x)*delta})
  Var = apply(apply((sweep(YHAT, c(1,3), apply(YHAT, c(1,3), mean), "-"))^2, c(1,3), mean), 1, function(x){sum(x)*delta})
  # Bias2 = mean((apply(apply(sweep(YHAT, c(1,3), y0, "-"), c(1,3), mean), 1, function(x){sum(x)*delta}))^2)
  # Var = mean(apply(apply((sweep(YHAT, c(1,3), apply(YHAT, 2, mean), "-"))^2, c(1,3), mean), 1, function(x){sum(x)*delta}))
  
  return(list(Bias2=Bias2, Var = Var))
}


f.mse.w_ram = function(m_est, true){
  
  RAND = length(m_est$beta); 
  ncurve = nrow(true[[1]]$bcoef)
  norder = length(m_est$theta[[1]])
  
  w0 = matrix(NA, ncurve, norder)
  WHAT = array(NA, dim = c(ncurve, RAND, norder))
  for(i in 1:ncurve){
    for(j in 1:RAND){
      WHAT[i,j,] = with(m_est, theta[[j]][i,])
    }
    w0[i,] = with(true[[1]], theta) # + Theta %*% ascore[i,]
  }
  
  Bias2 = apply(apply(WHAT, c(1,3), mean) - w0, 1, crossprod)
  Var = rowMeans(apply(sweep(WHAT, c(1,3), apply(WHAT, c(1,3), mean), "-"), c(1,2), crossprod))
  
  return(list(Bias2 = Bias2, Var = Var))
}

f.mse.w = function(m_est, true){
  
  RAND = length(m_est$beta); 
  ncurve = nrow(true[[1]]$bcoef)
  norder = length(m_est$theta[[1]])

  w0 = matrix(NA, ncurve, norder)
  WHAT = array(NA, dim = c(ncurve, RAND, norder))
  for(i in 1:ncurve){
    for(j in 1:RAND){
      WHAT[i,j,] = with(m_est, theta[[j]] + Theta[[j]] %*% ascore[[j]][i,])
    }
    w0[i,] = with(true[[1]], theta + Theta %*% ascore[i,])
  }
  
  Bias2 = apply(apply(WHAT, c(1,3), mean) - w0, 1, crossprod)
  Var = rowMeans(apply(sweep(WHAT, c(1,3), apply(WHAT, c(1,3), mean), "-"), c(1,2), crossprod))
  
  return(list(Bias2 = Bias2, Var = Var))
}


table.mse = function(m1, m2, m3){
  res = rbind(colMeans(cbind(m1$Bias2, m2$Bias2, m3$Bias2)), 
              colMeans(cbind(m1$Var, m2$Var, m3$Var)))
  colnames(res) = c("Ramsay", "Two-step", "Integrated")
  rownames(res) = c("Bias squared", "Variance")
  return(res)
}


boxplot.mse = function(m1, m2, m3, plotit = c("bias2", "var", "mse"), sigv){
  if ("bias2" %in% plotit){
    boxplot(log(c(m1$Bias2)), log(c(m2$Bias2)), log(c(m3$Bias2)), 
            outline = F, col = "grey", 
            main = sigv, #expression(paste(sigma, " = ", sigv)), 
            ylim = range(c(log(c(m1$Bias2)), log(c(m2$Bias2)), log(c(m3$Bias2)))), #cex.main = 1.5,
            ylab = expression(paste("log BIAS2","(",hat(m),")")))}
  if ("var" %in% plotit){
    boxplot(log(c(m1$Var)), log(c(m2$Var)), log(c(m3$Var)), 
            outline = F, col = "grey", 
            main = sigv, #expression(paste(sigma, " = ", sigv)), 
            ylim = range(c(log(c(m1$Var)), log(c(m2$Var)), log(c(m3$Var)))), #cex.main = 1.5,
            ylab = expression(paste("log VAR","(",hat(m),")")))}
  if ("mse" %in% plotit){
    boxplot(log(m1$Bias2+m2$Var), log(m2$Bias2+m2$Var), log(m3$Bias2+m3$Var), 
            outline = F, col = "grey", 
            main = sigv, #expression(paste(sigma, " = ", sigv)), 
            ylim = range(c(log(m1$Bias2+m2$Var), log(m2$Bias2+m2$Var), log(m3$Bias2+m3$Var))), #cex.main = 1.5,
            ylab = expression(paste("log MSE","(",hat(m),")")))}
  
}
