f.mse.y_ram = function(m_est, true){
  
  RAND = length(m_est); ncurve = nrow(true[[1]]$bcoef)
  tgrid = seq(0,1,l=100); delta = diff(tgrid)[1]
  oBsp = o.bsp.mat(tgrid, 10, 5, T)
  
  mse = matrix(NA, RAND, ncurve)
  for(j in 1:RAND){
    for(i in 1:ncurve){
      yhat = with(m_est, beta[[j]][i,1] + beta[[j]][i,2]*Hfn(tgrid, theta[[j]][i,], oBsp))
      y0 = with(true[[j]], bcoef[i,1] + bcoef[i,2]*Hfn(tgrid, theta, oBsp))
      mse[j, i] = sum((y0-yhat)^2)*delta
    }
  }
  
  return(c(mse))
}

f.mse.y = function(m_est, true){

    RAND = length(m_est); ncurve = nrow(true[[1]]$bcoef)
    tgrid = seq(0,1,l=100); delta = diff(tgrid)[1]
    oBsp = o.bsp.mat(tgrid, 10, 5, T)
    
    mse = matrix(NA, RAND, ncurve)
    for(j in 1:RAND){
      for(i in 1:ncurve){
        yhat = with(m_est, beta[[j]][i,1] + beta[[j]][i,2]*Hfn(tgrid, theta[[j]] + Theta[[j]] %*% ascore[[j]][i,], oBsp))
        y0 = with(true[[j]], bcoef[i,1] + bcoef[i,2]*Hfn(tgrid, theta + Theta %*% ascore[i,], oBsp))
        mse[j, i] = sum((y0-yhat)^2)*delta
        }
    }
    
  return(c(mse))
}


f.mse.w_ram = function(m_est, true){
  
  RAND = length(m_est); ncurve = nrow(true[[1]]$bcoef)
  
  mse = matrix(NA, RAND, ncurve)
  for(j in 1:RAND){
    for(i in 1:ncurve){
      what = with(m_est, theta[[j]][i,])
      w0 = with(true[[j]], theta)
      mse[j, i] = as.numeric(crossprod(what-w0))
    }
  }
  
  return(c(mse))
}

f.mse.w = function(m_est, true){
  
  RAND = length(m_est); ncurve = nrow(true[[1]]$bcoef)

  mse = matrix(NA, RAND, ncurve)
  for(j in 1:RAND){
    for(i in 1:ncurve){
      what = with(m_est, theta[[j]] + Theta[[j]] %*% ascore[[j]][i,])
      w0 = with(true[[j]], theta + Theta %*% ascore[i,])
      mse[j, i] = as.numeric(crossprod(what-w0))
    }
  }
  
  return(c(mse))
}

