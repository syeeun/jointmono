library(fda)
library(cvTools)

source("functions/fn_ftnlz.R")
source("functions/fn_cvfit_mfpca.R")
source("functions/fn_fit_mfpca.R")
source("functions/fn_ortho.R")
source("functions/fn_fig_mfpca_wind.R")

DAT = DAT0 = readRDS("application_windpower/windpower_data.rds")

## relative power production (0~1)*100 (%)
DAT0$yobs = (DAT$yobs)/max(DAT$yobs) 

## number of curves
M = max(DAT$id)
## extended grids to draw curves smoothly at boundaries
ttemp = seq(3, 13, l=100) 
## basis and integrated basis functions to draw curves
obsp = o.bsp.mat(ttemp, 10, 5, F)
oBsp = o.bsp.mat(ttemp, 10, 5, T)
## wind speed range
vlim = c(4,12)

#### Set initial parameters for each number of pc functions #### 
ini.par3 = list(theta = matrix(0, 10, 1), Theta = matrix(1:6/10, 10, 3),
                ascore = matrix(0, M, 3), bcoef = matrix(c(0,1), M, 2, byrow=T))

#### To determine the number of PC functions ####
### (Note: The inputs 'c_lmu' and 'c_lpc' are grids to run the cross-validation 
###        so that they can be changed for practical uses. 
###        We started with four PC functions and decrease to two. 
###        Each run takes a while to give an optimal grid.) 

## npc = 3
cv_3pc = cv.search(c_lmu = 10^{-(4:0)}, c_lpc = 10^{6:9}, K = 5, trng = range(ttemp),
                   dat = DAT0, ini = ini.par3, elim = 0.1, maxiter = 30)
round(cv_3pc$opt_cv_res, 5); cv_3pc$opt_lambda; # plot_cv_res(cv_3pc$opt_cv_res)
obj_mfpca_3 = fit_mfpca(DAT0, ini.par = ini.par3, trng = range(ttemp), norder = 5,
                        lmu = 10^{-4}, lpc = 10^{9}, elim = 0.05, maxiter = 20, printit = T)

## (a) scree_box.png
boxplot(obj_mfpca_3$est$ascore, xaxt="n", 
        xlab = "principal component functions", ylab = "principal component scores", col = "grey", pch = 16, cex = .6)
axis(1, at = 1:ncol(obj_mfpca_3$est$ascore), labels = c("1st", "2nd", "3rd"))
## (b) scree.png
plot(1:3, apply(obj_mfpca_3$est$ascore, 2, var), type = "o", xaxt = "n",
     xlab = "principal component function", ylab = "variance of principal component scores", pch = 4)
axis(1, at = c(1:3), labels = c("1st", "2nd", "3rd"))

#### Check convergence ####
obj_mfpca = obj_mfpca_3
track_fitted_mfpca(obj_mfpca, "mu")
track_fitted_mfpca(obj_mfpca, "fpc")
track_fitted_mfpca(obj_mfpca, "alpha")
track_fitted_mfpca(obj_mfpca, "beta")

#### Fitted mean (m) and relative curvature (w) curves ####
tgridtemp = seq(vlim[1]+.2, vlim[2]-.2, l = 1000)
## (a) mhat.png
plot(tgridtemp, Hfn(tgridtemp, with(obj_mfpca$est, theta), oBsp), 
     ylab = "m", cex.lab = 1.3, lwd = 2,
     xlab = "wind speed (m/s)", main = "", type = "l", xlim = vlim, ylim = c(0, 1))
abline(v=vlim, h = c(0,1), lty = 3)
## (b) what.png
plot(tgridtemp, Wfn(tgridtemp, with(obj_mfpca$est, theta), obsp), 
     ylab = "w", cex.lab = 1.3,lwd = 2,
     xlab = "wind speed (m/s)", main = "", type = "l", xlim = vlim, ylim = c(-1, 1))
abline(v=vlim, h = 0, lty = 3)

#### Fitted principal component functions ####
## (a) windpc1.png
tgridtemp = seq(vlim[1], vlim[2],l=100)
par(mar = c(5, 5, 4, 1) + 0.1)
plot(tgridtemp, Wfn(tgridtemp, with(obj_mfpca$est, Theta[,1]), obsp), 
     ylab = expression(f[1]), cex.lab = 1.3, lwd = 2, 
     xlab = "wind speed (m/s)", main = "", type = "l", xlim = c(4, 12), ylim = c(-0.6, 0.8))
abline(v = vlim, h = 0, lty = 3)
## (b) windpc2.png
plot(tgridtemp, Wfn(tgridtemp, with(obj_mfpca$est, Theta[,2]), obsp), 
     ylab = expression(f[2]), cex.lab = 1.3, lwd = 2,
     xlab = "wind speed (m/s)", main = "", type = "l", xlim = c(4, 12), ylim = c(-0.6, 0.8))
abline(v = vlim, h = 0, lty = 3)
## (b) windpc3.png
plot(tgridtemp, Wfn(tgridtemp, with(obj_mfpca$est, Theta[,3]), obsp), 
     ylab = expression(f[3]), cex.lab = 1.3, lwd = 2,
     xlab = "wind speed (m/s)", main = "", type = "l", xlim = c(4, 12), ylim = c(-0.6, 0.8))
abline(v = vlim, h = 0, lty = 3)


