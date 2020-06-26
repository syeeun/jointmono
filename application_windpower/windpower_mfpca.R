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
ini.par2 = list(theta = matrix(0, 10, 1), Theta = matrix(1:4/10, 10, 2),
                ascore = matrix(0, M, 2), bcoef = matrix(c(0,1), M, 2, byrow=T))

ini.par3 = list(theta = matrix(0, 10, 1), Theta = matrix(1:6/10, 10, 3),
                ascore = matrix(0, M, 3), bcoef = matrix(c(0,1), M, 2, byrow=T))

ini.par4 = list(theta = matrix(0, 10, 1), Theta = matrix(1:7/10, 10, 4),
                ascore = matrix(0, M, 4), bcoef = matrix(c(0,1), M, 2, byrow=T))


#### To determine the number of PC functions ####
### (Note: The inputs 'c_lmu' and 'c_lpc' are grids to run the cross-validation 
###        so that they can be changed for practical uses. 
###        We started with four PC functions and decrease to two. 
###        Each run takes a while to give an optimal grid.) 

## npc = 4
cv_4pc = cv.search(c_lmu = 10^{-(4:2)}, c_lpc = 10^{4:6}, K = 5, trng = range(ttemp),
                   dat = DAT0, ini = ini.par4, elim = 0.1, maxiter = 30)
round(cv_4pc$opt_cv_res, 5); cv_4pc$opt_lambda; # plot_cv_res(cv_4pc$opt_cv_res)
obj_mfpca_4 = fit_mfpca(DAT0, ini.par = ini.par4, trng = range(ttemp), norder = 5,
                        lmu = 10^{-3}, lpc = 10^{5}, elim = 0.1, maxiter = 20, printit = T)
## npc = 3
cv_3pc = cv.search(c_lmu = 10^{-(4:0)}, c_lpc = 10^{6:9}, K = 5, trng = range(ttemp),
                   dat = DAT0, ini = ini.par3, elim = 0.1, maxiter = 30)
round(cv_3pc$opt_cv_res, 5); cv_3pc$opt_lambda; # plot_cv_res(cv_3pc$opt_cv_res)
obj_mfpca_3 = fit_mfpca(DAT0, ini.par = ini.par3, trng = range(ttemp), norder = 5,
                        lmu = 1, lpc = 10^{9}, elim = 0.05, maxiter = 20, printit = T)
## npc = 2
cv_2pc = cv.search(c_lmu = 10^{-(6:4)}, c_lpc = 10^{4:6}, K = 5, trng = range(ttemp),
                   dat = DAT0, ini = ini.par2, elim = 0.1, maxiter = 30)
round(cv_2pc$opt_cv_res, 5); cv_2pc$opt_lambda; # plot_cv_res(cv_2pc$opt_cv_res)
obj_mfpca_2 = fit_mfpca(DAT0, ini.par = ini.par2, trng = range(ttemp), norder = 5,
                        lmu = 0.001, lpc = 10^{5}, elim = 0.05, maxiter = 20, printit = T)

#### PC scores for npc = 4 ####
### Figure 6 in the paper
## (a) scree_box.png
boxplot(obj_mfpca_4$est$ascore, xaxt="n", 
        xlab = "principal component functions", ylab = "principal component scores", col = "grey", pch = 16, cex = .6)
axis(1, at = 1:ncol(obj_mfpca_4$est$ascore), labels = c("1st", "2nd", "3rd", "4th"))
## (b) scree.png
plot(1:4, apply(obj_mfpca_4$est$ascore, 2, var), type = "o", xaxt = "n",
     xlab = "principal component function", ylab = "variance of principal component scores", pch = 4)
axis(1, at = c(1:4), labels = c("1st", "2nd", "3rd", "4th"))

#### Check convergence ####
obj_mfpca = obj_mfpca_3 #obj_mfpca_2 #obj_mfpca_4
track_fitted_mfpca(obj_mfpca, "mu")
track_fitted_mfpca(obj_mfpca, "fpc")
track_fitted_mfpca(obj_mfpca, "alpha")
track_fitted_mfpca(obj_mfpca, "beta")

#### Fitted mean (m) and relative curvature (w) curves ####
### Figure 7
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
### Figure 8
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


##### Examples of fitted individual curves (y) ####
### Figure 9
par(mar = c(5,4.5,4,1) + 0.1)
## (a) a1pos.png; Steepy start for poistive high a1 
plot.ym(obj_mfpca, 31, trng = c(4.2,11.8), yrnginput = c(0,1), cex.lab = 1.3,
        xlab = "wind speed (m/s)", ylab = "relative power curve")
## (b) a1neg.png; relatively low start for negative high a1
plot.ym(obj_mfpca, 70, trng = c(4.2,11.8), yrnginput = c(0,1), cex.lab = 1.3, 
        xlab = "wind speed (m/s)", ylab = "relative power curve")
## (c) a2pos.png; greater curvature before inflection for positive high a2
plot.ym(obj_mfpca, 63, trng = c(4.2,11.8), yrnginput = c(0,1), cex.lab = 1.3, 
        xlab = "wind speed (m/s)", ylab = "relative power curve")
## (d) a2neg.png; greater curvature after inflection for negative high a2 
plot.ym(obj_mfpca, 3, trng = c(4.2,11.8), yrnginput = c(0,1), cex.lab = 1.3, 
        xlab = "wind speed (m/s)", ylab = "relative power curve")



