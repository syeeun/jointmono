source("fn_ftnlz.R")
source("fn_ortho.R")

source("fn_fit_mfpca.R")
source("fn_fig_mfpca.R")

# library(fda)
# library(MASS)
#dat.sim = read.table("data-simul.text", header=T)
sim50 = fn_data_simul(50)
dat.sim = sim50$dat
true.par = sim50$par

#### INITIAL = TRUE ####
true.par = readRDS(file = "truepar.rds")
true.par$Theta = true.par$Theta#[,1,drop = F]
true.par$ascore = true.par$ascore#[,1,drop = F]

obj_mfpca0 = fit_mfpca(dat.sim, true.par, norder = 5, trng = c(0,1),
                       lmu = 5, lpc = 1, elim = 0.05, maxiter = 20)

#### INITIAL = from Ramsay + James ####
ini.par = readRDS(file = "inipar.rds")
ini.par$Theta = ini.par$Theta#[,1,drop = F]
ini.par$ascore = ini.par$ascore#[,1,drop = F]

obj_mfpca1 = fit_mfpca(dat.sim, ini.par, norder = 5, trng = c(0,1),
                       lmu = 10, lpc = 10, elim = 0.1, maxiter = 20)

#### INITIAL = arbitrary ####
ini.par2 = list(theta = matrix(1, 10, 1), Theta = matrix(1, 10, 2),
                ascore = matrix(0, 50, 2), bcoef = matrix(c(0,1), 50, 2, byrow=T))

obj_mfpca2 = fit_mfpca(dat.sim, ini.par2, norder = 5,
                       lmu = 10, lpc = 10, elim = 0.05, maxiter = 20, printit = T)


#### CHECK ####
obj_mfpca = obj_mfpca2

track_fitted_mfpca(obj_mfpca, "mu")
track_fitted_mfpca(obj_mfpca, "fpc")
track_fitted_mfpca(obj_mfpca, "alpha")
track_fitted_mfpca(obj_mfpca, "beta")

plot_fitted_mfpca(obj_mfpca, "mu", true.par = true.par); abline(v=c(0,1), lty=3)
plot_fitted_mfpca(obj_mfpca, "fpc", true.par = true.par); abline(v=c(0,1), lty=3)
plot_fitted_mfpca(obj_mfpca, "y", true.par = true.par); abline(v=c(0,1), lty=3)
plot_fitted_mfpca(obj_mfpca, "w", true.par = true.par); abline(v=c(0,1), lty=3)

plot_est_true_ini(obj_mfpca0, true.par, true.par)
plot_est_true_ini(obj_mfpca1, true.par, ini.par)
plot_est_true_ini(obj_mfpca2, true.par, ini.par2)






tempfit_ram = with(dat.sim, fit_ramsay(tobs[id==1], yobs[id==1], ini.theta = matrix(1,10,1), trng = c(0,1), 
                                       norder = 5, lmu = 0.005, elim = 0.01, maxiter = 10))
track_fitted_ramsay(tempfit_ram, "mu")
track_fitted_ramsay(tempfit_ram, "beta")
plot_fitted_ramsay(tempfit_ram, "y")

