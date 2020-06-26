library(fda)
library(MASS)

source("functions/fn_data_simul.R")
source("functions/fn_ftnlz.R")
source("functions/fn_ortho.R")
source("functions/fn_mse.R")

m1_01 = readRDS("simulation/simresult_compare//m1_01.rds")
m1_05 = readRDS("simulation/simresult_compare//m1_05.rds")
m1_10 = readRDS("simulation/simresult_compare//m1_10.rds")

m2_01 = readRDS("simulation/simresult_compare//m2_01.rds")
m2_05 = readRDS("simulation/simresult_compare//m2_05.rds")
m2_10 = readRDS("simulation/simresult_compare//m2_10.rds")

m3_01 = readRDS("simulation/simresult_compare//m3_01.rds")
m3_05 = readRDS("simulation/simresult_compare//m3_05.rds")
m3_10 = readRDS("simulation/simresult_compare//m3_10.rds")

true_par = list()
for(j in 1:100){
true_par[[j]] = fn_data_simul(50, seedn = 1234 + j - 1)$par
}

mse_m1_01 = f.mse.y_ram(m1_01, true_par)
mse_m1_05 = f.mse.y_ram(m1_05, true_par)
mse_m1_10 = f.mse.y_ram(m1_10, true_par)

mse_m2_01 = f.mse.y(m2_01, true_par)
mse_m2_05 = f.mse.y(m2_05, true_par)
mse_m2_10 = f.mse.y(m2_10, true_par)

mse_m3_01 = f.mse.y(m3_01, true_par)
mse_m3_05 = f.mse.y(m3_05, true_par)
mse_m3_10 = f.mse.y(m3_10, true_par)

msew_m1_01 = f.mse.w_ram(m1_01, true_par)
msew_m1_05 = f.mse.w_ram(m1_05, true_par)
msew_m1_10 = f.mse.w_ram(m1_10, true_par)

msew_m2_01 = f.mse.w(m2_01, true_par)
msew_m2_05 = f.mse.w(m2_05, true_par)
msew_m2_10 = f.mse.w(m2_10, true_par)

msew_m3_01 = f.mse.w(m3_01, true_par)
msew_m3_05 = f.mse.w(m3_05, true_par)
msew_m3_10 = f.mse.w(m3_10, true_par)

#### Figure 5 - top panel: simbox_w.png ####
par(mfrow = c(1,3), family = "", mar = c(5, 5, 4, 2)+0.1)
boxplot(log(msew_m1_01), log(msew_m2_01), log(msew_m3_01), ylim = c(-1,4), outline = F,
        main = expression(paste(sigma, " = ", 0.01)), cex.main = 1.5,
        names = c("Ramsay", "Two-step", "Integrated"), cex.axis = 1.2, 
        ylab = expression(paste("log MISE(",hat(w),")")), cex.lab = 1.5, col = "grey")
boxplot(log(msew_m1_05), log(msew_m2_05), log(msew_m3_05), ylim = c(-1,4), outline = F,  
        main = expression(paste(sigma, " = ", 0.05)), cex.main = 1.5,
        names = c("Ramsay", "Two-step", "Integrated"), cex.axis = 1.2, col = "grey")
boxplot(log(msew_m1_10), log(msew_m2_10), log(msew_m3_10), ylim = c(-1,4), outline = F, 
        main = expression(paste(sigma, " = ", 0.1)), cex.main = 1.5,
        names = c("Ramsay", "Two-step", "Integrated"), cex.axis = 1.2, col = "grey")

#### Figure 5 - bottom panel: simbox_y.png ####
par(mfrow = c(1,3), family = "", mar = c(5, 5, 4, 2)+0.1)
boxplot(log(mse_m1_01), log(mse_m2_01), log(mse_m3_01), ylim = c(-16,1), outline = F, 
        main = expression(paste(sigma, " = ", 0.01)), cex.main = 1.5,
        names = c("Ramsay", "Two-step", "Integrated"), cex.axis = 1.2, 
        ylab = expression(paste("log MISE(",hat(m),")")), cex.lab = 1.5, col = "grey")
boxplot(log(mse_m1_05), log(mse_m2_05), log(mse_m3_05), ylim = c(-16,1), outline = F, 
        main = expression(paste(sigma, " = ", 0.05)), cex.main = 1.5,
        names = c("Ramsay", "Two-step", "Integrated"), cex.axis = 1.2, col = "grey")
boxplot(log(mse_m1_10), log(mse_m2_10), log(mse_m3_10), ylim = c(-16,1), outline = F, 
        main = expression(paste(sigma, " = ", 0.1)), cex.main = 1.5,
        names = c("Ramsay", "Two-step", "Integrated"), cex.axis = 1.2, col = "grey")

