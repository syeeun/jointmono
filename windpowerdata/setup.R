dat.real = read.table("EXP.txt", header = T)
temp_day = as.factor(substr(dat.real$time, 1, 7))

idrm = ((dat.real$V > 12) | (dat.real$V < 4) | (as.integer(temp_day) %in% which(table(temp_day) < 10)))
dat.real = dat.real[!idrm,]

dat = with(dat.real, data.frame(trt = trt, 
                                mon = substr(time, 1, 7), day = substr(time, 1, 10), 
                                P = Ptest, V = V))
dat$id = as.numeric(as.factor(dat$day)) %/% 7 + 1

order_id = tapply(dat$V, dat$id, order)
DAT = NULL
for(k in 1:length(order_id)){
  DAT = rbind(DAT, dat[dat$id == k,][order_id[[k]],])
}

DAT = DAT[,c("trt", "id", "P", "V")]; colnames(DAT) = c("trt", "id", "yobs", "tobs" )

saveRDS(DAT, "windpower_data.rds")
