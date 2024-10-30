rm(list=ls())
# plot spectra of compound 23
pdf("CM2_compound23.pdf")
for(i in 1:10){
     file_name = paste0("23-",formatC(i,width=2,flag=0),".csv");
     data = read.csv(file_name,skip=2)
     mz = as.numeric(data[,2])
     ab = as.numeric(data[,3])
     ab = 100*ab/max(ab);
     plot(mz,ab,type="h",xlim=c(0,450),ylim=c(0,100),xlab = "m/z", ylab="rel.intensity",main=file_name)
}
dev.off()
