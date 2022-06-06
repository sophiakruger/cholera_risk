## Running eight cores
n_replicates <- 8
n_iters <- 500

## For loop range for each core
ranges <- cbind(round(seq(1,n_iters, n_iters/n_replicates)), round(seq(n_iters/n_replicates, n_iters,  n_iters/n_replicates)))
ranges <- matrix(c(c(1,62), c(63, 125), c(126, 188), c(189, 250), c(251, 312), c(313, 375), c(376,438), c(439, 500)), nrow=8, byrow=T)

## RF iterations script
r_script_name <- "rf_iters.R"

## Create folder for each replicate (core)
for (i in 1:n_replicates)
	{
	system(paste("mkdir Replicate", i,sep=""))
	setwd(paste("Replicate",i,sep=""))
	write.table(ranges[i,], "local_ranges.txt",col.names=F,row.names=F)
	system(paste("cp ../", r_script_name, " .", sep=""))
	setwd("..")
	}

## Run replicates
system(paste("cd Replicate1; Rscript ", r_script_name, " & ", paste(paste(paste("cd ../Replicate", 2:n_replicates, "; Rscript ", r_script_name, sep=""),collapse=" & "), "&& fg")))
