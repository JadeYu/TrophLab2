source("step_optimization.R")

nsp = 50
tdistr="geometric"
Dr = 0.9
gamma = 0.001 ##from 0.001 to 0.02


FW = one_food_web(nsp,tdistr,Dr,gamma)
plot_foodweb(FW)
chains <- get_chains(FW)

##Write the back and forth behavior in the analysis!!
#num_basal(FW) + num_predators(FW)
num_basal(FW)
num_top(FW)
log(length(chains))
connectance(FW)
Chain_length(chains,loop=T)
loops(chains)

sum(FW$sp_mat$theta*FW$sp_mat$APN)

Dr_seq <- runif(50,0.1,0.9)
lgamma_seq <- runif(50,-2,0)


Result = Analytical_LS(Dr_seq,lgamma_seq,nsp=30,tdistr="uniform",nrep=1,graph=NA)

Analytical_LS <- function(Dr_seq,lgamma_seq,nsp=30,tdistr="uniform",nrep=1,graph=NA){
	L <- S <- matrix(nrow=nrep,ncol=length(Dr_seq))
	for(t in 1:length(Dr_seq)){
		for(r in 1:nrep){
				FW <- one_food_web(nsp,tdistr,Dr=Dr_seq[t],gamma=10^lgamma_seq[t],disperse=c(0,0))
			L[r,t] <- get_links(FW)
			S[r,t] <- basal(FW) +  predators(FW)
		}
		print(paste("t =",t))	
	}
	mL = L/S
	C = L/S^2
	pmL = (1-Dr_seq)/Dr_seq/10^lgamma_seq
	zoom = 1
	if(!is.na(graph)){
		png(paste("new_graphs/",graph,sep=""),height=800,width=800)
		par(mfrow=c(2,2))
		zoom=2
	}
	plot(as.numeric(L)~rep(Dr_seq,each=nrep),xlab="Dr",ylab="L",cex=zoom,cex.lab=zoom,cex.axis=zoom)
	plot(as.numeric(L)~rep(10^lgamma_seq,each=nrep),xlab="gamma",ylab="L",cex=zoom,cex.lab=zoom,cex.axis=zoom)
	
	plot(as.numeric(S)~rep(Dr_seq,each=nrep),xlab="Dr",ylab="S",cex=zoom,cex.lab=zoom,cex.axis=zoom)
	plot(as.numeric(S)~rep(10^lgamma_seq,each=nrep),xlab="gamma",ylab="S",cex=zoom,cex.lab=zoom,cex.axis=zoom)
	

plot(as.numeric(mL)~rep(pmL,each=nrep),xlab="(1-Dr)/(gamma Dr)",ylab="L/S",cex=zoom,cex.lab=zoom,cex.axis=zoom)

plot(as.numeric(S)~rep(pmL,each=nrep),xlab="(1-Dr)/(gamma Dr)",ylab="S",cex=zoom,cex.lab=zoom,cex.axis=zoom)

plot(as.numeric(C)~rep(pmL,each=nrep),xlab="(1-Dr)/(gamma Dr)",ylab="connectance",cex=zoom,cex.lab=zoom,cex.axis=zoom)
	

if(!is.na(graph)){
	dev.off()
}
	list(L=L,S=S)
}


