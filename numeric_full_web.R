library(nleqslv)

thetas <- runif(30,1,1000)
Drs <- 0.2*log(thetas)- min(0.2*log(thetas))
R0 <- 100
tau <- 0.1

dynamic_full_web(thetas,Drs,R0,tau)

dynamic_full_web <- function(thetas,Drs,R0,tau){
	Rs <- c(rep(0,length(thetas)),R0)
	newRs <- one_step(Rs,thetas,Drs,tau)
	t <- 1
	while(sum((Rs-newRs)^2)>R0*10^-5){
		#if(!is.na(match(t,tseq))){
			#Y <- hist(newRs)
			#lines(Y$)
		#}
		Rs <- newRs
		newRs <- one_step(Rs,thetas,Drs,tau)
		t <- t+1
	}
	links <- matrix(unlist(lapply(newRs,one_node,thetas=thetas,Drs=Drs)),ncol=length(newRs))
	
	list(Rs = newRs,links=links)
}

one_step <- function(initRs,thetas,Drs,tau){
	
	links <- matrix(unlist(lapply(initRs,one_node,thetas=thetas,Drs=Drs)),ncol=length(initRs))
	newRs <- c()
	for(i in 1:length(thetas)){
		newRs[i] <- initRs[i]+ sum(links[i,])*tau - sum(links[,i]) ## input minus output; since initRs[i] = sum(links[,i]), equal to sum(links[i,])*tau 
	}
	c(newRs,initRs[length(initRs)])
}

one_node <- function(R,thetas,Drs){
	if(R==0){
		return(rep(0,length(thetas)))
	}
	init.C <- (R/sum(thetas^(Drs/(Drs-1))))^(mean(Drs)-1)
	C <- nleqslv(init.C,single_constraint,R=R,thetas=thetas,Drs=Drs)
	if(C[[3]]!=1){
		print(C[[4]])
	}
	one_allocation(C[[1]],thetas,Drs)
}

single_constraint <- function(C,R,thetas,Drs){
	sum(one_allocation(C,thetas,Drs)) - R
}

one_allocation <- function(C,thetas,Drs){
	C^(1/(Drs-1))*thetas^(Drs/(Drs-1))
}