#Calculate link width for a fully connected food web

tau_u <- 0.1

FW_full <- function(nsp,tdistr,Dr){##small disperse meaning linear, disperse = 0 meaning flat line
	Y <- try(try_FW_full(nsp,tdistr,Dr))
	while(class(Y)=="try-error"){
		Y <- try(try_FW_full(nsp,tdistr,Dr))	
	}
	Y
}

try_FW_full <- function(nsp,tdistr,Dr){
	community <- generate.pool(nsp,tdistr,Dr) 
	community0 <- community
	equilibrium <- FALSE
	step <- 0
	while(!equilibrium&step<100){
		community <- one_step_full(community)
		equilibrium <- judge_equi_full(community,community0)
		community0 <- community
		step <- step +1
		hist(community$sp_mat$BPN,xlab="abundance")
	}
	print(paste("step=",step))
	community <- get_after_abundance(community)#APN should be exactly half of BPN
	community <- flow_pattern(community)
	community
}

Trim_web <- function(link_mat,L){
	threshold <- sort(as.numeric(link_mat),decreasing=T)[L]
	link_mat[link_mat<threshold]<- 0 
	link_mat
}

one_step_full <- function(community){
	indice <- community$sp_mat$index
	for(index in indice){
		community <- all_link(community,index)
	}
	community
}

all_link <- function(community,index){
	preys <- community$sp_mat$index
	##Everything can prey on everything (including itself)
	SSN <- sum(unlist(lapply(preys,try_prey,community=community,index=index,tau_c=1)))
	community$link_mat[,index] <- 1
	community$sp_mat$BPN[index] <- SSN
	community
}

try_prey <- function(prey,community,index,tau_c){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	predators <- sp_mat$index[link_mat[prey,]>0]
	predators <- predators[predators!=index]
	if(length(predators)==0){##no other predators
		SSN <- sp_mat$BPN[prey]*sp_mat$theta[prey]/2*tau_c*tau_u/sp_mat$theta[index]
	}else{
		theta_seq <- c(sp_mat$theta[index]/tau_c,(sp_mat$theta/sp_mat$tau_c)[predators])
		Dr_seq <- sp_mat$Dr[c(index,predators)]
		R_prey <- sp_mat$theta[prey]*sp_mat$BPN[prey]/2*tau_u
		SSN <-solve.analytical(theta_seq,Dr_seq,R_prey,1)[1]
	}
	SSN
}


judge_equi_full <- function(community,community0,threshold=10^-5){
	N_seq <- community$sp_mat$BPN
	N0_seq <- community0$sp_mat$BPN
	C = sd(N_seq-N0_seq)/mean(N0_seq)<threshold
	C
}

flow_pattern <- function(community){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	for(sp in sp_mat$index[1:dim(link_mat)[2]]){
		preys <- sp_mat$index[link_mat[,sp]>0]
		link_mat[preys,sp] <- unlist(lapply(preys,try_prey,community=community,index=sp,tau_c=sp_mat$tau_c[sp]))*sp_mat$theta[sp]/tau_u
	}
	list(sp_mat=sp_mat,link_mat=link_mat)
}

pool.theta <- function(nsp,mu,distr){
	if(distr=="poisson"){
		theta_seq <- rpois(nsp,mu)+1
	}else if(distr=="geometric"){
		theta_seq <- rgeom(nsp,1/(mu+1))+1
	}else if(distr=="uniform"){
		theta_seq <- runif(nsp,0,2*mu)
	}else{
		print("no such bs distribution")
	}
	theta_seq
}

generate.pool <- function(nsp,tdistr,Dr,init="MERA"){
	sp_pool <- matrix(nrow=nsp,ncol=7)
	sp_pool[,1] <- 1:nsp
	sp_pool[,2] <- pool.theta(nsp,1,tdistr)
	sp_pool[,3] <- rep(Dr,nsp)
	sp_pool[,4] <- rep(NA,nsp)
	if(init=="MERA"){
		sp_pool[,5] <- solve.analytical(sp_pool[,2],sp_pool[,3],100,1)
	}else if(init=="uniform"){
		sp_pool[,5] <- runif(nsp,1,100)
	}else if(init=="poisson"){
		sp_pool[,5] <- rpois(nsp,50)
	}else if(init=="geometric"){
		sp_pool[,5] <- rgeom(nsp,0.05)+1
	}else{
		print("no such initial state defined")
		return()
	}
	sp_pool[,6] <- rep(0,nsp)
	sp_pool[,7] <- 1
	#sp_pool <- sp_pool[order(sp_pool[,5]),]##order competitiveness (from lowest to highest)
	colnames(sp_pool) <- c("index","theta","Dr","gamma","BPN","APN","tau_c")
	link_mat <- matrix(0,nrow=nsp,ncol=nsp)
	colnames(link_mat) <- rownames(link_mat) <- paste("Sp",1:nsp)
	list(sp_mat=data.frame(sp_pool),link_mat=link_mat)
}

update_links <- function(community,selection){
	community$link_mat[selection$added_prey,selection$index] <- 1
	if(length(selection$added_prey)>0){
		community$link_mat[-selection$added_prey,selection$index] <- 0
	}else{
		community$link_mat[,selection$index] <- 0
	}
	community$sp_mat$tau_c[selection$index] <- selection$tau_c
	community
}

get_after_abundance <- function(community){
	for(sp in community$sp_mat$index){
		community$sp_mat$APN[sp] <- SSN_after_predation(sp,community)
	}
	community
}

SSN_after_predation <- function(sp,community){##half the prey population and add the uncaptured prey back to its steady state abundance
	sp_mat = community$sp_mat
	link_mat = community$link_mat
	SSN <- sp_mat$BPN[sp]
	predators <- sp_mat$index[1:dim(link_mat)[2]][link_mat[sp,]>0]
	if(length(predators)>0){
		SSN <- 0.5*SSN
		for(i in 1:length(predators)){
			tau_c <- sp_mat$tau_c[predators[i]]
			SSN <- SSN+try_prey(sp,community,predators[i],tau_c)*sp_mat$theta[predators[i]]/sp_mat$theta[sp]/tau_u*(1-tau_c)/tau_c##the uncaptured amount from a given predator
		}
	}
	SSN
}
