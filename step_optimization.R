source("MERA_SSN.R")
source("anal_metrics.R")
source("Visualization.R")

##Global variables
tau_u = 0.1

##Rank species by competitiveness and start from the lowest (most likely to be basal)

##functions
one_food_web <- function(nsp,tdistr,Dr,gamma){##small disperse meaning linear, disperse = 0 meaning flat line
	Y <- try(try_one_food_web(nsp,tdistr,Dr,gamma))
	while(class(Y)=="try-error"){
		Y <- try(try_one_food_web(nsp,tdistr,Dr,gamma))	
	}
	Y
}

try_one_food_web <- function(nsp,tdistr,Dr,gamma){
	community <- generate.pool(nsp,tdistr,Dr,gamma) 
	scenario <-FWA_switch_cost(community,tau_u)
##omni specifies whether the species can consume the fundamental resource or not
	scenario
}

##functions################
FWA_switch_cost <- function(community,tau_u,threshold=10^-3){
	community0 <- community
	equilibrium <- FALSE
	step <- 0
	while(!equilibrium&step<100){
		community <- one_step(community,tau_u)
		equilibrium <- judge_equi(community,community0)
		community0 <- community
		step <- step +1
		print(paste("S=",num_basal(community) + num_predators(community),"L=",num_links(community)))
	}
	community <- get_after_abundance(community,tau_u)
	community <- flow_pattern(community,tau_u)
	community
}

one_step <- function(community,tau_u){
	indice <- community$sp_mat$index
	for(index in indice){
		predation <- select_prey(community,index,tau_u)
		basal <- basal_SSN(index,community)
		if(basal$SSN>=predation$SSN){
			selection <- basal
		}else{
			selection <- predation
		}
		#new_comm <- update_links(new_comm,selection)
		community <- update(community,selection)
	}
	#new_comm <- update_abundance(new_comm)
	community
}

basal_SSN <- function(index,community){
	indice <- community$sp_mat$index
	basals <- indice[is.basal(indice,community)]
	species <- c(index,basals[basals!=index])
	basal_mat <- community$sp_mat[species,]
	SSN <- solve.analytical(basal_mat$theta,basal_mat$Dr,100,1)[1]
	list(index=index,SSN=SSN,tau_c=1,added_prey=c())
}

select_prey <- function(community,index,tau_u){
	sp_mat <- community$sp_mat
	tau_c <- 1
	increase <- SSN <- 0
	added_prey <- c()
	predators <- sp_mat$index[community$link_mat[index,]>0]
	preys <- sp_mat$index[-c(index,predators)]
	#preys <- community$sp_mat$index[-index]
	#preys <- community$sp_mat$index[-index]
	##Remove species that prey on it to avoid direct loops
	SSNs <- numeric(length(preys))
	if(length(preys)>0){
		for(p in 1:length(preys)){
			SSNs[p] <- try_prey(preys[p],community,index,tau_c*exp(-sp_mat$gamma[preys[p]]*(sp_mat$Dr[preys[p]])),tau_u)
		}
		preys <- preys[order(SSNs,decreasing=T)]##determine priority 
	}
	increase <- sum(unlist(lapply(c(added_prey,preys[1]),try_prey,tau_c=tau_c*exp(-sp_mat$gamma[preys[1]]*(sp_mat$Dr[preys[1]])),community=community,index=index,tau_u=tau_u))) - SSN
	while(increase>0&&length(preys)>0){
		added_prey <- c(added_prey,preys[1])
		tau_c <- tau_c*exp(-sp_mat$gamma[preys[1]]*(sp_mat$Dr[preys[1]]))
		preys <- preys[-1]##remove the best prey	
		SSN <- SSN+increase
		if(length(preys)>0){
			increase <- sum(unlist(lapply(c(added_prey,preys[1]),try_prey,tau_c=tau_c*exp(-sp_mat$gamma[preys[1]]*(sp_mat$Dr[preys[1]])),community=community,index=index,tau_u=tau_u))) - SSN
		}
	}
	list(index=index,SSN=SSN,tau_c=tau_c,added_prey=added_prey)
}

try_prey <- function(prey,community,index,tau_c,tau_u){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	predators <- sp_mat$index[1:dim(link_mat)[2]][link_mat[prey,]>0]
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

SSN_before_predation <- function(sp,community,tau_u){##calculate SSN before a species is preyed on
	sp_mat= community$sp_mat
	link_mat = community$link_mat
	preys <- sp_mat$index[link_mat[,sp]>0]
	SSN <- 0
	if(length(preys)>0){
		SSN <- sum(unlist(lapply(preys,try_prey,community=community,index=sp,tau_c=sp_mat$tau_c[sp],tau_u=tau_u)))
	}
	SSN
}

SSN_after_predation <- function(sp,community,tau_u){##half the prey population and add the uncaptured prey back to its steady state abundance
	sp_mat = community$sp_mat
	link_mat = community$link_mat
	SSN <- sp_mat$BPN[sp]
	predators <- sp_mat$index[1:dim(link_mat)[2]][link_mat[sp,]>0]
	if(length(predators)>0){
		SSN <- 0.5*SSN
		for(i in 1:length(predators)){
			tau_c <- sp_mat$tau_c[predators[i]]
			SSN <- SSN+try_prey(sp,community,predators[i],tau_c,tau_u)*sp_mat$theta[predators[i]]/sp_mat$theta[sp]/tau_u*(1-tau_c)/tau_c##the uncaptured amount from a given predator
		}
	}
	SSN
}

get_after_abundance <- function(community,tau_u){
	for(sp in community$sp_mat$index){
		community$sp_mat$APN[sp] <- SSN_after_predation(sp,community,tau_u)
	}
	community
}

judge_equi <- function(community,community0,threshold=0.05){
	C1 = abs(num_links(community)-num_links(community0)) <= threshold * num_links(community0)
	#C2 = (num_basal(community)+num_predators(community))== (num_basal(community0)+num_predators(community0))
	#C3 = sd(community$sp_mat$BPN - community0$sp_mat$BPN)/mean(community0$sp_mat$BPN)<threshold
	#print(paste("C1",C1,"C2",C2))
	C1
}

update <- function(community,selection){
	community$link_mat[selection$added_prey,selection$index] <- 1
	if(length(selection$added_prey)>0){
		community$link_mat[-selection$added_prey,selection$index] <- 0
	}else{
		community$link_mat[,selection$index] <- 0
	}
	community$sp_mat$tau_c[selection$index] <- selection$tau_c
	community$sp_mat$BPN[selection$index] <- selection$SSN
	community
}


flow_pattern <- function(community,tau_u){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	for(sp in sp_mat$index[1:dim(link_mat)[2]]){
		preys <- sp_mat$index[link_mat[,sp]>0]
		link_mat[preys,sp] <- unlist(lapply(preys,try_prey,community=community,index=sp,tau_c=sp_mat$tau_c[sp],tau_u=tau_u))*sp_mat$theta[sp]/tau_u
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

generate.pool <- function(nsp,tdistr,Dr,gamma){
	sp_pool <- matrix(nrow=nsp,ncol=7)
	sp_pool[,1] <- 1:nsp
	sp_pool[,2] <- pool.theta(nsp,1,tdistr)
	sp_pool[,3] <- rep(Dr,nsp)
	sp_pool[,4] <- rep(gamma,nsp)
	sp_pool[,5] <- solve.analytical(sp_pool[,2],sp_pool[,3],100,1)
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

update_abundance <- function(community){
	##first calculate for basal species
	indice <- community$sp_mat$index
	basals <- indice[is.basal(indice,community$link_mat)]
	basal_mat <- community$sp_mat[basals,]
	community$sp_mat$BPN[basals] <- solve.analytical(basal_mat$theta,basal_mat$Dr,100,1)
	non_basals <- indice[-basals]
	match <- FALSE
	BPN0 <- BPN <- community$sp_mat$BPN
	while(!match){
		for(sp in non_basals){
			BPN[sp] <- SSN_before_predation(sp,community,tau_u)
		}
		match <- sd(BPN-BPN0)/mean(BPN0)<10^-3
		BPN0 <- BPN
	}
	community$sp_mat$BPN <- BPN
	community
}



