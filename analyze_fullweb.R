rank_R <- function(R_seq){
	R_seq <- sort(R_seq)
	pR <- R_seq/sum(R_seq)
	cdf <- rank(R_seq)/length(R_seq)
	list(cdf=cdf,pR=pR)
}

plot_Rdistr <- function(R_seq,R_comp){
	rr <- rank_R(R_seq)
	rrc <- rank_R(R_comp)
	pRrange <- range(c(rr$pR,rrc$pR))
	plot(rr$cdf~rr$pR,ylab="relative rank",xlab="resource proportion (in R0)",xlim=pRrange,pch=16)
	points(rrc$pR,rrc$cdf,col=2,pch=16)
	legend("bottomright",c("with interaction","without interaction"),pch=16,col=1:2)
}

Trim_web <- function(link_mat,threshold,ctype="links"){
	if(ctype=="links"){
		threshold <- sort(as.numeric(link_mat),decreasing=T)[threshold]
	}else if(ctype=="MVFF"){
		threshold <- sum(as.numeric(link_mat))*threshold
	}else{
		print("constraint type invalid")
		return()
	}
	link_mat[link_mat<threshold]<- 0 
	link_mat
}

is.basal <- function(sps,link_mat){
	diag(link_mat) = 0
	colSums(link_mat[,sps])==0&rowSums(link_mat[sps,])!=0
}

is.top <- function(sps,link_mat){
	diag(link_mat) = 0 ##excluding cannibals
	colSums(link_mat[,sps])!=0&rowSums(link_mat[sps,])==0
}

get_basals <- function(link_mat){
	indice <- 1:dim(link_mat)[1]
	diag(link_mat) = 0
	indice[colSums(link_mat)==0&rowSums(link_mat)!=0]
}

num_basal <- function(link_mat){## does not include cannibals
	length(get_basals(link_mat))
}

num_predator <- function(link_mat){##species that are feeding on at least one other species
	sum(colSums(link_mat)>0)
}

num_top <- function(link_mat){##species that are not fed on by any other species
	diag(link_mat) <- 0
	sum(colSums(link_mat)>0&rowSums(link_mat)==0)
}

num_cann <- function(link_mat){##number of cannibal species
	sum(diag(link_mat) > 0)
}

get_chains <- function(link_mat){
	basals <- get_basals(link_mat)
	chains <- list()
	for(b in basals){
		chain <- b
		chains <- c(chains,part2all(chain,link_mat))
	}
	chains
}

part2all <- function(chain,link_mat){ #recursive function: for a given partial chain (basal to an intermediate species), get all possible complete chains containing it
	indice <- 1:dim(link_mat)[1]
	sp <- chain[length(chain)]
	predators <- indice[link_mat[sp,]>0]
	if(length(predators)==0){## when it's already complete
		all_chains <- list(chain)
	}else{
		all_chains <- list()
		for(p in predators){
			if(!is.na(match(p, chain))){## if the new species is already in the chain, chain should end here
				new_chains <- list(c(chain,p))
			}else{
				chain = c(chain,p)
				new_chains <- part2all(chain,link_mat)
			}
			all_chains <- c(all_chains,new_chains)
		}
	}
	all_chains
}

Chain_length <- function(chains,loop){
	if(!loop){
		chains = chains[unlist(lapply(chains,Loop_length))==0]
	}
	CL = unlist(lapply(chains,length))
	list(mCL= mean(CL),sCL= sd(CL))
}

Loop_length <- function(chain){
	Loop_sp <- c()
	if(length(chain)==length(unique(chain))){
		LL = 0
	}else{
		sps <- unique(chain)
		sp <- sort(sps)[table(chain)>1]
		
		indice <- (1:length(chain))[chain==sp]
		Loop_sp <- chain[indice[1]:indice[2]]
		LL = indice[2]-indice[1]+1
	}
	list(LL,Loop_sp)
}

loops <- function(chains){
	lsp <- c()
	LL_seq <- numeric(length(chains))
	for(i in 1:length(chains)){
		Loop <- Loop_length(chains[[i]])
		LL_seq[i] <- Loop[[1]]
		lsp <- c(lsp,Loop[[2]])
	}
	N_loops <- sum(LL_seq>0)
	MLL <- mean(LL_seq[LL_seq>0])
	#list(N_loops = N_loops,nlsp=length(unique(lsp)),MLL = MLL)
	length(unique(lsp))
}
