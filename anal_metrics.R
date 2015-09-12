#calculate analytical metrics
is.connected <- function(index,community){
	rowSums(community$link_mat)!=0|colSums(community$link_mat)!=0
}

is.basal <- function(index,community){
	if(length(index)==1){
		sum(community$link_mat[,index])==0&sum(community$link_mat[index,])!=0
	}else{
		colSums(community$link_mat[,index])==0&rowSums(community$link_mat[index,])!=0
	}
}

is.top <- function(index,community){
	colSums(community$link_mat)>0&rowSums(community$link_mat)==0
}

num_links <- function(community){## not including links to the fundamental resource
	sum(community$link_mat[1:dim(community$link_mat)[2],]>0)
}

connectance <- function(community){
	totalS <- num_predators(community)+num_basal(community)
	C= num_links(community)/totalS^2
	C
}

num_basal <- function(community){## does not include cannibals
	sum(colSums(community$link_mat)==0&rowSums(community$link_mat)!=0)
}

num_predators <- function(community){##species that are feeding on at least one other species
	sum(colSums(community$link_mat)>0)
}

num_top <- function(community){##species that are not fed on by any other species
	diag(community$link_mat) <- 0
	sum(colSums(community$link_mat)>0&rowSums(community$link_mat)==0)
}

get_chains <- function(community){
	link_mat <- community$link_mat
	indice <- community$sp_mat$index
	basals <- indice[is.basal(indice,community)]
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
	mean(unlist(lapply(chains,length)))
}

Loop_length <- function(chain){
	if(length(chain)==length(unique(chain))){
		LL = 0
	}else{
		sps <- unique(chain)
		sp <- sort(sps)[table(chain)>1]
		
		indice <- (1:length(chain))[chain==sp]
		LL = indice[2]-indice[1]+1
	}
	LL
}

loops <- function(chains){
	LL_seq <- unlist(lapply(chains,Loop_length))
	N_loops <- sum(LL_seq>0)
	MLL <- mean(LL_seq[LL_seq>0])
	list(N_loops = N_loops,MLL = MLL)
}
