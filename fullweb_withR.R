

##compare webs with different Dr-theta relationship (random, positive or negative) scenarios

##For each scenario, compare resource rank distribution with when there is no interactions.

##For each scenario, get the top ranking links and compare with food web structure

##plot without the resource!!



library(nleqslv)

Ri2j <- function(C_i,Dr_j,theta_j){

	D <- 1/(Dr_j-1)

	(C_i*theta_j)^D * theta_j

}



in_flows <- function(C_seq,Dr_i,theta_i,tau_u){

	flows <- numeric(length(C_seq))

	for(s in 1:length(C_seq)){

		flows[s] <- Ri2j(C_seq[s],Dr_i,theta_i)*tau_u

	}

	sum(flows)

}



out_flows <- function(C_i,Dr_seq,theta_seq){

	flows <- numeric(length(Dr_seq))

	for(s in 1:length(Dr_seq)){

		flows[s] <- Ri2j(C_i,Dr_seq[s],theta_seq[s])

	}

	sum(flows)

}



R_eqs <- function(C_seq,theta_seq,Dr_seq,R0,tau_u){

	EQS <- numeric(length(C_seq))

	EQS[1] <- R0 - out_flows(C_seq[1],Dr_seq,theta_seq)#for fundamental resource

	for(i in 2:length(C_seq)){

			EQS[i] <- in_flows(C_seq,Dr_seq[i-1],theta_seq[i-1],tau_u) - out_flows(C_seq[i],Dr_seq,theta_seq)

	}

	EQS

}



solve_full_web <- function(theta_seq,Dr_seq,R0,tau_u){

	MERA_R <- theta_seq^(Dr_seq/(1-Dr_seq)) #

	C0 <- (R0*tau_u/sum(MERA_R))^(mean(Dr_seq)-1)

	R_seq <- R0*tau_u*MERA_R/sum(MERA_R)

	C_init <- c(C0,(R_seq/sum(MERA_R))^(mean(Dr_seq)-1))

	C.sol <- nleqslv(C_init,R_eqs,theta_seq=theta_seq,Dr_seq=Dr_seq,R0=R0,tau_u = tau_u)

	if(C.sol[[3]]!=1){

		print(C.sol[[4]])

	}

	C_seq <- C.sol[[1]]

	link_mat <- matrix(nrow=length(C_seq)-1,ncol=length(C_seq)-1)

	for(i in 1:(length(C_seq)-1)){

		for(j in 1:(length(C_seq)-1)){

			link_mat[i,j] <- Ri2j(C_seq[i+1],Dr_seq[j],theta_seq[j])

		}

	}

	R_seq <- rowSums(link_mat)

	#compare with colSums(link_mat)*tau_u

	list(link_mat=link_mat,R_seq=R_seq)

}



