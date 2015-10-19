Dr_seq <- runif(10,0,0.8)
theta_seq <- runif(10,0,1)
R0 <- 1
FFW <- solve_full_web(theta_seq,Dr_seq,R0)
##compare webs with different Dr-theta relationship (random, positive or negative)
##plot without the resource!!

##Incorporate basal-or-not by brutal binary selection (whether basal link > sum of other links); might be mechanistically hard to justify
##The conundrum for Dr distribution is still unsolved.

library(nleqslv)
Ri2j <- function(C_i,Dr_j,theta_j){
	D <- 1/(Dr_j-1)
	(C_i*theta_j)^D * theta_j
}


in_flows <- function(C_seq,Dr_i,theta_i){
	flows <- numeric(length(C_seq))
	for(s in 1:length(C_seq)){
		flows[s] <- Ri2j(C_seq[s],Dr_i,theta_i)
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

R_eqs <- function(C_seq,theta_seq,Dr_seq,R0){
	EQS <- numeric(length(C_seq))
	for(i in 1:(length(C_seq)-1)){
			EQS[i] <- in_flows(C_seq,Dr_seq[i],theta_seq[i]) - out_flows(C_seq[i],Dr_seq,theta_seq)
	}
	EQS[length(C_seq)] <- sum(all_links(C_seq,theta_seq,Dr_seq)$R_seq)-R0
	EQS
}

all_links <- function(C_seq,theta_seq,Dr_seq){
	link_mat <- matrix(nrow=length(C_seq),ncol=length(C_seq))
	for(i in 1:length(C_seq)){
		for(j in 1:length(C_seq)){
			link_mat[i,j] <- Ri2j(C_seq[i],Dr_seq[j],theta_seq[j])
		}
	}
	R_seq <- rowSums(link_mat)
	list(link_mat=link_mat,R_seq=R_seq)
}


solve_full_web <- function(theta_seq,Dr_seq,R0){
	MERA_R <- theta_seq^(Dr_seq/(1-Dr_seq)) #
	R_seq <- R0*MERA_R/sum(MERA_R)
	C_init <- (R_seq/sum(MERA_R))^(mean(Dr_seq)-1)
	C.sol <- nleqslv(C_init,R_eqs,theta_seq=theta_seq,Dr_seq=Dr_seq,R0=R0)
	if(C.sol[[3]]!=1){
		print(C.sol[[4]])
	}
	C_seq <- C.sol[[1]]
	all_links(C_seq,theta_seq,Dr_seq)
}