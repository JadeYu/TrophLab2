## MERA competition model
#library(nleqslv)

solve.sameDr <- function(theta_seq,Dr,R){
	base <- sum(theta_seq^(Dr/(Dr-1)))
	R* theta_seq^(1/(Dr-1))/base
}

solve.analytical <- function(theta_seq,Dr_seq,R,r){##given r solve for N+g; r=1 is the special case for equilibrium where N+g=Ne
	if(length(theta_seq)==1){
		sol <- R/theta_seq[1]
	}else if(sum(Dr_seq!=mean(Dr_seq))==0){
		sol <- solve.sameDr(theta_seq,Dr_seq[1],R)
	}else{
		r_seq <- c(r,rep(1,length(theta_seq)-length(r)))
		C <- get.C(theta_seq,Dr_seq,R,r_seq)
		sol <- analytical.sol(C,theta_seq,Dr_seq,r_seq)
	}
	sol
}

analytical.sol <- function(C,theta_seq,Dr_seq,r_seq){
	2/exp(1)*(C*theta_seq*r_seq^(0.5/theta_seq))^(1/(Dr_seq-1))
}

get.C <- function(theta_seq,Dr_seq,R,r_seq){
	init.C <- 1/mean(theta_seq)
	C <- nleqslv(init.C,R_constraint,theta_seq=theta_seq,Dr_seq=Dr_seq,R=R,r_seq=r_seq)
	if(C[[3]]!=1&&C[[3]]!=2){
		print(C[[4]])
	}
	C[[1]]
}

R_constraint <- function(C,theta_seq,Dr_seq,R,r_seq){
	R-sum(theta_seq*analytical.sol(C,theta_seq,Dr_seq,r_seq))
}

steady.state.mR <- function(Sp_mat,R_mat,combine=T){
	SSN_mat <- matrix(ncol=dim(R_mat)[1],nrow=dim(Sp_mat)[1])
	for(i in 1:dim(R_mat)[1]){
		SSN_mat[,i] <- solve.analytical(Sp_mat[,1]/R_mat[i,2],Sp_mat[,2],R_mat[i,1],1)
	}
	rownames(SSN_mat) <- paste("Sp",1:dim(SSN_mat)[1])
	colnames(SSN_mat) <- paste("R",1:dim(SSN_mat)[2])
	if(combine){
		SSN_mat <- rowSums(SSN_mat)
	}
	SSN_mat
}