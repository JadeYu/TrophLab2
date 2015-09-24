source("MERA_SSN.R")
source("fullweb_withR.R")
source("analyze_fullweb.R")

nsp <- 30
Dr_seq <- runif(nsp,0.1,0.8)

theta_seq <- runif(nsp,1,2)
theta_seq <- 2*Dr_seq

R0 <- 100

tau_u <- 0.1

FFW <- solve_full_web(theta_seq,Dr_seq,R0,tau_u)
R_inter <- FFW$R_seq

R_no_inter <- solve.analytical(theta_seq,Dr_seq,R0,1)*tau_u*theta_seq

plot_Rdistr(R_inter,R_no_inter)

connectance <- 0.3
trimmed_mat <- Trim_web(FFW$link_mat,L=connectance*nsp^2)
chains <- get_chains(trimmed_mat)


num_basal(trimmed_mat)
num_predator(trimmed_mat)
num_top(trimmed_mat)

Chain_length(chains,loop=T)
loops(chains)

num_cann(trimmed_mat)
