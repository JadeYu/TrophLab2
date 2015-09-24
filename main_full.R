source("MERA_SSN.R")
source("fullweb_withR.R")
source("analyze_fullweb.R")

nsp <- 30


mu <- 10
theta_seq <- rgeom(nsp,1/(mu+1))+1

Dr_seq <- rnorm(nsp,0.5,0.1)
Dr_seq <- theta_seq^0.5/max(theta_seq^0.5)*0.5+0.2

R0 <- 100

tau_u <- 0.1

FFW <- solve_full_web(theta_seq,Dr_seq,R0,tau_u)
R_inter <- FFW$R_seq

R_no_inter <- solve.analytical(theta_seq,Dr_seq,R0,1)*tau_u*theta_seq

plot_Rdistr(R_inter,R_no_inter)

connectance <- 0.3
trimmed_mat <- Trim_web(FFW$link_mat,L=connectance*nsp^2)

num_basal(trimmed_mat)
num_predator(trimmed_mat)
num_top(trimmed_mat)


chains <- get_chains(trimmed_mat)
Chain_length(chains,loop=T)
loops(chains)

num_cann(trimmed_mat)
