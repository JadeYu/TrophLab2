source("MERA_SSN.R")
source("fullweb_withR.R")
source("analyze_fullweb.R")
source("Visualization.R")

nsp <- 50
mu <- 10
theta_seq <- runif(nsp,1,2*mu-1)
Dr_seq <- rnorm(nsp,0.5,0.1)+0.2*(theta_seq-mu)/4/mu
#Dr_seq <- theta_seq^0.5/max(theta_seq^0.5)*0.5+0.2

R0 <- mu*nsp
tau_u <- 0.1

FFW <- solve_full_web(theta_seq,Dr_seq,R0,tau_u)

##compare resource distribution
R_inter <- FFW$R_seq
R_no_inter <- solve.analytical(theta_seq,Dr_seq,R0,1)*tau_u*theta_seq
plot_Rdistr(R_inter,R_no_inter)


##See shape of sub-web given constraints
connectance <- 0.05

trimmed_mat <- Trim_web(FFW$link_mat,connectance*nsp^2,"links")
trimmed_mat <- Trim_web(FFW$link_mat,0.001,"MVFF")

metrics <- web_summary(trimmed_mat,R_inter,verbose=T)

