library(igraph)


plot_foodweb <- function(link_mat,sizes){
	index <- 1:dim(link_mat)[1]
	colors = rainbow(3)[as.numeric(!is.basal(index,link_mat))+as.numeric(is.top(index,link_mat))+1]
	node_size <- sizes
	AJ = as.matrix(link_mat)
	net=graph.adjacency(AJ,mode="directed",weighted="width",diag=F)

	V(net)$color = colors
	V(net)$size = 10*node_size^0.5
	set.seed(2)

	par(bg="gray",mar=rep(0,4))

		plot.igraph(net,vertex.label=rep("",dim(link_mat)[1]),edge.arrow.size = E(net)$width^0.5/4,vertex.frame.color=0,vertex.color=colors,edge.color="black",edge.width=E(net)$width^0.2*4)
	legend("bottomright",c("basal","intermediate","top"),pch=16,col=rainbow(3))

}

web_metrics <- function(trimmed_mat){
	B <- num_basal(trimmed_mat)
 	Top <- num_top(trimmed_mat)
	I <- num_predator(trimmed_mat) - Top
	Total <- sum(B+Top+I)
	L <- sum(trimmed_mat>0)
	Cann <- num_cann(trimmed_mat)
	list(C = L/Total^2, Total = Total,Basal=B/Total,Top=Top/Total,Inter=I/Total,Cann=Cann/Total)
}

chain_metrics <- function(trimmed_mat){
	chains <- get_chains(trimmed_mat)
	chlen <- Chain_length(chains,loop=T)
	#nsp_loop <- loops(chains)
	#list(chlen= chlen,loop = nsp_loop)
	chlen
}

web_summary <- function(link_mat,R_seq,verbose=F){
	plot_foodweb(link_mat,R_seq)
	metrics <- web_metrics(link_mat)
	print(metrics)
	if(verbose){
		extra <- chain_metrics(link_mat)
		extra$loop <- extra$loop/metrics$Total
		print(extra)
		metrics <- c(metrics,extra)
	}
	metrics
}

	
