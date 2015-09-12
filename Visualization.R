library(igraph)


plot_foodweb <- function(community){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	colnames(link_mat) <- rownames(link_mat) <- c()
	colors = rainbow(3)[as.numeric(!is.basal(sp_mat$index,community))+as.numeric(is.top(sp_mat$index,community))+1]
	node_size <- sp_mat$theta*sp_mat$APN
	AJ = as.matrix(link_mat)
	net=graph.adjacency(AJ,mode="directed",weighted="width",diag=F)

	V(net)$color = colors
	V(net)$size = 10*node_size^0.2
	set.seed(1)

	par(bg="gray",mar=rep(0,4))

		plot.igraph(net,vertex.label=rep("",dim(link_mat)[1]),edge.arrow.size = E(net)$width^0.5/4,vertex.frame.color=0,vertex.color=colors,edge.color="black",edge.width=E(net)$width^0.5*3)
	legend("bottomright",c("basal","intermediate","top"),pch=16,col=rainbow(3))

}
