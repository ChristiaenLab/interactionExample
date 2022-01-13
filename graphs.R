library(igraph)

#' Writes an igraph object to a pdf after applying a color function to the edges.
#' @param g An igraph
#' @param out The output file
#' @param layout A plot layout function
#' @return NULL
plotfn <- function(g,out,layout=layout_nicely){
	require(circlize)

	colfn <- colorRamp2(c(max(E(g)$weight),0),c('red','white'))
	E(g)$edge.color <- colfn(E(g)$weight)

	lgd <- seq(round(max(E(g)$weight)),
		   0,length.out=6)

	#draw curved edges only if they're reciprocal
	tmp <- do.call(rbind,strsplit(as_ids(E(g)),'\\|'))

	curved <- sapply(1:nrow(tmp),function(x) { 
		y <- which(tmp[,1]==tmp[x,2]&tmp[,2]==tmp[x,1]) 
		if(length(y)>0) return(c(x,y))
	})

	curves <- curve_multiple(g)
	curves[unique(unlist(curved))] <- 0.3

	pdf(out)
	plot(g,vertex.shape='circle',
	     edge.curved=curves,
	     autocurve=T,
	     edge.arrow.size=0.5,
	     edge.color=colfn(E(g)$weight))
	legend('topleft',as.character(lgd),fill=colfn(lgd),title='weight')
	dev.off()
}

interactions <- read.csv('interactions.csv')
homologs <- read.csv('mmusculusHomologs.csv')
row.names(homologs) <- homologs$protein_external_id

interactions[,1] <- homologs[interactions[,1],'name']
interactions[,2] <- homologs[interactions[,2],'name']
interactions <- interactions[interactions[,1]!=interactions[,2],]

interactions <- interactions[!duplicated(interactions),]
interactions <- split(interactions,paste(interactions[,1],interactions[,2]))
interactions <- do.call(rbind,lapply(interactions,function(x) x[which.max(x[,3]),]))

g1 <- graph_from_data_frame(interactions,F)
E(g1)$weight <- interactions[,3]
plotfn(g1,'interactions.pdf')

odds <- read.csv('conditionPoisLog2OR.csv',row.names=1)
fdr <- read.csv('conditionPoisFDR.csv',row.names=1)

odds[odds<0] <- 0
odds[fdr>0.05] <- 0

g2 <- graph_from_adjacency_matrix(as.matrix(odds),'directed',T,F)

plotfn(g,'conditionNetwork.pdf')

