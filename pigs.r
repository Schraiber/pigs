require("phytools")
require("phangorn")

#STOLEN FROM http://blog.phytools.org/2012/02/mrca-for-set-of-taxa.html
oldest.mrca<-function(tree,tips){
   H<-nodeHeights(tree)
   X<-mrca(tree)
   n<-length(tips)
   nodes<-height<-vector(); k<-1
   for(i in 1:(n-1)) for(j in (i+1):n){
      nodes[k]<-X[tips[i],tips[j]]
      height[k]<-H[match(nodes[k],tree$edge[,1]),1]
      k<-k+1
   }
   z<-match(min(height),height)
   return(nodes[z])
}

nestedSubclades = function(tree,t1,t2) {
	mrca1 = oldest.mrca(tree,t1)
	mrca2 = oldest.mrca(tree,t2)
	mrcaFull = mrca(tree, full=TRUE)
	mrcaBoth = mrcaFull[mrca1,mrca2]
	if (mrca1!=mrcaBoth && mrca2!=mrcaBoth) {
		return(FALSE)
	} else {
		return(TRUE)
	}
	
}

testNest = function(dna,t1,t2,slide=dim(test)[2]) {
	dnaLength = dim(dna)[2]
	print(dnaLength)
	i = 0 
	nested = c()
	while (i*slide < dnaLength) {
		start = (i*slide+1)
		end = min((i+1)*slide,dnaLength)
		print(c(start,end))
		curDNA = dna[,start:end]
		curDist = dist.dna(curDNA)
		curUPGMA = upgma(curDist)
		nested = c(nested,nestedSubclades(curUPGMA,t1,t2))
		i = i + 1
	}
	return(nested)
}

testNestBulk = function(files,inds1,inds2,slide=1000) {
	res = list()
	for (i in 1:length(files)) {
		curDNA = read.dna(files[i],format="fasta")
		t1inds = sapply(inds1,grep,labels(curDNA))
		dim(t1inds) = NULL
		t1 = labels(curDNA)[t1inds]
		t2inds = sapply(inds2,grep,labels(curDNA))
		dim(t2inds) = NULL
		t2 = labels(curDNA)[t2inds]
		res[[i]] = testNest(curDNA,t1,t2,slide=slide)
	}
	return(res)
}

getUPGMAtree = function(file) {
	dna = read.dna(file,format="fasta")
	dist = dist.dna(dna)
	upgma(dist)
}

plotUPGMAcolored = function(file,inds,colors,show.tip.label=TRUE,exclude=NULL) {
	#inds is a LIST of different groups of individuals
	#colors is a VECTOR of colors, same length as inds
	tree = getUPGMAtree(file)
	#print("inferred UPGMA tree")
	if (!is.null(exclude)) {
		toExclude = sapply(exclude,grep,tree$tip.label)
		toExclude = unlist(toExclude)
		dim(toExclude) = NULL
		tree = drop.tip(tree,tree$tip.label[toExclude])
	}
	col = rep("black",length(tree$edge.length))
	for (i in 1:length(inds)) {
		curIndices = sapply(inds[[i]],grep,tree$tip.label)
		curIndices = unlist(curIndices)
		dim(curIndices)=NULL
		curEdges = match(curIndices,tree$edge[,2])
		col[curEdges]=colors[i]
	}
	#print("prepared colors")
	plot(tree,cex=.5,adj=.5,edge.color=col,edge.width=2,show.tip.label=show.tip.label)
	invisible(list(tree=tree,col=col))
}
