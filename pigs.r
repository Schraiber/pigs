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
	dnaLength = dim(test)[2]
	i = 0 
	nested = c()
	while (i*slide <= dnaLength) {
		print(i*slide)
		start = (i*slide+1)
		end = min((i+1)*slide,dnaLength)
		curDNA = dna[,start:end]
		curDist = dist.dna(curDNA)
		curUPGMA = upgma(curDist)
		nested = c(nested,nestedSubclades(curUPGMA,t1,t2))
		i = i + 1
	}
	return(nested)
}
