debug =0

test.propOB <-function(){
	numDescs = 15 # this can change

	p1 = prop_OB("SMI","C1CCCCC1")
	if(debug) print(p1)
	checkEquals(nrow(p1),1)
	checkEquals(ncol(p1),numDescs)
	#checkEquals(p1$atoms,6)
	checkEquals(p1$MW,84.15948)

	require(ChemmineR)
	data(sdfsample)
	n=5
	defs = Map(function(x) paste(x,collapse="\n"),as(as(sdfsample[1:n],"SDFstr"),"list"))
	p2 = prop_OB("SDF",paste(defs,"\n",sep="",collapse=""))

	if(debug) print(p2)
	checkEquals(nrow(p2),n)
	checkEquals(ncol(p2),numDescs)
	checkEquals(p2$MW[2],MW(sdfsample[2])[[1]])

}
test.fingerprintOB <-function(){

	f1 = fingerprint_OB("SMI","C1CCCCC1","FP2")
	checkEquals(which(f1==1)-1,c(260,384,429,441,670,984))
	#print(paste(f1[1,],collapse=""))
	
}


