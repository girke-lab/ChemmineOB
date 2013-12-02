debug =0

test.propOB <-function(){
	numDescs = 12 # this can change

	molRefs = forEachMol("SMI","C1CCCCC1",c,identity)

	p1 = prop_OB(molRefs)
	if(debug) print(p1)
	checkEquals(nrow(p1),1)
	checkEquals(ncol(p1),numDescs)
	#checkEquals(p1$atoms,6)
	checkEquals(p1$MW,84.15948)

	require(ChemmineR)
	data(sdfsample)
	n=5
	p2 = prop_OB(obmol(sdfsample[1:n]))

	if(debug) print(p2)
	checkEquals(nrow(p2),n)
	checkEquals(ncol(p2),numDescs)
	checkEquals(p2$MW[2],MW(sdfsample[2])[[1]])

}
test.fingerprintOB <-function(){

	molRefs = forEachMol("SMI","C1CCCCC1\ncc1ccc1",c,identity)
	f1 = fingerprint_OB(molRefs,"FP2")
	checkEquals(which(f1[1,]==1)-1,c(260,384,429,441,670,984))

	checkException(fingerprint_OB(molRefs,"badfingerprintname"))
}


