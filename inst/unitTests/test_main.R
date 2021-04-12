debug =0
library(ChemmineR)
data(sdfsample)

test.propOB <-function(){
	numDescs = 13 # this can change

	molRefs = forEachMol("SMI","C1CCCCC1",identity)

	p1 = prop_OB(molRefs)
	if(debug) print(p1)
	checkEquals(nrow(p1),1)
	checkEquals(ncol(p1),numDescs)
	#checkEquals(p1$atoms,6)
	checkEquals(p1$MW,84.15948)

	n=5
	p2 = prop_OB(obmol(sdfsample[1:n]))

	if(debug){message("probOB: "); print(p2)}
	checkEquals(nrow(p2),n)
	checkEquals(ncol(p2),numDescs)
	checkEquals(p2$MW[2],MW(sdfsample[2])[[1]])

}
test.fingerprintOB <-function(){

	molRefs = forEachMol("SMI","C1CCCCC1\ncc1ccc1",identity)
	f1 = fingerprint_OB(molRefs,"FP2")
	checkEquals(which(f1[1,]==1)-1,c(260,384,429,441,670,984))

	checkException(fingerprint_OB(molRefs,"badfingerprintname"))
}


test.smartsOB <- function(){
	
	message("start of smartsOB")
	molRefs = forEachMol("SMI","C1CCCCC1\ncc1ccc1",identity)
	molRefs = obmol(sdfsample)
	x = smartsSearch_OB(molRefs,"[CH3X4]")
	checkEquals(length(x[x>0]),80)

	message("before smartsSearch_OB")
	rotableBonds = smartsSearch_OB(molRefs[1:5],"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]",uniqueMatches=FALSE)
	message("rotable bonds:")
	print(rotableBonds)
	message("sdfid: ")
	print(sdfid(sdfsample[1:5]))
	message("checking equals...")
	checkEquals(as.vector(rotableBonds[1:5]),c(24,20,14,30,10))
}

test.exactMassOB <- function(){
	#molRefs = forEachMol("SMI","C1CCCCC1\ncc1ccc1",identity)
	#masses = exactMass_OB(molRefs)
	#print("masses:")
	#print(masses)
	m = exactMass_OB(obmol(sdfsample[1]))
	message("650001 mass: ",m)
	checkEqualsNumeric(m,456.200885,tolerance=0.000001)


}
test.canonicalize <- function(){

	sdfstrList=as(as(sdfsample[1],"SDFstr"),"list")
	sdfDef= paste(Map(function(x) paste(x,collapse="\n"),
						  sdfstrList),collapse="\n" )

	canSdfDef=convertFormat("SDF","SDF",sdfDef,options=data.frame(names=c("canonical","gen2d"),args=""))
	canSdf = read.SDFset(unlist(strsplit(canSdfDef,"\n",fixed=TRUE)))

	bb=bondblock(canSdf[[1]])

	checkEqualsNumeric(bb[1,1:3],c(2,3,1))
	checkEqualsNumeric(bb[2,1:3],c(2,4,1))
}
test.canonicalLabels <- function() {

	labels=canonicalNumbering_OB(obmol(sdfsample[[1]]))
	print("labels:")
	print(labels)

	# compare only the first 11 since there are two possible mappings
	# given the symetry of this compound and we can't predict which one
	# it will use.
	checkEqualsNumeric(labels[[1]][1:11],c(25,47,48,16,6,55,24,15,2,7,30,34,33,40,37,19,17,
										 18,12,9,23,22,29,5,46,27,49,50,3,4,8,57,58,32,38,
										 39,36,35,43,44,41,42,31,20,21,13,14,28,26,10,11,
										 45,51,52,53,54,1,56,59,60,61)[1:11])

}
test.formatConversions <- function(){
	sdfFile = tempfile()
	smiFile = tempfile()
	write.SDF(sdfsample[1],sdfFile)
	convertFormatFile("SDF","SMI",sdfFile,smiFile)
	smi = read.SMIset(smiFile)
	checkEquals(smi[[1]]@smiles,"O=C(NC1CCCC1)CN(c1cc2OCCOc2cc1)C(=O)CCC(=O)Nc1noc(c1)C")

	if(.Platform$OS.type != "windows"){
		pngFile = tempfile()
		convertToImage("SMI","PNG",smi[[1]]@smiles,pngFile)
		checkTrue(file.size(pngFile) > 0)
	}
}
test.writeMols <- function(){
	ref = obmol(sdfsample[1])
	sdfFile = tempfile()
	ChemmineOB:::writeMols(ref,sdfFile,"SDF")
	checkTrue(file.size(sdfFile) > 0)
}
