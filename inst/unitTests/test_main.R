debug =1

test.propOB <-function(){
	numDescs = 13 # this can change

	molRefs = forEachMol("SMI","C1CCCCC1",identity)

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
	
	molRefs = forEachMol("SMI","C1CCCCC1\ncc1ccc1",identity)
	library(ChemmineR)
	data(sdfsample)
	molRefs = obmol(sdfsample)
	x = smartsSearch_OB(molRefs,"[CH3X4]")
	checkEquals(length(x[x>0]),80)

	rotableBonds = smartsSearch_OB(molRefs[1:5],"[!$(*#*)&!D1]-!@[!$(*#*)&!D1]",uniqueMatches=FALSE)
	print(rotableBonds)
	print(sdfid(sdfsample[1:5]))
	checkEquals(as.vector(rotableBonds[1:5]),c(24,20,14,30,10))
}

test.exactMassOB <- function(){
	#molRefs = forEachMol("SMI","C1CCCCC1\ncc1ccc1",identity)
	#masses = exactMass_OB(molRefs)
	#print("masses:")
	#print(masses)
	library(ChemmineR)
	data(sdfsample)
	m = exactMass_OB(obmol(sdfsample[1]))
	message("650001 mass: ",m)
	checkEqualsNumeric(m,456.200885,tolerance=0.000001)


}
test.canonicalNumbers <- function(){
	library(ChemmineR)
	data(sdfsample)

	#print(bondblock(sdfsample[[1]]))
	sdfstrList=as(as(sdfsample[1],"SDFstr"),"list")
	sdfDef= paste(Map(function(x) paste(x,collapse="\n"),
						  sdfstrList),collapse="\n" )

	canSdfDef=convertFormat("SDF","SDF",sdfDef,options=data.frame(names=c("canonical","gen2d"),args=""))
	canSdf = read.SDFset(unlist(strsplit(canSdfDef,"\n",fixed=TRUE)))
	#print(bondblock(canSdf[[1]]))
	#write.SDF(canSdf,file="sdf1-renumbered.sdf")
	bb=bondblock(canSdf[[1]])
	checkEqualsNumeric(bb[1,1:3],c(2,3,1))
	checkEqualsNumeric(bb[2,1:3],c(2,4,1))


#	ref= obmol(sdfsample[1])
#	canonicalNumbering(ref)




	#ChemmineOB:::writeMols(ref,file="sdf1-renumbered-noH.sdf","SDF")

#
#	dir = tempdir()
#	f = file.path(dir,"sdfsample1.sdf")
#	write.SDF(sdfsample[1],file=f)
#
##	sdf = read.SDFset(f)[[1]]
##	ref= obmol(sdf)
##	ChemmineOB:::writeMols(ref,file="sdf1-orig.sdf","SDF")
##	ChemmineOB:::writeMols(ref,file="sdf1-orig.smi","SMI")
#
#	sdf = read.SDFset(f)[[1]]
#	ref= obmol(sdf)
#	canonicalNumbering(ref)
#	ChemmineOB:::writeMols(ref,file="sdf1-renumbered.sdf","SDF")
#	ChemmineOB:::writeMols(ref,file="sdf1-renumbered.smi","SMI")
#
#
#	#sdf <- sdfsample[1]
#	sdf = read.SDFset(f)[[1]]
#	#ref= obmol(sdf)
#
#	#sdf[["bondblock"]] <- bondblock(sdf)[nrow(bondblock(sdf)):1,]
#	#swap first two atoms
#	tempAtom = atomblock(sdf)[1,]
#	sdf[["atomblock"]][1,] = atomblock(sdf)[2,]
#	sdf[["atomblock"]][2,] = tempAtom
#
#	#renumber bond block to account for differnt atom numbers
#	sdf[["bondblock"]][ bondblock(sdf)[,1]==1,1] = 999
#	sdf[["bondblock"]][ bondblock(sdf)[,1]==2,1] = 1
#	sdf[["bondblock"]][ bondblock(sdf)[,1]==999,1] = 2
#
##	write.SDF(sdf,"sdf1-before.sdf")
##	sdf = as(sdf[[1]],"SDFset")[[1]]
#	#sdf = as(sdf,"SDFset")
##	write.SDF(sdf,"sdf1-after.sdf")
#
#	revref= obmol(sdf)
#	canonicalNumbering(revref)
#	ChemmineOB:::writeMols(revref,file="sdf1-rev-renumbered.smi","SMI")
#	ChemmineOB:::writeMols(revref,file="sdf1-rev-renumbered.sdf","SDF")
#	ChemmineOB:::writeMols(revref,file="sdf1-rev-renumbered.can","CAN")
##	print(atomblock(sdfsample[1]))
	
}
