debug=0
packageName = "ChemmineOB"

.onLoad <- function(libname,pkgname) {
	#The package library must be loaded globally (local=FALSE) in order for 
	# the openbabel plugins to load properly. Further, when calling
	# functions from this library, you must set PACKAGE=packageName in the .Call function
	# or it will not find any of the symbols.
  
  if(! .supportedPlatform()){
    warning("ChemmineOB is not fully supported on this platform. ", .Platform$OS," ", .Platform$r_arch)
  }

	library.dynam(pkgname,package=pkgname,lib.loc=libname,local=FALSE)

	if(.Platform$OS == "windows")
	 Sys.setenv(BABEL_DATADIR=system.file("openbabel_data",package=pkgname))

}
.supportedPlatform <- function(){
  !(.Platform$OS == "windows" && .Platform$r_arch == "i386")
}
convertFormat <- function(from,to,source,options=data.frame(names="gen2D",args="")){
	
	inStr = istreamFromString(source)
	outStr = ostreamToString()

	conv = OBConversion(inStr,outStr)

	if(!OBConversion_SetInAndOutFormats(conv,from,to))
		stop("failed to set 'from' and
			  'to' formats: ",from," ",to)

	for(i in seq(to=nrow(options),length=nrow(options))){
		OBConversion_AddOption(conv,as.character(options[i,1]),"GENOPTIONS",as.character(options[i,2]))
	}
	#OBConversion_AddOption(conv,"gen2D","GENOPTIONS")
	OBConversion_Convert(conv)

	stringFromOstream(outStr)
}

convertFormatFile <- function(from,to,fromFile,toFile, options=data.frame(names="gen2D",args="")){
	is= istreamFromFile(fromFile)
	os= ostreamToFile(toFile)

	conv = OBConversion(is,os)

	if(!OBConversion_SetInAndOutFormats(conv,from,to))
		stop("failed to set 'from' and 'to' formats: ",from," ",to)
	for(i in seq(to=nrow(options),length=nrow(options))){
		OBConversion_AddOption(conv,as.character(options[i,1]),"GENOPTIONS",as.character(options[i,2]))
	}
	#OBConversion_AddOption(conv,"gen2D","GENOPTIONS")
	OBConversion_Convert(conv)

	closeOfstream(os)
}
convertToImage <- function(from,to,source,toFile,options=data.frame(names="gen2D",args=""),
									out_options=data.frame(names=c("p"),args=c(300))){
	
	inStr = istreamFromString(source)
	os= ostreamToFile(toFile)

	conv = OBConversion(inStr,os)

	if(!OBConversion_SetInAndOutFormats(conv,from,to))
		stop("failed to set 'from' and
			  'to' formats: ",from," ",to)

	for(i in seq(to=nrow(options),length=nrow(options))){
		OBConversion_AddOption(conv,as.character(options[i,1]),"GENOPTIONS",as.character(options[i,2]))
	}
	for(i in seq(to=nrow(out_options),length=nrow(out_options))){
		OBConversion_AddOption(conv,as.character(out_options[i,1]),"OUTOPTIONS",as.character(out_options[i,2]))
	}
	#OBConversion_AddOption(conv,"gen2D","GENOPTIONS")
	OBConversion_Convert(conv)

	closeOfstream(os)
}
canonicalNumbering_OB <- function(obmolRefs) {

	if(length(obmolRefs)==1)
		obmolRefs=c(obmolRefs)

	Map(function(mol){
		.Call("canonicalReordering",mol@ref,PACKAGE="ChemmineOB")
	},obmolRefs)
	
}


prop_OB<- function(obmolRefs) {
	if(length(obmolRefs)==1)
		obmolRefs=c(obmolRefs)

	strDescrNames = c( "cansmi",
							  "cansmiNS",
							  "formula",
							  "title",
							  "InChI")
	descriptorNames = c(
							  "HBA1",
							  "HBA2",
							  "HBD",
							  "logP",
							  "MR",
							  "MW",
							  "nF",
							  "TPSA")

	Reduce(rbind,Map( function(mol){
			row=list()
			for(descName  in strDescrNames){
				desc = OBDescriptor_FindType(descName)
				result = stringp()
				OBDescriptor_GetStringValue(desc,mol,result$cast())
				row[[descName]] = result$value()
			}
			for(descName in descriptorNames){
				desc = OBDescriptor_FindType(descName)
				row[[descName]] =OBDescriptor_Predict(desc,mol)
			}
			#if(debug) print(row)
			as.data.frame(row)
	},obmolRefs))
}

fingerprint_OB <- function(obmolRefs, fingerprintName,reverse=FALSE){

	OBConversion()

	numBits = -1;
	fpHandle = OBFingerprint_FindFingerprint(fingerprintName)
	if(isNullPtr(fpHandle))
		stop("fingerprint ",fingerprintName," not found")

	Reduce(rbind,Map(function(mol){
		fp = numeric(1)

		fp=OBFingerprint_GetFingerprint(fpHandle,mol,fp)[[2]] # new fp arg is returned as second element of a list
		if(numBits == -1)
			numBits = length(fp) * 4 * 8
		row = unlist(Map(function(i){
					 r=OBFingerprint_GetBit(fpHandle,fp,i-1)[[1]]
					 if(r) 1 else 0
				},seq(1,numBits,length.out=numBits)))
		if(debug) print(row)

		# for some reason the fingerpints computed here are the 
		# reverse bit string of what obabel outputs. 
		if(reverse) 
			row = rev(row)
		row
	  },obmolRefs))
}
exactMass_OB <- function(obmolRefs){
	unlist(Map(function(mol){
		 OBMol_GetExactMass(mol)
	},obmolRefs))

}
forEachMol <- function(inFormat,inString,f,reduce=NULL){

	inStr = istreamFromString(inString)
	conv = OBConversion(inStr)

	if(!OBConversion_SetInFormat(conv,inFormat))
		stop("failed to set input format: ",inFormat)
	
	numMols = OBConversion_NumInputObjects(conv)
	x=Map(function(i){
		mol = OBMol()
		if(!OBConversion_Read(conv,mol))
			stop("failed to read ",numMols,"molecules from input")
		f(mol)
	},seq(1,numMols,length.out=numMols))

	if(is.null(reduce))
		x
	else
		Reduce(reduce,x)
}
writeMols <- function(obmolRefs,file,format) {

	if(length(obmolRefs)==1)
		obmolRefs=c(obmolRefs)
	conv = OBConversion()
	os= ostreamToFile(file)
	if(!OBConversion_SetOutFormat(conv,format))
		stop("failed to set output format to ",format)
	for(mol in obmolRefs){
		OBConversion_Write(conv,mol,os)
	}
	closeOfstream(os)
}
smartsSearch_OB<- function(obmolRefs,smartsPattern,uniqueMatches=TRUE){

	sp = OBSmartsPattern()
	if(!OBSmartsPattern_Init(sp,smartsPattern))
		stop("failed to parse smarts patter: ",smartsPattern)

	unlist(Map(function(mol){
		OBSmartsPattern_Match(sp,mol)
		if(uniqueMatches){
			umap = OBSmartsPattern_GetUMapList(sp)
			length(umap)
		}else
			OBSmartsPattern_NumMatches(sp)

   },obmolRefs))
}
