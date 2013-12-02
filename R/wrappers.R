debug=0
packageName = "ChemmineOB"

.onLoad <- function(libname,pkgname) {
	#The package library must be loaded globally (local=FALSE) in order for 
	# the openbabel plugins to load properly. Further, when calling
	# functions from this library, you must set PACKAGE=packageName in the .Call function
	# or it will not find any of the symbols.

	library.dynam(pkgname,package=pkgname,lib.loc=libname,local=FALSE)

	if(.Platform$OS == "windows")
	 Sys.setenv(BABEL_DATADIR=system.file("openbabel_data",package=pkgname))

}

convertFormat <- function(from,to,source){
	
	inStr = istreamFromString(source)
	outStr = ostreamToString()

	conv = OBConversion(inStr,outStr)

	if(!OBConversion_SetInAndOutFormats(conv,from,to))
		stop("failed to set 'from' and
			  'to' formats: ",from," ",to)
	OBConversion_AddOption(conv,"gen2D","GENOPTIONS")
	OBConversion_Convert(conv)

	stringFromOstream(outStr)
}

convertFormatFile <- function(from,to,fromFile,toFile){
	is= istreamFromFile(fromFile)
	os= ostreamToFile(toFile)

	conv = OBConversion(is,os)

	if(!OBConversion_SetInAndOutFormats(conv,from,to))
		stop("failed to set 'from' and 'to' formats: ",from," ",to)
	OBConversion_AddOption(conv,"gen2D","GENOPTIONS")
	OBConversion_Convert(conv)

	stringFromOstream(outStr)
}

prop_OB<- function(obmolRefs) {
	if(length(obmolRefs)==1)
		obmolRefs=c(obmolRefs)

	strDescrNames = c( "cansmi",
							  "cansmiNS",
							  "formula",
							  "title")
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

fingerprint_OB <- function(obmolRefs, fingerprintName){

	OBConversion()

	numBits = -1;
	fpHandle = OBFingerprint_FindFingerprint(fingerprintName)
	if(isNullPtr(fpHandle))
		stop("fingerprint ",fingerprintName," not found")

	Reduce(rbind,Map(function(mol){
		fp = vectorUnsignedInt(1)

		OBFingerprint_GetFingerprint(fpHandle,mol,fp)
		if(numBits == -1)
			numBits = vectorUnsignedInt_size(fp) * 4 * 8
		row = unlist(Map(function(i){
					 r=OBFingerprint_GetBit(fpHandle,fp,i-1)
					 if(r) 1 else 0
				},seq(1,numBits,length.out=numBits)))
		if(debug) print(row)
		row
	  },obmolRefs))
}
forEachMol <- function(inFormat,inString,reduce,f){

	inStr = istreamFromString(inString)
	conv = OBConversion(inStr)

	if(!OBConversion_SetInFormat(conv,inFormat))
		stop("failed to set input format: ",inFormat)
	
	numMols = OBConversion_NumInputObjects(conv)
	Reduce(reduce,Map(function(i){
		mol = OBMol()
		if(!OBConversion_Read(conv,mol))
			stop("failed to read ",numMols," from input")
		f(mol)
	},seq(1,numMols,length.out=numMols)))
}

