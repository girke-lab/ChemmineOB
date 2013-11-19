packageName = "ChemmineOB"

.onLoad <- function(libname,pkgname) {
	#The package library must be loaded globally (local=FALSE) in order for 
	# the openbabel plugins to load properly. Further, when calling
	# functions from this library, you must set PACKAGE=packageName in the .Call function
	# or it will not find any of the symbols.
#	message("pkgname: ",pkgname," libname: ",libname)

	library.dynam(pkgname,package=pkgname,lib.loc=libname,local=FALSE)

#	arch = R.version$arch
#	libPrefix = file.path(libname,pkgname,"libs")
#	sharedLibName = paste(pkgname,"so",sep=".")
#
#	candidates = c( file.path(libPrefix,sharedLibName),file.path(libPrefix,arch,sharedLibName))
#
#	lib = Find(file.exists,candidates)
#	if(is.null(lib))
#		stop("Could not find shared library, looked in these places: ",
#			  paste(candidates,collapse=","))
#
#
#	#dyn.load(file.path(libname,pkgname,"libs",paste(pkgname,"so",sep=".")),local=FALSE)
#	dyn.load(lib,local=FALSE)
	if(.Platform$OS == "windows")
	 Sys.setenv(BABEL_DATADIR=system.file("openbabel_data",package=pkgname))

}

convertFormat <- function(from,to,source){
	#.Call("ob_convert",as.character(from),as.character(to),as.character(source),PACKAGE=packageName)
	
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
	#.Call("ob_convert_file",as.character(from),as.character(to),as.character(fromFile),as.character(toFile),PACKAGE=packageName)
	is= istreamFromFile(fromFile)
	os= ostreamToFile(toFile)

	conv = OBConversion(is,os)

	if(!OBConversion_SetInAndOutFormats(conv,from,to))
		stop("failed to set 'from' and 'to' formats: ",from," ",to)
	OBConversion_AddOption(conv,"gen2D","GENOPTIONS")
	OBConversion_Convert(conv)

	stringFromOstream(outStr)

}

prop_OB<- function(from,source) {
#	descriptorNames = c("abonds", "atoms", "bonds", "dbonds", "HBA1", "HBA2", "HBD", "logP", "MR",
#							 "MW", "nF", "sbonds", "tbonds", "TPSA")
	#  "InChIKey",  not availabel on bioc servers
	strDescrNames = c( "cansmi",
							  "cansmiNS",
							  "formula",
							  "title")
	descriptorNames = c(
							  "HBA1",
							  "HBA2",
							  "HBD",
							  #"InChI",
							  "logP",
							  "MR",
							  "MW",
							  "nF",
							  #"s",
							  #"smarts",
							  "TPSA")
	values = obCall("propOB",as.character(from),as.character(source),descriptorNames,strDescrNames)




	#print(values)
	df = as.data.frame(values)
	#print(df)
	colnames(df) = c(descriptorNames,strDescrNames)
	df
}
obCall <-function(...)
	.Call(...,PACKAGE=packageName)

fingerprint_OB <- function(format,source, fingerprintName){
	obCall("fingerprintOB",as.character(format),as.character(source),as.character(fingerprintName))
}
