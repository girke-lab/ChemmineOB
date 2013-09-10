
packageName = "ChemmineOB"

.onLoad <- function(libname,pkgname) {
	#The package library must be loaded globally (local=FALSE) in order for 
	# the openbabel plugins to load properly. Further, when calling
	# functions from this library, you must set PACKAGE=packageName in the .Call function
	# or it will not find any of the symbols.

	dyn.load(file.path(libname,pkgname,"libs",paste(pkgname,"so",sep=".")),local=FALSE)

}

convertFormat <- function(from,to,source){
	.Call("ob_convert",as.character(from),as.character(to),as.character(source),PACKAGE=packageName)
}

convertFormatFile <- function(from,to,fromFile,toFile){
	.Call("ob_convert_file",as.character(from),as.character(to),as.character(fromFile),as.character(toFile),PACKAGE=packageName)
}

prop_OB<- function(from,source) {
	descriptorNames = c("abonds", "atoms", "bonds", "dbonds", "HBA1", "HBA2", "HBD", "logP", "MR",
							 "MW", "nF", "sbonds", "tbonds", "TPSA")
	values = .Call("propOB",as.character(from),as.character(source),PACKAGE=packageName)
	df = as.data.frame(values)
	colnames(df) = descriptorNames
	df
}

