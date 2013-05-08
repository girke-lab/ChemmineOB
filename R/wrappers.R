
convertFormat <- function(from,to,source){
	.Call("ob_convert",as.character(from),as.character(to),as.character(source))
}

convertFormatFile <- function(from,to,fromFile,toFile){
	.Call("ob_convert_file",as.character(from),as.character(to),as.character(fromFile),as.character(toFile))
}

