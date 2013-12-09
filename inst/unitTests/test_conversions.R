
debug=1

test.formatConversions <- function(){
 
	smiles = "OC(=O)C(Br)(Cl)N\ttest1\n"

	sdf=convertFormat("SMI","SDF",smiles)
	s2 = convertFormat("SDF","SMI",sdf)
	
	if(debug) {
		print(sdf)
		print(smiles)
		print(s2)
	}
	checkEquals(smiles,s2)

	checkException(convertFormat("BADFORMAT","SDF",smiles))

	t1=file.path(tempdir(),"t1.smi")
	# 650001
	cat("CC1=CC(=NO1)NC(=O)CCC(=O)N(CC(=O)NC2CCCC2)C3=CC4=C(C=C3)OCCO4",
		 file=t1)
	convertFormatFile("SMI","svg",t1,"test.svg")
	checkTrue(file.exists("test.svg"))
	checkTrue(file.info("test.svg")$size > 0)
	unlink("test.svg")
}
