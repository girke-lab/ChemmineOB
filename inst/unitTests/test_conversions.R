
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
}
