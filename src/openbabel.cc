


#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
#include <openbabel/mol.h>
#include <openbabel/oberror.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string>

using namespace std;
using namespace OpenBabel;

extern "C" {
	SEXP ob_convert_file(SEXP fromE,SEXP toE, SEXP sourceFileE,SEXP destinationFileE);
	SEXP ob_convert(SEXP fromE,SEXP toE, SEXP sourceStrE);
	SEXP propOB(SEXP fromFormatE, SEXP sourceStrE,SEXP descriptorNamesE);
}


SEXP ob_convert_file(SEXP fromE,SEXP toE, SEXP sourceFileE,SEXP destinationFileE)
{
	const char* from = CHAR(STRING_ELT(fromE,0));
	const char* to = CHAR(STRING_ELT(toE,0));

   ifstream ifs(CHAR(STRING_ELT(sourceFileE,0)));
   if(!ifs){
      error("cannot open smile file for input: %s",CHAR(STRING_ELT(sourceFileE,0)));
      return R_NilValue;
   }

   ofstream ofs(CHAR(STRING_ELT(destinationFileE,0)));
   if(!ofs){
      error("cannot open sdf file for output: %s",CHAR(STRING_ELT(destinationFileE,0)));
      return R_NilValue;
   }

   OpenBabel::OBConversion conv(&ifs,&ofs);
	if(!conv.SetInAndOutFormats(from,to)) {
      error("conversion from %s to %s is not available",from,to);
      return R_NilValue;
   }
	if(conv.FindFormat(from)==0 ){
		error("could not find format %s",from);
		return R_NilValue;
	}
	if(conv.FindFormat(to)==0){
		error("could not find format %s",to);
		return R_NilValue;
	}

   conv.AddOption("gen2D",OBConversion::GENOPTIONS);
	
   conv.Convert();

	return R_NilValue;

}
SEXP ob_convert(SEXP fromE,SEXP toE, SEXP sourceStrE)
{
	
	const char* from = CHAR(STRING_ELT(fromE,0));
	const char* to = CHAR(STRING_ELT(toE,0));

	istringstream ifs(CHAR(STRING_ELT(sourceStrE,0)));
	ostringstream ofs;


   OpenBabel::OBConversion conv(&ifs,&ofs);
	if(!conv.SetInAndOutFormats(from,to)) {
      error("conversion from %s to %s is not available",from,to);
      return R_NilValue;
   }
	if(conv.FindFormat(from)==0 ){
		error("could not find format %s",from);
		return R_NilValue;
	}
	if(conv.FindFormat(to)==0){
		error("could not find format %s",to);
		return R_NilValue;
	}

   conv.AddOption("gen2D",OBConversion::GENOPTIONS);
	
   conv.Convert();

	SEXP sdfString;
   PROTECT(sdfString = NEW_STRING(1));
   SET_STRING_ELT(sdfString,0,mkChar(ofs.str().c_str()));
   UNPROTECT(1);


   return sdfString;
}

SEXP propOB(SEXP fromFormatE, SEXP sourceStrE,SEXP descriptorNamesE)
{
	const int numDescriptors = length(descriptorNamesE);
	//const char* descriptorNames[] = 

	const char* from = CHAR(STRING_ELT(fromFormatE,0));
	istringstream ifs(CHAR(STRING_ELT(sourceStrE,0)));
   OpenBabel::OBConversion conv(&ifs);
	OBMol mol;
	if(conv.SetInFormat(from))
	{
		SEXP result;
		int numObjects = conv.NumInputObjects();
		//Rprintf("found %d compounds in input string\n",numObjects);
		PROTECT(result = allocMatrix(REALSXP,numObjects,numDescriptors));
		int count=0;
		while( conv.Read(&mol))
		{
			OBDescriptor* pDescr;
			for(int i=0; i < numDescriptors; i++){
				if(pDescr =OBDescriptor::FindType( CHAR(STRING_ELT(descriptorNamesE,i)) ) ){
					double val = pDescr->Predict(&mol);

					//REAL(result)[count*numDescriptors + i] = val;
					REAL(result)[i*numObjects +count] = val;

				} else{
					error("Could not find descriptor  %s",CHAR(STRING_ELT(descriptorNamesE,i)));
					UNPROTECT(1);
					return R_NilValue;
				}
			}
			count++;
		}

		UNPROTECT(1);
		return result;

	}else{
		error("Could not read given compound in format %s",from);
		return R_NilValue;
	}
}
