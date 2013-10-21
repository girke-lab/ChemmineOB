


#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
#include <openbabel/fingerprint.h>
#include <openbabel/mol.h>
#include <openbabel/oberror.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string>
#include <iomanip>

using namespace std;
using namespace OpenBabel;

extern "C" {
	SEXP ob_convert_file(SEXP fromE,SEXP toE, SEXP sourceFileE,SEXP destinationFileE);
	SEXP ob_convert(SEXP fromE,SEXP toE, SEXP sourceStrE);
	//SEXP propOB(SEXP fromFormatE, SEXP sourceStrE,SEXP descriptorNamesE);
	SEXP propOB(SEXP fromFormatE, SEXP sourceStrE,SEXP descriptorNamesE, SEXP strDescNamesE);
	SEXP fingerprintOB(SEXP fromFormatE, SEXP sourceStrE,SEXP fiingerprintNameE);
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

SEXP propOB(SEXP fromFormatE, SEXP sourceStrE,SEXP descriptorNamesE, SEXP strDescNamesE)
{
	const int numDescriptors = length(descriptorNamesE);
	const int numStrDescriptors = length(strDescNamesE);
	//const char* descriptorNames[] = 

	const char* from = CHAR(STRING_ELT(fromFormatE,0));
	istringstream ifs(CHAR(STRING_ELT(sourceStrE,0)));
   OpenBabel::OBConversion conv(&ifs);
	OBMol mol;
	if(conv.SetInFormat(from))
	{
		SEXP doubleResults,stringResults;
		int numObjects = conv.NumInputObjects();
		//Rprintf("found %d compounds in input string\n",numObjects);
		PROTECT(doubleResults= allocMatrix(REALSXP,numObjects,numDescriptors));
		PROTECT(stringResults= allocMatrix(STRSXP,numObjects,numStrDescriptors));
		int count=0;
		while( conv.Read(&mol))
		{
			OBDescriptor* pDescr;
			for(int i=0; i < numDescriptors; i++){
				if((pDescr =OBDescriptor::FindType( CHAR(STRING_ELT(descriptorNamesE,i)) ) )){
					REAL(doubleResults)[i*numObjects +count] = pDescr->Predict(&mol);
				} else{
					error("Could not find descriptor  %s",CHAR(STRING_ELT(descriptorNamesE,i)));
					UNPROTECT(1);
					return R_NilValue;
				}
			}
			for(int i=0; i < numStrDescriptors; i++){
				if((pDescr =OBDescriptor::FindType( CHAR(STRING_ELT(strDescNamesE,i)) ) )){
					string value;
					pDescr->GetStringValue(&mol,value);
					SET_STRING_ELT(stringResults,i*numObjects +count, mkChar(value.c_str()));
				} else{
					error("Could not find descriptor  %s",CHAR(STRING_ELT(descriptorNamesE,i)));
					UNPROTECT(1);
					return R_NilValue;
				}
			}
			count++;
		}

		SEXP result;
		PROTECT(result =  allocVector(VECSXP,2));
		SET_VECTOR_ELT(result,0,doubleResults);
		SET_VECTOR_ELT(result,1,stringResults);
		
		UNPROTECT(3);
		return result;

	}else{
		error("Could not read given compound in format %s",from);
		return R_NilValue;
	}
}
 
static inline unsigned short bswap_16(unsigned short x) {
      return (x>>8) | (x<<8);
}

static inline unsigned int bswap_32(unsigned int x) {
	  return (bswap_16(x&0xffff)<<16) | (bswap_16(x>>16));
}

char *binary (unsigned int v) {
	static char binstr[17] ;
	int i ;

	binstr[16] = '\0' ;
	for (i=0; i<16; i++) {
		binstr[15-i] = v & 1 ? '1' : '0' ;
		v = v / 2 ;
	}

	return binstr ;
}

SEXP fingerprintOB(SEXP fromFormatE, SEXP sourceStrE,SEXP fingerprintNameE)
{
	const char* from = CHAR(STRING_ELT(fromFormatE,0));
	istringstream ifs(CHAR(STRING_ELT(sourceStrE,0)));
   OpenBabel::OBConversion conv(&ifs);
	OBMol mol;
	int numBits=-1;//4096;

	OBFingerprint* pFP = OBFingerprint::FindFingerprint(CHAR(STRING_ELT(fingerprintNameE,0)));

	if(conv.SetInFormat(from))
	{
		SEXP result;
		int numObjects = conv.NumInputObjects();
		//Rprintf("found %d compounds in input string\n",numObjects);
		int count=0;
		while( conv.Read(&mol))
		{
			vector<unsigned int> fp;
			if(pFP!=0){
				pFP->GetFingerprint(&mol,fp);
				if(numBits == -1) { //not initalized yet
					numBits = fp.size()*sizeof(int)*8;
					//cout<<"real bits: "<<numBits<<endl;
					PROTECT(result = allocMatrix(REALSXP,numObjects,numBits));
				}

				for(int i=0; i< numBits; i++)
					REAL(result)[i*numObjects +count] = pFP->GetBit(fp,i)?1:0;

				//cout<<"bits at ";
				//for(int i=0; i< numBits;i++)
				//	if(pFP->GetBit(fp,i))
				//		cout<<i<<",";
				//cout<<endl;

			}else{
					error("Could not find fingerping %s",CHAR(STRING_ELT(fingerprintNameE,0)));
					UNPROTECT(1);
					return R_NilValue;
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
