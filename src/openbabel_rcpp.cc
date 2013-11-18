#include <Rcpp.h>
#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
#include <openbabel/fingerprint.h>
#include <openbabel/mol.h>
#include <openbabel/oberror.h>
#include <string>
#include <iomanip>


using namespace Rcpp;
using namespace std;
using namespace OpenBabel;

double normx(double x,double y) {
	return sqrt(x*x+y*y);
}
/*
RcppExport SEXP normx_wrapper ( SEXP x_ , SEXP y_) {
	// step 0: convert input to C ++ types
	double x = as < double >( x_), y = as < double >( y_);
	// step 1: call the underlying C ++ function 
	double res = normx ( x , y);
	// step 2: return the result as a SEXP
	return wrap ( res);
}
*/
RCPP_MODULE(test_mod){
	using namespace Rcpp;
	function("normx",&normx);
}

RCPP_EXPOSED_CLASS(OpenBabel::OBBase)
RCPP_EXPOSED_CLASS(OpenBabel::OBMol)
RCPP_EXPOSED_CLASS(OpenBabel::OBDescriptor)
RCPP_EXPOSED_CLASS(OpenBabel::OBConversion)

double GetStringValue(OBDescriptor* obd,OBBase* pOb,string &value){
	return obd->GetStringValue(pOb,value);
}
double Predict(OBDescriptor* obd,OBBase* pOb){
	return obd->Predict(pOb);
}
OBDescriptor* FindType(string type){
	return OBDescriptor::FindType(type.c_str());
}

RCPP_MODULE(openbabel) {
	using namespace Rcpp;
	using namespace std;
	using namespace OpenBabel;

	//class_<OBMol>("OBMol")
	//.constructor()
	//;
	class_<OBDescriptor>("OBDescriptor")
	;
	function("FindType",&FindType);
	function("Predict",&Predict);
	function("GetStringValue",&GetStringValue);

	class_<OBConversion>("OBConversion")
	.constructor<std::istream*,std::ostream*>()
	;


}
