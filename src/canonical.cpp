
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/graphsym.h>
#include <openbabel/canon.h>

#include <Rcpp.h>

using namespace Rcpp;
using namespace OpenBabel;

RcppExport SEXP canonicalReordering(SEXP pmolS){
   void *p = R_ExternalPtrAddr(pmolS);
   if(!p)
      Rf_error("bad pointer given, could not extract pointer");

   OBMol *pmol = reinterpret_cast<OBMol*>(p);
   
	if(!pmol)
      Rf_error("bad pointer given, could not cast to OBMol");

	std::vector<unsigned int> symmetry_classes;
	std::vector<unsigned int> canon_labels;

	OBGraphSym gs(pmol);
	gs.GetSymmetry(symmetry_classes);

	CanonicalLabels(pmol, symmetry_classes, canon_labels);

   NumericVector labels(canon_labels.begin(),canon_labels.end());


   return labels;
}
