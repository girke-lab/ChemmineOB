\name{prop_OB}
\alias{prop_OB}
\title{
	Properties from OpenBabel
}
\description{
	Generates the following descriptors: 
  "cansmi", "cansmiNS", "formula", "HBA1", "HBA2", "HBD", "InChI", "InChIKey", 
  "logP", "MR", "MW", "nF","title", "TPSA".
}
\usage{
prop_OB(obmolRefs)
}
\arguments{
  \item{obmolRefs}{
	  A list of OBMol references ( of class '_p_OpenBabel__OBMol') representing the molecules to
	  compute properites for. If you have your molecules in string format, you can create a list
	  of OBMol references using the \code{forEachMol} function, see the example.
	}
}
\value{
	Returns a data frame with the following OpenBabel descriptors:
  "cansmi", "cansmiNS", "formula", "HBA1", "HBA2", "HBD", "InChI",  
  "logP", "MR", "MW", "nF", "title", "TPSA". The "InChI" descriptor is not supported on windows currently, so 
  that column will always be blank. 
}
\author{
	Kevin Horan
}

\examples{
	\dontrun{ # remove when ubuntu 16.04 bug fixed
		molRefs = forEachMol("SMILES","C1CCCCC1\ttest-compound-name",identity)
		prop_OB(molRefs)
	}
}
