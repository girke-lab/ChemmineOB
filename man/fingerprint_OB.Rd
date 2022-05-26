\name{fingerprint_OB}
\alias{fingerprint_OB}
\title{
	Fingerprints from OpenBabel
}
\description{
	Generates fingerprints using OpenBabel. The compound format
	can be specified as anything supported by OpenBabel.
	The fingerprint name can also be specified.
}
\usage{
fingerprint_OB(obmolRefs, fingerprintName)
}
\arguments{
  \item{obmolRefs}{
	  A list of OBMol references ( of class '_p_OpenBabel__OBMol') representing the molecules to
	  compute properites for. If you have your molecules in string format, you can create a list
	  of OBMol references using the \code{forEachMol} function, see the example.
	}

  \item{fingerprintName}{
	  The name of the fingerprint to generate. A list of available names can be found with 
	  "obabel -L fingerprints". Currently that list is: "FP2", "FP3", "FP4", and "MACCS". Currently "MACCS"
	  is not supported on windows.
	}
}
\value{
	A matrix of binary values is returned. There is a row for each compound. The length of a 
	row is determined by the fingerprint specified.
}
\author{
	Kevin Horan
}

\examples{
\dontrun{
	molRefs = forEachMol("SMILES","C1CCCCC1\ttest-compound-name",identity)
	fingerprint_OB(molRefs,"FP3")
}
}
