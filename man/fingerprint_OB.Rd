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
fingerprint_OB(obmolRefs, fingerprintName,reverse=FALSE)
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
	\item{reverse}{
		Reverse the order of bits in each fingerprint. The fingerprints
		returned by this function are in reverse order to those returned by
		obabel when using the same fingerprint calculation. As long as all
		fingerprings being compared came from the same source (ie, all from
		obabel, or all from ChemmineOB) the results should be correct. But if
		you need to mix fingerprints from different sources, care muse be
		taken to make sure they are all in the same direction. To make
		ChemmineOB match obabel output, set \code{reverse=TRUE}.

		We cannot set this to TRUE by default at this point as that could
		break code written prior to this change.
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
