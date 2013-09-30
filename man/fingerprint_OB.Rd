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
fingerprint_OB(format, source, fingerprintName)
}
\arguments{
  \item{format}{
	  Format of string in \code{source}. This can be any OpenBabel format such as "SDF" or "SMILES".
	  A full list can be found by executing "obabel -L formats".
	}
  \item{source}{
	  The compounds to generate fingerprints for. The format should be exactly what would be in a file of 
	  the same format. Newlines can be represented with "\\n".
	}
  \item{fingerprintName}{
	  The name of the fingerprint to generate. A list of available names can be found with 
	  "obabel -L fingerprints".
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

	fingerprint_OB("SMILES","C1CCCCC1\ttest-compound-name","FP3")
}
