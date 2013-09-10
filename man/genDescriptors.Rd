\name{genDescriptors}
\alias{genDescriptors}
\title{
	Generate Descriptors
}
\description{
	Generated the following descriptors in the given order:
	'abonds', 'atoms', 'bonds', 'dbonds', 'HBA1', 'HBA2', 'HBD', 'logP', 'MR', 'MW', 'nF', 'sbonds', 'tbonds', 'TPSA'.
}
\usage{
genDescriptors(from, source)
}
\arguments{
  \item{from}{
	  Format of string in \code{source}. This can be any OpenBabel format such as "SDF" or "SMILES".
	}
  \item{source}{
	  The compound to generate descriptors for. The format should be exactly what would be in a file of 
	  the same format. Newlines can be represented with "\\n".
	}
}
\value{
	Returns a vector of the following OpenBabel descriptors:
	'abonds', 'atoms', 'bonds', 'dbonds', 'HBA1', 'HBA2', 'HBD', 'logP', 'MR', 'MW', 'nF', 'sbonds', 'tbonds', 'TPSA'.
}
\author{
	Kevin Horan
}

\examples{
	genDescriptors("SMILES","C1CCCCC1\ttest-compound-name")
}
