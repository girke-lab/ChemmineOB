\name{smartsSearch_OB}
\alias{smartsSearch_OB}
\title{
   SMARTS Search
}
\description{
   Returns the number of matches found for each
   compound given.
}
\usage{
smartsSearch_OB(obmolRefs, smartsPattern, uniqueMatches = TRUE)
}
\arguments{
  \item{obmolRefs}{
	  A list of OBMol references ( of class '_p_OpenBabel__OBMol') representing the molecules to
	  compute properites for. If you have your molecules in string format, you can create a list
	  of OBMol references using the \code{forEachMol} function, see the example.
   }
  \item{smartsPattern}{
      Any valid SMARTS pattern.
   }
  \item{uniqueMatches}{
     Should only unique matches be counted?
   }
}
\value{
   A vector of counts.
}
\author{
   Kevin Horan
}
\examples{
\dontrun{
	molRefs = forEachMol("SMILES","C1CCCCC1\ttest-compound-name",identity)
	smartsSearch_OB(molRefs,"[CH3X4]")
}
}
