\name{exactMass_OB}
\alias{exactMass_OB}
\title{
   Exact Mass (Monoisotopic Mass)
}
\description{
   Computes the exact mass of each compound given.
}
\usage{
exactMass_OB(obmolRefs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obmolRefs}{
	  A list of OBMol references ( of class '_p_OpenBabel__OBMol') representing the compounds.
   }
}
\value{
   A vector of mass values.
}
\author{
   Kevin Horan
}

\examples{
\dontrun{
	molRefs = forEachMol("SMILES","C1CCCCC1\ttest-compound-name",identity)
	exactMass_OB(molRefs)
}
}
\keyword{monoisotopic}
\keyword{mass}
\keyword{exact}
