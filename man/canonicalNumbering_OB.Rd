\name{canonicalNumbering_OB}
\alias{canonicalNumbering_OB}
\title{
   Canonical Numbering OB
}
\description{
   Computes a re-arrangement required to transform the atom numbering
   of the given compound into the canonical atom numbering. This
   function uses the OBGraphSym and CanonicalLabels classes of Open
   Babel to compute the re-arrangement.
}
\usage{
canonicalNumbering_OB(obmolRefs)
}
\arguments{
  \item{obmolRefs}{
	  A list of OBMol references ( of class '_p_OpenBabel__OBMol') representing the compounds.
   }
}
\value{
   A list of vectors of index values. Each item in the list corresponds to one of the given
 compounds. The values of a list item are the re-arrangement of the atoms.
 For example, if the value in item 1, column 1 is  25, that means that
 atom number 1 in the original compound should become atom number 25
 in the canonical version of that compound.
}
\author{
   Kevin Horan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\references{
   http://openbabel.org/api/2.3/canonical_code_algorithm.shtml
}

\seealso{
   \code{\link{convertFormat}} to return a new compound in the canonical format.

}
\examples{
\dontrun{
   library(ChemmineR)
   data(sdfsample)
   labels = canonicalNumbering_OB(obmol(sdfsample[[1]]))
}

}
\keyword{canonical}
\keyword{morgan}
