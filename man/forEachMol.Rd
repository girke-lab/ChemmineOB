\name{forEachMol}
\alias{forEachMol}
\title{
   For Each Mol
}
\description{
   Reads in molecules from the given string in the given
   format and calls function \code{f} on each molecule.
   The results are then combinded using the \code{reduce}
   function, if given.
}
\usage{
forEachMol(inFormat, inString, f, reduce)
}
\arguments{
  \item{inFormat}{
	  Format of string in \code{source}. This can be any OpenBabel format such as "SDF" or "SMILES".
	  A full list can be found by executing "obabel -L formats".
	}
  \item{inString}{
	  The compounds to generate fingerprints for. The format should be exactly what would be in a file of 
	  the same format. Newlines can be represented with "\\n".
	}

  \item{f}{
     A function taking one OBMol reference and possibly
     returning a result.
   }

  \item{reduce}{
     This function will be passed to the \code{Reduce}
     function along with the results of all the \code{f}
     calls. This can be used to combine the results. 
  }
}
\value{
   The result will be a List of return values from the
   \code{f} function if not \code{reduce} function was
   given. Otherwise it will be the result of the
   \code{reduce} function applied to the results of the
   \code{f} function.
}
\author{
   Kevin Horan
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
	molRefs = forEachMol("SMILES","C1CCCCC1\ttest-compound-name",
                        identity,c)
}
}
