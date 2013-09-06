\name{convertFormat}
\alias{convertFormat}
\title{
	Convert Formats
}
\description{
	Converts compound data from one format to another.
}
\usage{
convertFormat(from, to, source)
}
\arguments{
  \item{from}{
	  The format that \code{source} is in. This should be a string supported by OpenBabel.
	}
  \item{to}{
	  The format to convert \code{source} to. 
	}
  \item{source}{
	  The initial compound format, as a string. The format of the string should be identical to
	  the file format of the same name. Tabs and newlines may be represented with \\t and \\n,
	  repsectivly.
	}
}
\value{
	Returns the compound given in \code{source} in the format specified by \code{to}.
}
\references{
OpenBabel  http://openbabel.org
}
\author{
	Kevin Horan
}


\seealso{
	\code{\link{convertFormatFile}}
}
\examples{
	sdfStr = convertFormat("SMI","SDF","CC(=O)OC1=CC=CC=C1C(=O)O\ttest_name")
}
