\name{convertFormatFile}
\alias{convertFormatFile}
\title{
	Convert Format of Files
}
\description{
	Convert a file from one format to another
}
\usage{
convertFormatFile(from, to, fromFile, toFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{
	  The format that \code{fromFile} is in. This should be a string supported by OpenBabel.
	}
  \item{to}{
	  The format to convert \code{toFile} to. 
	}
  \item{fromFile}{
	  The name of the file to be converted
	}
  \item{toFile}{
	  The name of the new file to be created or overwritten
	}
}
\value{
	No value is returned. \code{toFile} will be created with the compound in the new format.
}
\references{
OpenBabel  http://openbabel.org
}
\author{
	Kevin Horan
}

\seealso{
	\code{\link{convertFormat}}
}
\examples{
	\dontrun{
		convertFormatFile("SMI","SDF","test.smiles","test.sdf")
	}
}
