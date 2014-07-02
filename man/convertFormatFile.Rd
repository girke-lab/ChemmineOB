\name{convertFormatFile}
\alias{convertFormatFile}
\title{
	Convert Format of Files
}
\description{
	Convert a file from one format to another
}
\usage{
convertFormatFile(from, to, fromFile, toFile,options=data.frame(names="gen2D",args=""))
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
 \item{options}{
      This is a data frame listing OpenBabel Options that should be
      applied while converting each compound. The "names" column
      should be the name of each option and the "args" column sould be
      the arguments to the corresponding option. For example, the
      default value shown above will cause 2D coordinates to be
      generated from the input compound and saved in the output
      compound, assuming the ouptut format supports it.

      A full list of available options can be found by executing
      "obabel -L ops" at the command line. The current list is:

      0xout    <file.xxx> Additional file output
      AddInIndex    Append input index to title
      AddPolarH    Adds hydrogen to polar atoms only
      canonical    Canonicalize the atom order
      energy    ForceField Energy Evaluation (not displayed in GUI)
      fillUC    <param> Fill the unit cell (strict or keepconnect)
      gen2D    Generate 2D coordinates
      gen3D    Generate 3D coordinates
      minimize    ForceField Energy Minimization (not displayed in GUI)
      partialcharge    <method> Calculate partial charges by specified method
      readconformer    Adjacent conformers combined into a single molecule
      s    Isomorphism filter(-s, -v options replacement)(not displayed in GUI)
      sort    <desc> Sort by descriptor(~desc for reverse)
      unique    [param] remove duplicates by descriptor;default inchi
      v    Isomorphism filter(-s, -v options replacement)(not displayed in GUI)

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
