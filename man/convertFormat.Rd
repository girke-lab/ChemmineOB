\name{convertFormat}
\alias{convertFormat}
\title{
	Convert Formats
}
\description{
	Converts compound data from one format to another. Some formats are not supported on windows.
	This currently includes "INCHI" and "INCHIKEY".
}
\usage{
convertFormat(from, to, source,options=data.frame(names="gen2D",args=""))
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
   #create an SDF from a SMILES string and put the atoms in the canonical order
\dontrun{
	sdfStr = convertFormat("SMI","SDF","CC(=O)OC1=CC=CC=C1C(=O)O\ttest_name",
                          options=data.frame(names="canonical",args=""))
}
}
