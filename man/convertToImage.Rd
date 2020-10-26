\name{convertToImage}
\alias{convertToImage}
\title{
	Convert To Image
}
\description{
	Converts compound data to an image format
}
\usage{

convertToImage(from,to,source,toFile,options=data.frame(names="gen2D",args=""),
									out_options=data.frame(names=c("p"),args=c(300)))
}
\arguments{
  \item{from}{
	  The format that \code{source} is in. This should be a string supported by OpenBabel.
	}
  \item{to}{
	  The  image format to convert \code{source} to. This should be either "PNG", "POV", or "SVG".
	}
  \item{source}{
	  The initial compound format, as a string. The format of the string should be identical to
	  the file format of the same name. Tabs and newlines may be represented with \\t and \\n,
	  repsectivly.
	}
	\item{toFile}{
		The name of the output file.
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
	\item{out_options}{
		Similar to the \code{options} argument, but these options apply to the output format.
		
		For the "PNG" format, possible options are: 

		"p" for the number of pixels on the square output file.

		For the "SVG" format, see the OpenBabel docs for a full list (https://openbabel.org/docs/dev/FileFormats/SVG_2D_depiction.html).
	}
}
\value{
	Generates a square image with the given compounds arranged in a grid format. The width (and height, which must be the same)
	can be given using the \code{out_options} with argument "p". Writes the resulting image to a file specified by \code{to_file}.
}
\references{
OpenBabel  http://openbabel.org
}
\author{
	Kevin Horan
}


\examples{
   # Convert the given smiles string to a PNG image with size 500x500.
	\dontrun{
		sdfStr = convertToImage("SMI","PNG","CC(=O)OC1=CC=CC=C1C(=O)O\ttest_name","out_image.png",
                          out_options=data.frame(names="p",args="500"))
	}
}
