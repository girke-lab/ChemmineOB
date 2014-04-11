all:
#	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit2html('ChemmineOB.Rmd',options=c('toc',markdown::markdownHTMLOptions(TRUE)))"
#	pandoc -s ChemmineOB.html -o ChemmineOB.pdf
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit('ChemmineOB.Rmd')"
	 pandoc --output ChemmineOB.html -H header.html --from markdown --smart --table-of-contents ChemmineOB.md; 

clean:
	rm -f *.md
