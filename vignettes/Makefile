all:
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit2html('ChemmineOB.Rmd',options=c('toc',markdown::markdownHTMLOptions(TRUE)))"
#	pandoc -s ChemmineOB.html -o ChemmineOB.pdf
clean:
	rm -f *.md
