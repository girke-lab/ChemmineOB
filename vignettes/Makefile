PKG_NAME=ChemmineOB

all:
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit('${PKG_NAME}.Rmd')"
	pandoc --self-contained --output ${PKG_NAME}.html -H header.html --from markdown --smart --table-of-contents --template html.template ${PKG_NAME}.md;

clean:
	rm -f *.md
