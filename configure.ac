
AC_INIT([ChemmineOB], [1.0])

# Checks for libraries.
AC_ARG_WITH([openbabel-include],AS_HELP_STRING([--with-openbabel-include],
				[The location of OpenBabel header files]),
				[openbabel_include_path=$withval])

AS_IF([ test -n "$openbabel_include_path" ],
			[ OPENBABEL_CFLAGS="-I${openbabel_include_path}"],
		[ test -n  "$OPEN_BABEL_INCDIR" ],
			[ OPENBABEL_CFLAGS="-I${OPEN_BABEL_INCDIR}"],
		[
			AC_LANG_PUSH([C++])
				 PKG_CHECK_MODULES([OPENBABEL],[openbabel-3.0],[],[A=2])
			AC_LANG_POP([C++])
		]
)

AC_ARG_WITH([openbabel-lib],AS_HELP_STRING([--with-openbabel-lib],
				[The location of OpenBabel libraries]),
				[openbabel_lib_path=$withval])

AS_IF([ test -n "$openbabel_lib_path" ],
			[ OPENBABEL_LIBS="-L${openbabel_lib_path} -lopenbabel"],
		[ test -n  "$OPEN_BABEL_LIBDIR" ],
			[ OPENBABEL_LIBS="-L${OPEN_BABEL_LIBDIR} -lopenbabel"],
		[
			AC_LANG_PUSH([C++])
				 PKG_CHECK_MODULES([OPENBABEL],[openbabel-3.1],[],[OPENBABEL_LIBS="-lopenbabel"])
			AC_LANG_POP([C++])
		]
)



dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
	echo "could not determine R_HOME"
	exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_PROG_CC
AC_PROG_CPP


AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
