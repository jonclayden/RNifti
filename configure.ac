#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([RNifti], [1.6.0], [code@clayden.org])

# Get compiler flags from R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

# Checks for libraries.
AC_SEARCH_LIBS([inflate], [z], [], [
    AC_CONFIG_COMMANDS([zlib], [cd src/zlib && ./configure --static])
    ZLIB_OBJECT="zlib/libz.a"
    ZLIB_LIBPATH="-Lzlib"
])

AC_SUBST(ZLIB_OBJECT)
AC_SUBST(ZLIB_LIBPATH)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
