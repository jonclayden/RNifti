#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([RNifti], [1.4.0], [code@clayden.org])

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
ZLIB_FLAGS="-Izlib"
ZLIB_OBJECTS=`ls -1 src/zlib/*.c | sed 's/^src\///' | sed 's/\\.c$/\\.o/' | tr '\n' ' '`
AC_SEARCH_LIBS([inflate], [z], [
    ZLIB_FLAGS=""
    ZLIB_OBJECTS=""
])

AC_SUBST(ZLIB_FLAGS)
AC_SUBST(ZLIB_OBJECTS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
