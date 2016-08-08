#ifndef _PRINT_H_
#define _PRINT_H_

#define R_USE_C99_IN_CXX

#include <R_ext/Print.h>

#define Rc_printf Rprintf
#define Rc_fprintf_stdout(...) Rprintf(__VA_ARGS__)
#define Rc_fprintf_stderr(...) REprintf(__VA_ARGS__)
#define Rc_fputs_stdout(str) Rprintf(str)
#define Rc_fputs_stderr(str) REprintf(str)
#define Rc_fputc_stderr(ch) REprintf("%c", ch)

#endif
