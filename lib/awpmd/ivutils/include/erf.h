# ifndef ERF_H
# define ERF_H

# if defined(_WIN32) && !defined(HAVE_ERF)

# ifdef __cplusplus
extern "C" {
# endif

double erf(double x);
double erfc(double x);

# ifdef __cplusplus
}
# endif

# endif

# endif
