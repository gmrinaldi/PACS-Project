#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP eval_FEM_fd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_mass_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_PDE_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_PDE_space_varying_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_FEM_stiff_matrix(SEXP, SEXP, SEXP, SEXP);
extern SEXP get_integration_points(SEXP, SEXP, SEXP, SEXP);
extern SEXP R_triangulate_native(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_Laplace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_PDE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regression_PDE_space_varying(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Smooth_FPCA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"eval_FEM_fd",                      (DL_FUNC) &eval_FEM_fd,                       8},
    {"get_FEM_mass_matrix",              (DL_FUNC) &get_FEM_mass_matrix,               4},
    {"get_FEM_PDE_matrix",               (DL_FUNC) &get_FEM_PDE_matrix,               17},
    {"get_FEM_PDE_space_varying_matrix", (DL_FUNC) &get_FEM_PDE_space_varying_matrix, 18},
    {"get_FEM_stiff_matrix",             (DL_FUNC) &get_FEM_stiff_matrix,              4},
    {"get_integration_points",           (DL_FUNC) &get_integration_points,            4},
    {"R_triangulate_native",             (DL_FUNC) &R_triangulate_native,              8},
    {"regression_Laplace",               (DL_FUNC) &regression_Laplace,               14},
    {"regression_PDE",                   (DL_FUNC) &regression_PDE,                   17},
    {"regression_PDE_space_varying",     (DL_FUNC) &regression_PDE_space_varying,     18},
    {"Smooth_FPCA",                      (DL_FUNC) &Smooth_FPCA,                      13},
    {NULL, NULL, 0}
};

void R_init_fdaPDE(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}