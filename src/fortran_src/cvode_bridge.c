#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sundials/sundials_context.h>

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int (*uclchem_rhs_fn)(double t, const double *y, double *ydot, void *user_data);
extern int uclchem_cvode_dense_jacobian(double t, void *y_ptr, void *data_ptr);
extern int uclchem_sparse_jacobian_nnz(void);
extern void uclchem_sparse_jacobian_pattern(int *colptr, int *rowind);
extern int uclchem_cvode_sparse_jacobian(double t, void *y_ptr, void *data_ptr);

struct uclchem_callback_data {
  uclchem_rhs_fn rhs;
  void *user_data;
};

struct uclchem_cvode_solver {
  void *cvode_mem;
  SUNContext sunctx;
  N_Vector yvec;
  N_Vector abstol_vec;
  N_Vector precond_diag;
  SUNMatrix A;
  SUNLinearSolver LS;
  struct uclchem_callback_data callbacks;
  sunindextype neq;
  sunindextype nnz;
  int linear_solver_kind;
  int initialized;
  int sparse_jac_current;
  realtype cached_gamma;
};

static struct uclchem_cvode_solver solver = {0};
enum {
  UCLCHEM_LINEAR_SPARSE_KLU = 0,
  UCLCHEM_LINEAR_DENSE = 1,
  UCLCHEM_LINEAR_SPGMR = 2
};

static void uclchem_append_cvode_log(const char *fmt, ...) {
  FILE *fp;
  va_list args;

  fp = fopen("cvode_c_debug.log", "a");
  if (fp == NULL) {
    return;
  }

  va_start(args, fmt);
  vfprintf(fp, fmt, args);
  va_end(args);
  fputc('\n', fp);
  fclose(fp);
}

static void uclchem_log_solver_diagnostics(const char *stage, int cvode_flag) {
  long int nst = 0;
  long int nfe = 0;
  long int netf = 0;
  long int nni = 0;
  long int nncf = 0;
  long int nje = 0;
  long int nfeLS = 0;
  long int nli = 0;
  long int nlcf = 0;
  long int npe = 0;
  long int nps = 0;
  long int njt = 0;
  long int njv = 0;
  long int last_lin_flag = 0;
  int qcur = 0;
  realtype tcur = 0.0;
  char *lin_flag_name = NULL;

  if (solver.cvode_mem == NULL) {
    return;
  }

  (void)CVodeGetNumSteps(solver.cvode_mem, &nst);
  (void)CVodeGetNumRhsEvals(solver.cvode_mem, &nfe);
  (void)CVodeGetNumErrTestFails(solver.cvode_mem, &netf);
  (void)CVodeGetNumNonlinSolvIters(solver.cvode_mem, &nni);
  (void)CVodeGetNumNonlinSolvConvFails(solver.cvode_mem, &nncf);
  (void)CVodeGetCurrentOrder(solver.cvode_mem, &qcur);
  (void)CVodeGetCurrentTime(solver.cvode_mem, &tcur);
  (void)CVodeGetLinSolveStats(solver.cvode_mem, &nje, &nfeLS, &nli, &nlcf, &npe, &nps, &njt, &njv);
  (void)CVodeGetLastLinFlag(solver.cvode_mem, &last_lin_flag);

  lin_flag_name = CVodeGetLinReturnFlagName(last_lin_flag);
  uclchem_append_cvode_log("C diagnostics: stage=%s cvode_flag=%d tcur=%.17g qcur=%d nst=%ld nfe=%ld netf=%ld nni=%ld nncf=%ld nje=%ld nfeLS=%ld nli=%ld nlcf=%ld npe=%ld nps=%ld njt=%ld njv=%ld last_lin_flag=%ld last_lin_flag_name=%s",
                           stage, cvode_flag, (double)tcur, qcur, nst, nfe, netf, nni, nncf, nje,
                           nfeLS, nli, nlcf, npe, nps, njt, njv, last_lin_flag,
                           (lin_flag_name != NULL) ? lin_flag_name : "(null)");

  if (solver.linear_solver_kind == UCLCHEM_LINEAR_SPARSE_KLU && solver.LS != NULL) {
    sunindextype klu_last_flag = SUNLinSolLastFlag_KLU(solver.LS);
    sun_klu_common *common = SUNLinSol_KLUGetCommon(solver.LS);
    uclchem_append_cvode_log("KLU diagnostics: last_flag=%ld", (long)klu_last_flag);
    if (common != NULL) {
      uclchem_append_cvode_log("KLU common: status=%d nrealloc=%d structural_rank=%d numerical_rank=%d singular_col=%d noffdiag=%d rcond=%.17g condest=%.17g rgrowth=%.17g",
                               common->status, common->nrealloc, common->structural_rank,
                               common->numerical_rank, common->singular_col, common->noffdiag,
                               common->rcond, common->condest, common->rgrowth);
    }
  } else if (solver.linear_solver_kind == UCLCHEM_LINEAR_SPGMR && solver.LS != NULL) {
    sunindextype spgmr_last_flag = SUNLinSolLastFlag_SPGMR(solver.LS);
    uclchem_append_cvode_log("SPGMR diagnostics: last_flag=%ld", (long)spgmr_last_flag);
  }
}

static int rhs_bridge(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  struct uclchem_callback_data *data = (struct uclchem_callback_data *)user_data;
  double *y_ptr = N_VGetArrayPointer_Serial(y);
  double *ydot_ptr = N_VGetArrayPointer_Serial(ydot);

  if (data == NULL || data->rhs == NULL || y_ptr == NULL || ydot_ptr == NULL) {
    return CV_RHSFUNC_FAIL;
  }

  return data->rhs((double)t, y_ptr, ydot_ptr, data->user_data);
}

static int jac_sparse_bridge(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  double *y_ptr = N_VGetArrayPointer_Serial(y);
  realtype *jac_data = SM_DATA_S(J);

  (void)fy;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;

  if (y_ptr == NULL || jac_data == NULL) {
    return CVLS_JACFUNC_UNRECVR;
  }

  return uclchem_cvode_sparse_jacobian((double)t, (void *)y_ptr, (void *)jac_data);
}

static int jac_dense_bridge(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  double *y_ptr = N_VGetArrayPointer_Serial(y);
  realtype *jac_data = SM_DATA_D(J);

  (void)fy;
  (void)user_data;
  (void)tmp1;
  (void)tmp2;
  (void)tmp3;

  if (y_ptr == NULL || jac_data == NULL) {
    return CVLS_JACFUNC_UNRECVR;
  }

  return uclchem_cvode_dense_jacobian((double)t, (void *)y_ptr, (void *)jac_data);
}

static const char *uclchem_linear_solver_name(int linear_solver_kind) {
  switch (linear_solver_kind) {
  case UCLCHEM_LINEAR_DENSE:
    return "dense";
  case UCLCHEM_LINEAR_SPGMR:
    return "spgmr";
  case UCLCHEM_LINEAR_SPARSE_KLU:
  default:
    return "sparse_klu";
  }
}

static int uclchem_solver_uses_sparse_jacobian(int linear_solver_kind) {
  return linear_solver_kind == UCLCHEM_LINEAR_SPARSE_KLU ||
         linear_solver_kind == UCLCHEM_LINEAR_SPGMR;
}

static int uclchem_linear_solver_kind_from_env(void) {
  const char *value = getenv("UCLCHEM_CVODE_LINEAR_SOLVER");
  if (value == NULL) {
    return UCLCHEM_LINEAR_SPARSE_KLU;
  }
  if (strcmp(value, "dense") == 0 || strcmp(value, "DENSE") == 0) {
    return UCLCHEM_LINEAR_DENSE;
  }
  if (strcmp(value, "spgmr") == 0 || strcmp(value, "SPGMR") == 0 ||
      strcmp(value, "gmres") == 0 || strcmp(value, "GMRES") == 0 ||
      strcmp(value, "iterative") == 0 || strcmp(value, "ITERATIVE") == 0) {
    return UCLCHEM_LINEAR_SPGMR;
  }
  return UCLCHEM_LINEAR_SPARSE_KLU;
}

static int uclchem_env_int(const char *name, int default_value) {
  const char *value = getenv(name);
  char *endptr = NULL;
  long parsed;

  if (value == NULL || *value == '\0') {
    return default_value;
  }

  parsed = strtol(value, &endptr, 10);
  if (endptr == value || (endptr != NULL && *endptr != '\0')) {
    return default_value;
  }

  return (int)parsed;
}

static double uclchem_env_double(const char *name, double default_value) {
  const char *value = getenv(name);
  char *endptr = NULL;
  double parsed;

  if (value == NULL || *value == '\0') {
    return default_value;
  }

  parsed = strtod(value, &endptr);
  if (endptr == value || (endptr != NULL && *endptr != '\0')) {
    return default_value;
  }

  return parsed;
}

static int uclchem_env_bool(const char *name, int default_value) {
  const char *value = getenv(name);

  if (value == NULL || *value == '\0') {
    return default_value;
  }

  if (strcmp(value, "1") == 0 || strcmp(value, "true") == 0 || strcmp(value, "TRUE") == 0 ||
      strcmp(value, "yes") == 0 || strcmp(value, "YES") == 0 || strcmp(value, "on") == 0 ||
      strcmp(value, "ON") == 0) {
    return 1;
  }

  if (strcmp(value, "0") == 0 || strcmp(value, "false") == 0 || strcmp(value, "FALSE") == 0 ||
      strcmp(value, "no") == 0 || strcmp(value, "NO") == 0 || strcmp(value, "off") == 0 ||
      strcmp(value, "OFF") == 0) {
    return 0;
  }

  return default_value;
}

static void uclchem_configure_klu(SUNLinearSolver LS) {
  int retval;
  int ordering;
  sun_klu_common *common;

  if (LS == NULL) {
    return;
  }

  ordering = uclchem_env_int("UCLCHEM_KLU_ORDERING", SUNKLU_ORDERING_DEFAULT);
  retval = SUNLinSol_KLUSetOrdering(LS, ordering);
  if (retval != SUNLS_SUCCESS) {
    uclchem_append_cvode_log("KLU config: failed to set ordering=%d retval=%d", ordering, retval);
    return;
  }

  common = SUNLinSol_KLUGetCommon(LS);
  if (common == NULL) {
    return;
  }

  common->scale = uclchem_env_int("UCLCHEM_KLU_SCALE", common->scale);
  common->btf = uclchem_env_int("UCLCHEM_KLU_BTF", common->btf);
  common->tol = uclchem_env_double("UCLCHEM_KLU_TOL", common->tol);

  uclchem_append_cvode_log("KLU config: ordering=%d scale=%d btf=%d tol=%.17g",
                           ordering, common->scale, common->btf, common->tol);
}

static void uclchem_configure_spgmr(SUNLinearSolver LS) {
  int retval;
  int max_restarts;
  int gstype;

  if (LS == NULL) {
    return;
  }

  max_restarts = uclchem_env_int("UCLCHEM_SPGMR_MAX_RESTARTS", 2);
  if (max_restarts >= 0) {
    retval = SUNLinSol_SPGMRSetMaxRestarts(LS, max_restarts);
    if (retval != SUNLS_SUCCESS) {
      uclchem_append_cvode_log("SPGMR config: failed to set max_restarts=%d retval=%d",
                               max_restarts, retval);
    }
  }

  gstype = uclchem_env_int("UCLCHEM_SPGMR_GSTYPE", SUNSPGMR_GSTYPE_DEFAULT);
  retval = SUNLinSol_SPGMRSetGSType(LS, gstype);
  if (retval != SUNLS_SUCCESS) {
    uclchem_append_cvode_log("SPGMR config: failed to set gstype=%d retval=%d", gstype, retval);
  }

  uclchem_append_cvode_log("SPGMR config: maxl=%d max_restarts=%d gstype=%d",
                           uclchem_env_int("UCLCHEM_SPGMR_MAXL", 30), max_restarts, gstype);
}

static int uclchem_fill_sparse_jacobian_values(realtype t, N_Vector y) {
  double *y_ptr;
  realtype *jac_data;

  if (solver.A == NULL) {
    return CVLS_MEM_NULL;
  }

  y_ptr = N_VGetArrayPointer_Serial(y);
  jac_data = SM_DATA_S(solver.A);
  if (y_ptr == NULL || jac_data == NULL) {
    return CVLS_JACFUNC_UNRECVR;
  }

  return uclchem_cvode_sparse_jacobian((double)t, (void *)y_ptr, (void *)jac_data);
}

static int uclchem_sparse_matvec(N_Vector v, N_Vector Jv) {
  realtype *jv_data;
  realtype *v_data;
  realtype *jac_data;
  sunindextype *index_ptrs;
  sunindextype *index_vals;
  sunindextype col;
  sunindextype idx;

  if (solver.A == NULL) {
    return CVLS_MEM_NULL;
  }

  jv_data = N_VGetArrayPointer_Serial(Jv);
  v_data = N_VGetArrayPointer_Serial(v);
  jac_data = SM_DATA_S(solver.A);
  index_ptrs = SM_INDEXPTRS_S(solver.A);
  index_vals = SM_INDEXVALS_S(solver.A);
  if (jv_data == NULL || v_data == NULL || jac_data == NULL || index_ptrs == NULL ||
      index_vals == NULL) {
    return CVLS_MEM_FAIL;
  }

  for (col = 0; col < solver.neq; ++col) {
    jv_data[col] = 0.0;
  }

  for (col = 0; col < solver.neq; ++col) {
    realtype vcol = v_data[col];
    if (vcol == 0.0) {
      continue;
    }

    for (idx = index_ptrs[col]; idx < index_ptrs[col + 1]; ++idx) {
      jv_data[index_vals[idx]] += jac_data[idx] * vcol;
    }
  }

  return 0;
}

static int uclchem_build_diag_preconditioner(realtype gamma) {
  realtype *diag_inv;
  realtype *jac_data;
  sunindextype *index_ptrs;
  sunindextype *index_vals;
  realtype diag_floor;
  sunindextype col;
  sunindextype idx;

  if (solver.A == NULL || solver.precond_diag == NULL) {
    return CVLS_MEM_NULL;
  }

  diag_inv = N_VGetArrayPointer_Serial(solver.precond_diag);
  jac_data = SM_DATA_S(solver.A);
  index_ptrs = SM_INDEXPTRS_S(solver.A);
  index_vals = SM_INDEXVALS_S(solver.A);
  if (diag_inv == NULL || jac_data == NULL || index_ptrs == NULL || index_vals == NULL) {
    return CVLS_MEM_FAIL;
  }

  diag_floor = uclchem_env_double("UCLCHEM_DIAG_PRECOND_FLOOR", 1.0e-30);

  for (col = 0; col < solver.neq; ++col) {
    realtype diag = 1.0;
    for (idx = index_ptrs[col]; idx < index_ptrs[col + 1]; ++idx) {
      if (index_vals[idx] == col) {
        diag -= gamma * jac_data[idx];
        break;
      }
    }

    if (fabs(diag) <= diag_floor) {
      diag_inv[col] = 1.0;
    } else {
      diag_inv[col] = 1.0 / diag;
    }
  }

  solver.cached_gamma = gamma;
  return 0;
}

static int jtimes_setup_bridge(realtype t, N_Vector y, N_Vector fy, void *user_data) {
  int retval;

  (void)fy;
  (void)user_data;

  retval = uclchem_fill_sparse_jacobian_values(t, y);
  if (retval == 0) {
    solver.sparse_jac_current = 1;
  }
  return retval;
}

static int jtimes_bridge(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy,
                         void *user_data, N_Vector tmp) {
  int retval;

  (void)fy;
  (void)user_data;
  (void)tmp;

  if (!solver.sparse_jac_current) {
    retval = uclchem_fill_sparse_jacobian_values(t, y);
    if (retval != 0) {
      return retval;
    }
    solver.sparse_jac_current = 1;
  }

  return uclchem_sparse_matvec(v, Jv);
}

static int prec_setup_bridge(realtype t, N_Vector y, N_Vector fy, booleantype jok,
                             booleantype *jcurPtr, realtype gamma, void *user_data) {
  int retval;

  (void)fy;
  (void)user_data;

  if (jcurPtr == NULL) {
    return CVLS_ILL_INPUT;
  }

  if (!jok || !solver.sparse_jac_current) {
    retval = uclchem_fill_sparse_jacobian_values(t, y);
    if (retval != 0) {
      return retval;
    }
    solver.sparse_jac_current = 1;
    *jcurPtr = SUNTRUE;
  } else {
    *jcurPtr = SUNFALSE;
  }

  return uclchem_build_diag_preconditioner(gamma);
}

static int prec_solve_bridge(realtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                             realtype gamma, realtype delta, int lr, void *user_data) {
  realtype *r_data;
  realtype *z_data;
  realtype *diag_inv;
  sunindextype i;

  (void)t;
  (void)y;
  (void)fy;
  (void)gamma;
  (void)delta;
  (void)lr;
  (void)user_data;

  if (solver.precond_diag == NULL) {
    return CVLS_PMEM_NULL;
  }

  r_data = N_VGetArrayPointer_Serial(r);
  z_data = N_VGetArrayPointer_Serial(z);
  diag_inv = N_VGetArrayPointer_Serial(solver.precond_diag);
  if (r_data == NULL || z_data == NULL || diag_inv == NULL) {
    return CVLS_MEM_FAIL;
  }

  for (i = 0; i < solver.neq; ++i) {
    z_data[i] = diag_inv[i] * r_data[i];
  }

  return 0;
}

static void uclchem_cvode_free_solver(void) {
  if (solver.cvode_mem != NULL) {
    CVodeFree(&solver.cvode_mem);
  }
  if (solver.LS != NULL) {
    SUNLinSolFree(solver.LS);
  }
  if (solver.A != NULL) {
    SUNMatDestroy(solver.A);
  }
  if (solver.precond_diag != NULL) {
    N_VDestroy(solver.precond_diag);
  }
  if (solver.abstol_vec != NULL) {
    N_VDestroy(solver.abstol_vec);
  }
  if (solver.yvec != NULL) {
    N_VDestroy(solver.yvec);
  }
  if (solver.sunctx != NULL) {
    SUNContext_Free(&solver.sunctx);
  }
  memset(&solver, 0, sizeof(solver));
}

int uclchem_cvode_init(int neq, uclchem_rhs_fn rhs) {
  int retval;
  int max_order;
  sunindextype *index_ptrs;
  sunindextype *index_vals;
  int nnz;

  if (neq <= 0 || rhs == NULL) {
    return CV_ILL_INPUT;
  }

  uclchem_cvode_free_solver();

  retval = SUNContext_Create(NULL, &solver.sunctx);
  if (retval != 0) {
    return CV_CONTEXT_ERR;
  }

  solver.neq = (sunindextype)neq;
  solver.linear_solver_kind = uclchem_linear_solver_kind_from_env();
  uclchem_append_cvode_log("CVODE init: linear_solver=%s", uclchem_linear_solver_name(solver.linear_solver_kind));

  if (uclchem_solver_uses_sparse_jacobian(solver.linear_solver_kind)) {
    nnz = uclchem_sparse_jacobian_nnz();
    if (nnz <= 0) {
      uclchem_cvode_free_solver();
      return CV_ILL_INPUT;
    }
    solver.nnz = (sunindextype)nnz;
  } else {
    solver.nnz = 0;
  }

  solver.yvec = N_VNew_Serial(solver.neq, solver.sunctx);
  solver.abstol_vec = N_VNew_Serial(solver.neq, solver.sunctx);
  if (solver.yvec == NULL || solver.abstol_vec == NULL) {
    uclchem_cvode_free_solver();
    return CV_MEM_FAIL;
  }

  if (uclchem_solver_uses_sparse_jacobian(solver.linear_solver_kind)) {
    solver.A = SUNSparseMatrix(solver.neq, solver.neq, solver.nnz, CSC_MAT, solver.sunctx);
    if (solver.A == NULL) {
      uclchem_cvode_free_solver();
      return CV_MEM_FAIL;
    }

    index_ptrs = SM_INDEXPTRS_S(solver.A);
    index_vals = SM_INDEXVALS_S(solver.A);
    if (index_ptrs == NULL || index_vals == NULL) {
      uclchem_cvode_free_solver();
      return CV_MEM_FAIL;
    }
    uclchem_sparse_jacobian_pattern((int *)index_ptrs, (int *)index_vals);

    if (solver.linear_solver_kind == UCLCHEM_LINEAR_SPARSE_KLU) {
      solver.LS = SUNLinSol_KLU(solver.yvec, solver.A, solver.sunctx);
      if (solver.LS == NULL) {
        uclchem_cvode_free_solver();
        return CV_MEM_FAIL;
      }

      uclchem_configure_klu(solver.LS);
    } else {
      int maxl = uclchem_env_int("UCLCHEM_SPGMR_MAXL", 30);
      if (maxl <= 0) {
        maxl = 30;
      }

      solver.LS = SUNLinSol_SPGMR(solver.yvec, SUN_PREC_LEFT, maxl, solver.sunctx);
      solver.precond_diag = N_VClone(solver.yvec);
      if (solver.LS == NULL || solver.precond_diag == NULL) {
        uclchem_cvode_free_solver();
        return CV_MEM_FAIL;
      }

      uclchem_configure_spgmr(solver.LS);
    }
  } else {
    solver.A = SUNDenseMatrix(solver.neq, solver.neq, solver.sunctx);
    if (solver.A == NULL) {
      uclchem_cvode_free_solver();
      return CV_MEM_FAIL;
    }
    solver.LS = SUNLinSol_Dense(solver.yvec, solver.A, solver.sunctx);
    if (solver.LS == NULL) {
      uclchem_cvode_free_solver();
      return CV_MEM_FAIL;
    }
  }

  solver.cvode_mem = CVodeCreate(CV_BDF, solver.sunctx);
  if (solver.cvode_mem == NULL) {
    uclchem_cvode_free_solver();
    return CV_MEM_FAIL;
  }

  max_order = uclchem_env_int("UCLCHEM_CVODE_MAX_ORDER", 5);
  if (max_order > 0) {
    retval = CVodeSetMaxOrd(solver.cvode_mem, max_order);
    if (retval != CV_SUCCESS) {
      uclchem_cvode_free_solver();
      return retval;
    }
    uclchem_append_cvode_log("CVODE config: max_order=%d", max_order);
  }

  solver.callbacks.rhs = rhs;
  solver.callbacks.user_data = NULL;
  solver.initialized = 0;
  solver.sparse_jac_current = 0;
  solver.cached_gamma = 0.0;

  return CV_SUCCESS;
}

int uclchem_cvode_integrate(double *y, double *t, double tout, double reltol, const double *abstol,
                            long int mxstep) {
  int retval;
  sunindextype i;
  double *y_data;
  double *abstol_data;

  if (solver.cvode_mem == NULL || solver.yvec == NULL || solver.abstol_vec == NULL ||
      solver.LS == NULL || (solver.linear_solver_kind != UCLCHEM_LINEAR_SPGMR && solver.A == NULL)) {
    return CV_NO_MALLOC;
  }
  if (y == NULL || t == NULL || abstol == NULL) {
    return CV_ILL_INPUT;
  }

  y_data = N_VGetArrayPointer_Serial(solver.yvec);
  abstol_data = N_VGetArrayPointer_Serial(solver.abstol_vec);
  if (y_data == NULL || abstol_data == NULL) {
    return CV_MEM_FAIL;
  }

  for (i = 0; i < solver.neq; ++i) {
    y_data[i] = y[i];
    abstol_data[i] = abstol[i];
  }
  solver.sparse_jac_current = 0;

  if (!solver.initialized) {
    retval = CVodeInit(solver.cvode_mem, rhs_bridge, (realtype)(*t), solver.yvec);
    if (retval != CV_SUCCESS) {
      return retval;
    }

    retval = CVodeSetUserData(solver.cvode_mem, &solver.callbacks);
    if (retval != CV_SUCCESS) {
      return retval;
    }

    retval = CVodeSetLinearSolver(
        solver.cvode_mem, solver.LS,
        (solver.linear_solver_kind == UCLCHEM_LINEAR_SPGMR) ? NULL : solver.A);
    if (retval != CV_SUCCESS) {
      return retval;
    }

    if (solver.linear_solver_kind == UCLCHEM_LINEAR_SPGMR) {
      retval = CVodeSetJacTimes(solver.cvode_mem, jtimes_setup_bridge, jtimes_bridge);
      if (retval != CV_SUCCESS) {
        return retval;
      }

      retval = CVodeSetPreconditioner(solver.cvode_mem, prec_setup_bridge, prec_solve_bridge);
      if (retval != CV_SUCCESS) {
        return retval;
      }
    } else {
      retval = CVodeSetJacFn(
          solver.cvode_mem,
          (solver.linear_solver_kind == UCLCHEM_LINEAR_SPARSE_KLU) ? jac_sparse_bridge : jac_dense_bridge);
      if (retval != CV_SUCCESS) {
        return retval;
      }
    }

    {
      int msbj = uclchem_env_int("UCLCHEM_CVODE_JAC_EVAL_FREQUENCY", -1);
      int linscale = uclchem_env_bool("UCLCHEM_CVODE_LINEAR_SOLUTION_SCALING", -1);
      double epslin = uclchem_env_double("UCLCHEM_CVODE_EPS_LIN", -1.0);

      if (msbj > 0) {
        retval = CVodeSetJacEvalFrequency(solver.cvode_mem, (long int)msbj);
        if (retval != CV_SUCCESS) {
          return retval;
        }
        uclchem_append_cvode_log("CVODE config: jac_eval_frequency=%d", msbj);
      }

      if (linscale >= 0) {
        retval = CVodeSetLinearSolutionScaling(solver.cvode_mem, linscale ? SUNTRUE : SUNFALSE);
        if (retval != CV_SUCCESS) {
          return retval;
        }
        uclchem_append_cvode_log("CVODE config: linear_solution_scaling=%d", linscale);
      }

      if (epslin > 0.0) {
        retval = CVodeSetEpsLin(solver.cvode_mem, (realtype)epslin);
        if (retval != CV_SUCCESS) {
          return retval;
        }
        uclchem_append_cvode_log("CVODE config: eps_lin=%.17g", epslin);
      }
    }

    solver.initialized = 1;
  } else {
    retval = CVodeReInit(solver.cvode_mem, (realtype)(*t), solver.yvec);
    if (retval != CV_SUCCESS) {
      return retval;
    }
  }

  retval = CVodeSVtolerances(solver.cvode_mem, (realtype)reltol, solver.abstol_vec);
  if (retval != CV_SUCCESS) {
    return retval;
  }

  if (mxstep > 0) {
    retval = CVodeSetMaxNumSteps(solver.cvode_mem, mxstep);
    if (retval != CV_SUCCESS) {
      return retval;
    }
  }

  retval = CVode(solver.cvode_mem, (realtype)tout, solver.yvec, (realtype *)t, CV_NORMAL);
  for (i = 0; i < solver.neq; ++i) {
    y[i] = y_data[i];
  }

  if (retval != CV_SUCCESS && retval != CV_TSTOP_RETURN) {
    uclchem_log_solver_diagnostics("CVode-return", retval);
  }

  return retval;
}

int uclchem_cvode_finalize(void) {
  uclchem_cvode_free_solver();
  return CV_SUCCESS;
}
