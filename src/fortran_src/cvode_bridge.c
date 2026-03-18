#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sundials/sundials_context.h>

#include <string.h>

typedef int (*uclchem_rhs_fn)(double t, const double *y, double *ydot, void *user_data);
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
  SUNMatrix A;
  SUNLinearSolver LS;
  struct uclchem_callback_data callbacks;
  sunindextype neq;
  sunindextype nnz;
  int initialized;
};

static struct uclchem_cvode_solver solver = {0};

static int rhs_bridge(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  struct uclchem_callback_data *data = (struct uclchem_callback_data *)user_data;
  double *y_ptr = N_VGetArrayPointer_Serial(y);
  double *ydot_ptr = N_VGetArrayPointer_Serial(ydot);

  if (data == NULL || data->rhs == NULL || y_ptr == NULL || ydot_ptr == NULL) {
    return CV_RHSFUNC_FAIL;
  }

  return data->rhs((double)t, y_ptr, ydot_ptr, data->user_data);
}

static int jac_bridge(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
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
  nnz = uclchem_sparse_jacobian_nnz();
  if (nnz <= 0) {
    uclchem_cvode_free_solver();
    return CV_ILL_INPUT;
  }
  solver.nnz = (sunindextype)nnz;

  solver.yvec = N_VNew_Serial(solver.neq, solver.sunctx);
  solver.abstol_vec = N_VNew_Serial(solver.neq, solver.sunctx);
  if (solver.yvec == NULL || solver.abstol_vec == NULL) {
    uclchem_cvode_free_solver();
    return CV_MEM_FAIL;
  }

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

  solver.LS = SUNLinSol_KLU(solver.yvec, solver.A, solver.sunctx);
  if (solver.LS == NULL) {
    uclchem_cvode_free_solver();
    return CV_MEM_FAIL;
  }

  retval = SUNLinSol_KLUSetOrdering(solver.LS, SUNKLU_ORDERING_DEFAULT);
  if (retval != SUNLS_SUCCESS) {
    uclchem_cvode_free_solver();
    return CV_LINIT_FAIL;
  }

  solver.cvode_mem = CVodeCreate(CV_BDF, solver.sunctx);
  if (solver.cvode_mem == NULL) {
    uclchem_cvode_free_solver();
    return CV_MEM_FAIL;
  }

  retval = CVodeSetMaxOrd(solver.cvode_mem, 3);
  if (retval != CV_SUCCESS) {
    uclchem_cvode_free_solver();
    return retval;
  }

  solver.callbacks.rhs = rhs;
  solver.callbacks.user_data = NULL;
  solver.initialized = 0;

  return CV_SUCCESS;
}

int uclchem_cvode_integrate(double *y, double *t, double tout, double reltol, const double *abstol,
                            long int mxstep) {
  int retval;
  sunindextype i;
  double *y_data;
  double *abstol_data;

  if (solver.cvode_mem == NULL || solver.yvec == NULL || solver.abstol_vec == NULL ||
      solver.A == NULL || solver.LS == NULL) {
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

  if (!solver.initialized) {
    retval = CVodeInit(solver.cvode_mem, rhs_bridge, (realtype)(*t), solver.yvec);
    if (retval != CV_SUCCESS) {
      return retval;
    }

    retval = CVodeSetUserData(solver.cvode_mem, &solver.callbacks);
    if (retval != CV_SUCCESS) {
      return retval;
    }

    retval = CVodeSetLinearSolver(solver.cvode_mem, solver.LS, solver.A);
    if (retval != CV_SUCCESS) {
      return retval;
    }

    retval = CVodeSetJacFn(solver.cvode_mem, jac_bridge);
    if (retval != CV_SUCCESS) {
      return retval;
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

  return retval;
}

int uclchem_cvode_finalize(void) {
  uclchem_cvode_free_solver();
  return CV_SUCCESS;
}
