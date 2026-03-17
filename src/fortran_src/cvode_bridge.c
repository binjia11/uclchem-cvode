#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_context.h>

typedef int (*uclchem_rhs_fn)(double t, const double *y, double *ydot, void *user_data);

struct uclchem_rhs_data {
  uclchem_rhs_fn rhs;
  void *user_data;
};

static int rhs_bridge(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  struct uclchem_rhs_data *data = (struct uclchem_rhs_data *)user_data;
  double *y_ptr = N_VGetArrayPointer_Serial(y);
  double *ydot_ptr = N_VGetArrayPointer_Serial(ydot);

  if (data == NULL || data->rhs == NULL || y_ptr == NULL || ydot_ptr == NULL) {
    return CV_RHSFUNC_FAIL;
  }

  return data->rhs((double)t, y_ptr, ydot_ptr, data->user_data);
}

int uclchem_cvode_integrate(int neq, double *y, double *t, double tout, double reltol,
                            const double *abstol, long int mxstep, uclchem_rhs_fn rhs) {
  int retval = CV_SUCCESS;
  void *cvode_mem = NULL;
  SUNContext sunctx = NULL;
  N_Vector yvec = NULL;
  N_Vector abstol_vec = NULL;
  SUNMatrix A = NULL;
  SUNLinearSolver LS = NULL;
  struct uclchem_rhs_data rhs_data;
  int i;

  if (neq <= 0 || y == NULL || t == NULL || abstol == NULL || rhs == NULL) {
    return CV_ILL_INPUT;
  }

  retval = SUNContext_Create(NULL, &sunctx);
  if (retval != 0) {
    return CV_CONTEXT_ERR;
  }

  yvec = N_VMake_Serial((sunindextype)neq, y, sunctx);
  abstol_vec = N_VNew_Serial((sunindextype)neq, sunctx);
  if (yvec == NULL || abstol_vec == NULL) {
    retval = CV_MEM_FAIL;
    goto cleanup;
  }
  for (i = 0; i < neq; ++i) {
    NV_Ith_S(abstol_vec, i) = abstol[i];
  }

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (cvode_mem == NULL) {
    retval = CV_MEM_FAIL;
    goto cleanup;
  }

  rhs_data.rhs = rhs;
  rhs_data.user_data = NULL;
  retval = CVodeInit(cvode_mem, rhs_bridge, (realtype)(*t), yvec);
  if (retval != CV_SUCCESS) {
    goto cleanup;
  }

  retval = CVodeSetUserData(cvode_mem, &rhs_data);
  if (retval != CV_SUCCESS) {
    goto cleanup;
  }

  retval = CVodeSVtolerances(cvode_mem, (realtype)reltol, abstol_vec);
  if (retval != CV_SUCCESS) {
    goto cleanup;
  }

  if (mxstep > 0) {
    retval = CVodeSetMaxNumSteps(cvode_mem, mxstep);
    if (retval != CV_SUCCESS) {
      goto cleanup;
    }
  }

  A = SUNDenseMatrix((sunindextype)neq, (sunindextype)neq, sunctx);
  if (A == NULL) {
    retval = CV_MEM_FAIL;
    goto cleanup;
  }
  LS = SUNLinSol_Dense(yvec, A, sunctx);
  if (LS == NULL) {
    retval = CV_MEM_FAIL;
    goto cleanup;
  }

  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if (retval != CV_SUCCESS) {
    goto cleanup;
  }

  retval = CVode(cvode_mem, (realtype)tout, yvec, (realtype *)t, CV_NORMAL);

cleanup:
  if (cvode_mem != NULL) {
    CVodeFree(&cvode_mem);
  }
  if (LS != NULL) {
    SUNLinSolFree(LS);
  }
  if (A != NULL) {
    SUNMatDestroy(A);
  }
  if (abstol_vec != NULL) {
    N_VDestroy(abstol_vec);
  }
  if (yvec != NULL) {
    N_VDestroy(yvec);
  }
  if (sunctx != NULL) {
    SUNContext_Free(&sunctx);
  }

  return retval;
}
