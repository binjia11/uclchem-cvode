SHELL = /bin/sh
################################################################################
#
# User Options
################################################################################
FC =gfortran
#FC=ifort
CC = cc

f2pyFC=gnu95 #this must match FC (ifort=intelem, gfortran=gnu95)

#Unforgiving debugging flags
#FFLAGS =-g -fbacktrace -Wall -fcheck=all -fPIC
#Fast optimizing flags
FFLAGS = -O3 -fPIC -ffree-line-length-0
CFLAGS = -O3 -fPIC
JACOBIAN_FFLAGS = -O0 -fPIC -ffree-line-length-0

SUNDIALS_PREFIX ?= $(CONDA_PREFIX)
SUNDIALS_CFLAGS =
SUNDIALS_LDFLAGS =
ifneq ($(strip $(SUNDIALS_PREFIX)),)
SUNDIALS_CFLAGS += -I$(SUNDIALS_PREFIX)/include
SUNDIALS_LDFLAGS += -L$(SUNDIALS_PREFIX)/lib
endif
SUNDIALS_LIBS = -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsolklu -lsundials_sunmatrixsparse -lsundials_generic -lklu -lamd -lcolamd -lbtf


###############################################################################
EXT_SUFFIX = $(shell python3-config --extension-suffix)


PHYSICS_OBS= physics-core.o cloud.o hotcore.o sputtering.o cshock.o jshock.o collapse.o postprocess.o
CHEMISTRY_OBS = odes.o photoreactions.o surfacereactions.o rates.o jacobian_sparse.o chemistry.o
CONSTANTS_OBS = constants.o defaultparameters.o f2py_constants.o
CVODE_OBS = cvode_bridge.o

##simple makefile for uclchem
##user must point the FC variable to their preferred fortran compiler
##builds ode solver, physics module and chemistry module before linking together for main

##User can also compile a python module which contains the subroutines in wrap.f90 as functions
##to do this run "python -m numpy.f2py -c --help-fcompiler" to find your fortran compiler
##and edit the python makef

##physics module selected by changing physics variable to chosen fortran file.
main: ${CVODE_OBS} ${CONSTANTS_OBS} network.o ${PHYSICS_OBS}  ${CHEMISTRY_OBS} io.o wrap.f90
	${FC} ${FFLAGS} ${SUNDIALS_LDFLAGS} -o ../../uclchem ${CONSTANTS_OBS} io.o ${PHYSICS_OBS} ${CVODE_OBS} ${CHEMISTRY_OBS} network.o wrap.f90 main.f90 ${SUNDIALS_LIBS}

# The order here is the compilation order and it's vital that modules are compiled after their dependencies
# We need to pass f2py_constants.f90 into wrap.f90 in order for nspec (and other parameters) to be in scope in wrap.
python: ${CONSTANTS_OBS} network.o ${PHYSICS_OBS} ${CVODE_OBS} ${CHEMISTRY_OBS} io.o wrap.f90
	python3 -m numpy.f2py -c --fcompiler=${f2pyFC} ${SUNDIALS_CFLAGS} ${SUNDIALS_LDFLAGS} ${CONSTANTS_OBS} network.o io.o ${PHYSICS_OBS} ${CVODE_OBS} ${CHEMISTRY_OBS} f2py_constants.f90 wrap.f90 ${SUNDIALS_LIBS} -m uclchemwrap
# Optional useful flags for debugging: --debug-capi, --debug, --f90flags="-fcheck=all"

chemistry.o: photoreactions.f90 surfacereactions.f90 rates.f90 odes.f90 jacobian_sparse.f90 chemistry.f90
	${FC} ${FFLAGS} -c chemistry.f90

jacobian_sparse.o: physics-core.f90 jacobian_sparse.f90
	${FC} ${JACOBIAN_FFLAGS} -c jacobian_sparse.f90

sputtering.o: surfacereactions.o sputtering.f90
	${FC} ${FFLAGS} -c sputtering.f90

%.o: %.f90 
	${FC} ${FFLAGS} -c $<

%.o: %.c
	${CC} ${CFLAGS} ${SUNDIALS_CFLAGS} -c $<

# This can de disabled for fast development, but it will cause major issues when you change chemical networks.
clean: 
	rm -f *.o *.mod *.so ../uclchem/*.so
