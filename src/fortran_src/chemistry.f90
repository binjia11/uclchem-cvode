! Chemistry module of UCL_CHEM.                                                               !
! Contains all the core machinery of the code, not really intended to be altered in standard  !
! use. Use a (custom) physics module to alter temp/density behaviour etc.                     !
!                                                                                             !
! chemistry module contains rates.f90, a series of subroutines to calculate all reaction rates!
! when updateChemistry is called from main, these rates are calculated, the ODEs are solved   !
! from currentTime to targetTime to get abundances at targetTime and then all abundances are  !
! written to the fullOutput file.                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE chemistry
USE constants
USE DEFAULTPARAMETERS
USE iso_c_binding, ONLY: c_int, c_long, c_double, c_ptr, c_funptr, c_funloc, c_f_pointer
!f2py INTEGER, parameter :: dp
USE physicscore, only: points, dstep, cloudsize, radfield, h2crprate, improvedH2CRPDissociation, &
& zeta, currentTime, targetTime, timeinyears, freefall, density, ion, densdot, gastemp, dusttemp, av
USE network
USE photoreactions
USE surfacereactions
USE jacobian, ONLY: jac_nnz, GETJACOBIAN_DENSE, GETJACOBIAN_SPARSE_VALUES, SPARSE_VALUES_TO_DENSE
use f2py_constants, only: nspec, nreac
USE postprocess_mod, only: lusecoldens,usepostprocess,tstep,lnh,lnh2,lnco,lnc
USE rates
USE odes
IMPLICIT NONE
    !f2py integer, intent(aux) :: points
    !These integers store the array index of important species and reactions, x is for ions    
    !loop counters    
    INTEGER :: i,j,l,writeCounter=0,loopCounter,failedIntegrationCounter
    INTEGER, PARAMETER :: maxLoops=10,maxConsecutiveFailures=10

    !Array to store reaction rates
    REAL(dp) :: rate(nreac)
    
    !ODE integrator state
    INTEGER :: ISTATE,NEQ
    REAL(dp), ALLOCATABLE :: abstol(:)
    !initial fractional elemental abudances and arrays to store abundances
    REAL(dp) :: h2col,cocol,ccol,h2colToCell,cocolToCell,ccolToCell
    REAL(dp), ALLOCATABLE :: abund(:,:)
    
    REAL(dp) :: MIN_ABUND = 1.0d-30 !Minimum abundance allowed

    INTEGER :: nion,ionlist(nspec)
    LOGICAL :: cvode_log_open = .FALSE.
    INTEGER, PARAMETER :: cvode_log_unit = 97
    CHARACTER(len=*), PARAMETER :: cvode_log_path = "cvode_solver.log"
    LOGICAL :: jac_debug_enabled = .FALSE.
    INTEGER :: jac_debug_max_calls = 0
    INTEGER :: jac_debug_call_count = 0
    PRIVATE :: update_state_dependent_rates, ensure_cvode_log_open, log_cvode_event, &
    & prepare_jacobian_state, configure_jacobian_debug, log_jacobian_comparison

    ! CVODE return codes from cvode/cvode.h
    INTEGER, PARAMETER :: CV_SUCCESS=0, CV_TSTOP_RETURN=1
    INTEGER, PARAMETER :: CV_TOO_MUCH_WORK=-1, CV_TOO_MUCH_ACC=-2
    INTEGER, PARAMETER :: CV_ERR_FAILURE=-3, CV_CONV_FAILURE=-4
    INTEGER, PARAMETER :: CV_LINIT_FAIL=-5, CV_LSETUP_FAIL=-6, CV_LSOLVE_FAIL=-7
    INTEGER, PARAMETER :: CV_RHSFUNC_FAIL=-8, CV_FIRST_RHSFUNC_ERR=-9
    INTEGER, PARAMETER :: CV_REPTD_RHSFUNC_ERR=-10, CV_UNREC_RHSFUNC_ERR=-11
    INTEGER, PARAMETER :: CV_NLS_INIT_FAIL=-13, CV_NLS_SETUP_FAIL=-14, CV_NLS_FAIL=-16
    INTEGER, PARAMETER :: CV_MEM_FAIL=-20, CV_MEM_NULL=-21, CV_ILL_INPUT=-22
    INTEGER, PARAMETER :: CV_NO_MALLOC=-23, CV_CONTEXT_ERR=-32

    INTERFACE
        INTEGER(c_int) FUNCTION uclchem_cvode_init(neq, rhs) &
                & BIND(C, name="uclchem_cvode_init")
            USE iso_c_binding, ONLY: c_int, c_funptr
            INTEGER(c_int), VALUE :: neq
            TYPE(c_funptr), VALUE :: rhs
        END FUNCTION uclchem_cvode_init

        INTEGER(c_int) FUNCTION uclchem_cvode_integrate(y, t, tout, reltol, abstol, mxstep) &
                & BIND(C, name="uclchem_cvode_integrate")
            USE iso_c_binding, ONLY: c_int, c_long, c_double, c_funptr
            REAL(c_double), INTENT(INOUT) :: y(*)
            REAL(c_double), INTENT(INOUT) :: t
            REAL(c_double), VALUE :: tout, reltol
            REAL(c_double), INTENT(IN) :: abstol(*)
            INTEGER(c_long), VALUE :: mxstep
        END FUNCTION uclchem_cvode_integrate

        INTEGER(c_int) FUNCTION uclchem_cvode_finalize() BIND(C, name="uclchem_cvode_finalize")
            USE iso_c_binding, ONLY: c_int
        END FUNCTION uclchem_cvode_finalize
    END INTERFACE
CONTAINS
    SUBROUTINE initializeChemistry(readAbunds)
        LOGICAL, INTENT(IN) :: readAbunds
        !f2py integer, intent(aux) :: points
        CHARACTER(len=32) :: linear_solver_env
        INTEGER :: env_length, env_status

        ! Sets variables at the start of every run.
        ! Since python module persists, it's not enough to set initial
        ! values in module definitions above. Reset here.
        NEQ=nspec+1
        IF (ALLOCATED(abund)) DEALLOCATE(abund,vdiff)
        ALLOCATE(abund(NEQ,points),vdiff(SIZE(iceList)))
        !Set abundances to initial elemental if not reading them in.
        IF (.NOT. readAbunds) THEN
            !ensure abund is initially zero
            abund= MIN_ABUND

            !Start by filling all metallicity scaling elements
            !neutral atoms  
            abund(no,:) = fo  
            abund(nn,:) = fn               
            abund(nmg,:) = fmg
            abund(np,:) = fp
            abund(nf,:) = ff
            !abund(nfe,:) = ffe
            abund(nna,:) = fna
            abund(nli,:) = fli
            abund(npah,:) = fpah
            !default to ions
            abund(nsx,:) = fs
            abund(nsix,:) = fsi                
            abund(nclx,:) = fcl 
            !Decide how much carbon is initiall ionized using parameters.f90
            SELECT CASE (ion)
                CASE(0)
                    abund(nc,:)=fc
                    abund(ncx,:)=1.d-10
                CASE(1)
                    abund(nc,:)=fc*0.5
                    abund(ncx,:)=fc*0.5
                CASE(2)
                    abund(nc,:)=1.d-10
                    abund(ncx,:)=fc
            END SELECT

            !isotopes
            abund(n18o,:) = f18o  
            abund(n15n,:) = f15n           
            abund(n13c,:) = f13c    

            abund(nelec,:)=abund(ncx,:)+abund(nsix,:)+abund(nsx,:)+abund(nclx,:)+abund(nmgx,:)

            abund=abund*metallicity

            !Total H nuclei is always 1 so put fh into H and whatever is left over in H2
            abund(nh,:) = fh
            abund(nh2,:) = 0.5*(1.0e0-fh) 
            abund(nd,:)=fd

            abund(nhe,:) = fhe  
        ENDIF
        abund(neq,:)=density  
        !Initial calculations of diffusion frequency for each species bound to grain
        !and other parameters required for diffusion reactions
        DO  i=lbound(iceList,1),ubound(iceList,1)
            j=iceList(i)
            vdiff(i)=VDIFF_PREFACTOR*bindingEnergy(i)/mass(j)
            vdiff(i)=dsqrt(vdiff(i))
        END DO

        ! get list of positive-charged species to conserve charge later
        nion = 0
        do i=1,nspec
           if (index(specname(i),'+') .ne. 0) then
              nion = nion + 1
              ionlist(nion) = i
           end if
        end do
        
        CALL ensure_cvode_log_open(reset=.TRUE.)
        CALL configure_jacobian_debug()
        CALL get_environment_variable("UCLCHEM_CVODE_LINEAR_SOLVER", linear_solver_env, &
        & length=env_length, status=env_status)
        IF (env_status == 0 .AND. env_length > 0 .AND. &
        & (TRIM(ADJUSTL(linear_solver_env(:env_length))) == "dense" .OR. &
        & TRIM(ADJUSTL(linear_solver_env(:env_length))) == "DENSE")) THEN
            CALL log_cvode_event("Initializing CVODE solver with dense linear solver")
        ELSE
            CALL log_cvode_event("Initializing CVODE solver with sparse KLU linear solver")
        END IF

        ISTATE = INT(uclchem_cvode_finalize(), KIND(ISTATE))
        ISTATE = INT(uclchem_cvode_init(INT(NEQ, c_int), c_funloc(cvode_rhs_c)), KIND(ISTATE))
        IF (ISTATE /= CV_SUCCESS) THEN
            WRITE(*,*) "Failed to initialize CVODE solver", ISTATE
            CALL log_cvode_event("Failed to initialize CVODE solver", ISTATE)
            STOP 1
        END IF
        CALL log_cvode_event("CVODE solver initialized", ISTATE)

        !ODE integrator settings
        ISTATE=1

        !set integration counts
        loopCounter=0
        failedIntegrationCounter=0

        IF (.NOT. ALLOCATED(abstol)) THEN
            ALLOCATE(abstol(NEQ))
        END IF
        !Set rates to zero to ensure they don't hold previous values or random ones if we don't set them in calculateReactionRates
        rate=0.0
        !We typically don't recalculate rates that only depend on temperature if the temp hasn't changed
        !use arbitrarily high value to make sure they are calculated at least once.
        lastTemp=99.0d99
    END SUBROUTINE initializeChemistry



    SUBROUTINE updateChemistry(successFlag)
    !Updates the abundances for the next time step, first updating chemical variables and reaction rates,
    !then by solving the ODE system to obtain new abundances.
    !Solving ODEs is complex so we have two checks to try to automatically overcome difficulties and end stalled models
    !Firstly, the integration subroutine is called up to maxLoops times whilst adjusting variables to help integration converge.
    !If it succeeds before maxLoops, we continue as normal, otherwise we'll call it a fail.
    !Secondly, we check for stalls caused by the the solver loop reducing the targetTime to overcome difficulties.
    !That reduction the possibility of the code "succeeding" by integrating tiny target times. We have a counter that resets each time
    !the code integrates to the planned targetTime rather than a reduced one. If the counter reaches maxConsecutiveFailures, we end the code.
        !f2py integer, intent(aux) :: points
        INTEGER, INTENT(OUT) :: successFlag
        real(dp) :: originalTargetTime !targetTime can be altered by integrator but we'd like to know if it was changed
        real(dp) :: surfaceCoverage


        !Integration can fail in a way that we can manage. Allow maxLoops tries before giving up.
        loopCounter=0
        successFlag=0
        originalTargetTime=targetTime
        DO WHILE((currentTime .lt. targetTime) .and. (loopCounter .lt. maxLoops)) 
            !allow option for dens to have been changed elsewhere.
            IF (.not. freefall) abund(nspec+1,dstep)=density(dstep)

            !First sum the total column density over all points further towards edge of cloud
            IF (dstep.gt.1) THEN
                h2ColToCell=(sum(abund(nh2,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
                coColToCell=(sum(abund(nco,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
                cColToCell=(sum(abund(nc,:dstep-1)*density(:dstep-1)))*(cloudSize/real(points))
            ELSE
                h2ColToCell=0.0
                coColToCell=0.0
                cColToCell=0.0
            ENDIF
            !then add half the column density of the current point to get average in this "cell"
            h2Col=h2ColToCell+0.5*abund(nh2,dstep)*density(dstep)*(cloudSize/real(points))
            coCol=coColToCell+0.5*abund(nco,dstep)*density(dstep)*(cloudSize/real(points))
            cCol=cColToCell+0.5*abund(nc,dstep)*density(dstep)*(cloudSize/real(points))

            ! Postprocessed tracers have column densities provided
            if (lusecoldens) then
               h2col = lnh2(tstep)
               cocol = lnco(tstep)
               ! ccol = lnc(dstep, tstep) ! TODO enable C column density support
               ccol = lnh(tstep) * abund(nc,dstep) ! No C column densities yet...
            end if

            !Reset surface and bulk values in case of integration error or sputtering
            abund(nBulk,dstep)=sum(abund(bulkList,dstep))
            abund(nSurface,dstep)=sum(abund(surfaceList,dstep))
            !recalculate coefficients for ice processes
            safeMantle=MAX(1d-30,abund(nSurface,dstep))
            safeBulk=MAX(1d-30,abund(nBulk,dstep))
            
            if (refractoryList(1) .gt. 0) safeBulk=safeBulk-SUM(abund(refractoryList,dstep))
            bulkLayersReciprocal=MIN(1.0,NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
            surfaceCoverage=bulkGainFromMantleBuildUp()

            
            CALL calculateReactionRates(abund,safeMantle, h2col, cocol, ccol, rate)

            !Integrate chemistry, and return fail if unrecoverable error was reached
            CALL integrateODESystem(successFlag)
            IF (successFlag .lt. 0) THEN
                write(*,*) "Integration failed, exiting"
                RETURN
            END IF



            !1.d-30 stops numbers getting too small for fortran.
            WHERE(abund<MIN_ABUND) abund=MIN_ABUND
            density(dstep)=abund(NEQ,dstep)
            loopCounter=loopCounter+1

            ! For postprocessing, force solver to try and reach original target time
            if (usepostprocess) targettime = originaltargettime
        END DO

        ! Postprocessing needs to reach next timestep whatever the cost
        if (.not. usepostprocess) then
        IF (loopCounter .eq. maxLoops) successFlag=INT_TOO_MANY_FAILS_ERROR

        !Since targetTime can be altered, eventually leading to "successful" integration we want to
        !check if integrator ever just reaches the planned target time. If it doesn't for many attempts,
        !we will call the run a failure. This stops the target being constantly reduced to tiny increments
        !so that the code all but stalls as the time is increased by seconds each integraiton.
        IF (ABS(originalTargetTime- targetTime) .lt. 0.001*originalTargetTime) THEN
            failedIntegrationCounter=0
        ELSE
            failedIntegrationCounter=failedIntegrationCounter+1
        END IF
        IF (failedIntegrationCounter .gt. maxConsecutiveFailures)&
             &successFlag=INT_TOO_MANY_FAILS_ERROR
        end if
    END SUBROUTINE updateChemistry

    SUBROUTINE integrateODESystem(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
        INTEGER(c_int) :: cvode_state
        INTEGER(c_long) :: mxstep_cvode
        successFlag=0

    !This subroutine calls CVODE until it can reach targetTime with acceptable errors (reltol/abstol)
        abstol=abstol_factor*abund(:,dstep) !absolute tolerances depend on value of abundance
        WHERE(abstol<abstol_min) abstol=abstol_min ! to a minimum degree
        mxstep_cvode = INT(MXSTEP, c_long)
        cvode_state = uclchem_cvode_integrate(abund(:,dstep), currentTime, targetTime, reltol, abstol, mxstep_cvode)
        ISTATE = INT(cvode_state, KIND(ISTATE))

        SELECT CASE(ISTATE)
            CASE(CV_TOO_MUCH_WORK)
                ! CVODE couldn't take enough steps before reaching targetTime
                write(*,*) "ISTATE -1: Reducing time step to ", (targetTime-currentTime)*0.1/SECONDS_PER_YEAR, "years"
                CALL log_cvode_event("CVODE too much work; reducing targetTime", ISTATE)
                targetTime=currentTime+(targetTime-currentTime)*0.1
            CASE(CV_TOO_MUCH_ACC)
                ! Tolerances are too small for machine precision
                write(*,*) "ISTATE -2: Tolerances too small"
                CALL log_cvode_event("CVODE tolerances too small; increasing abstol_factor", ISTATE)
                abstol_factor=abstol_factor*10.0
            CASE(CV_MEM_FAIL, CV_MEM_NULL, CV_ILL_INPUT, CV_NO_MALLOC, CV_CONTEXT_ERR)
                write(*,*) "CVODE found invalid inputs"
                write(*,*) "abstol:"
                write(*,*) abstol
                CALL log_cvode_event("CVODE invalid input or memory/context failure", ISTATE)
                successFlag=INT_UNRECOVERABLE_ERROR
                RETURN
            CASE(CV_ERR_FAILURE, CV_CONV_FAILURE, CV_LINIT_FAIL, CV_LSETUP_FAIL, CV_LSOLVE_FAIL, &
                    & CV_NLS_INIT_FAIL, CV_NLS_SETUP_FAIL, CV_NLS_FAIL)
                write(*,*) "ISTATE -4 - shortening step"
                CALL log_cvode_event("CVODE nonlinear/linear solve failure; shortening step", ISTATE)
                targetTime=currentTime+(targetTime-currentTime)*0.1
            CASE(CV_RHSFUNC_FAIL, CV_FIRST_RHSFUNC_ERR, CV_REPTD_RHSFUNC_ERR, CV_UNREC_RHSFUNC_ERR)
                write(*,*) "CVODE RHS failure at time", timeInYears,"years"
                CALL log_cvode_event("CVODE RHS failure", ISTATE)
                successFlag=INT_UNRECOVERABLE_ERROR
                RETURN
            CASE(CV_SUCCESS, CV_TSTOP_RETURN)
                CALL log_cvode_event("CVODE step completed", ISTATE)
                MXSTEP=10000
            CASE default
                IF (ISTATE .lt. 0) THEN
                    write(*,*) "CVODE returned error code", ISTATE
                    CALL log_cvode_event("CVODE returned unclassified negative code; shortening step", ISTATE)
                    targetTime=currentTime+(targetTime-currentTime)*0.1
                ELSE
                    CALL log_cvode_event("CVODE returned unclassified non-negative code", ISTATE)
                    MXSTEP=10000
                END IF
        END SELECT
    if (enforceChargeConservation) then
        ! REALLY ensure charge is always conserved (also after integrating)
        abund(nelec,dstep) = sum(abund(ionlist(1:nion),dstep))
    end if
    END SUBROUTINE integrateODESystem

    SUBROUTINE ensure_cvode_log_open(reset)
        LOGICAL, INTENT(IN), OPTIONAL :: reset
        LOGICAL :: do_reset
        INTEGER :: ios

        do_reset = .FALSE.
        IF (PRESENT(reset)) do_reset = reset

        IF (cvode_log_open .AND. do_reset) THEN
            CLOSE(cvode_log_unit)
            cvode_log_open = .FALSE.
        END IF

        IF (.NOT. cvode_log_open) THEN
            IF (do_reset) THEN
                OPEN(unit=cvode_log_unit, file=cvode_log_path, status="replace", action="write", iostat=ios)
            ELSE
                OPEN(unit=cvode_log_unit, file=cvode_log_path, status="unknown", action="write", &
                     & position="append", iostat=ios)
            END IF

            IF (ios == 0) THEN
                cvode_log_open = .TRUE.
            ELSE
                WRITE(*,*) "WARNING: unable to open CVODE log file", TRIM(cvode_log_path), "iostat=", ios
            END IF
        END IF
    END SUBROUTINE ensure_cvode_log_open

    SUBROUTINE log_cvode_event(message, code)
        CHARACTER(len=*), INTENT(IN) :: message
        INTEGER, INTENT(IN), OPTIONAL :: code

        CALL ensure_cvode_log_open()
        IF (.NOT. cvode_log_open) RETURN

        WRITE(cvode_log_unit,'(A)') REPEAT("-", 72)
        WRITE(cvode_log_unit,'(A)') TRIM(message)
        IF (PRESENT(code)) WRITE(cvode_log_unit,'(A,1X,I0)') "ISTATE:", code
        WRITE(cvode_log_unit,'(A,1X,I0)') "dstep:", dstep
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "currentTime:", currentTime
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "targetTime:", targetTime
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "deltaTime:", targetTime-currentTime
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "timeInYears:", timeInYears
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "reltol:", reltol
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "abstol_factor:", abstol_factor
        WRITE(cvode_log_unit,'(A,1X,I0)') "MXSTEP:", MXSTEP
        FLUSH(cvode_log_unit)
    END SUBROUTINE log_cvode_event

    SUBROUTINE configure_jacobian_debug()
        CHARACTER(len=32) :: env_value
        INTEGER :: env_length, env_status

        jac_debug_enabled = .FALSE.
        jac_debug_max_calls = 0
        jac_debug_call_count = 0

        CALL get_environment_variable("UCLCHEM_JAC_DEBUG", env_value, length=env_length, status=env_status)
        IF (env_status /= 0 .OR. env_length <= 0) RETURN

        SELECT CASE (TRIM(ADJUSTL(env_value(:env_length))))
            CASE ("1", "true", "TRUE", "yes", "YES", "on", "ON")
                jac_debug_enabled = .TRUE.
            CASE DEFAULT
                RETURN
        END SELECT

        jac_debug_max_calls = 8
        CALL get_environment_variable("UCLCHEM_JAC_DEBUG_MAX_CALLS", env_value, length=env_length, status=env_status)
        IF (env_status == 0 .AND. env_length > 0) THEN
            READ(env_value(:env_length), *, ERR=10) jac_debug_max_calls
10          CONTINUE
        END IF
    END SUBROUTINE configure_jacobian_debug

    SUBROUTINE prepare_jacobian_state(Y, D, totalSwapLocal)
        REAL(dp), INTENT(IN) :: Y(:), D
        REAL(dp), INTENT(OUT) :: totalSwapLocal

        CALL update_state_dependent_rates(Y, D)
        safeMantle = MAX(1d-30, Y(nSurface))
        safeBulk = MAX(1d-30, Y(nBulk))
        bulkLayersReciprocal = MIN(1.0_dp, NUM_SITES_PER_GRAIN / (GAS_DUST_DENSITY_RATIO * safeBulk))
        totalSwapLocal = GETTOTALSWAP(rate, Y, bulkLayersReciprocal)
    END SUBROUTINE prepare_jacobian_state

    SUBROUTINE log_jacobian_comparison(t, Y, dense_jac, sparse_values)
        REAL(dp), INTENT(IN) :: t, Y(:)
        REAL(dp), INTENT(IN) :: dense_jac(:,:), sparse_values(:)
        REAL(dp) :: sparse_dense(NEQ, NEQ)
        REAL(dp) :: abs_diff, rel_diff, denom, max_abs_diff, max_rel_diff
        REAL(dp) :: dense_value, sparse_value
        INTEGER :: row, col, worst_row, worst_col, mismatch_count

        IF (.NOT. jac_debug_enabled) RETURN
        IF (jac_debug_call_count >= jac_debug_max_calls) RETURN

        CALL SPARSE_VALUES_TO_DENSE(sparse_values, sparse_dense)

        max_abs_diff = 0.0_dp
        max_rel_diff = 0.0_dp
        mismatch_count = 0
        worst_row = 1
        worst_col = 1
        dense_value = dense_jac(1,1)
        sparse_value = sparse_dense(1,1)

        DO col = 1, NEQ
            DO row = 1, NEQ
                abs_diff = ABS(dense_jac(row, col) - sparse_dense(row, col))
                denom = MAX(ABS(dense_jac(row, col)), ABS(sparse_dense(row, col)), 1.0d-30)
                rel_diff = abs_diff / denom
                IF (abs_diff > 1.0d-12 .AND. rel_diff > 1.0d-10) mismatch_count = mismatch_count + 1
                IF (abs_diff > max_abs_diff) THEN
                    max_abs_diff = abs_diff
                    max_rel_diff = rel_diff
                    worst_row = row
                    worst_col = col
                    dense_value = dense_jac(row, col)
                    sparse_value = sparse_dense(row, col)
                END IF
            END DO
        END DO

        jac_debug_call_count = jac_debug_call_count + 1
        CALL ensure_cvode_log_open()
        IF (.NOT. cvode_log_open) RETURN

        WRITE(cvode_log_unit,'(A)') REPEAT("=", 72)
        WRITE(cvode_log_unit,'(A,1X,I0)') "Jacobian compare call:", jac_debug_call_count
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "callback_t:", t
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "state_density:", Y(NEQ)
        WRITE(cvode_log_unit,'(A,1X,I0)') "mismatch_count:", mismatch_count
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "max_abs_diff:", max_abs_diff
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "max_rel_diff:", max_rel_diff
        WRITE(cvode_log_unit,'(A,1X,I0,A,I0)') "worst_entry:", worst_row, ",", worst_col
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "dense_value:", dense_value
        WRITE(cvode_log_unit,'(A,1X,ES24.16)') "sparse_value:", sparse_value
        FLUSH(cvode_log_unit)
    END SUBROUTINE log_jacobian_comparison

    INTEGER(c_int) FUNCTION cvode_rhs_c(t, y_ptr, ydot_ptr, user_data) BIND(C)
        REAL(c_double), VALUE :: t
        TYPE(c_ptr), VALUE :: y_ptr, ydot_ptr, user_data
        REAL(c_double), POINTER :: y(:), ydot(:)
        REAL(dp) :: t_local

        CALL c_f_pointer(y_ptr, y, [NEQ])
        CALL c_f_pointer(ydot_ptr, ydot, [NEQ])
        t_local = REAL(t, dp)
        CALL F(NEQ, t_local, y, ydot)
        cvode_rhs_c = 0_c_int
    END FUNCTION cvode_rhs_c

    INTEGER(c_int) FUNCTION cvode_sparse_jac_c(t, y_ptr, data_ptr) BIND(C, name="uclchem_cvode_sparse_jacobian")
        REAL(c_double), VALUE :: t
        TYPE(c_ptr), VALUE :: y_ptr, data_ptr
        REAL(c_double), POINTER :: y(:), jac_data(:)
        REAL(dp) :: D, totalSwapLocal
        REAL(dp) :: dense_jac(NEQ, NEQ)

        CALL c_f_pointer(y_ptr, y, [nspec + 1])
        CALL c_f_pointer(data_ptr, jac_data, [jac_nnz])

        D = y(nspec + 1)
        CALL prepare_jacobian_state(y, D, totalSwapLocal)
        CALL GETJACOBIAN_SPARSE_VALUES(rate, y, safeMantle, D, bulkLayersReciprocal, totalSwapLocal, jac_data)
        IF (jac_debug_enabled .AND. jac_debug_call_count < jac_debug_max_calls) THEN
            CALL GETJACOBIAN_DENSE(rate, y, safeMantle, D, bulkLayersReciprocal, totalSwapLocal, dense_jac)
            CALL log_jacobian_comparison(REAL(t, dp), y, dense_jac, jac_data)
        END IF

        cvode_sparse_jac_c = 0_c_int
    END FUNCTION cvode_sparse_jac_c

    INTEGER(c_int) FUNCTION cvode_dense_jac_c(t, y_ptr, data_ptr) BIND(C, name="uclchem_cvode_dense_jacobian")
        REAL(c_double), VALUE :: t
        TYPE(c_ptr), VALUE :: y_ptr, data_ptr
        REAL(c_double), POINTER :: y(:), jac_dense(:,:)
        REAL(dp) :: D, totalSwapLocal
        REAL(dp) :: sparse_values(jac_nnz)

        CALL c_f_pointer(y_ptr, y, [nspec + 1])
        CALL c_f_pointer(data_ptr, jac_dense, [NEQ, NEQ])

        D = y(nspec + 1)
        CALL prepare_jacobian_state(y, D, totalSwapLocal)
        CALL GETJACOBIAN_DENSE(rate, y, safeMantle, D, bulkLayersReciprocal, totalSwapLocal, jac_dense)
        IF (jac_debug_enabled .AND. jac_debug_call_count < jac_debug_max_calls) THEN
            CALL GETJACOBIAN_SPARSE_VALUES(rate, y, safeMantle, D, bulkLayersReciprocal, totalSwapLocal, sparse_values)
            CALL log_jacobian_comparison(REAL(t, dp), y, jac_dense, sparse_values)
        END IF

        cvode_dense_jac_c = 0_c_int
    END FUNCTION cvode_dense_jac_c

    SUBROUTINE update_state_dependent_rates(Y, D)
        REAL(dp), INTENT(IN) :: Y(:), D

        if (.not. lusecoldens) then
            cocol=coColToCell+0.5*Y(nco)*D*(cloudSize/real(points))
            h2col=h2ColToCell+0.5*Y(nh2)*D*(cloudSize/real(points))
            rate(nR_H2_hv)=H2PhotoDissRate(h2Col,radField,av(dstep),turbVel)
            rate(nR_CO_hv)=COPhotoDissRate(h2Col,coCol,radField,av(dstep))
        end if
    END SUBROUTINE update_state_dependent_rates

    SUBROUTINE F (NEQUATIONS, T, Y, YDOT)
        USE ODES
        INTEGER, PARAMETER :: WP = KIND(1.0D0)
        INTEGER NEQUATIONS
        REAL(WP) T
        REAL(WP), DIMENSION(NEQUATIONS) :: Y, YDOT
        INTENT(IN)  :: NEQUATIONS, T, Y
        INTENT(OUT) :: YDOT
        REAL(dp) :: D,loss,prod
        REAL(dp) :: surfaceCoverage
        REAL(dp) :: phi,cgr(6),grec,denom
        integer :: ii
        !Set D to the gas density for use in the ODEs
        D=y(NEQ)
        ydot=0.0

        ! Column densities are fixed for postprocessing data, so don't do this bit
        !changing abundances of H2 and CO can causes oscillation since their rates depend on their abundances
        !recalculating rates as abundances are updated prevents that.
        !thus these are the only rates calculated each time the ODE system is called.
        CALL update_state_dependent_rates(Y, D)

        !recalculate coefficients for ice processes
        safeMantle=MAX(1d-30,Y(nSurface))
        safeBulk=MAX(1d-30,Y(nBulk))
        bulkLayersReciprocal=MIN(1.0,NUM_SITES_PER_GRAIN/(GAS_DUST_DENSITY_RATIO*safeBulk))
        surfaceCoverage=bulkGainFromMantleBuildUp()

        !The ODEs created by MakeRates go here, they are essentially sums of terms that look like k(1,2)*y(1)*y(2)*dens. Each species ODE is made up
        !of the reactions between it and every other species it reacts with.
        ! INCLUDE 'odes.f90'
        CALL GETYDOT(RATE, Y, bulkLayersReciprocal, surfaceCoverage, safeMantle,safeBulk, D, YDOT)
        ! get density change from physics module to send to DLSODE
        
        if (enforceChargeConservation) then 
            ydot(nelec) = sum(ydot(ionlist(1:nion)))
            ! ! replace electron ydot with sum of positive ion ydots to conserve charge
            ! ydot(nelec) = 0.
            ! prod = 0.
            ! loss = 0.
            ! ! Enforce the conservation of charge by summing the ydot of all positive ions
            ! do ii=1,nion
            ! if (ydot(ionlist(ii)) .ge. 0.) then
            !     prod = prod + ydot(ionlist(ii))
            ! else
            !     loss = loss + ydot(ionlist(ii))
            ! end if
            ! end do
            ! ydot(nelec) = prod + loss
        end if 

        ydot(NEQUATIONS)=densdot(y(NEQUATIONS))
    
    END SUBROUTINE F

END MODULE chemistry
