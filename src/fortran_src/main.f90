! 2022 UCLCHEM v3.0
! The canonical main file where the core code is written is actually found in wrap.f90
! main.f90 just provides a simple fortran interface to the core code so that a binary can
! be built and used directly from the command line.
PROGRAM uclchem

USE uclchemwrap, only: cloud,hot_core,cshock,postprocess
USE io, only: inputId
USE constants, only: dp
USE f2py_constants, only: nspec
IMPLICIT NONE
    CHARACTER (LEN=100):: modelType
    CHARACTER (LEN=100):: paramFile 
    CHARACTER (LEN=32):: model_arg1, model_arg2, model_arg3
    CHARACTER(:), ALLOCATABLE :: paramDict
    REAL(dp) :: abundances(nspec),dissipationResult
    CHARACTER(LEN=32) :: specname_out(nspec)
    INTEGER :: success,fileLength,model_index
    INTEGER, PARAMETER :: cli_timepoints=500, cli_gridpoints=1
    REAL(8) :: max_temp, vshock, timestep_factor, minimum_temp
    !Any subset of the parameters can be passed in a file on program start
    !see example.inp
    CALL GET_COMMAND_ARGUMENT(1, modelType)  
    CALL GET_COMMAND_ARGUMENT(2, paramFile)

    OPEN(unit=inputId, file=paramFile, action="read", &
    form="unformatted", access="stream")
    INQUIRE(unit=inputId, size=fileLength)

    ALLOCATE(character(fileLength) :: paramDict)
    READ(inputId) paramDict
    CLOSE(inputId)
    SELECT CASE(modelType)
    CASE("CLOUD")
        CALL cloud(dictionary=paramDict,outSpeciesIn="",returnArray=.false.,returnRates=.false.,&
        &givestartabund=.false.,timePoints=cli_timepoints,gridPoints=cli_gridpoints,&
        &abundance_out=abundances,specname_out=specname_out,successFlag=success)
    CASE("HOTCORE")
        !call hot_core (temp_indx,maxTemp,....)
        ! Read the strings from the CLI, then convert them to int/float
        CALL GET_COMMAND_ARGUMENT(3, model_arg1)  
        CALL GET_COMMAND_ARGUMENT(4, model_arg2)
        read(model_arg1, *) model_index
        read(model_arg2, *) max_temp
        CALL hot_core(temp_indx=model_index,max_temp=max_temp,dictionary=paramDict,outSpeciesIn="",&
        &returnArray=.false.,returnRates=.false.,givestartabund=.false.,&
        &timePoints=cli_timepoints,gridPoints=cli_gridpoints,&
        &abundance_out=abundances,specname_out=specname_out,successFlag=success)
    CASE("CSHOCK")
        !call cshock(vs,timestep_factor,minimum_temp,....)
        !sensible defaults in order: 20.0, 0.01d0, 10.0
        ! Read the strings from the CLI, then convert them to int/float
        CALL GET_COMMAND_ARGUMENT(3, model_arg1)  
        CALL GET_COMMAND_ARGUMENT(4, model_arg2)
        CALL GET_COMMAND_ARGUMENT(5, model_arg3)
        read(model_arg1, *) vshock
        read(model_arg2, *) timestep_factor
        read(model_arg3, *) minimum_temp
       CALL cshock(shock_vel=vshock,timestep_factor=timestep_factor,minimum_temperature=minimum_temp,&
       &dictionary=paramDict,outSpeciesIn="",returnArray=.false.,returnRates=.false.,&
       &givestartabund=.false.,timePoints=cli_timepoints,gridPoints=cli_gridpoints,&
       &abundance_out=abundances,dissipation_time=dissipationResult,specname_out=specname_out,&
       &successFlag=success)
    CASE("POSTPROCESS")
       WRITE(*,*) 'POSTPROCESS from Fortran CLI is not supported in this interface.'
       WRITE(*,*) 'Use Python API: uclchem.model.postprocess(...) with trajectory arrays.'
       STOP 2
    CASE default
        write(*,*) 'Model type not recognised'
        WRITE(*,*) 'Supported models are: CLOUD, CSHOCK, HOTCORE and POSTPROCESS'
        STOP
    END SELECT
END PROGRAM uclchem
