program emulator
    implicit none

    ! This standalone emulator:
    ! - Applies a prescribed total strain
    ! - Splits it into increments
    ! - Calls UMAT sequentially
    ! - Outputs stress–strain results to emulator_results.txt

    ! ============================================================
    ! IMPORTANT: the emulator is configured for a plane stress 
    ! calculation. Modify the total_stran and DSTRAN arrays for 
    ! a different modelling assumption. Please check the USER 
    ! INPUT section below to change the applied strain.
    ! ============================================================ 

    ! Declaration of variables as per UMAT's argument list
    DOUBLE PRECISION :: start, finish, hyd_strain, hyd_stress, equiv_stress, equiv_strain
    DOUBLE PRECISION :: dev_stress(6), dev_strain(6)
    INTEGER :: NDI, NSHR, NTENS, NSTATV, NPROPS, NOEL, NPT, LAYER, i, j
    INTEGER :: KSPT, JSTEP(4), KINC, iostat, count, n_div
    DOUBLE PRECISION :: SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, DTEMP, PNEWDT, CELENT, value
    DOUBLE PRECISION :: max_value, temp_value, max_strain_increment
    CHARACTER(LEN=80) :: CMNAME

    ! The dimensions of these arrays should match with what your UMAT subroutine expects
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: STRESS, STATEV, DDSDDT, DRPLDE, STRAN, DSTRAN
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DDSDDE, DROT
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TIME, PREDEF, DPRED, PROPS, COORDS, DFGRD0, DFGRD1
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: total_stran, temp_dstran
    
    ! Set values for integer and real variables
    CMNAME = 'UMAT'
    SSE = 0.0
    SPD = 0.0
    SCD = 0.0
    RPL = 0.0
    DRPLDT = 0.0
    TEMP = 0.0
    DTEMP = 0.0
    PNEWDT = 0.0
    CELENT = 0.0
    NOEL = 0
    NPT = 0
    LAYER = 0
    KSPT = 0
    KINC = 0

    ! ============================================================
    ! USER INPUT SECTION
    ! ============================================================

    ! Number of state variables
    NSTATV = 651

    ! Number of direct stress/strain tensor components
    NDI = 3   

    ! Number of shear stress/strain tensor components
    NSHR = 1

    ! Number of total components
    NTENS = NDI + NSHR

    ! Integrated size of all the vectors passed to the UMAT
    ! Example: in this case, 8 * 23 = 184
    NPROPS = 184

    ! Allocate and initialize the arrays
    ALLOCATE(STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
        STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), &
        DROT(3,3), DFGRD0(NTENS), DFGRD1(NTENS), temp_dstran(NTENS), total_stran(NTENS))

    ! Total applied strain (Voigt notation)
    ! This is the final strain that the material point will reach
    total_stran(1) = 0.005D0   ! epsilon_11
    total_stran(2) = 0.0D0     ! epsilon_22
    total_stran(3) = 0.0D0     ! epsilon_33
    total_stran(4) = 0.0D0     ! gamma_12

    ! Time increment per step
    DTIME = 0.01

    ! Maximum strain increment per step that will be used to 
    ! reach the applied strain values
    max_strain_increment = 4.0e-4
    ! ============================================================

    ! Get the start time
    CALL CPU_TIME(start)
    
    STRESS = 0.0
    STATEV = 0.0
    DDSDDE = 0.0
    DDSDDT = 0.0
    DRPLDE = 0.0
    STRAN = 0.0

    ! Find the maximum applied strain
    max_value = abs(total_stran(1))
    do i = 2, NTENS
        if (abs(total_stran(i)) > max_value) then
            max_value = abs(total_stran(i))
        end if
    end do

    ! Calculate the number of steps required to reach the applied strain
    if (max_value == 0.0D0) then
        n_div = 1
    else
        n_div = CEILING(max_value / max_strain_increment)
    end if

    DSTRAN(1) = total_stran(1) / n_div
    DSTRAN(2) = total_stran(2) / n_div
    DSTRAN(3) = total_stran(3) / n_div
    DSTRAN(4) = total_stran(4) / n_div
    
    TIME = 0.0
    PREDEF = 0.0
    DPRED = 0.0
    COORDS = 0.0
    DROT = 0.0
    DFGRD0 = 0.0
    DFGRD1 = 0.0

    ! Read values for the material properties
    count = 0
    open(unit=10, file='cp_params.csv', status='old', action='read')
    do
        read(10, *, iostat=iostat) value
        if (iostat /= 0) exit

        count = count + 1
        ! Process the value as needed
        PROPS(count) = value
    end do
    close(unit=10)
    
    ! Open a file to write the results
    open(unit=20, file='emulator_results.txt', status='REPLACE', action='WRITE')

    equiv_stress = 0.0
    equiv_strain = 0.0

    do i = 1, n_div
        ! Call the subroutine
        CALL UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
            RPL, DDSDDT, DRPLDE, DRPLDT, & 
            STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
            NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
            CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, JSTEP, KINC)

        STRAN = STRAN + DSTRAN 

        hyd_stress = sum(STRESS(1:NDI))/3.d0
        hyd_strain = sum(STRAN(1:NDI))/3.d0

        ! Calculate the deviatoric tensor
        dev_stress(1:NDI) = STRESS(1:NDI) - hyd_stress
        dev_stress(NDI+1:NDI+NSHR) = STRESS(NDI+1:NDI+NSHR)

        dev_strain(1:NDI) = STRAN(1:NDI) - hyd_strain
        dev_strain(NDI+1:NDI+NSHR) = STRAN(NDI+1:NDI+NSHR)

        ! Compute the equivalent Von Mises stress
        ! stressVM = sqrt(3/2 * dev_stress : dev_stress)
        ! strainVM = sqrt(2/3 * dev_strain : dev_strain)
        equiv_stress = sqrt(3.d0/2.d0 * (dev_stress(1)**2.d0 + &
                               dev_stress(2)**2.d0 + dev_stress(3)**2.d0 + &
                               2.d0*dev_stress(4)**2.d0 + 2.d0*dev_stress(5)**2.d0 + &    
                               2.d0*dev_stress(6)**2.d0))

        ! Compute the equivalent Von Mises strain
        equiv_strain = sqrt(2.d0/3.d0 * (dev_strain(1)**2.d0 + &
                                dev_strain(2)**2.d0 + dev_strain(3)**2.d0 + &
                                2.d0*dev_strain(4)**2.d0 + 2.d0*dev_strain(5)**2.d0 + &
                                2.d0*dev_strain(6)**2.d0))

        ! write(20, *) equiv_stress, equiv_strain
        write(20, *) STRAN(1), STRESS(1)
    end do

    close(unit=20)

    print *, 'UMAT executed successfully.'
    do i = 1, NTENS
        print *, 'STRESS(', i, ') = ', STRESS(i)
    end do
    do i = 1, NTENS
        print *, 'STRAIN(', i, ') = ', STRAN(i)
    end do
    
    ! Get the finish time
    CALL CPU_TIME(finish)

    ! Print the elapsed time
    PRINT *, 'Elapsed time: ', finish - start, ' seconds.'

    print *, 'Press any key to exit.'
    read(*,*)

    ! Clean up
    DEALLOCATE(STRESS, STATEV, DDSDDE, DDSDDT, DRPLDE, &
        STRAN, DSTRAN, TIME, PREDEF, DPRED, PROPS, COORDS, &
        DROT, DFGRD0, DFGRD1, temp_dstran, total_stran)
end program emulator
