program test_ludcmp
    implicit none
    integer, parameter :: n = 3
    double precision :: A(n,n), origA(n,n), det
    integer :: indx(n)
    integer :: i, j, k
    logical :: pass
    double precision :: L(n,n), U(n,n), LU(n,n), PA(n,n), temp(n)

    ! Test matrix
    origA = reshape([0d0, 2d0, 1d0, &
                     2d0, 2d0, 3d0, &
                     4d0, -3d0, 8d0], shape(origA))
    A = origA

    ! Call LUDCMP
    call ludcmp(A, n, n, indx, det, 0)  ! NOEL=0 for testing

    ! Extract L and U from packed A
    L = 0.0d0
    U = 0.0d0
    do i = 1, n
        L(i,i) = 1.0d0
        do j = 1, i-1
            L(i,j) = A(i,j)
        end do
        do j = i, n
            U(i,j) = A(i,j)
        end do
    end do

    ! Reconstruct PA by applying row swaps to origA
    PA = origA
    do j = n, 1, -1
        i = indx(j)
        if (i /= j) then
            temp = PA(j, :)
            PA(j, :) = PA(i, :)
            PA(i, :) = temp
        end if
    end do

    ! Compute L * U
    LU = matmul(L, U)

    ! Compare LU to P*A
    pass = .true.
    do i = 1, n
        do j = 1, n
            if (abs(LU(i,j) - PA(i,j)) > 1.0d-10) then
                pass = .false.
            end if
        end do
    end do

    if (pass) then
        print *, "LUDCMP test passed."
    else
        print *, "LUDCMP test FAILED."
        print *, "P * A:"
        print *, PA
        print *, "L * U:"
        print *, LU
        stop 1
    end if
end program
