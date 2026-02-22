program test_lubksb
    implicit none
    integer, parameter :: n = 3
    double precision :: A(n,n), B(n), X(n), det
    integer :: indx(n)
    logical :: pass
    integer :: i

    ! Test matrix A and known solution X = [1, 2, 3]
    A = reshape([2d0, 1d0, 1d0, &
                 4d0, -6d0, 0d0, &
                 -2d0, 7d0, 2d0], [n,n])

    B = matmul(A, [1d0, 2d0, 3d0])  ! Right-hand side: B = A * X_true

    ! Factor A
    call ludcmp(A, n, n, indx, det, 0)

    ! Solve A * X = B
    X = B
    call lubksb(A, n, n, indx, X)

    ! Check if X ≈ [1, 2, 3]
    pass = all(abs(X - [1d0, 2d0, 3d0]) < 1d-10)

    if (pass) then
        print *, "LUBKSB test passed."
    else
        print *, "LUBKSB test FAILED. Solution:"
        do i = 1, n
            print *, "X(", i, ") = ", X(i)
        end do
        stop 1
    end if
end program test_lubksb
