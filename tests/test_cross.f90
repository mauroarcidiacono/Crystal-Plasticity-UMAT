program test_cross
    implicit none
    double precision :: A(3), B(3), C(3,3), ANGLE
    double precision :: dot_ab, dot_ac, dot_bc
    double precision :: norm1, norm2, norm3
    double precision, parameter :: TOL = 1.0d-10
    logical :: pass
    integer :: i

    ! Test input: orthogonal unit vectors
    A = (/ 1.0d0, 0.0d0, 0.0d0 /)
    B = (/ 0.0d0, 1.0d0, 0.0d0 /)

    ! Call cross subroutine
    call CROSS(A, B, C, ANGLE)

    ! Compute norms of the three vectors
    norm1 = sqrt(sum(C(:,1)**2))
    norm2 = sqrt(sum(C(:,2)**2))
    norm3 = sqrt(sum(C(:,3)**2))

    ! Check orthogonality
    dot_ab = sum(C(:,1)*C(:,2))
    dot_ac = sum(C(:,1)*C(:,3))
    dot_bc = sum(C(:,2)*C(:,3))

    ! Check results
    pass = .true.
    if (abs(norm1 - 1.0d0) > TOL) pass = .false.
    if (abs(norm2 - 1.0d0) > TOL) pass = .false.
    if (abs(norm3 - 1.0d0) > TOL) pass = .false.
    if (abs(dot_ab) > TOL) pass = .false.
    if (abs(dot_ac) > TOL) pass = .false.
    if (abs(dot_bc) > TOL) pass = .false.
    if (abs(ANGLE - 1.57079632679d0) > TOL) pass = .false.

    ! Output
    if (pass) then
        print *, "CROSS test passed."
    else
        print *, "CROSS test FAILED."
        print *, "C(:,1) =", C(:,1)
        print *, "C(:,2) =", C(:,2)
        print *, "C(:,3) =", C(:,3)
        print *, "ANGLE =", ANGLE
        stop 1
    end if
end program
