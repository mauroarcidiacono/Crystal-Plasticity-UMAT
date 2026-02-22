      SUBROUTINE LUBKSB (A, N, NP, INDX, B)

!-----  Linear equation solver based on LU decomposition
!-----  Lower Upper Backwards Substitution
!
! Parameters
! A: The matrix containing LU decomposition.
! N: Number of rows in the solution vector and columns of A (practically, the size of our system).
! NP: Total number of rows and columns in the matrix A.
! INDX: Array that stores the row permutations due to pivoting.
! B: The right-hand side vector of the equation system. On exit, it contains the solution.

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NP,NP), INDX(N), B(N)

      ! Forward substitution - calculating LD = B
      ! Stores the solution of D in B.
      II=0
      DO I=1,N
         LL=INDX(I)   ! The original row location before any pivoting
         SUM=B(LL)    ! Accumulator for the matrix-vector multiplication
         B(LL)=B(I)

         IF (II.NE.0) THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            END DO
         ELSE IF (SUM.NE.0.) THEN
            II=I
         END IF

         B(I)=SUM
      END DO

      ! Backward substitution - calculating UX = D
      ! D was stored in B.
      DO I=N,1,-1
         SUM=B(I)

         IF (I.LT.N) THEN
            DO J=I+1,N
               SUM=SUM-A(I,J)*B(J)
            END DO
         END IF

         B(I)=SUM/A(I,I)
      END DO

      RETURN
      END