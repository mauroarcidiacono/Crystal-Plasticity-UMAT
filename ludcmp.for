      SUBROUTINE LUDCMP (A, N, NP, INDX, D, NOEL)

!-----  LU decomposition
!
! Parameters
! A: matrix to perform the LU decomposition. It will be used to store the matrices L and U.
! N: number of rows. It is used for pivoting.
! NP: number of rows and columns.
! INDX: stores the row index changes. INDX(J) will contain the row number of the pivot chosen for the column J.
! D: determinant sign of the matrix. Evaluates if the sign should be changed.

! Normalized pivoting is used in this subroutine.
! By scaling each row by the inverse of its largest absolute element, 
! we are essentially normalizing the rows such that we are choosing pivots based on 
! relative (scaled) sizes rather than absolute values. This helps to ensure that we are 
! pivoting based on the most numerically "significant" element in the column.


          IMPLICIT REAL*8 (A-H,O-Z)
          PARAMETER (NMAX=200, TINY=1.0E-20)
          INTEGER NOEL
          DIMENSION A(NP,NP), INDX(N), VV(NMAX)
    
          D=1.

          ! Checks for the maximum abs value each row.
          ! If all values = 0, then the matrix is singular and it has no inverse, then the program aborts (CALL XIT).
          ! Stores in VV the inverse of the max abs value of each row.
          DO I=1,N
             AAMAX=0.
    
             DO J=1,N
                IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
             END DO
             IF (AAMAX.EQ.0.) THEN
                WRITE(6,*) A
             END IF
             IF (AAMAX.EQ.0.) THEN 
                WRITE(6,*) 'Singular matrix.'
                WRITE(6,*) 'NOEL = ', NOEL
                CALL XIT
             END IF
             VV(I)=1./AAMAX
          END DO
          

          ! Main loop over column J for LU decomposition.
          DO J=1,N

             ! Loop over rows I, but only up to the current column J.
             ! Goal: Compute the upper triangle (U) of the LU decomposition.
             ! If there are not values on top of the diagonal, skip.
             DO I=1,J-1
                SUM=A(I,J)

                ! Loop over previous columns to adjust current element.
                DO K=1,I-1
                   SUM=SUM-A(I,K)*A(K,J)
                END DO
                
                ! Save U(I,J) in A(I,J)
                A(I,J)=SUM
             END DO
             AAMAX=0.


             ! Loop from row J down to the last row.
             ! Goal: Find the pivot element in column J.
             DO I=J,N
                SUM=A(I,J)

                ! Loop over previous columns to adjust current element.
                DO K=1,J-1
                   SUM=SUM-A(I,K)*A(K,J)
                END DO
    
                A(I,J)=SUM
                DUM=VV(I)*ABS(SUM)

                ! Check if the current element should be the pivot.
                ! IMAX is the best row to serve as pivot for the current column J.
                ! Division by the diagonal will be performed and in computational
                ! arithmetic using small numbers could lead to instability. 
                ! Therefore, partial pivoting is done to divide by the biggest number
                ! from the diagonal row index to bigger row indices (L matrix).
                IF (DUM.GE.AAMAX) THEN
                   IMAX=I
                   AAMAX=DUM
                END IF
             END DO
             
             ! If the pivot row isn't the current row, swap the rows.
             IF (J.NE.IMAX) THEN
                DO K=1,N
                   DUM=A(IMAX,K)
                   A(IMAX,K)=A(J,K)
                   A(J,K)=DUM
                END DO

                ! Change the sign of the determinant.  
                D=-D
                VV(IMAX)=VV(J)
             END IF
    
             INDX(J)=IMAX

             ! If pivot element is zero, set it to a very small value.
             IF (A(J,J).EQ.0.) A(J,J)=TINY

             ! For non-diagonal elements of the L matrix, scale the values in the current column.
             IF (J.NE.N) THEN
                DUM=1./A(J,J)
                DO I=J+1,N
                   A(I,J)=A(I,J)*DUM
                END DO
             END IF
    
          END DO
    
          RETURN
          END