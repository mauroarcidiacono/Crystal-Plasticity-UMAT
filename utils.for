           SUBROUTINE CROSS (A, B, C, ANGLE)

!-----  (1) normalize vectors A and B to unit vectors
!       (2) store A, B and A*B (cross product) in C

! C is a matrix that stores the normalized vectors A and B in the first two
! columns. The normalized vector C is the cross product of A and B and it
! is stored in the last column of C.
! ANGLE is the angle between the vectors A and B.

           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION A(3), B(3), C(3,3)

           SUM1=SQRT(A(1)**2+A(2)**2+A(3)**2)
           SUM2=SQRT(B(1)**2+B(2)**2+B(3)**2)

           IF (SUM1.EQ.0.) THEN
              WRITE (6,*) '***ERROR - first vector is zero'
              CALL XIT
           ELSE
              DO I=1,3
                 C(I,1)=A(I)/SUM1
              END DO
           END IF

           IF (SUM2.EQ.0.) THEN
              WRITE (6,*) '***ERROR - second vector is zero'
              CALL XIT
           ELSE
              DO I=1,3
                 C(I,2)=B(I)/SUM2
              END DO
           END IF

           ANGLE=0.
           DO I=1,3
              ANGLE=ANGLE+C(I,1)*C(I,2)
           END DO
           ANGLE=ACOS(ANGLE)

           C(1,3)=C(2,1)*C(3,2)-C(3,1)*C(2,2)
           C(2,3)=C(3,1)*C(1,2)-C(1,1)*C(3,2)
           C(3,3)=C(1,1)*C(2,2)-C(2,1)*C(1,2)
           SUM3=SQRT(C(1,3)**2+C(2,3)**2+C(3,3)**2)
           IF (SUM3.LT.1.E-8) THEN
              WRITE (6,*) 
     2           '***ERROR - first and second vectors are parallel'
               CALL XIT
           END IF

           RETURN
           END
