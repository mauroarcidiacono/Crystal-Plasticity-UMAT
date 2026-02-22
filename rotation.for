      SUBROUTINE ROTATION (PROP, ROTATE)

!-----  This subroutine calculates the rotation matrix, i.e. the 
!       direction cosines of crystal [100], [010] and [001] 
!       directions in global system.

!-----  The rotation matrix is stored in the array ROTATE.

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PROP(16), ROTATE(3,3), TERM1(3,3), TERM2(3,3), INDX(3) 

!-----  Subroutines:
!
!       CROSS  -- cross product of two vectors
!
!       LUDCMP -- LU decomposition
!
!       LUBKSB -- linear equation solver based on LU decomposition 
!                 method (must call LUDCMP first)


!-----  PROP -- constants characterizing the crystal orientation 
!               (INPUT)
!
!            PROP(1) - PROP(3) -- direction of the first vector in 
!                                 local cubic crystal system
!            PROP(4) - PROP(6) -- direction of the first vector in 
!                                 global system
!
!            PROP(9) - PROP(11)-- direction of the second vector in 
!                                 local cubic crystal system
!            PROP(12)- PROP(14)-- direction of the second vector in 
!                                 global system              
!
!-----  ROTATE -- rotation matrix (OUTPUT):
!
!            ROTATE(i,1) -- direction cosines of direction [1 0 0] in 
!                           local cubic crystal system
!            ROTATE(i,2) -- direction cosines of direction [0 1 0] in 
!                           local cubic crystal system
!            ROTATE(i,3) -- direction cosines of direction [0 0 1] in 
!                           local cubic crystal system

!-----  local matrix: TERM1
      CALL CROSS (PROP(1), PROP(9), TERM1, ANGLE1)

!-----  LU decomposition of TERM1
      CALL LUDCMP (TERM1, 3, 3, INDX, DCMP, NOEL)

!-----  inverse matrix of TERM1: TERM2
      DO J=1,3
         DO I=1,3
            IF (I.EQ.J) THEN
               TERM2(I,J)=1.
            ELSE
               TERM2(I,J)=0.
            END IF
         END DO
      END DO

      DO J=1,3
         CALL LUBKSB (TERM1, 3, 3, INDX, TERM2(1,J))
      END DO

!-----  global matrix: TERM1
      CALL CROSS (PROP(4), PROP(12), TERM1, ANGLE2)

!-----  Check: the angle between first and second vector in local and 
!     global systems must be the same.  The relative difference must be
!     less than 0.1%. 
!

      IF (ABS(ANGLE1/ANGLE2-1.).GT.0.001) THEN 
         WRITE (6,*) 
     2      '***ERROR - angles between two vectors are not the same'
            CALL XIT
      END IF

!-----  rotation matrix: ROTATE
!       TERM2 is the inverse of TERM1, which represents the global
!       to local transformation. Therefore, TERM2 is a linear map
!       of local -> global.
!       TERM1 performs the mapping from the first pair of vectors
!       to the second pair of vectors. 
!       To use a non-cubic crystal, a chage of basis would be
!       required. However, in this implementation that was done
!       directly in the noncubiccrystal subroutine.
      DO J=1,3
         DO I=1,3
            ROTATE(I,J)=0.
            DO K=1,3
               ROTATE(I,J)=ROTATE(I,J)+TERM1(I,K)*TERM2(K,J)
            END DO
         END DO
      END DO

      RETURN
      END