      SUBROUTINE NONCUBICCRYSTAL(SLPDIR, SLPNOR, ROTATE, NSLIP, NSLPTL)

! This subroutine can be used to adjust the UMAT for a noncubic crystal. 
! The current code will be accessed when PROPS(9).NE.0, meaning that the
! introduced properties are not for a cubic crystal.
!
! In this particular case, the code works for HCP but can be used for other
! noncubic systems. For that, the change of basis matrix (CBM) must be 
! adjusted. Then, the slip planes and directions must be changed according
! to the number of desired slip systems. The values must be expressed in a 
! three axes coordinate system that are oriented according to the crystal.
! The CBM will be used to transform from the crystal axes to the orthogonal
! basis. Finally, the values of NSLIP and NSLPT must be changed for the 
! specific case.

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION NSLIP(3), ROTATE(3,3), SLPDIR(3,150), SLPNOR(3,150),
     2             T1(3), T2(3), CBM(3,3), HSLPNOR(3,150),   
     3             HSLPDIR(3,150)

! This code was used to model Ti-64, which has an alpha phase with an HCP
! crystal structure.
! Lattice parameters
! c = 0.4673
! a = 0.2935

! Create the change of basis matrix.

      CBM(1,1) = 1
      CBM(1,2) = -0.5
      CBM(1,3) = 0

      CBM(2,1) = 0
      CBM(2,2) = SQRT(3.0)/2
      CBM(2,3) = 0

      CBM(3,1) = 0
      CBM(3,2) = 0
      CBM(3,3) = 1.592163543

! Define the HCP slip planes (normal vectors).

      ! Basal
      HSLPNOR(1,1) = 0
      HSLPNOR(2,1) = 0
      HSLPNOR(3,1) = 1

      HSLPNOR(1,2) = 0
      HSLPNOR(2,2) = 0
      HSLPNOR(3,2) = 1

      HSLPNOR(1,3) = 0
      HSLPNOR(2,3) = 0
      HSLPNOR(3,3) = 1

      ! Prismatic
      HSLPNOR(1,4) = 2
      HSLPNOR(2,4) = 1
      HSLPNOR(3,4) = 0

      HSLPNOR(1,5) = -1
      HSLPNOR(2,5) = -2
      HSLPNOR(3,5) = 0

      HSLPNOR(1,6) = 1
      HSLPNOR(2,6) = -1
      HSLPNOR(3,6) = 0

      ! First order pyramidal <a>
      HSLPNOR(1,7) = 2
      HSLPNOR(2,7) = 1
      HSLPNOR(3,7) = 0.5917195362

      HSLPNOR(1,8) = -1
      HSLPNOR(2,8) = -2
      HSLPNOR(3,8) = 0.5917195362

      HSLPNOR(1,9) = 1
      HSLPNOR(2,9) = -1
      HSLPNOR(3,9) = 0.5917195362

      HSLPNOR(1,10) = -2
      HSLPNOR(2,10) = -1 
      HSLPNOR(3,10) = 0.5917195362

      HSLPNOR(1,11) = 1
      HSLPNOR(2,11) = 2
      HSLPNOR(3,11) = 0.5917195362

      HSLPNOR(1,12) = -1
      HSLPNOR(2,12) = 1
      HSLPNOR(3,12) = 0.5917195362

! Define the HCP slip directions.

      ! Basal
      HSLPDIR(1,1) = 0
      HSLPDIR(2,1) = -1
      HSLPDIR(3,1) = 0

      HSLPDIR(1,2) = 1
      HSLPDIR(2,2) = 0
      HSLPDIR(3,2) = 0

      HSLPDIR(1,3) = 1
      HSLPDIR(2,3) = 1
      HSLPDIR(3,3) = 0

      ! Prismatic
      HSLPDIR(1,4) = 0
      HSLPDIR(2,4) = 1
      HSLPDIR(3,4) = 0

      HSLPDIR(1,5) = 1
      HSLPDIR(2,5) = 0
      HSLPDIR(3,5) = 0

      HSLPDIR(1,6) = 1
      HSLPDIR(2,6) = 1
      HSLPDIR(3,6) = 0

      ! First order pyramidal <a>
      HSLPDIR(1,7) = 0
      HSLPDIR(2,7) = 1
      HSLPDIR(3,7) = 0

      HSLPDIR(1,8) = 1
      HSLPDIR(2,8) = 0
      HSLPDIR(3,8) = 0

      HSLPDIR(1,9) = 1
      HSLPDIR(2,9) = 1
      HSLPDIR(3,9) = 0

      HSLPDIR(1,10) = 0
      HSLPDIR(2,10) = -1
      HSLPDIR(3,10) = 0

      HSLPDIR(1,11) = -1
      HSLPDIR(2,11) = 0
      HSLPDIR(3,11) = 0

      HSLPDIR(1,12) = -1
      HSLPDIR(2,12) = -1
      HSLPDIR(3,12) = 0


! Multiply the vectors by the change of basis matrix to obtain their coordinates in
! a orthogonal and orthonormal coordinate system.  

      DO J=1,12
         DO I=1,3
            T1(I)=0.
            DO K=1,3
               T1(I)=T1(I)+CBM(I,K)*HSLPDIR(K,J)
            END DO
         END DO
         DO I=1,3
            HSLPDIR(I,J)=T1(I)
         END DO

         DO I=1,3
            T2(I)=0.
            DO K=1,3
               T2(I)=T2(I)+CBM(I,K)*HSLPNOR(K,J)
            END DO
         END DO
         DO I=1,3
            HSLPNOR(I,J)=T2(I)
         END DO
      END DO 

! Rotate the slip directions and planes.

      DO J=1,12
         DO I=1,3
            T1(I)=0.
            DO K=1,3
               T1(I)=T1(I)+ROTATE(I,K)*HSLPDIR(K,J)
            END DO
         END DO
         DO I=1,3
            SLPDIR(I,J)=T1(I)
         END DO

         DO I=1,3
            T2(I)=0.
            DO K=1,3
               T2(I)=T2(I)+ROTATE(I,K)*HSLPNOR(K,J)
            END DO
         END DO
         DO I=1,3
            SLPNOR(I,J)=T2(I)
         END DO
      END DO

! Transform the SLPDIR and SLPNOR vectors to unit vectors.  

      DUMMY1 = 0
      DUMMY2 = 0
      
      DO J=1,12
         DO I=1,3
            DUMMY1 = DUMMY1 + SLPDIR(I,J)*SLPDIR(I,J)
            DUMMY2 = DUMMY2 + SLPNOR(I,J)*SLPNOR(I,J)
         END DO

         DUMMY1 = SQRT(DUMMY1)
         DUMMY2 = SQRT(DUMMY2)
         
         DO I=1,3
            SLPDIR(I,J) = SLPDIR(I,J)/DUMMY1
            SLPNOR(I,J) = SLPNOR(I,J)/DUMMY2
         END DO
         
         DUMMY1 = 0
         DUMMY2 = 0

      END DO

! Define the number of slip systems per set.

      NSLIP(1) = 3
      NSLIP(2) = 3
      NSLIP(3) = 6

! Define the total number of slip systems.
      
      NSLPTL = 12

      RETURN
      END