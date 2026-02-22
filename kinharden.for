      SUBROUTINE KINHARDEN (CKIN1, DKIN1, CKIN2, DKIN2,CKIN3, DKIN3,
     1                      DGAM, XKINSLP1, XKINSLP2, XKINSLP3,
     2                      DKINSLP1, DKINSLP2, DKINSLP3, NSLIP)

!-----  This subroutine calculates the kinematic hardening variable for 
!       each slip system based on the Armstrong-Frederick kinematic 
!       hardening rule.

!-----  Variables:
!
!     CKIN    -- kinematic hardening constant controlling saturation   
!                value of kinematic hardening variable (INPUT)
!     DKIN    -- kinematic hardening constant controlling rate at which     
!                saturation of kinematic hardening variable
!                is achieved (INPUT)
!     GAMDOT  -- current slip system plastic strain rate (INPUT)
!     XKINSLP -- current slip system kinematic hardening variable (INPUT)
!	DKINSLP -- current slip system increment of kinematic hardening 
!		     variable (OUTPUT)
!     NSLIP   -- number of slip systems in each set (INPUT)
!	DTIME   -- time increment

!-----  Use single precision on cray
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DGAM(NSLIP), XKINSLP1(NSLIP), XKINSLP2(NSLIP),  
     1          XKINSLP3(NSLIP), DKINSLP1(NSLIP), DKINSLP2(NSLIP),
     2          DKINSLP3(NSLIP) 

      DO I=1,NSLIP
         DKINSLP1(I)=CKIN1*DGAM(I)-DKIN1*XKINSLP1(I)*ABS(DGAM(I))
         DKINSLP2(I)=CKIN2*DGAM(I)-DKIN2*XKINSLP2(I)*ABS(DGAM(I))
         DKINSLP3(I)=CKIN3*DGAM(I)-DKIN3*XKINSLP3(I)*ABS(DGAM(I))
      END DO

      RETURN
      END