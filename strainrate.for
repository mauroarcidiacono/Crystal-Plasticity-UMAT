      SUBROUTINE STRAINRATE (GAMMA, TAUSLP, GSLIP, NSLIP, FSLIP, 
     2                       DFDXSP, PROP, XKINSLP1, XKINSLP2, XKINSLP3)

!-----  This subroutine calculates the shear strain-rate in each slip 
!     system for a rate-dependent single crystal.  The POWER LAW 
!     relation between shear strain-rate and resolved shear stress 
!     proposed by Hutchinson, Pan and Rice, is used here.

!-----  The power law exponents are assumed the same for all slip 
!     systems in each set, though they could be different from set to 
!     set, e.g. <110>{111} and <110>{100}.  The strain-rate coefficient
!     in front of the power law form are also assumed the same for all 
!     slip systems in each set. 

!-----  Users who want to use their own constitutive relation may 
!     change the function subprograms F and its derivative DFDX, 
!     where F is the strain hardening law, dGAMMA/dt = F(X), 
!     X=(TAUSLP-XKINSLP)/GSLIP.  The parameters characterizing F are passed into 
!     F and DFDX through array PROP.

!-----  Function subprograms:
!
!       F    -- User-supplied function subprogram which gives shear 
!               strain-rate for each slip system based on current 
!               values of resolved shear stress and current strength
!
!       DFDX -- User-supplied function subprogram dF/dX, where x is the
!               ratio of resolved shear stress over current strength

!-----  Variables:
!
!     GAMMA  -- shear strain in each slip system at the start of time 
!               step  (INPUT)
!     TAUSLP -- resolved shear stress in each slip system (INPUT)
!     GSLIP  -- current strength (INPUT)
!     NSLIP  -- number of slip systems in this set (INPUT)
!
!     FSLIP  -- current value of F for each slip system (OUTPUT)
!     DFDXSP -- current value of DFDX for each slip system (OUTPUT)
!
!     PROP   -- material constants characterizing the strain hardening 
!               law (INPUT)
!
!               For the current power law strain hardening law 
!               PROP(1) -- power law hardening exponent
!               PROP(1) = infinity corresponds to a rate-independent 
!               material
!               PROP(2) -- coefficient in front of power law hardening
!
!     XKINSLP   -- current slip system kinematic hardening variable (INPUT)

!-----  Use single precision on cray
!
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F, DFDX
      DIMENSION GAMMA(NSLIP), TAUSLP(NSLIP), GSLIP(NSLIP), 
     2          FSLIP(NSLIP), DFDXSP(NSLIP), PROP(8),
     3			XKINSLP1(NSLIP), XKINSLP2(NSLIP), XKINSLP3(NSLIP)

      DO I=1,NSLIP       		
         X=(TAUSLP(I)-XKINSLP1(I)-XKINSLP2(I)-XKINSLP3(I))/GSLIP(I)
         FSLIP(I)=F(X,PROP)
         DFDXSP(I)=DFDX(X,PROP)
   
      END DO

      RETURN
      END


!-----------------------------------


!-----  Use single precision on cray
!
           REAL*8 FUNCTION F(X,PROP)

!-----     User-supplied function subprogram which gives shear 
!        strain-rate for each slip system based on current values of 
!        resolved shear stress and current strength
!
!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)

           F=PROP(2)*(ABS(X))**PROP(1)*DSIGN(1.D0,X)

           RETURN
           END


!-----------------------------------


!-----  Use single precision on cray
!
           REAL*8 FUNCTION DFDX(X,PROP)

!-----     User-supplied function subprogram dF/dX, where x is the 
!        ratio of resolved shear stress over current strength

!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION PROP(8)
         
           DFDX=PROP(1)*PROP(2)*(ABS(X))**(PROP(1)-1.)

           RETURN
           END