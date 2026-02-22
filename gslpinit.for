      SUBROUTINE GSLPINIT (GSLIP0, NSLIP, NSLPTL, NSET, PROP)

!-----  This subroutine calculates the initial value of current 
!     strength for each slip system in a rate-dependent single crystal.
!     Two sets of initial values, proposed by Asaro, Pierce et al, and 
!     by Bassani, respectively, are used here.  Both sets assume that 
!     the initial values for all slip systems are the same (initially 
!     isotropic).

!-----  These initial values are assumed the same for all slip systems 
!     in each set, though they could be different from set to set, e.g.
!     <110>{111} and <110>{100}.

!-----  Users who want to use their own initial values may change the 
!     function subprogram GSLP0.  The parameters characterizing these 
!     initial values are passed into GSLP0 through array PROP.

!-----  Use single precision on cray

      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL GSLP0
      DIMENSION GSLIP0(NSLPTL), NSLIP(NSET), PROP(16,NSET)

!-----  Function subprograms:
!
!       GSLP0 -- User-supplied function subprogram given the initial 
!                value of current strength at initial state

!-----  Variables:
!
!     GSLIP0 -- initial value of current strength (OUTPUT)
!
!     NSLIP  -- number of slip systems in each set (INPUT)
!     NSLPTL -- total number of slip systems in all the sets (INPUT)
!     NSET   -- number of sets of slip systems (INPUT)
!
!     PROP   -- material constants characterizing the initial value of 
!               current strength (INPUT)
!
!               For Asaro, Pierce et al's law 
!               PROP(1,i) -- initial hardening modulus H0 in the ith 
!                            set of slip systems
!               PROP(2,i) -- saturation stress TAUs in the ith set of  
!                            slip systems
!               PROP(3,i) -- initial critical resolved shear stress 
!                            TAU0 in the ith set of slip systems
!
!               For Bassani's law 
!               PROP(1,i) -- initial hardening modulus H0 in the ith 
!                            set of slip systems
!               PROP(2,i) -- stage I stress TAUI in the ith set of  
!                            slip systems (or the breakthrough stress 
!                            where large plastic flow initiates)
!               PROP(3,i) -- initial critical resolved shear stress 
!                            TAU0 in the ith set of slip systems


      ID=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ID=ID+1
            GSLIP0(ID)=GSLP0(NSLPTL,NSET,NSLIP,PROP(1,I),ID,ISET)
         END DO
      END DO

      RETURN
      END


!----------------------------------


!-----  Use single precision on cray

           REAL*8 FUNCTION GSLP0(NSLPTL,NSET,NSLIP,PROP,ISLIP,ISET)

!-----     User-supplied function subprogram given the initial value of
!          current strength at initial state

!-----  Use single precision on cray

           IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION NSLIP(NSET), PROP(16)

           GSLP0=PROP(3)

           RETURN
           END