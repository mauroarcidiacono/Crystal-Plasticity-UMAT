      SUBROUTINE CORE(STRESS,STATEV,DDSDDE,
     1 STRAN,DSTRAN,DTIME, ND, NOEL,
     2 NDI,NSHR,NTENS,NSTATV,LPROPS,DROT)

      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER NDI, NSHR, NTENS, NSTATV, RMID

      REAL*8 STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS), DTIME,
     2 STRAN(NTENS),DSTRAN(NTENS),
     3 LPROPS(184),DROT(3,3)

      DIMENSION ISPDIR(3), ISPNOR(3), NSLIP(3), 
     2          SLPDIR(3,ND), SLPNOR(3,ND), SLPDEF(6,ND), 
     3          SLPSPN(3,ND), DSPDIR(3,ND), DSPNOR(3,ND), 
     4          DLOCAL(6,6), D(6,6), ROTD(6,6), ROTATE(3,3), 
     5          FSLIP(ND), DFDXSP(ND), DDEMSD(6,ND), 
     6          H(ND,ND), DDGDDE(ND,6), 
     7          DSTRES(6), DELATS(6), DSPIN(3), DVGRAD(3,3),
     8          DGAMMA(ND), DTAUSP(ND), DGSLIP(ND),
     9          WORKST(ND,ND), INDX(ND), TERM(3,3), TRM0(3,3), ITRM(3)


      DIMENSION FSLIP1(ND), STRES1(6), GAMMA1(ND), TAUSP1(ND),
     2          GSLP1(ND), SPNOR1(3,ND), SPDIR1(3,ND), DDSDE1(6,6),
     3          DSOLD(6), DGAMOD(ND), DTAUOD(ND), DGSPOD(ND), 
     4          DSPNRO(3,ND), DSPDRO(3,ND), 
     5          DHDGDG(ND,ND)


      DIMENSION PLDT(3,3), XKINSLP11(ND), XKINSLP12(ND), XKINSLP13(ND), 
     2          DKINSLP1(ND), DKINOD1(ND),
     3          DKINSLP2(ND), DKINOD2(ND), DKINSLP3(ND), DKINOD3(ND),
     4          ROTM1(3,3), ROTM2(3,3), HSLPNOR(3,ND), HSLPDIR(3,ND),
     5          T1(3), T2(3), CBM(3,3)


      !-----  Elastic matrix in local cubic crystal system: DLOCAL
      DO J=1,6
         DO I=1,6
            DLOCAL(I,J)=0.
         END DO
      END DO

      CHECK=0.
      DO J=10,21
         CHECK=CHECK+ABS(LPROPS(J))
      END DO

      IF (CHECK.EQ.0.) THEN
         DO J=4,9
            CHECK=CHECK+ABS(LPROPS(J))
         END DO

         IF (CHECK.EQ.0.) THEN

            IF (LPROPS(3).EQ.0.) THEN

!-----  Isotropic material
               GSHEAR=LPROPS(1)/2./(1.+LPROPS(2))
               E11=2.*GSHEAR*(1.-LPROPS(2))/(1.-2.*LPROPS(2))
               E12=2.*GSHEAR*LPROPS(2)/(1.-2.*LPROPS(2))

               DO J=1,3
                  DLOCAL(J,J)=E11

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=E12
                  END DO

                  DLOCAL(J+3,J+3)=GSHEAR
               END DO

            ELSE

!-----  Cubic material
               DO J=1,3
                  DLOCAL(J,J)=LPROPS(1)

                  DO I=1,3
                     IF (I.NE.J) DLOCAL(I,J)=LPROPS(2)
                  END DO

                  DLOCAL(J+3,J+3)=LPROPS(3)
               END DO

            END IF

         ELSE

!-----  Orthotropic material
            DLOCAL(1,1)=LPROPS(1)
            DLOCAL(1,2)=LPROPS(2)
            DLOCAL(2,1)=LPROPS(2)
            DLOCAL(2,2)=LPROPS(3)

            DLOCAL(1,3)=LPROPS(4)
            DLOCAL(3,1)=LPROPS(4)
            DLOCAL(2,3)=LPROPS(5)
            DLOCAL(3,2)=LPROPS(5)
            DLOCAL(3,3)=LPROPS(6)

            DLOCAL(4,4)=LPROPS(7)
            DLOCAL(5,5)=LPROPS(8)
            DLOCAL(6,6)=LPROPS(9)

         END IF

      ELSE

!-----  General anisotropic material
         ID=0
         DO J=1,6
            DO I=1,J
               ID=ID+1
               DLOCAL(I,J)=LPROPS(ID)
               DLOCAL(J,I)=DLOCAL(I,J)
            END DO
         END DO
      END IF

!-----  Assign the kinematic hardening parameters

      CKIN1 = LPROPS(161)
      DKIN1 = LPROPS(162)
      CKIN2 = LPROPS(169)
      DKIN2 = LPROPS(170)
      CKIN3 = LPROPS(177)
      DKIN3 = LPROPS(178)

!----- Rotation matrix: ROTATE, i.e. direction cosines of the local 
!      orthonormal basis vectors ([100], [010], [001]) expressed 
!      in the global coordinate system.

      CALL ROTATION (LPROPS(57), ROTATE)

!-----  Store the initial rotation matrix and calculate the current one.
!       The initial rotation matrix is stored from SDV482 to SDV490.
!       The current rotation matrix is stored from SDV491 to SDV499.
 
      RMID = 491

      IF (STATEV(1).EQ.0) THEN

         IRMID=482

         DO I=1,3
            DO J=1,3
               STATEV(IRMID)=ROTATE(I,J)
               IRMID=IRMID+1
               STATEV(RMID)=ROTATE(I,J)
               RMID=RMID+1
            END DO
         END DO

      ELSE

         DO I=1,3
            DO J=1,3
               ROTM1(I,J)=STATEV(RMID)
               RMID=RMID+1
            END DO
         END DO

         DO I=1,3
            DO J=1,3
               ROTM2(I,J)=0
               DO L=1,3
                  ROTM2(I,J)=ROTM2(I,J)+DROT(I,L)*ROTM1(L,J)
               END DO
            END DO
         END DO

         RMID=491

         DO I=1,3
            DO J=1,3
               STATEV(RMID)=ROTM2(I,J)
               RMID=RMID+1
            END DO
         END DO
      END IF

!-----  Rotation matrix: ROTD to transform local elastic matrix DLOCAL 
!       to global elastic matrix D

      DO J=1,3
         J1=1+J/3
         J2=2+J/2

         DO I=1,3
            I1=1+I/3
            I2=2+I/2

            ROTD(I,J)=ROTATE(I,J)**2
            ROTD(I,J+3)=2.*ROTATE(I,J1)*ROTATE(I,J2)
            ROTD(I+3,J)=ROTATE(I1,J)*ROTATE(I2,J)
            ROTD(I+3,J+3)=ROTATE(I1,J1)*ROTATE(I2,J2)+
     2                    ROTATE(I1,J2)*ROTATE(I2,J1)

         END DO
      END DO

!-----  Elastic matrix in global system: D
!     {D} = {ROTD} * {DLOCAL} * {ROTD}transpose

      DO J=1,6
         DO I=1,6
            D(I,J)=0.
         END DO
      END DO

      DO J=1,6
         DO I=1,J

            DO K=1,6
               DO L=1,6
                  D(I,J)=D(I,J)+DLOCAL(K,L)*ROTD(I,K)*ROTD(J,L)
               END DO
            END DO

            D(J,I)=D(I,J)

         END DO
      END DO

!-----  Total number of sets of slip systems: NSET
      NSET=NINT(LPROPS(25))
      IF (NSET.LT.1) THEN
         WRITE (6,*) '***ERROR - zero sets of slip systems'
         CALL XIT
      ELSE IF (NSET.GT.3) THEN
         WRITE (6,*) 
     2     '***ERROR - more than three sets of slip systems'
         CALL XIT
      END IF

!-----  Implicit integration parameter: THETA
      THETA=LPROPS(145)

!-----  Finite deformation ?
!-----  NLGEOM = 0,   small deformation theory
!       otherwise, theory of finite rotation and finite strain, Users 
!       must declare "NLGEOM" in the input file, at the *STEP card

      IF (LPROPS(146).EQ.0.) THEN
         NLGEOM=0
      ELSE
         NLGEOM=1
      END IF

!-----  Iteration?
!-----  ITRATN = 0, no iteration
!       otherwise, iteration (solving increments of stresses and 
!       solution dependent state variables)

      IF (LPROPS(153).EQ.0.) THEN
         ITRATN=0
      ELSE
         ITRATN=1
      END IF

      ITRMAX=NINT(LPROPS(154))
      GAMERR=LPROPS(155)

      NITRTN=-1

      DO I=1,NTENS
         DSOLD(I)=0.
      END DO

      DO J=1,ND
         DGAMOD(J)=0.
         DTAUOD(J)=0.
         DGSPOD(J)=0.
         DKINOD1(J)=0.
         DKINOD2(J)=0.
         DKINOD3(J)=0.
         DO I=1,3
            DSPNRO(I,J)=0.
            DSPDRO(I,J)=0.
         END DO
      END DO

!-----  Increment of spin associated with the material element: DSPIN
!     (only needed for finite rotation)

      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               TERM(I,J)=DROT(J,I) ! TERM = DROT^T
               TRM0(I,J)=DROT(J,I) ! TRM0 = DROT^T
            END DO

            TERM(J,J)=TERM(J,J)+1.D0 ! TERM = DROT^T + I
            TRM0(J,J)=TRM0(J,J)-1.D0 ! TRM0 = DROT^T - I
         END DO

         CALL LUDCMP (TERM, 3, 3, ITRM, DDCMP, NOEL)

         DO J=1,3
            CALL LUBKSB (TERM, 3, 3, ITRM, TRM0(1,J))
         END DO

         ! Axial vector of the skew-symmetric spin tensor
         DSPIN(1)=TRM0(2,1)-TRM0(1,2) ! w_x from skew part
         DSPIN(2)=TRM0(1,3)-TRM0(3,1) ! w_y
         DSPIN(3)=TRM0(3,2)-TRM0(2,3) ! w_z

      END IF

!-----  Increment of dilatational strain: DEV
      DEV=0.D0
      DO I=1,NDI
         DEV=DEV+DSTRAN(I)
      END DO

!-----  Iteration starts (only when iteration method is used)
1000  CONTINUE

!-----  Parameter NITRTN: number of iterations
!       NITRTN = 0 --- no-iteration solution

      NITRTN=NITRTN+1

!-----  Check whether the current stress state is the initial state
      IF (STATEV(1).EQ.0.) THEN

!-----  Initial state

!-----  Generating the following parameters and variables at initial 
!     state:
!          Total number of slip systems in all the sets NSLPTL
!          Number of slip systems in each set NSLIP
!          Unit vectors in initial slip directions SLPDIR
!          Unit normals to initial slip planes SLPNOR

         NSLPTL=0

!-----  Update of the UMAT to consider the HCP crystal structure.
!       The code introduces the basal, pyramidal and prismatic slip planes. 

         IF (LPROPS(9).NE.0) THEN   
            CALL NONCUBICCRYSTAL (SLPDIR, SLPNOR, ROTATE, NSLIP, 
     2                         NSLPTL)
         ELSE
            DO I=1,NSET
               ISPNOR(1)=NINT(LPROPS(25+8*I))
               ISPNOR(2)=NINT(LPROPS(26+8*I))
               ISPNOR(3)=NINT(LPROPS(27+8*I))

               ISPDIR(1)=NINT(LPROPS(28+8*I))
               ISPDIR(2)=NINT(LPROPS(29+8*I))
               ISPDIR(3)=NINT(LPROPS(30+8*I))

         CALL SLIPSYS (ISPDIR, ISPNOR, NSLIP(I), SLPDIR(1,NSLPTL+1), 
     2                    SLPNOR(1,NSLPTL+1), ROTATE)

               NSLPTL=NSLPTL+NSLIP(I)
            END DO
         END IF

         IF (ND.LT.NSLPTL) THEN
            WRITE (6,*) 
     2 '***ERROR - parameter ND chosen by the present user is less than
     3             the total number of slip systems NSLPTL'
            WRITE (6,*) ND
            WRITE (6,*) NSLPTL
            WRITE (6,*) NSLIP(1)
            CALL XIT
         END IF

!-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

!-----  Initial value of state variables: unit normal to a slip plane 
!       and unit vector in a slip direction

         STATEV(NSTATV)=FLOAT(NSLPTL)
         DO I=1,NSET
            STATEV(NSTATV-3*(NSLPTL+2)+I)=FLOAT(NSLIP(I))
            !STATEV(NSTATV-(5+NSLPTL)+I)=FLOAT(NSLIP(I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=SLPNOR(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=SLPDIR(I,J)
            END DO
         END DO

!-----  Initial value of the current strength for all slip systems

         CALL GSLPINIT (STATEV(1), NSLIP, NSLPTL, NSET, LPROPS(97))

!-----  Initial value of shear strain in slip systems
!-----  Initial value of cumulative shear strain in each slip systems

         DO I=1,NSLPTL
            STATEV(NSLPTL+I)=0.
            STATEV(9*NSLPTL+I)=0.
         END DO

         STATEV(500)=0.

!-----  Initial value of the resolved shear stress in slip systems
         DO I=1,NSLPTL
            TERM1=0.

            DO J=1,NTENS
               IF (J.LE.NDI) THEN
                  TERM1=TERM1+SLPDEF(J,I)*STRESS(J)
               ELSE
                  TERM1=TERM1+SLPDEF(J-NDI+3,I)*STRESS(J)
               END IF
            END DO

            STATEV(2*NSLPTL+I)=TERM1
         END DO

!-----  Initial value of the kinematic variable and its increment in slip systems
         DO I=1,NSLPTL
            STATEV(NSTATV-3*(NSLPTL+1)+I)=0.0
            DKINSLP1(I)=0.0
         END DO

         DO I=1,NSLPTL
            STATEV(NSTATV-2*(NSLPTL+1)+I)=0.0
            DKINSLP2(I)=0.0
         END DO

         DO I=1,NSLPTL
            STATEV(NSTATV-(NSLPTL+1)+I)=0.0
            DKINSLP3(I)=0.0
         END DO
		
      ELSE

!-----  Current stress state
!
!-----  Copying from the array of state variables STATVE the following
!          parameters and variables at current stress state:
!          Total number of slip systems in all the sets NSLPTL
!          Number of slip systems in each set NSLIP
!          Current slip directions SLPDIR
!          Normals to current slip planes SLPNOR
!
        
         NSLPTL=NINT(STATEV(NSTATV))
         DO I=1,NSET
            NSLIP(I)=NINT(STATEV(NSTATV-3*(NSLPTL+2)+I))
            !NSLIP(I)=NINT(STATEV(NSTATV-(5+NSLPTL)+I))
         END DO

         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               SLPNOR(I,J)=STATEV(IDNOR)

               IDDIR=IDDIR+1
               SLPDIR(I,J)=STATEV(IDDIR)
            END DO
         END DO

!-----  Slip deformation tensor: SLPDEF (Schmid factors)
         DO J=1,NSLPTL
            SLPDEF(1,J)=SLPDIR(1,J)*SLPNOR(1,J)
            SLPDEF(2,J)=SLPDIR(2,J)*SLPNOR(2,J)
            SLPDEF(3,J)=SLPDIR(3,J)*SLPNOR(3,J)
            SLPDEF(4,J)=SLPDIR(1,J)*SLPNOR(2,J)+SLPDIR(2,J)*SLPNOR(1,J)
            SLPDEF(5,J)=SLPDIR(1,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(1,J)
            SLPDEF(6,J)=SLPDIR(2,J)*SLPNOR(3,J)+SLPDIR(3,J)*SLPNOR(2,J)
         END DO

      END IF

!-----  Slip spin tensor: SLPSPN (only needed for finite rotation)
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            SLPSPN(1,J)=0.5*(SLPDIR(1,J)*SLPNOR(2,J)-
     2                       SLPDIR(2,J)*SLPNOR(1,J))
            SLPSPN(2,J)=0.5*(SLPDIR(3,J)*SLPNOR(1,J)-
     2                       SLPDIR(1,J)*SLPNOR(3,J))
            SLPSPN(3,J)=0.5*(SLPDIR(2,J)*SLPNOR(3,J)-
     2                       SLPDIR(3,J)*SLPNOR(2,J))
         END DO
      END IF

!-----  Double dot product of elastic moduli tensor with the slip 
!     deformation tensor (Schmid factors) plus, only for finite 
!     rotation, the dot product of slip spin tensor with the stress: 
!     DDEMSD
!
      DO J=1,NSLPTL
         DO I=1,6
            DDEMSD(I,J)=0.
            DO K=1,6
               DDEMSD(I,J)=DDEMSD(I,J)+D(K,I)*SLPDEF(K,J)
            END DO
         END DO
      END DO

      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL

            DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(1,J)*STRESS(1)
            DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(2,J)*STRESS(1)

            IF (NDI.GT.1) THEN
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(1,J)*STRESS(2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(3,J)*STRESS(2)
            END IF

            IF (NDI.GT.2) THEN
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(2,J)*STRESS(3)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(3,J)*STRESS(3)
            END IF

            IF (NSHR.GE.1) THEN
               DDEMSD(1,J)=DDEMSD(1,J)+SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(2,J)=DDEMSD(2,J)-SLPSPN(1,J)*STRESS(NDI+1)
               DDEMSD(5,J)=DDEMSD(5,J)-SLPSPN(3,J)*STRESS(NDI+1)
               DDEMSD(6,J)=DDEMSD(6,J)+SLPSPN(2,J)*STRESS(NDI+1)
            END IF

            IF (NSHR.GE.2) THEN
               DDEMSD(1,J)=DDEMSD(1,J)-SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(3,J)=DDEMSD(3,J)+SLPSPN(2,J)*STRESS(NDI+2)
               DDEMSD(4,J)=DDEMSD(4,J)+SLPSPN(3,J)*STRESS(NDI+2)
               DDEMSD(6,J)=DDEMSD(6,J)-SLPSPN(1,J)*STRESS(NDI+2)
            END IF

            IF (NSHR.EQ.3) THEN
               DDEMSD(2,J)=DDEMSD(2,J)+SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(3,J)=DDEMSD(3,J)-SLPSPN(3,J)*STRESS(NDI+3)
               DDEMSD(4,J)=DDEMSD(4,J)-SLPSPN(2,J)*STRESS(NDI+3)
               DDEMSD(5,J)=DDEMSD(5,J)+SLPSPN(1,J)*STRESS(NDI+3)
            END IF

         END DO
      END IF

!-----  Shear strain-rate in a slip system at the start of increment: 
!     FSLIP, and its derivative: DFDXSP
!	
      
      ID=1
      DO I=1,NSET
         IF (I.GT.1) ID=ID+NSLIP(I-1)
         CALL STRAINRATE (STATEV(NSLPTL+ID), STATEV(2*NSLPTL+ID), 
     2                    STATEV(ID), NSLIP(I), FSLIP(ID), DFDXSP(ID), 
     3                    LPROPS(65+8*I),STATEV(NSTATV-3*(NSLPTL+1)+ID),
     4   STATEV(NSTATV-2*(NSLPTL+1)+ID), STATEV(NSTATV-(NSLPTL+1)+ID))
      END DO
      
!-----  Self- and latent-hardening laws
  
       CALL LATENTHARDEN (STATEV(NSLPTL+1), STATEV(2*NSLPTL+1), 
     2                   STATEV(1), STATEV(9*NSLPTL+1),
     3                   STATEV(500), NSLIP, NSLPTL, 
     4                   NSET, H(1,1), LPROPS(97), ND)


!-----  LU decomposition to solve the increment of shear strain in a 
!     slip system
!
      TERM1=THETA*DTIME

      DO I=1,NSLPTL
         TAUSLP=STATEV(2*NSLPTL+I)
         GSLIP=STATEV(I)
         XKINSLP1=STATEV(NSTATV-3*(NSLPTL+1)+I)
         XKINSLP2=STATEV(NSTATV-2*(NSLPTL+1)+I)
         XKINSLP3=STATEV(NSTATV-(NSLPTL+1)+I)
         X=(TAUSLP-XKINSLP1-XKINSLP2-XKINSLP3)/GSLIP
         TERM2=TERM1*DFDXSP(I)/GSLIP
         TERM3=TERM1*X*DFDXSP(I)/GSLIP
         TERM4=-TERM1*DFDXSP(I)/GSLIP

         DO J=1,NSLPTL
            TERM5=0.
            DO K=1,6
               TERM5=TERM5+DDEMSD(K,I)*SLPDEF(K,J)
            END DO

            WORKST(I,J)=TERM2*TERM5+H(I,J)*TERM3*DSIGN(1.D0,FSLIP(J))

            IF (NITRTN.GT.0) WORKST(I,J)=WORKST(I,J)+TERM3*DHDGDG(I,J)

         END DO
         
         WORKST(I,I)=WORKST(I,I)+1.
     2		-TERM4*((CKIN1)-((DKIN1*XKINSLP1)*DSIGN(1.D0,FSLIP(I))))
     3      -TERM4*((CKIN2)-((DKIN2*XKINSLP2)*DSIGN(1.D0,FSLIP(I))))
     4      -TERM4*((CKIN3)-((DKIN3*XKINSLP3)*DSIGN(1.D0,FSLIP(I))))
      END DO
      
      CALL LUDCMP (WORKST, NSLPTL, ND, INDX, DDCMP, NOEL)

!-----  Increment of shear strain in a slip system: DGAMMA
      TERM1=THETA*DTIME
      DO I=1,NSLPTL

         IF (NITRTN.EQ.0) THEN
            TAUSLP=STATEV(2*NSLPTL+I)
            GSLIP=STATEV(I)
            XKINSLP1=STATEV(NSTATV-3*(NSLPTL+1)+I)
            XKINSLP2=STATEV(NSTATV-2*(NSLPTL+1)+I)
            XKINSLP3=STATEV(NSTATV-(NSLPTL+1)+I)
         	X=(TAUSLP-XKINSLP1-XKINSLP2-XKINSLP3)/GSLIP
            TERM2=TERM1*DFDXSP(I)/GSLIP

            DGAMMA(I)=0.
            DO J=1,NDI
               DGAMMA(I)=DGAMMA(I)+DDEMSD(J,I)*DSTRAN(J)
            END DO

            IF (NSHR.GT.0) THEN
               DO J=1,NSHR
                  DGAMMA(I)=DGAMMA(I)+DDEMSD(J+3,I)*DSTRAN(J+NDI)
               END DO
            END IF

            DGAMMA(I)=DGAMMA(I)*TERM2+FSLIP(I)*DTIME

         ELSE
            DGAMMA(I)=TERM1*(FSLIP(I)-FSLIP1(I))+FSLIP1(I)*DTIME
     2                -DGAMOD(I)

         END IF

      END DO

      CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DGAMMA)

      DO I=1,NSLPTL
         DGAMMA(I)=DGAMMA(I)+DGAMOD(I)
      END DO

!-----  Update the shear strain in a slip system: STATEV(NSLPTL+1) - 
!     STATEV(2*NSLPTL)
!
      DO I=1,NSLPTL
         STATEV(NSLPTL+I)=STATEV(NSLPTL+I)+DGAMMA(I)-DGAMOD(I)
      END DO 

!-----  Increment of current strength in a slip system: DGSLIP
      DO I=1,NSLPTL
         DGSLIP(I)=0.
         DO J=1,NSLPTL
            DGSLIP(I)=DGSLIP(I)+H(I,J)*ABS(DGAMMA(J))
         END DO
      END DO

!-----  Update the current strength in a slip system: STATEV(1) - 
!     STATEV(NSLPTL)
!
      DO I=1,NSLPTL
         STATEV(I)=STATEV(I)+DGSLIP(I)-DGSPOD(I)
      END DO

!-----  Increment of kinematic variable in a slip system: DKINSLP
!     

      ID=1
      DO I=1,NSET
         IF (I.GT.1) ID=ID+NSLIP(I-1)
         CALL KINHARDEN (LPROPS(161), LPROPS(162), LPROPS(169),   
     2      LPROPS(170), LPROPS(177), LPROPS(178),DGAMMA(ID),
     3	   STATEV(NSTATV-3*(NSLPTL+1)+ID),  
     4      STATEV(NSTATV-2*(NSLPTL+1)+ID), 
     5      STATEV(NSTATV-(NSLPTL+1)+ID), DKINSLP1(ID), 
     6      DKINSLP2(ID), DKINSLP3(ID), NSLIP(I)) 
     
      END DO

!-----  Update state variables for kinematic hardening variable
!
      DO I=1,NSLPTL
         STATEV(NSTATV-3*(NSLPTL+1)+I)=STATEV(NSTATV-3*(NSLPTL+1)+I)-
     2		DKINOD1(I)+DKINSLP1(I)
      END DO

      DO I=1,NSLPTL
         STATEV(NSTATV-2*(NSLPTL+1)+I)=STATEV(NSTATV-2*(NSLPTL+1)+I)-
     2		DKINOD2(I)+DKINSLP2(I)
      END DO

      DO I=1,NSLPTL
         STATEV(NSTATV-(NSLPTL+1)+I)=STATEV(NSTATV-(NSLPTL+1)+I)-
     2		DKINOD3(I)+DKINSLP3(I)
      END DO

!-----  Increment of strain associated with lattice stretching: DELATS
      DO J=1,6
         DELATS(J)=0.
      END DO

      DO J=1,3
         IF (J.LE.NDI) DELATS(J)=DSTRAN(J)
         DO I=1,NSLPTL
            DELATS(J)=DELATS(J)-SLPDEF(J,I)*DGAMMA(I)
         END DO
      END DO

      DO J=1,3
         IF (J.LE.NSHR) DELATS(J+3)=DSTRAN(J+NDI)
         DO I=1,NSLPTL
            DELATS(J+3)=DELATS(J+3)-SLPDEF(J+3,I)*DGAMMA(I)
         END DO
      END DO

!-----  Increment of deformation gradient associated with lattice 
!     stretching in the current state, i.e. the velocity gradient 
!     (associated with lattice stretching) times the increment of time:
!     DVGRAD (only needed for finite rotation)
!
      IF (NLGEOM.NE.0) THEN
         DO J=1,3
            DO I=1,3
               IF (I.EQ.J) THEN
                  DVGRAD(I,J)=DELATS(I)
               ELSE
                  DVGRAD(I,J)=DELATS(I+J+1)
               END IF
            END DO
         END DO

         DO J=1,3
            DO I=1,J
               IF (J.GT.I) THEN
                  IJ2=I+J-2
                  IF (MOD(IJ2,2).EQ.1) THEN
                     TERM1=1.
                  ELSE
                     TERM1=-1.
                  END IF

                  DVGRAD(I,J)=DVGRAD(I,J)+TERM1*DSPIN(IJ2)
                  DVGRAD(J,I)=DVGRAD(J,I)-TERM1*DSPIN(IJ2)

                  DO K=1,NSLPTL
                     DVGRAD(I,J)=DVGRAD(I,J)-TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                     DVGRAD(J,I)=DVGRAD(J,I)+TERM1*DGAMMA(K)*
     2                                       SLPSPN(IJ2,K)
                  END DO
               END IF

            END DO
         END DO

      END IF

!-----  Increment of resolved shear stress in a slip system: DTAUSP
      DO I=1,NSLPTL
         DTAUSP(I)=0.
         DO J=1,6
            DTAUSP(I)=DTAUSP(I)+DDEMSD(J,I)*DELATS(J)
         END DO
      END DO

!-----  Update the resolved shear stress in a slip system: 
!     STATEV(2*NSLPTL+1) - STATEV(3*NSLPTL)
!
      DO I=1,NSLPTL
         STATEV(2*NSLPTL+I)=STATEV(2*NSLPTL+I)+DTAUSP(I)-DTAUOD(I)
      END DO

!-----  Increment of stress: DSTRES
      IF (NLGEOM.EQ.0) THEN
         DO I=1,NTENS
            DSTRES(I)=0.
         END DO
      ELSE
         DO I=1,NTENS
            DSTRES(I)=-STRESS(I)*DEV
         END DO
      END IF

      DO I=1,NDI
         DO J=1,NDI
            DSTRES(I)=DSTRES(I)+D(I,J)*DSTRAN(J)
         END DO

         IF (NSHR.GT.0) THEN
            DO J=1,NSHR
               DSTRES(I)=DSTRES(I)+D(I,J+3)*DSTRAN(J+NDI)
            END DO
         END IF

         DO J=1,NSLPTL
            DSTRES(I)=DSTRES(I)-DDEMSD(I,J)*DGAMMA(J)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO I=1,NSHR

            DO J=1,NDI
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J)*DSTRAN(J)
            END DO

            DO J=1,NSHR
               DSTRES(I+NDI)=DSTRES(I+NDI)+D(I+3,J+3)*DSTRAN(J+NDI)
            END DO

            DO J=1,NSLPTL
               DSTRES(I+NDI)=DSTRES(I+NDI)-DDEMSD(I+3,J)*DGAMMA(J)
            END DO

         END DO
      END IF

!-----  Update the stress: STRESS
      DO I=1,NTENS
         STRESS(I)=STRESS(I)+DSTRES(I)-DSOLD(I)
      END DO

!-----  Increment of normal to a slip plane and a slip direction (only 
!    needed for finite rotation)
!
      IF (NLGEOM.NE.0) THEN
         DO J=1,NSLPTL
            DO I=1,3
               DSPNOR(I,J)=0.
               DSPDIR(I,J)=0.

               DO K=1,3
                  DSPNOR(I,J)=DSPNOR(I,J)-SLPNOR(K,J)*DVGRAD(K,I)
                  DSPDIR(I,J)=DSPDIR(I,J)+SLPDIR(K,J)*DVGRAD(I,K)
               END DO

            END DO
         END DO

!-----  Update the normal to a slip plane and a slip direction (only 
!     needed for finite rotation)
!
         IDNOR=3*NSLPTL
         IDDIR=6*NSLPTL
         DO J=1,NSLPTL
            DO I=1,3
               IDNOR=IDNOR+1
               STATEV(IDNOR)=STATEV(IDNOR)+DSPNOR(I,J)-DSPNRO(I,J)

               IDDIR=IDDIR+1
               STATEV(IDDIR)=STATEV(IDDIR)+DSPDIR(I,J)-DSPDRO(I,J)
            END DO
         END DO

      END IF

!-----  Derivative of shear strain increment in a slip system w.r.t. 
!     strain increment: DDGDDE
!
      TERM1=THETA*DTIME
      DO I=1,NTENS
         DO J=1,NSLPTL
            TAUSLP=STATEV(2*NSLPTL+J)
            GSLIP=STATEV(J)
            XKINSLP1=STATEV(NSTATV-3*(NSLPTL+1)+I)
            XKINSLP2=STATEV(NSTATV-2*(NSLPTL+1)+I)
            XKINSLP3=STATEV(NSTATV-(NSLPTL+1)+I)
         	X=(TAUSLP-XKINSLP1-XKINSLP2-XKINSLP3)/GSLIP
            TERM2=TERM1*DFDXSP(J)/GSLIP
            IF (I.LE.NDI) THEN
               DDGDDE(J,I)=TERM2*DDEMSD(I,J)
            ELSE
               DDGDDE(J,I)=TERM2*DDEMSD(I-NDI+3,J)
            END IF
         END DO

         CALL LUBKSB (WORKST, NSLPTL, ND, INDX, DDGDDE(1,I))

      END DO

!-----  Derivative of stress increment w.r.t. strain increment, i.e. 
!     Jacobian matrix
!
!-----  Jacobian matrix: elastic part
      DO J=1,NTENS
         DO I=1,NTENS
            DDSDDE(I,J)=0.
         END DO
      END DO

      DO J=1,NDI
         DO I=1,NDI
            DDSDDE(I,J)=D(I,J)
            IF (NLGEOM.NE.0) DDSDDE(I,J)=DDSDDE(I,J)-STRESS(I)
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR
            DO I=1,NSHR
               DDSDDE(I+NDI,J+NDI)=D(I+3,J+3)
            END DO

            DO I=1,NDI
               DDSDDE(I,J+NDI)=D(I,J+3)
               DDSDDE(J+NDI,I)=D(J+3,I)
               IF (NLGEOM.NE.0)
     2            DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-STRESS(J+NDI)
            END DO
         END DO
      END IF

!-----  Jacobian matrix: plastic part (slip)
      DO J=1,NDI
         DO I=1,NDI
            DO K=1,NSLPTL
               DDSDDE(I,J)=DDSDDE(I,J)-DDEMSD(I,K)*DDGDDE(K,J)
            END DO
         END DO
      END DO

      IF (NSHR.GT.0) THEN
         DO J=1,NSHR

            DO I=1,NSHR
               DO K=1,NSLPTL
                  DDSDDE(I+NDI,J+NDI)=DDSDDE(I+NDI,J+NDI)-
     2                                DDEMSD(I+3,K)*DDGDDE(K,J+NDI)
               END DO
            END DO

            DO I=1,NDI
               DO K=1,NSLPTL
                  DDSDDE(I,J+NDI)=DDSDDE(I,J+NDI)-
     2                            DDEMSD(I,K)*DDGDDE(K,J+NDI)
                  DDSDDE(J+NDI,I)=DDSDDE(J+NDI,I)-
     2                            DDEMSD(J+3,K)*DDGDDE(K,I)
               END DO
            END DO

         END DO
      END IF

      IF (ITRATN.NE.0) THEN
         DO J=1,NTENS
            DO I=1,NTENS
               DDSDDE(I,J)=DDSDDE(I,J)/(1.+DEV)
            END DO
         END DO
      END IF

!-----  Iteration ?
      IF (ITRATN.NE.0) THEN

!-----  Save solutions (without iteration):
!            Shear strain-rate in a slip system FSLIP1
!            Current strength in a slip system GSLP1
!            Shear strain in a slip system GAMMA1
!            Resolved shear stress in a slip system TAUSP1
!			    Kinematic variable in a slip system XKINSLP1
!            Normal to a slip plane SPNOR1
!            Slip direction SPDIR1
!            Stress STRES1
!            Jacobian matrix DDSDE1
!
         IF (NITRTN.EQ.0) THEN

            IDNOR=3*NSLPTL
            IDDIR=6*NSLPTL
            DO J=1,NSLPTL
               FSLIP1(J)=FSLIP(J)
               GSLP1(J)=STATEV(J)
               GAMMA1(J)=STATEV(NSLPTL+J)
               TAUSP1(J)=STATEV(2*NSLPTL+J)
               XKINSLP11(J)=STATEV(NSTATV-3*(NSLPTL+1)+J)
               XKINSLP12(J)=STATEV(NSTATV-2*(NSLPTL+1)+J)
               XKINSLP13(J)=STATEV(NSTATV-(NSLPTL+1)+J)
               DO I=1,3
                  IDNOR=IDNOR+1
                  SPNOR1(I,J)=STATEV(IDNOR)

                  IDDIR=IDDIR+1
                  SPDIR1(I,J)=STATEV(IDDIR)
               END DO
            END DO

            DO J=1,NTENS
               STRES1(J)=STRESS(J)
               DO I=1,NTENS
                  DDSDE1(I,J)=DDSDDE(I,J)
               END DO
            END DO

         END IF

!-----  Increments of stress DSOLD, and solution dependent state 
!     variables DGAMOD, DTAUOD, DGSPOD,DKINSLP, DSPNRO, DSPDRO 
!	  (for the next iteration)
!
         DO I=1,NTENS
            DSOLD(I)=DSTRES(I)
         END DO

         DO J=1,NSLPTL
            DGAMOD(J)=DGAMMA(J)
            DTAUOD(J)=DTAUSP(J)
            DGSPOD(J)=DGSLIP(J)
            DKINOD1(J)=DKINSLP1(J)
            DKINOD2(J)=DKINSLP2(J)
            DKINOD3(J)=DKINSLP3(J)
            DO I=1,3
               DSPNRO(I,J)=DSPNOR(I,J)
               DSPDRO(I,J)=DSPDIR(I,J)
            END DO
         END DO

!-----  Check if the iteration solution converges
         IDBACK=0
         ID=0
         DO I=1,NSET
            DO J=1,NSLIP(I)
               ID=ID+1
               TERM1=STATEV(2*NSLPTL+ID)-STATEV(NSTATV-3*(NSLPTL+1)+ID)-
     2       STATEV(NSTATV-2*(NSLPTL+1)+ID)-STATEV(NSTATV-(NSLPTL+1)+ID)
               X=TERM1/STATEV(ID)
               RESIDU=THETA*DTIME*F(X,LPROPS(65+8*I))+DTIME*(1.0-THETA)*
     2                FSLIP1(ID)-DGAMMA(ID)
               IF (ABS(RESIDU).GT.GAMERR) IDBACK=1
            END DO
         END DO

         IF (IDBACK.NE.0.AND.NITRTN.LT.ITRMAX) THEN
!-----  Iteration: arrays for iteration

            CALL ITERATION (STATEV(NSLPTL+1), STATEV(2*NSLPTL+1), 
     2                      STATEV(1), STATEV(9*NSLPTL+1), 
     3                      STATEV(500), NSLPTL, 
     4                      NSET, NSLIP, ND, LPROPS(97), DGAMOD,
     5                      DHDGDG)


            GO TO 1000

         ELSE IF (NITRTN.GE.ITRMAX) THEN
!-----  Solution not converge within maximum number of iteration (the 
!     solution without iteration will be used)
!
!           WRITE(6,*) 'The solution does not converge'
            DO J=1,NTENS
               STRESS(J)=STRES1(J)
               DO I=1,NTENS
                  DDSDDE(I,J)=DDSDE1(I,J)
               END DO
            END DO

            IDNOR=3*NSLPTL
            IDDIR=6*NSLPTL
            DO J=1,NSLPTL
               STATEV(J)=GSLP1(J)
               STATEV(NSLPTL+J)=GAMMA1(J)
               STATEV(2*NSLPTL+J)=TAUSP1(J)
               STATEV(NSTATV-3*(NSLPTL+1)+J)=XKINSLP11(J)
               STATEV(NSTATV-2*(NSLPTL+1)+J)=XKINSLP12(J)
               STATEV(NSTATV-(NSLPTL+1)+J)=XKINSLP13(J)

               DO I=1,3
                  IDNOR=IDNOR+1
                  STATEV(IDNOR)=SPNOR1(I,J)

                  IDDIR=IDDIR+1
                  STATEV(IDDIR)=SPDIR1(I,J)
               END DO
            END DO

         END IF

      END IF

!-----  Total cumulative shear strains on all slip systems (sum of the 
!       absolute values of shear strains in all slip systems)
!--     Total cumulative shear strains on each slip system (sum of the 
!       absolute values of shear strains in each individual slip system)

      DO I=1,NSLPTL

         STATEV(500)=STATEV(500)+ABS(DGAMMA(I))
         STATEV(9*NSLPTL+I)=STATEV(9*NSLPTL+I)+ABS(DGAMMA(I))

      END DO

!----- Effective plastic slip parameter calculated as a sum of slip rates 
!	    all slip systems integrated over time. 

	   PEFF=STATEV(501)

      DPEFF=0.0
      
	   DO I=1,3
      	DO J=1,3
       	PLDT(I,J)=0.0
        END DO
      END DO

	   DO I=1,NSLPTL
		   DO J=1,3
        	   DO K=1,3
        	      PLDT(J,K)=PLDT(J,K)+DGAMMA(I)*SLPDIR(J,I)*SLPNOR(K,I)
            END DO
         END DO
      END DO

	   DO I=1,3
         DO J=1,3
		      DPEFF=DPEFF+PLDT(I,J)*PLDT(I,J)	
         END DO 
	   END DO

	   DPEFF=SQRT(DPEFF*(2.0/3.0))
      STATEV(501)=PEFF+DPEFF

      
      RETURN
      END

!----------------------------------------------------------------------
! Subroutines
      
      include 'utils.for'
      include 'lubksb.for'
      include 'ludcmp.for'
      include 'kinharden.for'
      include 'iteration.for'
      include 'latentharden.for'
      include 'strainrate.for'
      include 'gslpinit.for'
      include 'slipsys.for'
      include 'rotation.for'
      include 'noncubiccrystal.for'