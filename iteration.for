      SUBROUTINE ITERATION (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, 
     2                      NSLPTL, NSET, NSLIP, ND, PROP, DGAMOD, 
     3                      DHDGDG)

!-----  This subroutine generates arrays for the Newton-Rhapson 
!     iteration method.

!-----  Users who want to use their own self- and latent-hardening law 
!     may change the function subprograms DHSELF (self hardening) and 
!     DHLATN (latent hardening).  The parameters characterizing these 
!     hardening laws are passed into DHSELF and DHLATN through array 
!     PROP.


!-----  Function subprograms:
!
!       DHSELF -- User-supplied function of the derivative of self-
!                 hardening moduli
!
!       DHLATN -- User-supplied function of the derivative of latent-
!                 hardening moduli

!-----  Variables:
!
!     GAMMA  -- shear strain in all slip systems at the start of time 
!               step  (INPUT)
!     TAUSLP -- resolved shear stress in all slip systems (INPUT)
!     GSLIP  -- current strength (INPUT)
!     GMSLTL -- total cumulative shear strains on each individual slip system 
!               (INPUT)
!     GAMTOL -- total cumulative shear strains over all slip systems 
!               (INPUT)
!     NSLPTL -- total number of slip systems in all the sets (INPUT)
!     NSET   -- number of sets of slip systems (INPUT)
!     NSLIP  -- number of slip systems in each set (INPUT)
!     ND     -- leading dimension of arrays defined in subroutine UMAT 
!               (INPUT) 
!
!     PROP   -- material constants characterizing the self- and latent-
!               hardening law (INPUT)
!
!               For the HYPER SECANT hardening law 
!               PROP(1,i) -- initial hardening modulus H0 in the ith 
!                            set of slip systems
!               PROP(2,i) -- saturation stress TAUs in the ith set of  
!                            slip systems
!               PROP(3,i) -- initial critical resolved shear stress 
!                            TAU0 in the ith set of slip systems
!               PROP(9,i) -- ratio of latent to self-hardening Q in the
!                            ith set of slip systems
!               PROP(10,i)-- ratio of latent-hardening from other sets 
!                            of slip systems to self-hardening in the 
!                            ith set of slip systems Q1
!
!               For Bassani's hardening law 
!               PROP(1,i) -- initial hardening modulus H0 in the ith 
!                            set of slip systems
!               PROP(2,i) -- stage I stress TAUI in the ith set of  
!                            slip systems (or the breakthrough stress 
!                            where large plastic flow initiates)
!               PROP(3,i) -- initial critical resolved shear stress 
!                            TAU0 in the ith set of slip systems
!               PROP(4,i) -- hardening modulus during easy glide Hs in 
!                            the ith set of slip systems
!               PROP(5,i) -- amount of slip Gamma0 after which a given 
!                            interaction between slip systems in the 
!                            ith set reaches peak strength
!               PROP(6,i) -- amount of slip Gamma0 after which a given 
!                            interaction between slip systems in the 
!                            ith set and jth set (i not equal j) 
!                            reaches peak strength
!               PROP(7,i) -- representing the magnitude of the strength
!                            of interaction in the ith set of slip 
!                            system
!               PROP(8,i) -- representing the magnitude of the strength
!                            of interaction between the ith set and jth
!                            set of system
!               PROP(9,i) -- ratio of latent to self-hardening Q in the
!                            ith set of slip systems
!               PROP(10,i)-- ratio of latent-hardening from other sets 
!                            of slip systems to self-hardening in the 
!                            ith set of slip systems Q1
!
!-----  Arrays for iteration:
!
!       DGAMOD (INPUT)
!
!       DHDGDG (OUTPUT)
!

!-----  Use single precision on cray
!
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL DHSELF, DHLATN

      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET), 
     3          DGAMOD(NSLPTL), DHDGDG(ND,NSLPTL)


      CHECK=0.
      DO I=1,NSET
         DO J=4,8
            CHECK=CHECK+ABS(PROP(J,I))
         END DO
      END DO

!-----  CHECK=0   --  HYPER SECANT hardening law
!       otherwise --  Bassani's hardening law

      ISELF=0
      DO I=1,NSET
         ISET=I
         DO J=1,NSLIP(I)
            ISELF=ISELF+1

            DO KDERIV=1,NSLPTL
               DHDGDG(ISELF,KDERIV)=0.

               DO LATENT=1,NSLPTL
                  IF (LATENT.EQ.ISELF) THEN

                     DHDG=DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
     3                           KDERIV)

                  ELSE

                     DHDG=DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP(1,I),CHECK,ISELF,ISET,
     3                           LATENT,KDERIV)

                  END IF

                  DHDGDG(ISELF,KDERIV)=DHDGDG(ISELF,KDERIV)+
     2                                 DHDG*ABS(DGAMOD(LATENT))
               END DO

            END DO
         END DO
      END DO

      RETURN
      END


!-----------------------------------


!-----  Use single precision on cray

           REAL*8 FUNCTION DHSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,
     3                            KDERIV)


!-----  User-supplied function of the derivative of self-hardening
!     moduli

!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)

           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), 
     2               NSLIP(NSET), PROP(16)


           IF (CHECK.EQ.0.) THEN

!-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
              DHSELF=-2.*PROP(1)*TERM2**2*TANH(TERM1)*TERM3

           ELSE

!-----  Bassani's hardening law

              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

              IF (KDERIV.EQ.ISELF) THEN
                 F=-2.*(PROP(1)-PROP(4))*TERM2**2*TANH(TERM1)*TERM3
                 ID=0
                 G=1.
                 DO I=1,NSET
                    IF (I.EQ.ISET) THEN
                       GAMMA0=PROP(5)
                       FAB=PROP(7)
                    ELSE
                       GAMMA0=PROP(6)
                       FAB=PROP(8)
                    END IF

                    DO J=1,NSLIP(I)
                       ID=ID+1

                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)

                    END DO
                 END DO

              ELSE
                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
                 ILOWER=0
                 IUPPER=NSLIP(1)
                 IF (ISET.GT.1) THEN
                    DO K=2,ISET
                       ILOWER=ILOWER+NSLIP(K-1)
                       IUPPER=IUPPER+NSLIP(K)
                    END DO
                 END IF

                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

                 TERM4=GMSLTL(KDERIV)/GAMMA0

                 TERM5=2.*EXP(-TERM4)/(1.+EXP(-2.*TERM4))
                 G=FAB/GAMMA0*TERM5**2

              END IF

              DHSELF=F*G

           END IF

           RETURN
           END


!-----------------------------------


!-----  Use single precision on cray

           REAL*8 FUNCTION DHLATN(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT,
     3                            KDERIV)


!-----  User-supplied function of the derivative of latent-hardening 
!     moduli

!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)

           DIMENSION GAMMA(NSLPTL), GMSLTL(NSLPTL), NSLIP(NSET), 
     2               PROP(16)

           ILOWER=0
           IUPPER=NSLIP(1)
           IF (ISET.GT.1) THEN
              DO K=2,ISET
                 ILOWER=ILOWER+NSLIP(K-1)
                 IUPPER=IUPPER+NSLIP(K)
              END DO
           END IF

           IF (LATENT.GT.ILOWER.AND.LATENT.LE.IUPPER) THEN
              Q=PROP(9)
           ELSE
              Q=PROP(10)
           END IF

           IF (CHECK.EQ.0.) THEN

!-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=PROP(1)/(PROP(2)-PROP(3))*DSIGN(1.D0,GAMMA(KDERIV))
              DHLATN=-2.*PROP(1)*TERM2**2*TANH(TERM1)*TERM3*Q

           ELSE

!-----  Bassani's hardening law

              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              TERM3=(PROP(1)-PROP(4))/(PROP(2)-PROP(3))

              IF (KDERIV.EQ.ISELF) THEN
                 F=-2.*(PROP(1)-PROP(4))*TERM2**2*TANH(TERM1)*TERM3
                 ID=0
                 G=1.
                 DO I=1,NSET
                    IF (I.EQ.ISET) THEN
                       GAMMA0=PROP(5)
                       FAB=PROP(7)
                    ELSE
                       GAMMA0=PROP(6)
                       FAB=PROP(8)
                    END IF

                    DO J=1,NSLIP(I)
                       ID=ID+1

                       IF (ID.NE.ISELF) G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)

                    END DO
                 END DO

              ELSE
                 F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)
                 ILOWER=0
                 IUPPER=NSLIP(1)
                 IF (ISET.GT.1) THEN
                    DO K=2,ISET
                       ILOWER=ILOWER+NSLIP(K-1)
                       IUPPER=IUPPER+NSLIP(K)
                    END DO
                 END IF

                 IF (KDERIV.GT.ILOWER.AND.KDERIV.LE.IUPPER) THEN
                    GAMMA0=PROP(5)
                    FAB=PROP(7)
                 ELSE
                    GAMMA0=PROP(6)
                    FAB=PROP(8)
                 END IF

                 TERM4=GMSLTL(KDERIV)/GAMMA0
                 TERM5=2.*EXP(-TERM4)/(1.+EXP(-2.*TERM4))
                 G=FAB/GAMMA0*TERM5**2

              END IF

              DHLATN=F*G*Q

           END IF

           RETURN
           END