      SUBROUTINE LATENTHARDEN (GAMMA, TAUSLP, GSLIP, GMSLTL, GAMTOL, 
     2                         NSLIP, NSLPTL, NSET, H, PROP, ND)

!-----  This subroutine calculates the current self- and latent-
!     hardening moduli for all slip systems in a rate-dependent single 
!     crystal.  Two kinds of hardening law are used here.  The first 
!     law, proposed by Asaro, and Pierce et al, assumes a HYPER SECANT 
!     relation between self- and latent-hardening moduli and overall 
!     shear strain.  The Bauschinger effect has been neglected.  The 
!     second is Bassani's hardening law, which gives an explicit 
!     expression of slip interactions between slip systems.  The 
!     classical three stage hardening for FCC single crystal could be 
!     simulated.

!-----  The hardening coefficients are assumed the same for all slip 
!     systems in each set, though they could be different from set to 
!     set, e.g. <110>{111} and <110>{100}.

!-----  Users who want to use their own self- and latent-hardening law 
!     may change the function subprograms HSELF (self hardening) and 
!     HLATNT (latent hardening).  The parameters characterizing these 
!     hardening laws are passed into HSELF and HLATNT through array 
!     PROP.


!-----  Function subprograms:
!
!       HSELF  -- User-supplied self-hardening function in a slip 
!                 system
!
!       HLATNT -- User-supplied latent-hardening function

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
!     NSLIP  -- number of slip systems in each set (INPUT)
!     NSLPTL -- total number of slip systems in all the sets (INPUT)
!     NSET   -- number of sets of slip systems (INPUT)
!
!     H      -- current value of self- and latent-hardening moduli 
!               (OUTPUT)
!               H(i,i) -- self-hardening modulus of the ith slip system
!                         (no sum over i)
!               H(i,j) -- latent-hardening molulus of the ith slip 
!                         system due to a slip in the jth slip system 
!                         (i not equal j)
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
!     ND     -- leading dimension of arrays defined in subroutine UMAT 
!               (INPUT) 


!-----  Use single precision on cray
!
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL HSELF, HLATNT

      DIMENSION GAMMA(NSLPTL), TAUSLP(NSLPTL), GMSLTL(NSLPTL),
     2          GSLIP(NSLPTL), NSLIP(NSET), PROP(16,NSET), 
     3          H(ND,NSLPTL)


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

            DO LATENT=1,NSLPTL
               IF (LATENT.EQ.ISELF) THEN

                  H(LATENT,ISELF)=HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,
     2                                  NSET,NSLIP,PROP(1,I),CHECK,
     3                                  ISELF,ISET)

               ELSE

                  H(LATENT,ISELF)=HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,
     2                                   NSET,NSLIP,PROP(1,I),CHECK,
     3                                   ISELF,ISET,LATENT)


               END IF
            END DO

         END DO
      END DO

      RETURN
      END


!-----------------------------------


!-----  Use single precision on cray

           REAL*8 FUNCTION HSELF(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                           NSLIP,PROP,CHECK,ISELF,ISET)


!-----     User-supplied self-hardening function in a slip system

!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)

           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
     2               GMSLTL(NSLPTL)


           IF (CHECK.EQ.0.) THEN

!-----  HYPER SECANT hardening law by Asaro, Pierce et al
              TERM1=PROP(1)*GAMTOL/(PROP(2)-PROP(3))
              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              HSELF=PROP(1)*TERM2**2

           ELSE

!-----  Bassani's hardening law

              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))

              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

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
                    IF (ID.NE.ISELF) THEN

		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)

		    END IF

                 END DO
              END DO

              HSELF=F*G

           END IF

           RETURN
           END


!-----------------------------------


!-----  Use single precision on cray

           REAL*8 FUNCTION HLATNT(GAMMA,GMSLTL,GAMTOL,NSLPTL,NSET,
     2                            NSLIP,PROP,CHECK,ISELF,ISET,LATENT)


!-----     User-supplied latent-hardening function

!-----  Use single precision on cray
!
           IMPLICIT REAL*8 (A-H,O-Z)

           DIMENSION GAMMA(NSLPTL), NSLIP(NSET), PROP(16),
     2               GMSLTL(NSLPTL)


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
              HLATNT=PROP(1)*TERM2**2*Q

           ELSE

!-----  Bassani's hardening law

              TERM1=(PROP(1)-PROP(4))*GMSLTL(ISELF)/(PROP(2)-PROP(3))

              TERM2=2.*EXP(-TERM1)/(1.+EXP(-2.*TERM1))
              F=(PROP(1)-PROP(4))*TERM2**2+PROP(4)

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
                    IF (ID.NE.ISELF) THEN

		       G=G+FAB*TANH(GMSLTL(ID)/GAMMA0)

		    END IF

                 END DO
              END DO

              HLATNT=F*G*Q

           END IF

           RETURN
           END