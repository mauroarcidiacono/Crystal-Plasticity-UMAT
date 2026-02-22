      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

! UMAT for Crystal Plasticity Modelling
! Last modified by: Mauro Francisco Arcidiacono

! This UMAT is an update from a previous UMAT originally developed by
! Yonggang Huang with the corrections introduced by J. W. Kysar. It
! was modified by M. F. Arcidiacono to include Modern FORTRAN features,
! kinematic hardening, the HCP crystal structure and to modularize
! the code. Multiple enhacements to improve readability were also
! included.  


!-----  Subroutines:
!
!       ROTATION     -- forming rotation matrix, i.e. the direction 
!                       cosines of cubic crystal [100], [010] and [001]
!                       directions in global system at the initial 
!                       state
!
!       SLIPSYS      -- calculating number of slip systems, unit 
!                       vectors in slip directions and unit normals to 
!                       slip planes in a cubic crystal at the initial 
!                       state
!
!       GSLPINIT     -- calculating initial value of current strengths 
!                       at initial state
!
!       STRAINRATE   -- based on current values of resolved shear 
!                       stresses and current strength, calculating 
!                       shear strain-rates in slip systems
!
!       LATENTHARDEN -- forming self- and latent-hardening matrix 
!
!       ITERATION    -- generating arrays for the Newton-Rhapson 
!                       iteration
!
!       LUDCMP       -- LU decomposition
!
!       LUBKSB       -- linear equation solver based on LU 
!                       decomposition method (must call LUDCMP first)
!
!       KINHARDEN	   -- calculates slip system increment of kinematic
!				hardening based on current kinematic hardening 
!				variable
!						 


!-----  Function subprogram:
!       F -- shear strain-rates in slip systems


!-----  Variables:
!
!       STRESS -- stresses (INPUT & OUTPUT)
!                 Cauchy stresses for finite deformation
!       STATEV -- solution dependent state variables (INPUT & OUTPUT)
!       DDSDDE -- Jacobian matrix (OUTPUT)

!-----  Variables passed in for information:
!
!       STRAN  -- strains
!                 logarithmic strain for finite deformation 
!                 (actually, integral of the symmetric part of velocity
!                  gradient with respect to time)
!       DSTRAN -- increments of strains
!       CMNAME -- name given in the *MATERIAL option
!       NDI    -- number of direct stress components
!       NSHR   -- number of engineering shear stress components
!       NTENS  -- NDI+NSHR
!       NSTATV -- number of solution dependent state variables (as 
!                 defined in the *DEPVAR option)
!       PROPS  -- material constants entered in the *USER MATERIAL 
!                 option
!       NPROPS -- number of material constants
!

!-----  This subroutine provides the plastic constitutive relation of 
!     single crystals for finite element code ABAQUS. The plastic slip
!     of single crystal obeys the Schmid law. The program gives the 
!     choice of small deformation theory and theory of finite rotation 
!     and finite strain.
!       The strain increment is composed of elastic part and plastic 
!     part. The elastic strain increment corresponds to lattice 
!     stretching, the plastic part is the sum over all slip systems of 
!     plastic slip. The shear strain increment for each slip system is
!     assumed a function of the ratio of corresponding resolved shear 
!     stress over current strength, and of the time step. The resolved
!     shear stress is the double product of stress tensor with the slip
!     deformation tensor (Schmid factor), and the increment of current 
!     strength is related to shear strain increments over all slip 
!     systems through self- and latent-hardening functions.

!-----  The implicit integration method proposed by Peirce, Shih and 
!     Needleman (1984) is used here. The subroutine provides an option
!     of iteration to solve stresses and solution dependent state 
!     variables within each increment.

!-----  The present program is for a single CUBIC crystal and a HCP crystal.  
!     However, this code can be generalized for other crystals (e.g. 
!     Tetragonal, Orthotropic, etc.).


!-----  Important notice:
!
!     (1) The number of state variables NSTATV must be larger than (or 
!         equal to) THIRTEEN (13) times the total number of slip systems
!         in all sets, NSLPTL, plus SEVEN (7)

!         NSTATV >= 13 * NSLPTL + 7

!         Denote s as a slip direction and m as normal to a slip plane.
!         Here (s,-m), (-s,m) and (-s,-m) are NOT considered 
!         independent of (s,m).  The number of slip systems in each set
!         could be either 6, 12, 24 or 48 for a cubic crystal, e.g. 12 
!         for {110}<111>.

!         Users who need more parameters to characterize the 
!         constitutive law of single crystal, e.g. the framework 
!         proposed by Zarka, should make NSTATV larger than (or equal 
!         to) the number of those parameters NPARMT plus thirteen times 
!         the total number of slip systems, NSLPTL, plus seven
!         
!         NSTATV >= NPARMT + 13 * NSLPTL + 7

!     (2) The tangent stiffness matrix in general is not symmetric if 
!         latent hardening is considered.  Users must declare "UNSYMM" 
!         in the input file, at the *USER MATERIAL card.


      PARAMETER (ND=150)
!  The parameter ND determines the dimensions of the arrays in 
!  this subroutine.  The current choice 150 is a upper bound for a 
!  cubic crystal with up to three sets of slip systems activated.  
!  Users may reduce the parameter ND to any number as long as larger
!  than or equal to the total number of slip systems in all sets.  
!  For example, if {110}<111> is the only set of slip system 
!  potentially activated, ND could be taken as twelve (12).  

      IMPLICIT REAL*8 (A-H,O-Z)

CHARACTER*80 CMNAME
      EXTERNAL F


      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)


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
     
      REAL*8 LPROPS(184), STRESS2(NTENS), DSTRAN2(NTENS) 
      REAL*8 total_stran(NTENS), max_value, temp_value
      REAL*8 S(6,6)
      INTEGER n_div

      
!-----  NSLIP  -- number of slip systems in each set
!-----  SLPDIR -- slip directions (unit vectors in the initial state)
!-----  SLPNOR -- normals to slip planes (unit normals in the initial 
!                 state)
!-----  SLPDEF -- slip deformation tensors (Schmid factors)
!                 SLPDEF(1,i) -- SLPDIR(1,i)*SLPNOR(1,i)
!                 SLPDEF(2,i) -- SLPDIR(2,i)*SLPNOR(2,i)
!                 SLPDEF(3,i) -- SLPDIR(3,i)*SLPNOR(3,i)
!                 SLPDEF(4,i) -- SLPDIR(1,i)*SLPNOR(2,i)+
!                                SLPDIR(2,i)*SLPNOR(1,i)
!                 SLPDEF(5,i) -- SLPDIR(1,i)*SLPNOR(3,i)+
!                                SLPDIR(3,i)*SLPNOR(1,i)
!                 SLPDEF(6,i) -- SLPDIR(2,i)*SLPNOR(3,i)+
!                                SLPDIR(3,i)*SLPNOR(2,i)
!                 where index i corresponds to the ith slip system
!-----  SLPSPN -- slip spin tensors (only needed for finite rotation)
!                 SLPSPN(1,i) -- [SLPDIR(1,i)*SLPNOR(2,i)-
!                                 SLPDIR(2,i)*SLPNOR(1,i)]/2
!                 SLPSPN(2,i) -- [SLPDIR(3,i)*SLPNOR(1,i)-
!                                 SLPDIR(1,i)*SLPNOR(3,i)]/2
!                 SLPSPN(3,i) -- [SLPDIR(2,i)*SLPNOR(3,i)-
!                                 SLPDIR(3,i)*SLPNOR(2,i)]/2
!                 where index i corresponds to the ith slip system
!-----  DSPDIR -- increments of slip directions
!-----  DSPNOR -- increments of normals to slip planes
!
!-----  DLOCAL -- elastic matrix in local cubic crystal system
!-----  D      -- elastic matrix in global system
!-----  ROTD   -- rotation matrix transforming DLOCAL to D
!
!-----  ROTATE -- rotation matrix, direction cosines of [100], [010] 
!                 and [001] of cubic crystal in global system
!
!-----  FSLIP  -- shear strain-rates in slip systems
!-----  DFDXSP -- derivatives of FSLIP w.r.t x=(TAUSLP-XKINSLP)/GSLIP, where 
!                 TAUSLP is the resolved shear stress, GSLIP is the 
!                 current strength and XKINSLP is the kinematic variable
!
!-----  DDEMSD -- double dot product of the elastic moduli tensor with 
!                 the slip deformation tensor plus, only for finite 
!                 rotation, the dot product of slip spin tensor with 
!                 the stress
!
!-----  H      -- self- and latent-hardening matrix
!                 H(i,i) -- self hardening modulus of the ith slip 
!                           system (no sum over i)
!                 H(i,j) -- latent hardening molulus of the ith slip 
!                           system due to a slip in the jth slip system
!                           (i not equal j)
!
!-----  DDGDDE -- derivatice of the shear strain increments in slip 
!                 systems w.r.t. the increment of strains
!
!-----  DSTRES -- Jaumann increments of stresses, i.e. corotational 
!                 stress-increments formed on axes spinning with the 
!                 material
!-----  DELATS -- strain-increments associated with lattice stretching
!                 DELATS(1) - DELATS(3) -- normal strain increments
!                 DELATS(4) - DELATS(6) -- engineering shear strain 
!                                          increments
!-----  DSPIN  -- spin-increments associated with the material element
!                 DSPIN(1) -- component 12 of the spin tensor
!                 DSPIN(2) -- component 31 of the spin tensor
!                 DSPIN(3) -- component 23 of the spin tensor
!
!-----  DVGRAD -- increments of deformation gradient in the current 
!                 state, i.e. velocity gradient times the increment of 
!                 time
!
!-----  DGAMMA -- increment of shear strains in slip systems
!-----  DTAUSP -- increment of resolved shear stresses in slip systems 
!-----  DGSLIP -- increment of current strengths in slip systems


!-----  Arrays for iteration:
!
!            FSLIP1, STRES1, GAMMA1, TAUSP1, GSLP1 , SPNOR1, SPDIR1, 
!            DDSDE1, DSOLD , DGAMOD, DTAUOD, DGSPOD, DSPNRO, DSPDRO,
!            DHDGDG


!-----  Solution dependent state variable STATEV:
!            Denote the number of total slip systems by NSLPTL, which 
!            will be calculated in this code.
!
!       Array STATEV:
!       1          - NSLPTL    :  current strength in slip systems
!       NSLPTL+1   - 2*NSLPTL  :  shear strain in slip systems
!       2*NSLPTL+1 - 3*NSLPTL  :  resolved shear stress in slip systems
!
!       3*NSLPTL+1 - 6*NSLPTL  :  current components of normals to slip
!                                 planes
!       6*NSLPTL+1 - 9*NSLPTL  :  current components of slip directions
!
!       9*NSLPTL+1 - 10*NSLPTL :  total cumulative shear strain on each 
!                                 slip system (sum of the absolute 
!                                 values of shear strains in each slip 
!                                 system individually)
!
!       10*NSLPTL+2 - NSTATV-3*(NSLPTL+2)-2  : additional parameters users may need 
!                                 		     to characterize the constitutive law 
!                                 		     of a single crystal (if there are 
!                                 		     any).
!
!       500                    : total cumulative shear strain on all 
!                                slip systems (sum of the absolute 
!                                values of shear strains in all slip 
!                                systems)
!
!       501                    :  effective plastic slip parameter
!
!       NSTATV-3*(NSLPTL+2)+1    :  number of slip systems in the 1st set
!       NSTATV-3*(NSLPTL+2)+2    :  number of slip systems in the 2nd set
!       NSTATV-3*(NSLPTL+2)+3    :  number of slip systems in the 3rd set
!
!	  NSTATV-3*(NSLPTL+1)+1 - NSTATV-1 : current kinematic variable in each 
!						       slip system (backstress)
!
!       NSTATV		             :  total number of slip systems in all 
!                                       sets
!		


!-----  Material constants PROPS:
!
!       PROPS(1) - PROPS(21) -- elastic constants for a general elastic
!                               anisotropic material
!
!            isotropic   : PROPS(i)=0  for  i>2
!                          PROPS(1) -- Young's modulus
!                          PROPS(2) -- Poisson's ratio
!
!            cubic       : PROPS(i)=0  for i>3
!                          PROPS(1) -- c11
!                          PROPS(2) -- c12
!                          PROPS(3) -- c44
!
!            orthotropic : PORPS(i)=0  for  i>9
!                          PROPS(1) - PROPS(9) are D1111, D1122, D2222,
!                          D1133, D2233, D3333, D1212, D1313, D2323, 
!                          respectively, which has the same definition 
!                          as ABAQUS for orthotropic materials
!                          (see *ELASTIC card)
!
!            anisotropic : PROPS(1) - PROPS(21) are D1111, D1122, 
!                          D2222, D1133, D2233, D3333, D1112, D2212, 
!                          D3312, D1212, D1113, D2213, D3313, D1213, 
!                          D1313, D1123, D2223, D3323, D1223, D1323, 
!                          D2323, respectively, which has the same 
!                          definition as ABAQUS for anisotropic 
!                          materials (see *ELASTIC card)
!
!
!       PROPS(25) - PROPS(56) -- parameters characterizing all slip 
!                                systems to be activated in a cubic 
!                                crystal
!
!            PROPS(25) -- number of sets of slip systems (maximum 3), 
!                         e.g. (110)[1-11] and (101)[11-1] are in the 
!                         same set of slip systems, (110)[1-11] and 
!                         (121)[1-11] belong to different sets of slip 
!                         systems
!                         (It must be a real number, e.g. 3., not 3 !)
!
!            PROPS(33) - PROPS(35) -- normal to a typical slip plane in
!                                     the first set of slip systems, 
!                                     e.g. (1 1 0)
!                                     (They must be real numbers, e.g. 
!                                      1. 1. 0., not 1 1 0 !)
!            PROPS(36) - PROPS(38) -- a typical slip direction in the 
!                                     first set of slip systems, e.g. 
!                                     [1 1 1]
!                                     (They must be real numbers, e.g. 
!                                      1. 1. 1., not 1 1 1 !)
!
!            PROPS(41) - PROPS(43) -- normal to a typical slip plane in
!                                     the second set of slip systems
!                                     (real numbers)
!            PROPS(44) - PROPS(46) -- a typical slip direction in the 
!                                     second set of slip systems
!                                     (real numbers)
!
!            PROPS(49) - PROPS(51) -- normal to a typical slip plane in
!                                     the third set of slip systems
!                                     (real numbers)
!            PROPS(52) - PROPS(54) -- a typical slip direction in the 
!                                     third set of slip systems
!                                     (real numbers)
!
!
!       PROPS(57) - PROPS(72) -- parameters characterizing the initial 
!                                orientation of a single crystal in 
!                                global system
!            The directions in global system and directions in local 
!            cubic crystal system of two nonparallel vectors are needed
!            to determine the crystal orientation.
!
!            PROPS(57) - PROPS(59) -- [p1 p2 p3], direction of first 
!                                     vector in local cubic crystal 
!                                     system, e.g. [1 1 0]
!                                     (They must be real numbers, e.g. 
!                                      1. 1. 0., not 1 1 0 !)
!            PROPS(60) - PROPS(62) -- [P1 P2 P3], direction of first 
!                                     vector in global system, e.g. 
!                                     [2. 1. 0.]
!                                     (It does not have to be a unit 
!                                      vector)
!
!            PROPS(65) - PROPS(67) -- direction of second vector in 
!                                     local cubic crystal system (real 
!                                     numbers)
!            PROPS(68) - PROPS(70) -- direction of second vector in 
!                                     global system
!
!
!       PROPS(73) - PROPS(96) -- parameters characterizing the visco-
!                                plastic constitutive law (shear 
!                                strain-rate vs. resolved shear 
!                                stress), e.g. a power-law relation
!
!            PROPS(73) - PROPS(80) -- parameters for the first set of 
!                                     slip systems
!            PROPS(81) - PROPS(88) -- parameters for the second set of 
!                                     slip systems
!            PROPS(89) - PROPS(96) -- parameters for the third set of 
!                                     slip systems
!
!
!       PROPS(97) - PROPS(144)-- parameters characterizing the self-
!                                and latent-hardening laws of slip 
!                                systems
!
!            PROPS(97) - PROPS(104)-- self-hardening parameters for the
!                                     first set of slip systems
!            PROPS(105)- PROPS(112)-- latent-hardening parameters for 
!                                     the first set of slip systems and
!                                     interaction with other sets of 
!                                     slip systems
!
!            PROPS(113)- PROPS(120)-- self-hardening parameters for the
!                                     second set of slip systems
!            PROPS(121)- PROPS(128)-- latent-hardening parameters for 
!                                     the second set of slip systems 
!                                     and interaction with other sets 
!                                     of slip systems
!
!            PROPS(129)- PROPS(136)-- self-hardening parameters for the
!                                     third set of slip systems
!            PROPS(137)- PROPS(144)-- latent-hardening parameters for 
!                                     the third set of slip systems and
!                                     interaction with other sets of
!                                     slip systems
!
!
!       PROPS(145)- PROPS(152)-- parameters characterizing forward time
!                                integration scheme and finite 
!                                deformation
!
!            PROPS(145) -- parameter theta controlling the implicit 
!                          integration, which is between 0 and 1
!                          0.  : explicit integration
!                          0.5 : recommended value
!                          1.  : fully implicit integration
!
!            PROPS(146) -- parameter NLGEOM controlling whether the 
!                          effect of finite rotation and finite strain 
!                          of crystal is considered,
!                          0.        : small deformation theory
!                          otherwise : theory of finite rotation and 
!                                      finite strain
!
!
!       PROPS(153)- PROPS(160)-- parameters characterizing iteration 
!                                method
!
!            PROPS(153) -- parameter ITRATN controlling whether the 
!                          iteration method is used, 
!                          0.        : no iteration
!                          otherwise : iteration
!
!            PROPS(154) -- maximum number of iteration ITRMAX 
!
!            PROPS(155) -- absolute error of shear strains in slip 
!                          systems GAMERR
!
!
!       PROPS(161)- PROPS(184)-- parameters characterizing kinematic 
!                                hardening
!
!            PROPS(161) -- parameter CKIN1 
!            PROPS(162) -- parameter DKIN1
!
!            PROPS(169) -- parameter CKIN2 
!            PROPS(170) -- parameter DKIN2
!
!            PROPS(177) -- parameter CKIN3 
!            PROPS(178) -- parameter DKIN3

      CALL CORE(STRESS,STATEV,DDSDDE,
     1 STRAN,DSTRAN,DTIME,ND,NOEL,
     2 NDI,NSHR,NTENS,NSTATV,PROPS,DROT)
      
      RETURN
      END

! Subroutines
      
      include 'core.for'
