 MODULE stepPotentials
!
! Purpose:
! To assign step potentials (squared wells) between pair of particles
! Main program and energy module call functions in this module to define potential energy
! interaction between particles i and j. This module does not decide which potential applies
! it just contains the collection of functions completing all cases.
!
! Record of revisions:
!      DATE            PROGRAMMER           DESCRIPTION OF CHANGE
!    =======          ============         =======================
!     ??              Agustí, JL             Original Code
!    25/08/11         Pedro                  Comments added
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC FACTE, BOND, SS, COUL, HFIL, HFOB, MAXSTEPS,   &
       stepPot, stepPotInt, getStepSSec, getStepCoul, &
       getStepHFil, getStepHFob
!
! Data dictionary and variable declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)  ! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6)   ! Single Precision real number definition processor independent
! Meaningful
INTEGER, PARAMETER :: MAXSTEPS = 3        ! Max number of steps of each potential
INTEGER, PARAMETER :: BOND=0              ! Code of kind of potential BONDED
INTEGER, PARAMETER :: SS=1                ! Code of kind of potential PSEUDO-BONDED SEC STRUCT
INTEGER, PARAMETER :: COUL=2              ! Code of kind of potential COULOMBIC
INTEGER, PARAMETER :: HFIL=3              ! Code of kind of potential HIDROPHILIC
INTEGER, PARAMETER :: HFOB=4              ! Code of kind of potential HIDROPHOBIC
REAL(SGL), PARAMETER :: FACTE = 4184._SGL ! Energy unit convertion from kcal to Joules (4184!)

TYPE stepPot
! Basic unit of step potential, radial distance from
! core and depth of the discontinuity
! ATT ORDER MATTERS
    REAL(SGL) :: r       ! Radial distance from particle center
                         ! where discontinuity occours
    REAL(SGL) :: e       ! Energy or depth of the well
END TYPE stepPot

TYPE stepPotInt
! Step potential interaction def and properties
    INTEGER nstep       ! Number of steps for this potential interaction
    INTEGER tipInt      ! Tipus of interaction (Bonded, coul, hfob, hfil ,SS)
    TYPE (stepPot) step(MAXSTEPS) ! Squared well components stepPot has
                                  ! MAXSTEPS maximum steps with (distance, energy)
    LOGICAL active      ! Control variable of step status
END TYPE stepPotInt
 

 CONTAINS
!===============================================================================
  FUNCTION getStepSSec (sigmago, ego, dist, active) RESULT (st)
 !
 ! Purpose:
 ! To assign pseudo-bonds interactions between particles in defined secondary structure
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE(stepPotInt) :: st                 ! Result: Step potential interaction
 REAL(SGL), INTENT(IN) :: sigmago       ! 2* sigmago is the amplitude of sq well
 REAL(SGL), INTENT(IN) :: ego           ! Depth (pot energy) of the well
 REAL(SGL), INTENT(IN) :: dist          ! Distance from core particle to sq well center
 LOGICAL active                         ! Control of th status of sqwell
    st%step(1)=stepPot(0._SGL, 0._SGL)
    st%step(2)=stepPot(0._SGL, 0._SGL)
    !
    st%nstep=2
    st%step(1)=stepPot((1.-SIGMAGO)*dist, -EGO*FACTE) ! ORDER MATTERS (r, e)
    st%step(2)=stepPot((1.+SIGMAGO)*dist, EGO*FACTE)  ! ORDER MATTERS (r, e)
    st%active= active
    st%tipInt=SS
 END FUNCTION getStepSSec
!===============================================================================
 FUNCTION getStepCoul (rvdwij, dpsint, dpsext, ecoul, hps, active) RESULT (st)
 !
 ! Purpose:
 ! To assign non-bonded coulombic interaction between two particles
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE(stepPotInt) :: st                  ! Result: Step potential interaction
 REAL(SGL), INTENT(IN) :: rvdwij         ! Sum of vdW radii of i j, hardcore repulsion
 REAL(SGL), INTENT(IN) :: dpsint         ! First energy discontinuity
 REAL(SGL), INTENT(IN) :: dpsext         ! Last energy discontinuity (most external)
 REAL(SGL), INTENT(IN) :: ecoul          ! Max energy of interaction
 REAL(SGL), INTENT(IN) :: hps            ! Shape modulator ( hps < 1) mimic 1/r effect
 LOGICAL active                          ! Control of th status of sqwell
    st%step(1)=stepPot(0._SGL, 0._SGL)
    st%step(2)=stepPot(0._SGL, 0._SGL)
    st%step(3)=stepPot(0._SGL, 0._SGL)
    !
    st%nstep=3
    st%step(1)=stepPot(rvdwij, -sign(3.,ecoul) * ecoul * FACTE) ! ORDER MATTERS (r, e)
    st%step(2)=stepPot(dpsint, -(1.-hps) * ecoul * FACTE) ! ORDER MATTERS (r, e)
    st%step(3)=stepPot(dpsext, -hps * ecoul * FACTE) ! ORDER MATTERS (r, e)
    st%active=active
    st%tipInt=COUL
 END FUNCTION getStepCoul
!===============================================================================
 FUNCTION getStepHFil (rvdwij, dhf, esolv, active) RESULT (st)
 !
 ! Purpose:
 ! To assign non-bonded hidrophilic interaction between two particles
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (stepPotInt) st                   ! Result: Step potential interaction
 REAL(SGL), INTENT(IN) :: rvdwij        ! Sum of vdW radii of i j, hardcore repulsion
 REAL(SGL), INTENT(IN) :: dhf           ! Last energy discontinuity (most external)
 REAL(SGL), INTENT(IN) :: esolv         ! Energy of interaction (depth sqwell)
 LOGICAL active
    st%step(1)=stepPot(0._SGL, 0._SGL)
    st%step(2)=stepPot(0._SGL, 0._SGL)
    !
    st%nstep=2
    st%step(1)=stepPot(rvdwij, -1.5*esolv*FACTE) ! ORDER MATTERS (r, e)
    st%step(2)=stepPot(dhf, -esolv*FACTE) ! ORDER MATTERS (r, e)
    st%active=active
    st%tipInt=HFIL
 END FUNCTION getStepHFil
!===============================================================================
 FUNCTION getStepHFob (rvdwij, dhf, esolv, active) RESULT (st)
 !
 ! Purpose:
 ! To assign non-bonded hidrophobic interaction between two particles
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (stepPotInt) st                   ! Result: Step potential interaction
 REAL(SGL), INTENT(IN) :: rvdwij        ! Sum of vdW radii of i j, hardcore repulsion
 REAL(SGL), INTENT(IN) :: dhf           ! Last energy discontinuity (most external)
 REAL(SGL), INTENT(IN) :: esolv         ! Energy of interaction (depth sqwell)
 LOGICAL active
    st%step(1)=stepPot(0._SGL, 0._SGL)
    st%step(2)=stepPot(0._SGL, 0._SGL)
    !
    st%nstep=2
    st%step(1)=stepPot(rvdwij, 3.0*esolv*FACTE)
    st%step(2)=stepPot(dhf, -esolv*FACTE)
    st%active=active
    st%tipInt=HFOB
 END FUNCTION getStepHFob
!===============================================================================
 END MODULE stepPotentials
