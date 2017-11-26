 MODULE paramSet
!
! Purpose:
! To define initial parameters, check that are valid ones, and tell
! the user what he/she have introduced as input.
! Parameters are initialized with default values, an they can be modified
! trough a namelist called input.
! Namelist is read after initialization allowing the tunning of parameters
! in a external text file.
!
! Record of revisions:
!      DATE            PROGRAMMER           DESCRIPTION OF CHANGE
!    =======          ============         =======================
!     ??              Agustí, JL             Original Code
!    24/08/11         Pedro                  Comments added
!
IMPLICIT NONE
!
! Data dictionary and variable declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6) ! Single Precision real number definition processor independent
! Physical Meaning
REAL(SGL), SAVE :: fvdw=1.2_SGL      ! Factor multiplying vdW energy term NO UNITS
REAL(SGL), SAVE :: fsolv=1.0_SGL     ! Factor multiplying solvation energy term  NO UNITS
REAL(SGL), SAVE :: eps=1.1_SGL       ! Factor multiplying electrostatic energy term NO UNITS
REAL(SGL), SAVE :: dpsext=7.0_SGL    ! Last sq-well discontinuity for electrostatics ARMSTRONGS
REAL(SGL), SAVE :: dpsint=5.0_SGL    ! Middle sq-well discontinuity for electrostatics ARMSTRONGS
REAL(SGL), SAVE :: dhf=5.0_SGL       ! Last sq-well discontinuity for vdW energy ARMSTRONGS
REAL(SGL), SAVE :: dcut=2.0_SGL      ! Security distance beyond last sq well to still
                                     ! compute possible interactions for a particle being in the
                                     ! center of a sphere radius= last sq well wall + dcut ARMSTRONGS
REAL(SGL), SAVE :: sigma=0.01_SGL    ! Half of the sq well amplidud -1 percent- (bonds)  ARMSTRONG
REAL(SGL), SAVE :: temp=300._SGL     ! Temperature of the system, KELVINS
REAL(SGL), SAVE :: rcutgo=10.0_SGL   ! Maximum distance to compute GO interactions (? CA models) ARMSTRONGS
REAL(SGL), SAVE :: sigmago=0.05_SGL   ! GO squared wells width allowance (?)  ARMSTRONGS
REAL(DBL), SAVE :: tmin=1.0E-22_DBL  ! Minimum time left to event below this time collision is not considered  SECONDS
REAL(SGL), SAVE :: ego=1.0_SGL       ! GO squared well depth (?ENERGY UNITS)
REAL(SGL), SAVE :: hps=0.4_SGL       ! Electrostatic sqwell modifier, it reduces depth of the outtermost
                                     ! part of sphere of interaction ( to mimic 1/r behaviour) NO UNITS
REAL(SGL), SAVE :: rnomax=4.1_SGL    ! Maximum distance between N and O in amide group (pseudo-bonds) ARMSTRONGS
REAL(SGL), SAVE :: rnomin=2.5_SGL    ! Minimum distance between N and O in amide group (pseudo-bonds) ARMSTRONGS
REAL(SGL), SAVE :: rcomax=5.0_SGL    ! Maximum distance between CA and O in amide group (pseudo-bonds) ARMSTRONGS
REAL(SGL), SAVE :: rcomin=3.2_SGL    ! Minimum distance between CA and O in amide group (pseudo-bonds) ARMSTRONGS
REAL(SGL), SAVE :: rcutcoul2=0._SGL  ! CutOff radius of coulomb terms squared ARMSTRONS^2
REAL(SGL), SAVE :: rcutsolv2=0_SGL   ! CutOff radius of solvation terms squared ARMSTRONS^2
LOGICAL, SAVE :: isolv=.TRUE.        ! LOGICAL                    1: Solvation implemented
                                     !                            0: Solvation Ignored NO UNITS
INTEGER, SAVE :: nbloc=1000          ! Number of time (tsnap) units the simulation should go for. NO UNITS
LOGICAL, SAVE :: idab=.TRUE.         ! LOGICAL: 1=TRUE 0=FALSE, To activate diedral angles restrictions between
                                     ! elements of defined secondary structure (i , i+1 )
LOGICAL, SAVE :: igoab=.TRUE.        ! LOGICAL: 1=TRUE 0=FALSE, To activate aditional restrictions in secondary
                                     ! structure elements , go potentials between i and i+4
LOGICAL, SAVE :: iwr=.TRUE.         ! LOGICAL: TRUE To activate energy writing
INTEGER, SAVE :: seed=2381           ! Seed of the number generator routine NO UNITS
INTEGER, SAVE :: tcalc=1             ! Simulation conditions NO UNITS
                                     !             E Const 0,  T Const 1, TBany 2
REAL(DBL), SAVE :: tsnap=1000._DBL   ! Time between each frame displayed. FEMTOSECONDS
REAL(DBL), SAVE :: tcorr=100._DBL    ! Time between each actualization of list
                                     ! possible events, neighbours ??     FEMTOSECONDS
REAL(DBL), SAVE :: tact=50._DBL      ! Time between correction of bonded distances
                                     ! ? SHAKE algorithm                  FEMTOSECONDS
! DIMS
REAL(DBL), SAVE :: trect=100._DBL    ! Time of rectification DIMS
INTEGER, SAVE :: idims=1             !
REAL(SGL), SAVE :: xbeta=0.5_SGL     !
REAL(SGL), SAVE :: sclim=0._SGL      !
REAL(SGL), SAVE :: RMSDlim=2.10_SGL  !
REAL(SGL), SAVE :: acceptance=0.30   !
REAL(SGL), SAVE :: acc_tolerance=0.05!
REAL(SGL), SAVE :: lower_acceptance=0.15!
REAL(SGL), SAVE :: rmsd_cutoff=3.0 !
REAL(SGL), SAVE :: xbeta_update = 5000 ! femtoseconds
REAL(SGL), SAVE :: Drmsd=0.02 !
REAL(SGL),SAVE :: MAX_TIME=1000 ! Maximum simulation times, picoseconds
REAL(SGL), SAVE :: NMCutoff=0.60
INTEGER, SAVE :: nevecs=10
REAL(DBL), SAVE :: timeMC_NM = 20.0 ! Picoseconds
REAL(SGL), SAVE :: time_in_eq = 30.0 ! Picoseconds
REAL(DBL), SAVE :: RMSD_STOP=1E10 ! Armstrongs
REAL(SGL), SAVE :: strict=1    !
INTEGER, SAVE :: NM_Rejections=6 !
LOGICAL, SAVE :: Superpose=.TRUE.
!REAL(SGL), SAVE :: OverlapCut=0.6 !

CONTAINS

 SUBROUTINE readInputParamSet (unit_i)
   !
   ! Purpose:
   ! To read the input parameter set from user. Reading froma a name list in
   ! external file unit_i
   !
   IMPLICIT NONE
   ! Data dictionary and variable declaration
   INTEGER, INTENT(IN) :: unit_i          ! External unit where parameters should be read
   INTEGER  :: istat                      ! Error handling auxiliar
!   
   NAMELIST /input/                                     &
                tsnap , sigma , temp,    &
                seed , rcutgo , sigmago , tmin, &
                dhf , dcut , ego , isolv , fsolv , fvdw &
                , eps , hps , dpsint , dpsext , idab ,  &
                igoab , iwr , rnomax , rnomin , rcomax  &
                , rcomin,acceptance,MAX_TIME,TCALC,&
                nevecs,timeMC_NM,time_in_eq,RMSD_STOP, &
                strict,NM_Rejections,Superpose

! Reading the name list
     READ(unit_i, INPUT, IOSTAT=istat)
     nlistcheck: IF ( istat /= 0) THEN
       ! Open failed stop
       WRITE(0,101)  istat
       101 FORMAT ('Reading Namelist failed IOSTAT =',I6)
       STOP
     ELSE
     ! reading was OK
     ! checking for proper time units
       trectcheck: IF (trect > tsnap) THEN
         WRITE(0,*) " trect cannot be larger than tsnap. Setting trect=tsnap"
         trect = tsnap
         END IF trectcheck
       tcorrcheck: IF (tcorr > tsnap ) THEN
       ! This makes no sense adjusting times
         WRITE(0,*) " tcorr cannot be larger than tsnap. Setting tcorr=tsnap"
         tcorr = tsnap
       END IF tcorrcheck
       tactcheck: IF (tact  > tcorr )  THEN
         ! This is a problem, adjusting times
         WRITE(0,*) " tact cannot be larger than tcorr. Setting tact=tcorr"
         tact = tcorr
       END IF tactcheck
       ! Setting rcut for coul and solv squared
       ! useful for?
       rcutcoul2 = (dpsext + dcut)**2
       rcutsolv2 = (dhf + dcut)**2

     END IF nlistcheck

 END SUBROUTINE readInputParamSet

 SUBROUTINE writeInputParamSet (unit_o)
 !
 ! Purpose:
 ! To write all the parameters introduced for the simulation.
 ! Notify user what the is the simulation about.
 !
 IMPLICIT NONE
 ! Data dictrionary and variable declaration
 INTEGER, INTENT(IN) :: unit_o           ! External unit where to print out parameters
   ! Copiem l'arxiu de parametres. Pendent format
   write (unit_o, *)
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | CALCULATION PARAMETERS                                   |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation settings                                      |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation Time (fs) (Nbloc x TSnap)       |",f12.3," |")') NBLOC * TSNAP / 1.e3
   write (unit_o, '(" | Output structure (fs)             | TSnap  |",f12.3," |")') TSNAP 
   write (unit_o, '(" | Update velocities (fs)            | Tcorr  |",f12.3," |")') TCORR
   write (unit_o, '(" | Update Lists, collision times (fs)| Tact   |",f12.3," |")') TACT
   write (unit_o, '(" | Min. accepted colision time (fs)  | TMin   |",f12.8," |")') TMIN*1.e15   
   write (unit_o, '(" | Temperature (K)                   | Temp   |",6X,f6.2, " |")') TEMP
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Well Potential Definitions                               |")')  
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Cov. Bond relative well width     | Sigma  |",7X,f5.2, " |")') SIGMA
   write (unit_o, '(" | N-O (min-max) (Angs)              | rno    |",f5.2,X,"-",f5.2, " |")') rnomin, rnomax
   write (unit_o, '(" | C-O (min-max) (Angs)              | rco    |",f5.2,X,"-",f5.2, " |")') rcomin, rcomin
   write (unit_o, '(" | CutOff Beta restrains (A)         | RCutGo |",7X,f5.2, " |")') RCUTGO
   write (unit_o, '(" | SSec relative well width          | SigmaGo|",7X,f5.2, " |")') SIGMAGO
   write (unit_o, '(" | SSec wells depth                  | EGo    |",7X,f5.2, " |")') EGO
   write (unit_o, '(" | Coulombic Int. wall (A)           | DPSInt |",7X,f5.2, " |")') DPSINT
   write (unit_o, '(" | Coulombic Ext. wall (A)           | DPSExt |",7X,f5.2, " |")') DPSEXT
   write (unit_o, '(" | 1/Dielectric                      | Eps    |",7X,f5.2, " |")') EPS
   write (unit_o, '(" | E. Fraction on 2nd Coul. well     | HPS    |",7X,f5.2, " |")') HPS   
   write (unit_o, '(" | Solvation Ext. wall (A)           | DHF    |",7X,f5.2, " |")') DHF   
   write (unit_o, '(" | Added Cutoff (A)                  | DCut   |",7X,f5.2, " |")') DCUT
   write (unit_o, '(" | Electrostatic cutoff (A) (DCUT + DPSExt)   |",7X,f5.2, " |")') sqrt(rcutcoul2)
   write (unit_o, '(" | Solvation cutoff (A)     (DCUT + DHF)      |",7X,f5.2, " |")') sqrt(rcutsolv2)
   write (unit_o, '(" | FSolv                                      |",7X,f5.2, " |")') FSOLV
   write (unit_o, '(" | FVdW                                       |",7X,f5.2, " |")') FVDW
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Other                                                    |")')  
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Random generator seed                      |",7X,i5  " |")') seed
   write (unit_o, '(" | IDAB, IGOAB, IWR, ISOLV                    |"L1,1X,L2,2X,2L2" |")') IDAB, IGOAB, IWR, ISOLV
   write (unit_o, '(" ------------------------------------------------------------")')
 END SUBROUTINE writeInputParamSet

 END MODULE paramSet
