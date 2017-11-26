MODULE energy
!
! Purpose:
! Set of tools regarding energy and sq-wells potential interactions.
!
! Record of revisions:
!      DATE            PROGRAMMER           DESCRIPTION OF CHANGE
!    =======          ============         =======================
!     ??              Agustí, JL             Original Code
!    08/09/11         Pedro                  Defined as module. Comments added.
!    21/09/11         Pedro                  Distribution of velocities centered at 0.
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC activateStepPot, thermalize, calcEpot, calcEkin!, &
!       MCcheck, distance,MCWeigth,readPDBtarg,readPDBtarg_CA, &
!       distanceCA
!
! Data dictionary and variable declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)  ! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6)   ! Single Precision real number definition processor independent
REAL, PARAMETER :: PI=3.141592654 ! PI number definition
!
CONTAINS
!===============================================================================
 SUBROUTINE activateStepPot(stepPts, r, rcutcoul2, rcutsolv2, natom, nblist, xsum)
 !
 ! Purpose:
 ! This subroutine checks for possibles interactions between particles.
 ! Finds possible interactions in local enviroment, and keep data of it.
 !
 USE stepPotentials
 USE geometryDP
 USE intList
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 INTEGER, INTENT(IN) :: natom           ! Total number of particles in the system
 REAL(SGL), INTENT(IN) :: rcutcoul2     ! Distance (squared) limit to consider Coulombic
                                        ! interactions (dpsext + dcut)**2
 REAL(SGL), INTENT(IN) :: rcutsolv2     ! Distance (squared) limit to consider Coulombic
                                        ! interactions (dhf + dcut)**2
 TYPE(pointDP), INTENT(IN) :: r(natom)  ! Coordinates of all particles
 TYPE(stepPotInt), INTENT(INOUT) :: stepPts(natom,natom) ! Squared-well potential interactions
                                                         ! between all particles
 TYPE(intpList), INTENT(INOUT) :: nblist(natom) ! Non-bonded list of all atoms
 REAL(SGL), INTENT(IN) :: xsum(natom,natom) ! Function of masses of i and j
                                            ! = ( mi + mj ) / ( 2 mi mj )
 ! Auxiliar block
 INTEGER :: i            ! Loop control index. i is a particle
 INTEGER :: j            ! Loop control index. i is a particle
 REAL(SGL) :: rij2       ! Squared distance between i and j.
! Make all sq-well potentials inactive
 WHERE (stepPts%active)
   stepPts%active = .FALSE.
 END WHERE
! Initialize local-enviroment counter
 nblist%nats = 0
! Iterate only upper diagonal of stepPts(i,j) because
! matrix is symmetric
 DO j = 2,natom
 DO i = 1,j-1
   rij2 = 0._SGL
   IF (stepPts(i,j)%tipInt > SS) rij2 = calcDist2DP(r(i), r(j)) ! SS is always active
   lookforsteppot: IF (stepPts(i,j)%tipInt == SS.OR.&
       (rij2 < rcutcoul2.AND.stepPts(i,j)%tipInt == COUL).OR. &
       (rij2 < rcutsolv2.AND.stepPts(i,j)%tipInt > COUL)) then
!  SS (always) , COUL and HFIL+HFOB are activated so keep intreaction data
       stepPts(i,j)%active = .TRUE.
! Keep track of the number of particles interacting in local enviroment.
       nblist(i)%nats = nblist(i)%nats + 1
       nblist(j)%nats = nblist(j)%nats + 1
! For this interaction found save interaction in intData derived data type structure
       nblist(i)%iData(nblist(i)%nats) = intData(j,nblist(j)%nats, stepPts(i,j),xsum(i,j),1.,0.)
       nblist(j)%iData(nblist(j)%nats) = intData(i,nblist(i)%nats, stepPts(i,j),xsum(i,j),1.,0.)
    END IF lookforsteppot
 END DO
 END DO
END SUBROUTINE activateStepPot
!===============================================================================
 SUBROUTINE thermalize(seed, iev, natom, TEMP, xmassa, v, xm, ekin)
 !
 ! Purpose:
 ! To assign random velocities to all particles. Impose that CM is not moving
 ! and rescale velocities to expected kynetic energy.
 !
 USE geometry
 USE random
 USE geometryDP
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 INTEGER, INTENT(IN) :: iev             !
 INTEGER, INTENT(IN) :: seed            ! Pseudo-random number generator seed
 INTEGER, INTENT(IN) :: natom           ! Total number of particles in the system
 REAL(SGL), INTENT(IN) :: xmassa        ! Total mass of the system
 REAL(SGL), INTENT(IN) :: temp          ! Temperature of the system
 REAL(SGL), INTENT(IN) :: xm(natom)     ! List of masses of all particles
 TYPE(point), INTENT(INOUT) :: v(natom) ! Velocities of all particles
 REAL(SGL), INTENT(OUT) :: ekin         ! Kynetic energy
 ! Auxiliar Block
 INTEGER :: i                           ! Loop control index
 INTEGER :: kk                          ! Random number generator ??
 REAL(SGL) :: sto                       ! ??
 !real calcEkin
 TYPE(point) :: vcm                     ! Velocity of CM
 REAL(DBL) :: fi                        ! Random number generator ??
!
! Initialize velocity of CM
   vcm = point(0., 0., 0.)
! Assing velocity v(x,y,z) of all atoms
  randomvel: DO i = 1,natom ! -1
    kk = seed + 3 * i + 1 + iev
    CALL ran1(fi, kk)
    v(i)%x = fi -DBLE(0.5_DBL)
    kk = seed + 3 * i + 2 + iev
    CALL ran1(fi, kk)
    v(i)%y = fi -DBLE(0.5_DBL)
    kk = seed + 3 * i + 3 + iev
    CALL ran1(fi, kk)
    v(i)%z = fi -DBLE(0.5_DBL)
! This should not be needed
! Compute velocity of center of mass
    vcm = vcm + xm(i) * v(i)
   END DO randomvel
! vcm = (1./(xmassa- xm(natom) )) * vcm
! Ensure vcm=0
!   v(natom)%x= -vcm%x*(xmassa-xm(natom))/xm(natom)
!   v(natom)%y= -vcm%y*(xmassa-xm(natom))/xm(natom)
!   v(natom)%z= -vcm%z*(xmassa-xm(natom))/xm(natom)
!
! vcm=point(0.,0.,0.)
! do i=1,natom!-1
! 	vcm = vcm + xm(i) * v(i)
! enddo
  vcm = (1./xmassa) * vcm
! write(0,*) v(natom)%x,v(natom)%y,v(natom)%z
! This should not be needed any longer
!   vcm = (1./xmassa) * vcm
! Atencion con este bloque!!
   DO i = 1,natom
      v(i) = v(i) - vcm
   END DO
!  vcm=point(0.,0.,0.)
!  do i=1,natom
! 	vcm = vcm + xm(i) * v(i)
 ! enddo
 ! vcm = (1./xmassa) * vcm

!
! write(0,*)"CM", vcm

 ! ??
   sto = sqrt(natom * TEMP / calcEkin(v, xm, natom))
   DO i = 1,natom
      v(i) =  sto * v(i)
   END DO

   ekin = natom * TEMP ! No cal recalcularla 
 END SUBROUTINE thermalize
!===============================================================================
PURE SUBROUTINE calcEpot (natom, r, stepPts, ego, epotgo, epotfis)
! PEDRO: pure attribute added
!
! Purpose:
! To compute potential energy of the system. Epotgo and epotfis are
! "different" kinds of potential energy. The first one without physical
! meaning just keep secondary structures. The second one represent the
! expected physico-chemical potential energy.
!
 USE geometryDP
 USE stepPotentials
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 INTEGER, INTENT(IN) :: natom           ! Total number of particles in the system
 TYPE(pointDP), INTENT(IN) :: r(natom)  ! Coordinates of all particles
 TYPE(stepPotInt), INTENT(IN) :: stepPts(natom,natom) ! Step-Potential interaction
                                                      ! between all particles
 REAL(SGL), INTENT(IN) :: EGO           ! Depth of GO sq-well
 REAL(SGL), intent(OUT) :: epotgo       ! Total GO potential energy of the system
 REAL(SGL), intent(OUT) :: epotfis      ! Total physical potential energy of the system
 ! Auxiliar block
 REAL(SGL) :: dist            ! Distance between i and j
 INTEGER :: i                 ! Loop control index. Particle i
 INTEGER :: j                 ! Loop control index. Particle j
 INTEGER :: k                 ! Index of the number of steps of sq-well interaction.
!PENDENT TREBALLAR AMB NBLIST
 ! Initialize total energy accumulators
 epotgo = 0.
 epotfis = 0.
! Iterate over all atoms
 DO j = 2,natom
 DO i = 1,j-1
   dist = sqrt(calcDist2DP(r(i), r(j)))
   IF (stepPts(i,j)%tipInt == SS) THEN
     ! Count GO energy term only if distance is appropiate
     IF (dist > stepPts(i,j)%step(1)%r.AND.dist < stepPts(i,j)%step(2)%r) epotgo = epotgo - ego
   ELSE IF (stepPts(i,j)%active) THEN
     ! Physical pot energy is active so sum its energy
     k = stepPts(i,j)%nstep
      overdiscontinuities:DO
       ! The following IF prevents us from summing wrong steps
       IF ( .NOT.((k >= 1).and. dist < stepPts(i,j)%step(k)%r )) EXIT
       epotfis = epotfis - stepPts(i,j)%step(k)%e / FACTE
       k = k - 1
      END DO overdiscontinuities
 !  changed above condition from k>1 to (k >= 1)
 !  this seems not to be needed any longer
    !  IF (dist < stepPts(i,j)%step(k)%r) then
    !  epotfis = epotfis - stepPts(i,j)%step(k)%e / FACTE
    !  write(0,*) "if raro"
    !  end if
   END IF
 END DO
 END DO

END SUBROUTINE calcEpot
!===============================================================================
PURE FUNCTION calcEkin (v, xm, natom) RESULT (ekin)
 !
 ! Purpose:
 ! To compute kynetic energy of the system
 !
 USE geometry
 !
 IMPLICIT NONE
 ! Data ditctionary and variable declaration
  INTEGER, INTENT(IN) :: natom          ! Total number of particles in the system
  REAL(SGL) :: ekin                     ! Result: total kynetic energy
  TYPE(point), INTENT(IN) :: v(natom)   ! Velocities of all atoms
  REAL(SGL), INTENT(IN) :: xm(natom)    ! List of masses of all atoms
  ! Auxiliar block
  REAL(SGL), PARAMETER :: a2 = 1.e-20   ! Armstrong to meters convertion
  INTEGER :: i                          ! Loop control index
  ekin = 0._SGL
  !
  DO i = 1,natom
     ekin = ekin + 0.5 * xm(i) * a2 * dot(v(i), v(i))
 END DO
 END FUNCTION calcEkin
!===============================================================================
END MODULE energy

