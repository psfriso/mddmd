 MODULE colisio
!
! Purpose:
!
! Record of revisions:
!      DATE            PROGRAMMER           DESCRIPTION OF CHANGE
!    =======          ============         =======================
!     ??              Agustí, JL             Original Code
!    01/09/11         Pedro                  Defined as module. Comments added.
!    01/09/11         Pedro                  Moved subroutines calcCM, CalcSingCosAng
!    05/09/11         Pedro                  TMIN upgraded to double precision
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC colisioBond, colisioNonBond, chgmom, &
        updateTCol, updateV, getMinTCol,     &
        nextCol, inici, updateXocPart
!
! Data dictionary and variable declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)  ! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6)   ! Single Precision real number definition processor independent
!
CONTAINS
!
 SUBROUTINE colisioBond(bpairs, blist, r, v, rb, drb, xm, xsum, natom)
 USE geometry
 USE geometryDP
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 INTEGER, INTENT(IN) :: natom            ! Total number of particles in the system
 INTEGER, INTENT(IN) :: bpairs           ! Total number of bonded pairs ( non-permited
                                         ! energetic change, inf sq well)
 INTEGER, INTENT(IN) :: blist(bpairs,2)  ! List of pairs of rigid bonded particles
 TYPE (pointDP), INTENT(IN) :: r(natom)  ! Coordinates of all particles
 REAL(SGL), INTENT(IN) :: rb(natom,natom)! Distance between centers of sq-wells of i-j
 REAL(SGL), INTENT(IN) :: drb(natom,natom) ! Half of the amplitude of the sq-well between i-j
 REAL(SGL), INTENT(IN) :: xm(natom)      ! List of masses of all the particles
 REAL(SGL), INTENT(IN) :: xsum(natom,natom) ! Function of masses of i and j
                                            ! = ( mi + mj ) / ( 2 mi mj )
 TYPE(point), INTENT(INOUT) :: v(natom)  ! Velocities of all partilces
! Auxiliar Block
   REAL(SGL) :: csang       ! Sign of Cos( angle(r(x,y,z),v(x,y,z) ) )
  !REAL calcSignCosAng
   REAL(SGL) :: rbmin       ! Lowest energy allowed distance between two particles
   REAL(SGL) :: rbmax       ! Maximum energy allowed distance between two particles
   REAL(SGL) :: rij         ! Distance (scalar) between i and j
   INTEGER :: i             ! Bonded particle i (of pair i-j )
   INTEGER :: j             ! Bonded particle j (of pair i-j )
   INTEGER :: np            ! Loop control index. Bonded pair number.
   !
  alloverbpairs: DO np = 1,bpairs
  ! Go over all bonded pairs
      i = blist(np,1)
      j = blist(np,2)
  ! Compute distances and bond properties
      csang = calcSignCosAng(r, v, i, j, natom)
      rij = calcDistDP(r(i), r(j))
      rbmin = rb(i,j) - drb(i,j)
      rbmax = rb(i,j) + drb(i,j)
  !
      violatingdist:IF((rij > rbmax .and. csang > 0.).or.(rij < rbmin.and.csang < 0.))  THEN
      ! Not-allowed position. Change their velocities. Collision.
         CALL chgmom(i, j, xm, xsum(i,j), r, v, natom)
      ENDIF violatingdist
   END DO alloverbpairs
 END SUBROUTINE colisioBond
!==============================================================================
!
 SUBROUTINE colisioNonBond(stepPts, temps, nblist, r, v, rhc, ind2, atom, xm, ierr, &
               TMIN, natom, isolv)
 USE geometry
 USE geometryDP
 USE intList
 USE stepPotentials
!
IMPLICIT NONE
! Data dictionary and variable declaration
   INTEGER, INTENT(IN) :: natom        ! Total number of particles in the system
   INTEGER, INTENT(IN) :: ind2(natom)  ! Residue number for all particles
   TYPE(pointDP), INTENT(IN) :: r(natom) ! Coordinates of all atoms
   REAL(SGL), INTENT(IN) :: rhc(natom) ! Hardcore radii
   REAL(SGL), INTENT(IN) :: xm(natom)  ! Masses
   REAL(DBL), INTENT(IN) :: TMIN       ! Minimum time to consider collision
   CHARACTER(len=4),INTENT(IN) :: atom(natom)  ! Atom type
   LOGICAL, INTENT(in) :: isolv        ! Switch Non-bonded interactions
!
   TYPE(stepPotInt), INTENT(INOUT) :: stepPts(natom, natom) ! Step potential interactions
   TYPE(intPList), INTENT(INOUT) :: nblist(natom)  ! Non bonded list for all particles
   INTEGER, INTENT(INOUT) :: ierr      ! Overlaps: errors control counter
   TYPE(point), INTENT(INOUT) :: v(natom) ! Velocities all particles
!
   REAL(SGL) :: csang         ! Sign of Cos(angle(r.v))
   REAL(SGL) :: dmin          ! Sum of hard-core radius (below dmin the
                              ! particles are in steric clash)
   REAL(SGL) :: rij           ! Squared distance between i and j
!   REAL(SGL) :: calcSignCosAng
   INTEGER :: i               ! Loop control index
   INTEGER :: j               ! Loop control index
   INTEGER :: k               ! Loop control index
   LOGICAL :: isalt           ! Logical variable to prevent particles non
                              ! bonded in amide bond to collide almost exclusively
   REAL(DBL) :: temps         ! Time evolved until collision (??)

! removes overlaps
   DO i = 1, natom-1
   DO k = 1, nblist(i)%nats
      j = nblist(i)%iData(k)%patnum
      IF (j > i) THEN ! Fem nomes la meitat perque les llistes son simetriques
         csang = calcSignCosAng(r, v, i, j, natom)
         dmin = (rhc(i) + rhc(j))**2
! evita els xocs dels molt adjacents
         isalt = (ind2(j).eq.ind2(i)+1).and.( &
                 (atom(i).eq.'C'.and.atom(j).eq.'CB').or. &
                 (atom(i).eq.'CB'.and.atom(j).eq.'N').or. &
                 (atom(i).eq.'O'.and.atom(j).eq.'CA'))
         IF (.NOT.isalt) THEN
         ! Only if atoms are NOT close neighboors
            rij = calcDist2DP(r(i), r(j))
            IF (rij < dmin .and. csang < 0) THEN
            ! They are overlapping, and getting closer
            ! change velocity is needed
               CALL chgmom(i, j, xm, nblist(i)%iData(k)%xsum, r, v, natom)
            ! This is an  " error ", count them
               ierr=ierr+1
            END IF
         END IF
      END IF
   END DO
   END DO
! Once removed overlaps, now it
! calculates colision times
   DO i = 1, natom-1
     DO k = 1, nblist(i)%nats
       j = nblist(i)%iData(k)%patnum
       IF (j > i) THEN ! ATT
        IF (stepPts(i,j)%tipInt == SS.or.isolv ) THEN
        ! Now update collision times for all particles
           CALL updateTCol (i, j, temps, r, v, nblist(i)%iData(k)%stepPt, TMIN, &
             nblist(i)%iData(k)%timp, nblist(i)%iData(k)%deltak, natom)
           nblist(j)%iData(nblist(i)%iData(k)%simp)%timp = nblist(i)%iData(k)%timp
           nblist(j)%iData(nblist(i)%iData(k)%simp)%deltak = nblist(i)%iData(k)%deltak
       END IF
      END IF
     END DO
   END DO
 END SUBROUTINE colisioNonBond
!
!===============================================================================
 SUBROUTINE colisioNonBondOLD(temps, nblist, r, v, rhc, ind2, &
                            atom, xm, ierr,TMIN, natom)
 USE geometry
 USE geometryDP
 USE intList
!
IMPLICIT NONE
!
   INTEGER, INTENT(IN) :: natom
   TYPE(intPList), INTENT(INOUT) :: nblist(natom)
   INTEGER, INTENT(IN) :: ind2(natom)
   TYPE (pointDP), INTENT(IN) :: r(natom)
   REAL(SGL), INTENT(IN) :: rhc(natom)
   REAL(SGL), INTENT(IN) :: xm(natom)
   REAL(DBL), INTENT(IN) :: TMIN
   CHARACTER(len=4), INTENT(IN) :: atom(natom)
!
   INTEGER, INTENT(INOUT) :: ierr
   TYPE(point), INTENT(INOUT) :: v(natom)
   
   REAL(SGL) :: csang
   REAL(SGL) :: dmin
   REAL(SGL) :: rij
   !REAL calcSignCosAng
   INTEGER :: i
   INTEGER :: j
   INTEGER :: k
   LOGICAL :: isalt
   REAL(DBL) :: temps
!
   DO i = 1, natom-1
   DO k = 1, nblist(i)%nats
      j = nblist(i)%iData(k)%patnum
      IF (j > i) THEN ! Fem nomes la meitat perque les llistes son simetriques
         csang = calcSignCosAng(r, v, i, j, natom)
         dmin = rhc(i) + rhc(j)
! evita els xocs dels molt adjacents !
         isalt = (ind2(j) == ind2(i)+1).AND.( &
                 (atom(i) == 'C'.AND.atom(j) == 'CB').OR. &
                 (atom(i) == 'CB'.AND.atom(j) == 'N').OR. &
                 (atom(i) == 'O'.AND.atom(j) == 'CA'))
         IF (.NOT.isalt) THEN
            rij = calcDist2DP(r(i), r(j))
            IF (rij < dmin.and.csang < 0) THEN
               CALL chgmom(i, j, xm, nblist(i)%iData(k)%xsum, r, v, natom)
               ierr=ierr+1
            END IF
         END IF
         CALL updateTCol (i, j, temps, r, v, nblist(i)%iData(k)%stepPt, TMIN, &
            nblist(i)%iData(k)%timp, nblist(i)%iData(k)%deltak, natom)
         nblist(j) % iData( nblist(i) % iData(k) %simp) %timp = nblist(i)%iData(k)%timp
         nblist(j)%iData(nblist(i)%iData(k)%simp)%deltak = nblist(i)%iData(k)%deltak
      END IF
   END DO
   END DO
 END SUBROUTINE colisioNonBondOLD
!===============================================================================
 PURE SUBROUTINE chgmom(mem1, mem2, xm, xsum, r, v, natom)
 !
 ! Purpose:
 ! To change velocity of particles mem1, mem2 by transfering linear momemtum,
 ! when potential energy change is forbbiden
 !
 USE geometry
 USE geometryDP
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
   INTEGER, INTENT(IN) :: natom         ! Total number of particles in the system
   REAL(SGL), INTENT(IN) :: xsum        ! Function of masses i, j ( mi+ mj) / 2mi mj
   REAL(SGL), INTENT(IN) :: xm(natom)   ! List of masses of all atoms
   TYPE(pointDP), INTENT(IN) :: r(natom)! Coordinates of all atoms
   INTEGER, INTENT(IN) :: mem1          ! First particle involved in collision
   INTEGER, INTENT(IN) :: mem2          ! Second particle involved in collsion
!
   TYPE(point), INTENT(INOUT) :: v(natom) ! Velocities of all particles
!
   TYPE(pointDP) :: dr                  ! Relative position (R3 vector)
   TYPE(point) :: dv                    ! Relative velocity (R3 vector)
   REAL(SGL) :: dp                      ! Function of transferred momentum
!
! Compute relative distance and relative velocity
   dr = r(mem2) - r(mem1)
   dv = v(mem2) - v(mem1)
! Velocity components will only change in the direction of dr
   dp = -(dr%x * dv%x + dr%y * dv%y + dr%z * dv%z) / dotDP(dr, dr) / xsum
! modul del moment transferit en la colisio dr dona la direccio
   v(mem1)%x = v(mem1)%x - dp / xm(mem1) * dr%x
   v(mem1)%y = v(mem1)%y - dp / xm(mem1) * dr%y
   v(mem1)%z = v(mem1)%z - dp / xm(mem1) * dr%z

   v(mem2)%x = v(mem2)%x + dp / xm(mem2) * dr%x
   v(mem2)%y = v(mem2)%y + dp / xm(mem2) * dr%y
   v(mem2)%z = v(mem2)%z + dp / xm(mem2) * dr%z
 END SUBROUTINE chgmom
!===============================================================================
 SUBROUTINE updateTCol(i ,j, temps, r, v, stepPt, TMIN, timp, deltak, natom)
! input r v de dues particules, output deltak i timp
 USE geometry
 USE geometryDP
 USE stepPotentials
 USE intList
 !
 ! Purpose:
 ! To update collision times.
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 INTEGER, INTENT(IN) :: natom           ! Total number of particles in the system
 INTEGER, INTENT(IN) :: i               ! Particle number possibly colliding with j (1-natom)
 INTEGER, INTENT(IN) :: j               ! Particle number possibly colliding with i (1-natom)
 TYPE(pointDP), INTENT(IN) :: r(natom)  ! Coordinates of all particles
 TYPE(point), INTENT(IN) :: v(natom)    ! Velocities of all particles
 TYPE(stepPotInt), INTENT(IN) :: stepPt ! Step-potential interaction between i and j
 REAL(DBL), INTENT(IN) :: temps         ! Time elapsed until event
 REAL(DBL), INTENT(IN) :: TMIN          ! Minimum time to consider collision
!
 REAL(DBL), INTENT(OUT) :: timp         ! Time of the first event
 REAL(SGL), INTENT(OUT) :: deltak       ! Energy tranferred during event
!
! Auxiliar block
 INTEGER :: k                           ! Loop control index. Number of steps in
                                        ! sq-well interaction
 REAL(DBL), DIMENSION(MAXSTEPS) :: argk ! Argument of sq root
 REAL(DBL), DIMENSION(MAXSTEPS,2)  :: tijs ! Collision times until events.
                                        ! we have to times per sq-well discontinuity
                                        ! as it can cross it in opposite directions
                                        ! without " bouncing" (esquema plano de agustí)
 REAL(SGL) :: rij                       ! Distance between i and j (scalar)
 REAL(SGL) :: rij2                      ! Squared distance between i and j (scalar)
 REAL(SGL) :: vij                       ! Module of relative velocity (scalar)
 REAL(SGL) :: vij2                      ! Squared module of relative velocity (scalar)
 REAL(SGL) :: a                         ! decomposition of grade 2 polinimial solution
 REAL(SGL) :: b                         ! decomposition of grade 2 polinimial solution
 REAL(SGL) :: dotrv                     ! Projection between relative velocity
                                        ! and relative position
 TYPE(pointDP) :: dr                    ! Relative postion (R3 distance vector)
 TYPE(point) :: dv                      ! Relative velocity (R3)
 INTEGER, DIMENSION(2) :: lc            ! Index ( integer, integer) where next collision is in the list
!
   dr = r(j) - r(i)
   dv = v(j) - v(i)
   rij2 = dotDP(dr,dr)
   rij = sqrt(rij2)
   vij2 = dot(dv,dv)
   vij = sqrt(vij2)
   dotrv = dr%x * dv%x + dr%y * dv%y + dr%z * dv%z ! dot function cannot yet
   ! be used because is not defined for dot(type ponint,type pointDP)
   tijs=0._DBL
  itsinteracting: IF (rij2.LE.stepPt%step(stepPt%nstep)%r**2 + dotrv**2/vij2) THEN
   ! Particles are close enough to interact with the outtermost discontinuity
   ! we shall look for this possible interaction
   ! We look over all k discontinuities
     DO k = 1,stepPt%nstep
     ! compute the argument of the squared root (solution of times of collision)
        argk(k) = stepPt%step(k)%r**2 - rij2 + dotrv**2/vij2
       collisionpossible: IF (argk(k) > 0._SGL) THEN
        ! if argument is possitive then we have two possible solutions (two possible collisions)
           a = -dotrv / vij2
           b = sqrt(argk(k)) / vij
        ! compute the time for this two solutions obtained, and record them
           tijs(k,1) = a - b
           tijs(k,2) = a + b
        END IF collisionpossible
     END DO
     ! now check wich one is the first to occour.
     lc = MINLOC(tijs, mask= tijs > TMIN)
  !   write(0,*) " little times" , count(tijs < TMIN)
    collisionfound: IF (lc(1) > 0)  THEN
     ! At least one possible collision appeared. take note of time
        timp = temps + tijs(lc(1),lc(2))
      energysign:  IF (lc(2) == 1) THEN
        ! This means that particle is entering in the well
        ! energy transfer conserve its sign
           deltak =  +stepPt%step(lc(1))%e
        ELSE
        ! if particle is leaving, energy transfer is opposite
        ! as defined
           deltak = -stepPt%step(lc(1))%e
        END IF energysign
     END IF collisionfound
   END IF itsinteracting
 END SUBROUTINE updateTCol
!===============================================================================
 PURE SUBROUTINE updateV (r, v, deltak, xm, xsum, mem1, mem2, natom)
 USE geometry
 USE geometryDP
 !
 ! Purpose:
 ! To change velocity of particles mem1, mem2 by transfering linear momemtum,
 ! when potential energy change is allowed
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
   INTEGER, INTENT(IN) :: natom         ! Total number of particles in the sistem
   INTEGER, INTENT(IN) :: mem1          ! Fisrt particle involved in  a collision
   INTEGER, INTENT(IN) :: mem2          ! Second particle involved in a collision
   TYPE(pointDP), INTENT(IN) :: r(natom)! Coordinates of all atoms
   REAL(SGL), INTENT(IN) :: deltak      ! Potential energy depth (sq-well)
   REAL(SGL), INTENT(IN) :: xm(natom)   ! List of masses of all atoms
   REAL(SGL), INTENT(IN) :: xsum        ! Function of masses i, j ( mi+ mj) / 2mi mj
   TYPE(point), INTENT(INOUT) :: v(natom) ! Velocities of all particles
!
   REAL(SGL), PARAMETER :: A2 = 1.e-20_SGL ! Armstrong2 to meter2
   TYPE(pointDP) :: dr                  ! Relative position
   TYPE(point) :: dv                    ! Relative velocity
   REAL(DBL) :: rij                     ! Distance (scalar) between i j
   REAL(SGL) :: vdmod                   ! Projected velocity on i-j coordinates axis
   REAL(SGL) ::  sto                    ! Argument of squared root
   REAL(SGL) :: dp                      ! Transferred linear momemtum
! Compute relative distance and velocity
   dr = r(mem2) - r(mem1)
   dv = v(mem2) - v(mem1)
! calculo vdmod, la projeccio de la diferencia de velocitats en l'eix que uneix les particules
   rij = sqrt(dotDP(dr,dr))
   vdmod = (dv%x * dr%x + dv%y * dr%y + dv%z * dr%z) / rij
! entra o surt d'un pou
   sto = vdmod**2 + 4. * deltak * xsum / A2
   IF (sto > 0) THEN
! vario la velocitat
      dp = -vdmod + sign(1.,vdmod) * sqrt(sto)
      dp = dp / 2. / xsum / rij
   ELSE
! les particules es queden atrapades al pou
      dp = -vdmod / xsum / rij
   END IF
   v(mem1)%x = v(mem1)%x - dp / xm(mem1) * dr%x
   v(mem1)%y = v(mem1)%y - dp / xm(mem1) * dr%y
   v(mem1)%z = v(mem1)%z - dp / xm(mem1) * dr%z

   v(mem2)%x = v(mem2)%x + dp / xm(mem2) * dr%x
   v(mem2)%y = v(mem2)%y + dp / xm(mem2) * dr%y
   v(mem2)%z = v(mem2)%z + dp / xm(mem2) * dr%z
 END SUBROUTINE updateV
!==============================================================================
 SUBROUTINE getMinTCol (m1, nblist, tpart, ipart, natom)
 USE intList
 !
 ! Purpose:
 ! To get the first collision per atom to occour. This subroutine creates
 ! ipart and tpart list by finding the smallest time of interaction
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
   INTEGER, INTENT(IN) :: m1            ! Particle number to look for its first collision (1-natom)
   INTEGER, INTENT(IN) :: natom         ! Total number of particles in the system
   TYPE (intpList), INTENT(IN) :: nblist(natom) ! Non bonded list for all atoms
   REAL(DBL), INTENT(OUT) :: tpart      ! Minimum collision time for m1 particle
   INTEGER, INTENT(OUT) :: ipart        ! Index of the particle that is going to interact
                                        ! first with m1 (referred to m1)
   INTEGER :: i                         ! Loop control index
!  actualitza la llista de temps
   tpart = 1.
   ipart = 0
   DO i = 1,nblist(m1)%nats
 ! For all nats particles interacting with m1 look for the smallest
 ! interaction time
      IF (nblist(m1)%iData(i)%timp < tpart ) THEN ! this time is smaller
      ! keep track of this particle and its time to collision
         tpart = nblist(m1)%iData(i)%timp
         ipart = i 
      END IF
   END DO
 END SUBROUTINE getMinTcol
!============================================================================================
! pure subroutine nextCol(mem1, mem2, np, tevent, ipart, tpart, natom,nblist)
 SUBROUTINE nextCol(mem1, mem2, np1, np2, tevent, ipart, tpart, natom,nblist)
 USE intList
 !
 ! Purpose:
 ! Check and retrieve the first collision to occour.
 ! In the list of first collision time per particle (tpart) look for the
 ! smallest collision time. Then keep the number of that particle and details
 ! of collision.
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
   INTEGER, INTENT(IN) :: natom         ! Total number of particles in the system
   INTEGER, INTENT(IN) :: ipart(natom)  ! List of particles. This is linked by location
                                        ! with tpart
   REAL(DBL), INTENT(IN) :: tpart(natom)! List of smallest time of collision per particle
                                        ! tpart and i part contain first collision per atom
   TYPE (intpList), INTENT(IN) :: nblist(natom) ! Non bonded list of all particles (list of all possible
                                                ! interactions)
   INTEGER, INTENT(OUT) :: mem1         ! First particle colliding (1-natom)
   INTEGER, INTENT(OUT) :: mem2         ! Second particle colliding (1-natom)
   REAL(DBL), INTENT(OUT) :: tevent     ! Time elapsed until collision
   INTEGER, INTENT(OUT) :: np1          ! First particle colliding referred to local mem1 enviroment
   INTEGER, INTENT(OUT) :: np2          ! Second particle colliding referred to local mem1 enviroment
 ! Get which is the partilce with smallest collision time
   mem1 = minloc(tpart,1)! ,mask= tpart > 1.e-22) ¿no debería aplicarse la mascara tambien?
   np1 = ipart(mem1)     ! point out which one is
   mem2 = nblist(mem1)%iData(np1)%patnum ! to whom is going to collide (1-natom)
   tevent = nblist(mem1)%iData(np1)%timp ! time elapsed until collision
   np2 = nblist(mem1)%iData(np1)%simp    ! to whom is going to collide, referred to mem1
! codigo viejo para revisar
!   mem2 = nblist(mem1)%iData(ipart(mem1))%patnum
!   tevent = nblist(mem1)%iData(ipart(mem1))%timp
!   np2 = nblist(mem1)%iData(ipart(mem1))%simp
 END SUBROUTINE nextCol
!============================================================================================
 SUBROUTINE inici (nblist, tpart, ipart, natom)
 USE intList
 !
 ! Purpose:
 ! To initilize list of first collision per particle.
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
   INTEGER, INTENT(IN) :: natom         ! Total number of particles in the system
   TYPE (intpList), INTENT(IN) :: nblist(natom) ! Non bonded list for all the particles
   REAL(DBL), INTENT(OUT) :: tpart(natom) ! List of time of first collision per particle
   INTEGER, INTENT(OUT) :: ipart(natom)   ! List of partilces of first collision per particle
   INTEGER :: i                         ! Loop control index
! Initialize lists
   tpart = 1.
   ipart = 0
! look smallest-time interaction per particle
   DO i = 1,natom
      CALL getMinTCol(i, nblist, tpart(i), ipart(i), natom)
   END DO
 END SUBROUTINE inici
!================================================================================
 SUBROUTINE updateXocPart(m1, nblist, temps, r, v, TMIN, natom)
 USE intList
 USE geometry
 USE geometryDP
 !
 ! Purpose:
 ! To update all collision times of particle m1 with all its local enviroment (nats)
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
    INTEGER, INTENT(IN) :: natom        ! Total number of particles in the system
    INTEGER, INTENT(IN) ::  m1          ! Particle whose local enviroment is being
                                        ! updated
    REAL(DBL), INTENT(IN) :: TMIN       ! Minimum time left to event below this time
                                        ! collision is not considered
    TYPE(pointDP), INTENT(IN) :: r(natom) ! Coordinates of all particles
    TYPE(point), INTENT(IN) :: v(natom)   ! Velocities of all partilcles
    TYPE(intpList), INTENT(INOUT) :: nblist(natom)! Non bonded list of all particles
                                                  !!OJO INOUT attributed added: previous IN
    REAL(DBL), INTENT(IN) :: temps      ! Time, which one?
    INTEGER :: i                        ! Loop control index
  !
    DO i = 1,nblist(m1)%nats
    ! Update collisions times in m1 local enviroment with particle i
    ! belonging to its.
       CALL updateTCol(m1, nblist(m1)%iData(i)%patnum, temps, r, v, nblist(m1)%iData(i)%stepPt, TMIN,&
           nblist(m1)%iData(i)%timp, nblist(m1)%iData(i)%deltak, natom)
    END DO
  !
 END SUBROUTINE updateXocPart

!===============================================================================   
!
! calc sing cos ang MOVED to geometryDP
!
!===============================================================================   
!
! calc CM MOVED to geometryDP
!
!===============================================================================   
END MODULE colisio
