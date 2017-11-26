 MODULE intList
 USE stepPotentials
!
! Purpose:
! This module contains the derived-type data structure to keep interaction information.
! Non bonded interactions are computed only within a distance of the particle i. This
! is periodically updated.
! Within that distance you can find particles that possibly will collide with particle i.
! Between this particles and i, interaction data is computed and stored in the intpList
! derived data type.
! You keep in this data type:
! number of particle that will interact with i (1-natom) in patnum
! number of particle interacting seen from other particle (j) (1-nblist(j)%nats) in simp
! type of sq well interaction in stepPt
! ???? in xsum
! time of the shortest interaction (timp), kin energy transfer in (deltak)
!
! Record of revisions:
!      DATE            PROGRAMMER           DESCRIPTION OF CHANGE
!    =======          ============         =======================
!     ??              Agustí, JL             Original Code
!    31/08/11         Pedro                  Comments added
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC intpList, allocateintPList, intData
!
! Data dictionary and variable declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)  ! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6)   ! Single Precision real number definition processor independent
! Meaningful
TYPE intData                      ! Single Interaction data storage
    INTEGER :: patnum             ! Particle number that could collide with i
    INTEGER :: simp               ! Particle number that could collide with i
                                  ! referred to j
    TYPE(stepPotInt) :: stepPt    ! sq well interaction with j
    REAL(SGL) :: xsum             ! sum of masses
    REAL(DBL) :: timp             ! time of the firt collsion to occour involving i
    REAL(SGL) :: deltak           ! kinetic energy transfer in collision
END TYPE intData
 
TYPE intpList                     ! All number of possible interactions
                                  ! are stored in this data type
    INTEGER nats                  ! total number of possible interactions
                                  ! involving i (counter)
    TYPE(intData), ALLOCATABLE :: iData(:) ! see above
END TYPE intpList
 
CONTAINS
 
FUNCTION allocateintPList(natom, ioerr) RESULT (pl)
 !
 ! Purpose:
 ! To allocate iData information for every particle (1-natom)
 ! for every particle we allocate natom possible interactions.
 ! This function will define nblist(i): non bonded list of particle i
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 INTEGER, INTENT(IN) :: natom
 INTEGER, INTENT(OUT) :: ioerr
 TYPE (intpList) pl
 ! This is called as nblist(i)=allocateintPList(natom, ieorr)
  ALLOCATE (pl%iData(natom), stat=ioerr)
 ! Reset counter for atom i nats is a scalar here.
  pl%nats=0
END FUNCTION allocateintPList
 
END MODULE intList