!
! $Date: 2011-04-01 14:38:02 $
! $Id: geometryDP.f,v 1.1 2011-04-01 14:38:02 gelpi Exp $
! $Revision: 1.1 $
!
 MODULE geometryDP
!
! Purpose:
! To dispose of functions to deal with R3 vectors defined as DOUBLE precision real vectors.
!
! Record of revisions:
!      DATE            PROGRAMMER           DESCRIPTION OF CHANGE
!    =======          ============         =======================
!     ??              Agustí, JL             Original Code
!    24/08/11         Pedro                  Comments added
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC operator(+), operator(*),operator(-), SPtoDP,  &
       DPtoSP, dotDP, dotDP2, crossDP, moduleDP,      &
       makeUnitDP, cosangDP, calcDistDP, calcDist2DP, &
       calcCM, calcSignCosAng
!
! Data dictionary and variables declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)  ! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6)   ! Single Precision real number definition processor independent

 TYPE, PUBLIC :: pointDP
   REAL(DBL) :: x      ! x cartesian coordinate, double presicion
   REAL(DBL) :: y      ! y cartesian coordinate, double presicion
   REAL(DBL) :: z      ! z cartesian coordinate, double presicion
 END TYPE pointDP

 INTERFACE operator (+)
  MODULE PROCEDURE sumavecDP
 END INTERFACE

 INTERFACE operator (-)
  MODULE PROCEDURE negvecDP
  MODULE PROCEDURE restavecDP
 END INTERFACE

 INTERFACE operator (*)
  MODULE PROCEDURE prodFactDP
  MODULE PROCEDURE prodFactDP2
  MODULE PROCEDURE dotDP
 END INTERFACE

CONTAINS

PURE FUNCTION SPtoDP (rsp) RESULT (rdp)
USE geometry
!
! Purpose:
! To switch from single precision R3 vector to double precision R3 vector
!
IMPLICIT NONE
! Data dictionary and variable declaration
 TYPE(point), INTENT(IN) :: rsp        ! Real single precision original point
 TYPE(pointDP) rdp                     ! Result: Real double precision point
 rdp%x = DBLE(rsp%x)
 rdp%y = DBLE(rsp%y)
 rdp%z = DBLE(rsp%z)
END FUNCTION SPtoDP

PURE FUNCTION DPtoSP (rdp) RESULT (rsp)
USE geometry
!
! Purpose:
! To switch from double precision R3 vector to single precision R3 vector
!
IMPLICIT NONE
! Data dictionary and variable declaration
 TYPE(point) rsp                       ! Result: Real single precision point
 TYPE(pointDP), INTENT(IN) :: rdp      ! Initial double precision point
 rsp%x = SNGL(rdp%x)
 rsp%y = SNGL(rdp%y)
 rsp%z = SNGL(rdp%z)
END FUNCTION DPtoSP

PURE FUNCTION sumavecDP (v1,v2)
!
! Purpose:
! To sum 2 double precision vectors. It shoulb be available only through + operator interface
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP):: sumavecDP             ! Result: sum of two vectors
TYPE (pointDP), INTENT(IN):: v1,v2     ! Initial vectors, v1, v2
sumavecDP%x = v1%x + v2%x
sumavecDP%y = v1%y + v2%y
sumavecDP%z = v1%z + v2%z
END FUNCTION sumavecDP

PURE FUNCTION prodFactDP (f,v)
!
! Purpose:
! To multiply by the left-hand side a double precision vector by a single precision scalar. It should
! be only available through the * operator interface
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP) :: prodFactDP           ! Result: multiplication of a vector by scalar
TYPE (pointDP), INTENT(IN):: v         ! Initial vector
REAL(SGL), INTENT(IN):: f              ! Initial scalar
prodFactDP%x = v%x * f
prodFactDP%y = v%y * f
prodFactDP%z = v%z * f
END FUNCTION prodFactDP

PURE FUNCTION prodFactDP2 (v,f)
!
! Purpose:
! To multiply by the right-hand side a double precision vector by a single precision scalar. It should
! be only available through the * operator interface
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP) :: prodFactDP2          ! Result: multiplication of a vector by scalar
TYPE (pointDP), INTENT(IN):: v         ! Initial vector
REAL(SGL), INTENT(IN):: f              ! Initial scalar
prodFactDP2%x = v%x * f
prodFactDP2%y = v%y * f
prodFactDP2%z = v%z * f
END FUNCTION prodFactDP2

PURE FUNCTION negvecDP (v)
!
! Purpose:
! To change sign to a double precision vector. It should be available only
! through the - operator's interface
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP) :: negvecDP             ! Result: opposite sign vector
TYPE (pointDP), INTENT(IN) :: v        ! Initial vector
negvecDP = prodFactDP (-1.0, v)
END FUNCTION negvecDP

PURE FUNCTION restaVecDP (v1,v2)
!
! Purpose:
! To substract to double precision vectors. It should be available only through
! the - operator's interface
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP) :: restaVecDP           ! Resulting vector from subtraction
TYPE (pointDP), INTENT(IN) :: v1, v2   ! Initial vectors
restaVecDP = v1 + (-v2)
END FUNCTION restaVecDP

PURE FUNCTION dotDP (v1,v2)
!
! Purpose:
! To compute the dot product between to double precision vectors.
! result is a single precision scalar.
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP), INTENT(IN):: v1,v2     ! Initial vectors
REAL(SGL) :: dotDP                     ! Result scalar, dot product
dotDP = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
END FUNCTION dotDP

PURE FUNCTION dotDP2 (v1,v2)
!
! Purpose:
! To compute the dot product between to double precision vectors.
! result is a DOUBLE precision scalar.
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP), INTENT(IN):: v1,v2     ! Initial double precision vectors
REAL(DBL) :: dotDP2                    ! Resulting dot product, double presicion
dotDP2 = v1%x * v2%x + v1%y * v2%y + v1%z * v2%z
END FUNCTION dotDP2

PURE FUNCTION crossDP(v1,v2)
!
! Purpose:
! To compute cross product between double precision vectors. Result is also a double presicion vector
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP) crossDP                 ! Result, cross product
TYPE (pointDP), INTENT(IN) :: v1,v2    ! Initial vectors
crossDP=pointDP(-v1%z*v2%y + v1%y*v2%z, v1%z*v2%x - v1%x*v2%z, -v1%y*v2%x + v1%x*v2%y)
END FUNCTION crossDP

PURE FUNCTION moduleDP (v)
!
! Purpose:
! To obtain the module of a double precision vector.
!
IMPLICIT NONE
! Data dictionary and variable declaration
REAL(SGL) :: moduleDP                  ! Result: module of double presicion vector
TYPE (pointDP), INTENT(IN) :: v        ! Initial vector
moduleDP = SQRT(v*v)
END FUNCTION moduleDP

PURE FUNCTION makeUnitDP (v)
!
! Purpose:
! To make a double presicion vector unitary
!
IMPLICIT NONE
! Data dictionary and variable declaration
TYPE (pointDP):: makeUnitDP            ! Result: unitary double precision vector
TYPE (pointDP), INTENT(IN) :: v        ! Initial vector
! ATT vector cannot be null vector
makeUnitDP = ( 1./ moduleDP(v) )*v
END FUNCTION makeUnitDP

PURE FUNCTION cosangDP (v1,v2)
!
! Purpose:
! To compute the cosinus of the angle between two double precision vectors
!
IMPLICIT NONE
! Data dictionary and variable declaration
REAL(SGL) :: cosangDP                  ! Result: COS of angle
TYPE (pointDP), INTENT(IN) :: v1,v2    ! Initial vectors
cosangDP= (v1 * v2) / moduleDP(v1) / moduleDP(v2)
END FUNCTION cosangDP

PURE FUNCTION calcDist2DP (p1,p2)
!
! Purpose:
! To compute the square of the distance between 2 doble precision vectors
!
IMPLICIT NONE
! Data dictionary and variable declaration
  REAL(SGL) :: calcDist2DP             ! Result: Squared distance UNITS
  TYPE (pointDP), INTENT(IN) :: p1,p2  ! Initial vectors
  TYPE (pointDP) :: v                  ! Auxiliar vector holding intial vector substraction
  v = restaVecDP(p1,p2)
  calcDist2DP = v * v
END FUNCTION calcDist2DP

PURE FUNCTION calcDistDP (p1,p2)
!
! Purpose:
! To compute the distance between 2 doble precision vectors
!
IMPLICIT NONE
! Data dictionary and variable declaration
  REAL :: calcDistDP                   ! Result: distance UNITS
  TYPE (pointDP), INTENT(IN) :: p1,p2  ! Initial vectors
  calcDistDP = SQRT( calcDist2DP (p1,p2) )
END FUNCTION calcDistDP

! ADDED

PURE FUNCTION calcSignCosAng(r, v, i, j, natom)
 USE geometry
!
! Purpose:
! To compute the sign of the cosinus of the angle between relative position (distance)
! and relative velocity
!
IMPLICIT NONE
! Data dictionary and variable declaration
 REAL(SGL) calcSignCosAng              ! Result: sign of the cosinus of the angle
 INTEGER, INTENT(IN) :: natom          ! Total number of particles in the system
 INTEGER, INTENT(IN) :: i              ! Index of particle i (1-natom)
 INTEGER, INTENT(IN) :: j              ! Index of particle j (1-natom)
 TYPE(pointDP), INTENT(IN) :: r(natom) ! Coordinates of all the particles of the system
 TYPE(point), INTENT(IN) :: v(natom)   ! Velocities of all the particles of the system
!
 REAL(SGL) :: dotv                     ! Dot product between dr and dv
 TYPE(pointDP) :: dr                   ! Difference in coordinates: R3 Vector DP. Distance
 TYPE(point) :: dv                     ! Difference in velocities: Relative velocity
!
 dr = r(j) - r(i)
 dv = v(j) - v(i)
 dotv = dr%x * dv%x + dr%y * dv%y + dr%z * dv%z
 calcSignCosAng = SIGN(1.,dotv) ! SIGN(A,B) returns the value A with sing of B
END FUNCTION calcSignCosAng

PURE FUNCTION calcCM (natom, r, xm) RESULT (rcm)
  USE geometry
!
! Purpose:
! To obtain the center of mass of a group of natom particles.
!
IMPLICIT NONE
! Data dictionary and variable declaration
 TYPE(point) rcm                       ! Result: coordinates of the CM
 INTEGER,INTENT(IN) :: natom           ! Total number of particles in the system
 TYPE(pointDP), INTENT(IN) :: r(natom) ! Cartesian coordinates of all particles
 REAL(SGL), INTENT(IN) :: xm(natom)    ! Mass of every particle
!
  INTEGER :: i                         ! Loop control index
  REAL(SGL) xmassa                     ! Total mass of the system
!
! Initialization of CM
  rcm = point(0.,0.,0.)
  xmassa = sum(xm)
!
  DO i = 1,natom
! explicit line code due to (type point + type pointDP) operation is not defined
    rcm%x = rcm%x + xm(i) * r(i)%x
    rcm%y = rcm%y + xm(i) * r(i)%y
    rcm%z = rcm%z + xm(i) * r(i)%z
  END DO
! CM
  rcm = (1./xmassa) * rcm
END FUNCTION calcCM

END MODULE geometryDP