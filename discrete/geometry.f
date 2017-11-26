!
! $Date: 2011-04-01 14:38:02 $
! $Id: geometry.f,v 1.1 2011-04-01 14:38:02 gelpi Exp $
! $Revision: 1.1 $
!
MODULE geometry
!
! Purpose:
! To dispose of functions to deal with R3 vectors defined as single precision real vectors.
!
! Record of revisions:
!      DATE            PROGRAMMER           DESCRIPTION OF CHANGE
!    =======          ============         =======================
!     ??              Agustí, JL             Original Code
!    23/08/11         Pedro                  Comments added
!    24/08/11         Pedro                  Prvate attribute
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC PI, cross, module, makeUnit, dot, cosang,     &
       operator(+), operator(*),operator(-),         &
       assignment (=) ! calcDist2,calcDist   !  Apparently not needed functions
       ! rotaAng!, pol2car, car2pol, dist2pol, rotquad, rotEuler
!
! DATA dictionary and variable declaration
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6) ! Single Precision real number definition (processor independent)
REAL, PARAMETER :: PI=3.1415926 ! PI number definition

TYPE, PUBLIC:: point
  REAL(SGL) :: x,y,z ! Cartesian coordinates of R3 point
END TYPE point

TYPE,PUBLIC :: pointPolar
  REAL(SGL) :: r,phi,the ! Polar coordinates of R3 point
END TYPE pointPolar

TYPE, PUBLIC :: pointInt
  INTEGER :: i  ! i coordinate of R3 integer point
  INTEGER :: j  ! j coordinate of R3 integer point
  INTEGER :: k  ! k coordinate of R3 integer point
END TYPE pointInt

INTERFACE operator (+)
 MODULE PROCEDURE sumavec
END INTERFACE

INTERFACE operator (-)
 MODULE PROCEDURE negvec
 MODULE PROCEDURE restavec
 MODULE PROCEDURE restaint
END INTERFACE

INTERFACE operator (*)
 MODULE PROCEDURE prodFact
 MODULE PROCEDURE prodFact2
 MODULE PROCEDURE dot
END INTERFACE

INTERFACE assignment (=)
 MODULE PROCEDURE topint
 MODULE PROCEDURE topreal
END INTERFACE

CONTAINS

PURE FUNCTION sumavec (v1,v2)
 !
 ! Purpose:
 ! To add two real-single-precision vectors. It should be only available through + interface operator
 !
  IMPLICIT NONE
 ! Data dictonary and variables declaration
  TYPE (point) :: sumavec
  TYPE (point), INTENT(IN):: v1,v2
  sumavec%x = v1%x + v2%x
  sumavec%y = v1%y + v2%y
  sumavec%z = v1%z + v2%z
END FUNCTION sumavec

PURE FUNCTION prodFact(f,v)
 !
 ! Purpose:
 ! To multiply to real-single-precision vector an scalar by the left-hand-side
 ! It should be only available through * interface operator
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point) :: prodFact      ! Result of th multiplication
 TYPE (point), INTENT(IN):: v  ! Inintial vector
 REAL(SGL), INTENT(IN):: f          ! Scalar to multiply
 prodFact%x = v%x * f
 prodFact%y = v%y * f
 prodFact%z = v%z * f
END FUNCTION prodFact

PURE FUNCTION prodFact2(v,f)
 !
 ! Purpose:
 ! To multiply to real-single-precision vector an scalar by the rigth-hand-side
 ! It should be only available through * interface operator
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point) :: prodFact2     ! Result of th multiplication
 TYPE (point), INTENT(IN):: v  ! Inintial vector
 REAL(SGL), INTENT(IN):: f          ! Scalar to multiply
 prodFact2%x = v%x * f
 prodFact2%y = v%y * f
 prodFact2%z = v%z * f
END FUNCTION prodFact2

PURE FUNCTION negvec(v)
 !
 ! Purpose:
 ! To change sing to R3 single-precision vector
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE(point) :: negvec        ! Result: opposite sign vector
 TYPE(point), INTENT(IN) :: v ! Initial vector
 negvec = prodFact(-1.0, v)
END FUNCTION negvec

PURE FUNCTION restaVec (v1,v2)
 !
 ! Purpose:
 ! To substract two single precision R3 vectors.
 ! It should be only available through - interface operator
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point) :: restaVec            ! Result: Substracted vectors
 TYPE (point), INTENT(IN) :: v1, v2  ! Initial 2 vectors
 restaVec = v1 + (-v2)
END FUNCTION restaVec

PURE FUNCTION dot (v1,v2)
 !
 ! Purpose:
 ! To compute dot product between vectors ( 3R, real-single-presicion)
 !
 TYPE (point), INTENT(IN):: v1,v2   ! Initial 2 vec
 REAL(SGL) :: dot                   ! Result: Real value, dot product of v1.v2
 dot = v1%x * v2%x  +  v1%y * v2%y  +  v1%z * v2%z
END FUNCTION dot

PURE FUNCTION cross(v1,v2)
 !
 ! Purpose:
 ! To compute cross product between vectors ( 3R, real-single-presicion)
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point) :: cross               ! Result: Cross product between v1, v2
 TYPE (point), INTENT(IN) :: v1,v2   ! Initial vectors v1, v2
 cross%x = -v1%z*v2%y + v1%y*v2%z
 cross%y =  v1%z*v2%x - v1%x*v2%z
 cross%z = -v1%y*v2%x + v1%x*v2%y
END FUNCTION cross

PURE FUNCTION module (v)
 !
 ! Purpose:
 ! To compute module of a single-precision vector
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 REAL(SGL) :: module            ! Result: Module of the vector
 TYPE (point), INTENT(IN) :: v  ! Initial vector
 module = SQRT(v*v)
END FUNCTION module

PURE FUNCTION makeUnit (v)
 !
 ! Purpose:
 ! To make a R3 single precision unit unitary
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point):: makeUnit        ! Result: Unitary R3 vector
 TYPE (point), INTENT(IN) :: v  ! Initial vector
 makeUnit = (1/module(v))*v
END FUNCTION makeUnit

PURE FUNCTION cosang (v1,v2)
 !
 !  Purpose:
 !  To compute the cosinus of the angle between to vectors
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 REAL(SGL) :: cosang                   ! Result: cosinus of the angle
 TYPE (point), INTENT(IN) :: v1,v2     ! Inital passed vectors
 cosang = (v1*v2) /module(v1) /module(v2)
END FUNCTION cosang


SUBROUTINE topint (a,p)
 !
 ! Purpose:
 ! ???
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (pointInt), INTENT(OUT) :: a    ! Result: add what it does
 TYPE (point), INTENT(IN) :: p        ! Original point
 a%i=INT(p%x)
 a%j=INT(p%y)
 a%k=INT(p%z)
END SUBROUTINE topint

SUBROUTINE topreal (a,p)
 !
 ! Purpose:
 ! ???
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point), intent (OUT) :: a   ! Result: add meaning
 TYPE (pointInt), intent (IN) :: p ! Original integer point
 a%x=REAL(p%i)
 a%y=REAL(p%j)
 a%z=REAL(p%k)
END SUBROUTINE topreal

PURE FUNCTION restaint (p,ip)
 !
 ! Purpose:
 ! ???
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point) restaint               ! Result: add meaning
 TYPE (point), INTENT(IN) :: p       ! Original real-single-precision point
 TYPE (pointInt), INTENT(IN) :: ip   ! Original integer point
 restaint%x = p%x - ip%i
 restaint%y = p%y - ip%j
 restaint%z = p%z - ip%k
END FUNCTION restaint

PURE FUNCTION calcDist2 (p1,p2)
 !
 ! Purpose:
 ! To compute the squared distance between 2 points p1, p2
 !
  IMPLICIT NONE
 ! Data dictonary and variables declaration
  REAL(SGL) :: calcDist2             ! Result squared distance
  TYPE (point), INTENT(IN) :: p1,p2  ! Original points
  TYPE (point) :: v                  ! Auxiliar variable to store distance
  v = restaVec(p1,p2)
  calcDist2 = v * v
END FUNCTION calcDist2

PURE FUNCTION calcDist (p1,p2)
 !
 ! Purpose:
 ! To compute the distance between 2 points p1, p2
 !
  IMPLICIT NONE
 ! Data dictonary and variables declaration
  REAL :: calcDist                      ! Result: Distance ??units
  TYPE (point), INTENT(IN) :: p1,p2     ! Initial points p1, p2
  calcDist = SQRT (calcDist2 (p1,p2))
END FUNCTION calcDist

PURE FUNCTION rotx (p,cint,sint)
 !
 ! Purpose:
 ! To compute the rotation arround x cartesian axes
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point), INTENT(IN) :: p         ! Original point
 TYPE (point) :: rotx                  ! Result : original point rotated arround x axes
 REAL(SGL), INTENT(IN) :: cint         ! ? COS of projected axes defining rotation
 REAL(SGL), INTENT(IN) :: sint         ! ? SIN of projected axes defining rotation
 rotx%x = p%x
 rotx%y = p%y * cint  -  p%z * sint
 rotx%z = p%y * sint  +  p%z * cint
END FUNCTION rotx

PURE FUNCTION roty (p,cint,sint)
 !
 ! Purpose:
 ! To compute the rotation arround y cartesian axes
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point), INTENT(IN) :: p         ! Original point
 TYPE (point) :: roty                  ! Result : original point rotated arround y axes
 REAL(SGL), INTENT(IN) :: cint         ! ? COS of projected axes defining rotation
 REAL(SGL), INTENT(IN) :: sint         ! ? SIN of projected axes defining rotation
 roty%y = p%y
 roty%x = p%x * cint  -  p%z * sint
 roty%z = p%x * sint  +  p%z * cint
END FUNCTION roty

PURE FUNCTION rotz (p,cint,sint)
 !
 ! Purpose:
 ! To compute the rotation arround z cartesian axes
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point), INTENT(IN) :: p         ! Original point
 TYPE (point) :: rotz                  ! Result : original point rotated arround z axes
 REAL(SGL), INTENT(IN) :: cint         ! ? COS of projected axes defining rotation
 REAL(SGL), INTENT(IN) :: sint         ! ? SIN of projected axes defining rotation
 rotz%x = p%x * cint  -  p%y * sint
 rotz%y = p%x * sint  +  p%y * cint
 rotz%z = p%z
END FUNCTION rotz

PURE FUNCTION rotphi (p,cint,sint)
 !
 ! Purpose:
 ! To compute the rotation arround z cartesian axes
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point), INTENT(IN) :: p         ! Original point
 TYPE (point) :: rotphi                ! Result : original point rotated arround z axes
 REAL(SGL), INTENT(IN) :: cint         ! ? COS of projected axes defining rotation
 REAL(SGL), INTENT(IN) :: sint         ! ? SIN of projected axes defining rotation
 rotphi = rotz (p,cint,sint)
END FUNCTION rotphi

PURE FUNCTION rotthe (p,cint,sint)
 !
 ! Purpose:
 ! To compute the rotation arround XY cartesian axes
 !
 IMPLICIT NONE
 ! Data dictonary and variables declaration
 TYPE (point), INTENT(IN) :: p         ! Original point
 TYPE (point) :: rotthe                ! Result : original point rotated arround xy axes
 REAL(SGL), INTENT(IN) :: cint         ! ? COS of projected axes defining rotation
 REAL(SGL), INTENT(IN) :: sint         ! ? SIN of projected axes defining rotation
 REAL(SGL) :: rxy                      ! Auxiliary variable
 rxy = sqrt( p%x**2 + p%y**2 )
 IF (rxy /= 0._SGL) THEN
   rotthe%x =p%x * cint  +  p%z * p%x / rxy * sint
 ELSE
   rotthe%x = p%x * cint
 END IF
!! ATT ¿por que estan separados?
 IF ( rxy /= 0._SGL ) THEN
   rotthe%y = p%y * cint  +  p%z * p%y / rxy * sint
 ELSE
   rotthe%y = p%y * cint
 END IF
 rotthe%z = p%z * cint  -  rxy * sint
END FUNCTION rotthe

PURE FUNCTION rota (p, func, cint, sint)
 !
 ! Purpose:
 ! To select which rotation apply
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (point) :: rota                    ! Result: Rotated point
 TYPE (point), INTENT(IN) :: p           ! Original point
 CHARACTER(len=3), INTENT(IN) :: func    ! Which rotation function applies
 REAL(SGL), INTENT(IN) :: cint           ! Projection of rotation (COS angle)
 REAL(SGL), INTENT(IN) :: sint           ! Projection of rotation (SIN angle)
 CHARACTER(len=4)  :: aux_file           ! Hold temp stop mssg
SELECT CASE (func)
 CASE('x')
   rota = rotx(p,cint,sint)
 CASE('y')
   rota = roty(p,cint,sint)
 CASE('z')
   rota = rotz(p,cint,sint)
 CASE('phi')
   rota = rotphi(p,cint,sint)
 CASE('the')
   rota = rotthe(p,cint,sint)
 CASE DEFAULT
   ! This should make program stop, but because of being pure function is not as easy as STOP
   WRITE(aux_file,*) "STOP"
END SELECT
END FUNCTION rota

PURE FUNCTION rotaAng (p,func,a)
 !
 !  Purpose:
 !  To rotate a point an angle a arroud any axis
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (point) :: rotaAng                 ! Result: rotated point
 TYPE (point), INTENT(IN) :: p           ! Original point
 CHARACTER(len=3), INTENT(IN) :: func    ! which rottating function applies
 REAL(SGL), INTENT(IN) :: a              ! Angle to rotate
 REAL(SGL) :: cint                       ! AUX variable, COS of rotating angle
 REAL(SGL) :: sint                       ! AUX variable, SIN of rotating angle
 cint = cos (a)
 sint = sin (a)
 rotaAng = rota (p,func,cint,sint)
END FUNCTION rotaAng

PURE FUNCTION pol2car (p)
 !
 !  Purpose:
 !  To swich from polar coordinates to cartesian coordinates
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (point):: pol2car                 ! Result cartesian coordinates
 TYPE (pointPolar), INTENT(IN):: p      ! Original polar-point
 pol2car%x = p%r * sin(p%the) * cos(p%phi)
 pol2car%y = p%r * sin(p%the) * sin(p%phi)
 pol2car%z = p%r * cos(p%the)
END FUNCTION pol2car

PURE FUNCTION car2pol (p)
 !
 !  Purpose:
 !  To swich from cartesian coordinates to polar coordinates
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (point), INTENT(IN):: p            ! Original point
 TYPE (pointPolar) :: car2pol            ! Result: point in polar coordinates
 car2pol%r = module(p)
 car2pol%phi = atan2( p%y, p%x )
 car2pol%the = atan2( sqrt( p%x**2 + p%y**2 ), p%z )
END FUNCTION car2pol

PURE FUNCTION dist2pol (p1,p2)
 !
 !  Purpose:
 !  To compute distance between to polar points
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (pointPolar), INTENT(IN) :: p1,p2  ! Original points
 REAL(SGL) :: dist2pol                   ! Result: distance between polar numbers
 dist2pol = calcDist2 ( pol2car(p1) , pol2car(p2) )
END FUNCTION dist2pol

PURE FUNCTION rotquat (p1,pv,sint,cint)
 !
 !  Purpose:
 !  ???
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (point), INTENT(IN) :: p1, pv     ! Initial points
 TYPE (point) :: rotquat                ! Result, add meaning
 REAL(SGL), INTENT(IN) :: cint          ! Projection of rotation (COS angle)
 REAL(SGL), INTENT(IN) :: sint          ! Projection of rotation (SIN angle)
 REAL(SGL) :: w                         ! Aux. ?
 REAL(SGL) :: x                         ! Aux. ?
 REAL(SGL) :: y                         ! Aux. ?
 REAL(SGL) :: z                         ! Aux. ?
 x = pv%x * sint
 y = pv%y * sint
 z = pv%z * sint
 w = cint
rotquat%x = p1%x *(1-2 *(y**2 +z **2))+p1%y * 2* (x*y+w*z)+p1%z*2*(x*z+w*y)
rotquat%y = p1%x *2 *(x*y + w*z) +p1%y* (1-2*(x**2+z**2))+p1%z*2*(y*z-w*x)
rotquat%z = p1%x *2 *(x*z - w*y) +p1%y* 2* (y*z+w*x)+p1%z*(1-2*(x**2+y**2))
END FUNCTION rotquat

PURE FUNCTION rotEuler (p1,s1,c1,s2,c2,s3,c3)
 !
 !  Purpose:
 !  ???
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (point), INTENT(IN) :: p1         ! Original point
 TYPE (point) :: rotEuler               ! Result add meaning
 REAL(SGL), INTENT(IN) :: s1            ! ?
 REAL(SGL), INTENT(IN) :: c1            ! ?
 REAL(SGL), INTENT(IN) :: s2            ! ?
 REAL(SGL), INTENT(IN) :: c2            ! ?
 REAL(SGL), INTENT(IN) :: s3            ! ?
 REAL(SGL), INTENT(IN) :: c3            ! ?
 rotEuler%x = p1%x*c3*c2*c1 -p1%x*s3*s1 +p1%y*c3*c2*s1 +p1%y*s3*c1 -p1%z*c3*s2
 rotEuler%y =-p1%x*s3*c2*c1 -p1%x*c3*s1 -p1%y*s3*c2*s1 +p1%y*c3*c1 +p1%z*s3*s2
 rotEuler%z = p1%x*s2*c1                +p1%y*s2 *s1              +p1%z*c2
END FUNCTION rotEuler

PURE FUNCTION rotEuler1 (s1,c1,s2,c2,s3,c3)
 !
 !  Purpose:
 !  ???
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (point) :: rotEuler1              ! Result add meaning
 REAL(SGL), INTENT(IN) :: s1            ! ?
 REAL(SGL), INTENT(IN) :: c1            ! ?
 REAL(SGL), INTENT(IN) :: s2            ! ?
 REAL(SGL), INTENT(IN) :: c2            ! ?
 REAL(SGL), INTENT(IN) :: s3            ! ?
 REAL(SGL), INTENT(IN) :: c3            ! ?
 rotEuler1 = rotEuler(point(1.0,0.0,0.0),s1,c1,s2,c2,s3,c3)
END FUNCTION rotEuler1

PURE FUNCTION rotaEix (r,n0,s1,c1)
 !
 !  Purpose:
 !  ???
 !
 IMPLICIT NONE
 ! Data dictionary and variable declaration
 TYPE (point) rotaEix                   ! Result: add meaning
 TYPE (point), INTENT(IN) :: r          ! Initial ?
 TYPE (point), INTENT(IN) :: n0         ! Initial ?
 TYPE (point) :: n                      ! Aux ?
 REAL(SGL), INTENT(IN) :: s1            ! ?
 REAL(SGL), INTENT(IN) :: c1            ! ?
 n=makeUnit(n0)
 rotaEix=(c1 * r + s1 * cross(r,n)) + (dot(n,r)*(1 - c1) * n)
END FUNCTION rotaEix

END MODULE geometry
