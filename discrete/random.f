MODULE random
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
PUBLIC ran1
!
CONTAINS
   SUBROUTINE ran1(z, idum)
   IMPLICIT NONE
   ! Data dictionary and variable declaration
   INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)! Double precision real number definition processor independent
   INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6) ! Single Precision real number definition processor independent
   ! RANDOM NUMBER GENERATORS
   INTEGER :: iff =0                   !AUX random num generator
   INTEGER, PARAMETER :: m1 = 259200   !AUX random num generator
   INTEGER, PARAMETER :: m2 = 134456   !AUX random num generator
   INTEGER, PARAMETER :: m3 = 243000   !AUX random num generator
   INTEGER, PARAMETER :: ia1 = 7141    !AUX random num generator
   INTEGER, PARAMETER :: ia2 = 8121    !AUX random num generator
   INTEGER, PARAMETER :: ia3 = 4561    !AUX random num generator
   INTEGER, SAVE :: ix1                !AUX random num generator
   INTEGER, SAVE :: ix2                !AUX random num generator
   INTEGER, SAVE :: ix3                !AUX random num generator
   INTEGER, SAVE :: j                  !AUX random num generator
   INTEGER, PARAMETER :: ic1 = 54773   !AUX random num generator
   INTEGER, PARAMETER :: ic2 = 28411   !AUX random num generator
   INTEGER, PARAMETER :: ic3 = 51349   !AUX random num generator
   INTEGER, INTENT(INOUT) :: idum      !AUX random num generator
   REAL(DBL), PARAMETER :: rm1 =1./m1  !AUX random num generator
   REAL(DBL), PARAMETER :: rm2 =1./m2  !AUX random num generator
   REAL(DBL), DIMENSION(97), SAVE :: r                !AUX random num generator
   REAL(DBL), INTENT(OUT) :: z         !Result: pseudo-random numer
!
!
   IF (idum < 0 .or. iff == 0) THEN
      iff = 1
      ix1 = mod(ic1 - idum , m1)
      ix1 = mod(ia1 * ix1 + ic1, m1)
      ix2 = mod(ix1, m2)
      ix1 = mod(ia1 * ix1 + ic1, m1)
      ix3 = mod(ix1, m3)
      DO j = 1, 97
         ix1 = mod(ia1 * ix1 + ic1, m1)
         ix2 = mod(ia2 * ix2 + ic2, m2)
         r(j) = (float(ix1) + float(ix2) * rm2) * rm1
      END DO
   END IF
!
!   except when initializing, this is were we start.
!
   ix1 = mod ( ia1 * ix1 + ic1 , m1 )
   ix2 = mod ( ia2 * ix2 + ic2 , m2 )
   ix3 = mod ( ia3 * ix3 + ic3 , m3 )
   j = 1 + (97 * ix3) / m3
!
!   ¿es un error? ¿debería el programa parar?
!   if (j.gt.97.or.j.lt.1) write(0,*) ' AAAUUGGHH !!!'
!
   z = r(j)
   r(j) = (float(ix1) + float(ix2) * rm2) * rm1
!
   END SUBROUTINE ran1
 END MODULE random
