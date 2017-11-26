MODULE dims_utils
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
!PUBLIC  distanceCA, MCCheck, MCWeigth, &
!       readPDBtarg, readPDBtarg_CA,makeunit_3N, &
!       get_coef,transition_def,NMalign, readPDBref_CA, &
!       read_evec,get_trans_vec,decide_dims,get_xbeta
PUBLIC distanceCA, MCCheck, MCWeigth, readPDBtarg, makeunit_3N, &
       transition_def, get_trans_vec, decide_dims, get_xbeta, &
       transition_NM,MCCheck_NM,lambda,overlap,transition_NM_targ



!
! Data dictionary and variable declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)  ! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6)   ! Single Precision real number definition processor independent
REAL, PARAMETER :: PI=3.141592654 ! PI number definition
!
CONTAINS

subroutine get_xbeta(xbeta,aux_acc,error,rmsd_cutoff,lower_acceptance,acceptance,acc_tolerance)
IMPLICIT NONE
   REAL(SGL),INTENT(INOUT) :: xbeta
   REAL(DBL), INTENT(IN) :: error
   REAL(SGL), INTENT(IN) :: rmsd_cutoff
   REAL(SGL), INTENT(IN) :: aux_acc
   REAL(SGL), INTENT(IN) :: lower_acceptance
   REAL(SGL), INTENT(IN) :: acceptance
   REAL(SGL), INTENT(IN) ::acc_tolerance

      IF (error < rmsd_cutoff ) THEN
         IF (aux_acc  < lower_acceptance) THEN
            xbeta=xbeta*1.10
         ENDIF
      ELSE
         IF (aux_acc  < acceptance-acc_tolerance) THEN
            xbeta=xbeta*1.10
         ENDIF
      ENDIF
      !
      IF( aux_acc > acceptance+acc_tolerance ) THEN
         xbeta =xbeta * 0.9
      ENDIF

 end subroutine get_xbeta
!===============================================================================
! function distance(r,rtarg,disttarg,w,ica,natom,nres) result (score)
! use geometryDP
!   integer, intent(IN) :: natom, nres,ica(nres)
!   type(pointDP), intent(IN) :: r(natom),rtarg(natom)
!  ! real disttarg(nres,nres),w(nres,nres)
!   real disttarg(nres), w(nres)
!   real score, rij
!   integer i,j
!   score=0.
!  ! do i=1, natom
!   do i=1,nres!-1
!     ! do j=i+1,nres
!     !    rij = calcDistDP(r(ica(i)),r(ica(j)))
!     !    score=score+w(i,j)*(rij-disttarg(i,j))**2
!     ! enddo
!   ! ESTO FUNCIONA
!     rij = calcDistDP(r(ica(i)),rtarg(ica(i)))
!     score=score+(rij)*w(i)!-disttarg(i))*(rij-disttarg(i))
!   !   rij= calcDistDP(r(i),rtarg(i))
!   !   score=score+w(i)*rij
!   enddo
!  ! write(0,*) score, rij
!  ! stop
!   return
! end function distance
!===============================================================================
!===============================================================================
 function distanceCA(r,rtarg,w,ica,natom,nres) result (score)
 use geometryDP
   integer, intent(IN) :: natom, nres,ica(nres)
   type(pointDP), intent(IN) :: r(natom),rtarg(natom)
   real w(nres)
   real score, rij
   integer :: i
   score=0.
  do i=1,nres
   rij = calcDistDP(r(ica(i)),rtarg(ica(i)))
   score=score+(rij)*w(i)
  enddo
 end function distanceCA
!===============================================================================

 function MCCheck (error,seed, xbeta, scoreprev, score) result (rej)
 ! If TRUE it makes trajectory rewind
 use random
   integer, intent(INOUT) :: seed
   real, intent(IN) :: xbeta
   real(DBL), intent(in) :: score,error,scoreprev
   logical rej
   real sto,dscore
   real*8 fi

dscore=scoreprev-score
   call ran1(fi, seed)
 sto=exp(1.0/(xbeta)**2 * sign(1.,dscore)*(dscore)**2 /error**2)

    rej = (sto.lt.fi)

 end function MCCheck
!===============================================================================
SUBROUTINE MCWeigth(nres,natom,r, rtarg, ica,w)
use geometryDP
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nres
TYPE(pointDP),DIMENSION(natom), INTENT(IN) :: r
TYPE(pointDP),DIMENSION(natom), INTENT(IN) :: rtarg
INTEGER,DIMENSION(nres), INTENT(IN) :: ica
REAL(SGL) , INTENT(INOUT):: w(nres)
INTEGER :: i
integer, INTENT(IN):: natom
REAL :: temp_distance

 do i=1,nres
  temp_distance = calcDistDP(rtarg(ica(i)) ,r(ica(i)) )
  IF (temp_distance > 1.5) THEN
      w(i)=temp_distance
      ELSE
       w(i)=0
      END IF
 enddo

END SUBROUTINE MCWeigth

SUBROUTINE readPDBtarg(natom,unit_targ, atom, res, rtarg)
USE geometryDP
IMPLICIT NONE
INTEGER, INTENT(IN) :: natom
INTEGER, INTENT(IN) :: unit_targ
CHARACTER(len=4), DIMENSION(:), INTENT(IN) :: res,atom
TYPE(pointDP), INTENT(INOUT), DIMENSION(natom) :: rtarg
!
character(len=4) :: c1,c2,c3
integer i, j, k
      do i = 1,natom
         read (unit_targ, '(A4,2X,I5,2X,A3,1X,A3,3X,I3,4X,3F8.3)') c1,j,c2,c3,k, rtarg(i)
         if(c2 == 'CA'.and. c2.ne.atom(i).or.c3.ne.res(i))then
            write(*,*)'Cadenes diferents 1!',  c2, atom(i), "res ", c3,res(i)
            ! stop 1
         endif
      enddo
 !
 END SUBROUTINE readPDBtarg
!
!!----------------------------------------------------------
!SUBROUTINE readPDBref_CA(nres,natom, unit_targ,ica, atom, res,rtarg)
!USE geometryDP
!IMPLICIT NONE
!INTEGER, INTENT(IN) :: nres,natom
!INTEGER, INTENT(IN) :: unit_targ
!INTEGER, DIMENSION(:), INTENT(IN) :: ica
!CHARACTER(len=4), DIMENSION(:), INTENT(IN) :: res,atom
!TYPE(pointDP), INTENT(INOUT), DIMENSION(natom) :: rtarg
!! REAL(DBL), INTENT(INOUT), DIMENSION(3, nres) :: coord2CA
!!
!character(len=4) :: c1,c2,c3
!integer i, j, k
!      do i = 1,nres
!         read (unit_targ, '(A4,2X,I5,2X,A3,1X,A3,3X,I3,4X,3F8.3)') c1,j,c2,c3,k,rtarg(ica(i))
!      !
!     if(c2 =='CA') then
!            if( c2.ne.atom(ica(i)).or.c3.ne.res(ica(i)))then
!               write(*,*)'Cadenes diferents 2!',c2,atom(ica(i)),c3,res(ica(i))
!              stop 1
!             endif
!       end if
!      !
!      enddo
! !
! END SUBROUTINE readPDBref_CA
!!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!SUBROUTINE readPDBtarg_CA(nres, unit_targ,ica, atom, res, rtarg, coord2CA)
!USE geometryDP
!IMPLICIT NONE
!INTEGER, INTENT(IN) :: nres
!INTEGER, INTENT(IN) :: unit_targ
!INTEGER, DIMENSION(:), INTENT(IN) :: ica
!CHARACTER(len=4), DIMENSION(:), INTENT(IN) :: res,atom
!TYPE(pointDP), INTENT(INOUT), DIMENSION(nres) :: rtarg
!REAL(DBL), INTENT(INOUT), DIMENSION(3, nres) :: coord2CA
!!
!character(len=4) :: c1,c2,c3
!integer i, j, k
!      do i = 1,nres
!         read (unit_targ, '(A4,2X,I5,2X,A3,1X,A3,3X,I3,4X,3F8.3)') c1,j,c2,c3,k, rtarg(i)
!      	!
!     if(c2 =='CA') then
!       !  coord2(1:3,i)=rtarg(i)
!         coord2CA(1,i)=rtarg(i)%x
!         coord2CA(2,i)=rtarg(i)%y
!         coord2CA(3,i)=rtarg(i)%z
!            if( c2.ne.atom(ica(i)).or.c3.ne.res(ica(i)))then
!               write(*,*)'Cadenes diferents 3!',ica(i),c3,res(ica(i))
!               stop 1
!             endif
!       end if
!      !
!      enddo
! !
! END SUBROUTINE readPDBtarg_CA
!!-----------------------------------------------------------------------
! SUBROUTINE NMalign(nres , evec , U , center1, center2)
! use geometryDP
! use geometry
! INTEGER, INTENT(IN):: nres
!! INTEGER, DIMENSION(nres), INTENT(IN) :: ica
! real(DBL), dimension(3) :: coord,evecaux
! real(DBL), dimension(3,3), intent(in) :: U
! real(DBL), dimension(3), intent(in) :: center1, center2
! TYPE(point), intent(inout),dimension(nres) :: evec
!!
! INTEGER i
!
!    do i=1, nres
!      coord(1)=evec(i)%x ; coord(2)=evec(i)%y ; coord(3)=evec(i)%z
!      evecaux(1:3)=matmul(U, coord(1:3)-center2)+center1! ojo
!      evec(i)%x=evecaux(1)
!      evec(i)%y=evecaux(2)
!      evec(i)%z=evecaux(3)
!    end do
!
!
!END SUBROUTINE NMalign
!!-----------------------------------------------------------------------

SUBROUTINE transition_def(nres,r,rtarg,ica,delta_rab)
USE geometryDP
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: nres
INTEGER, INTENT(IN), DIMENSION(nres) :: ica
TYPE(pointDP), INTENT(IN), DIMENSION(:) :: r
TYPE(pointDP), INTENT(IN), DIMENSION(:) :: rtarg
TYPE(pointDP), INTENT(OUT), DIMENSION(nres) :: delta_rab
INTEGER :: i
!
! Calcular el 3N vector de las diferentes estructuras
!
do i=1, nres
delta_rab(i)=rtarg(ica(i))-r(ica(i) )
enddo
call makeunit_3N(nres, delta_rab)
!
!
END SUBROUTINE transition_def
!----------------------------------------------------------

PURE FUNCTION get_coef(nres, vec_3N, evec) RESULT (coef)
USE geometry
USE geometryDP
IMPLICIT NONE
INTEGER, INTENT(IN) :: nres
TYPE(pointDP), INTENT(IN), DIMENSION(nres):: vec_3N
TYPE(pointDP), INTENT(IN), DIMENSION(nres):: evec
REAL(DBL) :: coef
INTEGER :: i

coef=0._DBL
DO i=1, nres
coef= coef+dotDP(vec_3N(i),evec(i))
END DO

END FUNCTION get_coef
!----------------------------------------------------------
!
SUBROUTINE makeunit_3NSP(N, vec_3N)
USE geometry
USE geometryDP
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
TYPE(point), INTENT(INOUT), DIMENSION(N) :: vec_3N
REAL(DBL) :: norma
INTEGER ::i
!
norma=0._DBL
DO i=1,N
  norma=norma+vec_3N(i)%x**2 +vec_3N(i)%y**2 + vec_3N(i)%z**2
ENDDO
norma=SQRT(norma)
DO i=1,N
  vec_3N(i)%x=vec_3N(i)%x/norma
  vec_3N(i)%y=vec_3N(i)%y/norma
  vec_3N(i)%z=vec_3N(i)%z/norma
END DO
END SUBROUTINE makeunit_3NSP

!--------------------------------------------------------------
!----------------------------------------------------------
!
SUBROUTINE makeunit_3N(N, vec_3N)
USE geometry
USE geometryDP
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
TYPE(pointDP), INTENT(INOUT), DIMENSION(N) :: vec_3N
REAL(DBL) :: norma
INTEGER ::i
!
norma=0._DBL
DO i=1,N
  norma=norma+vec_3N(i)%x**2 +vec_3N(i)%y**2 + vec_3N(i)%z**2
ENDDO
norma=SQRT(norma)
DO i=1,N
  vec_3N(i)%x=vec_3N(i)%x/norma
  vec_3N(i)%y=vec_3N(i)%y/norma
  vec_3N(i)%z=vec_3N(i)%z/norma
END DO
END SUBROUTINE makeunit_3N

!--------------------------------------------------------------

!SUBROUTINE  read_evec(nres,evec)
!use geometry
! INTEGER, INTENT(IN) :: unit_o
!INTEGER ioerr
!INTEGER ein_num
!INTEGER,INTENT(IN) :: nres
!TYPE(point), INTENT(INOUT), DIMENSION(nres):: evec
!INTEGER i
!real aux
!open(unit=53,file='/Users/psfriso/workspace/discreteNM/evec.dat', status='OLD')
!read(53,*) ein_num, aux
!  write(unit_o,*)ein_num,aux
!DO i=1,nres
!   read(53, *, iostat=ioerr) evec(i)%x,evec(i)%y,evec(i)%z
!   IF(ioerr /=0 ) then
!    write(0,*) "hola", ioerr
!    STOP 2
!    END IF
!END DO
!END SUBROUTINE read_evec
!

!----------------------------------------------------------
 SUBROUTINE get_trans_vec(nres,nevecs,evec,unit_o,delta_rab,acc_coef,trans_vec)
USE ls_rmsd
USE geometry
USE geometryDP
IMPLICIT NONE
!
INTEGER, INTENT(IN) ::nres!,natom
INTEGER, INTENT(IN) :: unit_o
TYPE(pointDP), DIMENSION(nres), INTENT(IN) :: delta_rab
!TYPE(pointDP), DIMENSION(natom), INTENT(IN) :: r
! INTEGER, DIMENSION(nres), INTENT(IN) :: ica
!CHARACTER(len=4), INTENT(IN) :: atom(natom), res(natom)
TYPE(pointDP), DIMENSION(nres), INTENT(OUT) :: trans_vec
INTEGER, INTENT(IN) :: nevecs
type(pointDP), dimension(nevecs, nres), INTENT(IN) :: evec
real(DBL), dimension(nevecs) :: coef
real(sgl),INTENT(OUT) ::acc_coef
integer i
integer j
!
   !==========================================================
   ! OBTENER VECTOR DE TRANSICION
   ! ENTRA nres, unit_evecs, unit_o, delta_rab, ica, atom, res
   ! SALE trans_vec
   !==========================================================


  DO i=1,nevecs
 !   call read_evec(nres,unit_o,evec(i,1:nres))
    coef(i)=get_coef(nres, delta_rab, evec(i,1:nres))
!    write(unit_o, '("coef ", f6.3)')coef(i)
  END DO
!
  WHERE (coef**2 < 0.0225)
  coef=0._SGL!
  END WHERE
!



!do i =1, nres
!write(*, '(3f10.5)') evec(1,i)
!enddo
acc_coef=0.0
  DO i=1,nevecs
    acc_coef=acc_coef+coef(i)**2
  ENDDO
!  write(unit_o, '("accumulate coef ", f10.6)')acc_coef
  trans_vec=pointDP(0._DBL,0._DBL,0._DBL)
  DO i=1,nres
    DO j=1,nevecs
      trans_vec(i)%x=trans_vec(i)%x+coef(j)*evec(j,i)%x
      trans_vec(i)%y=trans_vec(i)%y+coef(j)*evec(j,i)%y
      trans_vec(i)%z=trans_vec(i)%z+coef(j)*evec(j,i)%z
    END DO
  END DO
!
  call makeunit_3N(nres, trans_vec)
!
END SUBROUTINE get_trans_vec
!==========================================================

SUBROUTINE decide_dims(error,rmsdlim,unit_o)
IMPLICIT NONE
REAL(DBL), INTENT(IN) :: error
INTEGER, INTENT(IN) :: unit_o
REAL(SGL),INTENT(IN) :: rmsdlim
! This Block check that dims makes sense
    IF( error < rmsdlim ) then
       write(unit_o, '("Not appropiate resolution for DMD-DIMS dynamics, ",f6.2 ,1x,"<",1x,f6.2,/,"Attention !")')error,rmsdlim
     !   STOP 20
    endif
 !   write(unit_o,'("Initial rmsd CA",f8.2,/)')  error
END SUBROUTINE decide_dims

SUBROUTINE transition_NM(natom,nres,ica,unit_o,r,rtarg,acc_coef,v_trans)
use anm_laura
use geometry
use geometryDP
use paramSet
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: nres, natom
INTEGER, INTENT(IN), DIMENSION(nres) :: ica
INTEGER, INTENT(IN) :: unit_o
TYPE(pointDP), INTENT(IN), DIMENSION(natom) :: r
TYPE(pointDP), INTENT(IN), DIMENSION(natom) :: rtarg
TYPE(pointDP),INTENT(OUT), DIMENSION(nres) :: v_trans
REAL(SGL), INTENT(OUT) :: acc_coef
!
type(pointDP), dimension(:), allocatable :: delta_rab
TYPE(pointDP), DIMENSION(:,:), ALLOCATABLE ::evec
INTEGER :: i,ioerr
! REAL(SGL) :: aux_dot
!
allocate( evec(nevecs, nres), stat=ioerr )
call errorAllocmem(ioerr, 'evecs' )
allocate(delta_rab(nres), stat=ioerr)
call errorAllocmem(ioerr, 'delta coordinates')
!==========================================================
! OBTENER VECTOR DE TRANSICION (entra  r, rtarg : sale v_trans)
!
DO i=1,nevecs
evec(i,:)=pointDP(0._DBL,0._DBL,0._DBL)
ENDDO
!
CALL ANM(nres, r(ica(:)),nevecs, evec)
!
delta_rab=pointDP(0._DBL,0._DBL,0._DBL)
call transition_def(nres,r,rtarg,ica,delta_rab)
!
v_trans=pointDP(0._DBL,0._DBL,0._DBL)
call get_trans_vec(nres,nevecs,evec,unit_o,delta_rab,acc_coef,v_trans)

!aux_dot=0._SGL
!do i=1, nres
!   aux_dot= aux_dot+dotDP(delta_rab(i), v_trans(i))
!end do
! Write(unit_o, '("Overlap Transition Vector and transition is ", f7.5)')sqrt(acc_coef)
!!
!===============================================================================
END SUBROUTINE transition_NM

SUBROUTINE transition_NM_targ(natom,nres,ica,unit_o,r,rtarg,acc_coef,v_trans)
use anm_laura
use geometry
use geometryDP
use paramSet
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: nres, natom
INTEGER, INTENT(IN), DIMENSION(nres) :: ica
INTEGER, INTENT(IN) :: unit_o
TYPE(pointDP), INTENT(IN), DIMENSION(natom) :: r
TYPE(pointDP), INTENT(IN), DIMENSION(natom) :: rtarg
TYPE(pointDP),INTENT(OUT), DIMENSION(nres) :: v_trans
REAL(SGL), INTENT(OUT) :: acc_coef
!
type(pointDP), dimension(:), allocatable :: delta_rab
TYPE(pointDP), DIMENSION(:,:), ALLOCATABLE ::evec
INTEGER :: i,ioerr
! REAL(SGL) :: aux_dot
!
allocate( evec(nevecs, nres), stat=ioerr )
call errorAllocmem(ioerr, 'evecs' )
allocate(delta_rab(nres), stat=ioerr)
call errorAllocmem(ioerr, 'delta coordinates')
!==========================================================
! OBTENER VECTOR DE TRANSICION (entra  r, rtarg : sale v_trans)
!
DO i=1,nevecs
evec(i,:)=pointDP(0._DBL,0._DBL,0._DBL)
ENDDO
!
CALL ANM(nres, rtarg(ica(:)),nevecs, evec)
!
delta_rab=pointDP(0._DBL,0._DBL,0._DBL)
call transition_def(nres,r,rtarg,ica,delta_rab)
!
v_trans=pointDP(0._DBL,0._DBL,0._DBL)
call get_trans_vec(nres,nevecs,evec,unit_o,delta_rab,acc_coef,v_trans)

!aux_dot=0._SGL
!do i=1, nres
!   aux_dot= aux_dot+dotDP(delta_rab(i), v_trans(i))
!end do
! Write(unit_o, '("Overlap Transition Vector and transition is ", f7.5)')sqrt(acc_coef)
!!
!===============================================================================
END SUBROUTINE transition_NM_targ


!===============================================================================
function MCCheck_NM(seed,NMcut,proj) result (wto)
 ! If TRUE it makes trajectory rewind
 use random
   integer, intent(INOUT) :: seed
   real :: beta=10.
   real(SGL), intent(in) :: NMCut,proj
   logical wto
   real sto
   real*8 fi

   call ran1(fi, seed)
   sto=exp(beta * sign(1.,proj-NMcut) *(NMCut-proj)**2)

    wto = .not.(sto.lt.fi)
! write(*,*) "MC NM", STO, fi, wto, NMcut-proj,beta
 end function MCCheck_NM
!===============================================================================

FUNCTION lambda(error0,error) result (lmb)
use paramSet
!
IMPLICIT NONE
!
REAL(DBL) , INTENT(IN) :: error,error0
REAL(SGL) :: lmb
REAL(SGL) :: aux
!
aux= exp ( (DBLE(10.* RMSdlim) - 20._DBL*error )/error0)
lmb = 1.0/( 1.0+ 22026. *aux)
!
END FUNCTION lambda
! =========================================================
FUNCTION overlap(natom,nres,ica,r,rtarg,v_trans) RESULT (dot_prd)
use geometry
use geometryDP
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: nres,natom
INTEGER, INTENT(IN), DIMENSION(nres) :: ica
TYPE(pointDP), INTENT(IN), DIMENSION(natom) :: r
TYPE(pointDP), INTENT(IN), DIMENSION(natom) :: rtarg
TYPE(pointDP),INTENT(IN), DIMENSION(nres) :: v_trans
REAL(SGL) :: dot_prd
!
type(pointDP), dimension(nres) :: delta_rab
INTEGER :: i

delta_rab=pointDP(0._DBL,0._DBL,0._DBL)
call transition_def(nres,r,rtarg,ica,delta_rab)
!
dot_prd=0._SGL
do i=1, nres
   dot_prd= dot_prd+dotDP(delta_rab(i), v_trans(i))
end do
!
END FUNCTION overlap
!==========================================================
END MODULE dims_utils