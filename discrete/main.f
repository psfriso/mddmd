!===========================================================================
!  DIMS 0.1
!
!  Discrete Molecular Dynamics, DIMS. Based on Discrete 0.2.3
!
!   v. 0.1a    First public release.
!      Read Topology and coordinates prepared by auxiliary program
!        Restart available
!
!==========================================================================

program discrete
   use commLine
   use geometry
   use geometryDP
   use stepPotentials
   use intList
   use paramSet
   use random
   use energy
   use colisio
   use ls_rmsd
   use dims_utils
   use ANM_laura
   use least_sq_fit
!
!
! CommLine & I/O
!
   integer, parameter :: NFILES = 8

   integer unit_i, unit_o, unit_top, unit_r, unit_ener, unit_traj, unit_rst, unit_targ!,unit_toptarg

   type(commLineOption) :: files(NFILES) = (/&
   commLineOption("-i",    "param",          "formatted",   "old",     "Settings"),&
   commLineOption("-top",  "topology",       "unformatted", "old",     "Topology"),&
   commLineOption("-r",    "coordinates",    "unformatted", "old",     "Initial coordinates"),&
   commLineOption("-ener", "energy",         "formatted",   "unknown", "Energies"),&
   commLineOption("-traj", "trajectory.pdb", "formatted",   "unknown", "Trajectory (PDB)"),&
   commLineOption("-o",    "log",            "formatted",   "unknown", "Calculation Log"),&
   commLineOption("-rst",  "restart",        "unformatted", "unknown", "Restart coordinates"), &
   commLineOption("-targ",  "targpdb",       "formatted",   "unknown", "Target (PDB)")/)

! Coordinates and velocities
   type(pointDP), allocatable :: r(:)
   type(pointDP), allocatable :: rprev(:)
   type(pointDP), allocatable :: rtarg(:)
   type(point), allocatable :: v(:), rsp(:) ! Single precision version of r for binary I/O
   real, allocatable :: distat2(:,:),distat2targ(:,:), w(:)
   type(point)  rcm ! C of M

! Step Potentials
   type (stepPotInt), allocatable :: stepPts(:,:)

! Covalent bonds
   real, allocatable :: rb(:,:),drb(:,:)

! Structure
   character(len=4), allocatable :: atom(:), res(:), atp(:)
   integer, allocatable :: in(:),ica(:),ico(:), rnum(:)
   logical, allocatable :: istruct(:,:)                   ! true for SSec based interactions
   integer natom, nres, ncovpairs, nhbs, nhelix, nbeta
   integer, allocatable :: cov(:,:), hbs(:,:), helix(:,:), beta(:,:)

! ! Structure TARGET
!   character(len=4), allocatable :: atom2(:), res2(:), atp2(:)
!   integer, allocatable :: in2(:),ica2(:),ico2(:), rnum2(:)
!   logical, allocatable :: istructtarg(:,:)               ! true for SSec based interactions
!   integer natom2, nres2, ncovpairs2, nhbs2, nhelix2, nbeta2
!   integer, allocatable :: cov2(:,:), hbs2(:,:), helix2(:,:), beta2(:,:)
!   real :: ax1,ax2,ax4,ax5,ax6,ax7,ax8
!

! Potentials & energies
   real, allocatable :: evdw(:),rvdw(:),qq(:),gfree(:),vol(:),xlamb(:), xm(:), rhc(:), xsum(:,:)
   real potlk, rvdwij, xmassa, esolv, ecoul, ekin, epotfis, epotgo
   REAL(SGL) :: ekin0=0._SGL
   REAL(DBL) :: score=1._SGL
   REAL(DBL) :: scoreprev

! Interaction pairs
   integer bpairs
   integer, allocatable :: blist(:,:)
   type(intpList), allocatable :: nblist(:)

! Collisions
   integer :: ierr
   real*8, allocatable :: tpart(:)
   integer, allocatable :: ipart(:)
   integer :: mem1
   integer :: mem2
   integer :: ibloc=1
   integer :: iev
   integer :: npair1
   integer :: npair2
   integer :: idec=0
   integer :: iacc=0

! Time
   real*8 tacact, taccorr, tacrect, tacum, tevent, tevent1, temps, tactot
   real*8 tinit, tsetup, tfin

!
   integer i,j,k, ioerr,auxiliar
!
!
   real, parameter :: A = 1.e-10 ! M to Angs
! RMSd & superposing variables
 real(DBL), dimension(3,3) :: U               ! Rotation Matrix
 real(DBL), dimension(3) :: center1, center2  ! Center of initial and target structures
 real(DBL) :: error=10E10, error0             ! RMSd values RMSd(t), RMSd(t=0)
 integer :: accepted                          ! Accepted MC traj
 integer :: NonAccepted                       ! Wrong Name, Total proposed trajtories
 real :: aux_acc                              ! ratio of acceptance

! NORMAL MODE ANALISIS
! type(pointDP), dimension(:), allocatable :: delta_rab  ! TO MODULE DIMS
type(pointDP), dimension(:), allocatable :: trans_vec    ! 3N transition vector
type(pointDP), dimension(:), allocatable :: diff_frames  ! Difference between frames between NM rectification step
real :: aux_dotBIG                                       ! Internal product diff_frames and NM trans vec
CHARACTER(len=54), dimension(:,:), allocatable :: buffer_frames ! Buffer to keep frames until they are accepted
CHARACTER(len=120), dimension(:), allocatable :: buffer_out     ! Buffer to keep log until frames are accepted
CHARACTER(len=120), dimension(:), allocatable :: buffer_ener    ! Buffer to energy log frames until frames are accepted
INTEGER :: wiblocs                       ! Already writen frames
type(pointDP), dimension(:), allocatable :: r_prevNM ! Last NM accepted coordinates
INTEGER :: buffer_size                    ! Number of frames to store
REAL(DBL) :: time_prevNM,tpedro_prevNM    ! Total NM accepted time ; Xbeta last rect time
REAL(SGL) :: score_prevNM, xbeta_prevNM   ! Last accepted NM score , xbeta
!TYPE(pointDP), DIMENSION(:,:), ALLOCATABLE ::evec
LOGICAL:: NMcontrol               ! To activate NM rectification schema
integer:: nm_rejects              ! Number of NN rejection
REAL(SGL) :: acc_coef             ! Overlap between coordinates difference and NM-based transition vector
REAL(DBL) :: tpedro               ! Time to next XBETA rectification
REAL(DBL) :: tpedroNM             ! Time to next NM rectification
REAL(DBL) :: tpedroNM_prevNM      ! Last accepted time to next NM rectification
LOGICAL :: firstNM                ! First NM cycle exception active
REAL(SGL) :: solap                ! overlap between remaining transition and NM
LOGICAL :: NMtargOnce             !
LOGICAL :: stopping=.FALSE.       ! Whether to finish simulation or not
!
! Self-STOP
REAL(SGL), ALLOCATABLE, DIMENSION(:,:) :: error_vector   ! rmsd vs time vector
REAL(SGL) :: slope=0._SGL                                ! rmsd vs time slope
INTEGER :: error_vector_size=0                           ! number of rmsd point to use in linear regression
INTEGER :: errcont=0                                     ! index to write new rmsd value

   call cpu_time(tinit)
   call readCommLine(files, unit_i, unit_o, NFILES)
   call readInputParamSet(unit_i)
!
   call programHeader(unit_o)
   call printFileSummary(files, unit_o)
   call writeInputParamSet(unit_o)
!
   write (unit_o, '(/," Reading system from topology file ",/)')
! readTopology HO mantenim aqui per evitar problemes amb allocate
      unit_top = openFn(files, '-top')
      read(unit_top) nres, natom
      write (unit_o, '(" Residues:   ",i4," Atoms:  ",i6)') nres, natom
      allocate (in(nres), ica(nres), ico(nres), stat=ioerr)
      call errorAllocmem (ioerr, 'NRes')
      allocate (atom(natom), rnum(natom), res(natom), atp(natom), qq(natom),  gfree(natom), &
                vol(natom), evdw(natom), rvdw(natom), rhc(natom), xm(natom), xlamb(natom), stat=ioerr)
      call errorAllocmem (ioerr, 'NAtom')
!
      read (unit_top) (in(i), ica(i), ico(i), i=1,nres)

!
      read (unit_top) (atom(i), rnum(i), res(i), atp(i), i=1,natom)
!
      do i=1,natom
         read (unit_top) qq(i), gfree(i), vol(i), evdw(i), rvdw(i), rhc(i), xm(i)
      enddo
      xmassa = sum(xm)
      xlamb = 3.5
!
      read (unit_top) ncovpairs, nhbs, j
      write (unit_o, '(" Cov. Pairs:",i5," H. Bonds: ",i4)') ncovpairs, nhbs
!
      allocate (cov(ncovpairs,2), hbs(nhbs,2), stat=ioerr)
      call errorAllocmem (ioerr, 'Cov & Hbs')
!
      read (unit_top) (cov(i,1), cov(i,2), i=1,ncovpairs)
      read (unit_top) (hbs(i,1), hbs(i,2), i=1,nhbs)
      read (unit_top) nhelix
      write (unit_o, '(" Num. Helix:       ",i4)') nhelix
!
      allocate (helix(nhelix,2), stat=ioerr)
      call errorAllocmem (ioerr, 'Helix')
      read (unit_top) (helix(i,1), helix(i,2), i=1,nhelix)
!
      read (unit_top) nbeta
      write (unit_o, '(" Num. Beta Strands:",i4)') nbeta
!
      allocate (beta(nbeta,2), stat=ioerr)
      call errorAllocmem (ioerr, 'Beta')
!
      read (unit_top) (beta(i,1), beta(i,2), i=1,nbeta)
      close (unit_top)
! end readTopology
!

! readCoordinates
      unit_r = openFn(files, '-r')
      read (unit_r) j
      if (j.ne.natom) then
         write (0,*) "ERROR: coordinates and topology files do not match"
         stop 1
      endif
      allocate (r(natom), rsp(natom), rprev(natom), rtarg(natom), distat2(natom,natom), distat2targ(natom,natom),stat=ioerr)
      call errorAllocmem (ioerr, 'Coordinates')
      read (unit_r) (rsp(i)%x, rsp(i)%y, rsp(i)%z, i=1,natom)
      close (unit_r)
! end readCoordinates

! NORMAL MODE

allocate(trans_vec(nres), stat=ioerr)
call errorAllocmem(ioerr, 'Transition Vector')
allocate(diff_frames(nres), stat=ioerr)
call errorAllocmem(ioerr, 'diff frames')
allocate(r_prevNM(natom), stat=ioerr)
call errorAllocmem(ioerr, 'r prev NM')



   write (unit_o, '(" Initial coordinates read in ",/)')
!
   do i = 1,natom
      r(i) = SPtoDP(rsp(i))
   enddo
!
! Centrer CM at (0,0,0)
! Duplicate original coordinates (FOR DIMS ONLY) !OJO
   rcm = calcCM (natom, r, xm)
   do i = 1,natom
      r(i) = r(i) - SPtoDP(rcm)
      rprev(i) = r(i)
   enddo

!
!
   do j = 2, natom
   do i = 1, j-1
      distat2(i,j) = calcDist2DP(r(i), r(j))
      distat2(j,i) = distat2(i,j)
   enddo
   enddo
!
   allocate (v(natom), rb(natom,natom), drb(natom,natom), istruct(natom,natom),&
               stepPts(natom,natom), &
       xsum(natom,natom), stat=ioerr)!istructtarg(natom,natom),
   call errorAllocmem(ioerr, 'Setup/Natom')
!


timeMC_NM = 1.5*ACCEPTANCE*timeMC_NM*log(REAL(nres)/20.0)
time_in_eq=timeMC_NM
error_vector_size=FLOOR(1000*time_in_eq/TSNAP)
buffer_size=CEILING(1000*timeMC_NM/TSNAP)+1
allocate( error_vector(error_vector_size,2 ) , stat=ioerr)
call errorAllocmem(ioerr, 'error vector')
error_vector=0.

! Posem temps en segons

   TSNAP = TSNAP * 1.e-15
   TRECT = TRECT * 1.e-15
   TCORR = TCORR * 1.e-15
   TACT = TACT * 1.e-15
   xbeta_update=xbeta_update *1.e-15
!
! ============================================================================
! This is a DIMS dynamics  with complete PDB reference structure and
! normal mode analysis
! NOW read coordinates of aim structure
! OPEN unit
!
! NM ANALISIS
!
 unit_targ = openFn(files, '-targ')
 call readPDBtarg(natom, unit_targ, atom, res, rtarg)
 close(unit_targ)
!
allocate (w(nres),stat=ioerr)
call errorAllocmem(ioerr, 'MC weigths')
! OPEN UNITS
unit_traj = openFn(files, '-traj')
! Compute initial RMSd all atom
!  call rmsd(natom,r,rtarg,0, U,  center2, center1, error,.FALSE.)
!  write(unit_o,'("Initial rmsd all atom",f8.2,/)')  error
! Compute initial RMSd CA, this is the used one
   call rmsd(nres,  r(ica(:)),rtarg(ica(:)), 1, U, center2, center1, error0,.FALSE.)
! Use rmsd alignment to superpose structures ( r over rtarg)
   IF(Superpose) call pdbsuperpos(natom, U , center1, center2,r)
!   call writeSnapShot (unit_traj, 0, rtarg, atom, res, rnum, natom)

!==========================================================
! DECIDE DIMS
! call decide_dims(error0,rmsdlim,unit_o) ! Needs correction, rmsdlim changed
!==========================================================
call writeSnapShot (unit_traj, 0, r, atom, res, rnum, natom)!
call MCWeigth(nres,natom,r, rtarg, ica,w)

!!==========================================================
!! OBTENER VECTOR DE TRANSICION (entra  r, rtarg : sale trans_vec)
!!
trans_vec=pointDP(0._DBL,0._DBL,0._DBL)
CALL transition_NM(natom,nres,ica,unit_o,r,rtarg,acc_coef,trans_vec)
! Write(unit_o, '("Overlap Transition Vector and transition is ", f7.5)')sqrt(acc_coef)
!
!==========================================================

scoreprev=1.E10
iacc=0
idec=0
!write (unit_o, '(" Initial scoring: ",f10.2,/)') score0
allocate( buffer_frames(buffer_size,natom+2), stat=ioerr)
call errorAllocmem(ioerr, 'frames out')
allocate( buffer_out(buffer_size), stat=ioerr)
call errorAllocmem(ioerr, 'Info out')
allocate( buffer_ener(buffer_size), stat=ioerr)
call errorAllocmem(ioerr, 'Energy Out')
buffer_frames=""
buffer_out=""
buffer_ener=""

!   write (unit_o, '(/," Reading system from topology file Target",/)')
!! Leer la topologia final. Aqui es un desastre debe estar en un unico archivo de topologia
!! que transmita la informacion de la interseccion de restraints.
!! readTopology HO mantenim aqui per evitar problemes amb allocate
!      unit_toptarg = openFn(files, '-toptarg')
!      read(unit_toptarg) nres2, natom2
!      allocate (in2(nres), ica2(nres), ico2(nres), stat=ioerr)
!      call errorAllocmem (ioerr, 'NRes')
!      allocate (atom2(natom), rnum2(natom), res2(natom), atp2(natom), stat=ioerr)
!      call errorAllocmem (ioerr, 'NAtom')
!!
!      read (unit_toptarg) (in2(i), ica2(i), ico2(i), i=1,nres)
!
!!
!      read (unit_toptarg) (atom2(i), rnum2(i), res2(i), atp2(i), i=1,natom)
!!
!      do i=1,natom
!         read (unit_toptarg) ax1, ax2, ax4, ax4, ax5, ax6, ax7
!      enddo
!
!      read (unit_toptarg) ncovpairs2, nhbs2, ax8
!!
!      allocate (cov2(ncovpairs2,2), hbs2(nhbs2,2), stat=ioerr)
!      call errorAllocmem (ioerr, 'Cov & Hbs')
!!
!      read (unit_toptarg) (cov2(i,1), cov2(i,2), i=1,ncovpairs2)
!      read (unit_toptarg) (hbs2(i,1), hbs2(i,2), i=1,nhbs2)
!      read (unit_toptarg) nhelix2
!!      write (unit_o, '(" Num. Helix:       ",i4)') nhelix
!!
!      allocate (helix2(nhelix2,2), stat=ioerr)
!      call errorAllocmem (ioerr, 'Helix')
!      read (unit_toptarg) (helix2(i,1), helix2(i,2), i=1,nhelix2)
!!
!      read (unit_toptarg) nbeta2
!!      write (unit_o, '(" Num. Beta Strands:",i4)') nbeta
!!
!      allocate (beta2(nbeta2,2), stat=ioerr)
!      call errorAllocmem (ioerr, 'Beta')
!!
!      read (unit_toptarg) (beta2(i,1), beta2(i,2), i=1,nbeta2)
!      close (unit_toptarg)
!! end readTopology
!---------------------------------------------------------------------------

   do j = 2, natom
   do i = 1, j-1
      distat2targ(i,j) = calcDist2DP(rtarg(i), rtarg(j))
      distat2targ(j,i) = distat2targ(i,j)
   enddo
   enddo


   rb = 0.
   drb = 0.
   istruct = .false.
!   istructtarg=.false.
   stepPts%active = .false.
!
   do i = 1,ncovpairs
      rb(cov(i,1),cov(i,2)) = sqrt(distat2(cov(i,1),cov(i,2)))
      drb(cov(i,1),cov(i,2)) = SIGMA*rb(cov(i,1),cov(i,2))
   enddo

!! TARGET
!
!   do i = 1,nhbs2
!      if (atom2(hbs2(i,1)).eq.'N'.or.atom2(hbs2(i,2)).eq.'N') then
!         rb (hbs2(i,1),hbs2(i,2)) = 0.5 * (rnomax + rnomin)
!         drb(hbs2(i,1),hbs2(i,2)) = 1E10
!      else
!         rb (hbs2(i,1),hbs2(i,2)) = 0.5 * (rcomax + rcomin)
!         drb(hbs2(i,1),hbs2(i,2)) = 1E10
!      endif
!   enddo

! REAL
   do i = 1,nhbs
 !    IF( rb(hbs(i,1),hbs(i,2)) > 1e-20 ) THEN
         if (atom(hbs(i,1)).eq.'N'.or.atom(hbs(i,2)).eq.'N') then
            rb (hbs(i,1),hbs(i,2)) = 0.5 * (rnomax + rnomin)
            drb(hbs(i,1),hbs(i,2)) = 0.5 * (rnomax - rnomin)
         else
            rb (hbs(i,1),hbs(i,2)) = 0.5 * (rcomax + rcomin)
           drb(hbs(i,1),hbs(i,2)) = 0.5 * (rcomax - rcomin)
!     ENDIF
      endif
   enddo



! !  TARGET
!      do i = 1,ncovpairs
!      rb(cov(i,1),cov(i,2)) = sqrt(distat2(cov(i,1),cov(i,2)))
!      drb(cov(i,1),cov(i,2)) = SIGMA*rb(cov(i,1),cov(i,2))
!   enddo
!   do i = 1,nhbs
!      if (atom(hbs(i,1)).eq.'N'.or.atom(hbs(i,2)).eq.'N') then
!         rb (hbs(i,1),hbs(i,2)) = 0.5 * (rnomax + rnomin)
!         drb(hbs(i,1),hbs(i,2)) = 0.5 * (rnomax - rnomin)
!      else
!         rb (hbs(i,1),hbs(i,2)) = 0.5 * (rcomax + rcomin)
!         drb(hbs(i,1),hbs(i,2)) = 0.5 * (rcomax - rcomin)
!      endif
!   enddo


! fixa els angles diedrics dels residus que formen part d'un element d'estructura secundaria
   if(IDAB)then
      do i = 1,nhelix
         do j = helix(i,1)+1,helix(i,2)
            istruct(ico(j-1),ico(j)) = .true.
         enddo
         do j = helix(i,1),helix(i,2)-1
            istruct(in(j),in(j+1)) = .true.
            do k = j+1,helix(i,2)
               istruct(ica(j),ica(k)) = .true.
            enddo
         enddo
      enddo
      do i = 1,nbeta
         do j = beta(i,1)+1,beta(i,2)
            istruct(ico(j-1),ico(j)) = .true.
         enddo
         do j = beta(i,1),beta(i,2)-1
            istruct(in(j),in(j+1)) = .true.
            do k = j+1,beta(i,2)
               if (distat2(j,k).lt.RCUTGO**2) then
                  istruct(ica(j),ica(k)) = .true.
               endif
            enddo
         enddo
      enddo
   endif


!! TARGET
!   if(IDAB)then
!      do i = 1,nhelix2
!         do j = helix2(i,1)+1,helix2(i,2)
!            istructtarg(ico2(j-1),ico2(j)) = .true.
!         enddo
!         do j = helix2(i,1),helix2(i,2)-1
!            istructtarg(in2(j),in2(j+1)) = .true.
!            do k = j+1,helix2(i,2)
!               istructtarg(ica2(j),ica2(k)) = .true.
!            enddo
!         enddo
!      enddo
!      do i = 1,nbeta2
!         do j = beta2(i,1)+1,beta2(i,2)
!            istructtarg(ico2(j-1),ico2(j)) = .true.
!         enddo
!         do j = beta2(i,1),beta2(i,2)-1
!            istructtarg(in2(j),in2(j+1)) = .true.
!            do k = j+1,beta2(i,2)
!               if (distat2targ(j,k).lt.RCUTGO**2) then
!                  istructtarg(ica2(j),ica2(k)) = .true.
!               endif
!            enddo
!         enddo
!      enddo
!   endif

!
!!
stepPts%step(1)=stepPot(0._SGL, 0._SGL)
stepPts%step(2)=stepPot(0._SGL, 0._SGL)
!
   do j = 2, natom
   do i = 1, j-1
      if (istruct(i,j))then ! .and.istructtarg(i,j)) then
  !    write(unit_o,*)res(i)," ", rnum(i),"   ", res(j)," ",rnum(j)
      stepPts(i,j) = getStepSSec(SIGMAGO, EGO, sqrt(distat2(i,j)), .true.)
      END IF
   enddo
   enddo
!
   do j = 2,natom
   do i = 1, j-1
      if (.not.istruct(i,j).and.rb(i,j).lt.1.e-20)then
         esolv = -sqrt(evdw(i) * evdw(j)) * 0.5 * FVDW
         rvdwij = rvdw(i) + rvdw(j)
         potlk = gfree(i) * vol(j) / xlamb(i) + gfree(j) * vol(i) / xlamb(j)
         esolv = esolv - potlk / rvdwij**2 * 0.09 * 0.5 * FSOLV
         ecoul = qq(i)*qq(j) * EPS ! Si es constant dielectrica millor qq2/EPS!!
         if(abs(qq(i)*qq(j)).gt.1d-10.and.EPS.gt.1.d-10)then
            stepPts(i,j) = getStepCoul (rvdwij, DPSINT, DPSEXT, ecoul, HPS, .false.)
         else
            if(esolv.lt.1.d-10)then
               if(rnum(j).gt.(rnum(i)+1)) stepPts(i,j) = getStepHFob(rvdwij, DHF, esolv, .false.)
            else
               stepPts(i,j) = getStepHFil(rvdwij, DHF, esolv, .false.)
            endif
         endif
      endif
!
      xsum(i,j) = 0.5 * (1./xm(i) + 1./xm(j))
   enddo
   enddo
! precalcul parells
   bpairs = count(rb.gt.1.e-20)
   allocate(blist(bpairs,2), stat = ioerr)
   call errorAllocMem(ioerr, 'Bond pairs')
   bpairs=0
   do j = 2,natom
   do i = 1, j-1
      if (rb(i,j).gt.1e-20) then
  !       if (abs(rnum(i)-rnum(j))>3) write(unit_o,*)"BOND ",res(i)," ", rnum(i),"   ", res(j)," ",rnum(j)
         bpairs = bpairs + 1
         blist(bpairs,1) = i
         blist(bpairs,2) = j
         stepPts(i,j)%tipInt = BOND
         stepPts(i,j)%active = .false.
      endif
   enddo
   enddo
!
   allocate(nblist(natom), stat=ioerr)
   call errorAllocMem(ioerr, 'Non Bonded List 1')
   do i = 1,natom
        nblist(i) = allocateintPList(natom, ioerr)
        call errorAllocMem(ioerr, 'Non Bonded List')
   enddo
   call activateStepPot(stepPts, r, rcutcoul2, rcutsolv2, natom,nblist,xsum)
   write (unit_o, '(" Initial Pairs list")')
   write (unit_o, '(" ==================")')
   write (unit_o, '(" Total:              ",i9)') natom * (natom - 1) / 2
   write (unit_o, '(" Bonded:             ",i9)') bpairs
   write (unit_o, '(" Non Bonded:         ",i9)') count(stepPts%active)
   write (unit_o, '("   Secondary Struc.  ",i9)') count(stepPts%tipInt.eq.SS)
   write (unit_o, '("   Electrostatic     ",i9)') count(stepPts%tipInt.eq.COUL.and.stepPts%active)
   write (unit_o, '("   Hydrophilic       ",i9)') count(stepPts%tipInt.eq.HFil.and.stepPts%active)
   write (unit_o, '("   Hydrophobic       ",i9)') count(stepPts%tipInt.eq.HFob.and.stepPts%active)
!
   write (unit_o,'(/," System setup completed",/)')
!
! Deallocate, mantenim ICA!!!
   deallocate (in, ico, atp,qq, gfree, cov, hbs, helix, beta, distat2)
 !  DEALLOCATE(helix2,beta2,istructtarg,in2,ico2,ica2)
   deallocate (istruct)
!
   call cpu_time(tsetup)

! suma l'energia potencial de la conformacio inicial
!
   call calcEpot (natom, r, stepPts, EGO, epotgo, epotfis)
!

   unit_rst =  openFn(files, '-rst')
   unit_ener = openFn(files, '-ener')
!
!   call writeSnapshot (unit_traj, 0, r, atom, res, rnum, natom)
!
   iev=0
   call thermalize(seed, iev, natom, temp, xmassa, v, xm, ekin)
   write (unit_o, '(" Initial energy evaluation completed")')
   write (unit_o, '(" Epot (kcal/mol) = ",f10.3," Ekin (kcal/mol) = ",f10.3)') epotfis, ekin / FACTE

   if (iwr) write (unit_ener, '("# Initial energy ",4f10.3)') epotfis, epotgo, ekin/FACTE, &
   epotfis+epotgo+ekin/FACTE
   if (iwr) write (unit_ener, '("# time(ps), epotFis, epotGo, ekin0, ekin1, Etot, Temp, Rejected, Total")')

   write (unit_o,*)
   allocate (tpart(natom), ipart(natom), stat=ioerr)
   call errorAllocMem(ioerr, ' Colision times')
 !
temps = 0.
tactot = 0.
tpart = 1.
ipart = 0
accepted=0
nonaccepted=0
tpedro=0._DBL
tpedroNM=0._DBL
tpedroNM_prevNM=tpedroNM
wiblocs=0
r_prevNM=r
time_prevNM=0._DBL
score_prevNM=scoreprev
xbeta_prevNM=xbeta
tpedro_prevNM=tpedro
nm_rejects=0
firstNM=.TRUE.
solap=sqrt(acc_coef)
NMcontrol=(solap > NMCutoff)
NMtargOnce=.TRUE.
! WRITE (unit_o, '("ERROR VECTOR SIZE",1x,I3)') error_vector_size
stopping=.FALSE.
DRMSD=1E-4
! WRITE (unit_o,'("TOLERNACE IN RMSD ",f8.4)' ) Drmsd
RMSDlim=1.2*LOG10(REAL(nres)) ! Esto es provisional, tiene que ser dato de entrada para un
! usuario cualquiera para mi es mas comodo si es automatico asi no tengo que mirar cada estructura
! WRITE (unit_o,'("RMSD minimum to stop ",f8.4)' ) RMSDlim
RMSD_CUTOFF=1.5*RMSDlim
IF ( RMSD_STOP < 1.E10 ) RMSDlim=RMSD_STOP


! DIMS IS ACTIVATED INFO
      write (unit_o,'(" DIMS: Maxwell-Demon dynamics activated")')
      write (unit_o,'("     Re-scoring time (ps)     :", f8.2)') trect*1e12
      write (unit_o,'("     RMSD min to stop         :", f8.2)') RMSDlim
      write (unit_o,'("     Initial RMSd CA          :", f8.2)') error0
      write (unit_o,'("     Acceptance               :", f8.2)') acceptance
      write (unit_o,'("     Maximum Time (ps)        :", f9.1)') MAX_TIME
      write (unit_o,'("     Number of eigeinvectors  :",4x ,I2)') nevecs
      write (unit_o,'("     Normal Mode Control      :", 4x,L)') NMcontrol
      Write (unit_o, '("     Transition Vector Overlap:" ,4x,f7.5 )')acc_coef
      Write (unit_o, '("     Superposition            :" ,4x,L  )')Superpose
      write (unit_o,'("     NM time (ps)             :", f8.2)') timeMC_NM
      write (unit_o,'("     RMSD_CUTOFF              :", f8.2)') RMSD_CUTOFF
      write (unit_o,'("     Time Equlibrium          :", f8.2),/') time_in_eq

!----------------------------------------------------------------------------

     DO WHILE (.not.stopping)
      tacum = 0.
!----------------------------------------------------------------------------
      do while(tacum.lt.TSNAP)
         tacrect=0.
!----------------------------------------------------------------------------
         do while(tacrect.lt.TRECT)
            taccorr = 0.
            iev = 0
            ierr = 0
!----------------------------------------------------------------------------
            do while(taccorr.lt.TCORR)
               tacact = 0.
               call colisioBond(bpairs, blist, r, v, rb, drb, xm, xsum, natom)
               call activateStepPot(stepPts, r, rcutcoul2, rcutsolv2, natom, nblist, xsum)
               call colisioNonBond(stepPts, temps, nblist, r, v, rhc, rnum, atom, xm, ierr, &
                                   TMIN, natom, isolv)
               call inici (nblist, tpart, ipart, natom) !tpart(i) ipart(i) contenen la primera col per atom
! evolucio temporal
!------------------------------------------------------------------------------
               do while (tacact.lt.TACT)
! busca quina es la peopera colisio
                  tevent = 1.
                  call nextCol(mem1, mem2, npair1, npair2, tevent, ipart, tpart, natom, nblist)
                  tevent1 = tevent - temps

! translacio i variacio dels temps
                  do i=1,natom ! mantenim la versio inline degut a la barreja de tipus de real !!
                     r(i)%x = r(i)%x + v(i)%x * tevent1
                     r(i)%y = r(i)%y + v(i)%y * tevent1
                     r(i)%z = r(i)%z + v(i)%z * tevent1
                  enddo
!
                  tacact = tacact + tevent1
                  temps = tevent
!
                  call updateV(r, v, nblist(mem1)%iData(npair1)%deltak, xm, &
                     nblist(mem1)%iData(npair1)%xsum, mem1, mem2, natom)

                  iev = iev + 1
! ara calcula els temps de colisio per a les dues particules que han xocat
                  call updateXocPart(mem1, nblist, temps, r, v, TMIN, natom)
                  nblist(mem1)%iData(npair1)%timp = 1.
                  call getMinTCol (mem1, nblist, tpart(mem1), ipart(mem1), natom)
                  call updateXocPart(mem2, nblist, temps, r, v, TMIN, natom)
                  nblist(mem2)%iData(npair2)%timp = 1.
                  call getMinTCol (mem2, nblist, tpart(mem2), ipart(mem2), natom)
               enddo
! end do while (tacact.lt.tact)------------------------------------------------
               taccorr = taccorr + tacact
            enddo
! end do while(taccorr.lt.tcorr)-----------------------------------------------
            ekin0 = calcEkin(v,xm,natom)
            if (TCALC.eq.1) call thermalize (seed, iev, natom, temp, xmassa, v, xm, ekin)
            tacrect = tacrect + taccorr
         enddo
! end do while(tacrect.lt.trect)-----------------------------------------------
            idec=idec+1
            nonaccepted=nonaccepted+1
         if ( ibloc > 4) then ! Restraints do not apply up to 5 (arbitrary) time step to allow relaxation in case
         ! of a terrible structure
            score = distanceCA(r,rtarg,w,ica,natom,nres)
            auxiliar=seed+iev
            if (MCCheck(error,auxiliar, xbeta, scoreprev, score)) then
               r=rprev
            else
               rprev=r
               iacc=iacc+1
               scoreprev=score
               tacum=tacum+tacrect
               tactot=tactot+tacrect
               accepted=accepted+1
            endif
         else
               iacc=iacc+1
               tacum=tacum+tacrect
               tactot=tactot+tacrect
         endif
      enddo
! end DO while(tacum.lt.tsnap) -------------------------------------------------

! TIME IS OUT
IF (tactot*1.e12 > MAX_TIME ) stopping=.TRUE.

      ibloc=NINT(tactot/tsnap)
      call calcEpot (natom, r, stepPts, ego, epotgo, epotfis)

     ! Actualizacion de las coordenadas sobre las que se calcula el RMSd

     if (iwr) write (buffer_ener(ibloc-wiblocs),&
      '(1x,7f10.2,f5.2)') tactot*1.e12, epotfis, epotgo, ekin0/FACTE, ekin/FACTE, &
   epotfis+epotgo+ekin0/FACTE,ekin0/natom, iacc*1./idec

call rmsd(nres, r(ica(:)),rtarg(ica(:)), 1, U,  center2, center1, error,.FALSE.)
IF(Superpose) call pdbsuperpos(natom, U , center1, center2,r)

! SELF_STOP
errcont=errcont+1
IF (errcont > error_vector_size ) errcont=1
error_vector(errcont,2)=error
error_vector(errcont,1)=tactot*1e12
! END Bloque SELF-STOP
!
 write (buffer_frames(ibloc-wiblocs,1), '("MODEL",8X,I6)') ibloc
  do i = 1,natom
     write (buffer_frames(ibloc-wiblocs,i+1), '("ATOM",2X,I5,2X,A3,1X,A3,3X,I3,4X,3F8.3)') i, atom(i), res(i), rnum(i), r(i)
  enddo
  write (buffer_frames(ibloc-wiblocs,natom+2), '("ENDMDL")')

   aux_acc= accepted*1./nonaccepted

  IF((tactot-tpedro) > xbeta_update) THEN
      tpedro=tactot
      ! este if es un desatre, deberia ser una funcion continua con forma de escalon
      call get_xbeta(xbeta,aux_acc,error,rmsd_cutoff,lower_acceptance,acceptance,acc_tolerance)
      accepted=0
      NonAccepted=0
  ENDIF

 write ( buffer_out(ibloc-wiblocs),&
'("Time(ps):",f8.2," Ev: ", i4," Temp (K):" f8.2," Acc = ", f4.2," B= ",f4.2 ," RMSd =",f6.2," L= "f4.2," NM ="L1,1x,I1)')&
       tactot*1.e12, iev, ekin0/natom,iacc*1./idec,xbeta, error,solap,NMcontrol,nm_rejects


! DECIDE TO STOP
slope=lreg(error_vector_size,error_vector)
!write(unit_o, '("SLOPE ", f10.5)') slope
IF (abs(slope) < Drmsd.and..not.firstNM .AND. error < RMSDlim) stopping=.TRUE.
!
! Linea agregada para contrlolar mejor los calculos del servidor
IF (error < 0.75*RMSDlim ) NMControl=.FALSE.
!
!
! RMSD CUT_OFF
!! NM Monte Carlo
entraNM: IF( (tactot -tpedroNM)*1e12 > timeMC_NM .and..not.stopping) THEN

IF (abs(slope) < Drmsd.and..not.firstNM ) THEN
!     WRITE (unit_o, '("RECOMPUTING WEIGHTS")')
     call MCWeigth(nres,natom,r, rtarg, ica,w)
     scoreprev=1E10
     xbeta=0.9*xbeta
END IF

   solap= overlap(natom, nres,ica,r,rtarg,trans_vec)
!   write(unit_o,*) "Solap ", solap
!
    entraReCalVT:  IF (solap < NMCutoff .and.NMControl ) THEN
      Rmsdlejos: IF (error > rmsd_cutoff ) THEN
         CALL transition_NM(natom,nres,ica,unit_o,r,rtarg,acc_coef,trans_vec)
         IF(acc_coef < 0.3) NMControl=.FALSE. !Cambiar 0.3 por el cuadrado de NMCUT
      ELSE Rmsdlejos
         IF (NMtargOnce) CALL transition_NM_targ(natom,nres,ica,unit_o,r,rtarg,acc_coef,trans_vec)
         NMtargOnce=.FALSE.
         IF(acc_coef < 0.3) NMControl=.FALSE.
      END IF Rmsdlejos
   END IF entraReCalVT
!
   DO i=1,nres
     diff_frames(i)= r(ica(i)) -r_prevNM(ica(i))
   end do
   call makeunit_3N(nres,diff_frames)
      aux_dotBIG=0
   do i=1, nres
     aux_dotBIG= aux_dotBIG+dotDP(diff_frames(i), trans_vec(i))
   end do
!   write(unit_o,'(" *****  NEW MC  ",f8.4)')aux_dotBIG
entraCheckNM: IF ( MCCheck_NM(auxiliar,NMCutoff,aux_dotBIG) .or.firstNM .or..not.NMcontrol) THEN
!  Write(unit_o,'("ACEPTO: ",L1,1X,I1,1X,"PROJ ",f7.5)') NMcontrol,nm_rejects,aux_dotBIG
    wiblocs=ibloc
    do j=1,buffer_size
     do i=1,natom+2
       if(len_trim(buffer_frames(j,i)) < 3 ) CYCLE
       write(unit_traj,'(A54)') buffer_frames(j,i)
     enddo
    enddo
  !
  do i=1,buffer_size
    if(len_trim(buffer_out(i)) < 5 ) CYCLE
    write(unit_o, '(A120)')buffer_out(i)
   write(unit_ener, '(A120)')buffer_ener(i)
  end do
  buffer_out=""
  buffer_frames=""
  buffer_ener=""
  r_prevNM=r ! keep coordinates
  tpedro=tactot
  time_prevNM=tactot
  score_prevNM=scoreprev
  xbeta_prevNM=xbeta
  tpedro_prevNM=tpedro
  tpedroNM=tactot
  tpedroNM_prevNM=tpedroNM
!  nm_rejects=0
  firstNM=.FALSE.
 ELSE entraCheckNM
! Write(unit_o,'("NO ACEPTO: ",L1,1X,I1,1X,"PROJ ",f7.5)') NMcontrol,nm_rejects,aux_dotBIG
tactot=time_prevNM
r=r_prevNM
rprev=r_prevNM
scoreprev=score_prevNM
xbeta=xbeta_prevNM*0.9
buffer_out=""
buffer_frames=""
buffer_ener=""
tpedro=tpedro_prevNM
tpedroNM=tpedroNM_prevNM
xbeta_prevNM=xbeta_prevNM*0.9
nm_rejects=nm_rejects+1
errcont=0
END IF entraCheckNM

IF (nm_rejects >= NM_rejections) NMcontrol=.FALSE.

END IF entraNM

entraCheckStop: IF(stopping) THEN
  do j=1,buffer_size
    do i=1,natom+2
       if(len_trim(buffer_frames(j,i)) < 3 ) CYCLE
       write(unit_traj,'(A54)') buffer_frames(j,i)
    enddo
    if(len_trim(buffer_out(j)) < 5 ) CYCLE
       write(unit_o, '(A120)')buffer_out(j)
       write(unit_ener, '(A120)')buffer_ener(j)
  enddo
  if (tactot*1.e12 < MAX_TIME) THEN
  Write(unit_o,'(" **** RMSD converged *****","Current",f6.2,1x,"Fixed",1x,f6.2,3x,"SLOPE = ",f8.3,/)')error,rmsdlim,slope
  END IF
  EXIT
END IF entraCheckStop
!
enddo
! end DO iblock-----------------------------------------------------------------
   close(5)

!
   do i = 1,natom
      rsp(i) = DPtoSP (r(i))
   enddo
!
   write (unit_rst) natom
   write (unit_rst) (rsp(i), i = 1,natom)
   write (unit_o, '(/," Restart coordinates written",/)')
   call cpu_time(tfin)
   write (unit_o, '(" T I M I N G S ")')
   write (unit_o, '(" ============= ")')
   write (unit_o, '(" Setup: " f10.2," s")') tsetup-tinit
   write (unit_o, '(" Traj:  ",f10.2," s (",f10.2," ns/h)"),/') tfin-tsetup, temps*1.e15 / (tfin-tsetup) * 3.6 / 1000.
   end program discrete
!=============================================================================
