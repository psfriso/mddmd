!===============================================================================   
 subroutine readCommLine (files, unit_i, unit_o, NFILES)
 USE commLine
  integer, intent(IN) :: NFILES
  integer, intent(OUT) :: unit_i, unit_o
  type(commLineOption), intent(INOUT) :: files(NFILES) 
!
  call inputArgs(files)
  unit_i = openFn(files, '-i')
  if (fileName(files,'-o').ne.'log') then
        unit_o = openFn(files, '-o')
  else
        unit_o = 6
  endif
 end subroutine readCommLine
!===============================================================================   
 subroutine programHeader(unit_o)
   integer unit_o
   write (unit_o, *) "================================================="
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "=               DISCRETE  (0.2.3)               ="
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "=     A. Emperador, J. L. Gelpi, M.Orozco       =" 
   write (unit_o, *) "=                                               ="
   write (unit_o, *) "=                  (c) 2011                     =" 
   write (unit_o, *) "================================================="
   write (unit_o, *)
 end subroutine programHeader
!===============================================================================   
 subroutine errorAllocmem (ioerr, text)
  integer ioerr
  character(len=*) text
  if (ioerr.ne.0) then
    write (0, '("Error allocating memory",a30)') text
    stop 1
  endif
 end subroutine errorAllocmem
!===============================================================================   
 subroutine writeSnapshot (unit_traj, ibloc, r, atom, res, rnum, natom)
 use geometryDP
  integer, intent(IN) :: natom, unit_traj, ibloc
  type(pointDP), intent(IN) :: r(natom)
  character(len=4), intent(IN) :: atom(natom), res(natom)
  integer, intent(IN) :: rnum(natom)
  integer i
  write (unit_traj, '("MODEL",8X,I4)') ibloc
  do i = 1,natom
     write (unit_traj, '("ATOM",2X,I5,2X,A3,1X,A3,3X,I3,4X,3F8.3)') i, atom(i), res(i), rnum(i), r(i)
  enddo
  write (unit_traj, '("ENDMDL")')
 end subroutine writeSnapshot
!===============================================================================   

 subroutine writeSnapshotCA (unit_traj, ibloc, r, atom, res, rnum, natom,nres,ica)
 use geometryDP
  integer, intent(IN) :: natom, unit_traj, ibloc
  integer, intent(in), dimension(nres) :: ica
  integer, intent(in) :: nres
  type(pointDP), intent(IN) :: r(natom)
  character(len=4), intent(IN) :: atom(natom), res(natom)
  integer, intent(IN) :: rnum(natom)
  integer i
  write (unit_traj, '("MODEL",8X,I4)') ibloc
  do i = 1,nres
     write (unit_traj, '("ATOM",2X,I5,2X,A3,1X,A3,3X,I3,4X,3F8.3)') i, atom(ica(i)), &
             res(ica(i)), rnum(ica(i)), r(i)
  enddo
  write (unit_traj, '("ENDMDL")')
 end subroutine writeSnapshotCA
