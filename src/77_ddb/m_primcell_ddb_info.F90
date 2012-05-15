!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_primcell_ddb_info
!!
!! NAME
!! m_primcell_ddb_info
!!
!! FUNCTION
!! Module for container object passing cell information to SC phonon calculation in abinit 
!! Container type is defined, and destruction
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!! 
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_primcell_ddb_info

 use m_profiling

 use defs_basis

 implicit none

type primcell_ddb_info
! scalars
  integer :: brav,mpert,msym,natom,nrpt,nsym,ntypat
  integer :: dipdip
  real(dp) :: ucvol

! arrays
  integer , pointer :: indsym(:,:,:) ! indsym(4,nsym,natom)=label given by subroutine symatm
  integer , pointer :: symrec(:,:,:) ! (3,3,nsym) recip space symops
  integer , pointer :: symrel(:,:,:) ! (3,3,nsym) real  space symops
  integer , pointer :: typat(:)      ! typat(natom)=integer label of each type of atom (1,2,...)

  real(dp), pointer :: acell(:)      ! acell(3)
  real(dp), pointer :: amu(:)        ! amu(ntypat)=mass of the atoms (atomic mass unit)
  real(dp), pointer :: dielt(:,:)    ! dielt(3,3)=dielectric tensor
  real(dp), pointer :: dyewq0(:,:,:) ! dyewq0(3,3,natom)=Ewald part of the dynamical matrix, at q=0
  real(dp), pointer :: gmet(:,:)     ! gmet(3,3)
  real(dp), pointer :: gprim(:,:)    ! gprim(3,3)
  real(dp), pointer :: rcan(:,:)     ! rcan(3,natom)=atomic position in canonical coordinates
  real(dp), pointer :: rmet(:,:)     ! rmet(3,3)
  real(dp), pointer :: rprim(:,:)    ! rprim(3,3)=dimensionless primitive translations in real space
  real(dp), pointer :: rpt(:,:)      ! rpt(3,nrpt)=canonical coordinates of the R points in the unit cell
  real(dp), pointer :: trans(:,:)    ! trans(3,natom)=atomic translations : xred = rcan + trans
  real(dp), pointer :: wghatm(:,:,:) ! wghatm(natom,natom,nrpt)=weights assoc to a pair of atoms + R vector
  real(dp), pointer :: xred(:,:)     ! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
  real(dp), pointer :: zeff(:,:,:)   ! zeff(3,3,natom)=effective charge on each atom

end type primcell_ddb_info

! Now the subroutines in the module
contains
!!***

!!****f* m_primcell_ddb_info/init_primcell_ddb_info
!!
!! NAME
!! init_primcell_ddb_info
!!
!! FUNCTION
!!  init and fill primcell_ddb_info
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  pcell= structure to allocate and fill
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE
subroutine init_primcell_ddb_info (pcell,brav,dipdip,mpert,msym,natom,nrpt,nsym,ntypat,ucvol,&
&    indsym,symrec,symrel,typat,&
&    acell,amu,dielt,dyewq0,gmet,gprim,rcan,rmet,rprim,rpt,trans,wghatm,xred,zeff)
 
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_primcell_ddb_info'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(primcell_ddb_info), intent(inout) :: pcell

 integer, intent(in) :: brav,mpert,msym,natom,nrpt,nsym,ntypat,dipdip
 real(dp), intent(in) :: ucvol

 integer, intent(in) :: indsym(4,nsym,natom)
 integer, intent(in) :: symrec(3,3,nsym)
 integer, intent(in) :: symrel(3,3,nsym)
 integer, intent(in) :: typat(natom)

 real(dp), intent(in) :: acell(3)
 real(dp), intent(in) :: amu(ntypat)
 real(dp), intent(in) :: dielt(3,3)
 real(dp), intent(in) :: dyewq0(3,3,natom)
 real(dp), intent(in) :: gmet(3,3)
 real(dp), intent(in) :: gprim(3,3)
 real(dp), intent(in) :: rcan(3,natom)
 real(dp), intent(in) :: rmet(3,3)
 real(dp), intent(in) :: rprim(3,3)
 real(dp), intent(in) :: rpt(3,nrpt)
 real(dp), intent(in) :: trans(3,natom)
 real(dp), intent(in) :: wghatm(natom,natom,nrpt)
 real(dp), intent(in) :: xred(3,natom)
 real(dp), intent(in) :: zeff(3,3,natom)

!Local variables-------------------------------

! *************************************************************************

! init dimensions
  pcell%brav = brav
  pcell%mpert = mpert
  pcell%msym = msym
  pcell%natom = natom
  pcell%nrpt = nrpt
  pcell%nsym = nsym
  pcell%ntypat = ntypat

! init scalar reals
  pcell%ucvol = ucvol

! allocate int
  ABI_ALLOCATE(pcell%indsym,(4,nsym,natom))
  ABI_ALLOCATE(pcell%symrec,(3,3,nsym))
  ABI_ALLOCATE(pcell%symrel,(3,3,nsym))
  ABI_ALLOCATE(pcell%typat,(natom))

! allocate real
  ABI_ALLOCATE(pcell%acell,(3))
  ABI_ALLOCATE(pcell%amu,(ntypat))
  ABI_ALLOCATE(pcell%dielt,(3,3))
  ABI_ALLOCATE(pcell%dyewq0,(3,3,natom))
  ABI_ALLOCATE(pcell%gmet,(3,3))
  ABI_ALLOCATE(pcell%gprim,(3,3))
  ABI_ALLOCATE(pcell%rcan,(3,natom))
  ABI_ALLOCATE(pcell%rmet,(3,3))
  ABI_ALLOCATE(pcell%rprim,(3,3))
  ABI_ALLOCATE(pcell%rpt,(3,nrpt))
  ABI_ALLOCATE(pcell%trans,(3,natom))
  ABI_ALLOCATE(pcell%wghatm,(natom,natom,nrpt))
  ABI_ALLOCATE(pcell%xred,(3,natom))
  ABI_ALLOCATE(pcell%zeff,(3,3,natom))

! init int
  pcell%indsym = indsym
  pcell%symrec = symrec
  pcell%symrel = symrel
  pcell%typat = typat
  pcell%dipdip = dipdip

! init real
  pcell%acell = acell
  pcell%amu = amu
  pcell%dielt = dielt
  pcell%dyewq0 = dyewq0
  pcell%gmet = gmet
  pcell%gprim = gprim
  pcell%rcan = rcan
  pcell%rmet = rmet
  pcell%rprim = rprim
  pcell%rpt = rpt
  pcell%trans = trans
  pcell%wghatm = wghatm
  pcell%xred = xred
  pcell%zeff = zeff
  
 end subroutine init_primcell_ddb_info
!!***

!!****f* m_primcell_ddb_info/read_primcell_ddb_info
!!
!! NAME
!! read_primcell_ddb_info
!!
!! FUNCTION
!!  read in and fill primcell_ddb_info from the file name given in input
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filename= name of file to read in
!!
!! OUTPUT
!!  t_primcell_ddb_info= structure to allocate and fill
!!
!! NOTES
!!
!! PARENTS
!!      scphon
!!
!! CHILDREN
!!
!! SOURCE
subroutine read_primcell_ddb_info (filename,pcell)
 
 use defs_basis
 use m_io_tools, only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_primcell_ddb_info'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=fnlen), intent(in) :: filename
 type(primcell_ddb_info), intent(inout) :: pcell

!Local variables-------------------------------

 integer :: unit,istat
 character(len=fnlen):: trimfilename
 character(len=13):: buffer
 
! *************************************************************************

  unit=get_unit()
  trimfilename = trim(filename)
  open (unit=unit, file=trimfilename)

! read in dimensions
  read(unit,'(a,I6)') buffer, pcell%brav
  read(unit,'(a,I6)') buffer, pcell%dipdip
  read(unit,'(a,I6)') buffer, pcell%mpert
  read(unit,'(a,I6)') buffer, pcell%msym
  read(unit,'(a,I6)') buffer, pcell%natom
  read(unit,'(a,I6)') buffer, pcell%nrpt
  read(unit,'(a,I6)') buffer, pcell%nsym
  read(unit,'(a,I6)') buffer, pcell%ntypat

! read in scalar reals
  read(unit,'(a,E20.10)') buffer, pcell%ucvol

! allocate int
  ABI_ALLOCATE(pcell%indsym,(4,pcell%nsym,pcell%natom))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%indsym'
  ABI_ALLOCATE(pcell%symrec,(3,3,pcell%nsym))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%symrec'
  ABI_ALLOCATE(pcell%symrel,(3,3,pcell%nsym))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%symrel'
  ABI_ALLOCATE(pcell%typat,(pcell%natom))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%typat'

! allocate real
  ABI_ALLOCATE(pcell%acell,(3))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%acell'
  ABI_ALLOCATE(pcell%amu,(pcell%ntypat))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%amu'
  ABI_ALLOCATE(pcell%dielt,(3,3))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%dielt'
  ABI_ALLOCATE(pcell%dyewq0,(3,3,pcell%natom))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%dyewq0'
  ABI_ALLOCATE(pcell%gmet,(3,3))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%gmet'
  ABI_ALLOCATE(pcell%gprim,(3,3))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%gprim'
  ABI_ALLOCATE(pcell%rcan,(3,pcell%natom))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%rcan'
  ABI_ALLOCATE(pcell%rmet,(3,3))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%rmet'
  ABI_ALLOCATE(pcell%rprim,(3,3))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%rprim'
  ABI_ALLOCATE(pcell%rpt,(3,pcell%nrpt))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%rpt'
  ABI_ALLOCATE(pcell%trans,(3,pcell%natom))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%trans'
  ABI_ALLOCATE(pcell%wghatm,(pcell%natom,pcell%natom,pcell%nrpt))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%wghatm'
  ABI_ALLOCATE(pcell%xred,(3,pcell%natom))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%xred'
  ABI_ALLOCATE(pcell%zeff,(3,3,pcell%natom))
  istat = ABI_ALLOC_STAT
  if (istat /= 0) write(std_out,*) 'read_primcell_ddb_info problem allocating pcell%zeff'

! read in int
  read(unit,'(a)')
  read(unit,'(8I6)') pcell%indsym
  read(unit,'(a)')
  read(unit,'(8I6)') pcell%symrec
  read(unit,'(a)')
  read(unit,'(8I6)') pcell%symrel
  read(unit,'(a)')
  read(unit,'(8I6)') pcell%typat

! read in real
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%acell
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%amu
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%dielt
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%dyewq0
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%gmet
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%gprim
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%rcan
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%rmet
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%rprim
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%rpt
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%trans
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%wghatm
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%xred
  read(unit,'(a)')
  read(unit,'(3E20.10)') pcell%zeff

  close(unit)
  
 end subroutine read_primcell_ddb_info
!!***


!!****f* m_primcell_ddb_info/write_primcell_ddb_info
!!
!! NAME
!! write_primcell_ddb_info
!!
!! FUNCTION
!!  write out primcell_ddb_info to the file name given in input
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filename= name of file to read in
!!  t_primcell_ddb_info= structure to allocate and fill
!!
!! OUTPUT
!!   writes to file
!!
!! NOTES
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE
subroutine write_primcell_ddb_info (filename,pcell)
 
 use defs_basis
 use m_io_tools, only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_primcell_ddb_info'
!End of the abilint section

 implicit none

!Arguments ------------------------------------

 character(len=fnlen), intent(in) :: filename
 type(primcell_ddb_info), intent(in) :: pcell

!Local variables-------------------------------

 integer :: unit
 character(len=fnlen):: trimfilename

! *************************************************************************
  
  unit=get_unit()
  trimfilename = trim(filename)
  open (unit=unit, file=trim(filename), form='formatted',status='unknown')

! read out dimensions
  write(unit,'(a,I6)') 'pcell%brav   ', pcell%brav
  write(unit,'(a,I6)') 'pcell%dipdip ', pcell%dipdip
  write(unit,'(a,I6)') 'pcell%mpert  ', pcell%mpert
  write(unit,'(a,I6)') 'pcell%msym   ', pcell%msym
  write(unit,'(a,I6)') 'pcell%natom  ', pcell%natom
  write(unit,'(a,I6)') 'pcell%nrpt   ', pcell%nrpt
  write(unit,'(a,I6)') 'pcell%nsym   ', pcell%nsym
  write(unit,'(a,I6)') 'pcell%ntypat ', pcell%ntypat

! write out scalar reals
  write(unit,'(a,E20.10)') 'pcell%ucvol  ', pcell%ucvol

! write out int
  write(unit,'(a)') 'pcell%indsym'
  write(unit,'(8I6)') pcell%indsym
  write(unit,'(a)') 'pcell%symrec'
  write(unit,'(8I6)') pcell%symrec
  write(unit,'(a)') 'pcell%symrel'
  write(unit,'(8I6)') pcell%symrel
  write(unit,'(a)') 'pcell%typat'
  write(unit,'(8I6)') pcell%typat

! write out real
  write(unit,'(a)') 'pcell%acell'
  write(unit,'(3E20.10)') pcell%acell
  write(unit,'(a)') 'pcell%amu'
  write(unit,'(3E20.10)') pcell%amu
  write(unit,'(a)') 'pcell%dielt'
  write(unit,'(3E20.10)') pcell%dielt
  write(unit,'(a)') 'pcell%dyewq0'
  write(unit,'(3E20.10)') pcell%dyewq0
  write(unit,'(a)') 'pcell%gmet'
  write(unit,'(3E20.10)') pcell%gmet
  write(unit,'(a)') 'pcell%gprim'
  write(unit,'(3E20.10)') pcell%gprim
  write(unit,'(a)') 'pcell%rcan'
  write(unit,'(3E20.10)') pcell%rcan
  write(unit,'(a)') 'pcell%rmet'
  write(unit,'(3E20.10)') pcell%rmet
  write(unit,'(a)') 'pcell%rprim'
  write(unit,'(3E20.10)') pcell%rprim
  write(unit,'(a)') 'pcell%rpt'
  write(unit,'(3E20.10)') pcell%rpt
  write(unit,'(a)') 'pcell%trans'
  write(unit,'(3E20.10)') pcell%trans
  write(unit,'(a)') 'pcell%wghatm'
  write(unit,'(3E20.10)') pcell%wghatm
  write(unit,'(a)') 'pcell%xred'
  write(unit,'(3E20.10)') pcell%xred
  write(unit,'(a)') 'pcell%zeff'
  write(unit,'(3E20.10)') pcell%zeff

  close(unit)
  
 end subroutine write_primcell_ddb_info
!!***

!!****f* m_primcell_ddb_info/destroy_primcell_ddb_info
!!
!! NAME
!! destroy_primcell_ddb_info
!!
!! FUNCTION
!!  deallocate stuoff in primcell_ddb_info
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      anaddb,scphon
!!
!! CHILDREN
!!
!! SOURCE
subroutine destroy_primcell_ddb_info (pcell)
 
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_primcell_ddb_info'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(primcell_ddb_info), intent(inout) :: pcell

!Local variables-------------------------------

! *************************************************************************
  if (associated(pcell%indsym))  then
    ABI_DEALLOCATE(pcell%indsym)
  end if
  if (associated(pcell%symrec))  then
    ABI_DEALLOCATE(pcell%symrec)
  end if
  if (associated(pcell%symrel))  then
    ABI_DEALLOCATE(pcell%symrel)
  end if
  if (associated(pcell%typat ))  then
    ABI_DEALLOCATE(pcell%typat)
  end if
  if (associated(pcell%acell ))  then
    ABI_DEALLOCATE(pcell%acell)
  end if
  if (associated(pcell%amu   ))  then
    ABI_DEALLOCATE(pcell%amu)
  end if
  if (associated(pcell%dielt ))  then
    ABI_DEALLOCATE(pcell%dielt)
  end if
  if (associated(pcell%dyewq0))  then
    ABI_DEALLOCATE(pcell%dyewq0)
  end if
  if (associated(pcell%gmet  ))  then
    ABI_DEALLOCATE(pcell%gmet)
  end if
  if (associated(pcell%gprim ))  then
    ABI_DEALLOCATE(pcell%gprim)
  end if
  if (associated(pcell%rcan  ))  then
    ABI_DEALLOCATE(pcell%rcan)
  end if
  if (associated(pcell%rmet  ))  then
    ABI_DEALLOCATE(pcell%rmet)
  end if
  if (associated(pcell%rprim ))  then
    ABI_DEALLOCATE(pcell%rprim)
  end if
  if (associated(pcell%rpt   ))  then
    ABI_DEALLOCATE(pcell%rpt)
  end if
  if (associated(pcell%trans ))  then
    ABI_DEALLOCATE(pcell%trans)
  end if
  if (associated(pcell%wghatm))  then
    ABI_DEALLOCATE(pcell%wghatm)
  end if
  if (associated(pcell%xred  ))  then
    ABI_DEALLOCATE(pcell%xred)
  end if
  if (associated(pcell%zeff  ))  then
    ABI_DEALLOCATE(pcell%zeff)
  end if

 end subroutine destroy_primcell_ddb_info

end module m_primcell_ddb_info
!!***
