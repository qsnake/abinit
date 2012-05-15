!{\src2tex{textfont=tt}}
!!****f* ABINIT/readeig
!! NAME
!! readeig
!!
!! FUNCTION
!! Reading the - q point information: eigenvalues + eigenvectors
!! The information is in a specially formatted file written by anaddb
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ineig = file pointer
!!  natom = number of atoms in the unit cell
!!  nqpt = number of q points
!!
!! OUTPUT
!!  acell = lattice parameters
!!  rprim = orientation of the lattice parameters
!!  typat = atom type
!!  xred = atom fractional coordinates
!!  eigval = phonon eigenvalues
!!  eigvect = phonon eigenvectors
!!  qpoint = coordinates of the q points
!!
!! NOTES
!!  To be merged with invars7w in a later ABINIT version
!!
!! PARENTS
!!      lwf
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine readeig(acell,eigval,eigvect,ineig,natom,nqpt,qpoint,rprim,typat,xred)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'readeig'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ineig,natom,nqpt
!arrays
 integer,intent(out) :: typat(natom)
 real(dp),intent(out) :: acell(3),eigval(nqpt,3*natom)
 real(dp),intent(out) :: eigvect(nqpt,3*natom,natom,3,2),qpoint(nqpt,3)
 real(dp),intent(out) :: rprim(3,3),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,ii,iqpt,jatom,jj,matom,mqpt,nsym,ntypat
 character(len=150) :: tmpstr
 character(len=500) :: message
!arrays
 integer,allocatable :: symrel(:,:,:)

!******************************************************************
!BEGIN EXECUTABLE SECTION

!initialization of the matrices
!DEBUG
!write(std_out,*) 'readeig: enter'
!write(std_out,*) '  nqpt=',nqpt
!write(std_out,*) '  natom=',natom
!ENDDEBUG

 read(ineig,'(4i5)') matom,mqpt,nsym,ntypat
!write(std_out,*) matom,mqpt,nsym,ntypat

 if (nqpt .ne. mqpt) then
   write(message, '(a,a,a,a,i7,a,a,i7,a,a,a,a)' ) ch10,&
&   ' chkilwf : ERROR -',ch10,&
&   '  The number of the q points in the input file',nqpt,ch10,&
&   '  does not agree with the number of q points in the ifc-output file',mqpt,ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify the energy windows in the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if


 ABI_ALLOCATE(symrel,(3,3,nsym))

 read(ineig,'(a)') tmpstr
!write(std_out,*) tmpstr
 read(ineig,'(3f12.8)') acell(:)
!write(std_out,*) 'readeig: acell',acell(:)

 read(ineig,'(a)') tmpstr
!write(std_out,*) tmpstr
 read(ineig,'(3f12.8)') rprim(:,1)
!write(std_out,*) 'readeig: rprim',rprim(:,1)
 read(ineig,'(3f12.8)') rprim(:,2)
!write(std_out,*) 'readeig: rprim',rprim(:,2)
 read(ineig,'(3f12.8)') rprim(:,3)
!write(std_out,*) 'readeig: rprim',rprim(:,3)

 read(ineig,'(a)') tmpstr
!write(std_out,*) tmpstr
 read(ineig,'(i5)') matom
!write(std_out,*) '  matom=',matom
 read(ineig,'(a)') tmpstr
!write(std_out,*) tmpstr
 read(ineig,'(i5)') ntypat
!write(std_out,*) '  ntypat=',ntypat
 read(ineig,'(a)') tmpstr
!write(std_out,*) tmpstr
 do ii=1,natom
   read(ineig,'(i5,3f12.8)') typat(ii),(xred(jj,ii),jj=1,3)
!  write(std_out,*) 'readeig: xred',typat(ii),xred(:,ii)
 end do
 read(ineig,'(a)') tmpstr
!write(std_out,*) tmpstr
!write(std_out,*) 'special extra line, nsym=',nsym
 do ii=1,nsym
   read(ineig,'(9i4)') symrel(:,:,ii)
!  write(std_out,'(9i4)') symrel(:,:,ii)
 end do
 read(ineig,'(a)') tmpstr
!write(std_out,*) tmpstr
 read(ineig,'(i5)') mqpt
!write(std_out,*) 'readeig: mqpt',mqpt


 do iqpt=1,nqpt
!  write(std_out,*) 'readeig: qpoint no',iqpt
   read(ineig, '(a)') tmpstr
!  write(std_out,*) 'readeig: tmpstr',tmpstr
   read(ineig,'(3f9.5)' ) (qpoint(iqpt,ii),ii=1,3)
!  write(std_out,*) 'readeig: read qpoint',(qpoint(iqpt,ii),ii=1,3)
   read(ineig,'(a)') tmpstr
!  write(std_out,*) 'readeig: tmpstr',tmpstr

   do jj=1,3*natom,5
     if (3*natom-jj<5) then
       read(ineig, '(5d15.9)') (eigval(iqpt,ii),ii=jj,3*natom)
!      write(std_out,'(5d15.9)') (eigval(iqpt,ii),ii=jj,3*natom)
     else
       read(ineig, '(5d15.9)') (eigval(iqpt,ii),ii=jj,jj+4)
!      write(std_out,'(5d15.9)') (eigval(iqpt,ii),ii=jj,jj+4)
     end if
   end do

   do iatom=1,3*natom
     read(ineig,*)
     do jatom=1,natom
       read(ineig,'(i5,3f18.12)') jj,eigvect(iqpt,iatom,jatom,:,1)
!      write(std_out,'(a,i6,a,3f20.16)') ' atom ',jj,' real ',eigvect(iqpt,iatom,jatom,:,1)
       read(ineig,'(i5,3f18.12)') jj,eigvect(iqpt,iatom,jatom,:,2)
!      write(std_out,'(a,i6,a,3f20.16)') ' atom ',jj,' imag ',eigvect(iqpt,iatom,jatom,:,2)
     end do
   end do

 end do

!DEBUG
!write(std_out,*) 'readeig:exit'
!ENDDEBUG

end subroutine readeig
!!***
