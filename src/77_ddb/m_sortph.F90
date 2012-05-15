
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_sortph

 use m_profiling

use defs_basis

 implicit none

 private

 complex(dpc),allocatable :: eigvecLast(:,:)

 public :: end_sortph
 public :: sortph

contains

!{\src2tex{textfont=tt}}
!!****f* ABINIT/end_sortph
!! NAME
!! end_sortph
!!
!! FUNCTION
!! Deallocate array for sortph
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  Only deallocation
!!
!! NOTES
!!
!! PARENTS
!!      mkphbs,thm9
!!
!! CHILDREN
!!
!! SOURCE
subroutine end_sortph()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'end_sortph'
!End of the abilint section

  if (allocated(eigvecLast))  then
    ABI_DEALLOCATE(eigvecLast)
  end if
end subroutine end_sortph
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/sortph
!! NAME
!! sortph
!!
!! FUNCTION
!! Sort the energies in order to have fine phonon
!! dispersion curves
!! It is best not to include the gamma point in the list
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (FDortu,MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! MODIFIED
!! Takeshi Nishimatsu
!!
!! INPUTS
!!  eigvec(2*3*natom*3*natom)= contain
!!  the eigenvectors of the dynamical matrix.
!!  displ(2*3*natom*3*natom)= contain
!!   the displacements of atoms in cartesian coordinates.
!!   The first index means either the real or the imaginary part,
!!   The second index runs on the direction and the atoms displaced
!!   The third index runs on the modes.
!!  filnam=name of output files
!!   hacmm1,hartev,harthz,xkb= different conversion factors
!!  natom= number of atom
!!  phfrq(3*natom)= phonon frequencies in Hartree
!!  udispl=unit number for output of phonon eigendisplacements
!!  ufreq=unit number for output of phonon frequencies
!!
!! OUTPUT
!!  (only writing ?)
!!
!! NOTES
!! Called by one processor only
!!  FIXME: no one deallocates eigvecLast!
!!
!! PARENTS
!!      mkphbs,thm9
!!
!! CHILDREN
!!
!! SOURCE
subroutine sortph(eigvec,displ,filnam, natom,phfrq,udispl,ufreq)

use defs_basis
use defs_datatypes
use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sortph'
!End of the abilint section

implicit none

!Arguments -----------------------------------
!scalars
integer,intent(in) :: natom,udispl,ufreq
character(len=fnlen),intent(in) :: filnam
!arrays
real(dp),intent(in) :: eigvec(2,3,natom,3,natom)
real(dp),intent(in) :: displ(2*3*natom*3*natom)
real(dp),intent(in)    :: phfrq(3*natom)

!Local variables-------------------------------
!scalars
integer :: iatom,imode,j,idir1,idir2,ipert1,ipert2,i1,i2
character(len=fnlen) :: file_displ,file_freq
!arrays
integer :: i(1)
logical     ::               mask(3*natom)
real(dp)    ::           phfrqNew(3*natom)
complex(dpc) ::           displIn(3*natom,3*natom)
complex(dpc) ::           displNew(3*natom,3*natom)
complex(dpc) ::           eigvecIn(3*natom,3*natom)
complex(dpc) ::          eigvecNew(3*natom,3*natom)
complex(dpc) ::   transpose_eigvec(3*natom,3*natom)
real(dp)    ::     abs_similarity(3*natom,3*natom)  !|<displNew|displLast>|
character(len=20) :: fmt_phfrq

! *********************************************************************

!DEBUG
!write(std_out,*)' sortph : enter '
!call flush(6)
!ENDDEBUG

 do ipert2=1,natom
   do idir2=1,3
     i2=idir2+(ipert2-1)*3
     do ipert1=1,natom
       do idir1=1,3
         i1=idir1+(ipert1-1)*3
         eigvecIn(i1,i2)=cmplx(eigvec(1,idir1,ipert1,idir2,ipert2),eigvec(2,idir1,ipert1,idir2,ipert2))
         displIn(i1,i2)=cmplx(displ(1+2*(i1-1)+2*3*natom*(i2-1)),displ(2+2*(i1-1)+2*3*natom*(i2-1)))
       end do
     end do
   end do
 end do

 if(.not.allocated(eigvecLast)) then
   file_freq  = trim(filnam)//".freq" !---------------------------------------------------
   write(std_out,'(a,a)' )' sortph : opening file ',trim(file_freq)
   open(ufreq,FILE=trim(file_freq),STATUS='replace',ACCESS='sequential',ACTION='write')
   write(std_out,'(a,a,a)' )' sortph : file ',trim(file_freq),' opened '
   file_displ = trim(filnam)//".displ" !--------------------------------------------------
   write(std_out,'(a,a)' )' sortph : opening file ',trim(file_displ)
   open(udispl,FILE=trim(file_displ),STATUS='replace',ACCESS='sequential',ACTION='write')
   write(std_out,'(a,a,a)' )' sortph : file ',trim(file_displ),' opened '
   ABI_ALLOCATE(eigvecLast,(3*natom,3*natom))
   phfrqNew(:)   =  phfrq(:)
   displNew(:,:) =  displIn(:,:)
   eigvecNew(:,:) = eigvecIn(:,:)
 else
!  Avoid gfortran 4.2.1 bug, with which you CANNOT conjg(transpose(displ))
   transpose_eigvec = transpose(eigvecIn)
   abs_similarity = abs(matmul(conjg(transpose_eigvec),eigvecLast))
   mask(:) = .true.
   do j = 1, 3*natom
     i(:) = maxloc( abs_similarity(:,j), mask(:) )
     mask(i(1)) = .false.
     phfrqNew(j)   =    phfrq(i(1))
     displNew(:,j) =  displIn(:,i(1))
     eigvecNew(:,j) = eigvecIn(:,i(1))
   end do
 end if


!Write frequencies in a file
 write(fmt_phfrq,'(a,i3,a)') '(', 3*natom, 'e14.6)'
 write(ufreq,fmt_phfrq) (phfrqNew(j),j=1,3*natom)

!write displacements in a file
 do imode=1,3*natom
   do iatom=1,natom
     write(udispl,'(e14.6)') &
     sqrt(       displNew(3*(iatom-1)+1,imode)*   &
&     conjg(displNew(3*(iatom-1)+1,imode)) + &
&     displNew(3*(iatom-1)+2,imode)*   &
&     conjg(displNew(3*(iatom-1)+2,imode)) + &
&     displNew(3*(iatom-1)+3,imode)*   &
&     conjg(displNew(3*(iatom-1)+3,imode)) )
   end do
 end do

 eigvecLast(:,:) = eigvecNew(:,:)

!DEBUG
!write(std_out,*)' sortph : exit '
!call flush(6)
!ENDDEBUG

end subroutine sortph
!!***

end module m_sortph
!!***
