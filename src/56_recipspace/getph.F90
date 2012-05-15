!{\src2tex{textfont=tt}}
!!****f* ABINIT/getph
!!
!! NAME
!! getph
!!
!! FUNCTION
!! Compute three factors of one-dimensional structure factor phase
!! for input atomic coordinates, for all planewaves which fit in fft box.
!! The storage of these atomic factors is made according to the
!! values provided by the index table atindx. This will save time in nonlop.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  natom=number of atoms in cell.
!!  n1,n2,n3=dimensions of fft box (ngfft(3)).
!!  xred(3,natom)=reduced atomic coordinates.
!!
!! OUTPUT
!!  ph1d(2,(2*n1+1)*natom+(2*n2+1)*natom+(2*n3+1)*natom)=exp(2Pi i G.xred) for
!!   integer vector G with components ranging from -nj <= G <= nj.
!!   Real and imag given in usual Fortran convention.
!!
!! PARENTS
!!      afterscfloop,bethe_salpeter,extrapwf,gstate,loop3dte,loper3
!!      m_cprj_bspline,m_hamiltonian,m_wfs,partial_dos_fractions,prcref
!!      prcref_PMA,respfn,scfcv,screening,sigma,wfconv,wffile
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine getph(atindx,natom,n1,n2,n3,ph1d,xred)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getph'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,natom
!arrays
 integer,intent(in) :: atindx(natom)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ii
 real(dp) :: arg

! *************************************************************************

 do ia=1,natom

!  Store the phase factor of atom number ia in place atindx(ia)
   i1=(atindx(ia)-1)*(2*n1+1)
   i2=(atindx(ia)-1)*(2*n2+1)+natom*(2*n1+1)
   i3=(atindx(ia)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)

   do ii=1,2*n1+1
     arg=two_pi*dble(ii-1-n1)*xred(1,ia)
     ph1d(re,ii+i1)=cos(arg)
     ph1d(im,ii+i1)=sin(arg)
   end do

   do ii=1,2*n2+1
     arg=two_pi*dble(ii-1-n2)*xred(2,ia)
     ph1d(re,ii+i2)=cos(arg)
     ph1d(im,ii+i2)=sin(arg)
   end do

   do ii=1,2*n3+1
     arg=two_pi*dble(ii-1-n3)*xred(3,ia)
     ph1d(re,ii+i3)=cos(arg)
     ph1d(im,ii+i3)=sin(arg)
   end do

 end do

end subroutine getph
!!***
