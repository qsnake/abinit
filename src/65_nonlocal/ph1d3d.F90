!{\src2tex{textfont=tt}}
!!****f* ABINIT/ph1d3d
!! NAME
!! ph1d3d
!!
!! FUNCTION
!! Compute the three-dimensional phase factor $e^{i 2 \pi (k+G) cdot xred}$
!! from the three one-dimensional factors, the k point coordinates,
!! and the atom positions, for all planewaves which fit in the fft box.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (XG,DR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iatom, jatom= bounds of atom indices in ph1d for which
!!                                        ph3d has to be computed
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  ------- Taken away in beautification because unused MS -------
!!  kpt(3)=k point in terms of recip. translations
!!  --------------------------------------------------------------
!!  matblk= dimension of ph3d
!!  natom= dimension of ph1d
!!  npw=number of plane waves
!!  n1,n2,n3=dimensions of fft box (ngfft(3)).
!!  phkxred(2,natom)=phase factors exp(2 pi k.xred)
!!  ph1d(2,(2*n1+1)*natom+(2*n2+1)*natom+(2*n3+1)*natom)=exp(2Pi i G xred) for
!!   vectors (Gx,0,0), (0,Gy,0) and (0,0,Gz)
!!   with components ranging from -nj <= Gj <= nj
!!
!! OUTPUT
!!  ph3d(2,npw_k,matblk)=$e^{2 i \pi (k+G) cdot xred}$ for vectors (Gx,Gy,Gz),
!!   and for atoms in the range iatom to jatom with respect to ph1d
!!
!! PARENTS
!!      ctocprj,dens_in_sph,dyfnl3,eltfrnl3,energy,forstrnps,getcprj,ks_ddiago
!!      ladielmt,lavnl,m_cprj_bspline,m_shirley,m_wfs,nonlop_pl,nonlop_ylm
!!      nstpaw3,nstwf3,nstwf4,partial_dos_fractions,prctfvw1,prctfvw2,rhofermi3
!!      suscep_stat,vso_realspace_nonlop,vtorho,vtorho3,wfconv,wffile
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ph1d3d(iatom,jatom,kg_k,matblk,natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph1d3d'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatom,jatom,matblk,n1,n2,n3,natom,npw_k
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom) !,kpt(3)
 real(dp),intent(in) :: phkxred(2,natom)
 real(dp),intent(out) :: ph3d(2,npw_k,matblk)

!Local variables-------------------------------
!scalars
 integer :: i1,ia,iatblk,ig,shift1,shift2,shift3
 real(dp) :: ph12i,ph12r,ph1i,ph1r,ph2i,ph2r,ph3i,ph3r,phkxi,phkxr
 character(len=500) :: message
!arrays
 real(dp),allocatable :: ph1kxred(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' ph1d3d : enter '
!ENDDEBUG

 if(matblk-1 < jatom-iatom)then
   write(message, '(a,a,a,a,a,a,a,a,i4,a,a,i4,a,i4,a)' ) ch10,&
&   ' ph1d3d : BUG -',ch10,&
&   '  Input natom-1 must be larger or equal to jatom-iatom,',ch10,&
&   '  while their value is : ',ch10,&
&   '  natom-1 = ',natom-1,ch10,&
&   '  jatom=',jatom,', iatom=',iatom,'.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

 ABI_ALLOCATE(ph1kxred,(2,-n1:n1))

!ia runs from iatom to jatom
 do ia=iatom,jatom

!  iatblk runs from 1 to matblk
   iatblk=ia-iatom+1
   shift1=1+n1+(ia-1)*(2*n1+1)
   shift2=1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1)
   shift3=1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
!  Compute product of phkxred by phase for the first component of G vector
   phkxr=phkxred(1,ia)
   phkxi=phkxred(2,ia)
!  DEBUG (needed to compare with version prior to 2.0)
!  phkxr=1.0d0
!  phkxi=0.0d0
!  ENDDEBUG
   do i1=-n1,n1
     ph1kxred(1,i1)=ph1d(1,i1+shift1)*phkxr-ph1d(2,i1+shift1)*phkxi
     ph1kxred(2,i1)=ph1d(2,i1+shift1)*phkxr+ph1d(1,i1+shift1)*phkxi
   end do

!  Compute tri-dimensional phase factor
!  $OMP PARALLEL DO PRIVATE(ig,ph1i,ph1r,ph2i,ph2r,ph3i,ph3r,ph12i,ph12r)&
!  $OMP&SHARED(iatblk,kg_k,npw_k,ph1d,ph1kxred,ph3d,shift2,shift3)
   do ig=1,npw_k
     ph1r=ph1kxred(1,kg_k(1,ig))
     ph1i=ph1kxred(2,kg_k(1,ig))
     ph2r=ph1d(1,kg_k(2,ig)+shift2)
     ph2i=ph1d(2,kg_k(2,ig)+shift2)
     ph3r=ph1d(1,kg_k(3,ig)+shift3)
     ph3i=ph1d(2,kg_k(3,ig)+shift3)
     ph12r=ph1r*ph2r-ph1i*ph2i
     ph12i=ph1r*ph2i+ph1i*ph2r
     ph3d(1,ig,iatblk)=ph12r*ph3r-ph12i*ph3i
     ph3d(2,ig,iatblk)=ph12r*ph3i+ph12i*ph3r
   end do
!  $OMP END PARALLEL DO
 end do

 ABI_DEALLOCATE(ph1kxred)

!DEBUG
!write(std_out,*)' ph1d3d : exit '
!ENDDEBUG

end subroutine ph1d3d
!!***
