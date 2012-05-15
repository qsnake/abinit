!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_fsurf_1band
!! NAME
!! get_fsurf_1band
!!
!! FUNCTION
!! calculate Fermi surface in tetrahedra for 1 band 1 sppol
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtset=dataset structure for symmetry information
!! eigen_in(nkpt)=eigenenergies for each k point
!! fermie=Fermi energy
!! klatt=lattice vectors for full kpoint grid
!! kpt_fullbz=full BZ kpoints
!! mtetra= maximum number of tetrahedra
!! nfiner= sublattice of kpoints is nfiner times finer
!! nkpt=number of irreducible kpoints
!! nkpt_fullbz=number of kpoints in full brillouin zone
!! (THIS VARIABLE WAS UNUSED AND REMOVED DURING PEAUTIFICATION. PMA. ) ntetra=number of tetrahedra
!! tetra_full(4,2,mtetra)=for each irred tetrahedron, the list of k point vertices
!!     1 = irred and 1 = fullkpt
!! (THIS VARIABLE WAS UNUSED AND REMOVED DURING PEAUTIFICATION. PMA. )tetra_mult(mtetra)=
!!    for each irred tetrahedron, its multiplicity
!! tetra_wrap(3,4,mtetra) = integer to see if a given tetrahedron is wrapped
!!   around the edge of the Brillouin zone at +-0.5 (i.e. one apex is on the
!!   other side of the BZ), or even further into neighboring cells
!!   (tetra_wrap > 1). Real kpoint = kpoint_IBZ + tetra_wrap(:,isummit,itetra)
!! tolfermi=energy tolerance for surface, wrt fermie
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      matr3inv,sort_dp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine get_fsurf_1band(dtset,eigen_in,fermie,klatt,kpt_fullbz,&
&        mtetra,nfiner,nkpt_fullbz,tetra_full,tetra_wrap,tolfermi)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_fsurf_1band'
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mtetra,nfiner,nkpt_fullbz
 real(dp),intent(in) :: fermie,tolfermi
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: tetra_full(4,2,mtetra)
 integer,intent(in) :: tetra_wrap(3,4,mtetra)
 real(dp),intent(in) :: eigen_in(dtset%nkpt),klatt(3,3)
 real(dp),intent(in) :: kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
!scalars
 integer :: itetra,ii
!integer :: if1,if2,if3,ii,isummit,jsummit
 !real(dp) :: det1,det2,det3,det4,distka
 real(dp) :: fermidiff_coeff
!arrays
 integer :: ind_dum(4),ind_int(4),ind_wrp(3,4)
! real(dp) :: base_kpt(3)
 real(dp) :: eigen_1tetra(4),finer_klatt(3,3)
 real(dp) :: inv_klatt(3,3),k1(3),k2(3),k3(3),k4(3),normal_vect(3)
 real(dp) :: tmp1(3),tmpa(3),tmpb(3),tmpc(3)
 real(dp) :: tmpd(3)
! real(dp) :: tmp2(3),tmp3(3),tmp4(3),tmp_kpt(3),tmpka(3)
! *********************************************************************

!DEBUG
!write(std_out,*) ' get_fsurf_1band : enter'
!write(std_out,*) '   eigen_in ', eigen_in
!write(std_out,*) '   klatt', klatt(:,1)
!write(std_out,*) '        ', klatt(:,2)
!write(std_out,*) '        ', klatt(:,3)
!write(std_out,*) 'nfiner = ', nfiner
!ENDDEBUG

 finer_klatt(:,1) = klatt(:,1)/nfiner
 finer_klatt(:,2) = klatt(:,2)/nfiner
 finer_klatt(:,3) = klatt(:,3)/nfiner

 call matr3inv (klatt,inv_klatt)

!
!for each tetrahedron
!
 do itetra=1,mtetra
!  DEBUG
!  write(std_out,*) ' get_fsurf_1band : itetra ', itetra, ' / ', mtetra
!  ENDDEBUG
!  Here we need the original ordering to reference the correct irred kpoints
   eigen_1tetra(1) = eigen_in(tetra_full(1,1,itetra))
   eigen_1tetra(2) = eigen_in(tetra_full(2,1,itetra))
   eigen_1tetra(3) = eigen_in(tetra_full(3,1,itetra))
   eigen_1tetra(4) = eigen_in(tetra_full(4,1,itetra))
   ind_int(1) = 1
   ind_int(2) = 2
   ind_int(3) = 3
   ind_int(4) = 4
   call sort_dp(4,eigen_1tetra,ind_int,tol14)
   ind_dum(1) = tetra_full(ind_int(1),2,itetra)
   ind_dum(2) = tetra_full(ind_int(2),2,itetra)
   ind_dum(3) = tetra_full(ind_int(3),2,itetra)
   ind_dum(4) = tetra_full(ind_int(4),2,itetra)

   ind_wrp(:,1) = tetra_wrap(:,ind_int(1),itetra)
   ind_wrp(:,2) = tetra_wrap(:,ind_int(2),itetra)
   ind_wrp(:,3) = tetra_wrap(:,ind_int(3),itetra)
   ind_wrp(:,4) = tetra_wrap(:,ind_int(4),itetra)

!  DEBUG
!  for this tetra, draw edges
!  do isummit=1,4
!  do jsummit=isummit+1,4
!  write(std_out,*) ' TETRAS O ',  10.0*kpt_fullbz(:,ind_dum(isummit))
!  do if1=1,nfiner-1
!  write(std_out,*) ' TETRAS C ',  10.0*kpt_fullbz(:,ind_dum(isummit)) + &
!  &     10.0*if1/nfiner * (kpt_fullbz(:,ind_dum(jsummit)) - &
!  &                         kpt_fullbz(:,ind_dum(isummit)))
!  end do
!  end do
!  end do
!  ENDDEBUG

!  DEBUG
!  write(std_out,*) '   fermie, eigen_1tetra ', fermie, ' : ', eigen_1tetra(:)
!  ENDDEBUG

!  if tetrahedron is out of range, cycle
   if ( fermie+tolfermi < eigen_1tetra(1) .or. &
&   fermie-tolfermi > eigen_1tetra(4)) then
     cycle
   end if

!  
!  WARNING ! If e4-e1 < tolfermi then the whole tetrahedron is in the surface.
!  If the eigenvalues are really degenerate, there will be a FP overflow
!  


!  DEBUG
!  write(std_out,*) ' get_fsurf_1band : eigen_1tetra(1:4)',eigen_1tetra(1:4)
!  write(std_out,*) '   ind_dum = ', ind_dum
!  ENDDEBUG

!  
!  wrap coordinates if necessary to get connected tetrahedron
!  may be partially outside the first BZ !
!  
   k1(:) = kpt_fullbz(:,ind_dum(1))
   k2(:) = kpt_fullbz(:,ind_dum(2))
   k3(:) = kpt_fullbz(:,ind_dum(3))
   k4(:) = kpt_fullbz(:,ind_dum(4))
   do ii=1,3
     if (ind_wrp(ii,1) .ne. 0) then
       k1(ii) = k1(ii)+float(ind_wrp(ii,1))
     end if
     if (ind_wrp(ii,2) .ne. 0) then
       k2(ii) = k2(ii)+float(ind_wrp(ii,2))
     end if
     if (ind_wrp(ii,3) .ne. 0) then
       k3(ii) = k3(ii)+float(ind_wrp(ii,3))
     end if
     if (ind_wrp(ii,4) .ne. 0) then
       k4(ii) = k4(ii)+float(ind_wrp(ii,4))
     end if
   end do

!  DEBUG
!  write(std_out,*)  kpt_fullbz(:,ind_dum(1))
!  write(std_out,*)  kpt_fullbz(:,ind_dum(2))
!  write(std_out,*)  kpt_fullbz(:,ind_dum(3))
!  write(std_out,*)  kpt_fullbz(:,ind_dum(4))
!  write(std_out,*)  'N1 ', k1(:)
!  write(std_out,*)  'N2 ', k2(:)
!  write(std_out,*)  'N3 ', k3(:)
!  write(std_out,*)  'N4 ', k4(:)
!  ENDDEBUG

!  
!  find points of intersection between Fsurf and segments 1-2 1-3 1-4
!  NB: points may be outside tetrahedron, does matter !!!!
!  with these if else ifs should be ok
!  
!  fermie btw e1 and e234
   if (eigen_1tetra(2) - fermie > tol6 .and. &
&   fermie - eigen_1tetra(1) > tol6) then
     tmpa(:) = k1(:) + (k2(:)-k1(:)) &
&     * (fermie-eigen_1tetra(1)) / (eigen_1tetra(2)-eigen_1tetra(1))
     tmpb(:) = k1(:) + (k3(:)-k1(:)) &
&     * (fermie-eigen_1tetra(1)) / (eigen_1tetra(3)-eigen_1tetra(1))
     tmpc(:) = k1(:) + (k4(:)-k1(:)) &
&     * (fermie-eigen_1tetra(1)) / (eigen_1tetra(4)-eigen_1tetra(1))
!    fermie btw e12 and e34
   else if (eigen_1tetra(3) - fermie > tol6 .and. &
&     fermie - eigen_1tetra(2) > tol6) then
     tmpa(:) = k1(:) + (k3(:)-k1(:)) &
&     * (fermie-eigen_1tetra(1)) / (eigen_1tetra(3)-eigen_1tetra(1))
     tmpb(:) = k1(:) + (k4(:)-k1(:)) &
&     * (fermie-eigen_1tetra(1)) / (eigen_1tetra(4)-eigen_1tetra(1))
     tmpc(:) = k2(:) + (k3(:)-k2(:)) &
&     * (fermie-eigen_1tetra(2)) / (eigen_1tetra(3)-eigen_1tetra(2))
     tmpd(:) = k2(:) + (k2(:)-k4(:)) &
&     * (fermie-eigen_1tetra(4)) / (eigen_1tetra(4)-eigen_1tetra(2))
!    fermie btw e1=e2=e3 and e4
   else if (eigen_1tetra(4) - fermie > tol6 .and. &
&     fermie - eigen_1tetra(3) > tol6) then
     tmpa(:) = k4(:) + (k1(:)-k4(:)) &
&     * (fermie-eigen_1tetra(4)) / (eigen_1tetra(4)-eigen_1tetra(1))
     tmpb(:) = k4(:) + (k2(:)-k4(:)) &
&     * (fermie-eigen_1tetra(4)) / (eigen_1tetra(4)-eigen_1tetra(2))
     tmpc(:) = k4(:) + (k3(:)-k4(:)) &
&     * (fermie-eigen_1tetra(4)) / (eigen_1tetra(4)-eigen_1tetra(3))
!    all 4 degenerate eigenvalues: use 3 corners for plane
   else
     tmpa(:) = k2(:)
     tmpb(:) = k3(:)
     tmpc(:) = k4(:)
   end if

!  DEBUG
!  write(std_out,*) ' got 3 points of intersection btw aretes and Fsurf: '
!  write(std_out,*) tmpa(:)
!  write(std_out,*) tmpb(:)
!  write(std_out,*) tmpc(:)
!  write(std_out,*) tmpd(:)
!  ENDDEBUG

!  normal_vect is normal to fermisurface plane
   normal_vect(1) = (tmpb(2)-tmpa(2))*(tmpc(3)-tmpa(3)) &
&   - (tmpb(3)-tmpa(3))*(tmpc(2)-tmpa(2))
   normal_vect(2) = (tmpb(3)-tmpa(3))*(tmpc(1)-tmpa(1)) &
&   - (tmpb(1)-tmpa(1))*(tmpc(3)-tmpa(3))
   normal_vect(3) = (tmpb(1)-tmpa(1))*(tmpc(2)-tmpa(2)) &
&   - (tmpb(2)-tmpa(2))*(tmpc(1)-tmpa(1))
!  DEBUG
!  write(std_out,*) ' normal vector to plane ', normal_vect(:)
!  ENDDEBUG

   tmp1(:) = k1(:) - tmpa(:)
   fermidiff_coeff = (fermie-eigen_1tetra(1))&
&   / sqrt(tmp1(1)**2 + tmp1(2)**2 + tmp1(3)**2)

!  
!  print out homogeneous 2D grid on polyhedron
!  
!  write(std_out,*) 'fermi surface: homogeneous 2D grid '
!  do if2=0,nfiner
!  tmp1(:) = tmpb(:) + (nfiner-if2)*(tmpc(:)-tmpb(:))/nfiner
!  do if1=0,nfiner
!  write(std_out,*) 'FSURF ', tmpa(:) + if1*tmp1(:)/nfiner
!  end do
!  end do
!  write(std_out,*) 'FSURF ', tmpa(:) + 0.1*(tmpb(:)-tmpa(:)) + 0.1*(tmpc(:)-tmpa(:))
!  write(std_out,*) 'FSURF ', tmpa(:) + 0.1*(tmpb(:)-tmpa(:)) + 0.9*(tmpc(:)-tmpb(:))
!  write(std_out,*) 'FSURF ', tmpa(:) + 0.1*(tmpc(:)-tmpa(:)) + 0.9*(tmpb(:)-tmpc(:))
   write(std_out,*) 'FSURF ', 0.3*tmpa(:) + 0.7*0.5*(tmpb(:)+tmpc(:))

!  
!  find points on finer 3D grid which are closer than tolfermi to the surface
!  

!  do if3=0,nfiner-1
!  do if2=0,nfiner-1
!  do if1=0,nfiner-1
!  !DEBUG
!  write(std_out,*) 'if1,if2.if3 ',if1,if2,if3
!  !ENDDEBUG
!  
!  tmp_kpt(:) = base_kpt(:) + finer_klatt(:,1)*if1 &
!  &                             + finer_klatt(:,2)*if2 &
!  &                             + finer_klatt(:,3)*if3
!  !
!  !  find out if tmp_kpt is in tetrahedron
!  !
!  tmp1(:) = k1(:) - tmp_kpt(:)
!  tmp2(:) = k2(:) - tmp_kpt(:)
!  tmp3(:) = k3(:) - tmp_kpt(:)
!  tmp4(:) = k4(:) - tmp_kpt(:)
!  !DEBUG
!  write(std_out,*) 'vectors from test kpt to summits : '
!  write(std_out,*) tmp1(:)
!  write(std_out,*) tmp2(:)
!  write(std_out,*) tmp3(:)
!  write(std_out,*) tmp4(:)
!  !ENDDEBUG
!  det1 = tmp1(1)*(tmp2(2)*tmp3(3) - tmp2(3)*tmp3(2)) &
!  &         - tmp1(2)*(tmp2(1)*tmp3(3) - tmp2(3)*tmp3(1)) &
!  &         + tmp1(3)*(tmp2(1)*tmp3(2) - tmp2(2)*tmp3(1))
!  det2 = tmp1(1)*(tmp4(2)*tmp2(3) - tmp4(3)*tmp2(2)) &
!  &         - tmp1(2)*(tmp4(1)*tmp2(3) - tmp4(3)*tmp2(1)) &
!  &         + tmp1(3)*(tmp4(1)*tmp2(2) - tmp4(2)*tmp2(1))
!  det3 = tmp1(1)*(tmp3(2)*tmp4(3) - tmp3(3)*tmp4(2)) &
!  &         - tmp1(2)*(tmp3(1)*tmp4(3) - tmp3(3)*tmp4(1)) &
!  &         + tmp1(3)*(tmp3(1)*tmp4(2) - tmp3(2)*tmp4(1))
!  det4 = tmp2(1)*(tmp4(2)*tmp3(3) - tmp4(3)*tmp3(2)) &
!  &         - tmp2(2)*(tmp4(1)*tmp3(3) - tmp4(3)*tmp3(1)) &
!  &         + tmp2(3)*(tmp4(1)*tmp3(2) - tmp4(2)*tmp3(1))
!  !DEBUG
!  write(std_out,*) ' determinants of 4 tetrahedra : ', det1,det2,det3,det4
!  !ENDDEBUG
!  
!  !
!  !  calculate distance to Fermi plane in tetra
!  !
!  tmpka(:) = tmp_kpt(:) - tmpa(:)
!  distka = fermidiff_coeff*sqrt(tmpka(1)**2 + tmpka(2)**2 + tmpka(3)**2)
!  write(std_out,*) 'tmpka, distka, tolfermi ', tmpka(:), distka, tolfermi
!  
!  !
!  !  If point is inside tetrahedron
!  !
!  if ( (det1<tol6 .and. det2<tol6 .and. det3<tol6 .and. det4<tol6) &
!  &     .or.(det1>-tol6 .and. det2>-tol6 .and. det3>-tol6 .and. det4>-tol6) ) then
!  
!  write(std_out,*) 'found point inside tetrahedron : ', tmp_kpt(:)
!  !  if distance is smaller than tolfermi, then include point in fermisurface
!  if (distka < tolfermi) then
!  write(std_out,*) 'found point on fermi surface: ', distka
!  end if
!  end if
!  
!  end do
!  end do
!  end do
!  ! end do if1,if2,if3


 end do
!end do itetra

!DEBUG
!write(std_out,*)' get_fsurf_1band : exit '
!ENDDEBUG

end subroutine get_fsurf_1band
!!***
