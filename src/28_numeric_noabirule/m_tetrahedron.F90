!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_tetrahedron
!! NAME
!! m_tetrahedron
!!
!! FUNCTION
!!  module for tetrahedron interpolation of DOS and similar quantities
!!  depends on sort_dp and on m_kpt_rank
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2012 ABINIT group (MJV)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tetrahedron

 use m_profiling

implicit none

private

! tetrahedron geometry object

type, public :: t_tetrahedron
   integer :: ntetra 
   real(8) :: vv  ! volume of the tetrahedra
   integer,pointer :: tetra_full(:,:,:) !(4,2,ntetra)
   integer,pointer :: tetra_mult(:)     !(ntetra)
   integer,pointer :: tetra_wrap(:,:,:) !(3,4,ntetra)
end type t_tetrahedron

public :: init_tetra
public :: get_tetra_weight
public :: nullify_tetra
public :: destroy_tetra
!!***


contains

!{\src2tex{textfont=tt}}
!!****f* ABINIT/nullify_tetra
!! NAME
!! nullify_tetra
!!
!! FUNCTION
!! nullify tetrahedra pointers
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  tetra = nullified pointers
!!
!! PARENTS
!!      ep_fs_weights,m_bz_mesh,tetrahedron,thmeig
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

subroutine nullify_tetra (tetra)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_tetra'
!End of the abilint section

 implicit none
 type(t_tetrahedron), intent(out) :: tetra
 
 nullify(tetra%tetra_full)
 nullify(tetra%tetra_mult)
 nullify(tetra%tetra_wrap)

end subroutine nullify_tetra
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/destroy_tetra
!! NAME
!! destroy_tetra
!!
!! FUNCTION
!! deallocate tetrahedra pointers if needed
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
!!  tetra = deallocated pointers etc...
!!
!! PARENTS
!!      ep_fs_weights,m_bz_mesh,tetrahedron,thmeig
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

subroutine destroy_tetra (tetra)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_tetra'
!End of the abilint section

 implicit none
 type(t_tetrahedron), intent(inout) :: tetra

 if (associated(tetra%tetra_full))  then
   ABI_DEALLOCATE(tetra%tetra_full)
 end if
 if (associated(tetra%tetra_mult))  then
   ABI_DEALLOCATE(tetra%tetra_mult)
 end if
 if (associated(tetra%tetra_wrap))  then
   ABI_DEALLOCATE(tetra%tetra_wrap)
 end if

end subroutine destroy_tetra
!!***




!{\src2tex{textfont=tt}}
!!****f* ABINIT/init_tetra
!! NAME
!! init_tetra
!!
!! FUNCTION
!! get tetrahedra characterized by apexes
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  indkpt(nkpt_fullbz)=indexes of irred kpoints equivalent to kpt_fullbz
!!  gprimd(3,3) = reciprocal space vectors
!!  klatt(3,3)=reciprocal of lattice vectors for full kpoint grid
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!
!! OUTPUT
!!  tetrahedra%tetra_full(4,2,ntetra)=for each tetrahedron,
!!     the different instances of the tetrahedron (fullbz kpoints)
!!  tetrahedra%tetra_mult(ntetra) = store multiplicity of each irred tetrahedron
!!  tetrahedra%tetra_wrap(3,4,ntetra) = store flag to wrap tetrahedron summit into IBZ
!!  tetrahedra%ntetra = final number of irred tetrahedra (dimensions of tetra_* remain larger)
!!  tetrahedra%vv = tetrahedron volume divided by full BZ volume
!!
!! PARENTS
!!      ep_fs_weights,m_bz_mesh,m_phdos,tetrahedron,thmeig
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

subroutine init_tetra (indkpt,gprimd,klatt,kpt_fullbz,nkpt_fullbz,&
&                tetrahedra, ierr, errorstring)

 use m_kptrank

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_tetra'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_fullbz
!arrays
 integer,intent(in) :: indkpt(nkpt_fullbz)
 real(8),intent(in) :: gprimd(3,3),klatt(3,3),kpt_fullbz(3,nkpt_fullbz)

 integer, intent(out) :: ierr
 character(len=80), intent(out) :: errorstring
 type(t_tetrahedron),intent(out) :: tetrahedra

!Local variables-------------------------------
! 3 dimensions, 4 summits, and 6 tetrahedra / kpoint box
!scalars
 integer :: ialltetra,ikpt2,ikpt_full,isummit,itetra,jalltetra,jsummit
 integer :: symrankkpt
 integer :: mtetra
 real(8) :: itmp,shift1,shift2,shift3, rcvol
 type(kptrank_type) :: kptrank_t

!arrays
 integer :: tetra_shifts(3,4,6),tmptetra(4)
 real(8) :: k1(3),k2(3),k3(3)

 integer,allocatable :: tetra_full_(:,:,:) !(4,2,mtetra)
 integer,allocatable :: tetra_mult_(:)     !(ntetra)
 integer,allocatable :: tetra_wrap_(:,:,:) !(3,4,ntetra)

! *********************************************************************

 ierr = 0
 errorstring = ""

 mtetra = 6 * nkpt_fullbz
 ABI_ALLOCATE(tetra_full_,(4,2,mtetra))
 ABI_ALLOCATE(tetra_mult_,(mtetra))
 ABI_ALLOCATE(tetra_wrap_,(3,4,mtetra))

 tetra_mult_(:) = 1
 tetra_full_(:,:,:) = 0
 tetra_wrap_(:,:,:) = 0

 tetra_shifts(:,1,1) = (/0,0,0/)
 tetra_shifts(:,2,1) = (/0,1,0/)
 tetra_shifts(:,3,1) = (/0,1,1/)
 tetra_shifts(:,4,1) = (/1,1,0/)
 tetra_shifts(:,1,2) = (/0,0,0/)
 tetra_shifts(:,2,2) = (/0,1,1/)
 tetra_shifts(:,3,2) = (/1,1,0/)
 tetra_shifts(:,4,2) = (/1,1,1/)
 tetra_shifts(:,1,3) = (/0,0,0/)
 tetra_shifts(:,2,3) = (/1,0,0/)
 tetra_shifts(:,3,3) = (/1,1,0/)
 tetra_shifts(:,4,3) = (/1,1,1/)
 tetra_shifts(:,1,4) = (/0,0,0/)
 tetra_shifts(:,2,4) = (/0,0,1/)
 tetra_shifts(:,3,4) = (/1,0,0/)
 tetra_shifts(:,4,4) = (/1,1,1/)
 tetra_shifts(:,1,5) = (/0,0,1/)
 tetra_shifts(:,2,5) = (/1,0,0/)
 tetra_shifts(:,3,5) = (/1,0,1/)
 tetra_shifts(:,4,5) = (/1,1,1/)
 tetra_shifts(:,1,6) = (/0,0,0/)
 tetra_shifts(:,2,6) = (/0,0,1/)
 tetra_shifts(:,3,6) = (/0,1,1/)
 tetra_shifts(:,4,6) = (/1,1,1/)

!make full k-point rank arrays
 call mkkptrank (kpt_fullbz,nkpt_fullbz,kptrank_t)

 ialltetra = 1
 do ikpt_full=1,nkpt_fullbz
   do itetra=1,6
     do isummit=1,4

       k1(:) = kpt_fullbz(:,ikpt_full) &
&       + tetra_shifts(1,isummit,itetra)*klatt(:,1) &
&       + tetra_shifts(2,isummit,itetra)*klatt(:,2) &
&       + tetra_shifts(3,isummit,itetra)*klatt(:,3)

!      Find full kpoint which is summit isummit of tetrahedron itetra around full kpt ikpt_full !
       call get_rank_1kpt (k1,symrankkpt,kptrank_t)
       ikpt2 = kptrank_t%invrank(symrankkpt)
       if (ikpt2 < 1) stop 'Error in ranking k-points. Do you have a very dense k-grid?' 

!      Store irreducible kpoint equivalent to kpt_fullbz(:,ikpt2)
       tetra_full_(isummit,1,ialltetra) = indkpt(ikpt2)
       tetra_full_(isummit,2,ialltetra) = ikpt2
       if (shift1>0.5d0) then
         tetra_wrap_(1,isummit,ialltetra) = 1
       else if (shift1<-0.5d0) then
         tetra_wrap_(1,isummit,ialltetra) = -1
       end if
       if (shift2>0.5d0) then
         tetra_wrap_(2,isummit,ialltetra) = 1
       else if (shift2<-0.5d0) then
         tetra_wrap_(2,isummit,ialltetra) = -1
       end if
       if (shift3>0.5d0) then
         tetra_wrap_(3,isummit,ialltetra) = 1
       else if (shift3<-0.5d0) then
         tetra_wrap_(3,isummit,ialltetra) = -1
       end if

!      sort itetra summits
       do jsummit=isummit,2,-1
         if ( tetra_full_(jsummit,1,ialltetra) &
&         <  tetra_full_(jsummit-1,1,ialltetra) ) then
           itmp = tetra_full_(jsummit,1,ialltetra)
           tetra_full_(jsummit,1,ialltetra) = tetra_full_(jsummit-1,1,ialltetra)
           tetra_full_(jsummit-1,1,ialltetra) = itmp
           itmp = tetra_full_(jsummit,2,ialltetra)
           tetra_full_(jsummit,2,ialltetra) = tetra_full_(jsummit-1,2,ialltetra)
           tetra_full_(jsummit-1,2,ialltetra) = itmp
!          keep fullbz_kpt tetrahedra points in same order
           itmp = tetra_wrap_(1,jsummit,ialltetra)
           tetra_wrap_(1,jsummit,ialltetra) = tetra_wrap_(1,jsummit-1,ialltetra)
           tetra_wrap_(1,jsummit-1,ialltetra) = itmp
           itmp = tetra_wrap_(2,jsummit,ialltetra)
           tetra_wrap_(2,jsummit,ialltetra) = tetra_wrap_(2,jsummit-1,ialltetra)
           tetra_wrap_(2,jsummit-1,ialltetra) = itmp
           itmp = tetra_wrap_(1,jsummit,ialltetra)
           tetra_wrap_(3,jsummit,ialltetra) = tetra_wrap_(3,jsummit-1,ialltetra)
           tetra_wrap_(3,jsummit-1,ialltetra) = itmp
         end if
       end do !  loop jsummit

     end do !  loop isummit

     if (ialltetra > mtetra) then
       write (errorstring, '(3a,i6,a,i6)' ) &
&       ' init_tetra: BUG - ',&
&       '  ialltetra > mtetra ',&
&       '  ialltetra=  ',ialltetra,', mtetra=',mtetra
       ierr = 1
       return
     end if
     ialltetra = ialltetra+1
   end do !  loop itetra
 end do !  loop ikpt_full

 call destroy_kptrank (kptrank_t)

 rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
& -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
& +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

!Volume of all tetrahedra should be the same as that of tetra 1
!this is the volume of 1 tetrahedron, should be coherent with
!notation in Lehmann & Taut
 k1(:) = gprimd(:,1)*klatt(1,1) &
& +  gprimd(:,2)*klatt(2,1) &
& +  gprimd(:,3)*klatt(3,1)
 k2(:) = gprimd(:,1)*klatt(1,2) &
& +  gprimd(:,2)*klatt(2,2) &
& +  gprimd(:,3)*klatt(3,2)
 k3(:) = gprimd(:,1)*klatt(1,3) &
& +  gprimd(:,2)*klatt(2,3) &
& +  gprimd(:,3)*klatt(3,3)
 tetrahedra%vv  = abs (k1(1)*(k2(2)*k3(3)-k2(3)*k3(2)) &
& -k1(2)*(k2(1)*k3(3)-k2(3)*k3(1)) &
& +k1(3)*(k2(1)*k3(2)-k2(2)*k3(1))) / 6.d0 / rcvol

!
!eliminate equivalent tetrahedra by symmetry and account
!for them in multiplicity tetra_mult
!
 tetrahedra%ntetra = mtetra
!FIXME: could we replace this with a ranking algorithm to avoid the O(tetrahedra%ntetra^2) step? For example:
!get tetrahedron rank - problem too many combinations in principle = nkpt_irred^4 - only a few used in practice
!sort ranks and keep indices
!
!eliminate equal rank tetrahedra and accumulate multiplicity into first one
 do ialltetra=tetrahedra%ntetra,2,-1
   do jalltetra=1,ialltetra-1
!    check if tetra are equivalent
     if (tetra_full_(1,1,ialltetra) == tetra_full_(1,1,jalltetra) .and. &
&        tetra_full_(2,1,ialltetra) == tetra_full_(2,1,jalltetra) .and. &
&        tetra_full_(3,1,ialltetra) == tetra_full_(3,1,jalltetra) .and. &
&        tetra_full_(4,1,ialltetra) == tetra_full_(4,1,jalltetra) ) then
! if higher index tetrahedron already has multiplicity > 1 we had a problem before
       if ( tetra_mult_(ialltetra) > 1) then
         write(errorstring,'(a,I6,/,a,3I6)') 'found an equiv tetra with mult > 1 :', tetra_mult_(ialltetra), &
&         '  ialltetra,jalltetra,tetrahedra%ntetra = ',&
&         ialltetra,jalltetra,tetrahedra%ntetra
         ierr = 2
         return
       end if
!      accumulate multiplicity and positions into equiv tetrahedron jalltetra
       tetra_mult_(jalltetra) = tetra_mult_(jalltetra) + tetra_mult_(ialltetra)
       tetra_mult_(ialltetra) = 0
       tetrahedra%ntetra = tetrahedra%ntetra-1
       exit
     end if
   end do ! do jalltetra
 end do ! do ialltetra

!
!pack irred tetrahedra
!
 do ialltetra=mtetra,1,-1
   if (tetra_mult_(ialltetra) /= 0) then
!    look for an earlier place to put the irred tetrahedron
     do jalltetra=1,ialltetra-1
       if (tetra_mult_(jalltetra) == 0) then
!        
!        swap tetrahedrons
!        
         tmptetra(:) = tetra_full_(:,1,jalltetra)
         tetra_full_(:,1,jalltetra) = tetra_full_(:,1,ialltetra)
         tetra_full_(:,1,ialltetra) = tmptetra(:)
         tmptetra(:) = tetra_full_(:,2,jalltetra)
         tetra_full_(:,2,jalltetra) = tetra_full_(:,2,ialltetra)
         tetra_full_(:,2,ialltetra) = tmptetra(:)
!        
!        swap wrap flags
!        
         tmptetra(:) = tetra_wrap_(1,:,jalltetra)
         tetra_wrap_(1,:,jalltetra) = tetra_wrap_(1,:,ialltetra)
         tetra_wrap_(1,:,ialltetra) = tmptetra(:)
         tmptetra(:) = tetra_wrap_(2,:,jalltetra)
         tetra_wrap_(2,:,jalltetra) = tetra_wrap_(2,:,ialltetra)
         tetra_wrap_(2,:,ialltetra) = tmptetra(:)
         tmptetra(:) = tetra_wrap_(3,:,jalltetra)
         tetra_wrap_(3,:,jalltetra) = tetra_wrap_(3,:,ialltetra)
         tetra_wrap_(3,:,ialltetra) = tmptetra(:)
!        
!        swap multiplicities (tetra_mult_(jalltetra) was 0)
!        
         tetra_mult_(jalltetra) = tetra_mult_(ialltetra)
!  commented 10/2011 : where does this test come from?
!         if (tetra_mult_(ialltetra) > tetrahedra%ntetra) then
!           write(errorstring,'(a, 2I6)') 'problem : multiplicity > tetrahedra%ntetra ', &
!&           tetra_mult_(ialltetra),tetrahedra%ntetra
!           ierr = 3
!           return
!         end if
         tetra_mult_(ialltetra) = 0
         exit
       end if
     end do ! do jalltetra

   end if
 end do ! do ialltetra

 ! transfer to new arrays
 ABI_ALLOCATE(tetrahedra%tetra_full,(4,2,tetrahedra%ntetra))
 tetrahedra%tetra_full = tetra_full_(:,:,1:tetrahedra%ntetra)
 ABI_ALLOCATE(tetrahedra%tetra_mult,(tetrahedra%ntetra))
 tetrahedra%tetra_mult = tetra_mult_(1:tetrahedra%ntetra)
 ABI_ALLOCATE(tetrahedra%tetra_wrap,(3,4,tetrahedra%ntetra))
 tetrahedra%tetra_wrap = tetra_wrap_(:,:,1:tetrahedra%ntetra)

 ABI_DEALLOCATE(tetra_full_)
 ABI_DEALLOCATE(tetra_mult_)
 ABI_DEALLOCATE(tetra_wrap_)
end subroutine init_tetra
!!***

!!****f* ABINIT/get_tetra_weight
!! NAME
!! get_tetra_weight
!!
!! FUNCTION
!! calculate integration weights and their derivatives
!! from Blochl et al PRB 49 16223
!!
!! COPYRIGHT
!! Copyright (C) 2003-2012 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! eigen_in(nkpt)=eigenenergies for each k point
!! enemin=minimal energy for DOS
!! enemax=maximal energy for DOS
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! nkpt_fullbz=number of kpoints in full brillouin zone
!! max_occ=maximal occupation number (2 for nsppol=1, 1 for nsppol=2)
!!
!! tetrahedra%ntetra=number of tetrahedra
!! tetrahedra%tetra_full(4,2,ntetra)=for each irred tetrahedron, the list of k point vertices
!!   1 -> irred kpoint   2 -> fullkpt
!! tetrahedra%tetra_mult(ntetra)=for each irred tetrahedron, its multiplicity
!! tetrahedra%vv = ratio of volume of one tetrahedron in reciprocal space to full BZ volume
!!
!! OUTPUT
!!  tweight(nkpt,nene) = integration weights for each irred kpoint from all adjacent tetrahedra
!!  dtweightde(nkpt,nene) = derivative of tweight wrt energy
!!
!! PARENTS
!!      ep_fs_weights,m_ebands,m_phdos,tetrahedron,thmeig
!!
!! CHILDREN
!!      sort_dp
!!
!! SOURCE

subroutine get_tetra_weight(eigen_in,enemin,enemax,&
&             max_occ,nene,nkpt,tetrahedra,&
&            tweight,dtweightde)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_tetra_weight'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nene,nkpt
 type(t_tetrahedron), intent(in) :: tetrahedra
 real(8),intent(in) :: enemax,enemin,max_occ
!arrays
 real(8),intent(in) :: eigen_in(nkpt)

 real(8),intent(out) :: dtweightde(nkpt,nene),tweight(nkpt,nene)

!Local variables-------------------------------
!  needed for gaussian replacement of Dirac functions
!  the three coefficients of the DOS as quadratic form,
!    in the interval [eig(ikpt-1), eig(ikpt)]
!    for ikpt = 1 we add a point below eigen(1) which doesnt
!    contribute to the DOS in any tetrahedron
!scalars
 integer :: ieps,itetra,nn1,nn2,nn3,nn4
 real(8) :: cc,cc1,cc2,cc3,dcc1de,dcc2de,dcc3de,dccde,deltaene,eps
 real(8) :: epsilon21,epsilon31,epsilon32,epsilon41,epsilon42,epsilon43
 real(8) :: gau_prefactor,gau_width,gau_width2,inv_epsilon21,inv_epsilon31
 real(8) :: inv_epsilon32,inv_epsilon41,inv_epsilon42,inv_epsilon43
 real(8) :: deleps1, deleps2, deleps3, deleps4 
 real(8) :: invepsum, cc_pre, dccde_pre
 real(8) :: cc1_pre, cc2_pre, cc3_pre
 real(8) :: cc_tmp, dccde_tmp
 real(8) :: dcc1de_pre, dcc2de_pre, dcc3de_pre

 real(8) :: tmp,volconst,volconst_mult
 real(8) :: tol6 = 1.d-6, tol14 = 1.d-14, zero = 0.d0
 real(8) :: sqrtpi = 1.7724538509055159d0

!arrays
 real(8) :: tweight_tmp(nene,4)
 real(8) :: dtweightde_tmp(nene,4)

 integer :: ind_dum(4)
 real(8) :: eigen_1tetra(4), dtweightde_t(nene, nkpt), tweight_t(nene, nkpt)

! *********************************************************************

 tweight_t = zero
 dtweightde_t = zero

 volconst = tetrahedra%vv/4.d0
 if (nene <= 1) then
   stop 'get_tetra_weight : Error: nene must be at least 2'
 else
   deltaene = (enemax-enemin) / (nene-1)
 end if

!
!for each tetrahedron
!
 do itetra=1,tetrahedra%ntetra
   tweight_tmp = zero
   dtweightde_tmp = zero

   volconst_mult = max_occ*volconst*float(tetrahedra%tetra_mult(itetra))

!  Here we need the original ordering to reference the correct irred kpoints
   ind_dum(1) = tetrahedra%tetra_full(1,1,itetra)
   ind_dum(2) = tetrahedra%tetra_full(2,1,itetra)
   ind_dum(3) = tetrahedra%tetra_full(3,1,itetra)
   ind_dum(4) = tetrahedra%tetra_full(4,1,itetra)
   eigen_1tetra(1) = eigen_in(ind_dum(1))
   eigen_1tetra(2) = eigen_in(ind_dum(2))
   eigen_1tetra(3) = eigen_in(ind_dum(3))
   eigen_1tetra(4) = eigen_in(ind_dum(4))
   call sort_dp(4,eigen_1tetra,ind_dum,tol14)

!  
!  all notations are from Blochl PRB 49 16223 Appendix B
!  
   epsilon21 = eigen_1tetra(2)-eigen_1tetra(1)
   epsilon31 = eigen_1tetra(3)-eigen_1tetra(1)
   epsilon41 = eigen_1tetra(4)-eigen_1tetra(1)
   epsilon32 = eigen_1tetra(3)-eigen_1tetra(2)
   epsilon42 = eigen_1tetra(4)-eigen_1tetra(2)
   epsilon43 = eigen_1tetra(4)-eigen_1tetra(3)
   inv_epsilon21 = zero
   inv_epsilon31 = zero
   inv_epsilon41 = zero
   inv_epsilon32 = zero
   inv_epsilon42 = zero
   inv_epsilon43 = zero
   if (epsilon21 > tol6) then
     inv_epsilon21 = 1.d0 / epsilon21
   end if
   if (epsilon31 > tol6) then
     inv_epsilon31 = 1.d0 / epsilon31
   end if
   if (epsilon41 > tol6) then
     inv_epsilon41 = 1.d0 / epsilon41
   end if
   if (epsilon32 > tol6) then
     inv_epsilon32 = 1.d0 / epsilon32
   end if
   if (epsilon42 > tol6) then
     inv_epsilon42 = 1.d0 / epsilon42
   end if
   if (epsilon43 > tol6) then
     inv_epsilon43 = 1.d0 / epsilon43
   end if
   nn1 = int((eigen_1tetra(1)-enemin)/deltaene)+1
   nn2 = int((eigen_1tetra(2)-enemin)/deltaene)+1
   nn3 = int((eigen_1tetra(3)-enemin)/deltaene)+1
   nn4 = int((eigen_1tetra(4)-enemin)/deltaene)+1

   nn1 = max(1,nn1)
   nn1 = min(nn1,nene)
   nn2 = max(1,nn2)
   nn2 = min(nn2,nene)
   nn3 = max(1,nn3)
   nn3 = min(nn3,nene)
   nn4 = max(1,nn4)
   nn4 = min(nn4,nene)

   eps = enemin+nn1*deltaene
!  
!  interval enemin < eps < e1 nothing to do
!  
!  
!  interval e1 < eps < e2
!  
   deleps1 = eps-eigen_1tetra(1)
   cc_pre = volconst_mult*inv_epsilon21*inv_epsilon31*inv_epsilon41
   invepsum = inv_epsilon21+inv_epsilon31+inv_epsilon41
   dccde_pre = 3.d0*volconst_mult*inv_epsilon21*inv_epsilon31*inv_epsilon41
   do ieps=nn1+1,nn2
     cc = cc_pre * deleps1*deleps1*deleps1
     tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + cc*(4.d0-deleps1*invepsum)
     tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + cc*deleps1*inv_epsilon21
     tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + cc*deleps1*inv_epsilon31
     tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + cc*deleps1*inv_epsilon41

     dccde = dccde_pre * deleps1*deleps1
     dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) + &
&     dccde*(4.d0 - deleps1*invepsum) -cc*invepsum
     dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) + &
&     (dccde*deleps1 + cc) * inv_epsilon21
     dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) + &
&     (dccde*deleps1 + cc) * inv_epsilon31
     dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) + &
&     (dccde*deleps1 + cc) * inv_epsilon41

     deleps1 = deleps1 + deltaene
   end do
   eps = eps + (nn2-nn1)*deltaene
!  
!  interval e2 < eps < e3
!  
   deleps1 = eps-eigen_1tetra(1)
   deleps2 = eps-eigen_1tetra(2)
   deleps3 = eigen_1tetra(3)-eps
   deleps4 = eigen_1tetra(4)-eps

   cc1_pre = volconst_mult*inv_epsilon31*inv_epsilon41
   cc2_pre = volconst_mult*inv_epsilon41*inv_epsilon32*inv_epsilon31
   cc3_pre = volconst_mult*inv_epsilon42*inv_epsilon32*inv_epsilon41

   dcc1de_pre = 2.d0*cc1_pre
   dcc2de_pre = cc2_pre
   dcc3de_pre = cc3_pre
   do ieps=nn2+1,nn3
     cc1 = cc1_pre * deleps1*deleps1
     cc2 = cc2_pre * deleps1*deleps2*deleps3
     cc3 = cc3_pre * deleps2*deleps2*deleps4

     tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + &
&     cc1 + (cc1+cc2)*deleps3*inv_epsilon31 + &
&     (cc1+cc2+cc3)*deleps4*inv_epsilon41
     tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + &
&     cc1+cc2+cc3+(cc2+cc3)*deleps3*inv_epsilon32 +&
&     cc3*deleps4*inv_epsilon42
     tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + &
&     (cc1+cc2)*deleps1*inv_epsilon31 + &
&     (cc2+cc3)*deleps2*inv_epsilon32
     tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + &
&     (cc1+cc2+cc3)*deleps1*inv_epsilon41 + &
&     cc3*deleps2*inv_epsilon42

     dcc1de = dcc1de_pre * deleps1
     dcc2de = dcc2de_pre * (-deleps1*deleps2  +deleps1*deleps3  +deleps2*deleps3)
     dcc3de = dcc3de_pre * (2.d0*deleps2*deleps4  -deleps2*deleps2)

     dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) &
&     + dcc1de &
&     + ((dcc1de+dcc2de)*deleps3 -(cc1+cc2)) * inv_epsilon31 &
&     + ((dcc1de+dcc2de+dcc3de)*deleps4 -(cc1+cc2+cc3)) * inv_epsilon41
     dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) &
&     + dcc1de+dcc2de+dcc3de &
&     + ((dcc2de+dcc3de)*deleps3 -(cc2+cc3) ) * inv_epsilon32 &
&     + (dcc3de*deleps4  -cc3 ) * inv_epsilon42
     dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) &
&     + ((dcc1de+dcc2de)*deleps1 + (cc1+cc2) ) * inv_epsilon31 &
&     + ((dcc2de+dcc3de)*deleps2 + (cc2+cc3) ) * inv_epsilon32
     dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) &
&     + ((dcc1de+dcc2de+dcc3de)*deleps1 + (cc1+cc2+cc3) ) * inv_epsilon41 &
&     + (dcc3de*deleps2 + cc3) * inv_epsilon42

     deleps1 = deleps1 + deltaene
     deleps2 = deleps2 + deltaene
     deleps3 = deleps3 - deltaene
     deleps4 = deleps4 - deltaene
   end do
   eps = eps + (nn3-nn2)*deltaene
!  
!  interval e3 < eps < e4
!  
   deleps4 = eigen_1tetra(4)-eps
   cc_pre = volconst_mult*inv_epsilon41*inv_epsilon42*inv_epsilon43
   invepsum = inv_epsilon41+inv_epsilon42+inv_epsilon43
   dccde_pre = -3.d0*cc_pre
   do ieps=nn3+1,nn4
     cc = cc_pre * deleps4*deleps4*deleps4
     cc_tmp = cc * deleps4
     tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + volconst_mult - cc_tmp*inv_epsilon41
     tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + volconst_mult - cc_tmp*inv_epsilon42
     tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + volconst_mult - cc_tmp*inv_epsilon43
     tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + volconst_mult - cc*4.d0 + cc_tmp*invepsum

     dccde = dccde_pre * deleps4*deleps4
     dccde_tmp = -dccde*deleps4 + cc
     dtweightde_tmp(ieps,1) = dtweightde_tmp(ieps,1) + dccde_tmp * inv_epsilon41
     dtweightde_tmp(ieps,2) = dtweightde_tmp(ieps,2) + dccde_tmp * inv_epsilon42
     dtweightde_tmp(ieps,3) = dtweightde_tmp(ieps,3) + dccde_tmp * inv_epsilon43
     dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) - dccde*4.d0 - dccde_tmp*invepsum

     deleps4 = deleps4 - deltaene
   end do
   eps = eps + (nn4-nn3)*deltaene
!  
!  
!  interval e4 < eps < enemax
!  
   do ieps=nn4+1,nene
     tweight_tmp(ieps,1) = tweight_tmp(ieps,1) + volconst_mult
     tweight_tmp(ieps,2) = tweight_tmp(ieps,2) + volconst_mult
     tweight_tmp(ieps,3) = tweight_tmp(ieps,3) + volconst_mult
     tweight_tmp(ieps,4) = tweight_tmp(ieps,4) + volconst_mult
!    dtweightde unchanged by this tetrahedron
   end do

!  
!  if we have a fully degenerate tetrahedron,
!  1) the tweight is a Heaviside (step) function, which is correct above, but
!  2) the dtweightde should contain a Dirac function: add a Gaussian here
!  
   if (epsilon41 < tol6) then

!    to ensure the gaussian will integrate properly:
!    WARNING: this smearing could be problematic if too large
!    and doesnt integrate well if its too small
     gau_width = 10.0d0*deltaene
     gau_width2 = 1.0 / gau_width / gau_width
     gau_prefactor = volconst_mult / gau_width / sqrtpi
!    
!    average position since bracket for epsilon41 is relatively large
     cc = (eigen_1tetra(1)+eigen_1tetra(2)+eigen_1tetra(3)+eigen_1tetra(4))/4.d0
     eps = enemin
     do ieps=1,nene
       tmp = eps - cc
       dtweightde_tmp(ieps,4) = dtweightde_tmp(ieps,4) + &
&       gau_prefactor*exp(-tmp*tmp*gau_width2)
       eps = eps + deltaene
     end do
   end if
!  end degenerate tetrahedron if

! NOTE: the following blas calls are not working systematically, or do not give speed ups, strange...
!   call daxpy (nene, 1.d0,    tweight_tmp(:,1), 1,    tweight_t(:,ind_dum(1)), 1)
!   call daxpy (nene, 1.d0,    tweight_tmp(:,2), 1,    tweight_t(:,ind_dum(2)), 1)
!   call daxpy (nene, 1.d0,    tweight_tmp(:,3), 1,    tweight_t(:,ind_dum(3)), 1)
!   call daxpy (nene, 1.d0,    tweight_tmp(:,4), 1,    tweight_t(:,ind_dum(4)), 1)
!   call daxpy (nene, 1.d0, dtweightde_tmp(:,1), 1, dtweightde_t(:,ind_dum(1)), 1)
!   call daxpy (nene, 1.d0, dtweightde_tmp(:,2), 1, dtweightde_t(:,ind_dum(2)), 1)
!   call daxpy (nene, 1.d0, dtweightde_tmp(:,3), 1, dtweightde_t(:,ind_dum(3)), 1)
!   call daxpy (nene, 1.d0, dtweightde_tmp(:,4), 1, dtweightde_t(:,ind_dum(4)), 1)
   tweight_t(:,ind_dum(1)) = tweight_t(:,ind_dum(1)) + tweight_tmp(:,1)
   tweight_t(:,ind_dum(2)) = tweight_t(:,ind_dum(2)) + tweight_tmp(:,2)
   tweight_t(:,ind_dum(3)) = tweight_t(:,ind_dum(3)) + tweight_tmp(:,3)
   tweight_t(:,ind_dum(4)) = tweight_t(:,ind_dum(4)) + tweight_tmp(:,4)
   dtweightde_t(:,ind_dum(1)) = dtweightde_t(:,ind_dum(1)) + dtweightde_tmp(:,1)
   dtweightde_t(:,ind_dum(2)) = dtweightde_t(:,ind_dum(2)) + dtweightde_tmp(:,2)
   dtweightde_t(:,ind_dum(3)) = dtweightde_t(:,ind_dum(3)) + dtweightde_tmp(:,3)
   dtweightde_t(:,ind_dum(4)) = dtweightde_t(:,ind_dum(4)) + dtweightde_tmp(:,4)

 end do
!end do itetra

! transpose: otherwise the data access is crap and the code slows by an order of magnitude
 tweight    = transpose(tweight_t)
 dtweightde = transpose(dtweightde_t)

end subroutine get_tetra_weight
!!***



end module m_tetrahedron
!!***
