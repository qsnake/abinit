!{\src2tex{textfont=tt}}
!!****f* ABINIT/complete_gkk
!!
!! NAME
!! complete_gkk
!!
!! FUNCTION
!! Use the set of special q points calculated by the Monkhorst &
!! Pack Technique.
!! Check if all the informations for the q points are present in
!! the DDB to determine the elphon interaction matrices
!! Generate the gkk matrices of the set of q points which
!! samples homogeneously the entire Brillouin zone.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! elph_ds = datastructure for elphon information (mainly
!!      matrix elements and dimensions)
!!   elph_ds%k_phon%full2full = kpt_phon index mapping under symops
!! gkk_flag = flag for existence of matrix element
!! gprimd(3,3)=dimensionful primitive translations in reciprocal space
!! indsym = map of atoms by inverses of symrels
!! natom=number of atoms in unit cell
!! nsym=number of space group symmetries
!! qpttoqpt = qpoint index mapping under symops
!! rprimd(3,3)=dimensionful primitive translations in real space
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (recip space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! tnons(3,nsym)=nonsymmorphic translations associated to symrel
!!
!! OUTPUT
!! elph_ds%gkk_qpt = gkk matrices for all qpts on a full mesh
!!
!! TODO
!!
!! PARENTS
!!      get_all_gkq
!!
!! CHILDREN
!!      zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine complete_gkk(elph_ds,gkk_flag,&
&   gprimd,indsym,natom,nsym,qpttoqpt,rprimd,&
&   symrec,symrel)

 use m_profiling

 use defs_basis
 use defs_elphon

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'complete_gkk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: indsym(4,nsym,natom)
 integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt_full),symrec(3,3,nsym)
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol,elph_ds%nqpt_full)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ib1,ibranch,ieqqpt,ii
 integer :: iqpt,isppol,isym
 integer :: itim,jbranch,jj,kk,ll
 integer :: neqqpt,symikpt_phon
 integer :: iatom,ancestor_iatom

 real(dp),parameter :: tol=2.d-8
!arrays
 integer :: symmetrized_qpt(elph_ds%nqpt_full)
 real(dp) :: ss(3,3)
 real(dp) :: tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),allocatable :: gkk_qpt_new(:,:,:,:,:),gkk_qpt_tmp(:,:,:,:,:)

 real(dp) :: ss_allatoms(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: c_one(2), c_zero(2)


! *********************************************************************


 c_one = (/one,zero/)
 c_zero = (/zero,zero/)

!Generation of the gkk matrices relative to the q points
!of the set which samples the entire Brillouin zone

 symmetrized_qpt(:) = -1

 isppol=1

 ABI_ALLOCATE(gkk_qpt_new,(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol))
 ABI_ALLOCATE(gkk_qpt_tmp,(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nsppol))

 do iqpt=1,elph_ds%nqpt_full

!  Already symmetrized?
   if (symmetrized_qpt(iqpt) == 1) cycle

   gkk_qpt_new(:,:,:,:,:) = zero

!  loop over qpoints equivalent to iqpt
   neqqpt=0
!  do not use time reversal symmetry to complete the qpoints:
!  do not know what happens to the gamma matrices
   itim=1

!  do itim=1,2
   do isym=1,nsym
!    ieqqpt is sent onto iqpt by itim/isym
     ieqqpt = qpttoqpt(itim,isym,iqpt)


     if (gkk_flag(1,1,1,1,ieqqpt) == -1) cycle
!    if we have information on this qpt

!    iqpt is equivalent to ieqqpt: get it from file or memory
     if (elph_ds%gkqwrite == 0) then
       gkk_qpt_tmp(:,:,:,:,:) = elph_ds%gkk_qpt(:,:,:,:,:,ieqqpt)
     else
       do ikpt_phon=1, elph_ds%k_phon%nkpt
         read(elph_ds%unitgkq,REC=((ieqqpt-1)*elph_ds%k_phon%nkpt+ikpt_phon)) gkk_qpt_tmp(:,:,:,ikpt_phon,:)
       end do
     end if


     neqqpt=neqqpt+1

     if (elph_ds%ep_scalprod==1) then
       do ii=1,3
         do jj=1,3
           ss(ii,jj)=0.0_dp
           do kk=1,3
             do ll=1,3
               ss(ii,jj)=ss(ii,jj)+rprimd(ii,kk)*symrel(kk,ll,isym)*gprimd(ll,jj)
             end do
           end do
         end do
       end do
     else
       do ii=1,3
         do jj=1,3
           ss(ii,jj) = symrec(jj,ii,isym)
         end do
       end do
     end if

     ss_allatoms(:,:,:) = zero
     do iatom=1,natom
       ancestor_iatom = indsym(4,isym,iatom)
!      do jatom=1,natom
!      ancestor_jatom = indsym(4,isym,jatom)
       ss_allatoms(1,(ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3,&
&       (iatom-1)*3+1:         (iatom-1)*3+3) = ss(1:3,1:3)
!      end do
     end do


!    NOTE   ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*symrec(ll,kk,isym)

     do isppol=1,elph_ds%nsppol
       do ikpt_phon=1,elph_ds%k_phon%nkpt
!        symikpt_phon is sent onto ikpt_phon by itim/isym
         symikpt_phon=elph_ds%k_phon%full2full(itim,isym,ikpt_phon)

!        Do each element band1, band2 separately...
         do ib1=1,elph_ds%ngkkband*elph_ds%ngkkband

!          multiply by the ss matrices
           tmp_mat2(:,:,:) = zero
           tmp_mat(:,:,:) = reshape(gkk_qpt_tmp(:,ib1,:,ikpt_phon,isppol),&
&           (/2,elph_ds%nbranch,elph_ds%nbranch/))
           call ZGEMM ('N','N',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,&
&           c_one,ss_allatoms,elph_ds%nbranch,tmp_mat,elph_ds%nbranch,c_zero,&
&           tmp_mat2,elph_ds%nbranch)
           call ZGEMM ('N','T',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,&
&           c_one,tmp_mat2,elph_ds%nbranch,ss_allatoms,elph_ds%nbranch,c_zero,&
&           tmp_mat,elph_ds%nbranch)

!          add to gkk_qpt_new
           do ibranch =1,elph_ds%nbranch
             do jbranch =1,elph_ds%nbranch
               gkk_qpt_new(:,ib1,(jbranch-1)*elph_ds%nbranch+ibranch,symikpt_phon,isppol) = &
&               gkk_qpt_new(:,ib1,(jbranch-1)*elph_ds%nbranch+ibranch,symikpt_phon,isppol) + &
&               tmp_mat(:,jbranch,ibranch)
             end do
           end do

         end do ! end ib1 do
       end do ! end ikpt_phon do
     end do ! end isppol do

   end do ! end isym do

   if (neqqpt > 1) then
     write(std_out,*) ' found several equiv qpts and am symmetrizing them ', neqqpt
   end if

!  divide by number of equivalent qpts found
   gkk_qpt_new(:,:,:,:,:) = gkk_qpt_new(:,:,:,:,:)/neqqpt

!  copy the symmetrized version into all the equivalent qpoints, appropriately transformed
!  See above
   itim=1
!  do itim=1,2
   do isym=1,nsym
!    ieqqpt is sent onto iqpt by itim/isym
     ieqqpt = qpttoqpt(itim,isym,iqpt)

     if (symmetrized_qpt(ieqqpt) /= -1) cycle
     gkk_qpt_tmp(:,:,:,:,:) = zero

!    use symrec matrices to get inverse transform from isym^{-1}
     if (elph_ds%ep_scalprod==1) then
       do ii=1,3
         do jj=1,3
           ss(ii,jj)=0.0_dp
           do kk=1,3
             do ll=1,3
!              Use inverse of symop matrix here to get back to ieqqpt (inv+transpose is in symrec and in gprimd)
               ss(ii,jj)=ss(ii,jj)+rprimd(ii,kk)*symrec(ll,kk,isym)*gprimd(ll,jj)
             end do
           end do
         end do
       end do
     else
       do ii=1,3
         do jj=1,3
           ss(ii,jj) = symrel(ii,jj,isym)
         end do
       end do
     end if

     ss_allatoms(:,:,:) = zero
     do iatom=1,natom
       ancestor_iatom = indsym(4,isym,iatom)
!      do jatom=1,natom
!      ancestor_jatom = indsym(4,isym,jatom)
       ss_allatoms(1,(ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3,&
&       (iatom-1)*3+1:          (iatom-1)*3+3) = ss(1:3,1:3)
!      end do
     end do

!    ! Use inverse of symop matrix here to get back to ieqqpt
!    ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*symrel(kk,ll,isym)

     do isppol=1,elph_ds%nsppol
       do ikpt_phon=1,elph_ds%k_phon%nkpt
!        symikpt_phon is sent onto ikpt_phon by itim/isym
         symikpt_phon=elph_ds%k_phon%full2full(itim,isym,ikpt_phon)

         do ib1=1,elph_ds%ngkkband*elph_ds%ngkkband

!          multiply by the ss^{-1} matrices
           tmp_mat2(:,:,:) = zero
           tmp_mat(:,:,:) = reshape(gkk_qpt_new(:,ib1,:,ikpt_phon,isppol),&
&           (/2,elph_ds%nbranch,elph_ds%nbranch/))
           call ZGEMM ('N','N',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,&
&           c_one,ss_allatoms,elph_ds%nbranch,tmp_mat,elph_ds%nbranch,c_zero,&
&           tmp_mat2,elph_ds%nbranch)
           call ZGEMM ('N','T',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,&
&           c_one,tmp_mat2,elph_ds%nbranch,ss_allatoms,elph_ds%nbranch,c_zero,&
&           tmp_mat,elph_ds%nbranch)

           do ibranch =1,elph_ds%nbranch
             do jbranch =1,elph_ds%nbranch
               gkk_qpt_tmp(:,ib1,(jbranch-1)*elph_ds%nbranch+ibranch,symikpt_phon,isppol) =&
&               tmp_mat(:,jbranch,ibranch)
             end do
           end do

           if (gkk_flag (1,1,symikpt_phon,isppol,ieqqpt) == -1) then
             gkk_flag (:,:,symikpt_phon,isppol,ieqqpt) = 0
           end if

         end do ! end ib1 do
       end do ! end ikpt_phon do
     end do ! end isppol do


!    save symmetrized matrices for qpt ieqqpt
     if (elph_ds%gkqwrite == 0) then
       elph_ds%gkk_qpt(:,:,:,:,:,ieqqpt) = gkk_qpt_tmp(:,:,:,:,:)
     else
       do ikpt_phon=1, elph_ds%k_phon%nkpt
         write(elph_ds%unitgkq,REC=((ieqqpt-1)*elph_ds%k_phon%nkpt+ikpt_phon)) gkk_qpt_tmp(:,:,:,ikpt_phon,:)
       end do
     end if

     symmetrized_qpt(ieqqpt) = 1

   end do ! end isym do 
!  end do ! end itim do

 end do
!end iqpt do

end subroutine complete_gkk
!!***
