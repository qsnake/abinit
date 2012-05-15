!{\src2tex{textfont=tt}}
!!****f* ABINIT/integrate_gamma_alt
!!
!! NAME
!! integrate_gamma_alt
!!
!! FUNCTION
!! This routine integrates the electron phonon coupling matrix
!! over the kpoints on the fermi surface. A dependency on qpoint
!! remains for gamma_qpt
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!      elph_ds%qpt_full = qpoint coordinates
!!      elph_ds%nqptirred = number of irred qpoints
!!      elph_ds%qirredtofull = indexing of the GKK qpoints found
!!   elph_tr_ds = elphon transport datastructure
!!      elph_tr_ds
!!   nrpt = number of real space points for FT
!!
!! OUTPUT
!!   elph_ds = modified elph_ds%gamma_qpt with much larger size
!!   elph_tr_ds = modified elph_ds%gamma_qpt with much larger size
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      complete_gamma,ftgam,get_rank_1kpt,interpol3d,k_neighbors,leave_new
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine integrate_gamma_alt(elph_ds,elph_tr_ds,Cryst,gprim,kptrlatt,&
&   natom,nrpt,nsym,qpttoqpt,rpt,wghatm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_errors
 use m_kptrank

 use m_crystal,    only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'integrate_gamma_alt'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_77_ddb, except_this_one => integrate_gamma_alt
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nrpt
 integer, intent(in) :: natom, nsym
 
 type(crystal_structure),intent(in) :: Cryst
 type(elph_type),intent(inout) :: elph_ds
 type(elph_tr_type), intent(inout) :: elph_tr_ds
!arrays
 integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt_full)
 real(dp), intent(in) :: wghatm(natom,natom,nrpt), rpt(3,nrpt)
 real(dp), intent(in) :: gprim(3,3)
 integer, intent(in) :: kptrlatt(3,3)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ib1,ib2,ibeff,ierr,iqpt,iqpt_fullbz,isppol
 integer :: irec, qtor, rtoq, ibranch
 integer :: symrankkpt
 integer :: ikptfine,ikptqfine
 integer :: ineighbor
 integer :: counter_fine
 integer :: fib1, fib2
 integer :: itensor, icomp, jcomp

 real(dp) :: sd1, sd2, tolfs

!arrays
 integer :: kpt_phon_indices (8)
 real(dp) :: etaout, etain
 real(dp) :: qpt(3), rel_kpt(3)
 real(dp) :: kpt(3)
 real(dp) :: gam_now(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: elvelock(3), elvelockpq(3)

 character(len=500) :: message
 character(len=fnlen) :: fname

 real(dp),allocatable :: tmp_gkk(:,:,:,:,:)
 real(dp),allocatable :: tmp_gkk_1q(:,:,:,:,:)
 real(dp),allocatable :: tmp_gkr(:,:,:)

! *************************************************************************

 write (message,'(3a)')ch10,' entering integrate_gamma_alt ',ch10
 call wrtout(std_out,message,'COLL')

 tolfs = 1.e-11

 if (elph_ds%ep_keepbands == 0) then
   write (message,'(4a)')ch10,' ERROR: cannot call integrate_gamma_alt without all bands in memory',ch10,&
&   ' ACTION: set ep_keepbands = 1 in anaddb input file'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!allocations
 ABI_ALLOCATE(elph_ds%gamma_qpt,(2,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%k_fine%nkpt))
 ierr = ABI_ALLOC_STAT
 if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma_alt : ERROR- ',ch10,&
&   ' trying to allocate array elph_ds%gamma_qpt '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 elph_ds%gamma_qpt(:,:,:,:) = zero

!transport case
 if (elph_tr_ds%ifltransport==1) then
   ABI_ALLOCATE(elph_tr_ds%gamma_qpt_trin,(2,9,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%k_fine%nkpt))
   ierr = ABI_ALLOC_STAT
   if (ierr /= 0 ) then
     write (message,'(3a)')' integrate_gamma_alt : ERROR- ',ch10,&
&     ' trying to allocate array elph_tr_ds%gamma_qpt_trin '
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if 
   elph_tr_ds%gamma_qpt_trin = zero
   
   ABI_ALLOCATE(elph_tr_ds%gamma_qpt_trout,(2,9,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%k_fine%nkpt))
   ierr = ABI_ALLOC_STAT
   if (ierr /= 0 ) then
     write (message,'(3a)')' integrate_gamma_alt : ERROR- ',ch10,&
&     ' trying to allocate array elph_tr_ds%gamma_qpt_trout '
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if 
   elph_tr_ds%gamma_qpt_trout = zero
 end if 


!temp variables
 ABI_ALLOCATE(tmp_gkk_1q ,(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,8,elph_ds%nsppol))
 ierr = ABI_ALLOC_STAT
 
 ABI_CHECK(ierr==0,'trying to allocate tmp_gkk_1q')

 ABI_ALLOCATE(tmp_gkk ,(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqpt_full))
 ierr = ABI_ALLOC_STAT

 ABI_CHECK(ierr==0,' trying to allocate tmp_gkk')

 ABI_ALLOCATE(tmp_gkr ,(2,elph_ds%nbranch*elph_ds%nbranch,nrpt))
 ierr = ABI_ALLOC_STAT

 ABI_CHECK(ierr==0,' trying to allocate array tmp_gkr')

 if (elph_ds%gkqwrite == 0) then
   call wrtout(std_out,' integrate_gamma_alt : getting gamma matrices from memory','COLL')
 else if (elph_ds%gkqwrite == 1) then
   fname=trim(elph_ds%elph_base_name) // '_GKKQ'
   write (message,'(2a)')' integrate_gamma_alt : reading gamma matrices from file ',trim(fname)
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(a,i0)')' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if

 qtor = 1 ! q --> r
 rtoq = 0 ! r --> q

 counter_fine = 0

!FIXME: group fine kpt in batches with the same corners in sparse k-grid, to do only 1
!retrieval of gkk for each batch.
 do ikptfine=1,elph_ds%k_fine%nkpt
!  if we have found a good k-point on the FS
   if (sum(abs(elph_ds%k_fine%wtk(:,ikptfine,:))) < tolfs) cycle
   counter_fine = counter_fine+1

   call k_neighbors (elph_ds%k_fine%kpt(:,ikptfine), kptrlatt, elph_ds%k_phon%kptrank_t, &
&   rel_kpt, kpt_phon_indices)

!  FIXME: this read in of 8 matrices is done several times for fine kpt inside the same phon_kpt cell
!  more generally could keep 4 matrices when moving to the neighboring cell, but bookkeeping
!  will be a pain.
   tmp_gkk = zero
   do iqpt=1,elph_ds%nqptirred
     do ineighbor = 1, 8
       ikpt_phon = kpt_phon_indices(ineighbor)
       if (elph_ds%gkqwrite == 0) then
         tmp_gkk_1q(:,:,:,ineighbor,:) = elph_ds%gkk_qpt(:,:,:,ikpt_phon,:,iqpt)
       else if (elph_ds%gkqwrite == 1) then
         irec = (iqpt-1)*elph_ds%k_phon%nkpt+ikpt_phon
         read (elph_ds%unitgkq,REC=irec) tmp_gkk_1q(:,:,:,ineighbor,:)
       end if
     end do ! neighbor ikpt_phon

!    interpolate gkk matrix elements wrt k vector
!    FIXME: make a batch version of interpol3d for all matrix elements at once
     do isppol = 1, elph_ds%nsppol
       do ibranch = 1, elph_ds%nbranch*elph_ds%nbranch
         do ibeff = 1, elph_ds%ngkkband*elph_ds%ngkkband
           call interpol3d(rel_kpt,2,2,2,tmp_gkk(1,ibeff,ibranch,isppol,elph_ds%qirredtofull(iqpt)),&
&           tmp_gkk_1q(1,ibeff,ibranch,:,isppol))
           call interpol3d(rel_kpt,2,2,2,tmp_gkk(2,ibeff,ibranch,isppol,elph_ds%qirredtofull(iqpt)),&
&           tmp_gkk_1q(2,ibeff,ibranch,:,isppol))
         end do
       end do
     end do
   end do ! do iqpt
!  now tmpgkk is filled for present fine kpt and all irred qpoints on sparse grid


!  find spin and band for FS contributing state
   do isppol=1,elph_ds%nsppol
     do ib1=1,elph_ds%ngkkband
       if (abs(elph_ds%k_fine%wtk(ib1,ikptfine,isppol)) < tolfs) cycle
       sd1 = elph_ds%k_fine%wtk(ib1,ikptfine,isppol)

       if (elph_tr_ds%ifltransport==1) then
         fib1 = ib1 + elph_ds%minFSband - 1
         elvelock(:)=elph_tr_ds%el_veloc(ikptfine,fib1,:,isppol)
       end if

!      loop over k+q, n' point
!      FIXME: check if we need a second spin index here... looks like there is none in abinit matrix elements...
!      does this presume no SO operator in the perturbed potential? Otherwise the deltaV is spin diagonal...
       do ib2=1,elph_ds%ngkkband
         ibeff = ib2+(ib1-1)*elph_ds%ngkkband

!        do FT of gkq matrix to gkR
         if (elph_ds%symgkq ==1) then
!          FIXME: this does both sppol by default, which is a waste here inside the loop on isppol
           call complete_gamma(Cryst,elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqptirred,elph_ds%nqpt_full,&
&           elph_ds%ep_scalprod,elph_ds%qirredtofull,qpttoqpt,tmp_gkk(:,ibeff,:,:,:))
         end if

!        Now FT to real space too
         call ftgam(wghatm,tmp_gkk(:,ibeff,:,isppol,:),tmp_gkr(:,:,:),gprim,natom,&
&         elph_ds%nqpt_full,nrpt,qtor,rpt,elph_ds%qpt_full)

!        run _qpt_ over the dense grid of kpt_fine
!        FIXME: could at least use time reversal sym, and probably more, to do fewer loop iterations here
         do iqpt_fullbz = 1, elph_ds%k_fine%nkpt

!          find bloody index of        k             +             q
           kpt(:) = elph_ds%k_fine%kpt(:,ikptfine) + elph_ds%k_fine%kpt(:,iqpt_fullbz)
!          NB: this indexing is not the same order as the normal full qpts (elph_ds%k_fine%kpt(:,ikpt_fine) + qpt(:,iqpt_fullbz)
!          and it runs over all kpt_fine, not just the sparse input qpts.

!          which kpt is k+q among the full FS kpts
           call get_rank_1kpt (kpt,symrankkpt,elph_ds%k_fine%kptrank_t)
           ikptqfine = elph_ds%k_fine%kptrank_t%invrank(symrankkpt)

           if (elph_tr_ds%ifltransport==1) then
             fib2 = ib2 + elph_ds%minFSband - 1
             elvelockpq(:)=elph_tr_ds%el_veloc(ikptqfine,fib2,:,isppol)
           end if


           if (abs(elph_ds%k_fine%wtk(ib2,ikptqfine,isppol)) < tolfs) cycle
           sd2 = elph_ds%k_fine%wtk(ib2,ikptqfine,isppol)
           qpt = elph_ds%k_fine%kpt(:,iqpt_fullbz)

!          interpolate gkR to required gkq
           call ftgam(wghatm,gam_now,tmp_gkr(:,:,:),gprim,natom,1,nrpt,rtoq,rpt,qpt)

!          add to gamma(q) = sum over k and bands
           elph_ds%gamma_qpt(:,:,isppol,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,isppol,iqpt_fullbz) &
&           + gam_now*sd1*sd2

           if (elph_tr_ds%ifltransport==1) then
             do icomp = 1, 3
               do jcomp = 1, 3
                 itensor = (icomp-1)*3+jcomp
!                FIXME: could use symmetry i <-> j
                 
                 etain  = elvelock(icomp)*elvelockpq(jcomp)
                 etaout = elvelock(icomp)*elvelock(jcomp)

                 elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,iqpt_fullbz) = &
&                 elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,iqpt_fullbz) + &
&                 gam_now(:,:)*etain*sd1*sd2
                 
                 elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,iqpt_fullbz) = &
&                 elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,iqpt_fullbz) + &
&                 gam_now(:,:)*etaout*sd1*sd2
               end do ! jcomp
             end do ! icomp
           end if

         end do !iqpt_fullbz
       end do !ib2
     end do !ib1
   end do !isppol
 end do !ikpt_fine

 ABI_DEALLOCATE(tmp_gkk)
 ABI_DEALLOCATE(tmp_gkr)

 write (message,'(a,E20.10)') 'fraction of fine kpt used : ', dble(counter_fine) / elph_ds%k_fine%nkpt
 call wrtout(std_out,message,'COLL')

!need prefactor of 1/nkpt for each integration over 1 kpoint index.
!NOT INCLUDED IN elph_ds%k_fine%wtk
 elph_ds%gamma_qpt          = elph_ds%gamma_qpt          * elph_ds%occ_factor / elph_ds%k_fine%nkpt
 elph_tr_ds%gamma_qpt_trin  = elph_tr_ds%gamma_qpt_trin  * elph_ds%occ_factor / elph_ds%k_fine%nkpt
 elph_tr_ds%gamma_qpt_trout = elph_tr_ds%gamma_qpt_trout * elph_ds%occ_factor / elph_ds%k_fine%nkpt

!normalization in transport case
 if (elph_tr_ds%ifltransport==1) then
   do isppol = 1, elph_ds%nsppol
     do icomp = 1, 3
       do jcomp = 1, 3
         itensor = (icomp-1)*3+jcomp
         elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,:) = elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,:) / &
&         sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))
         elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,:) = elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,:) / &
&         sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))
       end do
     end do
   end do
 end if


 write (message,'(a,a)') ' integrate_gamma_alt : gamma matrices have been calculated',&
& 'for recip space and irred qpoints '
 call wrtout(std_out,message,'COLL')

end subroutine integrate_gamma_alt
!!***

