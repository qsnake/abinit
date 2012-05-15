!{\src2tex{textfont=tt}}
!!****f* ABINIT/integrate_gamma_tr
!!
!! NAME
!! integrate_gamma_tr
!!
!! FUNCTION
!! This routine integrates the TRANSPORT electron phonon coupling matrices
!! over the kpoints on the fermi surface. A dependency on qpoint
!! remains for gamma_qpt_in/out
!! Copied from integrate_gamma
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (JPC,MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!      elph_ds%qpt_full = qpoint coordinates
!!   FSfullpqtofull = mapping of k+q to k
!!   nrpt = number of real space points for FT
!!
!! OUTPUT
!!   elph_tr_ds%gamma_qpt_trout and created elph_tr_ds%gamma_rpt_trout
!!   elph_tr_ds%gamma_qpt_trin and created elph_tr_ds%gamma_rpt_trin
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine integrate_gamma_tr(elph_ds,FSfullpqtofull,nrpt,elph_tr_ds)

 use m_profiling

 use defs_basis
  use defs_datatypes
 use defs_abitypes
  use defs_elphon

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'integrate_gamma_tr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nrpt
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ierr,iqpt,iqpt_fullbz,isppol
 integer :: itensor, icomp, jcomp
 integer :: fib1, fib2
! integer :: ikpttemp
 real(dp) :: etain, etaout
 character(len=500) :: message
!arrays
 real(dp) :: elvelock(3), elvelockpq(3)
 real(dp),allocatable :: tmp_gkk(:,:,:,:)

! *************************************************************************

 ABI_ALLOCATE(elph_tr_ds%gamma_qpt_trin,(2,9,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqpt_full))
 ierr = ABI_ALLOC_STAT
 if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&   ' trying to allocate array elph_tr_ds%gamma_qpt_trin '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 elph_tr_ds%gamma_qpt_trin = zero

 ABI_ALLOCATE(elph_tr_ds%gamma_qpt_trout,(2,9,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqpt_full))
 ierr = ABI_ALLOC_STAT
 if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&   ' trying to allocate array elph_tr_ds%gamma_qpt_trout '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 elph_tr_ds%gamma_qpt_trout = zero

 ABI_ALLOCATE(elph_tr_ds%gamma_rpt_trout,(2,9,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,nrpt))
 ierr = ABI_ALLOC_STAT
 if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&   ' trying to allocate array elph_tr_ds%gamma_rpt_trout '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 elph_tr_ds%gamma_rpt_trout = zero

 ABI_ALLOCATE(elph_tr_ds%gamma_rpt_trin,(2,9,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,nrpt))
 ierr = ABI_ALLOC_STAT
 if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&   ' trying to allocate array elph_tr_ds%gamma_rpt_trin '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 elph_tr_ds%gamma_rpt_trin = zero


!information
 if (elph_ds%gkqwrite == 0) then
   write (message,'(a)')' integrate_gamma_tr : keeping gamma matrices in memory'
   call wrtout(std_out,message,'COLL')
 else if (elph_ds%gkqwrite == 1) then
   write (message,'(a)')' integrate_gamma_tr : reading gamma matrices from disk'
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(3a,i3)')' integrate_gamma_tr : BUG-',ch10,&
&   ' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!allocate temp variables
 ABI_ALLOCATE(tmp_gkk ,(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol))
 ierr = ABI_ALLOC_STAT
 if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&   ' trying to allocate array tmp_gkkout '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 do iqpt=1,elph_ds%nqptirred
   iqpt_fullbz = elph_ds%qirredtofull(iqpt)
   write(std_out,*)'iqpt, iqptfullbz  ',iqpt, iqpt_fullbz

   do ikpt_phon=1,elph_ds%k_phon%nkpt

     if (elph_ds%gkqwrite == 0) then
       tmp_gkk = elph_ds%gkk_qpt(:,:,:,ikpt_phon,:,iqpt)
     else if (elph_ds%gkqwrite == 1) then
       read(elph_ds%unitgkq,REC=((iqpt-1)*elph_ds%k_phon%nkpt+ikpt_phon)) tmp_gkk
     end if

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     do isppol=1,elph_ds%nsppol
       do ib1=1,elph_ds%ngkkband
         fib1=ib1+elph_ds%minFSband-1
         elvelock(:)=elph_tr_ds%el_veloc(ikpt_phon,fib1,:,isppol)

         do ib2=1,elph_ds%ngkkband
           ibeff=ib2+(ib1-1)*elph_ds%ngkkband
           fib2=ib2+elph_ds%minFSband-1
           elvelockpq(:)= elph_tr_ds%el_veloc(ikpt_phonq,fib2,:,isppol)


!          MJV 31/03/2009: Note that the following is valid for any geometry, not just cubic!
!          see eq 5 and 6 of prb 36 4103 (Al-Lehaibi et al 1987)
!          see also Allen PRB 17 3725
!          generalization to tensorial quantities is simple, by keeping the directional
!          references of velock and velockpq as indices.
           do icomp = 1, 3
             do jcomp = 1, 3
               itensor = (icomp-1)*3+jcomp
!              FIXME: could use symmetry i <-> j

               etain  = elvelock(icomp)*elvelockpq(jcomp)
               etaout = elvelock(icomp)*elvelock(jcomp)


               elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,iqpt_fullbz) = &
&               elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,iqpt_fullbz) + &
&               tmp_gkk(:,ibeff,:,isppol) &
&               *etain &
&               *elph_ds%gkk_intweight(ib1,ikpt_phon,isppol)*elph_ds%gkk_intweight(ib2,ikpt_phonq,isppol)
               
               elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,iqpt_fullbz) = &
&               elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,iqpt_fullbz) + &
&               tmp_gkk(:,ibeff,:,isppol) &
&               *etaout &
&               *elph_ds%gkk_intweight(ib1,ikpt_phon,isppol)*elph_ds%gkk_intweight(ib2,ikpt_phonq,isppol)
             end do
           end do
         end do
       end do

     end do ! isppol
   end do ! ik

 end do ! iq

 ABI_DEALLOCATE(tmp_gkk)

!
!normalize tensor with 1/sqrt(v_x**2 * v_y**2)
!
 do isppol=1, elph_ds%nsppol
   do icomp = 1, 3
     do jcomp = 1, 3
       itensor = (icomp-1)*3+jcomp
       elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,:) = elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,:) / &
&       sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))
       elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,:) = elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,:) / &
&       sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))
     end do 
   end do 
 end do ! isppol


!need prefactor of 1/nkpt for each integration over 1 kpoint index.
!NOT INCLUDED IN elph_ds%gkk_intweight
 elph_tr_ds%gamma_qpt_trout = elph_tr_ds%gamma_qpt_trout* elph_ds%occ_factor / elph_ds%k_phon%nkpt
 elph_tr_ds%gamma_qpt_trin  = elph_tr_ds%gamma_qpt_trin * elph_ds%occ_factor / elph_ds%k_phon%nkpt

 write (message,'(2a)')' integrate_gamma_tr : transport gamma matrices are calculated ',&
& ' in recip space and for irred qpoints'
 call wrtout(std_out,message,'COLL')

end subroutine integrate_gamma_tr
!!***
