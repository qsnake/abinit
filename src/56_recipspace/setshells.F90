!{\src2tex{textfont=tt}}
!!****f* ABINIT/setshells
!! NAME
!! setshells
!!
!! FUNCTION
!! Set consistently the number of shells, the number of plane-waves,
!! and the energy cut-off
!!
!! COPYRIGHT
!! Copyright (C) 2001-2012 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nsym=number of symmetry operations
!!  gmet(3,3)=metric tensor in reciprocal space
!!  gprimd(3,3)=dimensional primitive vectors in reciprocal space  
!!  symrel(3,3,nsym)=symmetry operations in real space
!!  tag=suffix to account for the different possibilities for these variables (npw, ecut or nsh ..)
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  ecut,npw,nsh=one of them is an input, the two others are output
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  npw=number of plane waves
!!  nsh=number of shells
!!
!! PARENTS
!!      invars2m,rdm,setup_screening,setup_sigma
!!
!! CHILDREN
!!      initmpi_seq,kpgsph,sort_dp,wrtout
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setshells(ecut,npw,nsh,nsym,gmet,gprimd,symrel,tag,ucvol)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setshells'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace, except_this_one => setshells
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(inout) :: npw,nsh
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: ecut
 character(len=*),intent(in) :: tag
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: exchn2n3d,ifound,ig,ii,ish,istat,isym,npw_found,npwave
 integer :: npwwrk,nsh_found,pad=50
 real(dp) :: ecut_found,ecut_trial,eps,scale=1.3_dp
 logical :: found
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer :: geq(3)
 integer,allocatable :: gvec(:,:),gvec_sh(:,:),insort(:),npw_sh(:)
 real(dp) :: gctr(3)
 real(dp),allocatable :: gnorm(:),gnorm_sh(:)

!******************************************************************

 DBG_ENTER("COLL")
!
!=== Check coherence of input variables ecut, npw, and nsh ===
!1-> one at least should be non-null
 if (npw==0.and.nsh==0.and.ecut==0) then
   write(msg,'(8a)')&
&   ' One of the three variables ecut',TRIM(tag),', npw',TRIM(tag),', or nsh',TRIM(tag),ch10,&
&   ' must be non-null. Returning.'
   MSG_COMMENT(msg)
   RETURN
 end if
!2-> one and only one should be non-null
 if (npw/=0.and.nsh/=0) then
   write(msg,'(6a)')&
&   ' Only one of the two variables npw',TRIM(tag),' and nsh',TRIM(tag),ch10,&
&   ' can be non-null. Modify the value of one of these in input file.'
   MSG_ERROR(msg)
 end if
 if (ecut>tol6.and.npw/=0) then
   write(msg,'(6a)')&
&   ' Only one of the two variables ecut',TRIM(tag),' and npw',TRIM(tag),ch10,&
&   ' can be non-null. Modify the value of one of these in input file.'
   MSG_ERROR(msg)
 end if
 if (ecut>tol6.and.nsh/=0) then
   write(msg,'(6a)')&
&   ' Only one of the two variables ecut',TRIM(tag),' and nsh',TRIM(tag),ch10,&
&   ' can be non-null Action : modify the value of one of these in input file.'
   MSG_ERROR(msg)
 end if
!
!=== Calculates an upper bound for npw ===
!* gctr is center of the g-vector sphere
 gctr(:)=(/zero,zero,zero/)
 if (ecut>tol6) then
!  The average number of plane-waves in the cutoff sphere is given by: 
!  npwave = (2*ecut)**(3/2)*ucvol/(6*pi**2)
!  The upper bound is calculated as npwwrk=int(scale * npwave) + pad
   npwave=NINT(ucvol*(two*ecut)**1.5_dp/(six*pi**2))
   npwwrk=NINT(DBLE(npwave)*scale)+pad
   ecut_trial=ecut
 else if (npw/=0) then
!  npw is given in the input
   npwwrk=NINT(DBLE(npw)*scale)+pad
   ecut_trial=(six*pi**2*npw/ucvol)**two_thirds/two
 else
!  If nsh is given in the input
   npwwrk=nsh*18+2*pad
   ecut_trial=(six*pi**2*nsh*18/ucvol)**two_thirds/two
 end if

 ABI_ALLOCATE(gvec,(3,npwwrk))
 istat = ABI_ALLOC_STAT
 ifound=0
 do while(ifound==0)
   write(msg,'(a,f8.2)')' setshells : ecut_trial = ',ecut_trial

   call wrtout(std_out,msg,'COLL')
   exchn2n3d=0 ! For the time being, no exchange of n2 and n3

   call initmpi_seq(MPI_enreg_seq)

   call kpgsph(ecut_trial,exchn2n3d,gmet,0,1,1,gvec,gctr,1,MPI_enreg_seq,npwwrk,npw_found)

   ABI_ALLOCATE(gnorm,(npw_found))
   ABI_ALLOCATE(insort,(npw_found))

   do ig=1,npw_found
     insort(ig)=ig
     gnorm(ig)=zero
     do ii=1,3
       gnorm(ig)=gnorm(ig)+(gvec(1,ig)*gprimd(ii,1)+&
&       gvec(2,ig)*gprimd(ii,2)+&
&       gvec(3,ig)*gprimd(ii,3))**2
     end do
   end do
   call sort_dp(npw_found,gnorm,insort,tol14)
   ABI_ALLOCATE(npw_sh,(npw_found))
   ABI_ALLOCATE(gnorm_sh,(npw_found))
   ABI_ALLOCATE(gvec_sh,(3,npw_found))
   npw_sh(:)=0
   gnorm_sh(:)=zero
   gvec_sh(:,:)=0
!  Count the number of shells:
!  (search for the G-vectors generating the others by symmetry)
   nsh_found=0

   do ig=1,npw_found
     eps=1.d-8*gnorm(ig)
     found=.FALSE.
     ish=1
     do while ((.not.found).and.(ish<=nsh_found))
       if (ABS(gnorm(ig)-gnorm_sh(ish))<=eps) then
         isym=1
         do while ((.not.found).and.(isym<=nsym))
           geq(:)=(symrel(1,:,isym)*gvec(1,insort(ig))+&
&           symrel(2,:,isym)*gvec(2,insort(ig))+&
&           symrel(3,:,isym)*gvec(3,insort(ig)))

           found=((geq(1)==gvec_sh(1,ish)).and.&
&           (geq(2)==gvec_sh(2,ish)).and.&
&           (geq(3)==gvec_sh(3,ish)))
           isym=isym+1
         end do
       end if
       ish=ish+1
     end do
     if (.not.found) then
       nsh_found=nsh_found+1
       gnorm_sh(nsh_found)=gnorm(ig)
       gvec_sh(:,nsh_found)=gvec(:,insort(ig))
       npw_sh(nsh_found)=1
     else
       ish=ish-1
       npw_sh(ish)=npw_sh(ish)+1
     end if
   end do

   ecut_found=two*pi**2*gnorm(npw_found)

   if(ecut>tol6) then
!    ecut is given in the input
     if (ecut_found<ecut-0.1) then
       write(msg,'(3a,e14.6,9a,e14.6,3a)')&
&       ' The value ecut',TRIM(tag),'=',ecut,' given in the input file leads to',ch10,&
&       ' the same values for nsh',TRIM(tag),' and npw',TRIM(tag),' as ecut',TRIM(tag),'=',ecut_found,ch10,&
&       ' This value will be adopted for the calculation.',ch10
       MSG_WARNING(msg)
     end if
     ifound=1
   else if (npw/=0) then
!    If npw is given in the input
     if (npw_found==npw) then
       ecut_found=two*pi**2*gnorm(npw_found)
       ifound=1
     else if (npw_found>npw) then
       npw_found=0
       nsh_found=0
       do while (npw_found<npw)
         nsh_found=nsh_found+1
         npw_found=npw_found+npw_sh(nsh_found)
       end do
!      check that the shell is closed
       if(npw_found>npw) then
!        shell not closed
         npw_found=npw_found-npw_sh(nsh_found)
         nsh_found=nsh_found-1
         do while (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001)
           npw_found=npw_found-npw_sh(nsh_found)
           nsh_found=nsh_found-1
         end do
         write(msg,'(3a,i6,5a,i6,3a)')&
&         ' The value npw',TRIM(tag),'=',npw,' given in the input file does not close the shell',ch10,&
&         ' The lower closed-shell is obtained for a value npw',TRIM(tag),'=',npw_found,ch10,&
&         ' This value will be adopted for the calculation.',ch10
         MSG_WARNING(msg)
       end if
       ecut_found=two*pi**2*gnorm(npw_found)
       ifound=1
     end if
   else if (nsh/=0) then
!    If nsh is given in the input
     if (nsh_found==nsh) then
       ecut_found=two*pi**2*gnorm(npw_found)
       ifound=1
     else if (nsh_found>nsh) then
       npw_found=0
       nsh_found=0
       do ish=1,nsh
         npw_found=npw_found+npw_sh(ish)
         nsh_found=nsh_found+1
       end do
       if (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001) then
         do while (ABS(gnorm_sh(nsh_found)-gnorm_sh(nsh_found+1))<0.000001)
           nsh_found=nsh_found+1
           npw_found=npw_found+npw_sh(nsh_found)
         end do
         write(msg,'(3a,i6,5a,i6,3a)')&
&         ' The value nsh',TRIM(tag),'=',nsh,' given in the input file corresponds to the same',ch10,&
&         ' cut-off energy as for closed-shell upto nsh',TRIM(tag),'=',nsh_found,ch10,&
&         ' This value will be adopted for the calculation.',ch10
         MSG_WARNING(msg)
       end if
       ecut_found=two*pi**2*gnorm(npw_found)
       ifound=1
     end if
   end if

   if (ifound==0) then
     ecut_trial=1.1*ecut_trial
     ABI_DEALLOCATE(gnorm)
     ABI_DEALLOCATE(gnorm_sh)
     ABI_DEALLOCATE(gvec_sh)
     ABI_DEALLOCATE(insort)
     ABI_DEALLOCATE(npw_sh)
   else
     ecut=ecut_found
     npw=npw_found
     nsh=nsh_found
   end if
 end do !while(ifound==0)

 ABI_DEALLOCATE(gnorm)
 ABI_DEALLOCATE(gnorm_sh)
 ABI_DEALLOCATE(gvec)
 ABI_DEALLOCATE(gvec_sh)
 ABI_DEALLOCATE(insort)
 ABI_DEALLOCATE(npw_sh)

 DBG_EXIT("COLL")
 
end subroutine setshells
!!***
