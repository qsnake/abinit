!{\src2tex{textfont=tt}}
!!****f* ABINIT/listkk
!! NAME
!! listkk
!!
!! FUNCTION
!! Given a list of nkpt1 initial k points kptns1 and a list of nkpt2
!! final k points kptns2, associates each final k pt with a "closest"
!! initial k point as determined by a (partly arbitrary) metric gmet.
!! Returns indirect indexing list indkk.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  kptns1(3,nkpt1)=list of initial k points (reduced coordinates)
!!  kptns2(3,nkpt2)=list of final k points
!!  nkpt1=number of initial k points
!!  nkpt2=number of final k points
!!  nsym=number of symmetry elements in space group
!!  sppoldbl=if 1, no spin-polarisation doubling
!!           if 2, spin-polarisation doubling using symafm
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symmat(3,3,nsym)=symmetry operations (symrel or symrec, depending on 
!!                   value of use_symrec
!!  timrev=1 if the use of time-reversal is allowed; 0 otherwise
!!  use_symrec :: if present and true, symmat assumed to be symrec, otherwise 
!!                assumed to be symrel
!!
!! OUTPUT
!!  dksqmax=maximal value of the norm**2 of the difference between
!!   a kpt2 vector and the closest k-point found from the kptns1 set,
!!   using symmetries.
!!  indkk(nkpt2*sppoldbl,6)=describe k point number of kpt1 that allows to
!!   generate wavefunctions closest to given kpt2
!!   if sppoldbl=2, use symafm to generate spin down wfs from spin up wfs
!!   indkk(:,1)=k point number of kptns1
!!   indkk(:,2)=symmetry operation to be applied to kpt1, to give kpt1a
!!    (if 0, means no symmetry operation, equivalent to identity )
!!   indkk(:,3:5)=shift in reciprocal space to be given to kpt1a,
!!    to give kpt1b, that is the closest to kpt2.
!!   indkk(:,6)=1 if time-reversal was used to generate kpt1a from kpt1, 0 otherwise
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  The tolerance tol12 aims at giving a machine-independent ordering.
!!  (this trick is used in bonds.f, listkk.f, prtrhomxmn.f and rsiaf9.f)
!!
!! PARENTS
!!      initberry,inwffil,mlwfovlp_qp,newsp
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine listkk(dksqmax,gmet,indkk,kptns1,kptns2,nkpt1,nkpt2,nsym,&
& sppoldbl,symafm,symmat,timrev,use_symrec)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'listkk'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt1,nkpt2,nsym,sppoldbl,timrev
 real(dp),intent(out) :: dksqmax
 logical,optional,intent(in) :: use_symrec
!arrays
 integer,intent(in) :: symafm(nsym),symmat(3,3,nsym)
 integer,intent(out) :: indkk(nkpt2*sppoldbl,6)
 real(dp),intent(in) :: gmet(3,3),kptns1(3,nkpt1),kptns2(3,nkpt2)

!Local variables-------------------------------
!scalars
 integer :: first_trial,ikpt1,ikpt2,isppol,isym,itimrev,jkpt1,jsym,jtime
 integer :: nsym_used,timrev_used,usesym
 real(dp) :: dksq,dksqmn
 character(len=500) :: message
!arrays
 integer :: dkint(3),jdkint(3)
 real(dp) :: dk(3),kpt1a(3)

! *************************************************************************

!DEBUG
!write(std_out,*)' listkk : enter'
!ENDDEBUG

 if(sppoldbl<1 .or. sppoldbl>2)then
   write(message, '(4a,i4,3a)' )ch10,&
&   ' listkk : BUG -',ch10,&
&   '  The value of sppoldbl is',sppoldbl,',',ch10,&
&   '  but it should be either 1 or 2.'
   call wrtout(std_out,message,'PERS')
   call leave_new('PERS')
 end if

!When usesym=0, the old way of converting the wavefunctions (without
!using the symmetries), is recovered.
 usesym=1

 nsym_used=nsym
 timrev_used=timrev
 if(usesym==0)nsym_used=1
 if(usesym==0)timrev_used=0

 dksqmax=zero
 do isppol=1,sppoldbl
   do ikpt2=1,nkpt2

     dksqmn=-one
     first_trial=1

     do itimrev=0,timrev_used
       do isym=1,nsym_used

!        Select magnetic characteristic of symmetries
         if(isppol==1 .and. symafm(isym)==-1)cycle
         if(isppol==2 .and. symafm(isym)==1)cycle

         do ikpt1=1,nkpt1

!          Compute symmetric point to kpt1
           if(usesym==1)then
!            original code only used transpose(symrel)
!            kpt1a(:)=symrel(1,:,isym)*kptns1(1,ikpt1)+&
!            &             symrel(2,:,isym)*kptns1(2,ikpt1)+&
!            &             symrel(3,:,isym)*kptns1(3,ikpt1)
             if (present(use_symrec)) then
               if (use_symrec) then
                 kpt1a(:) = MATMUL(symmat(:,:,isym),kptns1(:,ikpt1))
               else
                 kpt1a(:) = MATMUL(TRANSPOSE(symmat(:,:,isym)),kptns1(:,ikpt1))
               end if
             else
               kpt1a(:) = MATMUL(TRANSPOSE(symmat(:,:,isym)),kptns1(:,ikpt1))
             end if
             kpt1a(:)=(1-2*itimrev)*kpt1a(:)
           else
             kpt1a(:)=kptns1(:,ikpt1)
           end if

!          Compute difference with respect to kpt2, modulo a lattice vector
           dk(:)=kptns2(:,ikpt2)-kpt1a(:)
           if(usesym==1)then
!            The tolerance insure similar behaviour on different platforms
             dkint(:)=nint(dk(:)+tol12)
             dk(:)=dk(:)-dkint(:)
           else
             dkint(:)=0
           end if

!          Compute norm of the difference vector, and update kpt1 if better.
           dksq=gmet(1,1)*dk(1)**2+gmet(2,2)*dk(2)**2+&
&           gmet(3,3)*dk(3)**2+two*(gmet(2,1)*dk(2)*dk(1)+&
&           gmet(3,2)*dk(3)*dk(2)+gmet(3,1)*dk(3)*dk(1))
           if ( (ikpt1==1.and.first_trial==1).or.dksq+tol12<dksqmn) then
             dksqmn=dksq
             jkpt1=ikpt1
             jsym=isym
             jtime=itimrev
             jdkint(:)=dkint(:)
           end if

         end do ! ikpt1

         first_trial=0

       end do ! isym
     end do ! itimrev

     indkk(ikpt2+(isppol-1)*nkpt2,1)=jkpt1
     indkk(ikpt2+(isppol-1)*nkpt2,2)=jsym
     indkk(ikpt2+(isppol-1)*nkpt2,3:5)=jdkint(:)
     indkk(ikpt2+(isppol-1)*nkpt2,6)=jtime
     dksqmax=max(dksqmax,dksqmn)

     if(dksqmn<-tol12)then
       write(message, '(a,a,a,a,es16.6)' )ch10,&
&       ' listkk : BUG -',ch10,&
&       '  The minimum square of dk has negative norm: dksqmn=',dksqmn
       call wrtout(std_out,message,'PERS')
       call leave_new('PERS')
     end if

!    DEBUG
!    write(std_out,'(a,2i4,2x,6i3)' )' listkk: ikpt2,isppol,indkk(ikpt2+(isppol-1)*nkpt2,:)=',ikpt2,isppol,indkk(ikpt2+(isppol-1)*nkpt2,:)
!    ENDDEBUG

   end do ! ikpt2
 end do ! isppol

end subroutine listkk
!!***
