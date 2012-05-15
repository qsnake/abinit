!{\src2tex{textfont=tt}}
!!****f* ABINIT/occred
!! NAME
!! occred
!!
!! FUNCTION
!!  check occupations and redistribute them (?)
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2012 ABINIT group (NH)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  rdocc(nk,nb,nsppol) :: occupation numbers of natural orbitals
!!  nk :: number of k-points
!!  nb :: number of bands
!!  nsppol :: spin polarization (1 for spin-unpolarized, 2 for spin-polarized)
!!
!! OUTPUT
!!  rdocc(nk,nb,nsppol):: redistributed occupation numbers to fullfil N-representability
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      rdm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine occred(rdocc,nsppol,nk,nb)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'occred'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb,nk,nsppol
!arrays
 real(dp),intent(inout) :: rdocc(nk,nb,nsppol)

!Local variables-------------------------------
!scalars
 integer :: ib,ik,is
 real(dp) :: maxocc
 logical :: checkocc
 character(len=500) :: message
!arrays
 real(dp) :: excess(nk)

! *************************************************************************
 
!DEBUG
!write (std_out,*) ' occred : enter'
!ENDDEBUG

 if(nsppol==1) then
   maxocc=2.0
 else
   maxocc=1.0
 end if

 checkocc=.false.

 if(any(rdocc>maxocc)) then
   checkocc=.true.
   write(message,'(4a,F5.2)')ch10,&
&   ' checkocc: ',ch10,&
&   '  found occupation number larger than ', maxocc
   call wrtout(std_out,message,'COLL')
 end if

 if(any(rdocc<0d0)) then
   checkocc=.true.
   write(message,'(4a,F5.2)')ch10,&
&   ' checkocc: ',ch10,&
&   '  found occupation number smaller than ', 0d0
   call wrtout(std_out,message,'COLL')
 end if

 if(checkocc) then
   excess=0.0d0
   do is=1, nsppol
     do ib=1, nb
       do ik=1, nk
         if(rdocc(ik,ib,is)>maxocc) then
           excess(ik)=excess(ik)+(rdocc(ik,ib,is)-maxocc)
           rdocc(ik,ib,is)=maxocc
         end if
         if(rdocc(ik,ib,is)<0.0d0) then
           excess(ik)=excess(ik)+rdocc(ik,ib,is)
           rdocc(ik,ib,is)=0d0
         end if
       end do
     end do
     do ik=1, nk
       do while(excess(ik)>1d-11)
         do ib=1, nb
           if(rdocc(ik,ib,is)<maxocc) then
             rdocc(ik,ib,is)=rdocc(ik,ib,is)+excess(ik)
             excess(ik)=0d0
             if(rdocc(ik,ib,is)>maxocc) then      
               excess(ik)=rdocc(ik,ib,is)-maxocc
               rdocc(ik,ib,is)=maxocc             
             end if                                               
           end if                                                
         end do !loop over bands                                
       end do !excess loop
       do while(excess(ik)<-1d-11)                            
         do ib=1, nb                                        
           if(rdocc(ik,ib,is)>1d-11) then
             rdocc(ik,ib,is)=rdocc(ik,ib,is)+excess(ik)
             excess(ik)=0d0
             if(rdocc(ik,ib,is)<0d0) then
               excess(ik)=rdocc(ik,ib,is)
               rdocc(ik,ib,is)=0d0
             end if
           end if
         end do !loop over bands
       end do !excess loop
     end do !loop over kpoints
   end do !loop over spin
 end if !checkocc

 write(message,'(2a)')' rdm : checked occupation and redistributed ',ch10
 call wrtout(std_out,message,'COLL')



!DEBUG
!write (std_out,*) ' occred : exit'
!stop
!ENDDEBUG
 
end subroutine occred
!!***
