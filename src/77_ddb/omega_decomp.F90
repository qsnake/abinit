!{\src2tex{textfont=tt}}
!!****f* ABINIT/omega_decomp
!! NAME
!!  omega_decomp
!!
!! FUNCTION
!! Compute and return the eigenvalues (frequencies) of the short-range and 
!! long-range part of the dynamical matrix (as done in mkifc9).
!! See Europhys. Lett. 33 p.713 (1996) for details.
!! (included by U. Aschauer and EB)
!!  
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2012 ABINIT group (EB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!! mkifc9 has been modified accordingly to store and call omega_decomp 
!! (marked by "! OmegaSRLR:" in mkifc9.F90).
!!
!! PARENTS
!!      mkifc9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine omega_decomp(amu,natom,ntypat,typat,&
&   dynmatfl,dynmatsr,dynmatlr,iqpt,nqpt,eigenvec)

 use m_profiling
    
 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'omega_decomp'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 integer,intent(in) :: iqpt,nqpt
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat)
 real(dp),intent(inout) ::dynmatfl(2,3,natom,3,natom,nqpt)
 real(dp),intent(inout) ::dynmatsr(2,3,natom,3,natom,nqpt)
 real(dp),intent(inout) ::dynmatlr(2,3,natom,3,natom,nqpt)
 real(dp),intent(in)    ::eigenvec(2*3*natom*3*natom)

!Local variables -------------------------
!scalars
 integer :: i1,i2,idir1,idir2,imode,ipert1,ipert2,index1,index2
 real(dp),parameter :: break_symm=1.0d-12
 real(dp) :: fac
!arrays
 real(dp) :: nearidentity(3,3)
 real(dp) :: omegafl, omegasr, omegalr
 real(dp) :: sumfl,sumlr,sumsr,asr
! *********************************************************************

!write(ab_out,*)''
!write(std_out,*) 'SR/LR decomposition: enter for wavevector number :',iqpt
 
!apply asr (note the lr part is already asred by construction in mkifc9)
 do ipert1=1,natom
   do idir1=1,3
     do idir2=1,3
       asr=0.0d0
       do ipert2=1,natom
         asr=asr+dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)
       end do
       dynmatfl(1,idir1,ipert1,idir2,ipert1,iqpt)=dynmatfl(1,idir1,ipert1,idir2,ipert1,iqpt)-asr
       dynmatsr(1,idir1,ipert1,idir2,ipert1,iqpt)=dynmatsr(1,idir1,ipert1,idir2,ipert1,iqpt)-asr
     end do
   end do
 end do


!This slight breaking of the symmetry allows the
!results to be more portable between machines
 nearidentity(:,:)=1.0
 nearidentity(1,1)=1.0+break_symm
 nearidentity(3,3)=1.0-break_symm


!Include Mass
 do ipert1=1,natom
   do ipert2=1,natom
     
     fac=1.0d0/sqrt(amu(typat(ipert1))*amu(typat(ipert2)))/amu_emass
     
     do idir1=1,3
       do idir2=1,3
         
         dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)
         
         dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)

         dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&         dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)*&
&         fac*nearidentity(idir1,idir2)

!        This is to break slightly the translation invariance, and make
!        the automatic tests more portable

         if(ipert1==ipert2 .and. idir1==idir2)then
           
           dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatfl(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0

           dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatsr(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0

           dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)=&
&           dynmatlr(1,idir1,ipert1,idir2,ipert2,iqpt)+&
&           break_symm*natom/amu_emass/idir1*0.01d0

         end if
       end do
     end do
   end do
 end do


!Calculation of <eigvec|Dyn_tot,Dyn_SR,Dyn_LR|eigenvec>=omega**2

!write(ab_out,*)''
!write(ab_out,*)'==============================================================================='
 write(ab_out,*)''
 write(ab_out,*) 'Long-Range/Short-Range decomposed phonon freq. (cm-1)**2'
 write(ab_out,*) 'at wavevector number:',iqpt
 write(ab_out,*)''
 write(ab_out,'(a13,1x,a16,2x,a16,2x,a16)') ' Mode number.','tot**2','SR**2','LR**2'
 write(std_out,'(a13,1x,a16,2x,a16,2x,a16)') ' Mode number.','tot**2','SR**2','LR**2'
!write(ab_out,'(a12,2x,a10,2x,a10,2x,a10,2x,a16,2x,a16,2x,a16)') 'Mode number.','tot','SR','LR','tot**2','SR**2','LR**2'
!write(std_out,'(a12,2x,a10,2x,a10,2x,a10,2x,a16,2x,a16,2x,a16)') 'Mode number.','tot','SR','LR','tot**2','SR**2','LR**2'

 do imode=1,3*natom
   
   sumfl=0.0d0
   sumlr=0.0d0
   sumsr=0.0d0
   
   do ipert1=1,natom
     do ipert2=1,natom
       do i1=1,3
         do i2=1,3
           
           index1=i1+(ipert1-1)*3+3*natom*(imode-1)
           index2=i2+(ipert2-1)*3+3*natom*(imode-1)
           
           sumfl = sumfl + eigenvec(2*index1-1) * dynmatfl(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
           sumlr = sumlr + eigenvec(2*index1-1) * dynmatlr(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
           sumsr = sumsr + eigenvec(2*index1-1) * dynmatsr(1,i1,ipert1,i2,ipert2,iqpt) * eigenvec(2*index2-1)
         end do
       end do
     end do
   end do
   
   sumfl = sumfl * Ha_cmm1 * Ha_cmm1
   sumsr = sumsr * Ha_cmm1 * Ha_cmm1
   sumlr = sumlr * Ha_cmm1 * Ha_cmm1
   
!  Compute omega=sqrt(omega**2)
   if(sumfl>=1.0d-16)then
     omegafl=sqrt(sumfl)
   else if(sumfl>=-1.0d-16)then
     omegafl=zero
   else
     omegafl=-sqrt(-sumfl)
   end if

   if(sumsr>=1.0d-16)then
     omegasr=sqrt(sumsr)
   else if(sumsr>=-1.0d-16)then
     omegasr=zero
   else
     omegasr=-sqrt(-sumsr)
   end if

   if(sumlr>=1.0d-16)then
     omegalr=sqrt(sumlr)
   else if(sumlr>=-1.0d-16)then
     omegalr=zero
   else
     omegalr=-sqrt(-sumlr)
   end if

!  Output
   write(ab_out,'(i4,10x,sf16.4,2x,f16.4,2x,f16.4)') imode,sumfl,sumsr,sumlr
   write(std_out,'(i4,10x,sf16.4,2x,f16.4,2x,f16.4)') imode,sumfl,sumsr,sumlr
!  write(ab_out,'(i4,8x,f10.4,2x,f10.4,2x,f10.4,2x,sf16.6,2x,f16.6,2x,f16.6)') imode,omegafl,omegasr,omegalr,sumfl,sumsr,sumlr
!  write(std_out,'(i4,8x,f10.4,2x,f10.4,2x,f10.4,2x,sf16.6,2x,f16.6,2x,f16.6)') imode,omegafl,omegasr,omegalr,sumfl,sumsr,sumlr


 end do

!DEBUG
!write (std_out,*) ' omega_decomp : enter'
!ENDDEBUG

!DEBUG                                           ! to be uncommented, if needed
!if(option/=1 .and. option/=2 )then
!write(message,'(a,a,a,a,a,a,i6)') ch10,&
!&  ' omega_decomp: BUG -',ch10,&
!&  '  The argument option should be 1 or 2,',ch10,&
!&  '  however, option=',option
!call wrtout(std_out,message,'COLL')
!call leave_new('COLL')
!endif
!if(sizein<1)then
!write(message,'(a,a,a,a,a,a,i6)') ch10,&
!&  ' omega_decomp: BUG -',ch10,&
!&  '  The argument sizein should be a positive number,',ch10,&
!&  '  however, sizein=',sizein
!call wrtout(std_out,message,'COLL')
!call leave_new('COLL')
!endif
!ENDDEBUG


!DEBUG
!write (std_out,*) ' omega_decomp : exit'
!stop
!ENDDEBUG

end subroutine omega_decomp
!!***
