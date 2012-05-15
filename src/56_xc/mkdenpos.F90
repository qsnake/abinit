!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkdenpos
!! NAME
!! mkdenpos
!!
!! FUNCTION
!! Make a ground-state density positive everywhere :
!! when the density (or spin-density) is smaller than xc_denpos,
!! set it to the value of xc_denpos
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components (max. 2)
!!  option=0 if density rhonow is stored as (up,dn)
!!         1 if density rhonow is stored as (up+dn,up)
!!         Active only when nspden=2
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  Input/output
!!  iwarn=At input: iwarn=0 a warning will be printed when rho is negative
!!                  iwarn>0 no warning will be printed out
!!        At output: iwarn is increased by 1
!!  rhonow(nfft,nspden)=electron (spin)-density in real space,
!!     either on the unshifted grid (if ishift==0,
!!     then equal to rhor),or on the shifted grid
!!
!! NOTES
!!  At this stage, rhonow(:,1:nspden) contains the density in real space,
!!  on the unshifted or shifted grid. Now test for negative densities
!!  Note that, ignoring model core charge, as long as boxcut>=2
!!  the shifted density is derivable from the square of a Fourier
!!  interpolated charge density => CANNOT go < 0.
!!  However, actually can go < 0 to within machine precision;
!!  do not print useless warnings in this case, just fix it.
!!  Fourier interpolated core charge can go < 0 due to Gibbs
!!  oscillations; could avoid this by recomputing the model core
!!  charge at the new real space grid points (future work).
!!
!! PARENTS
!!      pawxc,pawxcm,pawxcmpositron,pawxcpositron,poslifetime,psolver_rhohxc
!!      rhohxc,rhohxcpositron
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkdenpos(iwarn,nfft,nspden,option,rhonow,xc_denpos)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkdenpos'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden,option
 integer,intent(inout) :: iwarn
 real(dp), intent(in) :: xc_denpos
!arrays
 real(dp),intent(inout) :: rhonow(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft,ispden,numneg
 real(dp) :: rhotmp,worst
 character(len=500) :: message
!arrays
 real(dp) :: rho(2)

! *************************************************************************

!DEBUG
 write(std_out,*)' mkdenpos : enter '
 write(std_out,*)' xc_denpos=',xc_denpos
!ENDDEBUG

 numneg=0
 worst=zero

 if(nspden==1)then

!  Non spin-polarized
!  $OMP PARALLEL DO PRIVATE(ifft,rhotmp) &
!  $OMP&REDUCTION(MIN:worst) &
!  $OMP&REDUCTION(+:numneg) &
!  $OMP&SHARED(nfft,rhonow)
   do ifft=1,nfft
     rhotmp=rhonow(ifft,1)
     if(rhotmp<xc_denpos)then
       if(rhotmp<-xc_denpos)then
!        This case is probably beyond machine precision considerations
         worst=min(worst,rhotmp)
         numneg=numneg+1
       end if
       rhonow(ifft,1)=xc_denpos
     end if
   end do
!  $OMP END PARALLEL DO
 else if (nspden==2) then

!  Spin-polarized

!  rhonow is stored as (up,dn)
   if (option==0) then

!    $OMP PARALLEL DO PRIVATE(ifft,ispden,rho,rhotmp) &
!    $OMP&REDUCTION(MIN:worst) &
!    $OMP&REDUCTION(+:numneg) &
!    $OMP&SHARED(nfft,nspden,rhonow)
     do ifft=1,nfft
!      For polarized case, rho(1) is spin-up density, rho(2) is spin-down density
       rho(1)=rhonow(ifft,1)
       rho(2)=rhonow(ifft,2)
       do ispden=1,nspden
         if (rho(ispden)<xc_denpos) then
           if (rho(ispden)<-xc_denpos) then
!            This case is probably beyond machine precision considerations
             worst=min(worst,rho(ispden))
             numneg=numneg+1
           end if
           rhonow(ifft,ispden)=xc_denpos
         end if
       end do
     end do
!    $OMP END PARALLEL DO

!    rhonow is stored as (up+dn,up)
   else if (option==1) then

!    $OMP PARALLEL DO PRIVATE(ifft,ispden,rho,rhotmp) &
!    $OMP&REDUCTION(MIN:worst) &
!    $OMP&REDUCTION(+:numneg) &
!    $OMP&SHARED(nfft,nspden,rhonow)
     do ifft=1,nfft
!      For polarized case, rho(1) is spin-up density, rho(2) is spin-down density
       rho(1)=rhonow(ifft,2)
       rho(2)=rhonow(ifft,1)-rho(1)
       do ispden=1,nspden
         if (rho(ispden)<xc_denpos) then
           if (rho(ispden)<-xc_denpos) then
!            This case is probably beyond machine precision considerations
             worst=min(worst,rho(ispden))
             numneg=numneg+1
           end if
           rho(ispden)=xc_denpos
           rhonow(ifft,1)=rho(1)+rho(2)
           rhonow(ifft,2)=rho(1)
         end if
       end do
     end do
!    $OMP END PARALLEL DO

   end if  ! option

 else

   write(message, '(4a)' ) ch10,&
&   ' mkdenpos : BUG -',ch10,&
&   '  nspden>2 not allowed !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')

!  End choice between non-spin polarized and spin-polarized.
 end if

 if (numneg>0) then
   if (iwarn==0) then
     write(message, '(a,a,a,a,i10,a,a,a,es10.2,a,e10.2,a,a,a,a)' ) ch10,&
&     ' mkdenpos : WARNING -',ch10,&
&     '  Density went too small (lower than xc_denpos) at',numneg,' points',ch10,&
&     '  and was set to xc_denpos=',xc_denpos,'.  Lowest was ',worst,'.',ch10,&
&     '  Likely due to too low boxcut or too low ecut for',&
&     ' pseudopotential core charge.'
     call wrtout(std_out,message,'COLL')
   end if
   iwarn=iwarn+1
 end if

!DEBUG
!write(std_out,'(a)') ' mkdenpos :  '
!write(std_out,'(a)')&
!& '   ir              rhonow(:,:,1:4)'
!do ir=1,nfft
!write(message,'(i5,a,4es14.6)')ir,' ',&
!&  rhonow(ir,1,1:4)
!call wrtout(std_out,message,'COLL')
!if(nspden==2)then
!write(message,'(a,2es14.6)')'               ',rhonow(ir,2,1:4)
!call wrtout(std_out,message,'COLL')
!end if
!end do
!ENDDEBUG

end subroutine mkdenpos
!!***
