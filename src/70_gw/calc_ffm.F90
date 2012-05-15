!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_ffm
!! NAME
!! calc_ffm
!!
!! FUNCTION
!! Calculate nth frequency moment of imaginary part of DM or its inverse
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (Rhaltaf,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nomega: total number of frequencies
!! npw: DM size
!! nq: total number of q vectors
!!
!!
!! OUTPUT
!! see side effects
!!
!! SIDE EFFECTS 
!! text file contains the nth frequency moment for each G, G '', q
!!
!! NOTES
!! RShaltaf
!! nfmidm : positive calculate the frequency moment for the DM
!! nfmidm : negative calcualte the frequency moment for the inverted
!! nfmidm : 0 calcualte the frequency moment for the full polarizibility
!! note that in a well converged results fmidm=1 and fmidm=-1 differ only in
!! minus sign
!! see M. Taut, J. Phys. C: Solid State Phys. 18 (1985) 2677-2690.
!!
!! PARENTS
!!
!! CHILDREN
!!      cmod_qpg,simpson_int,xginv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_ffm(epsm1,nq,npw,nomega,omega,gprimd,qq,gvec,nfmidm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

 use m_abilasi,    only : xginv
 use m_vcoul,      only : cmod_qpg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_ffm'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfmidm,nomega,npw,nq
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: gprimd(3,3),qq(3,nq)
 complex(gwpc),intent(in) :: epsm1(npw,npw,nomega,nq)
 complex(dpc),intent(in) :: omega(nomega)

!Local variables ------------------------------
!scalars
 integer :: ig,igp,iomega,iq
 real(dp) :: omega_delta
 complex(gwpc) :: freqm
!arrays
 real(dp) :: b1(3),b2(3),b3(3)
 real(dp),allocatable :: fun(:),qplusg(:),res(:)
 complex(gwpc),allocatable :: eps(:,:,:)

!************************************************************************
 omega_delta=omega(2)-omega(1)

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)

 ABI_ALLOCATE(res,(nomega))
 ABI_ALLOCATE(qplusg,(npw))
 ABI_ALLOCATE(fun,(nomega))
 ABI_ALLOCATE(eps,(npw,npw,nomega))

 do iq=1,nq
  eps(:,:,:)=epsm1(:,:,:,iq)
! if fmidm is positive calcualte the inverted DM for each frequency
  do iomega=1,nomega
   if(nfmidm>=0)then
    call xginv(eps(:,:,iomega),npw)
   end if
  end do
! if nfmidm=0, calcualte the full polarizibility 
  if(nfmidm==0)then
   call cmod_qpg(nq,iq,qq,npw,gvec,gprimd,qplusg)
  end if
  do ig=1,npw
   if(nfmidm==0)then
    eps(ig,ig,:)=eps(ig,ig,:)-1
   end if
   do igp=1,npw
    if(nfmidm==0)then
     eps(ig,igp,:)=qplusg(ig)*qplusg(igp)*eps(ig,igp,:)/(4*pi)
    end if  
    fun(:)=real(omega(:)**abs(nfmidm))*aimag(eps(ig,igp,:))
    call simpson_int(nomega,omega_delta,fun,res)
    freqm=res(nomega)
    write(18,'(1x,i3,3x,i3,3x,i3,3x,f16.11,2x,f16.11)')ig,igp,iq,real(freqm),aimag(freqm)
   end do !igp
  end do  !ig
 end do ! iq
 rewind(18)

 ABI_DEALLOCATE(res)
 ABI_DEALLOCATE(qplusg)
 ABI_DEALLOCATE(fun)

end subroutine calc_ffm
!!***
