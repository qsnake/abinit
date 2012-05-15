!{\src2tex{textfont=tt}}
!!****f* ABINIT/ewald3
!!
!! NAME
!! ewald3
!!
!! FUNCTION
!! Compute ewald contribution to the dynamical matrix, at a given
!! q wavevector.
!! Note : the q=0 part should be subtracted, by another call to
!! the present routine, with q=0. The present routine correspond
!! to the quantity C_bar defined in Eq.(24) or (27) in Phys. Rev. B 55, 10355 (1997). 
!! The two calls correspond to Eq.(23) of the same paper.
!! If q=0 is asked, sumg0 should be put to 0. Otherwise, it should
!! be put to 1.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (length units **-2)
!! natom=number of atoms in unit cell
!! qphon(3)=phonon wavevector (same system of coordinates as the
!!          reciprocal lattice vectors)
!! rmet(3,3)=metric tensor in real space (length units squared)
!! sumg0: if=1, the sum in reciprocal space must include g=0,
!!   if=0, this contribution must be skipped (q=0 singularity)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume in (whatever length scale units)**3
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! dyew(2,3,natom,3,natom)= Ewald part of the dynamical matrix,
!!    second energy derivative wrt xred(3,natom), Hartrees.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      derfc,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ewald3(dyew,gmet,natom,qphon,rmet,sumg0,typat,ucvol,xred,zion)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ewald3'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,sumg0
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gmet(3,3),qphon(3),rmet(3,3),xred(3,natom),zion(*)
 real(dp),intent(out) :: dyew(2,3,natom,3,natom)

!Local variables -------------------------
!nr, ng affect convergence of sums (nr=3,ng=5 is not good enough):
!scalars
 integer,parameter :: im=2,ng=10,nr=6,re=1
 integer :: ia,ib,ig1,ig2,ig3,ii,ir1,ir2,ir3,mu,nu
 real(dp) :: arg,arga,argb,c1i,c1r,da1,da2,da3,derfc_arg
 real(dp) :: direct,dot1,dot2,dot3,dotr1,dotr2,dotr3
 real(dp) :: eta,fac,gdot12,gdot13,gdot23,gsq,gsum,norm1
 real(dp) :: r1,r2,r3,rdot12,rdot13,rdot23,recip,reta
 real(dp) :: reta3m,rmagn,rsq,term,term1,term2
 real(dp) :: term3
 character(len=500) :: message
!arrays
 real(dp) :: gpq(3),rq(3)

! *********************************************************************

!Compute eta for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
 eta=pi*(dble(ng)/dble(nr))*sqrt(1.69_dp*recip/direct)

!Test Ewald s summation
!eta=1.2_dp*eta

!Sum terms over g space:
 fac=pi**2/eta
 gsum=0._dp
 da1=0._dp
 da2=0._dp
 da3=0._dp
 dyew(:,:,:,:,:)=0.0_dp
 do ig3=-ng,ng
   do ig2=-ng,ng
     do ig1=-ng,ng
       gpq(1)=dble(ig1)+qphon(1)
       gpq(2)=dble(ig2)+qphon(2)
       gpq(3)=dble(ig3)+qphon(3)
       gdot12=gmet(2,1)*gpq(1)*gpq(2)
       gdot13=gmet(3,1)*gpq(1)*gpq(3)
       gdot23=gmet(3,2)*gpq(2)*gpq(3)
       dot1=gmet(1,1)*gpq(1)**2+gdot12+gdot13
       dot2=gmet(2,2)*gpq(2)**2+gdot12+gdot23
       dot3=gmet(3,3)*gpq(3)**2+gdot13+gdot23
       gsq=dot1+dot2+dot3
!      Skip q=0:
       if (gsq<1.0d-20) then
         if (sumg0==1) then
           write(message, '(a,a,a,a,a,a,a)' )&
&           ' ewald3 : ERROR -',ch10,&
&           '  The phonon wavelength should not be zero : ',ch10,&
&           '  there are non-analytical terms that the code cannot handle.',&
&           ch10,'  Action : subtract this wavelength from the input.'
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if
       else
         arg=fac*gsq
!        Larger arg gives 0 contribution:
         if (arg <= 80._dp) then
           term=exp(-arg)/gsq
           do ia=1,natom
             arga=two_pi*(gpq(1)*xred(1,ia)+gpq(2)*xred(2,ia)+gpq(3)*xred(3,ia))
             do ib=1,ia
               argb=two_pi*(gpq(1)*xred(1,ib)+gpq(2)*xred(2,ib)+gpq(3)*xred(3,ib))
               arg=arga-argb
               c1r=cos(arg)*term
               c1i=sin(arg)*term

               do mu=1,3
                 do nu=1,mu
                   dyew(re,mu,ia,nu,ib)=dyew(re,mu,ia,nu,ib)+gpq(mu)*gpq(nu)*c1r
                   dyew(im,mu,ia,nu,ib)=dyew(im,mu,ia,nu,ib)+gpq(mu)*gpq(nu)*c1i
                 end do
               end do

             end do
           end do
         end if
!        Endif g/=0 :
       end if
!      End triple loop over G s:
     end do
   end do
 end do

!End G summation by accounting for some common factors.
!(for the charges:see end of routine)
 norm1=4.0_dp*pi/ucvol
 do ia=1,natom
   do ib=1,ia
     do mu=1,3
       do nu=1,mu
         dyew(:,mu,ia,nu,ib)=dyew(:,mu,ia,nu,ib)*norm1
       end do
     end do
   end do
 end do


!Do sums over real space:
 reta=sqrt(eta)
 reta3m=-eta*reta
 fac=4._dp/3.0_dp/sqrt(pi)
 do ir3=-nr,nr
   do ir2=-nr,nr
     do ir1=-nr,nr
       arg=-two_pi*(qphon(1)*ir1+qphon(2)*ir2+qphon(3)*ir3)
       c1r=cos(arg)*reta3m
       c1i=sin(arg)*reta3m
       do ia=1,natom
         do ib=1,ia
           r1=dble(ir1)+xred(1,ia)-xred(1,ib)
           r2=dble(ir2)+xred(2,ia)-xred(2,ib)
           r3=dble(ir3)+xred(3,ia)-xred(3,ib)
           rdot12=rmet(2,1)*r1*r2
           rdot13=rmet(3,1)*r1*r3
           rdot23=rmet(3,2)*r2*r3
           dotr1=rmet(1,1)*r1**2+rdot12+rdot13
           dotr2=rmet(2,2)*r2**2+rdot12+rdot23
           dotr3=rmet(3,3)*r3**2+rdot13+rdot23
           rsq=dotr1+dotr2+dotr3
           rmagn=sqrt(rsq)
!          Avoid zero denominators in term :
           if (rmagn>=1.0d-12) then
             arg=reta*rmagn
             term=0._dp
             if (arg<8.0_dp) then
!              Note: erfc(8) is about 1.1e-29,
!              so don t bother with larger arg.
!              Also: exp(-64) is about 1.6e-28,
!              so don t bother with larger arg**2 in exp.
               call derfc(derfc_arg,arg)
               term=derfc_arg/arg**3
               term1=2.0_dp/sqrt(pi)*exp(-arg**2)/arg**2
               term2=-(term+term1)
               term3=(3*term+term1*(3.0_dp+2.0_dp*arg**2))/rsq
               rq(1)=rmet(1,1)*r1+rmet(1,2)*r2+rmet(1,3)*r3
               rq(2)=rmet(2,1)*r1+rmet(2,2)*r2+rmet(2,3)*r3
               rq(3)=rmet(3,1)*r1+rmet(3,2)*r2+rmet(3,3)*r3
               do mu=1,3
                 do nu=1,mu
                   dyew(re,mu,ia,nu,ib)=dyew(re,mu,ia,nu,ib)+&
&                   c1r*(rq(mu)*rq(nu)*term3+rmet(mu,nu)*term2)
                   dyew(im,mu,ia,nu,ib)=dyew(im,mu,ia,nu,ib)+&
&                   c1i*(rq(mu)*rq(nu)*term3+rmet(mu,nu)*term2)
                 end do
               end do
             end if
           else
             if (ia/=ib)then
               write(message, '(a,a,a,a,a,a,a,i5,a,i5,a)' )&
&               ' ewald3 : ERROR -',ch10,&
&               '  The distance between two atoms vanishes.',ch10,&
&               '  This is not allowed.',ch10,&
&               '  Action : check the input for the atoms number',ia,' and',ib,'.'
               call wrtout(std_out,message,'COLL')
               call leave_new('COLL')
             else
               do mu=1,3
                 do nu=1,mu
                   dyew(re,mu,ia,nu,ib)=dyew(re,mu,ia,nu,ib)+&
&                   fac*reta3m*rmet(mu,nu)
                 end do
               end do
             end if
           end if
!          End loop over ib:
         end do
!        End loop over ia:
       end do
!      End triple loop over real space points:
     end do
   end do
 end do

!Take account of the charges
!write(std_out,*)' '
 do ia=1,natom
   do ib=1,ia
     do mu=1,3
       do nu=1,mu
         do ii=1,2
!          write(std_out,*)dyew(ii,mu,ia,nu,ib)
           dyew(ii,mu,ia,nu,ib)=dyew(ii,mu,ia,nu,ib)*&
&           zion(typat(ia))*zion(typat(ib))
         end do
       end do
     end do
   end do
 end do

!Symmetrize with respect to the directions
 do ia=1,natom
   do ib=1,ia
     do mu=1,3
       do nu=1,mu
         dyew(re,nu,ia,mu,ib)=dyew(re,mu,ia,nu,ib)
         dyew(im,nu,ia,mu,ib)=dyew(im,mu,ia,nu,ib)
       end do
     end do
   end do
 end do

!Fill the upper part of the matrix, with the hermitian conjugate
 do ia=1,natom
   do ib=1,ia
     do nu=1,3
       do mu=1,3
         dyew(re,mu,ib,nu,ia)=dyew(re,mu,ia,nu,ib)
         dyew(im,mu,ib,nu,ia)=-dyew(im,mu,ia,nu,ib)
       end do
     end do
   end do
 end do

end subroutine ewald3
!!***
