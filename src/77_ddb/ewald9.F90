!{\src2tex{textfont=tt}}
!!****f* ABINIT/ewald9
!!
!! NAME
!! ewald9
!!
!! FUNCTION
!! Compute ewald contribution to the dynamical matrix, at a given
!! q wavevector, including anisotropic dielectric tensor
!! and effective charges
!! See Phys. Rev. B 55, 10355 (1997), equations (71) to (75).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (DCA,JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell = lengths by which lattice vectors are multiplied
!! dielt(3,3)=dielectric tensor
!! gmet(3,3) = metric in reciprocal space.
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! natom=number of atoms in unit cell
!! qphon(3)=phonon wavevector (same system of coordinates as the
!!  reciprocal lattice vectors)
!! rmet = metric in real space
!! rprim(3,3)=dimensionless primitive translations in real space
!! sumg0: if=1, the sum in reciprocal space must include g=0,
!!  if=0, this contribution must be skipped (q=0 singularity)
!! ucvol=unit cell volume in (whatever length scale units)**3
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement
!!
!! OUTPUT
!! dyew(2,3,natom,3,natom)= Ewald part of the dynamical matrix,
!!  second energy derivative wrt xred(3,natom) in Hartrees
!!  (Denoted A-bar in the notes)
!!
!! NOTES
!! 1. The q=0 part should be subtracted, by another call to
!! the present routine, with q=0. The present routine correspond
!! to the quantity written A-bar in the explanatory notes.
!! If q=0 is asked, sumg0 should be put to 0. Otherwise, it should
!! be put to 1.
!! 2. Because this routine can be used many times in the
!! evaluation of phonons in ppddb9, it has been
!! optimized carefully. There is still possibility
!! for improvement, by using bloking on g s and r s !
!! 3. There can be small numerical variations due to the
!! fact that the input dielectric tensor is usually
!! not perfectly symmetric ....
!!
!! PARENTS
!!      gtdyn9,hybrid9,mkifc9
!!
!! CHILDREN
!!      derfc,leave_new,matr3inv,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ewald9(acell,dielt,dyew,gmet,gprim,natom,&
&                  qphon,rmet,rprim,sumg0,ucvol,xred,zeff)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ewald9'
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
 real(dp),intent(in) :: acell(3),dielt(3,3),gmet(3,3),gprim(3,3),qphon(3)
 real(dp),intent(in) :: rmet(3,3),rprim(3,3),xred(3,natom),zeff(3,3,natom)
 real(dp),intent(out) :: dyew(2,3,natom,3,natom)

!Local variables -------------------------
!scalars
 integer,parameter :: matom=1000,mr=10000
 integer :: i2,ia,ib,ig1,ig2,ig3,ii,ir,ir1,ir2,ir3,jj,mu,newg,newr,ng,nr,nu
 real(dp) :: arg1,arg2,arg3,arga,c123i,c123r,c23i,c23r,derfc_yy,detdlt
 real(dp) :: direct,eta,fac,fact1,fact2,fact3,gsq,invy,invy2,recip,reta,reta3
 real(dp) :: term1,term2,term3,term4,term5,y2,yy
! real(dp) :: tcpu,tcpui,twall,twalli 
 character(len=500) :: message
!arrays
 real(dp) :: c1i(2*mr+1),c1r(2*mr+1),c2i(2*mr+1),c2r(2*mr+1),c3i(2*mr+1)
 real(dp) :: c3r(2*mr+1),cosqxred(matom),gpq(3),gpqfac(3,3),gpqgpq(3,3)
 real(dp) :: invdlt(3,3),ircar(3),ircax(3),rr(3),sinqxred(matom)
 real(dp) :: xredcar(3,matom),xredcax(3,matom),xredicar(3),xredicax(3),xx(3)
 real(dp),allocatable :: dyewt(:,:,:,:,:)

! *********************************************************************

 if(matom<natom)then
   write(message, '(4a,i5,a,i5,a,a)' )ch10,&
&   ' ewald9 : ERROR -',ch10,&
&   '  matom=',matom,' is smaller than natom=',natom,ch10,&
&   '  Action : raise matom in ewald9.F90, and recompile the program.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 ABI_ALLOCATE(dyewt,(2,3,natom,3,natom))

!compute eta for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
 eta=pi*100.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)

 dyew(:,:,:,:,:)=0.0_dp
 dyewt(:,:,:,:,:)=0.0_dp

!Sum terms over g space:
 ng=0
 do
   ng=ng+1
   newg=0
   do ig3=-ng,ng
     do ig2=-ng,ng
       do ig1=-ng,ng
         if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng&
&         .or. ng==1 )then
!          
           gpq(1)=(ig1+qphon(1))*gprim(1,1)/acell(1)+(ig2+qphon(2))*&
&           gprim(1,2)/acell(2)+(ig3+qphon(3))*gprim(1,3)/acell(3)
           gpq(2)=(ig1+qphon(1))*gprim(2,1)/acell(1)+(ig2+qphon(2))*&
&           gprim(2,2)/acell(2)+(ig3+qphon(3))*gprim(2,3)/acell(3)
           gpq(3)=(ig1+qphon(1))*gprim(3,1)/acell(1)+(ig2+qphon(2))*&
&           gprim(3,2)/acell(2)+(ig3+qphon(3))*gprim(3,3)/acell(3)
           gsq=0.0_dp
           do jj=1,3
             do ii=1,3
               gpqgpq(ii,jj)=gpq(ii)*gpq(jj)
               gsq=gsq+gpqgpq(ii,jj)*dielt(ii,jj)
             end do
           end do

!          Skip q=0:
           if (gsq<1.0d-20) then
             if (sumg0==1) then
               write(message,'(a,a,a,a,a,a,a)' )&
&               ' ewald9 : ERROR -',ch10,&
&               '  The phonon wavelength should not be zero :',ch10,&
&               '  there are non-analytical terms that cannot be treated.',ch10,&
&               '  Action : subtract this wavelength from the input file.'
               call wrtout(std_out,message,'COLL')
               call leave_new('COLL')
             end if

           else

             arg1=(two_pi**2)*gsq/(4*eta)

!            Larger arg gives 0 contribution:
             if (arg1<=80._dp) then
               newg=1

!              Here calculate the term
               term1=exp(-arg1)/gsq
               do jj=1,3
                 do ii=1,3
                   gpqfac(ii,jj)=gpqgpq(ii,jj)*term1
                 end do
               end do

               do ia=1,natom
                 arga=two_pi*( (ig1+qphon(1))*xred(1,ia)&
&                 +(ig2+qphon(2))*xred(2,ia)&
&                 +(ig3+qphon(3))*xred(3,ia) )
                 cosqxred(ia)=cos(arga)
                 sinqxred(ia)=sin(arga)
               end do

!              First, the diagonal terms
               do nu=1,3
                 do ia=1,natom
                   do mu=nu,3
                     dyewt(1,mu,ia,nu,ia)=dyewt(1,mu,ia,nu,ia)+&
&                     gpqfac(mu,nu)
                   end do
                 end do
               end do

!              Then, the non-diagonal ones
               do ib=2,natom
                 do ia=1,ib-1
                   c123r=cosqxred(ia)*cosqxred(ib)+sinqxred(ia)*sinqxred(ib)
                   c123i=sinqxred(ia)*cosqxred(ib)-cosqxred(ia)*sinqxred(ib)
!                  The most inner loop
                   do nu=1,3
                     do mu=nu,3
                       dyewt(1,mu,ia,nu,ib)=dyewt(1,mu,ia,nu,ib)+&
&                       gpqfac(mu,nu)*c123r
                       dyewt(2,mu,ia,nu,ib)=dyewt(2,mu,ia,nu,ib)+&
&                       gpqfac(mu,nu)*c123i
                     end do
                   end do
                 end do
               end do

             end if
!            Endif g/=0 :
           end if
!          End triple summation over Gs:
         end if
       end do
     end do
   end do

!  Check if new shell must be calculated
   if(newg==0)exit

 end do

!Multiplies by common factor
 fact1=4.0_dp*pi/ucvol
 do ib=1,natom
   do ia=1,ib
     do nu=1,3
       do mu=nu,3
         dyewt(1,mu,ia,nu,ib)=dyewt(1,mu,ia,nu,ib)*fact1
         dyewt(2,mu,ia,nu,ib)=dyewt(2,mu,ia,nu,ib)*fact1
       end do
     end do
   end do
 end do

 reta=sqrt(eta)
 reta3=-eta*reta
 fac=4.0_dp/3.0_dp/sqrt(pi)
 fact2=2.0_dp/sqrt(pi)

!Calculating the inverse (transpose) of the dielectric tensor
 call matr3inv(dielt,invdlt)
!Calculating the determinant of the dielectric tensor
 detdlt=dielt(1,1)*dielt(2,2)*dielt(3,3)+dielt(1,3)*dielt(2,1)*&
& dielt(3,2)+dielt(1,2)*dielt(2,3)*dielt(3,1)-dielt(1,3)*&
& dielt(2,2)*dielt(3,1)-dielt(1,1)*dielt(2,3)*dielt(3,2)-&
& dielt(1,2)*dielt(2,1)*dielt(3,3)

 if(detdlt<tol6)then
   write(message, '(3a,es16.6,11a)' )&
&   ' ewald9 : ERROR -',ch10,&
&   '  The determinant of the dielectrix matrix, detdlt=',detdlt,' is smaller than 1.0d-6.',ch10,&
&   '  The use of the dipole-dipole model for interatomic force constants is not possible.',ch10,&
&   '  It is likely that you have not treated the electric field perturbations,',ch10,&
&   '  because you not are dealing with an insulator, so that',ch10,&
&   '  your dielectric matrix was simply set to zero in the Derivative DataBase.',ch10,&
&   '  Action : set the input variable dipdip to 0 .'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 fact3=reta3/sqrt(detdlt)

!Preparing the loop on real space
 do ia=1,natom
   do ii=1,3
     xredcar(ii,ia)=(xred(1,ia)*acell(1)*rprim(ii,1)+&
&     xred(2,ia)*acell(2)*rprim(ii,2)+&
&     xred(3,ia)*acell(3)*rprim(ii,3) )*reta
   end do
 end do
 do ia=1,natom
   do ii=1,3
     xredcax(ii,ia)= invdlt(1,ii)*xredcar(ii,ia)+&
&     invdlt(2,ii)*xredcar(ii,ia)+&
&     invdlt(3,ii)*xredcar(ii,ia)
   end do
 end do

!Prepare the evaluation of exp(iq*R)
 do ir=-mr,mr
   arg1=-two_pi*qphon(1)*ir
   arg2=-two_pi*qphon(2)*ir
   arg3=-two_pi*qphon(3)*ir
   c1r(ir+mr+1)=cos(arg1)
   c1i(ir+mr+1)=sin(arg1)
   c2r(ir+mr+1)=cos(arg2)
   c2i(ir+mr+1)=sin(arg2)
   c3r(ir+mr+1)=cos(arg3)
   c3i(ir+mr+1)=sin(arg3)
 end do

 do nr=1,mr
   newr=0

!  Begin big loop on real space vectors
   do ir3=-nr,nr
     do ir2=-nr,nr

!      Here, construct the cosine and sine of q*R for
!      components 2 and 3
       c23r = c2r(ir2+mr+1) * c3r(ir3+mr+1)&
&       - c2i(ir2+mr+1) * c3i(ir3+mr+1)
       c23i = c2i(ir2+mr+1) * c3r(ir3+mr+1)&
&       + c2r(ir2+mr+1) * c3i(ir3+mr+1)

!      Also multiplies by fact3, because it is a rather
!      economical place to do so
       c23r=c23r * fact3
       c23i=c23i * fact3

       do ir1=-nr,nr
         if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr&
&         .or. nr==1 )then

!          This is the real part and imaginary
!          part of the phase factor exp(iq*R)
           c123r = c1r(ir1+mr+1) * c23r - c1i(ir1+mr+1) * c23i
           c123i = c1i(ir1+mr+1) * c23r + c1r(ir1+mr+1) * c23i

           do ii=1,3
             ircar(ii)= ( ir1*acell(1)*rprim(ii,1)+&
&             ir2*acell(2)*rprim(ii,2)+&
&             ir3*acell(3)*rprim(ii,3) ) * reta
           end do
           do ii=1,3
             ircax(ii)=   invdlt(1,ii)*ircar(ii)+&
&             invdlt(2,ii)*ircar(ii)+&
&             invdlt(3,ii)*ircar(ii)
           end do

!          Here loops on atoms
           do ib=1,natom
             do ii=1,3
               xredicar(ii)=ircar(ii)-xredcar(ii,ib)
               xredicax(ii)=ircax(ii)-xredcax(ii,ib)
             end do
             do ia=1,ib
               do ii=1,3
                 rr(ii)=xredicar(ii)+xredcar(ii,ia)
                 xx(ii)=xredicax(ii)+xredcax(ii,ia)
               end do

               y2=rr(1)*xx(1)+rr(2)*xx(2)+rr(3)*xx(3)

!              The atoms should not be too far of each other
               if (y2<64.0_dp) then
!                Note: erfc(8) is about 1.1e-29,
!                so dont bother with larger y.
!                Also: exp(-64) is about 1.6e-28,
!                do dont bother with larger y**2 in exp.

!                Avoid zero denominators in term:
                 if (y2>=1.0d-24) then
                   newr=1
                   yy=sqrt(y2)
                   invy=1.0_dp/yy
                   invy2=invy**2
                   call derfc(derfc_yy,yy)
                   term2=derfc_yy*invy*invy2
                   term3=fact2*exp(-y2)*invy2
                   term4=-(term2+term3)
                   term5=(3*term2+term3*(3.0_dp+2.0_dp*y2))*invy2
                   do nu=1,3
                     do mu=nu,3
                       dyewt(1,mu,ia,nu,ib)=dyewt(1,mu,ia,nu,ib)+&
&                       c123r*(xx(nu)*xx(mu)*term5+term4*invdlt(nu,mu))
                       dyewt(2,mu,ia,nu,ib)=dyewt(2,mu,ia,nu,ib)+&
&                       c123i*(xx(nu)*xx(mu)*term5+term4*invdlt(nu,mu))
                     end do
                   end do
                 else
!                  If zero denominator, the atoms should be identical
                   if (ia/=ib)then
                     write(message, '(a,a,a,a,a,a,a,i5,a,i5,a)' )&
&                     ' ewald9 : ERROR -',ch10,&
&                     '  The distance between two atoms seem to vanish.',ch10,&
&                     '  This is not allowed.',ch10,&
&                     '  Action : check the input for the atoms number',ia,' and',ib,'.'
                     call wrtout(std_out,message,'COLL')
                     call leave_new('COLL')
                   else
!                    This is the correction when the atoms are identical
                     do nu=1,3
                       do mu=1,3
                         dyewt(1,mu,ia,nu,ib)=dyewt(1,mu,ia,nu,ib)+&
&                         fac*reta3*invdlt(nu,mu)/sqrt(detdlt)
                       end do
                     end do
                   end if

!                  End the condition for avoiding zero denominators
                 end if

!                End the condition of too large distance between atoms
               end if

!              End loop over ia and ib :
             end do
           end do
!          End triple loop over real space points:
         end if
       end do
     end do
   end do

!  Check if new shell must be calculated
   if(newr==0)exit
   if(newr==1 .and. nr==mr)then
     write(message,'(a,a,a)' )&
&     ' ewald9 : BUG ',ch10,&
&     '  mr is too small '
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

 end do

!DEBUG
!write(std_out,*)' ewald9 : quit with nr=',nr
!ENDDEBUG

!Now, symmetrizes
 do ib=1,natom-1
   do nu=1,3
     do ia=ib+1,natom
       do mu=nu,3
         dyewt(1,mu,ia,nu,ib)=dyewt(1,mu,ib,nu,ia)
         dyewt(2,mu,ia,nu,ib)=-dyewt(2,mu,ib,nu,ia)
       end do
     end do
   end do
 end do

 do ib=1,natom
   do nu=2,3
     do ia=1,natom
       do mu=1,nu-1
         dyewt(1,mu,ia,nu,ib)=dyewt(1,nu,ia,mu,ib)
         dyewt(2,mu,ia,nu,ib)=dyewt(2,nu,ia,mu,ib)
       end do
     end do
   end do
 end do

!Tests
!write(std_out,*)' ewald9 : take into account the effective charges '
!call timein(tcpu,twall)
!write(std_out,1000) tcpu-tcpui,twall-twalli
!
 do ib=1,natom
   do nu=1,3
     do ia=1,natom
       do mu=1,3
         do i2=1,2
           dyew(i2,mu,ia,nu,ib)=0.0_dp
           do jj=1,3
             do ii=1,3
               dyew(i2,mu,ia,nu,ib)=dyew(i2,mu,ia,nu,ib)+zeff(ii,mu,ia)*&
&               zeff(jj,nu,ib)*dyewt(i2,ii,ia,jj,ib)
             end do
           end do
         end do
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(dyewt)

end subroutine ewald9
!!***
