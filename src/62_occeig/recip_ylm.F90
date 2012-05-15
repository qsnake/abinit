!{\src2tex{textfont=tt}}
!!****f* ABINIT/recip_ylm
!! NAME
!! recip_ylm
!!
!! FUNCTION
!! Project input wavefunctions (real space) on to Ylm
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  bess_fit(mpw,nradintmax,ilang) = bessel functions for l=ilang splined
!!   with arguments $2 \pi |k+G| \Delta r$, for all G vectors in sphere
!!   and all points on radial grid.
!!  cgcband(2,npw_k)=wavefunction in recip space
!!  iatsph(natsph)=input variable iatsph, giving index of atoms around
!!    which ang mom projection has to be done
!!  istwfk= storage mode of cgcband
!!  mlang=maximum angular momentum
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms
!!  natsph=number of atoms around which ang mom projection has to be done (dimension of iatsph)
!!  npw_k=number of plane waves for kpt
!!  nradint(natsph)=number of points on radial real-space grid for a given atom iatsph
!!  nradintmax=dimension of rint array
!!  ntypat==number of types of atoms in cell.
!!  ph3d(2,npw_k,natom)=3-dim structure factors, for each atom and plane wave.
!!  prtsphere= if 1, print a complete analysis of the angular momenta in atomic spheres
!!  rint(nradintmax) = points on radial real-space grid for integration
!!  rmax(natom)=maximum radius for real space integration sphere
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  ylm(mpw,mlang*mlang)=real spherical harmonics for each G and k point
!!  znucl(ntypat)=gives the nuclear number for each type of atom
!!
!! OUTPUT
!!  sum_1ll_1atom(mlang,natsph)= projected scalars for each atom and ang. mom.
!!
!! NOTES
!!  ph3d atoms are ordered with atindx -> typat by typat
!!  The atoms have to be reordered !!
!!  ngfft is not used !!
!!
!! PARENTS
!!      partial_dos_fractions,wffile
!!
!! CHILDREN
!!      atmdata,dotprod_g,simpson_int
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine recip_ylm (bess_fit,cgcband,iatsph,istwfk,&
& nradint,nradintmax,mlang,mpi_enreg,mpw,natom,natsph,&
& npw_k,ntypat,ph3d,prtsphere,rint,rmax,sum_1ll_1atom,sum_1lm_1atom,typat,ucvol,ylm,znucl)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recip_ylm'
 use interfaces_32_util
 use interfaces_53_spacepar
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,mlang,mpw,natom,natsph,npw_k,nradintmax
 integer,intent(in) :: ntypat,prtsphere
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: iatsph(natsph),nradint(natsph)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: bess_fit(mpw,nradintmax,mlang),cgcband(2,npw_k)
 real(dp),intent(in) :: ph3d(2,npw_k,natom),rint(nradintmax)
 real(dp),intent(in) :: rmax(natom),ylm(mpw,mlang*mlang)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(out) :: sum_1ll_1atom(mlang,natsph)
 real(dp),intent(out) :: sum_1lm_1atom(mlang*mlang,natsph)

!Local variables-------------------------------
!scalars
 integer :: clmindex,iat,iatom,ilang,imylmind,ipw,ixint,ll,mm,option,reylmind
 real(dp) :: amu,doti,dotr,imil,invsqrt2,mmsign,monetom,rcov
 real(dp) :: reil,sum_all,tmp1,tmp2
 character(len=2) :: symbol
!arrays
 real(dp) :: sum_1atom(natsph),sum_1ll(mlang),sum_1lm(mlang*mlang)
 real(dp) :: tmppsia(2,npw_k),tmppsim(2,npw_k)
 real(dp),allocatable :: integ(:),psilmnorm(:),vect(:,:)

! *************************************************************************

 invsqrt2=one/sqrt(two)

!metric has been called in calling routine.

 sum_1lm_1atom(:,:) = zero
 sum_1ll_1atom(:,:) = zero

!Big loop  on all atoms
 do iat=1,natsph
   iatom = iatsph(iat)
   reil =  0 ; imil = -1
   ABI_ALLOCATE(integ,(nradint(iat)))
   ABI_ALLOCATE(psilmnorm,(nradint(iat)))

!  Temporary arrays for part of psi which depends only on iatom
   do ipw=1,npw_k
     tmppsia(1,ipw) = cgcband(1,ipw)*ph3d(1,ipw,iatom) &
&     -cgcband(2,ipw)*ph3d(2,ipw,iatom)
     tmppsia(2,ipw) = cgcband(1,ipw)*ph3d(2,ipw,iatom) &
&     +cgcband(2,ipw)*ph3d(1,ipw,iatom)
   end do

   do ilang=1,mlang
     ll=ilang-1
     tmp1 = -imil; tmp2 = reil
     reil =  tmp1; imil = tmp2
     monetom = (-one)**ll

     do mm=-ll,ll
       clmindex = (ll+1)**2-ll+mm
       reylmind = (ll+1)**2-ll+abs(mm)
       imylmind = (ll+1)**2-ll-abs(mm)
       mmsign = one
       if (mm < 0) then
         mmsign = -one
       end if

!      Temporary arrays for part of psi which doesnt depend on ixint
!      Take into account the fact that ylm are REAL spherical harmonics, see initylmg.f
       if(istwfk==1)then
         do ipw=1,npw_k
!          to get PDOS for real spherical harmonics, may be sufficient to multiply here by ylm instead of linear combination
!          tmppsim(1,ipw) = tmppsia(1,ipw)*ylm(ipw,clmindex)
!          tmppsim(2,ipw) = tmppsia(2,ipw)*ylm(ipw,clmindex)
           tmppsim(1,ipw) = invsqrt2* &
&           (        tmppsia(1,ipw)*ylm(ipw,reylmind)&
&           +mmsign*tmppsia(2,ipw)*ylm(ipw,imylmind) )
           tmppsim(2,ipw) = invsqrt2* &
&           (-mmsign*tmppsia(1,ipw)*ylm(ipw,imylmind)&
&           +tmppsia(2,ipw)*ylm(ipw,reylmind) )
         end do
       else
!        For time-reversal states, detailed treatment show that only the real or imaginary
!        part of tmppsia is needed here, depending on l .
!        TODO: check the invsqrt2 part! Could be incorrect if we go to real spherical harmonics
         if(mod(ll,2)/=1)then
           do ipw=1,npw_k
             tmppsim(1,ipw)= invsqrt2*       tmppsia(1,ipw)*ylm(ipw,reylmind)
             tmppsim(2,ipw)=-invsqrt2*mmsign*tmppsia(1,ipw)*ylm(ipw,imylmind)
           end do
         else
           do ipw=1,npw_k
             tmppsim(1,ipw)=invsqrt2*mmsign*tmppsia(2,ipw)*ylm(ipw,imylmind)
             tmppsim(2,ipw)=invsqrt2*       tmppsia(2,ipw)*ylm(ipw,reylmind)
           end do
         end if
       end if

       ABI_ALLOCATE(vect,(2,npw_k))
       do ixint=1,nradint(iat)
         vect(1,:)=bess_fit(1:npw_k,ixint,ilang)
         vect(2,:)=zero
         option=2
         call dotprod_g(dotr,doti,istwfk,mpi_enreg,npw_k,option,vect,tmppsim)

!        Multiply by 4 pi i^l and take norm
         psilmnorm(ixint) = (four_pi * rint(ixint))**2 * (dotr**2+doti**2)
       end do ! ixint
       ABI_DEALLOCATE(vect)

!      Integrate on rint: integrand is in psilmnorm
       call simpson_int(nradint(iat),rmax(iatom)/(nradint(iat)-1),psilmnorm,integ)
!      NOTE : could exploit full r dependency of integ
!      which is calculated in the call to simpson

       sum_1lm_1atom(clmindex,iat)=sum_1lm_1atom(clmindex,iat)+integ(nradint(iat))
       sum_1ll_1atom(ll+1,iat)    =sum_1ll_1atom(ll+1,iat)    +integ(nradint(iat))

       monetom = -monetom
     end do ! mm
   end do ! ilang
   ABI_DEALLOCATE(integ)
   ABI_DEALLOCATE(psilmnorm)
 end do ! iatom

!Normalize with unit cell volume
 sum_1lm_1atom(:,:) = sum_1lm_1atom(:,:) / ucvol
 sum_1ll_1atom(:,:) = sum_1ll_1atom(:,:) / ucvol

!Output
 if(prtsphere==1)then
   sum_1ll(:) = zero
   sum_1lm(:) = zero
   do iat=1,natsph
     sum_1atom(iat) = sum(sum_1lm_1atom(:,iat))
     sum_1ll(:)=sum_1ll(:)+sum_1ll_1atom(:,iat)
     sum_1lm(:)=sum_1lm(:)+sum_1lm_1atom(:,iat)
   end do
   sum_all = sum(sum_1atom)

   write(std_out,'(a)' ) ' Angular analysis '
   do iat=1,natsph
     call atmdata( amu,rcov,symbol,znucl(typat(iatsph(iat))) )
     write(std_out,'(a)' ) ' '
     write(std_out,'(a,i3,a,a,a,f10.6)' )&
&     ' Atom # ',iat, ' is  ',  symbol, &
&     ', in-sphere charge =',sum_1atom(iat)
     do ilang=1,mlang
       ll=ilang-1
       write(std_out,'(a,i1,a,f9.6,a,9f6.3)' )&
&       ' l=',ll,', charge=',sum_1ll_1atom(ll+1,iat),&
&       ', m=-l,l splitting:',sum_1lm_1atom(1+ll**2:(ll+1)**2,iat)
     end do ! ll
   end do ! iat
   write(std_out,'(a,a)') ch10,' Sum of angular contributions for all atomic spheres '
   do ll=0,mlang-1
     write(std_out,'(a,i1,a,f9.6,a,f9.6)' )&
&     ' l=',ll,', charge =',sum_1ll(ll+1),' proportion =',sum_1ll(ll+1)/sum_all
   end do
   write(std_out,'(a,a,f10.6)' ) ch10,' Total over all atoms and l=0 to 4 :',sum_all
   write(std_out,'(a)' ) ' '
 end if

end subroutine recip_ylm
!!***
