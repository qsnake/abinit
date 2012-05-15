!{\src2tex{textfont=tt}}
!!****f* ABINIT/vloca3
!! NAME
!! vloca3
!!
!! FUNCTION
!! Compute local part of 1st-order potential from the appropriate
!! atomic pseudopotential with structure and derivative factor.
!! In case of derivative with respect to k or
!! electric field perturbation, the 1st-order local potential
!! vanishes.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  cplex: if 1, real space 1-order functions on FFT grid
!!    are REAL, if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric (Bohr**-2)
!!  gsqcut=cutoff G**2 for included G s in fft box.
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of
!!    atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=number of the atom being displaced in the frozen-phonon
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=dimension of q grid for pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid.
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qgrid(mqgrid)=grid of q points from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  ucvol=unit cell volume (Bohr**3).
!!  vlspl(mqgrid,2,ntypat)=spline fit of q^2 V(q) for each type of atom.
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  vpsp1(cplex*nfft)=first-order local crystal pseudopotential in real space
!!    (including the minus sign, forgotten in the paper non-linear..
!!
!! PARENTS
!!      loop3dte,loper3,nstdy3
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vloca3(atindx,cplex,gmet,gsqcut,idir,ipert,&
& mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
& ntypat,n1,n2,n3,paral_kgb,ph1d,qgrid,qphon,ucvol,vlspl,vpsp1,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vloca3'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mqgrid,n1,n2,n3,natom,nfft,ntypat
 integer,intent(in) :: paral_kgb
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),qphon(3),vlspl(mqgrid,2,ntypat)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: vpsp1(cplex*nfft)

!Local variables -------------------------
!scalars
 integer :: i1,i2,i3,ia,ia1,iatom,id1,id2,id3,ig1,ig2,ig3,ii,ii1,im=2,isign
 integer :: itypat,jj,nri,re=1
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: aa,bb,cc,cutoff,dd,diff,dq,dq2div6,dqdiv6,dqm1,g1,g2,g3,gmag,gq1
 real(dp) :: gq2,gq3,gsq,gsquar,ph1,ph2,ph3,phi,phimag,phqim,phqre,phr
 real(dp) :: phre,qxred2pi,sfi,sfr,vion1,x1,x2,x3,xnorm,y1,y2,y3
 logical :: qeq0
!arrays
 integer :: ng(3)
 real(dp) :: gq(3)
 real(dp),allocatable :: work1(:,:)

! *********************************************************************

!Real and imaginary parts of phase--statment functions:
 phr(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 phi(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 ph1(nri,ig1,ia)=ph1d(nri,ig1+1+n1+(atindx(ia)-1)*(2*n1+1))
 ph2(nri,ig2,ia)=ph1d(nri,ig2+1+n2+(atindx(ia)-1)*(2*n2+1)+&
& natom*(2*n1+1))
 ph3(nri,ig3,ia)=ph1d(nri,ig3+1+n3+(atindx(ia)-1)*(2*n3+1)+&
& natom*(2*n1+1+2*n2+1))
 phre(ig1,ig2,ig3,ia)=phr(ph1(re,ig1,ia),ph1(im,ig1,ia),&
& ph2(re,ig2,ia),ph2(im,ig2,ia),ph3(re,ig3,ia),ph3(im,ig3,ia))
 phimag(ig1,ig2,ig3,ia)=phi(ph1(re,ig1,ia),ph1(im,ig1,ia),&
& ph2(re,ig2,ia),ph2(im,ig2,ia),ph3(re,ig3,ia),ph3(im,ig3,ia))
!
 gsq(g1,g2,g3)=g1*g1*gmet(1,1)+g2*g2*gmet(2,2)+&
& g3*g3*gmet(3,3)+2.0_dp*g1*g2*gmet(1,2)+&
& 2.0_dp*g2*g3*gmet(2,3)+2.0_dp*g3*g1*gmet(3,1)

 iatom=ipert

 if(iatom==natom+1 .or. iatom==natom+2 .or. iatom==natom+5)then

!  (In case of d/dk or an electric field)
   vpsp1(1:cplex*nfft)=zero

 else

!  (In case of a phonon perturbation)
   ABI_ALLOCATE(work1,(2,nfft))
   work1(1:2,1:nfft)=0.0_dp

   dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
   dqm1=1.0_dp/dq
   dqdiv6=dq/6.0_dp
   dq2div6=dq**2/6.0_dp
   cutoff=gsqcut*tolfix
   id1=n1/2+2
   id2=n2/2+2
   id3=n3/2+2


!  This is to allow q=0
   qeq0=.false.
   if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)qeq0=.true.

!  Determination of the atom type
   ia1=0
   itypat=0
   do ii=1,ntypat
     ia1=ia1+nattyp(ii)
     if(atindx(iatom)<=ia1.and.itypat==0)itypat=ii
   end do

!  Determination of phase qxred*
   qxred2pi=2.0_dp*pi*(qphon(1)*xred(1,iatom)+ &
&   qphon(2)*xred(2,iatom)+ &
&   qphon(3)*xred(3,iatom) )
   phqre=cos(qxred2pi)
   phqim=sin(qxred2pi)
   ii=0

   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     gq3=dble(ig3)+qphon(3)
     gq(3)=gq3
     do i2=1,n2
       if (((i2-1)/(n2/mpi_enreg%nproc_fft))==mpi_enreg%me_fft) then
         ig2=i2-(i2/id2)*n2-1
         gq2=dble(ig2)+qphon(2)
         gq(2)=gq2

!        Note the lower limit of the next loop
         ii1=1
         if(i3==1 .and. i2==1 .and. qeq0 .and. ig2==0 .and. ig3==0)then
           ii1=2
           ii=ii+1
         end if
         do i1=ii1,n1
           ig1=i1-(i1/id1)*n1-1
           gq1=dble(ig1)+qphon(1)
           gq(1)=gq1
           ii=ii+1
           gsquar=gsq(gq1,gq2,gq3)
!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Compute vion(G) for given type of atom
             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Evaluate spline fit from q^2 V(q) to get V(q):
!            (p. 86 Numerical Recipes, Press et al; NOTE error in book for sign
!            of "aa" term in derivative; also see splfit routine.
!            This bug fixed here 27 Jan 1992.)

             bb = diff*dqm1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*dq2div6
             dd = bb*(bb**2-1.0_dp)*dq2div6
             vion1 = (aa*vlspl(jj,1,itypat)+bb*vlspl(jj+1,1,itypat) + &
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat) ) &
&             / gsquar

!            Phase   G*xred  (complex conjugate) * -i *2pi*(g+q)*vion
             sfr=-phimag(ig1,ig2,ig3,iatom)*2.0_dp*pi*gq(idir)*vion1
             sfi=-phre(ig1,ig2,ig3,iatom)*2.0_dp*pi*gq(idir)*vion1

!            Phase   q*xred  (complex conjugate)
             work1(re,ii)=sfr*phqre+sfi*phqim
             work1(im,ii)=-sfr*phqim+sfi*phqre
           end if

         end do
       end if
     end do
   end do

   isign=1
   ng(1)=n1
   ng(2)=n2
   ng(3)=n3
   xnorm=1.0_dp/ucvol

!  Transform back to real space
   call fourdp(cplex,work1,vpsp1,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   vpsp1(1:cplex*nfft)=vpsp1(1:cplex*nfft)*xnorm

   ABI_DEALLOCATE(work1)

!  End the condition of non-electric-field
 end if

!DEBUG
!write(std_out,*)' vloca3 : ir,vpsp1(ir),vpsp1(ir+1)'
!do ii=1,nfft,42
!write(std_out,'(i5,2es16.6)' )ii,vpsp1(ii),vpsp1(ii+1)
!end do
!stop
!ENDDEBUG

end subroutine vloca3
!!***
