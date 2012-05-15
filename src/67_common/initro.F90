!{\src2tex{textfont=tt}}
!!****f* ABINIT/initro
!!
!! NAME
!! initro
!!
!! FUNCTION
!! Initialize the density using either:
!!  - a gaussian of adjustable decay length (norm-conserving psp)
!!  - PS atomic valence density from psp file (PAW)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA,XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! atindx(natom)=index table for atoms (see scfcv.f)
!! densty(ntypat,4)=parameters for initialisation of the density of each atom type
!! gmet(3,3)=reciprocal space metric (Bohr**-2)
!! gsqcut=cutoff G**2 for included G s in fft box (larger sphere).
!! izero=if 1, unbalanced components of rho(g) have to be set to zero
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=informations about mpi parallelization
!! mqgrid=number of grid pts in q array for n^AT(q) spline.
!! natom=number of atoms in cell.
!! nattyp(ntypat)=number of atoms of each type in cell.
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! ntypat=number of types of atoms in cell.
!! nspden=number of spin-density components
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase information for given atom coordinates.
!! qgrid(mqgrid)=q grid for spline atomic valence density n^AT(q) from 0 to qmax.
!! spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!! ucvol=unit cell volume (Bohr**3).
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! zion(ntypat)=charge on each type of atom (real number)
!! znucl(ntypat)=atomic number, for each type of atom
!!
!! OUTPUT
!! rhog(2,nfft)=initialized total density in reciprocal space
!! rhor(nfft,nspden)=initialized total density in real space.
!!         as well as spin-up part if spin-polarized
!!
!! PARENTS
!!      gstate,setup_positron
!!
!! CHILDREN
!!      atmlength,fourdp,leave_new,wrtout,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initro(atindx,densty,gmet,gsqcut,izero,mgfft,mpi_enreg,mqgrid,natom,nattyp,&
&  nfft,ngfft,nspden,ntypat,paral_kgb,pawtab,ph1d,qgrid,rhog,rhor,spinat,ucvol,usepaw,zion,znucl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initro'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,mgfft,mqgrid,natom,nfft,nspden,ntypat,paral_kgb
 integer,intent(in) :: usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(mpi_type) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: densty(ntypat,4),gmet(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),spinat(3,natom),zion(ntypat)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(out) :: rhog(2,nfft),rhor(nfft,nspden)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!The decay lengths should be optimized element by element, and even
!pseudopotential by pseudopotential.
!scalars
 integer,parameter :: im=2,natnum=110,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ig1,ig2,ig3,ii,ispden
 integer :: itypat,jj,jtemp,me_fft,n1,n2,n3,nproc_fft,nri
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: aa,alf2pi2,bb,cc,cutoff,dd,diff,dq,dq2div6,dqm1,fact,fact0,gmag
 real(dp) :: gsq,gsquar,ph1,ph2,ph3,phi,phimag,phr,phre,rhoat,sfi,sfr,x1,x2,x3
 real(dp) :: xnorm,y1,y2,y3
 character(len=500) :: message
!arrays
 real(dp),allocatable :: length(:),spinat_indx(:,:),work(:)
 logical,allocatable :: use_gaussian(:)

! *************************************************************************

!Real and imaginary parts of phase--statment functions:
 phr(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 phi(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 ph1(nri,ig1,ia)=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 ph2(nri,ig2,ia)=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+&
& natom*(2*n1+1))
 ph3(nri,ig3,ia)=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+&
& natom*(2*n1+1+2*n2+1))
 phre(ig1,ig2,ig3,ia)=phr(ph1(re,ig1,ia),ph1(im,ig1,ia),&
& ph2(re,ig2,ia),ph2(im,ig2,ia),ph3(re,ig3,ia),ph3(im,ig3,ia))
 phimag(ig1,ig2,ig3,ia)=phi(ph1(re,ig1,ia),ph1(im,ig1,ia),&
& ph2(re,ig2,ia),ph2(im,ig2,ia),ph3(re,ig3,ia),ph3(im,ig3,ia))

 gsq(i1,i2,i3)=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
& dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
& dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)

!DEBUG
!write(std_out,*)' initro : enter '
!ENDDEBUG

 if(nspden==4)then
   write(std_out,*)' initro : might work yet for nspden=4 (not checked)'
   write(std_out,*)'spinat',spinat(1:3,1:natom)
!  stop
 end if

!Check whether the values of spinat are acceptable
 if(nspden==2)then
   ia1=1
   do itypat=1,ntypat
!    ia1,ia2 sets range of loop over atoms:
     ia2=ia1+nattyp(itypat)-1
     do ia=ia1,ia2
       if( sqrt(spinat(1,ia)**2+spinat(2,ia)**2+spinat(3,ia)**2) &
&       > abs(zion(itypat))*(1.0_dp + epsilon(0.0_dp)) ) then
         write(message, '(a,a,a,a,i4,a,a,3es11.4,a,a,a,es11.4)' ) ch10,&
&         ' initro : WARNING - ',ch10,&
&         '  For atom number ',ia,ch10,&
&         '  input spinat=',spinat(:,ia),'  is larger, in magnitude,',ch10,&
&         '  than zion(ia)=',zion(itypat)
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
       end if
     end do
     ia1=ia2+1
   end do
 end if

 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)
 ABI_ALLOCATE(work,(nfft))
 ABI_ALLOCATE(spinat_indx,(3,natom))

!Transfer the spinat array to an array in which the atoms have the
!proper order, type by type.
 do ia=1,natom
   spinat_indx(:,atindx(ia))=spinat(:,ia)
 end do

!Compute the decay length of each type of atom
 ABI_ALLOCATE(length,(ntypat))
 ABI_ALLOCATE(use_gaussian,(ntypat))
 jtemp=0
 do itypat=1,ntypat

   use_gaussian(itypat)=.true.
   if (usepaw==1)then
     use_gaussian(itypat)=(pawtab(itypat)%usetvale==0)
   end if
   if (.not.use_gaussian(itypat)) jtemp=jtemp+1

   if (use_gaussian(itypat)) then
     call atmlength(densty(itypat,1),length(itypat),zion(itypat),znucl(itypat))
     write(message,'(a,i3,a,f12.4,a,a,a,f12.4,a,i3,a,es12.4,a)' )&
&     ' initro : for itypat=',itypat,', take decay length=',length(itypat),',',&
&     ch10,' initro : indeed, coreel=',znucl(itypat)-zion(itypat),', nval=',int(zion(itypat)),&
&     ' and densty=',densty(itypat,1),'.'
     call wrtout(std_out,message,'COLL')
   end if

 end do

 if (jtemp>0) then
   dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
   dqm1=1.0_dp/dq
   dq2div6=dq**2/6.0_dp
 end if

 cutoff=gsqcut*tolfix
 xnorm=1.0_dp/ucvol

 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 if(nspden /= 4) then

   do ispden=nspden,1,-1
!    This loop overs spins will actually be as follows :
!    ispden=2 for spin up
!    ispden=1 for total spin (also valid for non-spin-polarized calculations)
!    The reverse ispden order is chosen, in order to end up with
!    rhog containing the proper total density.

     rhog(:,:)=zero

     ia1=1
     do itypat=1,ntypat

       if (use_gaussian(itypat)) alf2pi2=(two_pi*length(itypat))**2

!      ia1,ia2 sets range of loop over atoms:
       ia2=ia1+nattyp(itypat)-1
       ii=0
       jtemp=0

       do i3=1,n3
         ig3=i3-(i3/id3)*n3-1
         do i2=1,n2
           ig2=i2-(i2/id2)*n2-1
!          if (me_fft==modulo(i2,nproc_fft)) then
           if (((i2-1)/(n2/nproc_fft))==me_fft) then
             do i1=1,n1

               ig1=i1-(i1/id1)*n1-1
               ii=ii+1
               gsquar=gsq(ig1,ig2,ig3)

!              Skip G**2 outside cutoff:
               if (gsquar<=cutoff) then

!                Assemble structure factor over all atoms of given type,
!                also taking into account the spin-charge on each atom:
                 sfr=zero;sfi=zero
                 if(ispden==1)then
                   do ia=ia1,ia2
                     sfr=sfr+phre(ig1,ig2,ig3,ia)
                     sfi=sfi-phimag(ig1,ig2,ig3,ia)
                   end do
                   if (use_gaussian(itypat)) then
                     sfr=sfr*zion(itypat)
                     sfi=sfi*zion(itypat)
                   end if
                 else
                   fact0=half;if (.not.use_gaussian(itypat)) fact0=half/zion(itypat)
                   do ia=ia1,ia2
!                    Here, take care only of the z component
                     fact=fact0*(zion(itypat)+spinat_indx(3,ia))
                     sfr=sfr+phre(ig1,ig2,ig3,ia)*fact
                     sfi=sfi-phimag(ig1,ig2,ig3,ia)*fact
                   end do
                 end if

!                Charge density integrating to one
                 if (use_gaussian(itypat)) then
                   rhoat=xnorm*exp(-gsquar*alf2pi2)
!                  Multiply structure factor times rhoat (atomic density in reciprocal space)
                   rhog(re,ii)=rhog(re,ii)+sfr*rhoat
                   rhog(im,ii)=rhog(im,ii)+sfi*rhoat        
                 else if(usepaw/=0)then
                   gmag=sqrt(gsquar)
                   jj=1+int(gmag*dqm1)
                   diff=gmag-qgrid(jj)
                   bb = diff*dqm1
                   aa = one-bb
                   cc = aa*(aa**2-one)*dq2div6
                   dd = bb*(bb**2-one)*dq2div6
                   rhoat=(aa*pawtab(itypat)%tvalespl(jj,1)+bb*pawtab(itypat)%tvalespl(jj+1,1)+&
&                   cc*pawtab(itypat)%tvalespl(jj,2)+dd*pawtab(itypat)%tvalespl(jj+1,2))&
&                   *xnorm
!                  Multiply structure factor times rhoat (atomic density in reciprocal space)
                   rhog(re,ii)=rhog(re,ii)+sfr*rhoat
                   rhog(im,ii)=rhog(im,ii)+sfi*rhoat
                 else
                   write(message,'(4a)')ch10,&
&                   ' initro : BUG -',ch10,&
&                   '  Initialization of density is non consistent.'
                   call wrtout(std_out,message,'COLL')
                   call leave_new('COLL')
                 end if

               else
                 jtemp=jtemp+1
               end if

!              End loop on i1
             end do
           end if
!          End loop on i2
         end do
!        End loop on i3
       end do
       ia1=ia2+1

!      End loop on type of atoms
     end do

!    Set contribution of unbalanced components to zero
     if (izero==1) call zerosym(rhog,2,mpi_enreg,n1,n2,n3)

!    Note, we end with ispden=1, so that rhog contains the total density
     call fourdp(1,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     rhor(:,ispden)=work(:)
!    End loop on spins
   end do

 else if(nspden==4) then
   do ispden=nspden,1,-1
!    This loop overs spins will actually be as follows :
!    ispden=2,3,4 for mx,my,mz
!    ispden=1 for total spin (also valid for non-spin-polarized calculations)
!    The reverse ispden order is chosen, in order to end up with
!    rhog containing the proper total density.

     rhog(:,:)=zero

     ia1=1
     do itypat=1,ntypat

       if (use_gaussian(itypat)) alf2pi2=(two_pi*length(itypat))**2

!      ia1,ia2 sets range of loop over atoms:
       ia2=ia1+nattyp(itypat)-1
       ii=0
       jtemp=0
       do i3=1,n3
         ig3=i3-(i3/id3)*n3-1
         do i2=1,n2
           ig2=i2-(i2/id2)*n2-1
!          if (me_fft==modulo(i2,nproc_fft)) then
           if (((i2-1)/(n2/nproc_fft))==me_fft) then
             do i1=1,n1

               ig1=i1-(i1/id1)*n1-1
               ii=ii+1
               gsquar=gsq(ig1,ig2,ig3)

!              Skip G**2 outside cutoff:
               if (gsquar<=cutoff) then

!                Assemble structure factor over all atoms of given type,
!                also taking into account the spin-charge on each atom:
                 sfr=zero;sfi=zero
                 if(ispden==1)then
                   do ia=ia1,ia2
                     sfr=sfr+phre(ig1,ig2,ig3,ia)
                     sfi=sfi-phimag(ig1,ig2,ig3,ia)
                   end do
                   if (use_gaussian(itypat)) then
                     sfr=sfr*zion(itypat)
                     sfi=sfi*zion(itypat)
                   end if
                 else
                   fact0=one;if (.not.use_gaussian(itypat)) fact0=one/zion(itypat)
                   do ia=ia1,ia2
!                    Here, take care of the components of m
                     fact=fact0*spinat_indx(ispden-1,ia)
                     sfr=sfr+phre(ig1,ig2,ig3,ia)*fact
                     sfi=sfi-phimag(ig1,ig2,ig3,ia)*fact
                   end do
                 end if

!                Charge density integrating to one
                 if (use_gaussian(itypat)) then
                   rhoat=xnorm*exp(-gsquar*alf2pi2)
                 else if(usepaw/=0)then
                   gmag=sqrt(gsquar)
                   jj=1+int(gmag*dqm1)
                   diff=gmag-qgrid(jj)
                   bb = diff*dqm1
                   aa = one-bb
                   cc = aa*(aa**2-one)*dq2div6
                   dd = bb*(bb**2-one)*dq2div6
                   rhoat=(aa*pawtab(itypat)%tvalespl(jj,1)+bb*pawtab(itypat)%tvalespl(jj+1,1)+&
&                   cc*pawtab(itypat)%tvalespl(jj,2)+dd*pawtab(itypat)%tvalespl(jj+1,2))&
&                   *xnorm
                 else
                   write(message,'(4a)')ch10,&
&                   ' initro : BUG -',ch10,&
&                   '  Initialization of density is non consistent.'
                   call wrtout(std_out,message,'COLL') 
                   call leave_new('COLL')
                 end if

!                Multiply structure factor times rhoat (atomic density in reciprocal space)
                 rhog(re,ii)=rhog(re,ii)+sfr*rhoat
                 rhog(im,ii)=rhog(im,ii)+sfi*rhoat
               else
                 jtemp=jtemp+1
               end if

!              End loop on i1
             end do
           end if
!          End loop on i2
         end do
!        End loop on i3
       end do
       ia1=ia2+1

!      End loop on type of atoms
     end do

!    Set contribution of unbalanced components to zero
     if (izero==1) call zerosym(rhog,2,mpi_enreg,n1,n2,n3)

!    Note, we end with ispden=1, so that rhog contains the total density
     call fourdp(1,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     rhor(:,ispden)=work(:)

!    End loop on spins
   end do
 end if
!DEBUG
!do ir=1,nfft,30
!do ir=1,nfft/(2*nproc_fft*20),24
!write(std_out,*)'rhor',ir,rhor(ir,:)
!end do
!stop
!ENDDEBUG

 ABI_DEALLOCATE(length)
 ABI_DEALLOCATE(use_gaussian)
 ABI_DEALLOCATE(spinat_indx)
 ABI_DEALLOCATE(work)

!DEBUG
!write(std_out,*)' initro : exit '
!ENDDEBUG

end subroutine initro
!!***
