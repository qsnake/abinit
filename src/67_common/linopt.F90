!{\src2tex{textfont=tt}}
!!****f* ABINIT/linopt
!! NAME
!! linopt
!!
!! FUNCTION
!! This routine compute optical frequency dependent dielectric function
!! for semiconductors
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (SSharma)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nspin=number of spins(integer)
!!  omega=crystal volume in au (real)
!!  nkpt=total number of kpoints (integer)
!!  wkpt(nkpt)=weights of kpoints (real)
!!  nsymcrys=number of crystal symmetry operations(integer)
!!  symcrys(3,3,nsymcrys)=symmetry operations in cartisian coordinates(real)
!!  nstval=total number of valence states(integer)
!!  occv(nstval,nspin,nkpt)=occupation number for each band(real)
!!  evalv(nstval,nspin,nkpt)=eigen value for each band in Ha(real)
!!  efermi=Fermi energy in Ha(real)
!!  pmat(nstval,nstval,nkpt,3,nspin)=momentum matrix elements in cartesian coordinates(complex)
!!  v1,v2=desired component of the dielectric function(integer) 1=x,2=y,3=z
!!  nmesh=desired number of energy mesh points(integer)
!!  de=desired step in energy(real); nmesh*de=maximum energy
!!  sc=scissors shift in Ha(real)
!!  brod=broadening in Ha(real)
!!  fnam=root for filename that will contain the output filename will be trim(fnam)//'-linopt.out'
!!
!! OUTPUT
!!  Dielectric function for semiconductors, on a desired energy mesh and for a desired
!!  direction of polarisation. The output is in a file named trim(fnam)//'-linopt.out' and contains
!!  Im(\epsilon_{v1v2}(\omega), Re(\epsilon_{v1v2}(\omega) and abs(\epsilon_{v1v2}(\omega).
!!  Comment:
!!  Right now the routine sums over the kpoints. In future linear tetrahedron method should be
!!  useful.
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine linopt(nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,occv,evalv,efermi,pmat, &
  v1,v2,nmesh,de,sc,brod,fnam)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'linopt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!no_abirules
integer, intent(in) :: nspin
real(dp), intent(in) :: omega
integer, intent(in) :: nkpt
real(dp), intent(in) :: wkpt(nkpt)
integer, intent(in) :: nsymcrys
real(dp), intent(in) :: symcrys(3,3,nsymcrys)
integer, intent(in) :: nstval
real(dp), intent(in) :: occv(nstval,nspin,nkpt)
real(dp), intent(in) :: evalv(nstval,nspin,nkpt)
real(dp), intent(in) :: efermi
complex(dpc), intent(in) :: pmat(nstval,nstval,nkpt,3,nspin)
integer, intent(in) :: v1
integer, intent(in) :: v2
integer, intent(in) :: nmesh
real(dp), intent(in) :: de
real(dp), intent(in) :: sc
real(dp), intent(in) :: brod
character(256), intent(in) :: fnam

!Local variables -------------------------
!no_abirules
integer :: isp
integer :: i,j,isym,lx,ly,ik
integer :: ist1,ist2,iw
real(dp) :: e1,e2,e12,deltav1v2
real(dp) :: ha2ev
real(dp) :: renorm_factor,emin,emax
real(dp) :: ene
complex(dpc) :: b11,b12
complex(dpc) :: ieta,w
character(256) :: fnam1
! local allocatable arrays
real(dp), allocatable :: s(:,:)
real(dp), allocatable :: sym(:,:)
complex(dpc), allocatable :: chi(:)
complex(dpc), allocatable :: eps(:)

! *********************************************************************

!fool proof:
!check polarisation
 if (v1.le.0.or.v2.le.0.or.v1.gt.3.or.v2.gt.3) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linopt:                         '
   write(std_out,*) '    the polarisation directions incorrect    '
   write(std_out,*) '    1=x and 2=y and 3=z                      '
   write(std_out,*) '---------------------------------------------'
   stop
 end if
!number of energy mesh points
 if (nmesh.le.0) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linopt:                         '
   write(std_out,*) '    number of energy mesh points incorrect   '
   write(std_out,*) '    number has to integer greater than 0     '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   stop
 end if
!step in energy
 if (de.le.0._dp) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linopt:                         '
   write(std_out,*) '    energy step is incorrect                 '
   write(std_out,*) '    number has to real greater than 0.0      '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   stop
 end if
!broadening
 if (brod.gt.0.009) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    ATTENTION: broadening is quite high      '
   write(std_out,*) '    ideally should be less than 0.005        '
   write(std_out,*) '---------------------------------------------'
 else if (brod.gt.0.015) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: broadening is too high   '
   write(std_out,*) '    ideally should be less than 0.005   '
   write(std_out,*) '----------------------------------------'
 end if
!fermi energy
 if(efermi<-1.0d4) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    ATTENTION: Fermi energy seems extremely  '
   write(std_out,*) '    low                                      '
   write(std_out,*) '---------------------------------------------'
 end if
!scissors operator
 if (sc.lt.0._dp) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linopt:                         '
   write(std_out,*) '    scissors shift is incorrect              '
   write(std_out,*) '    number has to real greater than 0.0      '
   write(std_out,*) '---------------------------------------------'
   stop
 end if
!fool proof end
!
!allocate local arrays
 ABI_ALLOCATE(chi,(nmesh))
 ABI_ALLOCATE(eps,(nmesh))
 ABI_ALLOCATE(s,(3,3))
 ABI_ALLOCATE(sym,(3,3))
 ieta=(0._dp,1._dp)*brod
 renorm_factor=1._dp/(omega*dble(nsymcrys))
 ha2ev=13.60569172*2._dp
!output file names
 fnam1=trim(fnam)//'-linopt.out'
!construct symmetrisation tensor
 sym(:,:)=0._dp
 do isym=1,nsymcrys
   s(:,:)=symcrys(:,:,isym)
   do i=1,3
     do j=1,3
       sym(i,j)=sym(i,j)+s(i,v1)*s(j,v2)
     end do
   end do
 end do
!calculate the energy window
 emin=0._dp
 emax=0._dp
 do ik=1,nkpt
   do isp=1,nspin
     do ist1=1,nstval
       emin=min(emin,evalv(ist1,isp,ik))
       emax=max(emax,evalv(ist1,isp,ik))
     end do
   end do
 end do
!start calculating linear optical response
 chi(:)=0._dp
 do ik=1,nkpt
   write(std_out,*) ik,'of',nkpt
   do isp=1,nspin
     do ist1=1,nstval
       e1=evalv(ist1,isp,ik)
!      if (e1.lt.efermi) then
!      do ist2=ist1,nstval
       do ist2=1,nstval
         e2=evalv(ist2,isp,ik)
!        if (e2.gt.efermi) then
         if (ist1.ne.ist2) then
!          scissors correction of momentum matrix
           e12=e1-e2-sc
           b11=0._dp
!          symmetrization of momentum matrix
           do lx=1,3
             do ly=1,3
               b11=b11+(sym(lx,ly)*pmat(ist1,ist2,ik,lx,isp)* &
               conjg(pmat(ist1,ist2,ik,ly,isp)))
             end do
           end do
           b12=b11*renorm_factor*(1._dp/(e12**2))
!          calculate on the desired energy grid
           do iw=2,nmesh
             w=(iw-1)*de+ieta
             chi(iw)=chi(iw)+(wkpt(ik)*(occv(ist1,isp,ik)-occv(ist2,isp,ik))* &
             (b12/(-e12-w)))
           end do
!          end loops over states
         end if
       end do
!      end if
     end do
!    end loop over spins
   end do
!  end loop over k-points
 end do

!open the output files
 open(92,file=fnam1,action='WRITE',form='FORMATTED')
!write the output
 write(92, '(a)' ) ' # Energy(eV)         Im(eps(w))'
 write(92, '(a,2i3,a)' )' #calculated the component:',v1,v2,'  of dielectric function'
 write(std_out,*) 'calculated the component:',v1,v2,'  of dielectric function'
 write(92, '(a,2es16.6)' ) ' #broadening:', real(ieta),aimag(ieta)
 write(std_out,*) ' with broadening:',ieta
 write(92, '(a,es16.6)' ) ' #scissors shift:',sc
 write(std_out,*) 'and scissors shift:',sc
 write(92, '(a,es16.6,a,es16.6,a)' ) ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
 write(std_out,*) 'energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
 eps(:)=0._dp
 deltav1v2=zero
 if(v1==v2)deltav1v2=one
 do iw=2,nmesh
   ene=(iw-1)*de
   ene=ene*ha2ev
   eps(iw)=deltav1v2+4._dp*pi*chi(iw)
   write(92, '(2es16.6)' ) ene,aimag(eps(iw))
 end do
 write(92,*)
 write(92,*)
 write(92, '(a)' ) ' # Energy(eV)         Re(eps(w))'
 do iw=2,nmesh
   ene=(iw-1)*de
   ene=ene*ha2ev
   write(92, '(2es16.6)' ) ene,dble(eps(iw))
 end do
 write(92,*)
 write(92,*)
 write(92, '(a)' )' # Energy(eV)         abs(eps(w))'
 do iw=2,nmesh
   ene=(iw-1)*de
   ene=ene*ha2ev
   write(92, '(2es16.6)' ) ene,abs(eps(iw))
 end do

!close output file
 close(92)
!deallocate local arrays
 ABI_DEALLOCATE(s)
 ABI_DEALLOCATE(sym)
 ABI_DEALLOCATE(chi)
 ABI_DEALLOCATE(eps)

 return

end subroutine linopt
!!***
