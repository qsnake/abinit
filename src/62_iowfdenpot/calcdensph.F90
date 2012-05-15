!{\src2tex{textfont=tt}}
!!****f* ABINIT/calcdensph
!! NAME
!! calcdensph
!!
!! FUNCTION
!! Compute and print integral of total density inside spheres around atoms.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  ntypat=number of atom types
!!  nunit=number of the unit for writing
!!  ratsph(ntypat)=radius of spheres around atoms
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type of each atom
!!  ucvol=unit cell volume in bohr**3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!      outscfcv,vtorho
!!
!! CHILDREN
!!      timab,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,nunit,ratsph,rhor,rprimd,typat,ucvol,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calcdensph'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,nunit
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: gmet(3,3),ratsph(ntypat),rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: ishift=5
 integer :: i1,i2,i3,iatom,ierr,ifft_local,ix,iy,iz,izmod,me_fft,n1,n1a,n1b,n2
 integer :: n2a,n2b,n3,n3a,n3b,nd3,nfftot,old_paral_level,spaceComm
 real(dp),parameter :: delta=0.99_dp
 real(dp) :: difx,dify,difz,r2,r2atsph,rr1,rr2,rr3,rx,ry,rz
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: intgden(:,:)

! *************************************************************************

 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nd3=n3/mpi_enreg%nproc_fft
 nfftot=n1*n2*n3;me_fft=mpi_enreg%me_fft
 ABI_ALLOCATE(intgden,(nspden,natom))
 intgden=zero

!Loop over atoms
!-------------------------------------------
 do iatom=1,natom

!  Define a "box" around the atom
   r2atsph=1.0000001_dp*ratsph(typat(iatom))**2
   rr1=sqrt(r2atsph*gmet(1,1))
   rr2=sqrt(r2atsph*gmet(2,2))
   rr3=sqrt(r2atsph*gmet(3,3))

   n1a=int((xred(1,iatom)-rr1+ishift)*n1+delta)-ishift*n1
   n1b=int((xred(1,iatom)+rr1+ishift)*n1      )-ishift*n1
   n2a=int((xred(2,iatom)-rr2+ishift)*n2+delta)-ishift*n2
   n2b=int((xred(2,iatom)+rr2+ishift)*n2      )-ishift*n2
   n3a=int((xred(3,iatom)-rr3+ishift)*n3+delta)-ishift*n3
   n3b=int((xred(3,iatom)+rr3+ishift)*n3      )-ishift*n3

   do i3=n3a,n3b
     iz=mod(i3+ishift*n3,n3)
     if(iz/nd3==me_fft) then
       izmod=modulo(iz,nd3)
       difz=dble(i3)/dble(n3)-xred(3,iatom)
       do i2=n2a,n2b
         iy=mod(i2+ishift*n2,n2)
         dify=dble(i2)/dble(n2)-xred(2,iatom)
         do i1=n1a,n1b
           ix=mod(i1+ishift*n1,n1)
           difx=dble(i1)/dble(n1)-xred(1,iatom)
           rx=difx*rprimd(1,1)+dify*rprimd(1,2)+difz*rprimd(1,3)
           ry=difx*rprimd(2,1)+dify*rprimd(2,2)+difz*rprimd(2,3)
           rz=difx*rprimd(3,1)+dify*rprimd(3,2)+difz*rprimd(3,3)
           r2=rx**2+ry**2+rz**2

!          Identify the fft indexes of the rectangular grid around the atom
           if (r2 <= r2atsph) then
             ifft_local=1+ix+n1*(iy+n2*izmod)
             if (nspden==1) then
!              intgden(1,iatom)= integral of total density
               intgden(1,iatom)=intgden(1,iatom)+rhor(ifft_local,1)
             else if (nspden==2) then
!              intgden(1,iatom)= integral of up density
!              intgden(1,iatom)= integral of dn density
               intgden(1,iatom)=intgden(1,iatom)+rhor(ifft_local,2)
               intgden(2,iatom)=intgden(2,iatom)+rhor(ifft_local,1)-rhor(ifft_local,2)
             else
!              intgden(1,iatom)= integral of total density
!              intgden(2,iatom)= integral of magnetization, x-component
!              intgden(3,iatom)= integral of magnetization, y-component
!              intgden(4,iatom)= integral of magnetization, z-component
               intgden(1,iatom)=intgden(1,iatom)+rhor(ifft_local,1)
               intgden(2,iatom)=intgden(2,iatom)+rhor(ifft_local,2)
               intgden(3,iatom)=intgden(3,iatom)+rhor(ifft_local,3)
               intgden(4,iatom)=intgden(4,iatom)+rhor(ifft_local,4)
             end if
           end if

         end do
       end do
     end if
   end do

   intgden(:,iatom)=intgden(:,iatom)*ucvol/dble(nfftot)

!  End loop over atoms
!  -------------------------------------------
 end do

!MPI parallelization
 if(mpi_enreg%paral_compil_fft==1)then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_fft)
   call timab(48,1,tsec)
   call xsum_mpi(intgden,spaceComm,ierr)
   call timab(48,2,tsec)
   mpi_enreg%paral_level=old_paral_level
 end if

!Printing
 write(message, '(4a)' ) ch10,&
& ' Integrated total density in atomic spheres:',ch10,&
& ' -------------------------------------------'
 call wrtout(nunit,message,'COLL')
 if (nspden==1) then
   write(message, '(a)' ) ' Atom  Sphere_radius  Integrated_density'
   call wrtout(nunit,message,'COLL')
   do iatom=1,natom
     write(message, '(i5,f15.5,f20.8)' ) iatom,ratsph(typat(iatom)),intgden(1,iatom)
     call wrtout(nunit,message,'COLL')
   end do
 else if(nspden==2) then
   write(message, '(a)' ) ' Atom  Sphere radius  Integrated_up_density  Integrated_dn_density  Total(up+dn)   Diff(up-dn)'
   call wrtout(nunit,message,'COLL')
   do iatom=1,natom
     write(message, '(i5,f15.5,2f23.8,2f14.8)' ) iatom,ratsph(typat(iatom)),intgden(1,iatom),intgden(2,iatom),&
&     (intgden(1,iatom)+intgden(2,iatom)),(intgden(1,iatom)-intgden(2,iatom))
     call wrtout(nunit,message,'COLL')
   end do
   write(message, '(3a)' ) ' Note: Diff(up-dn) can be considered as a rough ',&
&   'approximation of a local magnetic moment.'
   call wrtout(nunit,message,'COLL')
 else if(nspden==4) then
   write(message, '(a)' ) ' Atom  Sphere radius  Total_density Integrated_x_magnetiz Integrated_y_magnetiz Integrated_z_magnetiz'
   call wrtout(nunit,message,'COLL')
   do iatom=1,natom
     write(message, '(i5,f14.5,f14.8,3f20.8)' ) iatom,ratsph(typat(iatom)),(intgden(ix,iatom),ix=1,4)
     call wrtout(nunit,message,'COLL')
   end do
 end if
 write(message, '(a)' ) ch10
 call wrtout(nunit,message,'COLL')

 ABI_DEALLOCATE(intgden)

end subroutine calcdensph

!!***
