!{\src2tex{textfont=tt}}
!!****f* ABINIT/sphere_fft
!! NAME
!! sphere_fft
!!
!! FUNCTION
!! Array cg is defined in sphere with npw points. Insert cg inside box
!! of n1*n2*n3 points to define array cfft for fft box.
!! corresponds to given element in cg.  rest of cfft is filled with 0 s.
!!
!! iflag=1==>insert cg into cfft.
!! iflag=2==>insert cg into cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!! iflag=-1==> extract cg from cfft.
!! iflag=-2==> extract cg from cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!!  (WARNING : iflag=-2 cannot use symmetry operations)
!!
!! There is also the possibility to apply a symmetry operation,
!! as well as to make a shift in reciprocal space, or to multiply
!! by a constant factor, in the case iflag=-1.
!! Multiplication by a constant factor is also possible in the case iflag=-2.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, AR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cg(2,npw)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to do in //
!! npw=number of G vectors in basis at this k point
!! cfft(2,n4,n5,n6) = fft box
!! n1,n2,n3=physical dimension of the box (cfft)
!! n4,n5,n6=memory dimension of cfft
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! mpi_enreg=informations about MPI parallelization
!! nd2proc TO BE DESCRIBED SB 090831
!! iflag=option parameter. Possible values: -1, -2, 1, 2 ; this is used only in debug option
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! iflag=1 and 2, insert cg(input) into cfft(output)
!! iflag=-1 and -2, extract cg(output) from cfft(input)
!!
!! NOTES
!! cg and cfft are assumed to be of type COMPLEX, although this routine treats
!! them as real of twice the length to avoid nonstandard complex*16.
!!
!! WARNING
!! NO CHECK is DONE over iflag.
!!
!! TODO
!! Order arguments
!!
!! PARENTS
!!      fourwf
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine sphere_fft(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,kg_k,&
& mpi_enreg,nd2proc)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sphere_fft'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,n4,n5,nd2proc,ndat,npw
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(inout) :: cfft(2,n4,n5,nd2proc*ndat),cg(2,npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,idat,ipw,nproc_fft
!arrays

! *************************************************************************

!DEBUG
!write(std_out,*)' sphere : iflag=',iflag
!ENDDEBUG

!!Insert cg into cfft with extra 0 s around outside:
!!$OMP PARALLEL DO PRIVATE(i1,i2,i3) SHARED(cfft,ndat,n1,n2,n3)
!!do i2=1,nd2proc*ndat
!!do i3=1,n3
!!do i1=1,n1
!!cfft(1,i1,i3,i2)=0.0d0
!!cfft(2,i1,i3,i2)=0.0d0
!!end do
!!end do
!!end do
!!$OMP END PARALLEL DO
 nproc_fft = mpi_enreg%nproc_fft 

 cfft(:,:,:,:)=zero

!$OMP PARALLEL DO PRIVATE(i1,i2,i3,idat,ipw) SHARED(cfft,cg,kg_k,ndat,n1,n2,n3,npw,nproc_fft,nd2proc)
!write(std_out,*)'In sphere fft,i1,i2,i3,ipw,cfft=cg'
 do ipw=1,npw
   i1=kg_k(1,ipw); if(i1<0)i1=i1+n1 ; i1=i1+1
   i2=kg_k(2,ipw); if(i2<0)i2=i2+n2 ; i2=i2+1
   i3=kg_k(3,ipw); if(i3<0)i3=i3+n3 ; i3=i3+1
   do idat=1,ndat
     cfft(1,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))=cg(1,ipw+npw*(idat-1))
     cfft(2,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))=cg(2,ipw+npw*(idat-1))
   end do
!  write(std_out,'(4i3,2e24.12)')kg_k(:,ipw),ipw,cg(:,ipw)
 end do
!$OMP END PARALLEL DO

end subroutine sphere_fft
!!***

!!****f* ABINIT/sphere_fft1
!! NAME
!! sphere_fft1
!!
!! FUNCTION
!! Array cg is defined in sphere with npw points. Insert cg inside box
!! of n1*n2*n3 points to define array cfft for fft box.
!! corresponds to given element in cg.  rest of cfft is filled with 0 s.
!!
!! iflag=1==>insert cg into cfft.
!! iflag=2==>insert cg into cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!! iflag=-1==> extract cg from cfft.
!! iflag=-2==> extract cg from cfft, where the second and third dimension
!! have been switched (needed for new 2002 SGoedecker FFT)
!!  (WARNING : iflag=-2 cannot use symmetry operations)
!!
!! There is also the possibility to apply a symmetry operation,
!! as well as to make a shift in reciprocal space, or to multiply
!! by a constant factor, in the case iflag=-1.
!! Multiplication by a constant factor is also possible in the case iflag=-2.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, AR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! cg(2,npw)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to do in //
!! npw=number of G vectors in basis at this k point
!! cfft(2,n4,n5,n6) = fft box
!! n1,n2,n3=physical dimension of the box (cfft)
!! n4,n5,n6=memory dimension of cfft
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! mpi_enreg=informations about MPI parallelization
!! nd2proc TO BE DESCRIBED SB 090831
!! iflag=option parameter. Possible values: -1, -2, 1, 2
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! iflag=1 and 2, insert cg(input) into cfft(output)
!! iflag=-1 and -2, extract cg(output) from cfft(input)
!!
!! NOTES
!! cg and cfft are assumed to be of type COMPLEX, although this routine treats
!! them as real of twice the length to avoid nonstandard complex*16.
!!
!! WARNING
!! NO CHECK is DONE over iflag.
!!
!! TODO
!! Order arguments
!!
!! PARENTS
!!      sg_fourwf
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sphere_fft1(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,n6,kg_k,&
& mpi_enreg,nd2proc)

 use m_profiling

 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sphere_fft1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,n4,n5,n6,nd2proc,ndat,npw
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(inout) :: cfft(2,n4,n5,n6),cg(2,npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,idat,ipw,nproc_fft
!arrays

! *************************************************************************

!DEBUG
!write(std_out,*)' sphere : iflag=',iflag
!ENDDEBUG

!Insert cg into cfft with extra 0 s around outside:
!!$OMP PARALLEL DO PRIVATE(i1,i2,i3) SHARED(cfft,ndat,n1,n2,n3)
!!do i2=1,nd2proc*ndat
!!do i3=1,n3
!!do i1=1,n1
!!cfft(1,i1,i3,i2)=0.0d0
!!cfft(2,i1,i3,i2)=0.0d0
!!end do
!!end do
!!end do
!!$OMP END PARALLEL DO
 nproc_fft = mpi_enreg%nproc_fft

 cfft(:,:,:,:)=zero
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,idat,ipw) SHARED(cfft,cg,kg_k,ndat,n1,n2,n3,npw,nproc_fft,nd2proc)
!write(std_out,*)'In sphere fft,i1,i2,i3,ipw,cfft=cg'
 do ipw=1,npw
   i1=kg_k(1,ipw); if(i1<0)i1=i1+n1 ; i1=i1+1
   i2=kg_k(2,ipw); if(i2<0)i2=i2+n2 ; i2=i2+1
   i3=kg_k(3,ipw); if(i3<0)i3=i3+n3 ; i3=i3+1
   do idat=1,ndat
     cfft(1,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))=cg(1,ipw+npw*(idat-1))
     cfft(2,i1,i3,(i2-1)/nproc_fft +1+nd2proc*(idat-1))=cg(2,ipw+npw*(idat-1))
   end do
!  write(std_out,'(4i3,2e24.12)')kg_k(:,ipw),ipw,cg(:,ipw)
!  write(std_out,'(5i3,2e24.12)')i1,i2,(i2-1)/mpi_enreg%nproc_fft +1,i3,ipw,cg(:,ipw)
 end do
!$OMP END PARALLEL DO

!call leave_new("COLL")

end subroutine sphere_fft1
!!***
