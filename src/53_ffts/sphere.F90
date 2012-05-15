!{\src2tex{textfont=tt}}
!!****f* ABINIT/sphere
!! NAME
!! sphere
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
!! cg(2,npw*ndat)= contains values for npw G vectors in basis sphere
!! ndat=number of FFT to do in //
!! npw=number of G vectors in basis at this k point
!! cfft(2,n4,n5,n6) = fft box
!! n1,n2,n3=physical dimension of the box (cfft)
!! n4,n5,n6=memory dimension of cfft
!! kg_k(3,npw)=integer coordinates of G vectors in basis sphere
!! istwf_k=option parameter that describes the storage of wfs
!! iflag=option parameter. Possible values: -1, -2, 1, 2
!! mpi_enreg=Information about MPI parallelization.
!! shiftg(3)=The shift in reciprocal space.
!! symm(3,3)=symmetry operation in reciprocal space to be applied.
!! xnorm=Normalization factor.
!!
!! SIDE EFFECTS
!! Input/Output
!! iflag=1 and 2, insert cg(input) into cfft(output)
!! iflag=-1 and -2, extract cg(output) from cfft(input)
!!
!! NOTES
!! cg and cfft are assumed to be of type COMPLEX, although this routine treats
!! them as real of twice the length to avoid nonstandard complex*16.
!! If istwf_k differs from 1, then special storage modes must be taken
!! into account, for symmetric wavefunctions coming from k=(0 0 0) or other
!! special k points.
!!
!! TODO
!! Order arguments
!!
!! PARENTS
!!      fftw3_fourwf,fourwf,m_io_kss,m_wfs,sg_fourwf,wfconv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine sphere(cg,ndat,npw,cfft,n1,n2,n3,n4,n5,n6,kg_k,istwf_k,iflag,mpi_enreg,shiftg,symm,xnorm)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sphere'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iflag,istwf_k,n1,n2,n3,n4,n5,n6,ndat,npw
 real(dp),intent(in) :: xnorm
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw),shiftg(3),symm(3,3)
 real(dp),intent(inout) :: cfft(2,n4,n5,n6*ndat),cg(2,npw*ndat)

!Local variables-------------------------------
!scalars
 integer :: i1,i1inv,i2,i2inv,i3,i3inv,id1,id2,id3,idat,ipw
 integer :: j1,j2,j3,l1,l2,l3,npwmin,use_symmetry
 character(len=500) :: msg
!arrays
 integer :: identity(3,3)
 integer,allocatable :: i1inver(:),i2inver(:),i3inver(:)

! *************************************************************************

!In the case of special k-points, invariant under time-reversal,
!but not Gamma, initialize the inverse coordinates
!Remember indeed that
!u_k(G) = u_{k+G0}(G-G0); u_{-k}(G) = u_k(G)^*
!and therefore:
!u_{G0/2}(G) = u_{G0/2}(-G-G0)^*.

 if(istwf_k>=2)then
   ABI_ALLOCATE(i1inver,(n1))
   ABI_ALLOCATE(i2inver,(n2))
   ABI_ALLOCATE(i3inver,(n3))
   if(istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8)then
     i1inver(1)=1
     do i1=2,n1
       i1inver(i1)=n1+2-i1
     end do
   else
     do i1=1,n1
       i1inver(i1)=n1+1-i1
     end do
   end if
   if(istwf_k>=2 .and. istwf_k<=5)then
     i2inver(1)=1
     do i2=2,n2
       i2inver(i2)=n2+2-i2
     end do
   else
     do i2=1,n2
       i2inver(i2)=n2+1-i2
     end do
   end if
   if(istwf_k==2 .or. istwf_k==3 .or. istwf_k==6 .or. istwf_k==7)then
     i3inver(1)=1
     do i3=2,n3
       i3inver(i3)=n3+2-i3
     end do
   else
     do i3=1,n3
       i3inver(i3)=n3+1-i3
     end do
   end if
 end if

 if (iflag==1 .or. iflag==2) then

   if(iflag==1)then ! Insert cg into cfft with extra 0 s around outside:
!    $omp parallel do private(i1,i2,i3) shared(cfft,ndat,n1,n2,n3)
     do i3=1,n3*ndat
       do i2=1,n2
         do i1=1,n1
           cfft(1,i1,i2,i3)=0.0d0
           cfft(2,i1,i2,i3)=0.0d0
         end do
       end do
     end do
   end if

   if(iflag==2)then ! Insert cg into cfft with extra 0 s around outside:
!    $omp parallel do private(i1,i2,i3) shared(cfft,ndat,n1,n2,n3)
     do i2=1,n2*ndat
       do i3=1,n3
         do i1=1,n1
           cfft(1,i1,i3,i2)=0.0d0
           cfft(2,i1,i3,i2)=0.0d0
         end do
       end do
     end do
   end if

!  Take care of each plane wave, and complete cfft if needed

   if(istwf_k==1)then

     if(iflag==1)then
!      $omp parallel do private(i1,i2,i3,idat,ipw) shared(cfft,cg,kg_k,ndat,npw)
       do ipw=1,npw
         i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
         i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
         i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
         do idat=1,ndat
           cfft(1,i1,i2,i3+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i2,i3+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

     if(iflag==2)then
!      $omp parallel do private(i1,i2,i3,idat,ipw) shared(cfft,cg,kg_k,ndat,npw)
       do ipw=1,npw
         i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
         i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
         i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
         do idat=1,ndat
           cfft(1,i1,i3,i2+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i3,i2+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
         end do
       end do
     end if

   else if(istwf_k>=2)then

     npwmin=1
     if(istwf_k==2 .and. mpi_enreg%me_g0==1)then ! If gamma point, then cfft must be completed
       do idat=1,ndat
         cfft(1,1,1,1+n6*(idat-1))=cg(1,1+npw*(idat-1))
         cfft(2,1,1,1)=0.0d0
       end do
       npwmin=2
     end if

     if(iflag==1)then
!      $OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipw) &
!      $OMP&SHARED(cfft,cg,i1inver,i2inver,i3inver,kg_k,npw,npwmin)
       do ipw=npwmin,npw

         i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
         i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
         i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1
!        Construct the coordinates of -k-G
         i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)
         do idat=1,ndat
           cfft(1,i1,i2,i3+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i2,i3+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
           cfft(1,i1inv,i2inv,i3inv+n6*(idat-1))= cg(1,ipw+npw*(idat-1))
           cfft(2,i1inv,i2inv,i3inv+n6*(idat-1))=-cg(2,ipw+npw*(idat-1))
         end do
!        
       end do
     end if

     if(iflag==2)then
!      $OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipw) &
!      $OMP&SHARED(cfft,cg,i1inver,i2inver,i3inver,kg_k,npw,npwmin)
       do ipw=npwmin,npw

         i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
         i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
         i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

!        Construct the coordinates of -k-G
         i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)
         do idat=1,ndat
           cfft(1,i1,i3,i2+n6*(idat-1))=cg(1,ipw+npw*(idat-1))
           cfft(2,i1,i3,i2+n6*(idat-1))=cg(2,ipw+npw*(idat-1))
           cfft(1,i1inv,i3inv,i2inv+n6*(idat-1))= cg(1,ipw+npw*(idat-1))
           cfft(2,i1inv,i3inv,i2inv+n6*(idat-1))=-cg(2,ipw+npw*(idat-1))
         end do
!        
       end do
     end if

   end if

 else if (iflag==-1 .or. iflag==-2) then

   use_symmetry=0
   identity(:,:)=0
   identity(1,1)=1 ; identity(2,2)=1 ; identity(3,3)=1
   if(sum((symm(:,:)-identity(:,:))**2)/=0)use_symmetry=1
   if(sum(shiftg(:)**2)/=0)use_symmetry=1
!  Extract cg from cfft, ignoring components outside range of cg:
   if(istwf_k==1)then

     if(use_symmetry==0)then
       if(iflag==-1)then
!        $omp parallel do private(i1,i2,i3,ipw) shared(cfft,cg,kg_k,npw)
         do ipw=1,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           cg(1,ipw)=cfft(1,i1,i2,i3)*xnorm
           cg(2,ipw)=cfft(2,i1,i2,i3)*xnorm
         end do
       else
!        $omp parallel do private(i1,i2,i3,ipw) SHARED(cfft,cg,kg_k,npw)
         do ipw=1,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

           cg(1,ipw)=cfft(1,i1,i3,i2)*xnorm
           cg(2,ipw)=cfft(2,i1,i3,i2)*xnorm
         end do
       end if
     else
!      $OMP PARALLEL DO PRIVATE(i1,i2,i3,j1,j2,j3,l1,l2,l3,ipw)&
!      $OMP&SHARED(cfft,cg,id1,id2,id3,kg_k,npw,n1,n2,n3,shiftg,symm)
       do ipw=1,npw
         l1=kg_k(1,ipw)+shiftg(1)
         l2=kg_k(2,ipw)+shiftg(2)
         l3=kg_k(3,ipw)+shiftg(3)
         j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
         j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
         j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
         if(j1<0)j1=j1+n1; i1=j1+1
         if(j2<0)j2=j2+n2; i2=j2+1
         if(j3<0)j3=j3+n3; i3=j3+1
         cg(1,ipw)=cfft(1,i1,i2,i3)*xnorm
         cg(2,ipw)=cfft(2,i1,i2,i3)*xnorm
       end do
     end if

   else if(istwf_k>=2)then

     npwmin=1
     if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
!      Extract cg from cfft, in a way that project on a wavefunction with time-reversal symmetry
       cg(1,1)=cfft(1,1,1,1)*xnorm
       cg(2,1)=0.0d0
       npwmin=2
     end if

     if(use_symmetry==0)then
       if(iflag==-1)then
!        $OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipw) &
!        $OMP&SHARED(cfft,cg,i1inver,i2inver,i3inver,kg_k,npw,npwmin)
         do ipw=npwmin,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

!          Construct the coordinates of -k-G
           i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

!          Here the time-reversal symmetry is used to project from cfft
           cg(1,ipw)=(cfft(1,i1,i2,i3) + cfft(1,i1inv,i2inv,i3inv))*0.5d0*xnorm
           cg(2,ipw)=(cfft(2,i1,i2,i3) - cfft(2,i1inv,i2inv,i3inv))*0.5d0*xnorm
         end do
       else
!        $OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipw) &
!        $OMP&SHARED(cfft,cg,i1inver,i2inver,i3inver,kg_k,npw,npwmin)
         do ipw=npwmin,npw
           i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
           i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
           i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

!          Construct the coordinates of -k-G
           i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

!          Here the time-reversal symmetry is used to project from cfft
           cg(1,ipw)=(cfft(1,i1,i3,i2) + cfft(1,i1inv,i3inv,i2inv))*0.5d0*xnorm
           cg(2,ipw)=(cfft(2,i1,i3,i2) - cfft(2,i1inv,i3inv,i2inv))*0.5d0*xnorm
         end do
       end if
     else

       id1=n1/2+2
       id2=n2/2+2
       id3=n3/2+2

!      $OMP PARALLEL DO PRIVATE(i1,i1inv,i2,i2inv,i3,i3inv,ipw) &
!      $OMP&PRIVATE(j1,j2,j3,l1,l2,l3)&
!      $OMP&SHARED(cfft,cg,id1,id2,id3,kg_k)&
!      $OMP&SHARED(i1inver,i2inver,i3inver,npw,npwmin,n1,n2,n3,shiftg,symm)
       do ipw=npwmin,npw

         i1=kg_k(1,ipw); if(i1<0)i1=i1+n1; i1=i1+1
         i2=kg_k(2,ipw); if(i2<0)i2=i2+n2; i2=i2+1
         i3=kg_k(3,ipw); if(i3<0)i3=i3+n3; i3=i3+1

         i1inv=i1inver(i1) ; i2inv=i2inver(i2) ; i3inv=i3inver(i3)

         l1=kg_k(1,ipw)+shiftg(1)
         l2=kg_k(2,ipw)+shiftg(2)
         l3=kg_k(3,ipw)+shiftg(3)
         j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
         j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
         j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
         if(j1<0)j1=j1+n1 ; i1=j1+1
         if(j2<0)j2=j2+n2 ; i2=j2+1
         if(j3<0)j3=j3+n3 ; i3=j3+1

!        Construct the coordinates of -k-G
         l1=i1inv-(i1inv/id1)*n1-1+shiftg(1)
         l2=i2inv-(i2inv/id2)*n2-1+shiftg(2)
         l3=i3inv-(i3inv/id3)*n3-1+shiftg(3)
         j1=symm(1,1)*l1+symm(1,2)*l2+symm(1,3)*l3
         j2=symm(2,1)*l1+symm(2,2)*l2+symm(2,3)*l3
         j3=symm(3,1)*l1+symm(3,2)*l2+symm(3,3)*l3
         if(j1<0)j1=j1+n1 ; i1inv=j1+1
         if(j2<0)j2=j2+n2 ; i2inv=j2+1
         if(j3<0)j3=j3+n3 ; i3inv=j3+1

!        Here the time-reversal symmetry is used to project from cfft
         cg(1,ipw)=(cfft(1,i1,i2,i3) + cfft(1,i1inv,i2inv,i3inv))*0.5d0*xnorm
         cg(2,ipw)=(cfft(2,i1,i2,i3) - cfft(2,i1inv,i2inv,i3inv))*0.5d0*xnorm
       end do
!      $OMP END PARALLEL DO
     end if

   end if

 else
   write(msg,'(a,i0,a)')'  iflag=',iflag,' not acceptable.'
   MSG_PERS_BUG(msg)
 end if

 if(istwf_k>=2)then
   ABI_DEALLOCATE(i1inver)
   ABI_DEALLOCATE(i2inver)
   ABI_DEALLOCATE(i3inver)
 end if

!DEBUG
!if(iflag==-1)then
!write(std_out,*)'sphere : cg(:,:)'
!do ipw=1,npw
!write(std_out,'(4i6,2es16.6)') ipw,kg_k(1:3,ipw),cg(1:2,ipw)
!end do
!end if
!ENDDEBUG

end subroutine sphere
!!***
