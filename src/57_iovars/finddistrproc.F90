!{\src2tex{textfont=tt}}
!!****f* ABINIT/finddistrproc
!! NAME
!! finddistrproc
!!
!! FUNCTION
!! Find a distribution of processor to fill npkpt, npband, npfft and bandpp, knowing nproc
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!  mband=maximum number of bands.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mpi_enreg%keywp(6,250)= transient storage to print nproc,npkpt,npspinor,npband,npfft, and bandpp.
!!  mpi_enreg%trialproc(2)= transient storage for the number of processors and the error flag.
!!  dtset%npkpt= number of processors for parallelisation on k points
!!  dtset%npband=number of processors for parallelisation on bands
!!  dtset%npfft=number of processors for parallelisation on fft grid
!!  dtset%bandpp= internal parameter for lobpcg parallelisation algorithm
!!
!! PARENTS
!!      invars1
!!
!! CHILDREN
!!      initmpi_world,leave_new,sort_int,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine finddistrproc(dtset,mband,mpi_enreg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'finddistrproc'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_28_numeric_noabirule
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: mband
 type(dataset_type),intent(inout) :: dtset
!arrays

!Local variables-------------------------------
!scalars
 integer :: icount,i1,i3,i4,i5,ii,in,na,ncount,ncount1,ncount2,ncount3,ncount4,nleft
 integer :: npband,npfft,npkpt,nproc,nproc1,nprocmin,nspin
 character(len=500) :: message
!arrays
 integer,allocatable :: bandppndiv(:),ndiv(:),distp(:,:),distp1(:,:),distp2(:),iperm(:),work(:),work1(:),work2(:)
 real(dp),allocatable :: weight(:),rwork(:)

!******************************************************************
!DEBUG
!write(std_out,*)' finddistrproc : enter'
!flush(6)
!stop
!ENDDEBUG

!Nothing to do if the input file contains one of the parallelism keywords
 if (((dtset%npkpt/=1.or.dtset%npband/=1.or.dtset%npfft/=1.or. &
& dtset%npspinor/=1.or.dtset%bandpp/=1).and.dtset%paral_kgb>=0).or.&
& (mpi_enreg%nproc==1.and.dtset%paral_kgb==1).or.dtset%paral_kgb==0) return


 nprocmin=2
 if(mpi_enreg%paral_compil==1)then
   if (dtset%paral_kgb==0) return
   if(dtset%paral_kgb <0) then
     nproc=-dtset%paral_kgb;if (dtset%nimage>1) nproc=nproc/dtset%npimage
!    nprocmin=max(2,nproc-100)
   else
     nproc=mpi_enreg%nproc;if (dtset%nimage>1) nproc=nproc/dtset%npimage
     nprocmin=nproc
   end if
 else
   if (dtset%paral_kgb>=0) return
   nproc=-dtset%paral_kgb;if (dtset%nimage>1) nproc=nproc/dtset%npimage
!  nprocmin=max(2,nproc-100)
 end if
 mpi_enreg%trialproc(1)=nproc

!Spin treatment (could be improved)
 nspin=max(dtset%nsppol,dtset%nspinor)

!Compute bandpp and the list if divisors of mband
 ABI_ALLOCATE(ndiv,((mband/2)+1))
 ABI_ALLOCATE(bandppndiv,((mband/2)+1))
 bandppndiv=1
 ndiv=1
 icount=0
 do ii=1,mband/2
   na=mband/ii
   if(na*ii==mband)then
     icount=icount+1
     ndiv(icount)=na
     if(mband/(2*na)*(2*na)==mband)bandppndiv(icount)=2
     if(mband/(4*na)*(4*na)==mband)bandppndiv(icount)=4
   end if
 end do
 ncount=icount+1

!calculate all the possible npband and npfft
 ABI_ALLOCATE(distp,(3,ncount*(nproc-nprocmin+1)))
 distp=0
 icount=0
 do in=nprocmin,nproc
   npband=1;npfft=1;npkpt=1
!  if(2*dtset%nkpt*nspin>=in) then
   if(dtset%nkpt*nspin>=in) then
     npkpt=in
   else
     npkpt=dtset%nkpt*nspin
     if((dtset%paral_kgb==1).and.(in/npkpt*npkpt/=in)) npkpt=1
     nleft=in/npkpt
     i1=0
     do ii=1,ncount
!      if((mband/(bandpp*ndiv(ii))*(bandpp*ndiv(ii)))/=mband) cycle
       if(ndiv(ii)>nleft)cycle
       npband=ndiv(ii)
       if (npband==1) cycle
       i1=1
       if(dtset%use_gpu_cuda==1) then
         npfft=1
       else
         npfft=nleft/npband
       end if
       if (npfft>npband) cycle
       nproc1=npkpt*npband*npfft
       if (nproc1==0) cycle
       if((dtset%paral_kgb==1).and.(npband*npfft*npkpt/=nproc)) cycle
       icount=icount+1
       distp(1,icount)=npband;distp(2,icount)=npfft;distp(3,icount)=bandppndiv(ii)
     end do
   end if
 end do
 ncount1=icount
 if(icount==0) then
   distp(1,1)=npband;distp(2,1)=npfft;distp(3,1)=1
   ncount1=1
 end if

!Sort the possible npband by ascending order
 ABI_ALLOCATE(distp1,(3,ncount1))
 ABI_ALLOCATE(iperm,(ncount1))
 ABI_ALLOCATE(work,(ncount1))
 ABI_ALLOCATE(work1,(ncount1))
 do ii=1,ncount1
!  write(std_out,*)'coucou',distp(1,ii),distp(2,ii),distp(3,ii)
   iperm(ii)=ii
 end do
 call sort_int(ncount1,distp(1,1:ncount1),iperm)
 work(1:ncount1)=distp(2,1:ncount1)
 work1(1:ncount1)=distp(3,1:ncount1)
 do ii=1,ncount1
   distp(2,ii)=work(iperm(ii))
   distp(3,ii)=work1(iperm(ii))
 end do

!Eliminate the same distributions
 if (ncount1>=2) then
   i3=1;i4=1;work1(1)=0
   do icount=2,ncount1
     if ((distp(1,icount)==distp(1,icount-1))) then
       i4=i4+1
     else
       i3=i3+1
       work1(i3)=i4
       i4=1
     end if
   end do
   if ((distp(1,ncount1)/=distp(1,ncount1-1))) then
     i3=i3+1
     work1(i3)=1
   else
     i3=i3+1
     work1(i3)=i4
   end if
   ncount2=i3
   icount=0
   do ii=2,ncount2
     do i1=1,work1(ii)
       iperm(i1)=i1
     end do
     call sort_int(work1(ii),distp(2,icount+1:icount+work1(ii)),iperm(1:work1(ii)))
     work(1:work1(ii))=distp(3,1+icount:icount+work1(ii))
     do i5=1,work1(ii)
       distp(3,1+icount:icount+work1(ii))=work(iperm(i5))
     end do
     icount=icount+work1(ii)
   end do
   distp1(:,1)=distp(:,1)
   i3=1
   do icount=2,ncount1
     if ((distp(1,icount)/=distp1(1,i3)).or.(distp(2,icount)/=distp1(2,i3))) then
       i3=i3+1
       distp1(:,i3)=distp(:,icount)
     end if
   end do
   ncount4=i3
 else
   ncount4=1
   distp1(:,1)=distp(:,1)
 end if

 ABI_ALLOCATE(distp2,(ncount4))
 do icount=1,ncount4
   distp2(icount)=npkpt*distp1(1,icount)*distp1(2,icount)
 end do
 do i1=1,ncount4
   iperm(i1)=i1
 end do
 call sort_int(ncount4,distp2,iperm(1:ncount4))
 ABI_ALLOCATE(work2,(ncount4))
 work(1:ncount4)=distp1(1,1:ncount4)
 work1(1:ncount4)=distp1(2,1:ncount4)
 work2(1:ncount4)=distp1(3,1:ncount4)
 do i1=1,ncount4
   distp1(1,i1)=work(iperm(i1))
   distp1(2,i1)=work1(iperm(i1))
   distp1(3,i1)=work2(iperm(i1))
 end do
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(work1)
 ABI_DEALLOCATE(work2)

 ABI_ALLOCATE(weight,(ncount4))
 do icount=1,ncount4
   weight(icount)=distp1(1,icount)/distp1(2,icount)/4.d0
   if ((distp1(1,icount)<=50).and.(distp1(2,icount)==1)) weight(icount)=1
!  write(std_out,*),distp1(1,icount),distp1(2,icount),distp1(3,icount)
 end do

!Store and print the results
 if(ncount4>250) then
   write(message,'(a,a,a,a)' )ch10,&
&   'WARNING in finddistrproc: more than 250 possible choices for nproc',ch10,&
&   'Only the first 250 ones are printed'
 end if
 if(dtset%paral_kgb<0)then
   dtset%npkpt=mpi_enreg%nproc
   ncount3=min(ncount4,250)
   mpi_enreg%trialproc(2)=1
   do icount=ncount3,1,-1
     mpi_enreg%keywp(1,ncount3-icount+1)=npkpt*distp1(1,icount)*distp1(2,icount)
     if (dtset%nspinor==2) then
       if (mod(npkpt,2)==0) then
         mpi_enreg%keywp(2,ncount3-icount+1)=npkpt/2
         mpi_enreg%keywp(3,ncount3-icount+1)=2
       else
         mpi_enreg%keywp(1,ncount3-icount+1)=0
       end if
     else
       mpi_enreg%keywp(2,ncount3-icount+1)=npkpt
       mpi_enreg%keywp(3,ncount3-icount+1)=1
     end if
     mpi_enreg%keywp(4,ncount3-icount+1)=distp1(1,icount)
     mpi_enreg%keywp(5,ncount3-icount+1)=distp1(2,icount)
     mpi_enreg%keywp(6,ncount3-icount+1)=distp1(3,icount)
   end do
   mpi_enreg%trialproc(3)=mpi_enreg%nproc
 elseif(dtset%paral_kgb==1)then
   if (ncount4>1) then
     ABI_ALLOCATE(rwork,(ncount4))
     do ii=1,ncount4
       rwork(ii)=abs(weight(ii)-1.d0)
     end do
     do ii=1,ncount4
       if (abs(rwork(ii)-minval(rwork))<1d-8) i1=ii
     end do
     ABI_DEALLOCATE(rwork)
   else
     i1=1
   end if
   if (distp1(1,i1)*distp1(2,i1)/=distp1(1,ncount4)*distp1(2,ncount4)) i1=ncount4
   if (dtset%nspinor==2) then
     dtset%npkpt=npkpt/2
     dtset%npspinor=2-mod(npkpt,2)
   else
     dtset%npkpt=npkpt
     dtset%npspinor=1
   end if
   dtset%npband=distp1(1,i1)
   dtset%npfft=distp1(2,i1)
   dtset%bandpp=distp1(3,i1)
!  if(distp1(1,i1)*distp1(2,i1)*npkpt/=mpi_enreg%nproc) dtset%npkpt=mpi_enreg%nproc/distp1(1,i1)/distp1(2,i1)+1
   mpi_enreg%keywp(1,1)=npkpt*distp1(1,i1)*distp1(2,i1)
   mpi_enreg%keywp(2,1)=npkpt
   mpi_enreg%keywp(4,1)=distp1(1,i1)
   mpi_enreg%keywp(5,1)=distp1(2,i1)
   mpi_enreg%keywp(6,1)=distp1(3,i1)

   if(mpi_enreg%keywp(1,1)>mpi_enreg%nproc)then
     write(message,'(a,a)' )ch10,&
&     'ERROR in finddistrproc: nproc must be equal to npkpt*npband*npfft'
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   else
     if (mpi_enreg%trialproc(3)==0)then
       call initmpi_world(mpi_enreg,mpi_enreg%keywp(1,1)*dtset%npimage)
     else
       if (mpi_enreg%keywp(1,1)*dtset%npimage/=mpi_enreg%nproc) then
         write(message,'(a,a)' )ch10,&
&         'ERROR in finddistrproc: nproc must not be changed between two datasets'
         call wrtout(std_out,  message,'COLL')
         call leave_new('COLL')
       end if
     end if
   end if
   mpi_enreg%trialproc(3)=1
 end if

 ABI_DEALLOCATE(bandppndiv)
 ABI_DEALLOCATE(ndiv)
 ABI_DEALLOCATE(distp)
 ABI_DEALLOCATE(distp1)
 ABI_DEALLOCATE(distp2)
 ABI_DEALLOCATE(iperm)
 ABI_DEALLOCATE(weight)

!write(std_out,*)' finddistrproc : out'
!flush(6)
!DEBUG
!stop
!ENDDEBUG

end subroutine finddistrproc


!!***
