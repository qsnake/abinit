!{\src2tex{textfont=tt}}
!!****f* ABINIT/wfread
!! NAME
!! wfread
!!
!! FUNCTION
!! returns wave function in real space and energy eigen values
!! (this file is just a shortened version of wffile.F90)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (JB)
!! this file is distributed under the terms of the
!! gnu general public license, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! Needs an unformatted wave function from abinit.
!! exchn2n3d=if 1, n2 and n3 are exchanged
!! csppol = spin polarization
!! cbandpick = bandindex for the wf
!! ckpt = kpoint index for the wf
!! ecut = effective ecut (ecut*dilatmx**2)
!! headform=format of the wf file
!! istwfk = input variable indicating the storage option of each k-point
!! kpt = array of k-points coordinates
!! nband = input variable nband(nkpt*nsppol) !! MC 090902: the definition below in not consistent with what we find in defs_datatypes.F90
!! nbands = size of e_kpt
!! nkpt = number of k-points
!! npwarr = array holding npw for each k point
!! nr1,nr2,nr3 = grid size (nr1 x nr2 x nr3 = filrho dimension)
!! nspinor = number of spinorial components of the wavefunctions
!! nsppol = number of spin polarization
!! paral_kgb = parallization option, it is set to 0 in the parent subroutine
!! rprim = orientation of the unit cell axes
!!
!! OUTPUT
!!  cwave0 = wave function (corresponding to cbandpick ckpt csppol) in real space
!!  e_kpt = all the energy eigen values at the particular kpoint
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      localorb_S
!!
!! CHILDREN
!!      fourwf,getkpgnorm,hdr_skip,kpgio,metric,rwwf,sphereboundary
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wfread(cwave0,e_kpt,exchn2n3d,csppol,cbandpick,ckpt,&
     & ecut,headform,istwfk,kpt,nband,nbands,nkpt,npwarr,&
     & nr1,nr2,nr3,nspinor,nsppol,paral_kgb,rprimd)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfread'
 use interfaces_42_geometry
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cbandpick,ckpt,csppol,exchn2n3d,headform,nbands,nkpt,nr1
 integer,intent(in) :: nr2,nr3,nspinor,nsppol,paral_kgb
 real(dp),intent(in) :: ecut
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt),npwarr(nkpt)
 real(dp),intent(in) :: kpt(3,nkpt),rprimd(3,3)
 real(dp),intent(out) :: e_kpt(nbands)
!no_abirules
 complex(dp),intent(out) :: cwave0(nr1,nr2,nr3)

!Local variables-------------------------------
       character(*), parameter :: inputfile='cut.in'
!scalars
 integer,save :: tim_fourwf=0,tim_rwwf=0
 integer :: cband,cgshift,cplex,cspinor,formeig,i
 integer :: ierr
 integer :: ii1,ikpt,ioffkg,iout
 integer :: isppol,j,k,l,m1
! integer :: ix,iy,iz
 integer :: mband,mcg,mgfft,mkmem,mpw,n,n4,n5,n6,nband_disk
 integer :: npw_k
 integer :: oldcband,oldckpt,oldcspinor,oldcsppol,option
 integer :: unkg
 real(dp) :: energy
 real(dp) :: re1,re2
 real(dp) :: rnorm,tpi,ucvol
 real(dp) :: weight
 character(len=4) :: mode_paral
 character(len=fnlen) :: kgnam
 type(mpi_type) :: mpi_enreg
 type(wffile_type) :: wff
!arrays
 integer :: ngfft(18)
 integer,allocatable :: gbound(:,:),kg(:,:),kg_dum(:,:),kg_k(:,:)
 integer,allocatable :: npwarr1(:),npwtot1(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: cg(:,:),cgcband(:,:),denpot(:,:,:),eigen1(:)
 real(dp),allocatable :: fofgout(:,:),fofr(:,:,:,:),kpgnorm(:)
 real(dp),allocatable :: occ1(:)

! *************************************************************************

 cspinor = 1 ! only scalar wave functions.

!begin executable section
 mpi_enreg%paralbd=0

 formeig=0
 oldckpt=0
 oldcband=0
 oldcsppol=0
 oldcspinor=0

 iout=-1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

 ABI_ALLOCATE(kg_dum,(3,0))

!#############################################################################

 tpi = two_pi

 mband=maxval(nband(1:nkpt))
 mpw=maxval(npwarr)
 mcg=mpw*nspinor*mband

 ABI_ALLOCATE(cg,(2,mcg))
 ABI_ALLOCATE(eigen1,((2*mband)**formeig*mband))
 ABI_ALLOCATE(occ1,(mband))

!==========================================================================
!necessary procedures for the fft subroutine:
!==========================================================================
 mpi_enreg%paralbd=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%me=0
 mpi_enreg%me_fft=0
 mpi_enreg%paral_fft=0
 mpi_enreg%paral_compil_kpt=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%fft_option_lob=1
 mpi_enreg%mode_para="n"
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0
 mpi_enreg%paral_spin = 0

 ngfft(1)=nr1
 ngfft(2)=nr2
 ngfft(3)=nr3

 if (mod(nr1,2)==0)then
   ngfft(4)=nr1+1
 else
   ngfft(4)=nr1
 end if
 if(mod(nr2,2)==0)then
   ngfft(5)=nr2+1
 else
   ngfft(5)=nr2
 end if
 ngfft(6)=nr3
 ngfft(7)=111
 ngfft(8)=256

 mode_paral='pers'
 mkmem=nkpt
 mgfft=maxval(ngfft(1:3))

 ABI_ALLOCATE(npwarr1,(nkpt))
 ABI_ALLOCATE(kg,(3,mpw*mkmem))
 ABI_ALLOCATE(npwtot1,(nkpt))
!create positions index for pw
 call kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kgnam,kpt,mkmem,nband,nkpt,&
 mode_paral,mpi_enreg,mpw,npwarr1,npwtot1,nsppol,unkg)

!additional allocation:
 n4=       ngfft(4)
 n5=       ngfft(5)
 n6=       ngfft(6)

 cwave0 = 0.0
 e_kpt = 0.0
!-------------------------------------------------------------------
!reading wave function from "_wfk" file :
!-------------------------------------------------------------------
 cband = cbandpick
 cg = 0.0
 eigen1 = 0.0
 occ1 = 0.0
 rewind(19)

 wff%unwff=19
 wff%accesswff=IO_MODE_FORTRAN
 call hdr_skip(wff,ierr)

!call hdr_skip(19)

 do isppol=1,csppol
!  write(std_out,*)'ckpt',ckpt
   do ikpt=1,nkpt

     if(isppol==csppol .and. ikpt==ckpt)then
       option=1
!      write(std_out,*)'option',option,cband
     else
       option=-1
     end if
     call rwwf(cg,eigen1,formeig,headform,0,ikpt,isppol,kg_dum,&
&     mband,mcg,mpi_enreg,nband(ikpt),nband_disk,&
&     npwarr(ikpt),nspinor,occ1,option,0,tim_rwwf,wff)

     if(isppol==csppol.and.ikpt==ckpt)then
       re1 = 0.0
       re2 = 0.0
       do ii1 = 1,mband
         energy=eigen1(ii1)
         e_kpt(ii1) = energy
       end do
     end if

     if(option==1)exit       ! when the target wf has been read,
!    exit the wf file reading
   end do
   if(option==1)exit
 end do

 ioffkg=0
 do ikpt=1,ckpt-1
   ioffkg=ioffkg+npwarr1(ikpt)
 end do
 npw_k=npwarr(ckpt)

 ABI_ALLOCATE(gbound,(2*mgfft+8,2))
 ABI_ALLOCATE(kg_k,(3,npw_k))
 ABI_ALLOCATE(kpgnorm,(npw_k))

 kg_k(:,1:npw_k)=kg(:,1+ioffkg:npw_k+ioffkg)
 call getkpgnorm(gprimd,kpt(:,ckpt),kg_k,kpgnorm,npw_k)
 call sphereboundary(gbound,istwfk(ckpt),kg_k,mgfft,npw_k)
 n4=ngfft(4)
 n5=ngfft(5)
 n6=ngfft(6)
!cplex=0
 cplex=1
 cgshift=(cband-1)*npw_k*nspinor + (cspinor-1)*npw_k

 ABI_ALLOCATE(cgcband,(2,npw_k))
 ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
 ABI_ALLOCATE(fofgout,(2,npw_k))
 ABI_ALLOCATE(fofr,(2,n4,n5,n6))

 cgcband(:,1:npw_k)=cg(:,cgshift+1:cgshift+npw_k)

 call fourwf(cplex,denpot,cgcband,fofgout,fofr,&
& gbound,gbound,&
& istwfk(ckpt),kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,&
& npw_k,n4,n5,n6,0,paral_kgb,tim_fourwf,weight,weight)

!write(std_out,*)n4,n5,n6
!write(std_out,*)nr1,nr2,nr3
!rnorm = 0.0
!do iz=1,nr3
!do iy=1,nr2
!do ix=1,nr1
!cwave1(ix,iy,iz)=cmplx(fofr(1,ix,iy,iz),fofr(2,ix,iy,iz))
!rnorm = rnorm + conjg(cwave1(ix,iy,iz))*cwave1(ix,iy,iz)
!end do
!end do
!end do
!Swaping :
 rnorm = 0.0
 do k=1,nr3
   if(mod(real(nr3),2.0).eq.0.0)then
     if(k < ((nr3/2)+1))n=k+(nr3/2)
     if(k > (nr3/2))n=k-(nr3/2)
   else
     if(k < (nr3-1)/2+1)n=k+(nr3+1)/2
     if(k > (nr3-1)/2)n=k-(nr3-1)/2
   end if
   do j=1,nr2
     if(mod(real(nr2),2.0).eq.0.0)then
       if(j < ((nr2/2)+1))m1=j+(nr2/2)
       if(j > (nr2/2))m1=j-(nr2/2)
     else
       if(j < (nr2-1)/2+1)m1=j+(nr2+1)/2
       if(j > (nr2-1)/2)m1=j-(nr2-1)/2
     end if
     do i=1,nr1
       if(mod(real(nr1),2.0).eq.0.0)then
         if(i < ((nr1/2)+1))l=i+(nr1/2)
         if(i > (nr1/2))l=i-(nr1/2)
       else
         if(i < (nr1-1)/2+1)l=i+(nr1+1)/2
         if(i > (nr1-1)/2)l=i-(nr1-1)/2
       end if
!      cwave0(i,j,k)=cwave1(l,m1,n)
       cwave0(i,j,k)=cmplx(fofr(1,l,m1,n),fofr(2,l,m1,n))
       rnorm = rnorm + conjg(cwave0(i,j,k))*cwave0(i,j,k)
     end do
   end do
 end do

 ABI_DEALLOCATE(cgcband)
 ABI_DEALLOCATE(denpot)
 ABI_DEALLOCATE(fofgout)
 ABI_DEALLOCATE(fofr)

 ABI_DEALLOCATE(gbound)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kpgnorm)

 cwave0 = cwave0/sqrt(rnorm)
 return

end subroutine wfread
!!***
