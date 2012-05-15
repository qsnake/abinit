!{\src2tex{textfont=tt}}
!!****f* ABINIT/vso_realspace_nonlop
!! NAME
!!   vso_realspace_nonlop
!!
!! FUNCTION
!!
!!  Calculate real space (non local) values of the SO part of the
!!   pseudopotentials, from calls to nonlop, then Fourier transforming
!!   As of 10/2008 this routine is probably useless, as the nonlop
!!   call does not allow one to extract, e.g. L x pauli matrices, to
!!   get the factor in front of the momentum operator (ie. the effective A
!!   field)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (Mver)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf,leave_new,mkffnl,nonlop,ph1d3d,sphereboundary,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine vso_realspace_nonlop(atindx,atindx1,dtfil,dtset,gmet,gprimd,hdr,kg,&
       & mpi_enreg,nattyp,ph1d,position_op,psps,rmet,ucvol,vso_realspace_nl,ylm,ylmgr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vso_realspace_nonlop'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_53_ffts
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments -------------------------------

   type(hdr_type),intent(inout) :: hdr
   type(dataset_type),intent(in) :: dtset
   type(pseudopotential_type),intent(in) :: psps
   type(MPI_type),intent(inout) :: mpi_enreg
   type(datafiles_type),intent(in) :: dtfil

   real(dp),intent(in) :: ucvol

   real(dp),intent(in) :: rmet(3,3)
   real(dp),intent(in) :: gprimd(3,3),gmet(3,3)
   real(dp),intent(in) :: position_op(3,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3))

   integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(dtset%ntypat)
   integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)

   real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)

   real(dp),intent(out) :: vso_realspace_nl(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),&
       & dtset%nspinor,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor)

   real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
   real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)

!Local variables -------------------------

 ! variables for ph3d mkffnl and company
 ! dummy variables for nonlop
 ! real variables for nonlop
   integer :: choice,cplex,cpopt_dummy,dimenl1,dimenl2,dimffnl,fft_option
   integer :: idir_dummy,ider,ikg, ikpt, npw,ia,iatom,signs,only_SO,paw_opt_dummy
   integer :: matblk,ispinor,ispinorp,igp
   integer :: spcur_unit,iost
   integer :: i1,i2,i3,i1p,i2p,i3p,irealsp,irealsp_p
   integer,allocatable :: gbound(:,:),kg_k(:,:)

   real(dp) :: lambda_dummy,arg
   real(dp),allocatable :: dummy_denpot(:,:,:)
   real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),vso_realrecip(:,:,:,:,:)
   real(dp),allocatable :: vectin_ft(:,:),vectin(:,:),svectout_dummy(:,:)
   real(dp),allocatable :: vectout(:,:),vectout_ft(:,:,:,:),sij_dummy(:,:)
   real(dp),allocatable :: enlout_dummy(:),ffnl(:,:,:,:)
   real(dp),allocatable :: dummy_fofgout(:,:),kpg_dummy(:,:)

   character(len=fnlen) :: filnam
   character(len=500) :: message

   type(cprj_type),allocatable :: cprjin_dummy(:)

! *********************************************************************

 if (mpi_enreg%paral_spin==1) then
   write (message,'(2a)')' spin_current: ERROR-',ch10,&
&   'unable to treat parallelization over spinorial components !'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!variables for nonlop
 choice = 1 ! NL energy contribution, not derivatives
 signs = 2 ! get function of G instead of contracted KS matrix element
!only_SO 1 gets the full SO potential  (V_SO L.S) (G,s,G',s')
!only_SO 2 gets a partial SO potential (V_SO   S) (G,s,G',s') then FT wrt G,G'
 only_SO = 2

 cplex=2
 fft_option = 0 ! just do direct fft

 cpopt_dummy = -1
 idir_dummy = 0 ! should not be used
 lambda_dummy = zero
 paw_opt_dummy=0

!allocate stuff for nonlop that does not depend on npw/kpt
 ABI_ALLOCATE(cprjin_dummy,(dtset%natom*((cpopt_dummy+3)/3)))
 ABI_ALLOCATE(sij_dummy,(dimenl1,dtset%ntypat*((paw_opt_dummy+1)/3)))
 ABI_ALLOCATE(enlout_dummy,(1))
 ABI_ALLOCATE(dummy_denpot,(cplex*dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6)))
 ABI_ALLOCATE(gbound,(2*dtset%mgfft+8,2))

!dimensions for ffnl and nonlop
 dimenl1 = psps%dimekb
 dimenl2 = dtset%ntypat
 dimffnl=1
 matblk=dtset%natom

!choose which kpt we will use to get V_SO (closest to Gamma probably best)
 ikg=0
 do ikpt=1,dtset%nkpt
   if ( sum(abs(dtset%kpt(:,ikpt))) < tol10) exit
   ikg=ikg+hdr%npwarr(ikpt)
 end do
 write(std_out,*) 'Found Gamma to be ikpt ', ikpt, dtset%kpt(:,ikpt)
 write(std_out,*) ' ikg = ', ikg

 npw = hdr%npwarr(ikpt)

 ABI_ALLOCATE(kg_k,(3,npw))
 kg_k = kg(:,ikg+1:ikg+npw)

!rebuild phkxred
 ABI_ALLOCATE(phkxred,(2,dtset%natom))
 do ia=1,dtset%natom
   iatom=atindx(ia)
   arg=two_pi*(dtset%kpt(1,ikpt)*hdr%xred(1,ia)&
&   +dtset%kpt(2,ikpt)*hdr%xred(2,ia)&
&   +dtset%kpt(3,ikpt)*hdr%xred(3,ia))
   phkxred(1,iatom)=cos(arg)
   phkxred(2,iatom)=sin(arg)
 end do

!rebuild ph3d
 ABI_ALLOCATE(ph3d,(2,npw,matblk))
 call ph1d3d(1,dtset%natom,kg_k,matblk,dtset%natom,npw,&
& dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),&
& phkxred,ph1d,ph3d)


!rebuild ffnl
 ider=0
 ABI_ALLOCATE(ffnl,(npw,dimffnl,psps%lmnmax,dtset%ntypat))
 ABI_ALLOCATE(kpg_dummy,(npw,0))
 call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
& gmet,gprimd,ider,idir_dummy,psps%indlmn,kg_k,&
& kpg_dummy,dtset%kpt(:,ikpt),psps%lmnmax,&
& psps%lnmax,psps%mpsang,psps%mqgrid_ff,0,&
& npw,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
& psps%usepaw,psps%useylm,ylm,ylmgr)

!get gbound
 call sphereboundary(gbound,dtset%istwfk(ikpt),kg_k,dtset%mgfft,npw)

!allocations for nonlop
 ABI_ALLOCATE(vectin ,(2,dtset%nspinor*npw))
 ABI_ALLOCATE(vectout,(2,dtset%nspinor*npw))

 ABI_ALLOCATE(svectout_dummy,(2,dtset%nspinor*npw*(paw_opt_dummy/3)))
 ABI_ALLOCATE(vectin_ft,(2,npw))
 ABI_ALLOCATE(vectout_ft,(2,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6)))
 ABI_ALLOCATE(vso_realrecip,(2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%nspinor,npw,dtset%nspinor))

!for each spinorial component
 do ispinorp=1,dtset%nspinor
!  for each planewave G',
   do igp=1,npw
!    make wavefunction with only that component
!    probably to be changed: loop over ks states and call nonlop with them
!    eventually premultiplying with r????
!    
!    Aaaaah maybe not: want full spatial
!    dependency and nonlop gives you a projected quantity summed over the G of the
!    KS state
!    
!    This is actually a barbaric way of extracting the so potential 1 GG' pair
!    at a time
     vectin = zero
     vectin(1,(ispinorp-1)*npw+igp) = one

!    and call nonlop -> get <G|V_SO|G'> for all G
!    added flag to not calculate scalar relativistic term, only SO
     call nonlop(atindx1,choice,cpopt_dummy,cprjin_dummy,dimenl1,dimenl2,dimffnl,dimffnl,&
&     psps%ekb,enlout_dummy,ffnl,ffnl,gmet,gprimd,idir_dummy,psps%indlmn,dtset%istwfk(ikpt),&
&     kg_k,kg_k,kpg_dummy,kpg_dummy,dtset%kpt(:,ikpt),dtset%kpt(:,ikpt),&
&     lambda_dummy,psps%lmnmax,matblk,dtset%mgfft,&
&     mpi_enreg,psps%mpsang,psps%mpssoang,dtset%natom,nattyp,dtset%ngfft,0,0,dtset%nloalg,&
&     1,npw,npw,dtset%nspinor,dtset%nspinor,dtset%ntypat,only_SO,paw_opt_dummy,phkxred,&
&     phkxred,ph1d,ph3d,ph3d,signs,sij_dummy,svectout_dummy,&
&     0,ucvol,psps%useylm,vectin,vectout,use_gpu_cuda=dtset%use_gpu_cuda)


!    FT wrt G, one spinorial component of vectout at a time
     do ispinor=1,dtset%nspinor
       vectin_ft = vectout(:,(ispinor-1)*npw+1:(ispinor)*npw)

       call fourwf(cplex,dummy_denpot,vectin_ft,dummy_fofgout,&
&       vectout_ft,gbound,gbound,&
&       hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&       npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&       fft_option,dtset%paral_kgb,0,one,one,use_gpu_cuda=dtset%use_gpu_cuda)

       vso_realrecip(:,:,ispinor,igp,ispinorp)=&
&       reshape(vectout_ft(:,1:dtset%ngfft(1),1:dtset%ngfft(2),1:dtset%ngfft(3)),&
&       (/2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)/))
     end do  ! ispinor
   end do  ! igp

 end do ! ispinorp

 ABI_DEALLOCATE(kpg_dummy)
 ABI_DEALLOCATE(svectout_dummy)

!FT wrt Gprim
 do ispinor=1,dtset%nspinor
   do irealsp=1,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
     do ispinorp=1,dtset%nspinor
       vectin_ft = vso_realrecip(:,irealsp,ispinor,1:npw,ispinorp)

       call fourwf(cplex,dummy_denpot,vectin_ft,dummy_fofgout,&
&       vectout_ft,gbound,gbound,&
&       hdr%istwfk(ikpt),kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw,&
&       npw,dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&       fft_option,dtset%paral_kgb,0,one,one,use_gpu_cuda=dtset%use_gpu_cuda)

       vso_realspace_nl(:,irealsp,ispinor,:,ispinorp) = &
&       reshape(vectout_ft(:,1:dtset%ngfft(1),1:dtset%ngfft(2),1:dtset%ngfft(3)),&
&       (/2,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)/))
     end do
   end do
 end do
 ABI_DEALLOCATE(vso_realrecip)
 ABI_DEALLOCATE(vectout_ft)

!DEBUG check symmetric quality of vso_realspace_nl
!do ispinor=1,dtset%nspinor
!do irealsp=1,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
!do ispinorp=1,dtset%nspinor
!do irealsp_p=1,dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
!write (1666,'(2E16.10)') vso_realspace_nl(:,irealsp,ispinor,irealsp_p,ispinorp) &
!&                      - vso_realspace_nl(:,irealsp_p,ispinorp,irealsp,ispinor)
!end do
!end do
!end do
!end do
!ENDDEBUG

 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(vectin)
 ABI_DEALLOCATE(vectout)
 ABI_DEALLOCATE(ffnl)
 ABI_DEALLOCATE(phkxred)
 ABI_DEALLOCATE(ph3d)

 ABI_DEALLOCATE(cprjin_dummy)
 ABI_DEALLOCATE(sij_dummy)
 ABI_DEALLOCATE(enlout_dummy)
 ABI_DEALLOCATE(dummy_denpot)
 ABI_DEALLOCATE(gbound)

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!output SO potential (non local) for each pair of real space points
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
 filnam=trim(dtfil%filnam_ds(4))//"_VSO_rrp"
 spcur_unit=200
 open (file=filnam,unit=spcur_unit,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(2a)')' spin_current: ERROR- opening file ',trim(filnam)
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!print header
 write (spcur_unit,'(a)') &
& '#  SO operator (nonlocal) as a function of real space point rprime, for fixed r'
 write (spcur_unit,'(a,3(I5,1x))') &
& '#  fft grid is ', dtset%ngfft(1), dtset%ngfft(2),   dtset%ngfft(3)
 write (spcur_unit,'(a)') &
& '#  cart xprime * cart yprime * cart zprime ***   V_SO '
!
!NOTE: have chosen actual dims of grid (n123) instead of fft box, for which n45
!may be different - forced to be odd for FT
!
 i1=1
 i2=1
 i3=1
 write (spcur_unit,'(a,3(E12.5,1x))') &
& '# position of first r point for V_SO(r,rprime): ', &
& position_op(:,i1,i2,i3)
!look at a given spinorial component of V_SO matrix:
 ispinor=1
 ispinorp=2

!do i3=1,dtset%ngfft(3)
!do i2=1,dtset%ngfft(2)
!do i1=1,dtset%ngfft(1)

 irealsp = i1 + (i2-1)*dtset%ngfft(1) + (i3-1)*dtset%ngfft(2)*dtset%ngfft(1)
 do i3p=1,dtset%ngfft(3)
   do i2p=1,dtset%ngfft(2)
     do i1p=1,dtset%ngfft(1)
       irealsp_p = i1p + (i2p-1)*dtset%ngfft(1) + (i3p-1)*dtset%ngfft(2)*dtset%ngfft(1)
!      write (spcur_unit,'(3(E12.5,1x),3x,3(E12.5,1x),3x,2(E20.10,1x))')&
       write (spcur_unit,'(3(E12.5,1x),3x,2(E20.10,1x))')&
&       position_op(:,i1p,i2p,i3p), &
&       vso_realspace_nl(:,irealsp,ispinor,irealsp_p,ispinorp)
     end do
   end do
 end do

!end do
!end do
!end do

 close (spcur_unit)

end subroutine vso_realspace_nonlop
!!***
