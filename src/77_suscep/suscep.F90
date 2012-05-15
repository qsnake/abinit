!{\src2tex{textfont=tt}}
!!****f* ABINIT/suscep
!! NAME
!! suscep
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of the polarisability
!! within the random phase approximation (RPA)
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (XG,GMR,MF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mband =maximum number of bands
!!  mpi_enreg=informations about MPI parallelization
!!  mpw   =maximum number of planewaves in basis sphere (large number)
!!  natom =number of atoms in unit cell
!!  nkpt  =number of k points
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=number of channels for spin-polarization (1 or 2)
!!  nsym=number of symmetry elements in space group
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (no direct output : results written)
!!
!! SIDE EFFECTS
!!  mkmem =maximum number of k points which can fit in core memory
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      distrb2,getfreqsus,getmpw,getng,hdr_clean,inwffil3,ioarr,kpgio,metric
!!      mkrdim,newocc,setsym,sphereboundary,status,suscep_dyn,suscep_stat,timab
!!      wffclose,wrtout,xcacfd,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine suscep(dtfil,dtset,iexit,&
& mband,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,&
& nspden,nspinor,nsppol,nsym,occ,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_wffile

 use m_header,          only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'suscep'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_77_suscep, except_this_one => suscep
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iexit,mband,mpw,natom,nfft,nkpt,nspden,nsppol,nsym
 integer,intent(inout) :: mkmem,nspinor
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=20
 integer :: accessfil,dielop,fformr,ierr,ifreq,ii,ipw1,ipw2,isp1,isp2
 integer :: lmax_diel,master,mcg,mcprj,me,mgfftdiel,my_nspinor,nband_mx
 integer :: neglect_pawhat,nfftdiel,nfreqsus
 integer :: npwdiel,nspinor_,prtvol,rdwr,rdwrpaw,susopt,usetimerev
 real(dp) :: diecut,ecutsus,entropy,etotal,fermie,ucvol
 character(len=500) :: message
 character(len=fnlen) :: kgnam
 type(hdr_type) :: hdr
 type(pawang_type) :: pawang
 type(wffile_type) :: wff1
 type(wvl_internal_type) :: wvl
!arrays
 integer :: ngfftdiel(18),npwarr_diel(1),npwtot_diel(1)
 integer,allocatable :: atindx1_dum(:),dimcprj(:),gbound_diel(:,:)
 integer,allocatable :: indsym(:,:,:),irrzondiel(:,:,:),kg(:,:),kg_diel(:,:)
 integer,allocatable :: nband_dum(:),npwarr(:),npwtot(:),symrec(:,:,:)
 real(dp) :: dielar(7),gmet(3,3),gprimd(3,3),kpt_diel(3),rmet(3,3),rprimd(3,3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cg(:,:),dielinv(:,:,:,:,:),doccde(:),eigen(:),freq(:)
 real(dp),allocatable :: ph1ddiel(:,:),phnonsdiel(:,:,:),rhor(:,:)
 real(dp),allocatable :: susd_non_dyn(:,:,:,:),susmat(:,:,:,:,:)
 real(dp),allocatable :: susmat_dyn(:,:,:,:,:,:),wght_freq(:),ylmdiel(:,:)
 type(cprj_type),allocatable :: cprj(:,:)
 type(pawrhoij_type),allocatable :: rhoij_dum(:)
 type(pawtab_type),allocatable :: pawtab(:)

!***********************************************************************

 call timab(84,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'enter         ')

 mpi_enreg%paralbd=0
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%paral_fft=0
 mpi_enreg%nproc_spin=1
 mpi_enreg%paral_spin=0
 mpi_enreg%paral_level=2

!
!If dtset%accesswff == 2 set all array outputs to netcdf format
!
 accessfil = 0
 if (dtset%accesswff == IO_MODE_NETCDF) then
   accessfil = 1
 end if
 if (dtset%accesswff == IO_MODE_ETSF) then
   accessfil = 3
 end if

 master=0
!Init me
 call xme_init(mpi_enreg,me)

 if(mpi_enreg%paral_compil_kpt==1)then
   ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt,mband,nsppol))
   call distrb2(mband, dtset%nband, nkpt, nsppol, mpi_enreg)
 end if

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' suscep : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

!Loop input variables
 nfreqsus=dtset%nfreqsus

 dielar(1)=dtset%diecut
 dielar(2)=dtset%dielng
 dielar(3)=dtset%diemac
 dielar(4)=dtset%diemix
 dielar(5)=dtset%diegap
 dielar(6)=dtset%dielam
 dielar(7)=dtset%diemixmag

!Impose mkmem=0 to read cg from disk
 mkmem=0

 call status(0,dtfil%filstat,iexit,level,'call inwffil3 ')

!Compute different geometric tensor, as well as ucvol, from rprimd
 call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Get diecut, and the fft grid to be used for the susceptibility computation
 diecut=abs(dielar(1))
 if( dielar(1)<0.0_dp )then
   ecutsus= dtset%ecut
 else
   ecutsus= ( sqrt( dtset%ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
 end if

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 ngfftdiel(1:3)=0 ; ngfftdiel(7)=101 ; ngfftdiel(8:18)=dtset%ngfft(8:18)
 call getng(dtset%boxcutmin,ecutsus,gmet,mpi_enreg%me_fft,mgfftdiel,nfftdiel,ngfftdiel,&
& mpi_enreg%nproc_fft,nsym,mpi_enreg%fft_option_lob,mpi_enreg%paral_fft,dtset%symrel,&
& use_gpu_cuda=dtset%use_gpu_cuda)

!Compute the size of the dielectric matrix
 kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
 call getmpw(diecut,dtset%exchn2n3d,gmet,(/1/),kpt_diel,&
& mpi_enreg,npwdiel,1)

!Now, performs allocation
 mcg=mpw*my_nspinor*mband*mkmem*nsppol
 ABI_ALLOCATE(cg,(2,mcg))
 ABI_ALLOCATE(eigen,(mband*nkpt*nsppol))
 ABI_ALLOCATE(kg,(3,mpw*mkmem))
 ABI_ALLOCATE(kg_diel,(3,npwdiel))
 ABI_ALLOCATE(npwarr,(nkpt))
 ABI_ALLOCATE(npwtot,(nkpt))
 ABI_ALLOCATE(gbound_diel,(2*mgfftdiel+8,2))
 ABI_ALLOCATE(irrzondiel,(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4)))
 ABI_ALLOCATE(phnonsdiel,(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4)))
 ABI_ALLOCATE(nband_dum,(nsppol))
 mcprj=0

!Then, initialize and compute the values of different arrays
 call kpgio(dtset%ecut,dtset%exchn2n3d,gmet,dtset%istwfk,kg,dtfil%fnametmp_kgs,dtset%kptns,mkmem,&
& dtset%nband,nkpt,'PERS',mpi_enreg,mpw,npwarr,npwtot,nsppol,dtfil%unkg)
!This kpgio call for going from the suscep FFT grid to the diel sphere
!Note : kgnam is dummy, npwarr_diel is dummy, npwtot_diel is dummy, nband_dum is dummy
 nband_dum(:) = 1
 call kpgio(diecut,dtset%exchn2n3d,gmet,(/1/),kg_diel,kgnam,&
& kpt_diel,1,nband_dum,1,'COLL',mpi_enreg,npwdiel,npwarr_diel,npwtot_diel,&
& nsppol,tmp_unit)

 call sphereboundary(gbound_diel,1,kg_diel,mgfftdiel,npwdiel)

 ABI_ALLOCATE(indsym,(4,nsym,natom))
 ABI_ALLOCATE(symrec,(3,3,nsym))
 if (nsym>1) then
   call setsym(indsym,irrzondiel,dtset%iscf,natom,&
&   nfftdiel,ngfftdiel,nspden,nsppol,nsym,phnonsdiel,&
&   dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)
 end if

!Read eigenvalues from the wavefunction file
!Also, initialize wff1 and hdr
 eigen(:)=0.0_dp
!mpi_enreg%paralbd=0
 call inwffil3(dtset,eigen,hdr,dtset%istwfk,mband,mpi_enreg,dtset%nband,&
& nkpt,npwarr,nsppol,prtvol,wff1,dtfil%unwff1,dtfil%fnamewffk)

!Compute new occupation numbers if needed
 ABI_ALLOCATE(doccde,(mband*nkpt*nsppol))
 if(dtset%occopt>=3.and.dtset%occopt<=8) then
   call status(0,dtfil%filstat,iexit,level,'call newocc   ')
   call newocc(doccde,eigen,entropy,fermie,dtset%fixmom,mband,dtset%nband,&
&   dtset%nelect,nkpt,dtset%nspinor,nsppol,occ,dtset%occopt,prtvol,dtset%stmbias,&
&   dtset%tphysel,dtset%tsmear,dtset%wtk)
 end if

 dielop=2 ! Immediate computation of dielectric matrix
 if(nfreqsus==0) then

!  Perform allocations
   ABI_ALLOCATE(dielinv,(2,npwdiel,nspden,npwdiel,nspden))
   ABI_ALLOCATE(susmat,(2,npwdiel,nspden,npwdiel,nspden))
   susmat(:,:,:,:,:)=0._dp

!  Compute the static susceptibility matrix
   lmax_diel=0
   ABI_ALLOCATE(atindx1_dum,(dtset%natom))
   usetimerev=1;neglect_pawhat=0
   call suscep_stat(atindx1_dum,cg,cprj,&
&   dielar,dimcprj,doccde,eigen,gbound_diel,gprimd,&
&   irrzondiel,dtset%istwfk,kg,kg_diel,lmax_diel,&
&   mband,mcg,mcprj,mgfftdiel,mkmem,mpi_enreg,mpw,natom,&
&   dtset%nband,neglect_pawhat,nfftdiel,ngfftdiel,&
&   nkpt,npwarr,npwdiel,nspden,nspinor,nsppol,nsym,dtset%ntypat,&
&   occ,dtset%occopt,pawang,pawtab,phnonsdiel,ph1ddiel,rprimd,&
&   susmat,dtset%symafm,dtset%symrel,dtset%tnons,dtset%typat,ucvol,&
&   dtfil%unkg,0,dtfil%unpaw,0,usetimerev,wff1,dtset%wtk,ylmdiel)
   ABI_DEALLOCATE(atindx1_dum)

!  Print the susceptibility matrix
   do isp1=1,nspden
     do isp2=1,nspden
       write(std_out,'(5x,a,2i2)') 'Susceptibility matrix for spins=',isp1,isp2
       write(std_out,'(9x,a,13x,a,10x,a,10x,a)') "g","g'","real","imag"
       do ipw1=1,10
         do ipw2=ipw1,10
           write(std_out,'(2x,3i4,2x,3i4,2x,f12.8,2x,f12.8)') &
&           kg_diel(1:3,ipw1),kg_diel(1:3,ipw2),&
&           susmat(1,ipw1,isp1,ipw2,isp2),susmat(2,ipw1,isp1,ipw2,isp2)
         end do
       end do
     end do
   end do

!  Perform deallocations
   ABI_DEALLOCATE(dielinv)
   ABI_DEALLOCATE(susmat)
   ABI_DEALLOCATE(nband_dum)

 else if(nfreqsus > 0) then

!  Perform allocations
   ABI_ALLOCATE(freq,(nfreqsus))
   ABI_ALLOCATE(rhor,(nfft,nspden))
   ABI_ALLOCATE(wght_freq,(nfreqsus))

!  Perform initializations
   freq(:)=0.0_dp
   rhor(:,:)=0._dp
   wght_freq(:)=0.0_dp

!  Read in the density, needed for ALDA kernel
   call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
!  Read rho(r) from a disk file
!  Unit numbers and file name for the _KGS file
   rdwr=1;rdwrpaw=0
!  Note : etotal is read here, and might serve in the tddft routine.
   fformr=52
   call ioarr(accessfil,rhor, dtset, etotal,fformr,dtfil%fildensin,hdr, mpi_enreg, &
&   nfft,rhoij_dum,rdwr,rdwrpaw,wvl)
   call status(0,dtfil%filstat,iexit,level,'call fourdp   ')

!  DEBUG
!  leave in for MF to check the density
!  dummy=0._dp
!  do ii=1,nfft
!  dummy=dummy+rhor(ii,1)
!  end do
!  write(std_out,*) '%suscep: nfft=',nfft
!  write(std_out,*) '%suscep: dummy=',dummy*ucvol/dble(nfft)
!  write(std_out,*) '%suscep: ucvol=',ucvol
!  call flush(6)
!  Compute up+down rho(G) by fft
!  allocate(work(nfft))
!  work(:)=rhor(:,1)
!  call fourdp(1,rhog,work,-1,mpi_enreg,nfft,ngfft,0)
!  deallocate(work)
!  ENDDEBUG

!  Create frequency grid and weights, 2 stands for preassigned grid
   call getfreqsus(freq,wght_freq,nfreqsus,dtset%optfreqsus,dtset%freqsuslo,dtset%freqsusin)

!  DEBUG
!  write(std_out,*)' suscep : after getfreqsus '
!  call flush(6)
!  ENDDEBUG

   nspinor_=dtset%nspinor
   call xcacfd(dielar,dtfil,dtset,eigen,freq,gbound_diel,gmet,&
&   gprimd,irrzondiel,kg,kg_diel,mband,mgfftdiel,mkmem,&
&   mpi_enreg,mpw,nfft,nfftdiel,nfreqsus,dtset%ngfft,&
&   ngfftdiel,nkpt,npwarr,npwdiel,nspden,nspinor_,nsppol,&
&   nsym,occ,phnonsdiel,rhor,rprimd,ucvol,wff1,wght_freq)

!  Perform deallocations
   ABI_DEALLOCATE(freq)
   ABI_DEALLOCATE(rhor)
   ABI_DEALLOCATE(wght_freq)

 else

!  Leave intact for testing purposes
   nfreqsus=abs(nfreqsus)

!  Perform allocations
   ABI_ALLOCATE(freq,(nfreqsus))
   ABI_ALLOCATE(susd_non_dyn,(2,npwdiel,nspden,nfreqsus))
   ABI_ALLOCATE(susmat_dyn,(2,npwdiel,nspden,npwdiel,nspden,nfreqsus))

!  Perform initializations
   freq(:)=0.0_dp
   susmat_dyn(:,:,:,:,:,:)=0.0_dp

!  Create a linear frequency grid
   freq(1)=dtset%freqsuslo
   do ifreq=2,nfreqsus
     freq(ifreq)=freq(ifreq-1)+dtset%freqsusin
   end do

!  DEBUG
!  write(std_out,*) '%suscep: nband_mx=', nband_mx
!  ENDDEBUG

!  Compute the dynamical susceptibility matrices
   nspinor_=dtset%nspinor
   call suscep_dyn(dielar,dtset,&
&   eigen,freq,gbound_diel,gprimd,&
&   irrzondiel,dtset%istwfk,kg,kg_diel,&
&   mband,mgfftdiel,mkmem,mpi_enreg,mpw,dtset%nband,nband_mx,nfftdiel,nfreqsus,&
&   ngfftdiel,nkpt,npwarr,npwdiel,nspden,nspinor_,nsppol,nsym,&
&   occ,dtset%occopt,phnonsdiel,rprimd,&
&   susopt,susd_non_dyn,susmat_dyn,dtset%symafm,&
&   dtset%symrel,dtset%tnons,ucvol,dtfil%unkg,wff1,dtset%wtk)

!  Print the dynamical susceptibility matrices
   do ifreq=1,nfreqsus
     write(std_out,'(/,2x,a,f12.8,a)') '---Susceptibility matrices for frequency=',freq(ifreq),'i'
     do isp1=1,nspden
       do isp2=1,nspden
         write(std_out,'(5x,a,2i2)') 'Susceptibility matrix for spins=',isp1,isp2
         write(std_out,'(9x,a,13x,a,10x,a,10x,a)') "g","g'","real","imag"
         do ipw1=1,10
           do ipw2=ipw1,10
             write(std_out,'(2x,3i4,2x,3i4,2x,f12.8,2x,f12.8)') &
&             kg_diel(1:3,ipw1),kg_diel(1:3,ipw2),&
&             susmat_dyn(1,ipw1,isp1,ipw2,isp2,ifreq),susmat_dyn(2,ipw1,isp1,ipw2,isp2,ifreq)
           end do
         end do
       end do
     end do
   end do

!  Perform deallocations
   ABI_DEALLOCATE(freq)
   ABI_DEALLOCATE(susd_non_dyn)
   ABI_DEALLOCATE(susmat_dyn)

 end if !condition nfreqsus

!Performs deallocations
 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(gbound_diel)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(irrzondiel)
 ABI_DEALLOCATE(kg_diel)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(npwtot)
 ABI_DEALLOCATE(phnonsdiel)
 ABI_DEALLOCATE(symrec)

 if(mkmem==0)then
!  Sequential case
   if(mpi_enreg%paral_compil_kpt==0)then
!    Unit dtfil%unkg was opened in kpgio
     close (unit=dtfil%unkg,status='delete')

!    Parallel case
   else if(mpi_enreg%paral_compil_kpt==1)then

!    All procs close the file dtfil%unkg
     close(unit=dtfil%unkg)
     if(mpi_enreg%me==0)then
!      only proc 0 delete the file
       open(unit=dtfil%unkg,file=dtfil%fnametmp_kgs,form='unformatted',status='unknown')
       close(unit=dtfil%unkg,status='delete')
     end if

   end if
 end if

 call WffClose(wff1,ierr)

!Clean the header
 call hdr_clean(hdr)

 if(mpi_enreg%paral_compil_kpt==1)then
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 end if

 write(message, '(a,a)' ) ch10,' suscep : exiting '
 call wrtout(std_out,message,'COLL')

 call status(0,dtfil%filstat,iexit,level,'exit          ')
 call timab(84,2,tsec)

end subroutine suscep
!!***
