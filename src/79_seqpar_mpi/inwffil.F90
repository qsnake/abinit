!{\src2tex{textfont=tt}}
!!****f* ABINIT/inwffil
!! NAME
!! inwffil
!!
!! FUNCTION
!! Do initialization of wavefunction files.
!! Also call other relevant routines for this initialisation
!!    (initialization of wavefunctions from scratch or from file,
!!     translations of wavefunctions, ...)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, AR, MB, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ask_accurate= if 1, the wavefunctions and eigenvalues must be
!!    accurate, that is, they must come from a k point that is
!!    symmetric of the needed k point, with a very small tolerance,
!!    the disk file contained sufficient bands to initialize all of them,
!!    the spinor and spin-polarisation characteristics must be identical
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=effective kinetic energy planewave cutoff (hartree), beyond
!!    which the coefficients of plane waves are zero
!!  ecut_eff=effective kinetic energy planewave cutoff (hartree), needed
!!    to generate the sphere of plane wave
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  formeig=explained above
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ireadwf=option parameter described above for wf initialization
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!    to be initialized here.
!!  kg(3,mpw*mkmem)=dimensionless coords of G vecs in basis sphere at k point
!!  kptns(3,nkpt)=reduced coords of k points
!!  localrdwf=(for parallel case) if 1, the wffnm  file is local to each machine
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mkmem=maximum number of k-points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwarr(nkpt)=array holding npw for each k point.
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in space group
!!  occ(mband*nkpt*nsppol)=occupations (from disk or left at their initial value)
!!  optorth= 1 if the WFS have to be orthogonalized; 0 otherwise
!!  prtvol=control print volume and debugging
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!  unkg=unit number for storage of basis sphere data: stores indirect
!!   indexing array and integer coordinates for all planewaves in basis
!!   sphere for each k point being considered
!!  unwff1,unwfnow=
!!            unit numbers for files wffnm and wft1nm.
!!  wffnm=name (character data) of file for input wavefunctions.
!!  wft1nm= name (character data) of file used for temporary wf storage.
!!
!! OUTPUT
!!  wff1  = structure information for files wffnm .
!!  wffnow= structure information for wf file wft1nm
!!  if ground state format (formeig=0):
!!    eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), (Ha)
!!  if respfn format (formeig=1):
!!    eigen(2*mband*mband*nkpt*nsppol)=matrix of eigenvalues
!!                                     (input or init to large number), (Ha)
!! Conditional output (returned if mkmem/=0):
!!  cg(2,mcg)=complex wf array
!!    be careful : an array of size cg(2,npw*nspinor), as used
!!    in the response function code, is not enough !
!!  wvl <type(wvl_data)>=all wavelets data.
!!
!! NOTES
!! Detailed description :
!!  Initialize unit wff1%unwff for input of wf data if ireadwf=1
!!  Opens file on unit wffnow%unwff
!!   if the storage on disk is needed (mkmem==0)
!!  Initializes wf data on wffnow%unwff, by calling the appropriate routine.
!!
!! formeig option (format of the eigenvalues and occupations) :
!!   0 => ground-state format (initialisation of
!!        eigenvectors with random numbers, vector of eigenvalues,
!!        occupations are present)
!!   1 => respfn format (initialisation of
!!        eigenvectors with 0 s, hermitian matrix of eigenvalues)
!!
!! ireadwf options:
!!   0 => initialize with random numbers or 0 s
!!   1 => read from disk file wff1, initializing higher bands
!!        with random numbers or 0 s if not provided in disk file
!!
!! The wavefunctions after this initialisation are stored in unit wffnow%unwff
!!
!! WARNINGS
!! The symmetry operations are used to translate the data from one
!! k point to another, symmetric, k point.
!! They can be completely different from the symmetry operations
!! contained on the disk file. No check is performed between the two sets.
!!
!! Occupations will not be modified nor output,
!!  in the present status of this routine.
!!
!! If ground state format (formeig=0) occ(mband*nkpt*nsppol) was output.
!! NOT OUTPUT NOW !
!!
!! PARENTS
!!      gstate,loop3dte,loper3,nonlinear,respfn
!!
!! CHILDREN
!!      copy_mpi_enreg,destroy_mpi_enreg,handle_ncerr,hdr_check,hdr_clean
!!      hdr_io,hdr_io_etsf,hdr_io_netcdf,ini_wf_netcdf,leave_new,listkk
!!      matr3inv,newkpt,prep_kpgio,rhoij_copy,timab,wffkg,wffopen,wfsinp,wrtout
!!      wvl_wfsinp_disk,wvl_wfsinp_scratch,xcomm_init,xcomm_self,xdefineoff
!!      xmaster_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inwffil(ask_accurate,cg,dtset,ecut,ecut_eff,eigen,exchn2n3d,&
&           formeig,gmet,hdr,ireadwf,istwfk,kg,kptns,localrdwf,mband,&
&           mcg,mkmem,mpi_enreg,mpw,nband,ngfft,nkpt,npwarr,&
&           nsppol,nsym,occ,optorth,rprimd,symafm,symrel,tnons,unkg,wff1,&
&           wffnow,unwff1,unwfnow,wffnm,wft1nm, wvl)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_wffile
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_header,          only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inwffil'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_61_ionetcdf
 use interfaces_66_wfs
 use interfaces_67_common
 use interfaces_79_seqpar_mpi, except_this_one => inwffil
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: ask_accurate,exchn2n3d,formeig,ireadwf,localrdwf,mband,mcg,mkmem,mpw
 integer, intent(in) :: nkpt,nsppol,nsym,optorth,unkg,unwff1,unwfnow
 real(dp), intent(in) :: ecut,ecut_eff
 character(len=fnlen), intent(in) :: wffnm,wft1nm
 type(MPI_type),intent(inout),target :: mpi_enreg
 type(dataset_type), intent(in) :: dtset
 type(hdr_type), intent(inout) :: hdr
 type(wffile_type), intent(out) :: wff1,wffnow
 type(wvl_data), intent(inout) :: wvl
 integer, intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),ngfft(18)
 integer, intent(in) :: npwarr(nkpt),symafm(nsym),symrel(3,3,nsym)
 integer, intent(in),target :: nband(nkpt*nsppol)
 real(dp), intent(out),target :: cg(2,mcg),eigen((2*mband)**formeig*mband*nkpt*nsppol)
 real(dp), intent(in) :: gmet(3,3),kptns(3,nkpt),rprimd(3,3),tnons(3,nsym)
 real(dp), intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
! integer,save :: count=0 ! appears commented out below
! integer :: natt,ndim,nvar,uid ! appears commented out below
 integer :: accesswff,accurate,ceksp,debug,doorth,fform,fform_dum,fill
 integer :: headform0,iband,ibg,ibg0,icg,icg0,icgsft,ieigsft,ierr,ii
 integer :: ikassoc,ikpt,ikpt0,ikptsp,ikptsp0,imax,increase_nkassoc,isppol,isppol0
 integer :: master,mband0,mband0_rd,mband_eff,mcg_disk,me,me0,mkmem0,mpw0
 integer :: my_nspinor,my_nspinor0,nband_k,nband0_k
 integer :: nkassoc,nkpt0,npw,npw0,nspinor0,nspinor_eff,nsppol0,nsppol_eff,nsppol2nspinor
 integer :: rdwr,restart,restartpaw,spaceComm,spaceComm_io,sppoldbl,sppoldbl_eff,squeeze
 real(dp) :: dksqmax,ecut0
 character(len=500) :: message
 type(hdr_type) :: hdr0
 integer :: ngfft0(18)
 integer,allocatable :: indkk0(:,:),indx(:),istwfk0(:),kg0(:,:)
 integer,allocatable :: nband0_rd(:),npwarr0(:),npwi(:),npwtot0(:)
 integer,allocatable,target :: indkk(:,:),nband0(:)
 integer, pointer :: indkk_eff(:,:),nband_eff(:)
 logical,allocatable :: my_kpt(:)
 real(dp) :: gmet0(3,3),gprim0(3,3),rprim0(3,3),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),kptns0(:,:)
 real(dp),pointer :: cg_eff(:,:),eigen_eff(:)
 type(MPI_type),pointer :: mpi_enreg0

#if defined HAVE_TRIO_NETCDF
 integer :: ncid_hdr,ncerr
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' inwffil: enter'
!count=count+1
!write(std_out,*)' count=',count
!if(count==13)stop
!stop
!ENDDEBUG

!Keep track of total time spent in inwffil
 call timab(710,1,tsec)
 call timab(711,1,tsec)

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart_4d)
 call xcomm_self(spaceComm_io)
 if (mpi_enreg%mode_para=='b') spaceComm_io= mpi_enreg%commcart_3d

!Init me
 call xme_init(mpi_enreg,me,option_comm=2)
!Init master
 call xmaster_init(mpi_enreg,master)

!Check the validity of formeig
 if(formeig/=0.and.formeig/=1)then
   write(message, '(a,a,a,a,i12,a)' ) ch10,&
&   ' inwffil: BUG -',ch10,&
&   '  formeig=',formeig,' , but the only allowed values are 0 or 1.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 ngfft0(:)=ngfft(:)

!Default value for headform0 (will be needed later, to read wf blocks)
 headform0=0

!If the input data are on disk, determine the kind of restart
 if (ireadwf==1)then

   accesswff=dtset%accesswff
   if(localrdwf==0)accesswff=IO_MODE_FORTRAN_MASTER ! This is in case the wff file must be read by only the master proc

!  DEBUG
!  write(std_out,*)' inwffil : will open file ',trim(wffnm)
!  write(std_out,*)' inwffil : on unit', unwff1
!  call flush(6)
!  ENDDEBUG
   call WffOpen(accesswff,spaceComm,wffnm,ierr,wff1,master,me,unwff1,spaceComm_io)

!  DEBUG
!  write(std_out,*)' inwffil : before hdr_io '
!  call flush(6)
!  ENDDEBUG

#if defined HAVE_TRIO_NETCDF
!  DEBUG
!  if (accesswff == 2) then
!  ncid_hdr = wff1%unwff
!  ncerr = nf90_Inquire(ncid=ncid_hdr,nDimensions=ndim,nVariables=nvar,nAttributes=natt,unlimitedDimId=uid)
!  call handle_ncerr(ncerr, " general Inquire ")
!  write(std_out,*) 'inwffil : found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
!  end if
!  ENDDEBUG
#endif

!  Initialize hdr0 (sent to all procs), thanks to reading of wff1
   rdwr=1
   if ( ANY(wff1%accesswff == (/IO_MODE_FORTRAN_MASTER, IO_MODE_FORTRAN, IO_MODE_MPI/) )) then
     call hdr_io(fform_dum,hdr0,rdwr,wff1)
#if defined HAVE_TRIO_NETCDF
   else if (wff1%accesswff == IO_MODE_NETCDF) then
     call hdr_io_netcdf(fform_dum,hdr0,rdwr,wff1)
#endif
#if defined HAVE_TRIO_ETSF_IO
   else if (wff1%accesswff == IO_MODE_ETSF) then
     call hdr_io_etsf(fform_dum, hdr0, rdwr, wff1%unwff)
#endif
   end if

   write(std_out, '(a)')' inwffil : before hdr_check '
   write(message, '(a,a)' )&
&   ' inwffil : examine the header of disk file ',trim(wffnm)
   call wrtout(std_out,message,'COLL')


!  Check hdr0 versus hdr
!  (and from now on ignore header consistency and write new info
!  to header for each file)
   if (dtset%usewvl == 0) then
!    wait for plane waves.
     fform=2
   else
!    wait for wavelets.
     fform = 200
   end if
!  call hdr_check(fform,fform_dum,hdr,hdr0,'COLL',restart)
   call hdr_check(fform,fform_dum,hdr,hdr0,'PERS',restart,restartpaw)

   nkpt0=hdr0%nkpt
   nsppol0=hdr0%nsppol
   headform0=hdr0%headform


   write(message, '(a,a)' )&
&   '-inwffil : will read wavefunctions from disk file ',trim(wffnm)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

 else

   restart=1 ; restartpaw=0

!  Fill some data concerning an hypothetical file to be read
!  This is to allow the safe use of same routines than with ireadwf==1.
   nkpt0=nkpt ; nsppol0=nsppol

 end if ! end ireadwf

 sppoldbl=1
 if(minval(symafm(:))==-1)then
   if(nsppol0==1 .and. nsppol==2)sppoldbl=2
 end if

 ABI_ALLOCATE(indkk,(nkpt*sppoldbl,6))
 ABI_ALLOCATE(istwfk0,(nkpt0))
 ABI_ALLOCATE(kptns0,(3,nkpt0))
 ABI_ALLOCATE(nband0,(nkpt0*nsppol0))
 ABI_ALLOCATE(npwarr0,(nkpt0))

 if(restart==2)then ! restart with translations

   ecut0=hdr0%ecut_eff
   istwfk0(1:nkpt0)=hdr0%istwfk(1:nkpt0)
   kptns0(1:3,1:nkpt0)=hdr0%kptns(1:3,1:nkpt0)
   nband0(1:nkpt0*nsppol0)=hdr0%nband(1:nkpt0*nsppol0)
   ngfft0(1:3)=hdr0%ngfft(1:3)
   npwarr0(1:nkpt0)=hdr0%npwarr(1:nkpt0)
   nspinor0=hdr0%nspinor
   rprim0(:,:)=hdr0%rprimd(:,:)
   mpw0=maxval(npwarr0(:))

!  Compute reciprocal space metric gmet for unit cell of disk wf
   call matr3inv(rprim0,gprim0)
   do ii=1,3
     gmet0(:,ii)=gprim0(1,:)*gprim0(1,ii)+&
&     gprim0(2,:)*gprim0(2,ii)+&
&     gprim0(3,:)*gprim0(3,ii)
   end do

   if (mpi_enreg%mode_para=='b') then
     ABI_ALLOCATE(mpi_enreg0,)
     call copy_mpi_enreg(mpi_enreg,mpi_enreg0,0)
     ABI_ALLOCATE(kg0,(3,mpw0*nkpt0))
     ABI_ALLOCATE(npwtot0,(nkpt0))
     message="tmpfil"
     call prep_kpgio(dtset%accesswff,ecut0,dtset%exchn2n3d,gmet0,istwfk0,kg0,kptns0,&
&     message,dtset%mgfft,nkpt0,'PERS',mpi_enreg0,mpw0,nband0,nkpt0,&
&     npwarr0,npwtot0,nsppol0,tmp_unit)
     ABI_DEALLOCATE(kg0)
     ABI_DEALLOCATE(npwtot0)
   else
     mpi_enreg0 => mpi_enreg
   end if

!  At this stage, the header of the file wff1i%unwff is read, and
!  the pointer is ready to read the first wavefunction block.

!  Compute k points from input file closest to the output file
   call listkk(dksqmax,gmet0,indkk,kptns0,kptns,nkpt0,nkpt,nsym,sppoldbl,&
&   symafm,symrel,1)

 else if (restart==1) then ! direct restart

!  Fill variables that must be the same, as determined by hdr_check.f
!  This is to allow the safe use of the same routines than with restart==2.
   nspinor0=dtset%nspinor
   ecut0=ecut_eff
   gmet0(:,:)=gmet(:,:)
   istwfk0(:)=istwfk(:)
   kptns0(:,:)=kptns(:,:)
   npwarr0(:)=npwarr(:)
   mpw0=mpw

   do isppol=1,sppoldbl
     do ikpt=1,nkpt
       indkk(ikpt+(isppol-1)*nkpt,1)=ikpt
       indkk(ikpt+(isppol-1)*nkpt,2:6)=0
     end do
   end do
   dksqmax=0.0_dp

!  The treatment of nband0 asks for some care
   if(ireadwf==0)then
     nband0(:)=0
   else
     nband0(1:nkpt0*nsppol0)=hdr0%nband(1:nkpt0*nsppol0)
   end if

   mpi_enreg0 => mpi_enreg

 else
   mpi_enreg0 => mpi_enreg
 end if

 call xme_init(mpi_enreg0,me0,option_comm=2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Before hdr_clean:
!If restartpaw==1, store hdr0%pawrhoij in hdr%pawrhoij; else if restartpaw==0,
!hdr%pawrhoij(:)has been initialized in hdr_init.
 if(restartpaw==1) then
   call rhoij_copy(hdr0%pawrhoij,hdr%pawrhoij)
 end if

 call timab(711,2,tsec)
 call timab(712,1,tsec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!At this stage, all the relevant information from the header of the disk file,
!has been exploited, and stored in variables, on all processors.
!It is also contained in hdr0
!(on all processors, except if restart=1 and localrdwf=0,
!in which case it is only on the master)
!These information might be changed later, while processing the
!wavefunction data, and converting it. The variable hdr0 might be kept
!for further checking, or reference, or debugging, but at present,
!it is simpler to close it. The other header, hdr, will be used
!for the new file, if any.

 if(ask_accurate==1)then

!  Check whether the accuracy requirements might be fulfilled
   if(ireadwf==0)then
     write(message, '(a,a,a,a,a, a,a,a,a,a, a,a)' ) ch10, &
&     ' inwffil: ERROR ',ch10,&
&     '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&     '  present calculation. It was asked that the wavefunctions be accurate,',ch10,&
&     '  but they were not even read.',ch10,&
&     '  Action: use a wf file, with ireadwf/=0.'
     call wrtout(std_out,  message,'PERS')
     call leave_new('PERS')
   end if
   if(dksqmax>tol12)then
     write(message, '(12a,es16.6,4a)' ) ch10, &
&     ' inwffil: ERROR ',ch10,&
&     '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&     '  present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&     '  at least one of the k points could not be generated from a symmetrical one.',ch10,&
&     '  dksqmax=',dksqmax,ch10,&
&     '  Action: check your wf file and k point input variables',ch10,&
&     '    (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
     call wrtout(std_out,  message,'PERS')
     call leave_new('PERS')
   end if
   if(dtset%nspinor/=nspinor0)then
     write(message, '(a,a,a,a,a, a,a,a,a,a, a,a,2i5,a,a)' ) ch10, &
&     ' inwffil: ERROR ',ch10,&
&     '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&     '  present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&     '  nspinor differs in the file from the actual nspinor.',ch10,&
&     '  nspinor,nspinor0=',dtset%nspinor,nspinor0,ch10,&
&     '  Action: check your wf file, and nspinor input variables.'
     call wrtout(std_out,  message,'PERS')
     call leave_new('PERS')
   end if
   if((nsppol>nsppol0 .and. sppoldbl==1) .or. nsppol<nsppol0 ) then
     write(message, '(a,a,a,a,a, a,a,a,a,a, a,a,3i5,a,a)' ) ch10, &
&     ' inwffil: ERROR ',ch10,&
&     '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&     '  present calculation. It was asked that the wavefunctions be accurate, but',ch10,&
&     '  the nsppol variables do not match in the file and in the actual calculation',ch10,&
&     '  nsppol,nsppol,sppoldbl=',dtset%nspinor,nspinor0,sppoldbl,ch10,&
&     '  Action: check your wf file, and nsppol input variables.'
     call wrtout(std_out,  message,'PERS')
     call leave_new('PERS')
   end if

!  Now, check the number of bands
   accurate=1
   do isppol=1,nsppol
     do ikpt=1,nkpt
       ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
       ikptsp =ikpt +(isppol-1)*nkpt
       ikptsp0=ikpt0+(isppol-1)*(2-sppoldbl)*nkpt0
       if(nband0(ikptsp0)<nband(ikptsp))accurate=0
     end do
   end do
   if(accurate==0)then
     write(message, '(a,a,a,a,a, a,a,a,a,a, a,a)' ) ch10, &
&     ' inwffil: ERROR ',ch10,&
&     '  The file ',trim(wffnm),' cannot be used to start the ',ch10,&
&     '  present calculation. It was asked that the wavefunctions be accurate,',ch10,&
&     '  but the number of bands differ in the file and in the actual calculation.',ch10,&
&     '  Action: use a wf file with the correct characteristics.'
     call wrtout(std_out,  message,'PERS')
     call leave_new('PERS')
   end if

 end if

!Flag: do we need to translate WF to (from) spinors ?
 nsppol2nspinor=0
 if (nsppol0==2.and.dtset%nspinor==2) nsppol2nspinor=+1
 if (nspinor0==2.and.nsppol==2) nsppol2nspinor=-1

!Take into account parallism over spinors
 my_nspinor =max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 my_nspinor0=max(1,nspinor0/mpi_enreg0%nproc_spin)

!Not all bands might be read, if not needed to fill the wavefunctions
 mband0=maxval(nband0(1:nkpt0*nsppol0))
 mband0_rd=min(mband0,(mband/dtset%nspinor)*nspinor0)

!Open the wf file on unit wffnow%unwff if mkmem==0
!Also allocate temporary array, and write header at top of wffnow%unwff

 if (mkmem==0) then

#if defined HAVE_TRIO_NETCDF
   if(dtset%accesswff==IO_MODE_NETCDF) then
     accesswff=IO_MODE_NETCDF
!    Create empty netCDF file
     ncerr = nf90_create(path=trim(wft1nm), cmode=NF90_CLOBBER, ncid=ncid_hdr)
     call handle_ncerr(ncerr," inwffil create netcdf wavefunction file")
     ncerr = nf90_close(ncid_hdr)
     call handle_ncerr(ncerr," inwffil close netcdf wavefunction file")
   else if (dtset%accesswff == IO_MODE_ETSF) then
     write (std_out,*) "FIXME: ETSF I/O support in inwffil"
   end if
#endif

   write(message, '(a,i4,a,a)' ) &
&   ' inwffil about to open unit',unwfnow,' for file=',trim(wft1nm)
   call wrtout(std_out,  message,'PERS')
   call WffOpen(dtset%accesswff,spaceComm,wft1nm,ierr,wffnow,master,me,unwfnow)

   rdwr=2 ; fform=2
   if (wffnow%accesswff /= IO_MODE_NETCDF) then
     call hdr_io(fform,hdr,rdwr,wffnow)
#if defined HAVE_TRIO_NETCDF
   else if (wffnow%accesswff == IO_MODE_NETCDF) then
     call hdr_io_netcdf(fform,hdr,rdwr,wffnow)

     call ini_wf_netcdf(mpw,wffnow%unwff,formeig)
   else if (wffnow%accesswff == IO_MODE_ETSF) then
     write (std_out,*) "FIXME: ETSF I/O support in inwffil"
#endif
   end if


!  Note that the cg_disk array is allocated when cg uses no space (mkmem==0)
!  This dimensioning allows to treat all cases (note: not compatible with parallelization over spinors)
   mcg_disk=max(mpw0*my_nspinor0*mband0_rd,mpw*my_nspinor*mband)
   ABI_ALLOCATE(cg_disk,(2,mcg_disk))

!  Define offsets, in case of MPI I/O
   call WffKg(wffnow,1)   ! option optkg in wfsinp
   call xdefineOff(formeig,wffnow,mpi_enreg,hdr%nband,hdr%npwarr,hdr%nspinor,hdr%nsppol,hdr%nkpt)
 end if


!****************************************************************************
!If needed, transfer the input wf from disk to core memory
!(in the parallel case, it allows to change localrdwf=0 in localrdwf=1)

 mkmem0=0

 if(mpi_enreg%paral_compil_kpt == 1 .or. mpi_enreg%paral_compil_fft == 1) then
   if(localrdwf==0 .and. mkmem==0)then
     write(message,'(a,a,a,a)')ch10,&
&     ' inwffil: BUG -',ch10,&
&     '  localrdwf==0 and mkmem==0 are not allowed together (yet)'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

 call timab(712,2,tsec)

!Here, treat reading wavefunctions with mkmem/=0, first step
 if(ireadwf==1 .and. mkmem/=0)then

   call timab(713,1,tsec)

!  if(restart==1 .and. ireadwf==1 .and. mkmem/=0)then

!  Compute table of k point associations. Make a trial choice for nkassoc.
   nkassoc=(nkpt/nkpt0+1)*2
   ABI_ALLOCATE(indkk0,(nkpt0,nkassoc))
!  Infinite loops are allowed in F90
   do
     indkk0(:,:)=0
     increase_nkassoc=0
     do ikpt=1,nkpt*sppoldbl
       ikpt0=indkk(ikpt,1)
       do ikassoc=1,nkassoc
         if(indkk0(ikpt0,ikassoc)==0)then
           indkk0(ikpt0,ikassoc)=ikpt
           exit
         end if
         if(nkassoc==ikassoc)increase_nkassoc=1
       end do
       if(increase_nkassoc==1)then
         ABI_DEALLOCATE(indkk0)
         nkassoc=2*nkassoc
         ABI_ALLOCATE(indkk0,(nkpt0,nkassoc))
         exit
       end if
     end do
     if(increase_nkassoc==0)exit
   end do

!  DEBUG
!  write(std_out,*)' inwffil: indkk0, nkassoc=',nkassoc
!  do ikpt0=1,nkpt0
!  write(std_out,*)' ikpt0,indkk0(ikpt0,1)=',ikpt0,indkk0(ikpt0,1)
!  end do
!  ENDDEBUG

!  DEBUG
!  write(std_out,*)' inwffil : indkk(:,1)=',indkk(:,1)
!  write(std_out,*)' inwffil : sppoldbl=',sppoldbl
!  ENDDEBUG

!  To treat the case (nsppol0=2,nspinor0=1)<->(nsppol=1,nspinor=2),
!  apply the following trick:
!  1- We call wfsinp with fake arguments (nsppol_eff and nspinor_eff)
!  2- We transform collinear polarized WF into spinors
!  or  spinors into collinear polarized WF
   if (nsppol2nspinor/=0.and.mkmem==0.and.dtset%usewvl==0) then
     write(message, '(10a)') ch10,&
&     ' inwffil : ERROR -',ch10,&
&     '  When mkmem=0, the wavefunction translator is unable',ch10,&
&     '  to interchange spin-polarized wfs and spinor wfs.',ch10,&
&     '  Action : use a non-spin-polarized wf to start a spinor wf,',ch10,&
&     '           and a non-spinor wf to start a spin-polarized wf.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  === Fake arguments definition for wfsinp
   if (nsppol2nspinor==0.or.dtset%usewvl/=0) then
     indkk_eff => indkk
     nband_eff => nband
     eigen_eff => eigen
     cg_eff => cg
     nspinor_eff=dtset%nspinor;nsppol_eff=nsppol;sppoldbl_eff=sppoldbl
     mband_eff=maxval(nband_eff(1:nkpt*nsppol_eff))
   else if (nsppol2nspinor==1.and.mkmem/=0) then
     nsppol_eff=2;nspinor_eff=1;sppoldbl_eff=1
     ABI_ALLOCATE(indkk_eff,(nkpt*sppoldbl_eff,6))
     ABI_ALLOCATE(nband_eff,(nkpt*nsppol_eff))
     indkk_eff(1:nkpt,1:6)   =indkk(1:nkpt,1:6)
     nband_eff(1:nkpt)       =nband(1:nkpt)/2
     nband_eff(1+nkpt:2*nkpt)=nband(1:nkpt)/2
     mband_eff=maxval(nband_eff(1:nkpt*nsppol_eff))
     eigen_eff => eigen
     cg_eff => cg
   else if (nsppol2nspinor==-1.and.mkmem/=0) then
!    WARNING: MT 07072011 -> this is memory consuming
!    A copy a spinorial WF (and eigenvalues) is temporary kept in memory;
!    But the case (nspinor=2 => nsppol=2) might be rare
!    and only useful for testing purposes.
!    => print a warning for the user
!    NOTE: in that case (nsppol=2), parallelization over spinors is not activated
     write(message, '(8a)') ch10,&
&     ' inwffil : WARNING -',ch10,&
&     '  In the case of spinor WF read from disk and converted into',ch10,&
&     '  spin-polarized non-spinor WF, the WF translator is memory',ch10,&
&     '  consuming (a copy a spinor WF is temporary stored in memory).'
     call wrtout(std_out,message,'COLL')
     nsppol_eff=1;nspinor_eff=2;sppoldbl_eff=1
     ABI_ALLOCATE(indkk_eff,(nkpt*sppoldbl_eff,6))
     ABI_ALLOCATE(nband_eff,(nkpt*nsppol_eff))
     indkk_eff(1:nkpt,1:6)=indkk(1:nkpt,1:6)
     nband_eff(1:nkpt)    =2*nband(1:nkpt)
     mband_eff=maxval(nband_eff(1:nkpt*nsppol_eff))
     ABI_ALLOCATE(eigen_eff,((2*mband_eff)**formeig*mband_eff*nkpt*nsppol_eff))
     ABI_ALLOCATE(cg_eff,(2,mpw0*nspinor_eff*mband_eff*mkmem*nsppol_eff))
   end if

!  === nband0 argument definition for wfsinp
   squeeze=0
   if(mkmem/=0)then
     ABI_ALLOCATE(nband0_rd,(nkpt0*nsppol0))
     nband0_rd(:)=0
     do isppol=1,nsppol_eff
       do ikpt=1,nkpt
         ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
         isppol0=min(isppol,nsppol0)
         ikptsp =ikpt +(isppol -1)*nkpt
         ikptsp0=ikpt0+(isppol0-1)*(2-sppoldbl)*nkpt0
         nband0_k=min(nband0(ikptsp0),(nband_eff(ikptsp)/nspinor_eff)*nspinor0)
         nband0_rd(ikptsp0)=max(nband0_rd(ikptsp0),nband0_k)
         npw0=npwarr0(ikpt0)
         npw =npwarr (ikpt)
         if(npw0*nspinor0*nband0_k > npw*nspinor_eff*nband_eff(ikptsp))squeeze=1
       end do
     end do
     if(squeeze==1)then
       mcg_disk=mpw0*my_nspinor0*mband0_rd
       ABI_ALLOCATE(cg_disk,(2,mcg_disk))
     else
       if(mpi_enreg0%paral_compil_kpt == 1 .or. mpi_enreg0%paral_compil_fft == 1)then
         if(localrdwf==0)then
           mcg_disk=mpw0*my_nspinor0*mband0_rd
           ABI_ALLOCATE(cg_disk,(2,mcg_disk))
         end if
       end if
     end if
   end if

   call timab(713,2,tsec)
   call timab(714,1,tsec)

!  === call to wfsinp
   if (dtset%usewvl == 0) then
     call wfsinp(cg_eff,cg_disk,ecut,ecut0,ecut_eff,eigen,&
&     exchn2n3d,formeig,gmet,gmet0,headform0,&
&     indkk_eff,indkk0,istwfk,istwfk0,kptns,kptns0,localrdwf,&
&     mband_eff,mband0_rd,mcg,mcg_disk,mkmem,mpi_enreg,mpi_enreg0,mpw,mpw0,&
&     nband_eff,nband0_rd,ngfft,nkassoc,nkpt,nkpt0,npwarr,npwarr0,nspinor_eff,nspinor0,&
&     nsppol_eff,nsppol0,nsym,occ,optorth,dtset%prtvol,restart,rprimd,sppoldbl_eff,squeeze,&
&     symrel,tnons,wff1,wffnow)
     if (nsppol2nspinor/=0)  then
       ABI_DEALLOCATE(indkk_eff)
       ABI_DEALLOCATE(nband_eff)
     end if
   else
!    Read wavefunctions from file.
     call wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, 1, &
&     hdr%rprimd, wff1, wvl%wfs, wvl%descr, hdr%xred)
   end if

   call timab(714,2,tsec)
   call timab(715,1,tsec)

!  Now, update xyz0 variables, for use in newkpt
   nband0(:)=nband0_rd(:)

!  If squeeze, the conversion was done in wfsinp, so no conversion left.
   if(squeeze==1)then
     ecut0=ecut_eff
     gmet0(:,:)=gmet(:,:)
     ABI_DEALLOCATE(kptns0)
     ABI_DEALLOCATE(istwfk0)
     ABI_DEALLOCATE(nband0)
     ABI_DEALLOCATE(npwarr0)
     ABI_ALLOCATE(kptns0,(3,nkpt))
     ABI_ALLOCATE(istwfk0,(nkpt))
     ABI_ALLOCATE(nband0,(nkpt*nsppol))
     ABI_ALLOCATE(npwarr0,(nkpt))
     kptns0(:,:)=kptns(:,:)
     istwfk0(:)=istwfk(:)
     npwarr0(:)=npwarr(:)
     nband0(:)=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         ikpt0=indkk(ikpt+(isppol-1)*(sppoldbl-1)*nkpt,1)
         isppol0=min(isppol,nsppol0)
         ikptsp =ikpt +(isppol -1)*nkpt
         ikptsp0=ikpt0+(isppol0-1)*(sppoldbl-1)*nkpt0
         nband0(ikptsp)=(nband0_rd(ikptsp0)/nspinor0)*dtset%nspinor
#if defined FC_FUJITSU
!        DEBUG by MM for VPP
!        write(std_out,*) 'nband0(',ikptsp,'), nband0_rd(',ikptsp0,')=',&
!        &     nband0(ikptsp),nband0_rd(ikptsp0)
!        write(std_out,*) 'ikpt=',ikpt,' ikpt0=',ikpt0
!        write(std_out,*) 'isppol=',isppol,' isppol0=',isppol0
         write(std_out,*) 'ikptsp=',ikptsp,' ikptsp0=',ikptsp0         !this four lines
         write(std_out,*) 'nspinor=',dtset%nspinor,' nspinor0=',nspinor0     ! necessary
         write(std_out,*) 'nband0   (',ikptsp,')=',nband0(ikptsp)      !  for VPP
         write(std_out,*) 'nband0_rd(',ikptsp0,')=',nband0_rd(ikptsp0) !don''t comment!
!        END DEBUG by MM
#endif
       end do
     end do
     do ikpt=1,nkpt
       indkk(ikpt,1)=ikpt
       indkk(ikpt,2:6)=0
     end do
!    This transfer must come after the nband0 transfer
     nspinor0=dtset%nspinor
     nkpt0=nkpt
     nsppol0=nsppol
   end if ! end squeeze == 1

!  The input wavefunctions have been transferred from disk to core memory
   mkmem0=mkmem

   ABI_DEALLOCATE(indkk0)
   ABI_DEALLOCATE(nband0_rd)
   if(squeeze==1)ABI_DEALLOCATE(cg_disk)
   if(mpi_enreg%paral_compil_fft == 1 .or. mpi_enreg%paral_compil_kpt == 1 )then
     if(localrdwf==0)ABI_DEALLOCATE(cg_disk)
   end if

   call timab(715,2,tsec)

 else !ireadwf == 0
   if (dtset%usewvl == 1) then

     call timab(714,1,tsec)
!    Compute wavefunctions from input guess.
     call wvl_wfsinp_scratch(dtset, mpi_enreg, hdr%rprimd, wvl, hdr%xred)
     call timab(714,2,tsec)
   end if
 end if

 call timab(716,1,tsec)

!=== Eventual conversion of WF into (from) spinors
 if (dtset%usewvl==0) then

!  ***** No conversion (standard case) ****
   if (nsppol2nspinor==0) then
     nspinor_eff=nspinor0;nsppol_eff=nsppol0;sppoldbl_eff=sppoldbl
     indkk_eff => indkk
     nband_eff => nband0

!    ***** Conversion from collinear to spinorial WF ****
   else if (nsppol2nspinor==1.and.mkmem/=0) then
!    Translate the WF and eigenvalues from nsppol=2 to nspinor=2
!    This is tricky (because we do not want to create a temporary array for cg)
     nsppol_eff=1;nspinor_eff=2;sppoldbl_eff=1
     ABI_ALLOCATE(indkk_eff,(nkpt*sppoldbl_eff,6))
     ABI_ALLOCATE(nband_eff,(nkpt0*nsppol_eff))
     indkk_eff(1:nkpt,1:6)=indkk(1:nkpt,1:6)
     nband_eff(1:nkpt0)=2*nband0(1:nkpt0)
!    Compute some shifts from isspol0=1 to isppol0=2
     imax=0;icgsft=0;ieigsft=0
     ABI_ALLOCATE(my_kpt,(nkpt0))
     do ikpt0=1,nkpt0
       nband0_k=nband0(ikpt0);nband_k=nband(ikpt0)
       my_kpt(ikpt0)=(minval(abs(mpi_enreg0%proc_distrb(ikpt0,1:nband_k,1)-me))==0)
       ieigsft=ieigsft+(2*nband0_k)**formeig*nband0_k
       if(my_kpt(ikpt0)) then
         imax=imax+nband0_k;icgsft=icgsft+nband0_k*npwarr0(ikpt0)
       end if
     end do
!    --- First version: no parallelization over spinors
     if (mpi_enreg0%paral_spin==0) then
!      Compute some useful indexes
       ABI_ALLOCATE(indx,(2*imax))
       ABI_ALLOCATE(npwi,(imax))
       ii=0;icg=0
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             ii=ii+1;npwi(ii)=npw0
             indx(2*ii-1)=icg+mpw0;indx(2*ii)=icg+2*mpw0
             icg=icg+4*mpw0
           end do
         end if
       end do
!      Expand WF in cg (try to use the whole array)
       ii=nsppol0*imax;icg0=nsppol0*icgsft
       do isppol=nsppol0,1,-1
         do ikpt0=nkpt0,1,-1
           if(my_kpt(ikpt0)) then
             nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
             do iband=nband0_k,1,-1
               icg0=icg0-npw0;if (indx(ii)<icg0) stop "Bug: Unable to read WF !"
               cg(:,indx(ii)+1:indx(ii)+npw0)=cg(:,icg0+1:icg0+npw0)
               ii=ii-1
             end do
           end if
         end do
       end do
!      Convert polarized WF into spinors
       ii=1
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             npw0=npwi(ii)
             cg(:,indx(2*ii-1)-mpw0+1:indx(2*ii-1)-mpw0+npw0)=cg(:,indx(ii)+1:indx(ii)+npw0)
             cg(:,indx(2*ii  )+mpw0+1:indx(2*ii  )+mpw0+npw0)=cg(:,indx(ii+imax)+1:indx(ii+imax)+npw0)
             ii=ii+1
           end do
         end if
       end do
!      Compress new cg array (from mpw to npw) and cancel zero-components
       icg0=0;icg=0
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             cg(:,icg0       +1:icg0+  npw0)=cg(:,icg+1:icg+npw0)
             cg(:,icg0+  npw0+1:icg0+2*npw0)=zero
             cg(:,icg0+2*npw0+1:icg0+3*npw0)=zero
             cg(:,icg0+3*npw0+1:icg0+4*npw0)=cg(:,icg+3*mpw0+1:icg+3*mpw0+npw0)
             icg0=icg0+4*npw0;icg=icg+4*mpw0
           end do
         end if
       end do
!      --- Second version: parallelization over spinors
     else
!      Compute some useful indexes
       ABI_ALLOCATE(indx,(imax))
       ABI_ALLOCATE(npwi,(imax))
       ii=0;icg=0
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             ii=ii+1;npwi(ii)=npw0
             indx(ii)=icg+mpi_enreg0%me_spin*mpw0
             icg=icg+2*mpw0
           end do
         end if
       end do
!      Expand WF in cg
       ii=(mpi_enreg0%me_spin+1)*imax;icg0=(mpi_enreg0%me_spin+1)*icgsft
       do ikpt0=nkpt0,1,-1
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=nband0_k,1,-1
             icg0=icg0-npw0;if (indx(ii)<icg0) stop "Bug: Unable to read WF !"
             cg(:,indx(ii)+1:indx(ii)+npw0)=cg(:,icg0+1:icg0+npw0)
             ii=ii-1
           end do
         end if
       end do
!      Compress new cg array (from mpw to npw) and cancel zero-components
       icg0=0;icg=0
       do ikpt0=1,nkpt0
         if(my_kpt(ikpt0)) then
           nband0_k=nband0(ikpt0);npw0=npwarr0(ikpt0)
           do iband=1,nband0_k
             if (mpi_enreg0%me_spin==0) then
               cg(:,icg0     +1:icg0+  npw0)=cg(:,icg+1:icg+npw0)
               cg(:,icg0+npw0+1:icg0+2*npw0)=zero
             else
               cg(:,icg0     +1:icg0+  npw0)=zero
               cg(:,icg0+npw0+1:icg0+2*npw0)=cg(:,icg+mpw0+1:icg+mpw0+npw0)
             end if
             icg0=icg0+2*npw0;icg=icg+2*mpw0
           end do
         end if
       end do
     end if
!    Translate eigenvalues
     ibg0=2*ieigsft;ibg=2*ieigsft
     do ikpt0=nkpt0,1,-1
       nband0_k=nband0(ikpt0)
       ibg0=ibg0-  nband0_k*(2*nband0_k)**formeig
       ibg =ibg -2*nband0_k*(2*nband0_k)**formeig
       if(my_kpt(ikpt0)) then
         do iband=nband0_k*(2*nband0_k)**formeig,1,-1
           eigen(2*iband-1+ibg)=eigen(iband+ibg0-ieigsft)
           eigen(2*iband  +ibg)=eigen(iband+ibg0)
         end do
       end if
     end do
     ABI_DEALLOCATE(indx)
     ABI_DEALLOCATE(npwi)
     ABI_DEALLOCATE(my_kpt)

!    ***** Conversion from spinorial to collinear WF ****
   else if (nsppol2nspinor==-1.and.mkmem/=0) then
!    In that case parallelization over spinors is never activated
     nsppol_eff=2;nspinor_eff=1;sppoldbl_eff=1
     ABI_ALLOCATE(indkk_eff,(nkpt*sppoldbl_eff,6))
     ABI_ALLOCATE(nband_eff,(nkpt0*nsppol_eff))
     indkk_eff(1:nkpt,1:6)=indkk(1:nkpt,1:6)
     nband_eff(1:nkpt0)        =nband0(1:nkpt0)/2
     nband_eff(1+nkpt0:2*nkpt0)=nband0(1:nkpt0)/2
!    Compute shifts from isspol0=1 to isppol0=2
     icgsft=0;ieigsft=0
     do ikpt0=1,nkpt0
       nband0_k=nband0(ikpt0);nband_k=nband(ikpt0)
       ieigsft=ieigsft+(nband0_k/2)*(nband0_k)**formeig
       if(minval(abs(mpi_enreg0%proc_distrb(ikpt0,1:nband_k,1)-me))==0) &
&       icgsft=icgsft+(nband0_k/2)*npwarr0(ikpt0)
     end do
!    Translate the WF and eigenvalues from nspinor=2 to nsppol=2
     icg0=0;icg=0;ibg=0
     do ikpt0=1,nkpt0
       nband0_k=nband0(ikpt0);nband_k=nband(ikpt0);npw0=npwarr0(ikpt0)
       if(minval(abs(mpi_enreg0%proc_distrb(ikpt0,1:nband_k,1)-me))==0) then
         do iband=1,nband0_k/2
           do ii=1,npw0
             cg(:,ii+icg)       =cg_eff(:,ii+icg0)
             cg(:,ii+icg+icgsft)=cg_eff(:,ii+icg0+3*npw0)
           end do
           icg0=icg0+4*npw0;icg=icg+npw0
         end do
         do iband=(nband0_k/2)*(nband0_k)**formeig,1,-1
           eigen(iband+ibg)        =eigen_eff(2*iband-1+2*ibg)
           eigen(iband+ibg+ieigsft)=eigen_eff(2*iband  +2*ibg)
!          occ(iband+ibg)        =occ_eff(2*iband-1+2*ibg)
!          occ(iband+ibg+ieigsft)=occ_eff(2*iband  +2*ibg)
         end do
       end if
       ibg=ibg+(nband0_k/2)*(nband0_k)**formeig
     end do
     ABI_DEALLOCATE(cg_eff)
     ABI_DEALLOCATE(eigen_eff)

   else
     write(message, '(10a)') ch10,&
&     ' inwffil : BUG -',ch10,&
&     '  unable to interchange nsppol and nspinor when mkmem=0 !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if
 end if

!Clean hdr0
 if (ireadwf==1)then
   if( restart==2 .or. localrdwf==1 .or. master==me)then
     call hdr_clean(hdr0)
   end if
 end if

 call timab(716,2,tsec)
 call timab(717,1,tsec)


!****************************************************************************
!Now, treat translation of wavefunctions if wavefunctions are planewaves

 ceksp=0 ; debug=0 ; doorth=1 ; fill=1

 if (dtset%usewvl == 0) then
!  DEBUG
!  write(std_out,*)' inwffil : will call newkpt, me= ',me,wff1%offwff,wffnow%offwff
!  write(std_out,*)' inwffil : will call newkpt, me= ',me
!  write(std_out,*)' inwffil : restart=',restart
!  if(restart==2)stop
!  ENDDEBUG
   if(mkmem/=0)then
     call newkpt(ceksp,cg,debug,ecut0,ecut,ecut_eff,eigen,exchn2n3d,&
&     fill,formeig,gmet0,gmet,headform0,indkk_eff,&
&     ab_out,ireadwf,istwfk0,istwfk,kg,kptns0,kptns,&
&     mband,mcg,mkmem0,mkmem,mpi_enreg0,mpi_enreg,&
&     mpw0,mpw,nband_eff,nband,ngfft0,nkpt0,nkpt,npwarr0,npwarr,&
&     nspinor_eff,dtset%nspinor,nsppol_eff,nsppol,nsym,occ,optorth,dtset%prtvol,&
&     restart,rprimd,sppoldbl_eff,symrel,tnons,unkg,wff1,wffnow)
   else
     call newkpt(ceksp,cg_disk,debug,ecut0,ecut,ecut_eff,eigen,exchn2n3d,&
&     fill,formeig,gmet0,gmet,headform0,indkk_eff,&
&     ab_out,ireadwf,istwfk0,istwfk,kg,kptns0,kptns,&
&     mband,mcg_disk,mkmem0,mkmem,mpi_enreg0,mpi_enreg,&
&     mpw0,mpw,nband_eff,nband,ngfft0,nkpt0,nkpt,npwarr0,npwarr,&
&     nspinor_eff,dtset%nspinor,nsppol_eff,nsppol,nsym,occ,optorth,dtset%prtvol,&
&     restart,rprimd,sppoldbl_eff,symrel,tnons,unkg,wff1,wffnow)
   end if
   if (nsppol2nspinor/=0)  then
     ABI_DEALLOCATE(nband_eff)
   end if

 end if ! dtset%usewvl == 0

!****************************************************************************

 ABI_DEALLOCATE(indkk)
 ABI_DEALLOCATE(istwfk0)
 ABI_DEALLOCATE(kptns0)
 ABI_DEALLOCATE(nband0)
 ABI_DEALLOCATE(npwarr0)
 if (restart==2.and.mpi_enreg%mode_para=='b') then
   call destroy_mpi_enreg(mpi_enreg0)
   ABI_DEALLOCATE(mpi_enreg0)
 else
   nullify(mpi_enreg0)
 end if

!If disk file was created, tell its name, also deallocate the big array
 if (mkmem==0) then
   write(message, '(a,a)' ) &
&   ' inwffil: copy wf on disk file ',trim(wft1nm)
   call wrtout(std_out,  message,'PERS')
   ABI_DEALLOCATE(cg_disk)
 end if

 write(message,*)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 call timab(717,2,tsec)
 call timab(710,2,tsec)

!DEBUG
!write(std_out,*)' inwffil: exit, write eigen '
!do ikpt=1,nkpt
!do iband=1,mband
!write(std_out,*)'ikpt,iband,eigen',ikpt,iband,eigen(iband+(ikpt-1)*mband)
!end do
!end do
!if(formeig==1)stop
!ENDDEBUG

end subroutine inwffil
!!***
