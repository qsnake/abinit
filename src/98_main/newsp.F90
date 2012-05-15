!{\src2tex{textfont=tt}}
!!****p* ABINIT/newsp
!! NAME
!! newsp
!!
!! FUNCTION
!! This program is used to take wavefunction data at a set of k-points, a
!! lattice constant, and an energy cutoff and write them out at a
!! new set of k-points, lattice constant, and energy cutoff.  Any
!! or all of these parameters is allowed to vary.  Even if no
!! parameters vary, this program may be used to re-order a data
!! set (i.e. change the order of ikpt and kptns), or to change the
!! wavefunction storage mode (istwfk).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, ZL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!! ABINIT coding rules are NOT followed strictly in this program.
!! Most of NEWSP functionalities are now included in ABINIT (from version 1.9).
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,bstruct_clean,bstruct_init,destroy_mpi_enreg,getmpw
!!      getng,handle_ncerr,hdr_clean,hdr_init,hdr_io,hdr_io_netcdf,hdr_update
!!      herald,importcml,indefo1,ini_wf_netcdf,instrng,intagm,inupper,invars1
!!      invars2,isfile,kpgsph,leave_new,listkk,matr3inv,metric,mpi_comm_rank
!!      mpi_comm_size,newkpt,nullify_mpi_enreg,wrtout,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program newsp

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_wffile
 use m_build_info
#if defined HAVE_MPI2
 use mpi
#endif
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif

 use m_header,  only : hdr_init, hdr_clean
 use m_ebands,  only : bstruct_init, bstruct_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'newsp'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_42_parser
 use interfaces_47_xml
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
 use interfaces_57_iovars
 use interfaces_59_io_mpi
 use interfaces_61_ionetcdf
 use interfaces_66_wfs
 use interfaces_67_common
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
!
! Set maximum array dimensions **********************************************
!mpsang=highest angular momentum allowed for consideration plus one
!msppol=maximum number of spin-polarization (default=2)
!msym=maximum number of symmetry operations in space group
 integer,parameter :: msppol=2,msym=48
!Fortran units for disk files
 integer,parameter :: iout=6
 integer :: wfinp=10,wfout=11
 integer :: fform=2
 integer :: bant0,bantot,bantot2,ceksp2,debug,doorth,fform1,fill,formeig,gscase
 integer :: headform,iband,ierr,ii,ikpt,ikpt1,ikpt2,ipsp,ireadwf,isppol
 integer :: itypat
 integer :: jdtset,lenstr,marr,mband1,mband2,mband_rd,mband_upper,mcg2,me_fft
 integer :: mgfft2
 integer :: mkmem,mkmem2,mpw,mpw2,mpw2_tmp,mu,natom,natom2,nbks1,nbks2
 integer :: nfft2,nimage2
 integer :: nkpt,nkpt2,nproc_fft,npsp,npsp2,nspden,nspinor,nsppol,nsppol2,nsym
 integer :: nsym2
 integer :: ntypat,ntypat2,nu,option,optorth,paral_fft,prtvol,rdwr
 integer :: restart,sppoldbl
 integer :: tread,unkg2,usepaw
 integer :: ncid_hdr_in,ncid_hdr_out,ncerr,ndim,nvar,natt,uid
 integer :: bravais2(11),intarr(1),ngfft(18),ngfft2(18)
 integer,allocatable :: indkk(:,:),istwfk(:),istwfk2(:),kg2(:,:),nband(:)
 integer,allocatable :: nband2(:),npwarr(:),npwarr2(:),symafm(:),symafmtmp(:)
 integer,allocatable :: symrel(:,:,:),symreltmp(:,:,:)
 real(dp) :: boxcutmin2,dksqmax,ecut,ecut2,ucvol2,zion_max2
 real(dp) :: acell(3),dprarr(1),gmet(3,3),gmet2(3,3)
 real(dp) :: gprimd(3,3),gprimd2(3,3),rmet2(3,3),rprimd(3,3),rprimd2(3,3)
 real(dp),allocatable :: cg2(:,:),doccde(:),eigen(:),kptns(:,:),occ(:),occ2(:)
 real(dp),allocatable :: tnonstmp(:,:),xred(:,:),zionpsp2(:)
 character(len=24) :: codename
 character(len=fnlen) :: f1,f2,file2
 character(len=strlen) :: string,string_raw
 character(len=30) :: token
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg
 type(dataset_type) :: dtset2
 type(hdr_type) :: hdr,hdr2
 type(bandstructure_type) :: bstruct
 type(pseudopotential_type) :: psps
 type(wffile_type) :: wffinp,wffout
 type(wvl_internal_type) :: wvl
 type(pawrhoij_type), allocatable :: pawrhoij(:)
 type(pawtab_type), allocatable :: pawtab(:,:)
 type(pspheader_type),allocatable :: pspheads(:)

! ************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

 call xmpi_init()

!Initialize MPI : one should write a separate routine -init_mpi_enreg-
!for doing that !!

!Default for sequential use
 call nullify_mpi_enreg(mpi_enreg)
 mpi_enreg%world_comm=0
 mpi_enreg%world_group=0
 mpi_enreg%me=0
 mpi_enreg%nproc=1
 mpi_enreg%nproc_atom=1
 mpi_enreg%num_group_fft = 0 ! in some cases not initialized but referenced in xdef_comm.F90
 mpi_enreg%paral_compil=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%mode_para="n"
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0
 mpi_enreg%paral_spin = 0
!MG080916 If we want to avoid MPI preprocessing options, %proc_distr should be always allocated and
!set to mpi_enreg%me. In such a way we can safely test its value inside loops parallelized over k-points
!For the time being, do not remove this line since it is needed in outkss.F90.
!nullify(mpi_enreg%proc_distrb)
!nullify(mpi_enreg%bandfft_kpt,mpi_enreg%tab_kpt_distrib)

!Initialize MPI
#if defined HAVE_MPI

 mpi_enreg%world_comm=xmpi_world
 mpi_enreg%world_group=MPI_GROUP_NULL
 call MPI_COMM_RANK(xmpi_world,mpi_enreg%me,ierr)
 call MPI_COMM_SIZE(xmpi_world,mpi_enreg%nproc,ierr)
!write(std_out,*)' abinit : nproc,me=',mpi_enreg%nproc,mpi_enreg%me
 mpi_enreg%paral_compil=1
#endif

!Signal MPI I/O compilation has been activated
#if defined HAVE_MPI_IO
 mpi_enreg%paral_compil_mpio=1
 if(mpi_enreg%paral_compil==0)then
   write(message,'(6a)') ch10,&
&   ' abinit : ERROR -',ch10,&
&   '  In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
&   '  Action : recompile your code with different CPP flags.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
#endif

!Initialize paral_compil_kpt, actually always equal to paral_compil
!(paral_compil_kpt should be suppressed after big cleaning)
 mpi_enreg%paral_compil_kpt=0
 if(mpi_enreg%paral_compil==1) mpi_enreg%paral_compil_kpt=1

!Additional setup specific to newsp
 mpi_enreg%paral_compil_kpt=0

!Other values of mpi_enreg are dataset dependent, and should NOT be initialized
!inside newsp.F90.


 codename='NEWSP '//repeat(' ',18)
 call herald(codename,abinit_version,std_out)
!YP: calling dump_config() makes tests fail => commented
!call dump_config()

!Set up for debugging and set up for orthonormalized output
 debug=1
 doorth=1

!************************************************************************
!Take care of wf input file

 write(std_out,'(a)' ) ' Enter name of unformatted wf data input file'
 read (05, '(a)' ) f1
 write(std_out,'(a,a)' ) ' wfinp=',trim(f1)

!DEBUG
!write(std_out,*)' newsp : enter '
!stop
!ENDDEBUG

!Check that old wf file exists
 call isfile(f1,'old')

 wffinp%accesswff = IO_MODE_FORTRAN
 wffout%accesswff = IO_MODE_FORTRAN

#if defined HAVE_TRIO_NETCDF
!test if file is a netcdf file
 ncerr = nf90_open(path=f1, mode=NF90_NOWRITE, ncid=ncid_hdr_in)
 if (ncerr == NF90_NOERR) then
!  NOTE: the choice could be made here to use any combination of
!  netcdf or binary input or output.
   wffinp%accesswff = IO_MODE_NETCDF
   wffout%accesswff = IO_MODE_NETCDF
!  here we modify wfinp to save the netCDF identifier in it
   write(std_out,*) ' newsp : open a netCDF file ', trim(f1), ncid_hdr_in
   wfinp = ncid_hdr_in
   ncerr = nf90_Inquire(ncid=ncid_hdr_in,nDimensions=ndim,nVariables=nvar,&
&   nAttributes=natt,unlimitedDimId=uid)
   call handle_ncerr(ncerr, " general Inquire ")
   write(std_out,*) 'newsp : input found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
 end if
#endif

 wffinp%unwff=wfinp

!Open wf file for input wf
 if (wffinp%accesswff /= IO_MODE_NETCDF) then
   open (unit=wfinp,file=f1,form='unformatted',status='old')
 end if

!Read the header from the file
 rdwr=1
 if (wffinp%accesswff/= IO_MODE_NETCDF) then
   call hdr_io(fform1,hdr,rdwr,wfinp)
#if defined HAVE_TRIO_NETCDF
 else if (wffinp%accesswff==IO_MODE_NETCDF) then
   call hdr_io_netcdf(fform1,hdr,rdwr,wfinp)
 else if (wffinp%accesswff == IO_MODE_ETSF) then
   write (std_out,*) "FIXME: ETSF I/O support in newsp"
#endif
 end if

!Echo the header of the file
 rdwr=4
 call hdr_io(fform1,hdr,rdwr,6)

!Format 1 or 2 are allowed for input wavefunction
 if ( (fform1+1)/2 /= (fform+1)/2 ) then
   write(message, '(a,a,a,a,i10,a,i10,a,a)' ) ch10,&
&   ' newsp: ERROR -',ch10,&
&   '  Input wf file fform=',fform1,' should equal coded fform=',fform,ch10,&
&   '  Action : check that wf file is correct.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 headform=hdr%headform
 bantot=hdr%bantot
 natom=hdr%natom
 ngfft(1:3)=hdr%ngfft(1:3)
 nkpt=hdr%nkpt
 nspden=hdr%nspden
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 nsym=hdr%nsym
 ntypat=hdr%ntypat
 npsp=ntypat
 acell(:)=one
 ecut=hdr%ecut_eff
 rprimd(:,:)=hdr%rprimd(:,:)

!Make sure nsppol does not exceed msppol
 if (nsppol>msppol) then
   write(message, '(a,a,a,a,i4,a,a,a)' ) ch10,&
&   ' newsp: ERROR -',ch10,&
&   '  Wf file nsppol=',nsppol,' exceeds msppol=2.',ch10,&
&   '  Action : the wf file must be corrupted, change it.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Make sure nsym does not exceed msym
 if (nsym>msym) then
   write(message, '(a,a,a,a,i8,a,i8,a,a)' ) ch10,&
&   ' newsp: BUG -',ch10,&
&   '  Wf file nsym=',nsym,' exceeds dimensioned msym=',msym,ch10,&
&   '  msym must be raised in newsp, and the code recompiled.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 ABI_ALLOCATE( istwfk,(nkpt))
 ABI_ALLOCATE(kptns,(3,nkpt))
 ABI_ALLOCATE(nband,(nkpt*msppol))
 ABI_ALLOCATE( npwarr,(nkpt))
 ABI_ALLOCATE(occ,(bantot))
 ABI_ALLOCATE( symafm,(msym))
 ABI_ALLOCATE(symrel,(3,3,msym))
 ABI_ALLOCATE( xred,(3,natom))

 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 npwarr(1:nkpt)=hdr%npwarr(1:nkpt)
 symrel(:,:,1:nsym)=hdr%symrel(:,:,1:nsym)
 istwfk(1:nkpt)=hdr%istwfk(1:nkpt)
 kptns(:,1:nkpt)=hdr%kptns(:,1:nkpt)
 occ(1:bantot)=hdr%occ(1:bantot)
 xred(:,:)=hdr%xred(:,:)

!Get mpw, as the maximum value of npwarr(:)
 mpw=maxval(npwarr(:))

!DEBUG
!write(std_out,*)' newsp : after read'
!write(std_out,*)' newsp : nsym=',nsym
!stop
!ENDDEBUG

!************************************************************************
!Take care of input file for new data set

 write(std_out,'(a)' ) ' Enter name of formatted input file for new data set'
 read (05, '(a)' ) file2
 write(std_out,'(a,a)' ) ' input=',trim(file2)
!Check that second input file exists
 call isfile(file2,'old')

!Read the data from input file,
!Really need: ceksp2, ecut2, istwfk2
!nband2, ngfft2, ntypat2, occ2,
!Note : all variables in the next three calls and the allocate statement
!have 2 as index, EXCEPT : abinit_version,iout,lenstr,string

!strlen from defs_basis module
 option=1
 call instrng (file2,lenstr,option,strlen,string)

!Copy original file, without change of case
 string_raw=string

!To make case-insensitive, map characters of string to upper case:
 call inupper(string(1:lenstr))

!Might import data from CML file(s)
 call importcml(lenstr,string_raw,string,strlen)

 ntypat2=1 ; marr=1
 token = 'ntypat'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) ntypat2=intarr(1)
!Check that ntypat2 is greater than 0
 if (ntypat2<=0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' newsp : ERROR -',ch10,&
&   '  Input ntypat must be > 0, but was ',ntypat2,ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify ntypat in the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!Note that the default npsp is ntypat
 npsp2=ntypat2 ; marr=1
 token = 'npsp'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) npsp2=intarr(1)
!Check that npsp is greater than 0
 if (npsp2<=0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' abinit : ERROR -',ch10,&
&   '  Input npsp must be > 0, but was ',npsp2,ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify npsp in the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

 ABI_ALLOCATE(zionpsp2,(npsp2))
 zionpsp2(1:npsp2)=zero
 zion_max2=zero

!In newsp, msym is a parameter.
 ABI_ALLOCATE(symafmtmp,(msym))
 ABI_ALLOCATE(symreltmp,(3,3,msym))
 ABI_ALLOCATE(tnonstmp,(3,msym))

 token = 'natom'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
!Might initialize natom from CML file
 if(tread==0)then
   token = '_natom'
   call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 end if
 if(tread==1)then
   natom2=intarr(1)
 else
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' abinit: ERROR -',ch10,&
&   '  Input natom must be defined, but was absent.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!Defaults of variables related to invars0
 dtset2%maxnsym=384
 dtset2%natom=natom2;dtset2%npsp=npsp2;dtset2%ntypat=ntypat2
 dtset2%useria=0;dtset2%userib=0
 dtset2%useric=0;dtset2%userid=0;dtset2%userie=0
 dtset2%userra=zero;dtset2%userrb=zero
 dtset2%userrc=zero;dtset2%userrd=zero;dtset2%userre=zero
 dtset2%usewvl = 0;dtset2%macro_uj=0
 dtset2%paral_kgb=0;dtset2%paral_rf=0;dtset2%ngroup_rf=1
 dtset2%istatr=49;dtset2%istatshft=1
 dtset2%nimage=1;dtset2%use_gpu_cuda=0

 ABI_ALLOCATE(dtset2%acell_orig,(3,1))
 ABI_ALLOCATE(dtset2%algalch,(ntypat2))
 ABI_ALLOCATE(dtset2%amu,(ntypat2))
 ABI_ALLOCATE(dtset2%corecs,(ntypat2))
 ABI_ALLOCATE(dtset2%densty,(ntypat2,4))
 ABI_ALLOCATE(dtset2%dynimage,(1))
 ABI_ALLOCATE(dtset2%iatfix,(3,natom2))
 ABI_ALLOCATE(dtset2%jpawu,(ntypat2))
 ABI_ALLOCATE(dtset2%kberry,(3,20))
 ABI_ALLOCATE(dtset2%lexexch,(ntypat2))
 ABI_ALLOCATE(dtset2%lpawu,(ntypat2))
 ABI_ALLOCATE(dtset2%mixalch,(npsp2,ntypat2))
 ABI_ALLOCATE(dtset2%pimass,(ntypat2))
 ABI_ALLOCATE(dtset2%ptcharge,(ntypat2))
 ABI_ALLOCATE(dtset2%quadmom,(ntypat2))
 ABI_ALLOCATE(dtset2%ratsph,(ntypat2))
 ABI_ALLOCATE(dtset2%rprim_orig,(3,3,1))
 ABI_ALLOCATE(dtset2%rprimd_orig,(3,3,1))
 ABI_ALLOCATE(dtset2%so_psp,(npsp2))
 ABI_ALLOCATE(dtset2%spinat,(3,natom2))
 ABI_ALLOCATE(dtset2%shiftk,(3,8))
 ABI_ALLOCATE(dtset2%typat,(natom2))
 ABI_ALLOCATE(dtset2%upawu,(ntypat2))
 ABI_ALLOCATE(dtset2%vel_orig,(3,natom2,1))
 ABI_ALLOCATE(dtset2%xred_orig,(3,natom2,1))
 ABI_ALLOCATE(dtset2%ziontypat,(ntypat2))
 ABI_ALLOCATE(dtset2%znucl,(npsp2))

!Defaults of variables related to invars1
 call indefo1(dtset2)

 dtset2%charge=zero
 dtset2%boxcutmin=two

!This one if not set in indefo1, but should likely be ...
 dtset2%nshiftk=1

!Note that this data is transferred directly from header of disk file
!It requires npsp=npsp2
 dtset2%znucl(:)=hdr%znuclpsp(:)

 jdtset=0

 symafmtmp(:)=1
 symreltmp(:,:,:)=0
 do ii=1,3
   symreltmp(ii,ii,:)=1
 end do
 tnonstmp(:,:)=0.0d0

 call invars1(bravais2,dtset2,iout,jdtset,lenstr,&
& mband_upper,mpi_enreg,msym,string,&
& symafmtmp,symreltmp,tnonstmp,zion_max2)

 nsppol2=dtset2%nsppol
 nimage2=dtset2%nimage
 nsym2=dtset2%nsym
 nkpt2=dtset2%nkpt
 ABI_ALLOCATE(istwfk2,(nkpt2))
 ABI_ALLOCATE(nband2,(nkpt2*nsppol2))
 ABI_ALLOCATE(dtset2%istwfk,(nkpt2))
 ABI_ALLOCATE(dtset2%occ_orig,(mband_upper*nkpt2*nsppol2))
 ABI_ALLOCATE(dtset2%kpt,(3,nkpt2))
 ABI_ALLOCATE(dtset2%kptns,(3,nkpt2))
 ABI_ALLOCATE(dtset2%wtk,(nkpt2))
 ABI_ALLOCATE(dtset2%nband,(nkpt2*nsppol2))
 ABI_ALLOCATE(dtset2%bdgw,(2,dtset2%nkptgw,nsppol2))
 ABI_ALLOCATE(dtset2%kptgw,(3,dtset2%nkptgw))
 if(associated(dtset2%dynimage))  then
   ABI_DEALLOCATE(dtset2%dynimage)
 end if
 ABI_ALLOCATE(dtset2%dynimage,(nimage2))
 ABI_ALLOCATE(dtset2%symafm,(nsym2))
 ABI_ALLOCATE(dtset2%symrel,(3,3,nsym2))
 ABI_ALLOCATE(dtset2%tnons,(3,nsym2))
 dtset2%symafm(:)=symafmtmp(1:nsym2)
 dtset2%symrel(:,:,:)=symreltmp(:,:,1:nsym2)
 dtset2%tnons(:,:)=tnonstmp(:,1:nsym2)

!Really need: ceksp2, ecut2
!ngfft2, ntypat2, occ2

!Here set defaults
 ceksp2=0
 dtset2%kptrlen=zero
 ecut2=-1.0d0
 dtset2%istwfk(:)=0
 dtset2%kpt(:,:)=zero
 dtset2%kptnrm=one
 dtset2%nband(:)=0
 dtset2%nspden=1
 dtset2%nspinor=1
 dtset2%occ_orig(:)=zero
 dtset2%wtk(:)=zero
 dtset2%dynimage(:)=1

!These default values should come from indefo.f
!except for iscf, that is set to 0 here

 dtset2%dilatmx=one
 dtset2%enunit=0
 dtset2%exchn2n3d=0
 dtset2%ionmov=0
 dtset2%intxc=0
 dtset2%iprcch=2
 dtset2%iprcel=0
 dtset2%iprcfc=0
 dtset2%irdwfk=0
 dtset2%iscf=0
 dtset2%isecur=0
 dtset2%ixc=1
 dtset2%nqpt=0
 dtset2%restartxf=0
 dtset2%optcell=0
 dtset2%irdwfq=0
 dtset2%ird1wf=0
 dtset2%irdddk=0
 dtset2%kptopt=0
 dtset2%chkexit=2
 dtset2%ikhxc=0
 dtset2%nbdbuf=0
 dtset2%localrdwf=1
 dtset2%nberry=1
 dtset2%bdberry(1)=0
 dtset2%bdberry(2)=0
 dtset2%bdberry(3)=0
 dtset2%bdberry(4)=0
 dtset2%delayperm=0
 dtset2%signperm=1
 dtset2%nbandkss=0
 dtset2%npwkss=-1
 dtset2%ntypalch=0
 dtset2%berryopt=0
 dtset2%wfoptalg=0
 dtset2%nbdblock=1
 dtset2%kssform=0
 dtset2%usedmatpu=0
 dtset2%usepaw=0
 dtset2%usepawu=0
 dtset2%usewvl=0
 dtset2%useylm=0
 dtset2%td_mexcit=0
 dtset2%prtdensph=0

 dtset2%pawujat=1
 dtset2%pawujrad=20.0_dp
 dtset2%pawujv=0.1_dp/Ha_eV

 dtset2%rfasr=0
 dtset2%rfatpol(1)=1
 dtset2%rfatpol(2)=1
 dtset2%rfdir(1)=0
 dtset2%rfdir(2)=0
 dtset2%rfdir(3)=0
 dtset2%rfddk=0
 dtset2%rfelfd=0
 dtset2%rfmeth=1
 dtset2%rfphon=0
 dtset2%rfstrs=0
 dtset2%rfuser=0
 dtset2%rf1elfd=0
 dtset2%rf1phon=0
 dtset2%rf2elfd=0
 dtset2%rf2phon=0
 dtset2%rf3elfd=0
 dtset2%rf3phon=0

 dtset2%qptn(:)=zero

!DEBUG
!write(std_out,*)' newsp : string '
!write(std_out,*)string(1:lenstr)
!ENDDEBUG
!DEBUG
!write(std_out,*)' newsp : before invars2'
!write(std_out,*)' nkpt2=',nkpt2
!stop
!write(std_out,*)' newsp : before invars2 '
!write(std_out,*)' newsp : dtset2%nband=',dtset2%nband
!ENDDEBUG

 usepaw=0
!Default values for FFT sequential execution
 paral_fft=0
 me_fft=0
 nproc_fft=1

 ABI_ALLOCATE(pspheads,(npsp2))
 call invars2(bravais2,dtset2,iout,jdtset,lenstr,mband_upper,msym,&
& npsp2,string,usepaw,zionpsp2)
 ABI_DEALLOCATE(pspheads)

 istwfk2(:)=dtset2%istwfk(:)
 ecut2=dtset2%ecut
 boxcutmin2=dtset2%boxcutmin
 ngfft2(:)=dtset2%ngfft(:)

!DEBUG
 write(std_out,*)' newsp : after invars2 '
 write(std_out,*)' newsp : nkpt2=',nkpt2
 write(std_out,*)' newsp : dtset2%nband=',dtset2%nband
!ENDDEBUG

 ABI_DEALLOCATE(symafmtmp)
 ABI_DEALLOCATE(symreltmp)
 ABI_DEALLOCATE(tnonstmp)

!Compute mgfft2,mpw2_tmp,nfft2 for this data set
 rprimd2(:,:)=dtset2%rprimd_orig(:,:,1)
 call metric(gmet2,gprimd2,-1,rmet2,rprimd2,ucvol2)
 call getng(boxcutmin2,ecut2,gmet2,me_fft,mgfft2,nfft2,ngfft2,nproc_fft,nsym2,1,paral_fft,dtset2%symrel,&
& use_gpu_cuda=dtset2%use_gpu_cuda)
 dtset2%ngfft(:)=ngfft2(:)
 mpi_enreg%nproc_fft=nproc_fft; ngfft(10)=nproc_fft
 mpi_enreg%me_fft=me_fft;ngfft(11)=me_fft
 mpi_enreg%fft_option_lob=1
 call getmpw(ecut2,dtset2%exchn2n3d,gmet2,istwfk2,dtset2%kptns,&
& mpi_enreg,mpw2_tmp,nkpt2)
 nband2(1:nkpt2*nsppol2)=dtset2%nband(1:nkpt2*nsppol2)
 mband2=maxval(nband2(1:nkpt2*nsppol2))
 dtset2%mband=mband2

 ABI_ALLOCATE(occ2,(mband2*nkpt2*nsppol2))
 occ2(1:mband2*nkpt2*nsppol2)=dtset2%occ_orig(1:mband2*nkpt2*nsppol2)

!NOTE: At this point, much more checking could be done to
!guarantee some consistency between the new input file and the
!old wf file (e.g. same number of atoms, etc.).  I am NOT doing
!this checking right now (23 Jun 1993).  Most of the data from
!the original wf file (e.g. natom) is being passed along to the newly
!created file.

!Check same choice of spin-polarization
!(But I think it should work ...)
 if (nsppol/=nsppol2) then
   write(message, '(a,a,a,a,i6,a,i6,a,a,a,a)' ) ch10,&
&   ' newsp: ERROR -',ch10,&
&   '  Wf disk file nsppol=',nsppol,&
&   '  not equal input file nsppol=',nsppol2,ch10,&
&   '  These must be same for newsp to run.',ch10,&
&   '  Action : correct one of these files.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!************************************************************************
!At last, need name of wf output file, then initialize its header.

 write(std_out,'(a)' ) ' Enter name of unformatted wf data output file'
 read(05, '(a)' ) f2
 write(std_out,'(a,a)' ) ' wfout=',trim(f2)
!Check that new file does NOT already exist
 call isfile(f2,'new')

 write(std_out,'(a,3i6)' ) ' wf    file ng',ngfft(1:3)
 write(std_out,'(a,3i6)' ) ' input file ng',ngfft2(1:3)

 write(std_out,'(a,1p,e12.4)' ) ' wf    file ecut',ecut
 write(std_out,'(a,1p,e12.4)' ) ' input file ecut',ecut2

 write(std_out,'(a,1p,3e15.7)' ) ' wf    file acell',acell(1:3)
 write(std_out,'(a,1p,3e15.7)' ) ' input file acell',dtset2%acell_orig(1:3,1)

!Compute reciprocal space metric gmet from original unit cell
 call matr3inv(rprimd,gprimd)
 do nu=1,3
   do mu=1,3
     gmet(mu,nu)=gprimd(1,mu)*gprimd(1,nu)+&
&     gprimd(2,mu)*gprimd(2,nu)+&
&     gprimd(3,mu)*gprimd(3,nu)
   end do
 end do

 sppoldbl=1
 if(minval(dtset2%symafm(:))==-1)then
   if(nsppol==1 .and. nsppol2==2)sppoldbl=2
 end if
 ABI_ALLOCATE(indkk,(nkpt2*sppoldbl,6))

!Find nearest old k point to each new k point
 call listkk(dksqmax,gmet,indkk,kptns,dtset2%kptns,nkpt,nkpt2,nsym2,sppoldbl,&
& dtset2%symafm,dtset2%symrel,1)

 if (debug>0 .and. sppoldbl==1) then
!  For new wf file, compute overall total number of bands, bantot
   write(std_out,*)' original bantot=',bantot
   bantot2=0
   do isppol=1,nsppol2
     do ikpt2=1,nkpt2
       bantot2=bantot2+nband2(indkk(ikpt2,1)+(isppol-1)*nkpt2)
     end do
   end do
   write(std_out,*)' final    bantot=',bantot2
 end if

!DEBUG
!write(std_out,*)' newsp : before kpgsph, nkpt2=',nkpt2
!stop
!ENDDEBUG

!Compute npw for each k point of new wf file, then compute actual mpw2
 ABI_ALLOCATE(kg2,(3,mpw2_tmp))
 ABI_ALLOCATE(npwarr2,(nkpt2))

 nullify(mpi_enreg%bandfft_kpt)
 mpi_enreg%flag_ind_kg_mpi_to_seq = 0;
 if( (mpi_enreg%mode_para         =='b') .and. &
& (mpi_enreg%paral_compil_mpio == 1)  .and. &
& (wffout%accesswff==IO_MODE_MPI .or. wffinp%accesswff==IO_MODE_MPI) ) then
   mpi_enreg%flag_ind_kg_mpi_to_seq = 1
   ABI_ALLOCATE(mpi_enreg%bandfft_kpt,(nkpt2))

   do ikpt=1,nkpt2
     nullify(mpi_enreg%bandfft_kpt(ikpt)%ind_kg_mpi_to_seq)
   end do
 end if
 do ikpt2=1,nkpt2
   call kpgsph(ecut2,dtset2%exchn2n3d,gmet2,0,ikpt2,istwfk2(ikpt2),kg2,dtset2%kptns(:,ikpt2),&
&   1,mpi_enreg,mpw2_tmp,npwarr2(ikpt2))
 end do
 mpw2=maxval(npwarr2(:))
 ABI_DEALLOCATE(kg2)
 ABI_ALLOCATE(kg2,(3,mpw2))

!DEBUG
!write(std_out,*)' newsp : after kpgsph '
!write(std_out,*)' nsppol,nkpt,nkpt2',nsppol,nkpt,nkpt2
!write(std_out,*)' nband(1:nkpt)',nband(1:nkpt)
!write(std_out,*)' nband2(1:nkpt2)',nband2(1:nkpt2)
!write(std_out,*)' occ(1:nband(1)*nkpt)',occ(1:nband(1)*nkpt)
!stop
!ENDDEBUG

!Copy over occupancies, band by band, for new k points
!(also gives value of bantot2)
!--NOTE that OLD occupancies are used here, NOT the presumably
!desired new occupancies as requested in the new input file

 bantot2=0
!Warning : should be adapted to sppoldbl=2
 do isppol=1,nsppol
   do ikpt2=1,nkpt2
!    Count over skipped bands from original set
     bant0=0
     if (indkk(ikpt2,1)>1) then
       do ikpt1=1,indkk(ikpt2,1)-1
         bant0=bant0+nband(ikpt1+(isppol-1)*nkpt)
       end do
     end if
!    Determine nband2 for each k, spin--this may not exceed
!    nband1 for the same spin and the associated k
!    nbks1 is the nband for the same spin and associated k
!    nbks2 is the requested nband from formatted input file
     nbks1=nband(indkk(ikpt2,1)+(isppol-1)*nkpt)
     nbks2=nband2(ikpt2+(isppol-1)*nkpt2)
!    If number of bands is being increased, print warning and
!    reset nband2 to only available number
!    WARNING : does not take into account nspinor, while it should
     if (nbks2>nbks1) then
       write(std_out,'(a,i2,i5,2i6,a,i6)' ) &
&       ' newsp: isppol,ikpt2,nband1,nband2=',&
&       isppol,ikpt2,nbks1,nbks2,'=> lower nband2 to',nbks1
       nband2(ikpt2+(isppol-1)*nkpt2)=nbks1
       nbks2=nbks1
     end if
!    Define occupancies for new wf file using appropriate
!    nband for each k,spin
     do iband=1,nbks2
       occ2(iband+bantot2)=occ(iband+bant0)
     end do
     bantot2=bantot2+nbks2
   end do
 end do

!After all checks out, write first part of new header for new wf file

!The following are from the original wf file, NOT from the new
!input file:
!abinit_version, fform
!natom,nsppol,ntypat
!occ(shifted)
!All data from psp headers (of course), and xred.
!The following are modified from original wf file:
!bantot (now reflects new set of k points),ngfft,nkpt,
!nband(shifted and possibly made smaller), ecut,rprimd,kptns,
!Note especially that no account is taken of any consistency in
!the use of symmetry.

!Set up the psps part of the header (many other infos form
!the psps datastructure, but only these matters here)
 psps%ntypat=ntypat
 psps%npsp=npsp
 psps%usepaw=0
 ABI_ALLOCATE(psps%pspcod,(npsp))
 ABI_ALLOCATE(psps%pspdat,(npsp))
 ABI_ALLOCATE(psps%pspso,(npsp))
 ABI_ALLOCATE(psps%pspxc,(npsp))
 ABI_ALLOCATE(psps%title,(npsp))
 ABI_ALLOCATE(psps%znuclpsp,(npsp))
 ABI_ALLOCATE(psps%znucltypat,(ntypat))
 ABI_ALLOCATE(psps%zionpsp,(npsp))
 do ipsp=1,npsp
   psps%pspcod(ipsp)=hdr%pspcod(ipsp)
   psps%pspdat(ipsp)=hdr%pspdat(ipsp)
   psps%pspso(ipsp)=hdr%pspso(ipsp)
   psps%pspxc(ipsp)=hdr%pspxc(ipsp)
   psps%title(ipsp)=hdr%title(ipsp)
   psps%znuclpsp(ipsp)=hdr%znuclpsp(ipsp)
   psps%zionpsp(ipsp)=hdr%zionpsp(ipsp)
 end do
 do itypat=1,ntypat
   psps%znucltypat(itypat)=hdr%znucltypat(itypat)
 end do

!Initialize band structure datatype
 ABI_ALLOCATE(doccde,(bantot2))
 ABI_ALLOCATE(eigen,(bantot2))
 doccde(:)=zero ; eigen(:)=zero


!this quantities are not initialized but are not needed here
 dtset2%nelect=zero
 dtset2%tphysel=zero
 dtset2%tsmear=0.04
 dtset2%occopt=1

 call bstruct_init(bantot2,bstruct,dtset2%nelect,doccde,eigen,dtset2%istwfk,dtset2%kptns,&
& nband2,nkpt2,npwarr2,nsppol2,dtset2%nspinor,dtset2%tphysel,dtset2%tsmear,dtset2%occopt,occ2,dtset2%wtk)

 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)

!DEBUG
!write(std_out,*)' newsp: before hdr_init'
!ENDDEBUG

!Set up the band structure information in hdr2
 gscase=0
 call hdr_init(bstruct,abinit_version,dtset2,hdr2,pawtab,gscase,psps,wvl)

!DEBUG
!write(std_out,*)' newsp: after hdr_init'
!ENDDEBUG


!Clean band structure datatype
 call bstruct_clean(bstruct)

!Update header, with evolving variables
!WARNING : The number of atom in the disk file and the input file
!cannot change.
 ABI_ALLOCATE(pawrhoij,(natom*usepaw))
 call hdr_update(bantot2,hdr%etot,hdr%fermie,hdr2,natom,&
& hdr%residm,rprimd2,occ2,mpi_enreg,pawrhoij,usepaw,xred)
 ABI_DEALLOCATE(pawrhoij)

#if defined HAVE_TRIO_NETCDF
 if(wffout%accesswff==IO_MODE_NETCDF) then
!  Create empty netCDF file
   ncerr = nf90_create(path=f2, cmode=NF90_CLOBBER, ncid=ncid_hdr_out)
   call handle_ncerr(ncerr," create netcdf wavefunction file")
   ncerr = nf90_close(ncid_hdr_out)
   call handle_ncerr(ncerr," close netcdf wavefunction file")
 else if (wffout%accesswff == IO_MODE_ETSF) then
   write (std_out,*) "FIXME: ETSF I/O support in newsp"
 end if
#endif

!Open wf file for output wf
 if (wffout%accesswff/=IO_MODE_NETCDF) then
   open (unit=wfout,file=f2,form='unformatted',status='unknown')
#if defined HAVE_TRIO_NETCDF
 else if (wffout%accesswff==IO_MODE_NETCDF) then
!  here we modify wfout to save the netCDF identifier in it
   ncerr = nf90_open(path=f2, mode=NF90_WRITE, ncid=ncid_hdr_out)
   call handle_ncerr(ncerr," newsp : open netcdf wavefunction file")
   write(std_out,*) ' newsp : open a netCDF file ', trim(f2), ncid_hdr_out
   wfout = ncid_hdr_out
 else if (wffout%accesswff == IO_MODE_ETSF) then
   write (std_out,*) "FIXME: ETSF I/O support in newsp"
#endif
 end if

!Only now can we initialize wffout%unwff, in the NetCDF case.
 wffout%unwff=wfout

 rdwr=2
 if (wffout%accesswff/=IO_MODE_NETCDF) then
   call hdr_io(fform,hdr2,rdwr,wfout)
#if defined HAVE_TRIO_NETCDF
 else if (wffout%accesswff==IO_MODE_NETCDF) then
   call hdr_io_netcdf(fform,hdr2,rdwr,wfout)

   call ini_wf_netcdf(mpw2,ncid_hdr_out,gscase)
   ncerr = nf90_Inquire(ncid=ncid_hdr_out,nDimensions=ndim,nVariables=nvar,&
&   nAttributes=natt,unlimitedDimId=uid)
   call handle_ncerr(ncerr, " general Inquire ")
   write(std_out,*) 'newsp : output found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
 else if (wffout%accesswff == IO_MODE_ETSF) then
   write (std_out,*) "FIXME: ETSF I/O support in newsp"
#endif
 end if

!Echo header of future disk file
 write(std_out,*)' Echo the header of the new file, ',trim(f2)

 rdwr=4
 call hdr_io(fform,hdr2,rdwr,6)

 call hdr_clean(hdr2)
 ABI_DEALLOCATE(psps%pspcod)
 ABI_DEALLOCATE(psps%pspdat)
 ABI_DEALLOCATE(psps%pspso)
 ABI_DEALLOCATE(psps%pspxc)
 ABI_DEALLOCATE(psps%title)
 ABI_DEALLOCATE(psps%znuclpsp)
 ABI_DEALLOCATE(psps%zionpsp)
 ABI_DEALLOCATE(psps%znucltypat)

!write(unit=wfout) hdr%etot
!Header data is now complete for new file
!For old file, one is ready to read the first wf block

!nband  tells how many bands are to be read (how many expected)
!nband2 tells how many are to be written

!************************************************************************

 if (debug>0) then
   write(std_out,*)'newsp: Calling newkpt with :'
   write(std_out,*)' mband2,mpw,mpw2,nkpt,nkpt2='
   write(std_out,'(6i5)' ) mband2,mpw,mpw2,nkpt,nkpt2
   write(std_out,*)' nspinor,nspinor2,nsppol,nsppol2='
   write(std_out,'(6i5)' ) nspinor,dtset2%nspinor,nsppol,nsppol2
   write(std_out,'(15i4)')istwfk(1:nkpt)
   write(std_out,'(15i4)')istwfk2(1:nkpt2)
 end if

 fill=0 ; formeig=0 ; ireadwf=1 ; mkmem=0 ; mkmem2=0
 prtvol=dtset2%prtvol ; mpi_enreg%paralbd=0 ; restart=2

 mband1=maxval(nband(1:nkpt*nsppol))
 mband_rd=min(mband1,(mband2/dtset2%nspinor)*nspinor)
 mcg2=max(mpw*nspinor*mband_rd,mpw2*mband2*dtset2%nspinor)
 ABI_ALLOCATE(cg2,(2,mcg2))
 ABI_ALLOCATE(eigen,(mband2*nkpt2*nsppol))

!DEBUG
 write(std_out,*)' newsp : before newkpt, stop'
!stop
!ENDDEBUG
 optorth=1
 call newkpt (ceksp2,cg2,debug,ecut,ecut2,ecut2,eigen,dtset2%exchn2n3d,fill,formeig,&
& gmet,gmet2,headform,indkk,iout,ireadwf,istwfk,istwfk2,kg2,kptns,dtset2%kptns,&
& mband2,mcg2,mkmem,mkmem2,mpi_enreg,mpi_enreg,mpw,mpw2,&
& nband,nband2,ngfft,nkpt,nkpt2,npwarr,npwarr2,nspinor,dtset2%nspinor,&
& nsppol,nsppol2,nsym2,occ2,optorth,prtvol,restart,rprimd2,sppoldbl,&
& dtset2%symrel,dtset2%tnons,unkg2,wffinp,wffout)

!-------------------------------------------------------------------

 ABI_DEALLOCATE(istwfk)
 ABI_DEALLOCATE(kptns)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(zionpsp2)
 ABI_DEALLOCATE(istwfk2)
 ABI_DEALLOCATE(nband2)
 ABI_DEALLOCATE(occ2)
 ABI_DEALLOCATE(indkk)
 ABI_DEALLOCATE(npwarr2)
 ABI_DEALLOCATE(kg2)
 ABI_DEALLOCATE(cg2)
 ABI_DEALLOCATE(eigen)

 ABI_DEALLOCATE(dtset2%acell_orig)
 ABI_DEALLOCATE(dtset2%algalch)
 ABI_DEALLOCATE(dtset2%amu)
 ABI_DEALLOCATE(dtset2%corecs)
 ABI_DEALLOCATE(dtset2%densty)
 ABI_DEALLOCATE(dtset2%iatfix)
 ABI_DEALLOCATE(dtset2%jpawu)
 ABI_DEALLOCATE(dtset2%kberry)
 ABI_DEALLOCATE(dtset2%lexexch)
 ABI_DEALLOCATE(dtset2%lpawu)
 ABI_DEALLOCATE(dtset2%mixalch)
 ABI_DEALLOCATE(dtset2%pimass)
 ABI_DEALLOCATE(dtset2%ptcharge)
 ABI_DEALLOCATE(dtset2%quadmom)
 ABI_DEALLOCATE(dtset2%ratsph)
 ABI_DEALLOCATE(dtset2%rprim_orig)
 ABI_DEALLOCATE(dtset2%rprimd_orig)
 ABI_DEALLOCATE(dtset2%so_psp)
 ABI_DEALLOCATE(dtset2%spinat)
 ABI_DEALLOCATE(dtset2%shiftk)
 ABI_DEALLOCATE(dtset2%typat)
 ABI_DEALLOCATE(dtset2%upawu)
 ABI_DEALLOCATE(dtset2%vel_orig)
 ABI_DEALLOCATE(dtset2%xred_orig)
 ABI_DEALLOCATE(dtset2%ziontypat)
 ABI_DEALLOCATE(dtset2%znucl)

 ABI_DEALLOCATE(dtset2%bdgw)
 ABI_DEALLOCATE(dtset2%kptgw)
 ABI_DEALLOCATE(dtset2%istwfk)
 ABI_DEALLOCATE(dtset2%kpt)
 ABI_DEALLOCATE(dtset2%kptns)
 ABI_DEALLOCATE(dtset2%symafm)
 ABI_DEALLOCATE(dtset2%symrel)
 ABI_DEALLOCATE(dtset2%tnons)
 ABI_DEALLOCATE(dtset2%wtk)
 ABI_DEALLOCATE(dtset2%occ_orig)
 ABI_DEALLOCATE(dtset2%nband)
 ABI_DEALLOCATE(dtset2%dynimage)

 if (associated(mpi_enreg%bandfft_kpt)) then
   do ikpt=1,nkpt2
     if (associated(mpi_enreg%bandfft_kpt(ikpt)%ind_kg_mpi_to_seq)) then
       ABI_DEALLOCATE(mpi_enreg%bandfft_kpt(ikpt)%ind_kg_mpi_to_seq)
       nullify(mpi_enreg%bandfft_kpt(ikpt)%ind_kg_mpi_to_seq)
     end if
   end do
   ABI_DEALLOCATE(mpi_enreg%bandfft_kpt)
   nullify(mpi_enreg%bandfft_kpt)
 end if

 call hdr_clean(hdr)

#if defined HAVE_TRIO_NETCDF
 if (wffinp%accesswff==IO_MODE_NETCDF) then
   ncerr = nf90_close (ncid=ncid_hdr_in)
   call handle_ncerr(ncerr," newsp : final close netcdf input wavefunction file")
 else if (wffinp%accesswff == IO_MODE_ETSF) then
   write (std_out,*) "FIXME: ETSF I/O support in newsp"
 end if
 if (wffout%accesswff==IO_MODE_NETCDF) then
   ncerr = nf90_close (ncid=ncid_hdr_out)
   call handle_ncerr(ncerr," newsp : final close netcdf output wavefunction file")
 else if (wffout%accesswff == IO_MODE_ETSF) then
   write (std_out,*) "FIXME: ETSF I/O support in newsp"
 end if
#endif

 write(std_out,*)' newsp:  program ends normally'

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program newsp
!!***
