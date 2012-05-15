!{\src2tex{textfont=tt}}
!!****p* ABINIT/anaddb
!! NAME
!! anaddb
!!
!! FUNCTION
!! Main routine for analysis of the interatomic force constants and associated
!! properties.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG,DCA,JCC,CL,XW)
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
!! WARNING
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,anaddb_dtset_clean,anaddb_dtset_nullify,asria_calc
!!      asria_corr,asrprs,chneu9,create_ddb_blk,destroy_ddb_blk,destroy_ifc
!!      destroy_mpi_enreg,destroy_phondos,destroy_primcell_ddb_info,diel9,dtchi
!!      dtech9,elast9,electrooptic,elphon,gtblk9,gtdyn9,herald,init9
!!      init_primcell_ddb_info,inprep8,instr9,instrng,inupper,invars9,isfile
!!      mkherm,mkifc9,mkphbs,mkphdos,mkrdim,mpi_comm_rank,mpi_comm_size
!!      nullify_ddb_blk,nullify_mpi_enreg,outvars9,phfrq3,piezo9,print_phondos
!!      print_phondos_debye,prtph3,ramansus,rdddb9,refineblk,relaxpol,symph3
!!      thm9,thmeig,timein,write_primcell_ddb_info,wrtout,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program anaddb

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_build_info
 use m_xmpi
 use m_errors
 use m_ifc
 use m_ddb_blk
 use m_phdos
 use m_primcell_ddb_info
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools,     only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anaddb'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_42_parser
 use interfaces_51_manage_mpi
 use interfaces_72_response
 use interfaces_77_ddb
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
! Set array dimensions
!  msym =maximum number of symmetry elements in space group
 integer :: msym
!Define input and output unit numbers (some are defined in defs_basis -all should be there ...):
! FIXME: these should not be reserved unit numbers!
 integer,parameter :: ddbun=2
!Define unit number for the files that can be analysed with band2eps
 integer,parameter :: udispl=19,ufreq=18
 integer :: iodyn
 integer :: dimekb,dims
 integer :: iatom,iblok,iblok_stress,ibloknl,idir,ii,index
 integer :: ierr,iphl2,lenstr,lmnmax
 integer :: mband,mtyp,mpert,msize,natom,nblok,nblok2
 integer :: nkpt,nph2l,nsym,ntypat
 integer :: nunit,occopt,option,rftyp
 integer :: usepaw,vrsddb
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 integer, allocatable :: symq(:,:,:)
 integer, allocatable :: symrec(:,:,:)
 integer, allocatable :: symrel(:,:,:)
 integer,allocatable :: d2flg(:),indsym(:,:),typat(:)
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: etotal,factor,tcpu,tcpui,twall
 real(dp) :: twalli,ucvol
 real(dp) :: acell(3),dielt(3,3),dielt_rlx(3,3),elast(6,6),epsinf(3,3),gmet(3,3),gprim(3,3)
 real(dp) :: pel(3)
 real(dp) :: piezo(6,3),qphnrm(3),qphon(3,3),rmet(3,3),rprim(3,3)
 real(dp) :: rprimd(3,3),strten(6)
 real(dp), allocatable :: tnons(:,:)
 real(dp),allocatable :: amu(:)
 real(dp),allocatable :: d2asr(:,:,:,:,:),d2cart(:,:),dchide(:,:,:)
 real(dp),allocatable :: dchidt(:,:,:,:),displ(:),dyewq0(:),eigval(:,:)
 real(dp),allocatable :: eigvec(:,:,:,:,:),fact_oscstr(:,:,:),instrain(:,:)
 real(dp),allocatable :: fred(:,:),lst(:),phfrq(:)
 real(dp),allocatable :: rcan(:,:),rsus(:,:,:),trans(:,:)
 real(dp),allocatable :: singular(:),uinvers(:,:), vtinvers(:,:)
 real(dp),allocatable :: xcart(:),xred(:,:),zeff(:,:,:),zion(:)
 real(dp),allocatable :: d2asr_res(:,:,:,:,:)
 character(len=24) :: codename
 character(len=strlen) :: string
 character(len=fnlen) :: filnam(7),elph_base_name
 character(len=fnlen) :: tmpfilename
 character(len=fnlen) :: scphon_filename
 character(len=500) :: message
 type(anaddb_dataset_type) :: anaddb_dtset
 type(phonon_dos_type) :: phonon_dos
 type(primcell_ddb_info) :: pcell
 ! here only since needed for call to rwwf
 type(MPI_type) :: mpi_enreg

 type(ifc_type) :: ifc_obj
 type(ddb_blk_type), pointer :: ddb_blk

!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()

!Initialize MPI : one should write a separate routine -init_mpi_enreg-
!for doing that !!
 call nullify_mpi_enreg(mpi_enreg)

!Default for sequential use
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

 ABI_ALLOCATE(ddb_blk,)
 call nullify_ddb_blk(ddb_blk)

!MG080916 If we want to avoid MPI preprocessing options, %proc_distr should be always allocated and
!set to mpi_enreg%me. In such a way we can safely test its value inside loops parallelized over k-points
!For the time being, do not remove this line since it is needed in outkss.F90.
!nullify(mpi_enreg%proc_distrb)
!nullify(mpi_enreg%bandfft_kpt,mpi_enreg%tab_kpt_distrib)

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
   write(message,'(3a)')&
&   '  In order to use MPI_IO, you must compile with the MPI flag ',ch10,&
&   '  Action : recompile your code with different CPP flags.'
   MSG_ERROR(message)
 end if
#endif

!Initialize paral_compil_kpt, actually always equal to paral_compil
!(paral_compil_kpt should be suppressed after big cleaning)
 mpi_enreg%paral_compil_kpt=0
 if(mpi_enreg%paral_compil==1) mpi_enreg%paral_compil_kpt=1

!Other values of mpi_enreg are dataset dependent, and should NOT be initialized
!inside anaddb.F90.

!Initialisation of the timing
 call timein(tcpui,twalli)

 codename='ANADDB'//repeat(' ',18)
 call herald(codename,abinit_version,std_out)
!YP: calling dump_config() makes tests fail => commented
!call dump_config()

 call anaddb_dtset_nullify(anaddb_dtset)

!Initialise the code : write heading, and read names of files.
 call init9(filnam)

!******************************************************************

 write(message, '(a,a,a,a)' )&
& ch10,ch10,' Read the input file',ch10
 call wrtout(std_out,message,'COLL')

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')

!Must read natom from the DDB before being able to allocate
!some arrays needed for invars9
 vrsddb=100401
 call inprep8(dimekb,filnam(3),lmnmax,mband,mtyp,msym,&
& natom,nblok,nkpt,ntypat,ddbun,usepaw,vrsddb)

!allocate arrays depending on msym (which is actually fixed to nsym inside inprep8)
 ABI_ALLOCATE(symq,(4,2,msym))
 ABI_ALLOCATE(symrec,(3,3,msym))
 ABI_ALLOCATE(symrel,(3,3,msym))
 ABI_ALLOCATE(tnons,(3,msym))

 mpert=natom+6
 msize=3*mpert*3*mpert
 if(mtyp==3)msize=msize*3*mpert

!Read the input file, and store the information in a long string of characters
!strlen from defs_basis module
 option=1
 call instrng (filnam(1),lenstr,option,strlen,string)

!To make case-insensitive, map characters to upper case:
 call inupper(string(1:lenstr))

 write(std_out,*) 'will read the inputs completely'

!Read the inputs
 call invars9 (anaddb_dtset,lenstr,natom,qtol,string)

 nph2l=anaddb_dtset%nph2l
 ABI_ALLOCATE(lst,(nph2l))

 write(std_out,*) 'read the inputs completely'

!Echo the inputs to console
 nunit=6
 call outvars9 (anaddb_dtset,nunit)

!Open output files iodyn and ab_out (might change its name if needed)
!MJV 1/2010 : now output file is open, but filnam(2) continues unmodified
!so the other output files are overwritten instead of accumulating.
 tmpfilename = filnam(2)
 call isfile(tmpfilename,'new')
 open (unit=ab_out,file=tmpfilename,form='formatted',status='new')
 rewind (unit=ab_out)
 call herald(codename,abinit_version,ab_out)
!MJV : standardize file name for DOS output
 if (anaddb_dtset%eivec==3) then
   tmpfilename = trim(filnam(2))//"_LWF"
   call isfile(tmpfilename,'new')
   iodyn = get_unit()
   open (unit=iodyn,file=tmpfilename,form='formatted',status='new')
   rewind (unit=iodyn)
 end if

!Echo the inputs to long printout
 call outvars9 (anaddb_dtset,ab_out)

!******************************************************************

!Read the DDB information, also perform some checks, and symmetrize partially the DDB

 write(message, '(a,a)' ) &
& ' read the DDB information and perform some checks',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a,a)' )&
& '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 call create_ddb_blk(msize, nblok, ddb_blk)

 ABI_ALLOCATE(amu,(ntypat))
 ABI_ALLOCATE(d2cart,(2,msize))
 ABI_ALLOCATE(indsym,(4,msym*natom))
 ABI_ALLOCATE(typat,(natom))
 ABI_ALLOCATE(xcart,(3*natom))
 ABI_ALLOCATE(xred,(3,natom))
 ABI_ALLOCATE(zion,(ntypat))
 ABI_ALLOCATE(instrain,(3*natom,6))

!DEBUG
 write(std_out,*)' anaddb : call rdddb9, read the file with name :'
 write(std_out,*)'  ',trim(filnam(3))
!ENDDEBUG

 call rdddb9(acell,anaddb_dtset%atifc,amu,ddb_blk,&
& ddbun,dimekb,filnam(3),gmet,gprim,indsym,ab_out,&
& lmnmax,mband,mpert,msize,msym,&
& anaddb_dtset%natifc,natom,nkpt,nsym,ntypat,&
& occopt,rmet,rprim,symq,symrec,symrel,&
& tnons,typat,ucvol,usepaw,xcart,xred,zion)

!Now the whole DDB is in central memory, contained in the
!array ddb_blk%val(2,msize,nblok).
!The information on it is contained in the four arrays
!ddb_blk%flg(msize,nblok) : blok flag for each element
!ddb_blk%qpt(9,nblok)  : blok wavevector (unnormalized)
!ddb_blk%nrm(3,nblok)  : blok wavevector normalization
!ddb_blk%typ(nblok)    : blok type

 if(anaddb_dtset%brav/=1 .and. abs(abs(rprim(1,2))-half)>tol10)then
!  Renormalize rprim to possibly satisfy the constraint abs(rprim(1,2))=half when brav/=1
   if(abs(rprim(1,2))<tol6)then
     write(message, '(a,i6,7a)' )&
&     '  The input DDB value of brav is',anaddb_dtset%brav,',',ch10,&
&     '  and the one of rprim(1,2) is zero.',ch10,&
&     '  These are incompatible',ch10,&
&     '  Action : check the value of brav and rprim(1,2) in your DDB.'
     MSG_ERROR(message)
   end if
   factor=abs(rprim(1,2))*two
   acell(:)=acell(:)*factor
   rprim(:,:)=rprim(:,:)/factor
   gprim(:,:)=gprim(:,:)*factor
 end if

 call mkrdim(acell,rprim,rprimd)

 ABI_ALLOCATE(displ,(2*3*natom*3*natom))
 ABI_ALLOCATE(dyewq0,(3*3*natom))
 ABI_ALLOCATE(d2asr,(2,3,natom,3,natom))
 ABI_ALLOCATE(eigval,(3,natom))
 ABI_ALLOCATE(eigvec,(2,3,natom,3,natom))
 ABI_ALLOCATE(phfrq,(3*natom))
 ABI_ALLOCATE(rcan,(3,natom))
 ABI_ALLOCATE(trans,(3,natom))
 ABI_ALLOCATE(zeff,(3,3,natom))

!**********************************************************************
!**********************************************************************

!Acoustic Sum Rule
 d2asr = zero

!In case the interatomic forces are not calculated, the
!ASR-correction (d2asr) has to be determined here from
!the Dynamical matrix at Gamma.
 if(anaddb_dtset%ifcflag==0 .or. &
& anaddb_dtset%instrflag/=0 .or. &
& anaddb_dtset%elaflag/=0)then

!  Find the Gamma block in the DDB (no need for E-field entries)
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0
   rfphon(1:2)=1
   rfelfd(:)=0
   rfstrs(:)=0
   rftyp=anaddb_dtset%rfmeth

   call gtblk9(ddb_blk,iblok,mpert,natom,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

   if(anaddb_dtset%asr==1 .or. anaddb_dtset%asr==2) then
     if (iblok /=0) then
       call asria_calc(anaddb_dtset%asr,d2asr,ddb_blk%val(:,:,iblok),mpert,natom)
     end if
   end if

!  Rotational invariance for 1D and 0D systems

   if(anaddb_dtset%asr==3 .or. anaddb_dtset%asr==4) then
     dims=3*natom*(3*natom-1)/2
     ABI_ALLOCATE(uinvers,(1:dims,1:dims))
     ABI_ALLOCATE(vtinvers,(1:dims,1:dims))
     ABI_ALLOCATE(singular,(1:dims))
     uinvers=0d0
     vtinvers=0d0
     singular=0d0
     if (iblok /= 0) then
       call asrprs(anaddb_dtset%asr,1,3,uinvers,vtinvers,singular,ddb_blk%val(:,:,iblok),mpert,natom,xcart)
     end if
   end if
   if(anaddb_dtset%asr==5) then
     if (iblok /=0 ) then
!      d2cart is a temp variable here
       d2cart = ddb_blk%val(:,:,iblok)
!      calculate diagonal correction
       call asria_calc(2,d2asr,d2cart,mpert,natom)
!      apply diagonal correction
       call asria_corr(2,d2asr,d2cart,mpert,natom)
!      hermitianize
       call mkherm(d2cart,3*mpert)
!      remove remaining ASR rupture due to Hermitianization
       ABI_ALLOCATE(d2asr_res,(2,3,natom,3,natom))
       call asria_calc(anaddb_dtset%asr,d2asr_res,d2cart,mpert,natom)
!      full correction is sum of both
       d2asr = d2asr + d2asr_res
       ABI_DEALLOCATE(d2asr_res)
     else
       d2asr = zero
     end if
   end if

 end if

!**********************************************************************

!Dielectric Tensor and Effective Charges

!Look for the Gamma Blok in the DDB
 qphon(:,1)=0.0d0
 qphnrm(1)=0.0d0
 rfphon(1:2)=1
 rfelfd(1:2)=2
 rfstrs(1:2)=0
 rftyp=anaddb_dtset%rfmeth

 call gtblk9(ddb_blk,iblok,mpert,natom,&
& qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

!Compute effective charges and dielectric tensor only if the
!Gamma-blok was found in the DDB or the occupation is metallic
!In case it was not found, iblok = 0

 zeff(:,:,:)=0.0d0
 dielt(:,:)=0.0d0
 dielt(1,1)=1.0d0 ; dielt(2,2)=1.0d0 ; dielt(3,3)=1.0d0

 if ((iblok/=0)) then

   write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Dielectric Tensor and Effective Charges ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   write(message, '(a,i6)' )' The Gamma block is : ',iblok
   call wrtout(std_out,message,'COLL')

!  Make the imaginary part of the Gamma block vanish
   write(message, '(a,a,a,a,a)'  ) ch10,&
&   ' anaddb : Zero the imaginary part of the Dynamical Matrix at Gamma,',ch10,&
&   '   and impose the ASR on the effective charges ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  Impose the charge neutrality on the effective charges
!  and eventually select some parts of the effective charges
   call chneu9(anaddb_dtset%chneut,ddb_blk%val(:,:,iblok),mpert,natom,ntypat,&
&   anaddb_dtset%selectz,typat,zion)

!  Extraction of the dielectric tensor and the effective charges
   call dtech9(ddb_blk%val,dielt,iblok,mpert,natom,ddb_blk%nblok,zeff)

 end if    ! iblok not found

!**********************************************************************

!Structural response at fixed polarization
 if (anaddb_dtset%polflag == 1) then

   ABI_ALLOCATE(d2flg,(msize))

   if(iblok/=0)then

!    Save the second-order derivatives
     d2cart(1:2,1:msize) = ddb_blk%val(1:2,1:msize,iblok)
     d2flg(1:msize) = ddb_blk%flg(1:msize,iblok)

   else ! the gamma blok has not been found

     if(anaddb_dtset%relaxat==0 .and. anaddb_dtset%relaxstr==0)then

!      The gamma blok is not needed
       d2cart(1:2,1:msize)=zero
       d2flg(1:msize)=1

     else ! There is a problem !

       write(message, '(7a)' )&
&       '  The dynamical matrix at Gamma is needed, in order to perform ',ch10,&
&       "  relaxation at constant polarisation (Na Sai's method)",ch10,&
&       '  However, this was not found in the DDB.',ch10,&
&       '  Action : complete your DDB with the dynamical matrix at Gamma.'
       MSG_ERROR(message)
     end if

   end if ! iblok not found

!  Read the block with the total energy
   qphon(:,:) = zero
   qphnrm(:) = zero
   rfphon(:) = 0
   rfelfd(:) = 0
   rfstrs(:) = 0
   rftyp = 0
   call gtblk9(ddb_blk,iblok,mpert,natom,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
   etotal = ddb_blk%val(1,1,iblok)

!  Read the block with the gradients
   ABI_ALLOCATE(fred,(3,natom))
   rftyp = 4
   rfelfd(:) = 2
   if (anaddb_dtset%relaxat == 1) rfphon(:) = 1
   if (anaddb_dtset%relaxstr == 1) rfstrs(:) = 3
   call gtblk9(ddb_blk,iblok,mpert,natom,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

   if (anaddb_dtset%relaxat == 1) then
     index = 0
     do iatom = 1, natom
       do idir = 1, 3
         index = index + 1
         fred(idir,iatom) = ddb_blk%val(1,index,iblok)
       end do
     end do
   end if

   pel(1:3) = ddb_blk%val(1,3*natom+4:3*natom+6,iblok)

   if (anaddb_dtset%relaxstr == 1) then
     index = 3*natom + 6
     do ii = 1, 6
       index = index + 1
       strten(ii) = ddb_blk%val(1,index,iblok)
     end do
   end if

   call relaxpol(d2flg,d2cart,etotal,fred,anaddb_dtset%iatfix,&
&   indsym,ab_out,anaddb_dtset%istrfix,&
&   mpert,msize,msym,anaddb_dtset%natfix,natom,&
&   anaddb_dtset%nstrfix,nsym,ntypat,pel,&
&   anaddb_dtset%relaxat,anaddb_dtset%relaxstr,&
&   rprimd,strten,symrel,anaddb_dtset%targetpol,typat,ucvol,xcart,xred,zion)

   ABI_DEALLOCATE(fred)
   ABI_DEALLOCATE(d2flg)

 end if

!***************************************************************************

!Compute non-linear optical susceptibilities and
!First-order change in the linear dielectric susceptibility
!induced by an atomic displacement

 if (anaddb_dtset%nlflag > 0) then

   qphon(:,:) = 0_dp
   qphnrm(:)  = 1_dp
   rfphon(1)  = 1 ; rfphon(2:3) = 0
   rfelfd(:)  = 2
   rfstrs(:)  = 0
   rftyp = 3

   call gtblk9(ddb_blk,iblok,mpert,natom,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

   ibloknl = iblok
   ABI_ALLOCATE(dchide,(3,3,3))
   ABI_ALLOCATE(dchidt,(natom,3,3,3))

   call dtchi(ddb_blk%val(:,:,ibloknl),dchide,dchidt,mpert,natom,anaddb_dtset%ramansr)

 end if ! nlflag

!**********************************************************************
!**********************************************************************

!Interatomic Forces Calculation

!DEBUG
!write(std_out,*)' anaddb : before ifcflag check, ifcflg=',anaddb_dtset%ifcflag,anaddb_dtset%thmflag
!stop
!ENDDEBUG
!if the ifc_obj is not calculated, we need to initialize it to 0 for call to mkphbs etc...
 if (anaddb_dtset%ifcflag ==0 ) then
   ifc_obj%nrpt = 0
   ABI_ALLOCATE(ifc_obj%atmfrc,(2,3,natom,3,natom,ifc_obj%nrpt))
   ABI_ALLOCATE(ifc_obj%rpt,(3,ifc_obj%nrpt))
   ABI_ALLOCATE(ifc_obj%wghatm,(natom,natom,ifc_obj%nrpt))
   ifc_obj%atmfrc = zero
   ifc_obj%rpt = zero
   ifc_obj%wghatm = zero

!  ifc to be calculated for interpolation
 else

   if (anaddb_dtset%qrefine > 1) then
!    if we are using refinement scheme, this modifies the content of ddb_blk%* arrays
     call refineblk(acell,amu,anaddb_dtset,ddb_blk,&
&     dielt,gmet,gprim,indsym,ab_out,&
&     mpert,msym,natom,nsym,ntypat,rmet,rprim,&
&     symrec,symrel,tcpui,twalli,typat,&
&     ucvol,xred,zeff)
   end if

   write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the interatomic forces ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call mkifc9(acell,amu,anaddb_dtset,ddb_blk,&
&   dielt,dyewq0,gmet,gprim,ifc_obj,indsym,ab_out,&
&   mpert,msym,natom,anaddb_dtset%ngqpt,nsym,ntypat,rcan,rmet,rprim,&
&   symrec,symrel,tcpui,trans,twalli,typat,&
&   ucvol,xred,zeff)


   if (anaddb_dtset%outscphon == 1) then
     call init_primcell_ddb_info (pcell,anaddb_dtset%brav,anaddb_dtset%dipdip,mpert,msym,natom,ifc_obj%nrpt,nsym,ntypat,ucvol,&
&     indsym,symrec,symrel,typat,&
&     acell,amu,dielt,dyewq0,gmet,gprim,rcan,rmet,rprim,ifc_obj%rpt,trans,ifc_obj%wghatm,xred,zeff)
     scphon_filename=trim(filnam(2))//"_PCINFO"
     scphon_filename=trim(scphon_filename)
     call write_primcell_ddb_info (scphon_filename,pcell)
     call destroy_primcell_ddb_info (pcell)
   end if

   call wrtout(std_out,' anaddb    : end of the IFC section ','COLL')
 end if

 write(std_out,*)' anaddb : after ifcflag check, ifcflg,thmflag,elphflag,prtdos='
 write(std_out,*)anaddb_dtset%ifcflag,anaddb_dtset%thmflag,anaddb_dtset%elphflag,anaddb_dtset%prtdos

!**********************************************************************
!**********************************************************************

!Short-Range/Long-Range decomposition of the phonon frequencies
!if (anaddb_dtset%prtsrlr == 1) then
!call wrtout(std_out,' anaddb    : start of the SR/LR decomposition ','COLL')
!end if



!**********************************************************************
!**********************************************************************

!Electron-phonon section
 if (anaddb_dtset%elphflag == 1) then
   call elphon(anaddb_dtset,filnam,acell,amu,ifc_obj%atmfrc,dielt,dyewq0,gmet,&
&   gprim,indsym,mpert,mpi_enreg,natom,ifc_obj%nrpt,nsym,ntypat,&
&   rcan,rmet,rprim,ifc_obj%rpt,symrec,symrel,tnons,trans,typat,ucvol,ifc_obj%wghatm,xred,zeff)
 end if


!**********************************************************************
!**********************************************************************

!Phonon density of states calculation, Start if interatomic forces have been calculated
 if (anaddb_dtset%ifcflag==1.and.(anaddb_dtset%prtdos==1.or.anaddb_dtset%prtdos==2)) then
   write(message,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of phonon density of states ',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   call mkphdos(phonon_dos,anaddb_dtset%prtdos,anaddb_dtset%dosdeltae,anaddb_dtset%dossmear,&
&   anaddb_dtset%dipdip,anaddb_dtset%symdynmat,&
&   acell,amu,anaddb_dtset,ifc_obj%atmfrc,dielt,dyewq0,gmet,gprim,indsym,&
&   mpert,msym,natom,ifc_obj%nrpt,nsym,ntypat,rmet,rprim,ifc_obj%rpt,symrec,symrel,&
&   trans,typat,ucvol,ifc_obj%wghatm,xred,zeff)

   call print_phondos(phonon_dos,"PHDOS")

   call print_phondos_debye(phonon_dos, ucvol)

   call destroy_phondos(phonon_dos)
 end if

!Phonon density of states and thermodynamical properties calculation
!Start if interatomic forces and thermal flags are on
 if(anaddb_dtset%ifcflag==1 .and. (anaddb_dtset%thmflag==1 .or. anaddb_dtset%thmflag==2)) then

   write(message, '(a,(80a),a,a,a,a,a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of phonon density of states, ',ch10,&
&   '    thermodynamical properties, ',ch10,&
&   '    and Debye-Waller factors.',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if(anaddb_dtset%thmflag==1) then
     call thm9(acell,amu,anaddb_dtset,ifc_obj%atmfrc,dielt,displ,dyewq0,d2cart,&
&     eigval,eigvec,gmet,gprim,indsym,ab_out,mpert,msym,natom,&
&     ifc_obj%nrpt,nsym,ntypat,filnam(2),phfrq,rmet,rprim,ifc_obj%rpt,symrec,symrel,tcpui,&
&     trans,twalli,typat,ucvol,ifc_obj%wghatm,xred,zeff)
   else if (anaddb_dtset%thmflag==2) then
     write(message, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&     ch10,' Entering thm9 routine with thmflag=2 ',ch10
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     call thm9(acell,amu,anaddb_dtset,ifc_obj%atmfrc,dielt,displ,dyewq0,d2cart,&
&     eigval,eigvec,gmet,gprim,indsym,ab_out,mpert,msym,natom,&
&     ifc_obj%nrpt,nsym,ntypat,filnam(2),phfrq,rmet,rprim,ifc_obj%rpt,symrec,symrel,tcpui,&
&     trans,twalli,typat,ucvol,ifc_obj%wghatm,xred,zeff, anaddb_dtset%thmflag, udispl, ufreq)

   end if

 end if

!**********************************************************************
!**********************************************************************

!Now treat the first list of vectors (without non-analyticities)
 call mkphbs(acell,amu,anaddb_dtset,ifc_obj%atmfrc,ddb_blk,&
& d2asr,dielt,dyewq0,filnam(2),gmet,gprim,indsym,iodyn,&
& mpert,msize,msym,natom,ifc_obj%nrpt,nsym,ntypat,&
& qtol,rmet,rprim,ifc_obj%rpt,singular,symrel,tcpui,  &
& trans,twalli,typat,ucvol,uinvers,vtinvers,ifc_obj%wghatm,xred,zeff)


!***********************************************************************
!***********************************************************************
!Test thmeig
 if(anaddb_dtset%thmflag>=3 .and. anaddb_dtset%thmflag<=8) then

!  DEBUG
   write(std_out,*)' anaddb : call rdddb9, read the second-order electron-phonon file with name :'
   write(std_out,*)'  ',trim(filnam(5))
!  ENDDEBUG

!  Obtain the number of bloks contained in this file.
   call inprep8(dimekb,filnam(5),lmnmax,mband,mtyp,msym,&
&   natom,nblok2,nkpt,ntypat,ddbun,usepaw,vrsddb)

!  DEBUG
!  write(std_out,*)' anaddb : nblok=',nblok
!  write(std_out,*)' anaddb : nblok2=',nblok2
!  ENDDEBUG

   write(std_out,*)'Entering thmeig: '
   elph_base_name=trim(filnam(2))//"_ep"
   call thmeig(anaddb_dtset%a2fsmear,acell,amu,anaddb_dtset,d2asr,&
&   elph_base_name,mband,mpert,msize,natom,nkpt,anaddb_dtset%ntemper,&
&   ntypat,rprim,anaddb_dtset%telphint,anaddb_dtset%temperinc,&
&   anaddb_dtset%tempermin,anaddb_dtset%thmflag,typat,xred,&
&   ddb_blk,ddbun,dimekb,filnam(5),ab_out,& !new
&  lmnmax,msym,nblok2,nsym,occopt,symrel,tnons,usepaw,zion,& !new
&  symrec,anaddb_dtset%natifc,gmet,gprim,indsym,rmet,anaddb_dtset%atifc,& !new
&  ucvol,xcart) !new
 end if

!**********************************************************************

!DEBUG
!write(std_out,*)' anaddb : nph2l,dieflag=',nph2l,anaddb_dtset%dieflag
!stop
!ENDDEBUG
!Now treat the second list of vectors (only at the Gamma point,
!but can include non-analyticities), as well as the
!frequency-dependent dielectric tensor

 if (anaddb_dtset%nlflag > 0)  then
   ABI_ALLOCATE(rsus,(3*natom,3,3))
 end if
 ABI_ALLOCATE(fact_oscstr,(2,3,3*natom))

 if( nph2l/=0 .or. anaddb_dtset%dieflag==1 )then

   write(message, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&   ch10,' Treat the second list of vectors ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  Before examining every direction or the dielectric tensor,
!  generates the dynamical matrix at gamma
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0

!  Generation of the dynamical matrix in cartesian coordinates
   if(anaddb_dtset%ifcflag==1)then

!    Get d2cart using the interatomic forces and the
!    long-range coulomb interaction through Ewald summation
     call gtdyn9(acell,ifc_obj%atmfrc,dielt,anaddb_dtset%dipdip,&
&     dyewq0,d2cart,gmet,gprim,mpert,natom,&
&     ifc_obj%nrpt,qphnrm(1),qphon,rmet,rprim,ifc_obj%rpt,&
&     trans,ucvol,ifc_obj%wghatm,xred,zeff)

   else if(anaddb_dtset%ifcflag==0)then

!    Look after the information in the DDB
     rfphon(1:2)=1
     rfelfd(1:2)=2
     rfstrs(1:2)=0
     rftyp=anaddb_dtset%rfmeth
     call gtblk9(ddb_blk,iblok,mpert,natom,&
&     qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

!    Copy the dynamical matrix in d2cart
     d2cart(:,1:msize)=ddb_blk%val(:,:,iblok)

!    Eventually impose the acoustic sum rule
     call asria_corr(anaddb_dtset%asr,d2asr,d2cart,mpert,natom)

!    end of the generation of the dynamical matrix at gamma.
   end if

   if(nph2l/=0)then

!    Examine every wavevector of this list
     do iphl2=1,nph2l

!      Initialisation of the phonon wavevector
       qphon(:,1)=anaddb_dtset%qph2l(:,iphl2)
       qphnrm(1)=anaddb_dtset%qnrml2(iphl2)

!      Calculation of the eigenvectors and eigenvalues
!      of the dynamical matrix
       call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&       mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,rprimd,anaddb_dtset%symdynmat,&
&       symrel,typat,ucvol)

!      Write the phonon frequencies
       call prtph3(displ,anaddb_dtset%eivec,anaddb_dtset%enunit,-1,ab_out,natom,phfrq,qphnrm(1),qphon)

!      Determine the symmetries of the phonon modes at Gamma
       if(sum(abs(qphon(:,1)))<qtol)then
         call symph3(ab_out,acell,eigvec,indsym,natom,nsym,phfrq,rprim,symrel)
       end if

!      Write Raman susceptibilities
       if (anaddb_dtset%nlflag == 1) then
         call ramansus(d2cart,dchide,dchidt,displ,mpert,&
&         natom,phfrq,qphon,qphnrm(1),rsus,ucvol)
       end if

!      Prepare the evaluation of the Lyddane-Sachs-Teller relation
       if(anaddb_dtset%dieflag==1 .and. natom>1)then
         lst(iphl2)=zero
!        The fourth mode should have positive frequency, otherwise,
!        there is an instability, and the LST relationship should not
!        be evaluated
         if(phfrq(4)>tol6)then
           do ii=4,3*natom
             lst(iphl2)=lst(iphl2)+2*log(phfrq(ii))
           end do
         end if
       end if

     end do ! iphl2
   end if ! nph2l/=0

!  The frequency-dependent dielectric tensor (and oscillator strength).
   if (anaddb_dtset%dieflag==1)then

     write(message, '(a,a,a,a,a,a)' )&
&     ' anaddb : the frequency-dependent dielectric tensor (and also once more',&
&     ch10,' the phonons at gamma - without non-analytic part )',ch10,&
&     ch10,' The frequency-dependent dielectric tensor'
     call wrtout(std_out,message,'COLL')

!    Initialisation of the phonon wavevector
     qphon(:,1)=0.0d0
     qphnrm(1)=0.0d0

!    Calculation of the eigenvectors and eigenvalues
!    of the dynamical matrix
     call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&     mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,&
&     rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol)

!    Write the phonon frequencies (not to ab_out, however)
     call prtph3(displ,0,anaddb_dtset%enunit,-1,-1,natom,phfrq,qphnrm(1),qphon)

!    Evaluation of the oscillator strengths and frequency-dependent
!    dielectric tensor.
     call diel9(amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
&     ab_out,lst,mpert,natom,nph2l,ntypat,phfrq,qtol,typat,ucvol)

!    DEBUG
!    write(std_out,*)' anaddb : after diel9, dielt_rlx(:,:)=',dielt_rlx(:,:)
!    ENDDEBUG

   end if

!  If the electronic dielectric tensor only is needed...
   if (anaddb_dtset%dieflag==2.or.anaddb_dtset%dieflag==3&
&   .or. anaddb_dtset%dieflag==4)then

!    Everything is already in place...
     call diel9(amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
&     ab_out,lst,mpert,natom,nph2l,ntypat,phfrq,qtol,typat,ucvol)

   end if

!  End the condition of either nph2l/=0  or  dieflag==1
 end if

!**********************************************************************

!In case nph2l was equal to 0, the electronic dielectric tensor
!has to be computed independently.

 if( anaddb_dtset%dieflag==2 .and. anaddb_dtset%nph2l==0 )then

   write(message, '(a)' )&
&   ' anaddb : nph2l=0, so compute the electronic dielectric tensor independently'
   call wrtout(std_out,message,'COLL')

!  Look after the second derivative matrix at gamma in the DDB
!  Note that the information on the dielectric tensor is completely
!  independent of the interatomic force constant calculation
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0
   rfphon(1:2)=0
   rfelfd(1:2)=2
   rfstrs(1:2)=0
   rftyp=anaddb_dtset%rfmeth
   call gtblk9(ddb_blk,iblok,mpert,natom,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

   d2cart(:,1:msize)=ddb_blk%val(:,:,iblok)

!  Print the electronic dielectric tensor
   call diel9(amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
&   ab_out,lst,mpert,natom,nph2l,ntypat,phfrq,qtol,typat,ucvol)

!  DEBUG
!  write(std_out,*)' anaddb : after third diel9, dielt_rlx(:,:)=',dielt_rlx(:,:)
!  ENDDEBUG

 end if

!**********************************************************************

!Compute the electrooptic tensor

 if (anaddb_dtset%nlflag == 1) then

!  In case dieflag = 2, recompute phonon frequencies and
!  eigenvectors without non-analyticity
   if (anaddb_dtset%dieflag == 2) then
     qphon(:,1)=0.0d0
     qphnrm(1)=0.0d0
     call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&     mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,&
&     rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol)
   end if

   rsus(:,:,:) = 0_dp
   call ramansus(d2cart,dchide,dchidt,displ,mpert,&
&   natom,phfrq(1),qphon,qphnrm(1),rsus,ucvol)

   call electrooptic(dchide,anaddb_dtset%dieflag,epsinf,&
&   fact_oscstr,natom,phfrq,anaddb_dtset%prtmbm,rsus,ucvol)

 end if  ! condition on nlflag and dieflag

 ABI_DEALLOCATE(fact_oscstr)
 if (anaddb_dtset%nlflag > 0)  then
   ABI_DEALLOCATE(dchide)
   ABI_DEALLOCATE(dchidt)
   ABI_DEALLOCATE(rsus)
 end if

!**********************************************************************

!here treating the internal strain tensors at Gamma point
 if(anaddb_dtset%instrflag/=0)then

   write(message, '(a,a,(80a),a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the internal-strain  tensor',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message,'(a,f11.3,a,f11.3,a)')&
&   '-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'

   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  allocate(instrain(3*natom,6))
   if(anaddb_dtset%instrflag==1)then
     write(message,'(a)' )&
&     ' anaddb : instrflag=1, so extract the internal strain constant from the 2DTE'
     call wrtout(std_out,message,'COLL')

!    looking after the no. of blok that contians internal strain tensor
     qphon(:,1)=0.0d0
     qphnrm(1)=0.0d0
     rfphon(1:2)=0
     rfelfd(1:2)=0
     rfstrs(1:2)=3
     rftyp=anaddb_dtset%rfmeth
     call gtblk9(ddb_blk,iblok,mpert,natom,&
&     qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
!    then print the internal stain tensor
!    write(ab_out,'(/,a,i6)')'iblok is',iblok
     call instr9(anaddb_dtset%asr,ddb_blk%val,d2asr,iblok,instrain,ab_out,mpert,natom,ddb_blk%nblok)
   end if
 end if
!end the part for internal strain

!**********************************************************************

!here treating the elastic tensors at Gamma Point
 if(anaddb_dtset%elaflag/=0)then
   write(message, '(a,a,(80a),a,a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the elastic and compliances tensor (Voigt notation)',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message,'(a,f11.3,a,f11.3,a)')&
&   '-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'

   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')


   if(anaddb_dtset%elaflag==1 .or.anaddb_dtset%elaflag==2&
&   .or. anaddb_dtset%elaflag==3 .or.anaddb_dtset%elaflag==4&
&   .or. anaddb_dtset%elaflag==5)then
     write(message,'(a)' )&
&     ' anaddb : so extract the elastic constant from the 2DTE'
     call wrtout(std_out,message,'COLL')

!    look after the blok no. that contains the stress tensor
     qphon(:,1)=0.0d0
     qphnrm(1)=0.0d0
     rfphon(1:2)=0
     rfelfd(1:2)=0
     rfstrs(1:2)=0
     rftyp=4
     call gtblk9(ddb_blk,iblok,mpert,natom,&
&     qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
     iblok_stress=iblok

!    DEBUG
!    check the iblok number containing first order derivative
!    write(std_out,'(/,a,/)')'iblok_stress number'
!    write(std_out,'(i)')iblok_stress
!    ENDDEBUG

!    look after the blok no.iblok that contains the elastic tensor
     qphon(:,1)=0.0d0
     qphnrm(1)=0.0d0
     rfphon(1:2)=0
     rfelfd(1:2)=0
     rfstrs(1:2)=3
!    for both diagonal and shear parts
     rftyp=anaddb_dtset%rfmeth
     call gtblk9(ddb_blk,iblok,mpert,natom,&
&     qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
!    write(ab_out,'(/,a,i6)')'iblok is',iblok
!    print the elastic tensor
     call elast9(anaddb_dtset,ddb_blk%val,d2asr,elast,iblok,iblok_stress,instrain,ab_out,mpert,&
&     natom,ddb_blk%nblok,ucvol)
   end if
 end if
!ending the part for elastic tensors

!**********************************************************************

!here treating the piezoelectric tensor at Gamma Point
 if(anaddb_dtset%piezoflag/=0 .or. anaddb_dtset%dieflag==4&
& .or. anaddb_dtset%elaflag==4)then
   write(message, '(a,a,(80a),a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the tensor related to piezoelectric effetc',ch10,&
&   '  (Elastic indices in Voigt notation)',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   call timein(tcpu,twall)
   write(message,'(a,f11.3,a,f11.3,a)')&
&   '-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'

   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if(anaddb_dtset%piezoflag==1 .or.anaddb_dtset%piezoflag==2&
&   .or.anaddb_dtset%piezoflag==3 .or. anaddb_dtset%piezoflag==4&
&   .or.anaddb_dtset%piezoflag==5 .or. anaddb_dtset%piezoflag==6&
&   .or.anaddb_dtset%piezoflag==7 .or. anaddb_dtset%dieflag==4&
&   .or.anaddb_dtset%elaflag==4)then
     write(message,'(a)' )&
&     ' anaddb : extract the piezoelectric constant from the 2DTE'
     call wrtout(std_out,message,'COLL')
!    looking for the gamma point block
     qphon(:,1)=0.0d0
     qphnrm(1)=0.0d0
     rfphon(1:2)=0
     rfelfd(1:2)=0
     rfstrs(1:2)=3
!    for both diagonal and shear parts
     rftyp=anaddb_dtset%rfmeth

     call gtblk9(ddb_blk,iblok,mpert,natom,&
&     qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
!    write(ab_out,'(/,a,i6)')'iblok is',iblok
!    then print out the piezoelectric constants

!    DEBUG
!    write(std_out,*)' anaddb : before piezo9, dielt_rlx(:,:)=',dielt_rlx(:,:)
!    ENDDEBUG

     call piezo9(anaddb_dtset,ddb_blk%val,dielt_rlx,elast,iblok,instrain,ab_out,mpert,&
&     natom,ddb_blk%nblok,piezo,ucvol)
   end if
 end if

!**********************************************************************

 call anaddb_dtset_clean(anaddb_dtset)
 ABI_DEALLOCATE(symq)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(tnons)

 ABI_DEALLOCATE(amu)
 ABI_DEALLOCATE(displ)
 ABI_DEALLOCATE(dyewq0)
 ABI_DEALLOCATE(d2asr)
 ABI_DEALLOCATE(d2cart)
 ABI_DEALLOCATE(eigval)
 ABI_DEALLOCATE(eigvec)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(instrain)
 ABI_DEALLOCATE(lst)
 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(rcan)
 ABI_DEALLOCATE(trans)
 ABI_DEALLOCATE(typat)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(zeff)
 ABI_DEALLOCATE(zion)

 call destroy_ifc(ifc_obj)

 call destroy_ddb_blk(ddb_blk)
 ABI_DEALLOCATE(ddb_blk)

 call timein(tcpu,twall)
 write(message, '(a,(80a),a,a,a,f11.3,a,f11.3,a,a,a,a)' ) ch10,&
& ('=',ii=1,80),ch10,ch10,&
& '+Total cpu time',tcpu-tcpui,&
& '  and wall time',twall-twalli,' sec',ch10,ch10,&
& ' anaddb : the run completed succesfully.'
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 close(ab_out)

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program anaddb
!!***


