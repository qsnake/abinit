!{\src2tex{textfont=tt}}
!!****p* ABINIT/lwf
!! NAME
!! lwf
!!
!! FUNCTION
!! Main routine for the generation of the Lattice Wannier Functions
!! from disentangled phonon bands.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group
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
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,atmdata,bldlwf,chkilwf,destroy_mpi_enreg,instrng,intagm
!!      inupper,invars7w,leave_new,mpi_comm_rank,mpi_comm_size
!!      nullify_mpi_enreg,overlap_ph,readeig,rwwan,secinit,shellin,wanvec
!!      wrtout,xmpi_end,xmpi_init,zhpev,zmnbld
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program lwf

 use defs_basis
 use defs_abitypes
 use m_xmpi
#if defined HAVE_MPI2
 use mpi
#endif

!no_abirules
! Global Variables:
!! allerr= maximum allowed mismatch between initial quess and the frozen states
!!         it is expressed as ratio from the frozen window width
!! alpha= real parameter for the kmixing of Zmn_2 and Zmn_1
!! atmass(natom)= atomic masses, finally in sqrt form
!! decflg= flag for forced stop of the minimzation of Omega_I at its minimum
!! eigval(nqpt,3*natom)=array containing the eigenvalues, in cm-1
!! eigvect(nqpt,3*natom,natom,3,2)=array containing the eigenvectors:
!!  its components are: q-point, band, atom, displacement, re/im part
!! eps(nqpt,nwnn,maxqsize)= eigenvalues of the Zmn matrix
!! enwdmax,enwdmin = global energy window boundaries
!! f_,g_,z_ subsp = indexes of the bands in the free,global,frozen window
!! frozflg = flag for frozen window: 1 use; 0 no use of frozen states
!! g_, f_, z_ subspace (nqpt,3*natom)= index of the Global Free froZen bands
!! grdsize(3) = size of the grid of q points = limit of the shells of the LWF
!! ingss(nwnn,natom,3,2)=initial guess of the W states, user-input
!! indeps(nqpt,nwnn)= index for the maximum value of the Zmn eigenvalues
!! indwnz(nqpt,nwnn)= index of the W states in the frozen window
!! irwfl (+1/-1) = flag for reading previously generated Wannier bands
!!     +1 = write the lambdas
!!     -1 = reads the lambdas
!! lambda(nqpt,nwnn,maxqsize,2)= overlaps (real) over all states of the initial LWF
!! maxeps= (real) maximum value of the Z eigenvalues
!! maxqsize = maximum number of bands in the global window
!! mmnkb(nqpt,6,maxqsize,maxqsize,2)= overlaps between the global states
!! mqpt= number of maximally allowed q points (parameter)
!! natom= number of atoms per unit cell
!! nqpt= number of q points in the whole BZ
!! nstom= number of steps for the minimization of Omega_I
!! nwnn=no of wannier functions to be generated
!! nwnz(nqpt)= number of Wannier states in the frozen window
!! prtvol = debugging flag
!! qneigh(nqpt,6)= q neighbours
!! qpoint(nqpt,3)= list of the q points in fractional coordinates
!! qsize(nqpt,3)= number of the bands in the energy window : G - Z - F
!! rcenter(3)= center of the Wannier functions
!! subwdmax,subwdmin = frozen energy window boundaries
!! tcpui, twalli = for the time routine
!! tolomi= real, TOLerance in OMegaI
!! ut1(nqpt,nwnn,maxqsize,maxqsize,2)= eigenvectors of the Zmn matrix
!! wannval(nqpt,nwnn)= "eigenvalue" of the LWF
!! wannvect(nqpt,nwnn,natom,3,2)= "eigenvector" of the LWF (=displacement)
!! xred(natom,3)= atomic reduced coordinates
!! znucl(natom)= nuclear Z, real(dp)
!! Zmn,Zmn1,Zmn2(nqpt,maxqsize,maxqsize,2)=Z matrix

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lwf'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_parser
 use interfaces_51_manage_mpi
 use interfaces_77_lwf
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
! Set array dimensions
!  mqpt =maximum number of q wavevectors in the list to be computed
 integer,parameter :: mqpt=100000
!Define input and output unit numbers:
 integer,parameter :: ineig=3,iout=7,iwf=8!,ininp=4
!Integer scalar
 integer :: aa,contwr,decflg,frozflg,gband,i1,i2,iatom,iband,ierr,ii,iqpt,irwfl
 integer :: iwann,jdtset,lenstr,marr,maxqsize,natom,nqpt,nstom!,mnatom,mnqpt
 integer :: nwnn,option,prtvol,step,tao,tread,trialq,zband!,oklocalmode
!Integer arrays
 integer :: grdsize(3)
 integer,allocatable :: f_subsp(:,:),g_subsp(:,:),indeps(:,:),indwnz(:,:)
 integer,allocatable :: nwnz(:),qneigh(:,:),qsize(:,:),typat(:),z_subsp(:,:)
!Real scalars
 real(dp) :: allerr,alpha,enwdmax,enwdmin,omega
 real(dp) :: rcov,subwdmax,subwdmin,testnorm,tolomi
!Real arrays
 real(dp) :: acell(3),localqmode(3),rcenter(3),rprim(3,3)
 real(dp),allocatable :: Omega_I(:),Zmn(:,:,:,:),Zmn1(:,:,:,:),Zmn2(:,:,:,:)
 real(dp),allocatable :: atmass(:),eigval(:,:),eigvect(:,:,:,:,:),eigwan(:)
 real(dp),allocatable :: eps1(:),ingss(:,:,:,:),lambda(:,:,:,:),matrx(:,:)
 real(dp),allocatable :: mmnkb(:,:,:,:,:),qpoint(:,:),ut1(:,:,:),wannval(:,:)
 real(dp),allocatable :: wannvect(:,:,:,:,:),xred(:,:),zhpev1(:,:),zhpev2(:)
 real(dp),allocatable :: znucl(:)
!Dummy arguments for subroutine 'intagm' to parse input file
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)
!Characters
 character(len=5) :: symbol
 character(len=30) :: token
 character(len=500) :: message
 character(len=fnlen) :: filnam(4)
 character(len=strlen) :: string
!Other variables
 type(MPI_type) :: mpi_enreg

!******************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

!Initialize MPI : one should write a separate routine -init_mpi_enreg-
!for doing that !!
 call xmpi_init()

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

!Other values of mpi_enreg are dataset dependent, and should NOT be initialized
!inside lwf.F90.


!!Initialisation of the timing
!! call timein(tcpui,twalli)

!!Read the file names
 write(std_out,*)' Give name for      formatted input file : '
 read(05, '(a)' ) filnam(1)
 write(std_out,'(1x,a)' ) trim(filnam(1))
 write(std_out,*)' Give name for     formatted dynamical input file : '
 read(05, '(a)' ) filnam(2)
 write(std_out,'(1x,a)' ) trim(filnam(2))
 write(std_out,*)' Give name for output file : '
 read(05, '(a)' ) filnam(3)
 write(std_out,'(1x,a)' ) trim(filnam(3))
 write(std_out,*)' Give name for Wannier eigenvector output file : '
 read(05, '(a)' ) filnam(4)
 write(std_out,'(1x,a)' ) trim(filnam(4))

!!open the files
!open (unit=ininp,file=filnam(1),form='formatted',status='old')  ! This file is opened in instrng
 open (unit=ineig,file=filnam(2),form='formatted',status='old')
 open (unit=iout,file=filnam(3),form='formatted')
 open (unit=iwf,file=filnam(4),form='formatted')

!******************************************************************
!!read the parameter input file
!first we have to get natom, nqpt and nwnn for the correct allocation of further arrays
!Compute the maximum size of arrays intarr and dprarr
 marr=4*mqpt
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))
 option=1
 call instrng (filnam(1),lenstr,option,strlen,string)

!To make case-insensitive, map characters to upper case:
 call inupper(string(1:lenstr))

 natom=0
 jdtset=1
 token = 'natom'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) natom=intarr(1)
 if(natom<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' lwf-main : ERROR -',ch10,&
&   '  natom is',natom,', which is lower than 1 .',ch10,&
&   '  Action : correct natom in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 nqpt=0
 jdtset=1
 token = 'nqpt'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) nqpt=intarr(1)
 if(natom<0)then
   write(message, '(a,a,a,i8,a,a,a)' )&
&   ' lwf-main : ERROR -',ch10,&
&   '  natom is',nqpt,', which is lower than 1 .',ch10,&
&   '  Action : correct nqpt in your input file.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 nwnn=0
 token = 'nwnn'
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) nwnn=intarr(1)
!check the existence of Wannier states
 if (nwnn<1) then
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' lwf  : ERROR -',ch10,&
&   '  The number of Wannier states is zero.',ch10,&
&   '  This is not allowed.  ',ch10,&
&   '  Action : modify the number of Wannier states in the input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!DEBUG
!write(std_out,*)' lwf main: before allocations'
!write(std_out,*)' natom=',natom
!write(std_out,*)' nqpt =',nqpt
!write(std_out,*)' nwnn =',nwnn
!mnatom=natom
!mnqpt=nqpt
!ENDDEBUG
!**************************************************************************
!allocation of the main variables

 ABI_ALLOCATE(atmass,(natom))
 ABI_ALLOCATE(eigval,(nqpt,3*natom))
 ABI_ALLOCATE(eigvect,(nqpt,3*natom,natom,3,2))
 ABI_ALLOCATE(eigwan,(nwnn))
 ABI_ALLOCATE(f_subsp,(nqpt,3*natom))
 ABI_ALLOCATE(g_subsp,(nqpt,3*natom))
 ABI_ALLOCATE(indeps,(nqpt,nwnn))
 ABI_ALLOCATE(indwnz,(nqpt,nwnn))
 ABI_ALLOCATE(ingss,(nwnn,natom,3,2))
 ABI_ALLOCATE(nwnz,(nqpt))
 ABI_ALLOCATE(qneigh,(nqpt,6))
 ABI_ALLOCATE(qpoint,(nqpt,3))
 ABI_ALLOCATE(qsize,(nqpt,3))
 ABI_ALLOCATE(typat,(natom))
 ABI_ALLOCATE(wannval,(nqpt,nwnn))
 ABI_ALLOCATE(wannvect,(nqpt,nwnn,natom,3,2))
 ABI_ALLOCATE(xred,(natom,3))
 ABI_ALLOCATE(z_subsp,(nqpt,3*natom))
 ABI_ALLOCATE(znucl,(natom))

!DEBUG
!write(std_out,*) ' lwf main: before invars'
!ENDDEBUG

!********************************************************************
!!read completely the input file
 call invars7w(allerr,alpha,decflg,enwdmax,enwdmin,frozflg,grdsize,ingss,irwfl,lenstr,&
& localqmode,mqpt,natom,nstom,nwnn,prtvol,rcenter,&
& string,subwdmax,subwdmin,tolomi,trialq,znucl)

!checks the input variables
 call chkilwf(alpha,decflg,enwdmax,enwdmin,frozflg,grdsize,irwfl,&
& nqpt,subwdmax,subwdmin)


 if (prtvol==-1) then     !after read completely the input file
   write(std_out,*) 'STOP requested by user in lwf.f after reading the input file'
!  stop
   call leave_new('COLL')
 end if

 ABI_ALLOCATE(Omega_I,(nstom+1))

!********************************************************************
!get atomic masses
 do tao=1,natom
   call atmdata(atmass(tao),rcov,symbol,znucl(tao))
   atmass(tao)=dble(sqrt(atmass(tao)*1822.88851))
 end do

 do iwann=1,nwnn               !loop over wannier states
   do tao=1,natom               !loop over atoms
     do aa=1,3                   !loop over displacements
       ingss(iwann,tao,aa,1)=atmass(tao)*ingss(iwann,tao,aa,1)
     end do                       !displacements
   end do                        !atoms
 end do                         !wannier states
!end tests for input

!********************************************************************
!!read the phonon input file
 call readeig(acell,eigval,eigvect,ineig,natom,nqpt,qpoint,rprim,typat,xred)

 if (prtvol==-2) then     !after read the phonon input file
   write(std_out,*) 'STOP requested by user in lwf.f after reading the phonon input file'
   stop
 end if

!DEBUG
!write(std_out,*)' lwf main: after readeig'
!write(std_out,*)'   nqpt=',nqpt
!ENDDEBUG

 if (irwfl==1) then   ! the eigenvectors of the Wannier interpolated states are not available
!  *************************************************
!  !compute the subspaces at each q point
!  delimitation of the bands in the windows
!  checks for consistency
   maxqsize=0
   do iqpt=1,nqpt                       ! loop over q points
!    DEBUG
!    write(std_out,*) 'debug, qpoint no.:',iqpt
!    ENDDEBUG

     gband=0
     zband=0
     qsize(iqpt,1)=0
     qsize(iqpt,2)=0
     qsize(iqpt,3)=0
     f_subsp(iqpt,:)=0
     g_subsp(iqpt,:)=0
     z_subsp(iqpt,:)=0
     nwnz(iqpt)=0
     do iatom=1,3*natom                  ! loop over all bands
!      renormalize displacement with sqrt of atomic masses
       do tao=1,natom                     ! loop over all atoms
         do aa=1,3                         ! loop over displacements, x,y,z
           eigvect(iqpt,iatom,tao,aa,1)=atmass(tao)*eigvect(iqpt,iatom,tao,aa,1)
           eigvect(iqpt,iatom,tao,aa,2)=atmass(tao)*eigvect(iqpt,iatom,tao,aa,2)
         end do                             !displacements
       end do                              !atoms

       if (eigval(iqpt,iatom).ge.enwdmin .and. eigval(iqpt,iatom).le.enwdmax) then
         gband=gband+1
         g_subsp(iqpt,gband)=iatom
         qsize(iqpt,1)=gband
         if (frozflg==1) then              ! frozen window state is defined
           if (eigval(iqpt,iatom).ge.subwdmin .and. eigval(iqpt,iatom).le.subwdmax) then
             zband=zband+1
             if (zband>nwnn) then            ! check the max. no. of frozen bands
               write(message, '(a,a,a,a,a,a,a,a,3f9.5,a,a,a,a)' ) ch10,&
&               ' lwf  : ERROR -',ch10,&
&               '  The number of frozen-states is larger ',ch10,&
&               '  than the number of Wannier functions ',ch10,&
&               '  at the q-point: ',(qpoint(iqpt,ii),ii=1,3),ch10,&
&               '  This is not allowed.  ',ch10,&
&               '  Action : modify the frozen-states energy window in the input file.'
               call wrtout(std_out,  message,'COLL')
               call leave_new('COLL')
             end if
             z_subsp(iqpt,gband)=iatom
             nwnz(iqpt)=nwnz(iqpt)+1
             qsize(iqpt,2)=qsize(iqpt,2)+1
           else
             f_subsp(iqpt,gband)=iatom
             qsize(iqpt,3)=qsize(iqpt,3)+1
           end if
         else                              ! frozen window state is not defined
           nwnz(iqpt)=0
           qsize(iqpt,2)=0
           qsize(iqpt,3)=qsize(iqpt,1)
           f_subsp(iqpt,gband)=iatom
         end if
       end if
     end do                                !all bands
     if (gband>maxqsize) maxqsize=gband

!    DEBUG
!    write(std_out,'(a,a,3f10.6,a,a,i4,a,a,i4,a,a,i4,a)') ch10,' number of bands at q point: ',qpoint(iqpt,:),ch10,&
!    ' global states = ',qsize(iqpt,1),ch10,&
!    ' frozen states = ',qsize(iqpt,2),ch10,&
!    ' free   states = ',qsize(iqpt,3)
!    write(std_out,*) 'index of bands: global - free - frozen'
!    write(std_out,'(30i5)') g_subsp(iqpt,1:maxqsize)
!    write(std_out,'(30i5)') f_subsp(iqpt,1:maxqsize)
!    write(std_out,'(30i5)') z_subsp(iqpt,1:maxqsize)
!    ENDDEBUG

     if (qsize(iqpt,1)<nwnn) then
       write(message, '(a,a,a,a,a,a,a,a,3f9.5,a,a,a,a)' ) ch10,&
&       ' lwf  : ERROR -',ch10,&
&       '  The number of global-states is smaller ',ch10,&
&       '  than the number of Wannier functions ',ch10,&
&       '  at the q-point: ',(qpoint(iqpt,ii),ii=1,3),ch10,&
&       '  This is not allowed.  ',ch10,&
&       '  Action : modify the global-states energy window in the input file.'
       call wrtout(std_out,  message,'COLL')
       call leave_new('COLL')
     end if

   end do                                !q points

   if (prtvol==-3) then   !after assignement of the number of bands at each q point
     write(std_out,*) 'STOP requested by user in lwf.f after assignement of the bands'
     stop    !after overlap
   end if

!  DEBUG
!  write(std_out,*) 'debug, lwf main: before setup the shells of each q point'
!  ENDDEBUG

!  **************************************************
!  setup the shells of each q point
   call shellin(acell,nqpt,qneigh,qpoint,rprim)

   if (prtvol==-4) then   !after assignement of the shells of each q point
     write(std_out,*) 'STOP requested by user in lwf.f after assignement of the Q+B shells'
     stop    !after overlap
   end if

!  DEBUG
!  write(std_out,*) 'debug, lwf main: before compute initial trial subspace'
!  ENDDEBUG


!  !compute initial trial subspace
!  first step : read
!  if trialq = 1 then the initial guess is a phonon
!  else initial guess is a user defined trial function
!  and it was read from the input file

!  if (trialq==1) then            ! initial quess is a phonon
!  oklocalmode=0
!  do iqpt=1,nqpt                ! loop over q points
!  if (abs(qpoint(iqpt,1)-localqmode(1))+abs(qpoint(iqpt,2)-localqmode(2))+&
!  &abs(qpoint(iqpt,3)-localqmode(3)) < 0.0000001 ) then
!  if (qsize(iqpt)<nwnn) then
!  write(message, '(a,a,a,a,i12,a,a,a,a,i12,a,a,a)' ) ch10,&
!  &  ' lwf  : ERROR -',ch10,&
!  &  '  Input local mode must contain all the ',nwnn,ch10,&
!  &  ' bands appearing in the frozen window, while it contains ',ch10,&
!  &  ' only ',qsize(iqpt),'  This is not allowed.  ',ch10,&
!  &  '  Action : modify local mode or energy windows in the input file.'
!  call wrtout(std_out,  message,'COLL')
!  call leave_new('COLL')
!  else
!  oklocalmode=1
!  do ii=1,nwnn                ! loop over W states
!  ingss(ii,:,:,:)=eigvect(iqpt,z_subsp(iqpt,ii),:,:,:)
!  end do                       !W states
!  exit
!  end if
!  end if
!  end do                          !q points
!  
!  if (oklocalmode==0) then
!  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
!  &  ' lwf  : ERROR -',ch10,&
!  &  '  Input local mode is not contained in the q-points grid.',ch10,&
!  &  '  This is not allowed.  ',ch10,&
!  &  '  Action : modify local mode in the input file.'
!  call wrtout(std_out,  message,'COLL')
!  call leave_new('COLL')
!  end if
!  
!  DEBUG
!  else
!  write(std_out,*) 'Initial guess is user-defined'
!  ENDDEBUG

!  end if                         !initial guess is completely determined here

!  DEBUG
!  write(std_out,*) 'debug, lwf main: before second step'
!  write(std_out,*) ' lwf : start to compute initial guess at different q points'
!  ENDDEBUG

!  *********************************************************************
!  second step : compute initial guess at different q points
   ABI_ALLOCATE(lambda,(nqpt,nwnn,maxqsize,2))

!  NOTE: do loop over q points in the lwf instead of secinit
!  for the sake of input in secinit simiplity + reduce some of
!  the memory needs

   do iqpt=1,nqpt
     call secinit(eigval,eigvect,f_subsp,g_subsp,z_subsp,&
&     qsize(iqpt,1),ingss,iqpt,lambda,maxqsize,natom,nqpt,nwnn,nwnz,qpoint)

!    DEBUG
!    write(std_out,*) 'after secinit. lambda is:'
!    do aa=1,nwnn
!    do ii=1,qsize(iqpt,1)
!    write(std_out,'(a,2f12.7)') ' LAMBDA ',lambda(iqpt,aa,ii,1),lambda(iqpt,aa,ii,2)
!    end do
!    end do
!    ENDDEBUG

   end do

   if (prtvol==-5) then        !after secinit
     write(std_out,*) 'STOP requested by user in lwf.f after secinit'
     stop    !after overlap
   end if

!  !main loop for minimizing Omega_I
!  ***************************************************************************
!  ! build the overlaps between neighbours
   ABI_ALLOCATE(mmnkb,(nqpt,6,maxqsize,maxqsize,2))
   call overlap_ph(eigvect,g_subsp,maxqsize,mmnkb,natom,nqpt,qneigh,qsize)

   if (prtvol==-6) then   !after overlap
     write(std_out,*) 'STOP requested by user in lwf.f after overlap'
     stop    !after overlap
   end if

!  *************************************************************************
!  ! build Zmn
   ABI_ALLOCATE(Zmn1,(nqpt,maxqsize,maxqsize,2))
   ABI_ALLOCATE(Zmn2,(nqpt,maxqsize,maxqsize,2))
   ABI_ALLOCATE(Zmn,(nqpt,maxqsize,maxqsize,2))

   Zmn=zero
   Zmn1=zero
   Zmn2=zero

!  *****
!  compute Z_0
   call zmnbld(f_subsp,lambda,maxqsize,mmnkb,natom,nqpt,nwnn,nwnz,qneigh,qsize,Zmn1)

   if (prtvol==-7) then       !after compute Z_0
     write(std_out,*) 'STOP requested by user in lwf.f after computation of Z_0'
     stop    !after overlap
   end if


!  ***** deal with the first step

   step=1
   Omega_I(:)=zero
   omega=zero
   do iqpt=1,nqpt               ! loop over q points

!    write(std_out,*)
!    write(std_out,*)
!    write(std_out,*) 'current Q point for the diagonalization of Zmn1',iqpt,'qsize, G - Z - F',qsize(iqpt,:)

!    ***** diagonalize Zmn1(iqpt)
     if (nwnz(iqpt)<nwnn) then    ! if there are non-frozen states
       ierr=0
       ii=1
       ABI_ALLOCATE(eps1,(qsize(iqpt,3)))
       ABI_ALLOCATE(ut1,(2,qsize(iqpt,3),qsize(iqpt,3)))
       ut1=zero
       eps1=zero
       ABI_ALLOCATE(matrx,(2,qsize(iqpt,3)*(qsize(iqpt,3)+1)/2))
       do i2=1,qsize(iqpt,3)
         do i1=1,i2
           matrx(1,ii)=Zmn1(iqpt,i2,i1,1)
           matrx(2,ii)=Zmn1(iqpt,i2,i1,2)
           ii=ii+1
         end do
       end do
       ABI_ALLOCATE(zhpev1,(2,2*qsize(iqpt,3)-1))
       ABI_ALLOCATE(zhpev2,(3*qsize(iqpt,3)-2))
       call zhpev('V','U',qsize(iqpt,3),matrx,eps1,ut1,qsize(iqpt,3),zhpev1,zhpev2,ierr)
       ABI_DEALLOCATE(matrx)
       ABI_DEALLOCATE(zhpev1)
       ABI_DEALLOCATE(zhpev2)


!      write(std_out,*) 'number of frozen wannier states',nwnz(iqpt)
!      write(std_out,'(a,i4,3f8.4,a,i4)') 'eigenvalues Zmn1 at q point',iqpt,qpoint(iqpt,:),' f_bands ',qsize(iqpt,3)
!      do ii=1,qsize(iqpt,3)
!      write(std_out,'(a,f12.8,a,i4)') 'eps1 eigenval=',eps1(ii),' at q point no',iqpt
!      do iband=1,qsize(iqpt,3)
!      write(std_out,'(a,3f22.18)') '     eivec =',ut1(1,iband,ii),ut1(2,iband,ii),ut1(1,iband,ii)*ut1(1,iband,ii)+ut1(2,iband,ii)*ut1(2,iband,ii)
!      end do
!      end do

!      ***** assign lambdas
       do ii=nwnz(iqpt)+1,nwnn              ! loop over free Wannier states
         omega=omega+eps1(qsize(iqpt,3)-ii+1+nwnz(iqpt))
         testnorm=zero
         aa=0
!        write(std_out,*) 'ii = ',ii,' corresponding eps',eps1(qsize(iqpt,3)-ii+1+nwnz(iqpt))
!        write(std_out,*) 'qsize(iqpt,1)',qsize(iqpt,1)
         do iband=1,qsize(iqpt,1)
           if (f_subsp(iqpt,iband)>0) then
             aa=aa+1
             lambda(iqpt,ii,iband,1)=ut1(1,aa,qsize(iqpt,3)-ii+1+nwnz(iqpt))
             lambda(iqpt,ii,iband,2)=ut1(2,aa,qsize(iqpt,3)-ii+1+nwnz(iqpt))

             testnorm=testnorm+lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1) &
&             +lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2)
           else
             lambda(iqpt,ii,iband,1)=zero
             lambda(iqpt,ii,iband,2)=zero
           end if
!          write(std_out,'(a,i4,a,3f12.8)') 'aa=',aa,' lambda',lambda(iqpt,ii,iband,1),lambda(iqpt,ii,iband,2),&
!          lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1)+lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2)
         end do                               !iband
         if (abs(testnorm-1.0d0)>1.0d-8) write(std_out,'(a,f12.8)') 'due to zmnbld, ERROR with lambda norm',testnorm
       end do                                !free Wannier states

       ABI_DEALLOCATE(eps1)
       ABI_DEALLOCATE(ut1)
     end if                        !if there are non-frozen states

!    ***** write lambdas
!    DEBUG
!    write(std_out,*) ' after updating lambda:'
!    do ii=1,nwnn
!    write(std_out,*) 'current lambdas for the Wannier state no ',ii
!    do iband=1,qsize(iqpt,1)
!    write(std_out,'(a,3f12.8)') 'lambda',lambda(iqpt,ii,iband,1),lambda(iqpt,ii,iband,2),&
!    lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1)+lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2)
!    end do
!    end do
!    ENDDEBUG

   end do                        !q points

!  Omega_I(1)=nwnn/6-omega
   Omega_I(1)=-omega
   write(std_out,*) 'at step',step,'Omega_I(step) is',Omega_I(step)
   write(iout,'(a,i5,a,f16.13)') 'at step',step,'Omega_I(step) is',Omega_I(step)

!  ***** update wannier information
   wannval=zero
   do iqpt=1,nqpt
!    do ii=nwnz(iqpt)+1,nwnn              ! loop over Wannier states
     do ii=1,nwnn              ! loop over Wannier states
       do iband=1,qsize(iqpt,1)            ! loop over global bands
         wannval(iqpt,ii)=wannval(iqpt,ii)+(lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1)+&
&         lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2))*eigval(iqpt,g_subsp(iqpt,iband))
       end do                               !global bands
     end do                                !Wannier states
   end do
   do iqpt=1,nqpt
!    write(std_out,*) 'qpoint, Wanneig',qpoint(iqpt,:),wannval(iqpt,:)
     write(std_out,'(i5,3f10.5,a,3f14.7)') iqpt,qpoint(iqpt,:),' 1st Omega, Wannier eigenvalue ',wannval(iqpt,:)

   end do

   if (prtvol==-8) then     !after first Omega step
     write(std_out,*) 'STOP requested by user in lwf.f after computation of first Omega step'
     stop    !after overlap
   end if


!  ******************************************************************************
!  iterative procedure of minimizing Omega_I
   contwr=0
   do
     contwr=contwr+1
     step=step+1
!    write(std_out,*) 'current step',step
     omega=0
     wannval=zero

!    ***** get Zmn_2
!    Zmn=a*Zmn_2+(1-a)*Zmn_1
     call zmnbld(f_subsp,lambda,maxqsize,mmnkb,natom,nqpt,nwnn,nwnz,qneigh,qsize,Zmn2)
     Zmn(:,:,:,:)=alpha*Zmn2(:,:,:,:)+(1-alpha)*Zmn1(:,:,:,:)

     do iqpt=1,nqpt               ! loop over q points
!      write(std_out,*) 'current Q point for the diagonalization of Zmn1',iqpt
!      write(std_out,*) 'qsize, G - Z - F',qsize(iqpt,:)
!      write(std_out,*) 'at current point',iqpt,' there are',nwnz(iqpt),' frozen LWFs'


!      ***** diagonalize Zmn(iqpt)
       if (nwnz(iqpt)<nwnn) then   !diagonalize if there are non-frozen LWF
         ierr=0
         ii=1
         ABI_ALLOCATE(eps1,(qsize(iqpt,3)))
         ABI_ALLOCATE(ut1,(2,qsize(iqpt,3),qsize(iqpt,3)))
         ut1=zero
         eps1=zero
         ABI_ALLOCATE(matrx,(2,qsize(iqpt,3)*(qsize(iqpt,3)+1)/2))
         ut1=zero
         do i2=1,qsize(iqpt,3)
           do i1=1,i2
             matrx(1,ii)=Zmn(iqpt,i2,i1,1)
             matrx(2,ii)=Zmn(iqpt,i2,i1,2)
             ii=ii+1
           end do
         end do
         ABI_ALLOCATE(zhpev1,(2,2*qsize(iqpt,3)-1))
         ABI_ALLOCATE(zhpev2,(3*qsize(iqpt,3)-2))
         call zhpev('V','U',qsize(iqpt,3),matrx,eps1,ut1,qsize(iqpt,3),zhpev1,zhpev2,ierr)
         ABI_DEALLOCATE(matrx)
         ABI_DEALLOCATE(zhpev1)
         ABI_DEALLOCATE(zhpev2)

!        ***** assign lambdas
         do ii=nwnz(iqpt)+1,nwnn              ! loop over free Wannier states
           omega=omega+eps1(qsize(iqpt,3)-ii+1+nwnz(iqpt))
           testnorm=zero
           aa=0
!          write(std_out,*) 'ii = ',ii,' corresponding eps',eps1(qsize(iqpt,3)-ii+1+nwnz(iqpt))
!          write(std_out,*) 'qsize(iqpt,1)',qsize(iqpt,1)
           do iband=1,qsize(iqpt,1)
!            write(std_out,*) 'iband',iband,f_subsp(iqpt,iband)
             if (f_subsp(iqpt,iband)>0) then
               aa=aa+1
!              write(std_out,'(a,4i5)') 'in loop: aa,qsize(3),nwnz,iterator:',aa,qsize(iqpt,3),nwnz(iqpt),qsize(iqpt,3)-ii+1+nwnz(iqpt)
               lambda(iqpt,ii,iband,1)=ut1(1,aa,qsize(iqpt,3)-ii+1+nwnz(iqpt))
               lambda(iqpt,ii,iband,2)=ut1(2,aa,qsize(iqpt,3)-ii+1+nwnz(iqpt))
!              write(std_out,'(a,i4,a,3f12.8)') 'in loop, aa=',aa,' lambda',lambda(iqpt,ii,iband,1),lambda(iqpt,ii,iband,2),&
!              lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1)+lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2)
               testnorm=testnorm+lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1)+lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2)
             else
               lambda(iqpt,ii,iband,1)=zero
               lambda(iqpt,ii,iband,2)=zero
             end if
!            write(std_out,'(a,i4,a,3f12.8)') 'aa=',aa,' lambda',lambda(iqpt,ii,iband,1),lambda(iqpt,ii,iband,2),&
!            lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1)+lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2)
           end do                               !iband
           if (abs(testnorm-1.0d0)>1.0d-8) write(std_out,'(a,f12.8)') 'due to zmnbld, ERROR with lambda norm',testnorm
         end do                                !free Wannier states

         ABI_DEALLOCATE(eps1)
         ABI_DEALLOCATE(ut1)
       end if                        !if there are non-frozen states

!      ***** write lambdas
!      DEBUG
!      write(std_out,*) ' after updating lambda:'
!      do ii=1,nwnn
!      write(std_out,*) 'current lambdas for the Wannier state no ',ii
!      do iband=1,qsize(iqpt,1)
!      write(std_out,'(a,3f12.8)') 'lambda',lambda(iqpt,ii,iband,1),lambda(iqpt,ii,iband,2),&
!      lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1)+lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2)
!      end do
!      end do
!      ENDDEBUG

     end do                        !q points

!    ***** update wannier information
     do iqpt=1,nqpt
       do ii=1,nwnn              ! loop over Wannier states
         do iband=1,qsize(iqpt,1)            ! loop over global bands
           wannval(iqpt,ii)=wannval(iqpt,ii)+(lambda(iqpt,ii,iband,1)*lambda(iqpt,ii,iband,1)+&
&           lambda(iqpt,ii,iband,2)*lambda(iqpt,ii,iband,2))*eigval(iqpt,g_subsp(iqpt,iband))
         end do                               !global bands
!        write(std_out,'(a,i5,f14.7)') ' wannval ',iqpt,wannval(iqpt,ii)
       end do                                !Wannier states
     end do                        !q points

!    do iqpt=1,nqpt
!    write(std_out,'(a,i10,i5,3f14.7)') ' Wanneig ',step,iqpt,wannval(iqpt,:)
!    end do

!    ***** check convergence of Omega_I

     Omega_I(step)=nwnn/6-omega
     write(std_out,*) 'at step',step-1,'Omega_I(step) is',Omega_I(step)
     write(iout,*) 'at step',step-1,'Omega_I(step) is',Omega_I(step)

     do iqpt=1,nqpt
       write(iout,'(a,i10,i5,3f14.7)') ' Wanneig ',step,iqpt,wannval(iqpt,:)
     end do


     if (contwr==20) then
       contwr=0
       do iqpt=1,nqpt
         write(std_out,'(a,i10,i5,3f10.6,3f14.7)') 'Wanneig',step,iqpt,qpoint(iqpt,:),wannval(iqpt,:)
       end do
     end if
     if (step==nstom+1) exit
     if (decflg==1 .and. step>2) then
       if (Omega_I(step)>Omega_I(step-1)) then
         write(std_out,*) 'forced out from Omega minimum'
         exit
       end if
     end if
     if (abs( (Omega_I(step)-Omega_I(step-1)) / Omega_I(step-1) )<tolomi) exit
     Zmn1=Zmn

   end do                         !main loop

   do iqpt=1,nqpt
     write(std_out,'(i5,3f10.5,a,3f14.7)') iqpt,qpoint(iqpt,:),' Wannier state eigenvalue ',wannval(iqpt,:)
     write(iout,'(i5,3f10.5,a,3f14.7)') iqpt,qpoint(iqpt,:),' Wannier state eigenvalue ',wannval(iqpt,:)
   end do

   if (prtvol==-9) then     !after minimization of Omega_I
     write(std_out,*) 'STOP requested by user in lwf.f after minimization of Omega_I'
     stop
   end if


!  ***************************************************************************
!  compute the Wannier eigenvectors
   call wanvec(atmass,eigvect,g_subsp,iwf,lambda,maxqsize,natom,nqpt,nwnn,qpoint,qsize,wannvect)


 else     ! the eigenvectors of the Wannier states are available
   call rwwan(-1,iwf,natom,nqpt,nwnn,qpoint,wannvect)
 end if


!***************************************************************************
!build the lattice Wannier functions
 call bldlwf(grdsize,iout,natom,nqpt,nwnn,qpoint,rcenter,wannvect)

!***************************************************************************
!check addition in the center

!***************************************************************************
!end of everything

 write(std_out,*) ' End of LWF'

 if (irwfl==1) then
   ABI_DEALLOCATE(atmass)
   ABI_DEALLOCATE(dprarr)
   ABI_DEALLOCATE(eigval)
   ABI_DEALLOCATE(eigvect)
   ABI_DEALLOCATE(f_subsp)
   ABI_DEALLOCATE(indeps)
   ABI_DEALLOCATE(indwnz)
   ABI_DEALLOCATE(ingss)
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(g_subsp)
   ABI_DEALLOCATE(mmnkb)
   ABI_DEALLOCATE(nwnz)
   ABI_DEALLOCATE(Omega_I)
   ABI_DEALLOCATE(qneigh)
   ABI_DEALLOCATE(qpoint)
   ABI_DEALLOCATE(qsize)
   ABI_DEALLOCATE(wannvect)
   ABI_DEALLOCATE(z_subsp)
   ABI_DEALLOCATE(Zmn)
   ABI_DEALLOCATE(Zmn1)
   ABI_DEALLOCATE(Zmn2)
 else
   ABI_DEALLOCATE(qpoint)
   ABI_DEALLOCATE(wannvect)
 end if

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program lwf
!!***
