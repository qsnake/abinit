!{\src2tex{textfont=tt}}
!!****p* ABINIT/cut3d
!! NAME
!! cut3d
!!
!! FUNCTION
!! Main routine for the analysis of the density and potential files,
!! as well as other files with the ABINIT header.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (GMR, RC, LSI, XG, NCJ, JFB, MCote, LPizzagalli)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  (main program)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,date_and_time,destroy_mpi_enreg,handle_ncerr,hdr_clean
!!      hdr_io,hdr_io_etsf,herald,hirsh,leave_new,lineint,localorb_s
!!      mpi_comm_rank,mpi_comm_size,nullify_mpi_enreg,planeint,pointint,rrho
!!      rtau,volumeint,wffile,wffopen,wrtout,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program cut3d

 use defs_basis
#if defined FC_NAG
 use f90_unix
#endif
#if defined HAVE_TRIO_NETCDF
 use netcdf
#endif
 use defs_datatypes
 use defs_abitypes
 use m_wffile
 use m_build_info
 use m_xmpi
#if defined HAVE_MPI2
 use mpi
#endif

 use m_header,          only : hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cut3d'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
 use interfaces_83_cut3d
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------

!Local variables-------------------------------
  ! natom = number of atoms in the unit cell
  ! nr1,nr2,nr3 = grid size (nr1 x nr2 x nr3 = filrho dimension)
  ! ntypat = number of atom types
  ! ucvol = unit cell volume (> 0)
  ! densfileformat = flag for the format of the density file:
  !         0 = ASCII
  !         1 = binary
  !         2 = binary (ETSF)
  ! denval = density value exported by interpol, to be wrote in the output file
  ! filrho = name of the density file (ASCII or binary)
  ! filtau = name of the atomic position file (Xmol format)
  ! acell = unit cell parameters
  ! rprim = orientation of the unit cell axes
  character(len=*), parameter :: INPUTfile='cut.in'
  !For NetCDF *******************************************************
#if defined HAVE_TRIO_NETCDF
  character(len=3), parameter :: monnam(12)=(/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)
#endif
 character(len=1) :: outputchar,blank=' '
!scalars
 integer :: atomicnumVarID,atomposiVarID,dataVarID,dd,densfileformat,exchn2n3d
 integer :: fform0,grid1VarID,grid2VarID,grid3VarID,gridshift1,gridshift2
 integer :: gridshift3,gridsize1DimID,gridsize2DimID,gridsize3DimID,i1,i2,i3
 integer :: iatom,ierr,ifiles,igrid,ii,ii1,ii2,ii3,index,iprompt,ir1,ir2,ir3,ispden
 integer :: itask,jfiles,latDimID,latticevecVarID,mfiles,mm,natom,nbatomDimID,ncid,nfiles,nr1,nr2
 integer :: nr3,nr1_stored,nr2_stored,nr3_stored,nrws,nspden,nspden_stored
 integer :: ntypat,originVarID,paral_kgb=0,posDimID,rdwr,ncerr
 integer :: titlechoice,unitfi,yyyy
 real(dp) :: dotdenpot,maxmz,normz,sumdenpot,ucvol
 real(dp) :: xm,xnow,xp,ym,ynow,yp,zm,znow,zp
 logical :: filexist
 character(len=10) :: strtime
 character(len=11) :: stridate
 character(len=24) :: codename
 character(len=5) :: strzone
 character(len=50) :: chr_inputfname
 character(len=500) :: filetitle
 character(len=8) :: strdat
 character(len=fnlen) :: filetsf,filnam,filrho,filrho_tmp,filtau
 type(hdr_type) :: hdr
 type(MPI_type) :: mpi_enreg
 type(wffile_type) :: wff
!arrays
 integer :: values(8)
 integer, allocatable :: isdenpot(:)
 real(dp) :: acell(3),gridwavefun1(3,2),gridwavefun2(3,2),gridwavefun3(3,2)
 real(dp) :: originatt(3,3),rprim(3,3),rprimd(3,3),shift_tau(3)!,value(2)
 real(dp) :: xcart2(3)
 real(dp),allocatable :: grid(:,:,:),grid_full(:,:,:,:),grid_full_stored(:,:,:,:,:)
 real(dp),allocatable :: gridtt(:,:,:),nuclz(:)
 real(dp),allocatable :: tau2(:,:),xcart(:,:),xred(:,:)
 real(dp),allocatable :: rhomacu(:,:)

 real(dp),allocatable :: gridmz(:,:,:),gridmy(:,:,:),gridmx(:,:,:)
 character(len=fnlen), allocatable :: filrho_stored(:)
#if defined HAVE_MPI_IO
 character(len=500) :: message
#endif

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
!inside cut3d.F90.


 codename='CUT3D '//repeat(' ',18)
 call herald(codename,abinit_version,std_out)
!YP: calling dump_config() makes tests fail => commented
!call dump_config()

!BIG LOOP on files
 mfiles=10
 ABI_ALLOCATE(isdenpot,(mfiles))
 isdenpot=0
 ABI_ALLOCATE(filrho_stored,(mfiles))

 do ifiles=1,mfiles

!  Get name of density file
   write(std_out,*)
   write(std_out,*) ' What is the name of the 3D function (density, potential or wavef) file ?'
   read(5,'(a)')filrho
   filrho_tmp=adjustl(filrho)
   do ii=1,len_trim(filrho_tmp)
     if(filrho_tmp(ii:ii)==blank)then
       filrho=trim(filrho_tmp(1:ii-1))
       exit
     end if
   end do
   write(std_out,*) ' => Your 3D function file is : ',trim(filrho)
   write(std_out,*)

!  Checking the existence of data file
   inquire (file=filrho,exist=filexist)
   if (.NOT. filexist) then
     write(std_out,*) 'Error, missing data file: ',filrho
     stop
   end if

!  Get its type
   write(std_out,*) ' Does this file contain formatted 3D ASCII data (=0)  '
   write(std_out,*) '  or unformatted binary header + 3D data        (=1) ?'
   write(std_out,*) '  or ETSF binary                                (=2) ?'
   read(5,*)densfileformat

!  Treat the different cases : formatted or unformatted
   if(densfileformat==1 .or. densfileformat == 2)then

     if (densfileformat == 1) then
       write(std_out,*) ' 1 => Your file contains unformatted binary header + 3D data'
     else
       write(std_out,*) ' 2 => Your file contains ETSF data'
     end if
     write(std_out,*) ' The information it contains should be sufficient.'

     if (densfileformat == 1) then
       open(unit=19,file=filrho,form='unformatted',status='old')
       write(std_out,'(a,a,a,i4)' )&
&       '  cut3d : read file ',trim(filrho),' from unit number 19.'
       write(std_out,*)
!      Read the header
       rdwr=1 ; unitfi=19
       call hdr_io(fform0,hdr,rdwr,unitfi)
       wff%unwff=19
       wff%accesswff=IO_MODE_FORTRAN
     else
!      We remove -etsf.nc from the file name.
       write(filetsf,"(A)") filrho(1:len(trim(filrho)) - 8)
!      Note that the MPI information are dummy, to avoid errors when consistency checks inside wffopen
       call WffOpen(3, mpi_enreg%world_comm, filetsf, ierr, wff, mpi_enreg%me, mpi_enreg%me, 19)
       write(std_out,'(a,a,a,i4)' )&
&       '  cut3d : read file ',trim(filrho),'.'
       write(std_out,*)
!      Read the header
       rdwr=1 ; unitfi=wff%unwff
       call hdr_io_etsf(fform0,hdr,rdwr,unitfi)
     end if

!    Echo part of the header
     rdwr=4 ; unitfi=6
     call hdr_io(fform0,hdr,rdwr,unitfi)

     natom=hdr%natom
     nr1=hdr%ngfft(1)
     nr2=hdr%ngfft(2)
     nr3=hdr%ngfft(3)
     nspden=hdr%nspden
     ntypat=hdr%ntypat
     rprimd(:,:)=hdr%rprimd(:,:)

!    Need to know natom in order to allocate xcart
     ABI_ALLOCATE(xcart,(3,natom))
     ABI_ALLOCATE(xred,(3,natom))
     xred(:,:)=hdr%xred(:,:)
     do iatom=1,natom
       xcart(:,iatom)=rprimd(:,1)*xred(1,iatom)+&
&       rprimd(:,2)*xred(2,iatom)+&
&       rprimd(:,3)*xred(3,iatom)
     end do

     if(fform0>50)then
       ispden=0
       if(nspden/=1)then
         write(std_out,'(a)' )' '
         write(std_out,'(a)' )' * This file contains more than one spin component,'
         write(std_out,'(a,i3,a)' )'  (indeed, nspden=',nspden,' )'
         write(std_out,'(a)' )'  Some of the tasks that you will define later will concern all spin components.'
         write(std_out,'(a)' )'  Others tasks might require you to have chosen among the following :'
       end if
       if(nspden==2)then
         write(std_out,'(a)' )'   ispden= 0 ==> Total density'
         write(std_out,'(a)' )'   ispden= 1 ==> spin-up density'
         write(std_out,'(a)' )'   ispden= 2 ==> spin-down density'
         write(std_out,'(a)' )'   ispden= 3 ==> spin-polarization (or magnetization) density'
         write(std_out,'(a)' )'                 spin up - spin down difference.'
       end if
       if(nspden==4)then
         write(std_out,'(a)' )'   ispden= 0 ==> Total density'
         write(std_out,'(a)' )'   ispden= 1 ==> magnetization in the x direction'
         write(std_out,'(a)' )'   ispden= 2 ==> magnetization in the y direction'
         write(std_out,'(a)' )'   ispden= 3 ==> magnetization in the z direction'
         write(std_out,'(a)' )'   ispden= 4 might be used to plot the magnetization (3D) in the XCrysDen format,'
       end if
       if(nspden/=1)then
         write(std_out,*)'  Please define ispden :'
         read(5,*)ispden
         write(std_out,'(a,i3)' )' You entered ispden=',ispden
       end if
     end if

   else if(densfileformat==0)then

     write(std_out,*) ' 0 => Your file contains formatted 3D ASCII data'
     write(std_out,*) ' The complementary information is taken from ',trim(INPUTfile)

!    Checking the existence of input file
     inquire (file=INPUTfile,exist=filexist)
     if (.NOT. filexist) then
       write(std_out,*) 'Error, missing input file ',INPUTfile
       stop
     end if

!    Read in the input file INPUTfile
     write(std_out,*)
     write(std_out,*) 'READING FROM FILE ',INPUTfile

     open(32,file=INPUTfile,status='old')
     read(32,'(a)') filtau

!    Checking the existence of atomic position file
     inquire (file=filtau,exist=filexist)
     if (.NOT. filexist) then
       write(std_out,*) 'Error, missing atomic position file ',trim(filtau)
       stop
     end if
     write(std_out,*) 'atomic position file (Xmol format):',trim(filtau)

!    Read cell parameters
     read(32,*) acell(1),acell(2),acell(3)
     do ii=1,3
       if (acell(ii) <= 0.0) then
         write(std_out,*) 'Error, invalid value for acell(',ii,'): ',acell(ii)
         stop
       end if
     end do
     read(32,*) rprim

     do ii=1,3
       rprimd(:,ii)=rprim(:,ii)*acell(ii)
     end do

!    FFT grid, number of atoms, number of atom types
     read(32,*) nr1,nr2,nr3
     read(32,*) natom
     read(32,*) ntypat

     close(32)

!    Need to know natom in order to allocate xcart
     ABI_ALLOCATE(xcart,(3,natom))

     call rtau(filtau,xcart,natom,ntypat)

     write(std_out,*)

!    By default there is only one spin component, and one works with a total density
     nspden=1 ; ispden=0 ; fform0=52
     wff%unwff=19

     open(unit=wff%unwff,file=filrho,form='formatted',status='old')

   else

     write(std_out,*)'Error, value for density file format is invalid: ',densfileformat
     stop

   end if

   write(std_out,*)
   write(std_out,*) '==========================================================='
   write(std_out,*)

!  Echo the value of different input parameters
   write(std_out,*)'ECHO important input variables ...'
   write(std_out,*)
   write(std_out,*) ' Dimensional primitive vectors (ABINIT equivalent : rprimd):'
   write(std_out,'(3es16.6)' ) rprimd(1:3,1)
   write(std_out,'(3es16.6)' ) rprimd(1:3,2)
   write(std_out,'(3es16.6)' ) rprimd(1:3,3)

!  Here test the non-collinearity of rprim vectors
   ucvol=rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
&   rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
&   rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))

   if (abs(ucvol)<1.0d-12) then
     write(std_out,*)' At least two rprimd(,) vectors are collinear'
     write(std_out,*)' Please check the input rprim and acell, or rprimd.'
     write(std_out,*)' The program will stop'
     stop
   end if

   write(std_out,'(a,3i5)' ) '  Grid density (ABINIT equivalent : ngfft): ',nr1,nr2,nr3
   write(std_out,*) ' Number of atoms       :',natom
   write(std_out,*) ' Number of atomic types:',ntypat

   write(std_out,*)
   write(std_out,*) '  #    Atomic positions (cartesian coordinates - Bohr)'
   do iatom=1,natom
     write(std_out,'(i4,3es16.6)' )iatom,xcart(1:3,iatom)
   end do
   write(std_out,*)

!  ------------------------------------------------------------------------
!  Branching : either WF file, or DEN/POT file.

   if((densfileformat==1 .or. densfileformat == 2) .and. fform0<50)then
     write(std_out,*)' This file is a WF file. '
     isdenpot(ifiles)=0
     iprompt = 0 ! this needs to be initialized, as it is used after the loop on files...

     exchn2n3d=0

     write(std_out,*)" If you want to analyze one wavefunction,                   type  0 "
     write(std_out,*)" If you want to construct Wannier-type Localized Orbitals,  type  2 "
     read(*,*)ii1
     write(std_out,*)" You typed ",ii1

     if(ii1==0)then
       call wffile(hdr%ecut_eff,exchn2n3d,hdr%headform,hdr%istwfk,hdr%kptns,&
&       natom,hdr%nband,hdr%nkpt,hdr%npwarr,&
&       nr1,nr2,nr3,hdr%nspinor,hdr%nsppol,ntypat,paral_kgb,rprimd,xcart,hdr%typat,wff,hdr%znucltypat)
     else if(ii1==2)then
!      Read the name of the input file name :
       read(5,'(a)')chr_inputfname
       call localorb_S(chr_inputfname,hdr%ecut_eff,exchn2n3d,hdr%headform,hdr%istwfk,hdr%kptns,&
&       natom,hdr%nband,hdr%nkpt,hdr%npwarr,&
&       nr1,nr2,nr3,hdr%nspinor,hdr%nsppol,ntypat,paral_kgb,rprimd,xcart,hdr%typat,hdr%znucltypat)
       write(std_out,*)" "
       write(std_out,*)" ###################################################### "
       write(std_out,*)" "
       write(std_out,*)" Localized orbital files fort.1**1 for spin up "
       write(std_out,*)"                     and fort.1**2 for spin dn written."
       write(std_out,*)" "
       write(std_out,*)" ###################################################### "
     else
       write(std_out,*)" Option ",ii1," is not allowed => stop "
     end if

     call hdr_clean(hdr)

!    -------------------------------------------------------------------------
   else ! This is a DEN/POT file

!    This should become a subroutine
     write(std_out,*)' This file is a Density or Potential file '
     isdenpot(ifiles)=1

!    In the wavelet case (with isolated boundary conditions), ngfft is buffered.
     if (hdr%usewvl == 1) then
       nr1 = nr1 - 31
       nr2 = nr2 - 31
       nr3 = nr3 - 31
     end if

!    Read the function on the 3D grid
     ABI_ALLOCATE(grid,(nr1,nr2,nr3))
     ABI_ALLOCATE(grid_full,(nr1,nr2,nr3,nspden))
     ABI_ALLOCATE(gridtt,(nr1,nr2,nr3))
     ABI_ALLOCATE(gridmx,(nr1,nr2,nr3))
     ABI_ALLOCATE(gridmy,(nr1,nr2,nr3))
     ABI_ALLOCATE(gridmz,(nr1,nr2,nr3))
     call rrho(densfileformat,grid_full,nr1,nr2,nr3,nspden,wff)

!    Do not forget that the first sub-array of a density file is the total density,
!    while the first sub-array of a potential file is the spin-up potential
     if(fform0==51 .or. fform0==52)then   ! Density case

!      gridtt= grid --> Total density or potential.
!      gridmx= grid --> spin-Up density, or magnetization density in X direction.
!      gridmy= grid --> spin-Down density, or magnetization density in Y direction.
!      gridmz= grid --> spin-polarization density (Magnetization),
!      or magnetization density in Z direction.
       gridtt(:,:,:)=grid_full(:,:,:,1)
       if(nspden==2)then
         gridmx(:,:,:)=grid_full(:,:,:,2)
         gridmy(:,:,:)=grid_full(:,:,:,1)-grid_full(:,:,:,2)
         gridmz(:,:,:)=-grid_full(:,:,:,1)+two*grid_full(:,:,:,2)
       else if(nspden==4)then
         gridmx(:,:,:)=grid_full(:,:,:,2)
         gridmy(:,:,:)=grid_full(:,:,:,3)
         gridmz(:,:,:)=grid_full(:,:,:,4)
       end if

       if(nspden==1)then
         grid(:,:,:)=grid_full(:,:,:,1)
       else
         if(ispden==0)then
           grid(:,:,:)=gridtt(:,:,:)
         else if(ispden==1)then
           grid(:,:,:)=gridmx(:,:,:)
         else if(ispden==2)then
           grid(:,:,:)=gridmy(:,:,:)
         else if(ispden==3)then
           grid(:,:,:)=gridmz(:,:,:)
!          if(ispden==0)then
!          grid(:,:,:)=grid_full(:,:,:,1)
!          else if(ispden==1)then
!          grid(:,:,:)=grid_full(:,:,:,2)
!          else if(ispden==2)then
!          grid(:,:,:)=grid_full(:,:,:,1)-grid_full(:,:,:,2)
!          else if(ispden==-1)then
!          grid(:,:,:)=-grid_full(:,:,:,1)+two*grid_full(:,:,:,2)
         else if(ispden==4)then
           write(std_out,*) ' '
         else
           write(std_out,*)' Error: bad ispden value'
           stop
         end if
       end if

     else if(fform0==101 .or. fform0==102)then    ! Potential case
       if(ispden==0)then
         grid(:,:,:)=grid_full(:,:,:,1)
       else if(ispden==1 .or. ispden==2)then
         grid(:,:,:)=grid_full(:,:,:,ispden)
       else
         write(std_out,*)' Error: bad ispden value'
         stop
       end if
       gridtt = grid
     end if

     write(std_out,*)
     write(std_out,*) ' 3D function was read. Ready for further treatment.'
     write(std_out,*)
     write(std_out,*) '==========================================================='
     write(std_out,*)

!    ------------------------------------------------------------------------

!    At this moment all the input is done
!    The code knows the geometry of the system,
!    and the data file (electron density, potential, etc).
!    It will further calculate the electron density by interpolation in
!    a point, along a line or in a plane.

     do
       do
         write(std_out,*) ' What is your choice ? Type:'
         write(std_out,*) '  0 => exit'
         write(std_out,*) '  1 => point  (interpolation of data for a single point)'
         write(std_out,*) '  2 => line   (interpolation of data along a line)'
         write(std_out,*) '  3 => plane  (interpolation of data in a plane)'
         write(std_out,*) '  4 => volume (interpolation of data in a volume)'
         write(std_out,*) '  5 => 3D formatted data (output the bare 3D data - one column)'
         write(std_out,*) '  6 => 3D indexed data (bare 3D data, preceeded by 3D index)'
         write(std_out,*) '  7 => 3D Molekel formatted data '
         write(std_out,*) '  8 => 3D data with coordinates (tecplot ASCII format)'
         write(std_out,*) '  9 => output .xsf file for XCrysDen'
         write(std_out,*) ' 10 => output .dx file for OpenDx'
         write(std_out,*) ' 11 => compute atomic charge using the Hirshfeld method'
         write(std_out,*) ' 12 => NetCDF file'
         write(std_out,*) ' 14 => Gaussian/cube wavefunction module'
         read(*,*) itask
         write(std_out,'(a,a,i2,a)' ) ch10,' Your choice is ',itask,ch10

         if( 5<=itask .and. itask<=10 .or. itask==12 .or. itask==14 )then
           write(std_out,*) ch10,'  Enter the name of an output file:'
           read(*,*) filnam
           write(std_out,*) '  The name of your file is : ',trim(filnam)
         end if

         select case(itask)

           case(1)            ! point calculation
             call pointint(gridtt,gridmx,gridmy,gridmz,nr1,nr2,nr3,nspden,rprimd)
             exit

           case(2)            ! line calculation
             call lineint(gridtt,gridmx,gridmy,gridmz,nr1,nr2,nr3,nspden,rprimd)
             exit

           case(3)            ! plane calculation
             call planeint(gridtt,gridmx,gridmy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,xcart)
             exit

           case(4)            ! volume calculation
             write(std_out,*) ' Enter volume calculation'
             call volumeint(gridtt,gridmx,gridmy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,xcart)
             exit

           case(5)
!            Rewrite the data on a formatted file, just in one (or four) column(s)
             open(unit=31,file=trim(filnam),status='unknown')
             if(nspden==1)then
               do i3=1,nr3
                 do i2=1,nr2
                   do i1=1,nr1
                     write(31,'(4(es22.12))') grid(i1,i2,i3)
                   end do
                 end do
               end do
             else
               do i3=1,nr3
                 do i2=1,nr2
                   do i1=1,nr1
                     write(31,'(4(es22.12))') gridtt(i1,i2,i3), gridmx(i1,i2,i3), gridmy(i1,i2,i3), gridmz(i1,i2,i3)
                   end do
                 end do
               end do
             end if
             close(31)
             exit

           case(6)
!            Rewrite the data on a formatted file, 3D index + density
             open(unit=31,file=trim(filnam),status='unknown')
             if(nspden==1)then
               write(31,*)'   i1    i2    i3      data '
               do i3=1,nr3
                 do i2=1,nr2
                   do i1=1,nr1
                     write(31,'(3i6,4(es24.14))') i1,i2,i3,grid(i1,i2,i3)
                   end do
                 end do
               end do
             else
               if(nspden==2)then
                 write(31,*)'   i1    i2    i3     non-spin-polarized spin up  spin down  difference  '
               else if(nspden==4)then
                 write(31,*)'   i1    i2    i3     non-spin-polarized   x       y      z   '
               end if
               do i3=1,nr3
                 do i2=1,nr2
                   do i1=1,nr1
                     write(31,'(3i6,4(es24.14))') i1,i2,i3,gridtt(i1,i2,i3),gridmx(i1,i2,i3),gridmy(i1,i2,i3),gridmz(i1,i2,i3)
                   end do
                 end do
               end do
             end if ! nspden
             close(31)
             exit

           case(7)
             open(unit=31,file=trim(filnam),form='unformatted')
             xm=0 ; xp=rprimd(1,1)*Bohr_Ang
             ym=0 ; yp=rprimd(2,2)*Bohr_Ang
             zm=0 ; zp=rprimd(3,3)*Bohr_Ang
             write(std_out,'(/,a,/)' )&
&             ' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
             write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
             write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1,nr2,nr3
             write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', nr1*nr2*nr3
             write(31) xm,xp,ym,yp,zm,zp,nr1,nr2,nr3
             ABI_ALLOCATE(rhomacu,(nr1,nr2))
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   rhomacu(i1,i2)=grid(i1,i2,i3)
                 end do
               end do
               write(31) rhomacu(:,:)
             end do
             close(31)
             exit

           case (8)
             open(unit=31,file=trim(filnam),form='formatted')
             write(std_out,'(/,a,/)' )&
&             ' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
             write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
             write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1,nr2,nr3
             write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', nr1*nr2*nr3
             write(31,'(a)') 'TITLE = "  " '
             write(31,'(a)') 'VARIABLES = "X"  "Y"  "Z" (all three in Angstrom)  "DENSITY or POTENTIAL" (atomic units) '
             write(31,'(3(a,i6),a)') 'ZONE I=',nr1, ' J=', nr2, ' K=', nr3, ' F=POINT'
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   xnow = rprimd(1,1)*(i1-1)/nr1 + rprimd(1,2)*(i2-1)/nr2 + rprimd(1,3)*(i3-1)/nr3
                   ynow = rprimd(2,1)*(i1-1)/nr1 + rprimd(2,2)*(i2-1)/nr2 + rprimd(2,3)*(i3-1)/nr3
                   znow = rprimd(3,1)*(i1-1)/nr1 + rprimd(3,2)*(i2-1)/nr2 + rprimd(3,3)*(i3-1)/nr3
                   write(31,'(4es22.15)') Bohr_Ang*xnow, Bohr_Ang*ynow, Bohr_Ang*znow, grid (i1,i2,i3)
                 end do
               end do
             end do
             close(31)
             exit

           case (9)
             open(unit=31,file=trim(filnam),form='formatted')
             xm=0 ; xp=rprimd(1,1)*Bohr_Ang
             ym=0 ; yp=rprimd(2,2)*Bohr_Ang
             zm=0 ; zp=rprimd(3,3)*Bohr_Ang
             write(std_out,'(/,a,/)' )&
&             ' Extremas (x,y,z) of the cube in which the molecule is placed, in Angstroms'
             write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
             write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1+1,nr2+1,nr3+1
             write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', (nr1+1)*(nr2+1)*(nr3+1)
             write(std_out,*) '  znucl = ', hdr%znucltypat, ' type = ', hdr%typat, ' ntypat = ', ntypat

             gridshift1 = 0
             gridshift2 = 0
             gridshift3 = 0
             write(std_out,*) 'Do you want to shift the grid along the x,y or z axis (y/n)?'
             write(std_out,*)
             shift_tau(:) = zero
             read (*,*) outputchar
             if (outputchar == 'y' .or. outputchar == 'Y') then
               write(std_out,*) 'Give the three shifts (x,y,z < ',nr1,nr2,nr3,') :'
               write(std_out,*)
               read (*,*) gridshift1, gridshift2, gridshift3
               shift_tau(:) = gridshift1*rprimd(:,1)/(nr1+1) + gridshift2*rprimd(:,2)/(nr2+1) + gridshift3*rprimd(:,3)/(nr3+1)
             end if
!            
!            Generate translated coordinates to match density shift
!            
             ABI_ALLOCATE(tau2,(3,natom))
             do iatom = 1,natom
               tau2(:,iatom) = xcart(:,iatom) - shift_tau(:)
             end do
!            ################################################################### (LD)
!            Option only available for "xcrysden" format as documented at the beginning
             if (ispden==4) then
!              It is necessary to know previously how many atoms will be used.
!              in order to plot the necessary magnetization arrows only.
               write(std_out,*)'Is it possible to decrease the number of arrows in order to improve the'
               write(std_out,*)'visualization in the screen, and decrease the size of the xcrysden output file.'
               write(std_out,*)'How many arrows would you like to skip? 0 = take all. 1 = skip every other point...'
               read (*,*) nrws
               nrws=nrws+1
               index=natom
               maxmz=0.0
               do i1=1,nr1,nrws
                 do i2=1,nr2,nrws
                   do i3=1,nr3,nrws
                     normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                     if(normz > maxmz) maxmz=normz
                   end do
                 end do
               end do
               if(abs(maxmz)<tol10)then
                 write(std_out,*)' At least, one of the components must differ from zero.',&
&                 ' This is not the case. Stopping now.'
                 stop
               end if
               do i1=1,nr1,nrws
                 do i2=1,nr2,nrws
                   do i3=1,nr3,nrws
                     normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                     if(0.1*maxmz <= normz) index=index+1
                   end do
                 end do
               end do

               write(31,'(1X,A)') 'CRYSTAL'
               write(31,'(1X,A)') 'PRIMVEC'
               do i1 = 1,3
                 write(31,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
               end do
               write(31,'(1X,A)') 'PRIMCOORD'
               write(31,*) index, '1'
!              
!              write out atom types and positions
!              
               do iatom = 1,natom
                 write(31,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),&
&                 Bohr_Ang*tau2(1:3,iatom)
               end do
!              
!              write out magnetization vectors.
!              xcrysden consider these as X (dummy) atoms.
!              
               do i1=1,nr1,nrws
                 do i2=1,nr2,nrws
                   do i3=1,nr3,nrws
                     normz=gridmx(i1,i2,i3)**2+gridmy(i1,i2,i3)**2+gridmz(i1,i2,i3)**2
                     if(0.1*maxmz <= normz) then
                       xcart2 = matmul (rprimd, (/(i1-one)/nr1, (i2-one)/nr2, (i3-one)/nr3/))
                       write(31,'(A,1X,6(ES17.10,2X))')'X',&
                       Bohr_Ang*(xcart2(1)-shift_tau(1)),&
                       Bohr_Ang*(xcart2(2)-shift_tau(2)),&
                       Bohr_Ang*(xcart2(3)-shift_tau(3)),&
                       gridmx(i1,i2,i3),&
                       gridmy(i1,i2,i3),&
                       gridmz(i1,i2,i3)
                     end if
                   end do
                 end do
               end do
             else
!              ################################################################### (LD)
!              
!              normal case: output density or potential (scalar field)
!              
               write(31,'(1X,A)')  'DIM-GROUP'
               write(31,*) '3  1'
               write(31,'(1X,A)') 'PRIMVEC'
               do i1 = 1,3
                 write(31,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
               end do
               write(31,'(1X,A)') 'PRIMCOORD'
               write(31,*) natom, ' 1'
               do iatom = 1,natom
                 write(31,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),&
&                 Bohr_Ang*tau2(1:3,iatom)
               end do
               write(31,'(1X,A)') 'ATOMS'
               do iatom = 1,natom
                 write(31,'(i9,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),&
&                 Bohr_Ang*tau2(1:3,iatom)
               end do
!              write(31,'(1X,A)') 'FRAMES'
               write(31,'(1X,A)') 'BEGIN_BLOCK_DATAGRID3D'
               write(31,*) 'datagrids'
               write(31,'(1X,A)') 'DATAGRID_3D_DENSITY'
               write(31,*) nr1+1,nr2+1,nr3+1
               write(31,*) '0.0 0.0 0.0 '
               do i1 = 1,3
                 write(31,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(i2,i1), i2=1,3)
               end do

               index = 0
               do ir3=gridshift3+1,nr3+1
                 ii3=mod(ir3-1,nr3) + 1
                 do ir2=gridshift2+1,nr2+1
                   ii2=mod(ir2-1,nr2) + 1
                   do ir1=gridshift1+1,nr1+1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                   do ir1=1,gridshift1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                 end do
                 do ir2=1,gridshift2
                   ii2=mod(ir2-1,nr2) + 1
                   do ir1=gridshift1+1,nr1+1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                   do ir1=1,gridshift1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                 end do
               end do
               do ir3=1,gridshift3
                 ii3=mod(ir3-1,nr3) + 1
                 do ir2=gridshift2+1,nr2+1
                   ii2=mod(ir2-1,nr2) + 1
                   do ir1=gridshift1+1,nr1+1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                   do ir1=1,gridshift1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                 end do
                 do ir2=1,gridshift2
                   ii2=mod(ir2-1,nr2) + 1
                   do ir1=gridshift1+1,nr1+1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                   do ir1=1,gridshift1
                     ii1=mod(ir1-1,nr1) + 1
                     write(31,'(e20.5,2x)',ADVANCE='NO') grid(ii1,ii2,ii3)
                     index = index+1
                     if (mod (index,6) == 0) write (31,*)
                   end do
                 end do
               end do
               write (31,*)
               write(31,'(1X,A)') 'END_DATAGRID_3D'
               write(31,'(1X,A)') 'END_BLOCK_DATAGRID3D'

             end if

             close(31)
             exit

           case (10)      ! formatted for OpenDX
             open(unit=31,file=trim(filnam)//'.dx',form='formatted')
             open(unit=32,file=trim(filnam)//'.xyz',form='formatted')
             write(31, '(a,2x,3i5)' )'object 1 class gridpositions counts',nr3,nr2,nr1
             write(31, '(a)' )'origin 0 0 0'
             write(31,'(a,3(es17.10,2x))')'delta ',(Bohr_Ang*rprimd(i1,3)/nr3, i1=1,3)
             write(31,'(a,3(es17.10,2x))')'delta ',(Bohr_Ang*rprimd(i1,2)/nr2, i1=1,3)
             write(31,'(a,3(es17.10,2x))')'delta ',(Bohr_Ang*rprimd(i1,1)/nr1, i1=1,3)
             write(31, '(a,2x,3i5)' )'object 2 class gridconnections counts',nr3,nr2,nr1
             write(31, '(a)' )'attribute "element type" string "cubes"'
             write(31, '(a)' )'attribute "ref" string "positions"'
             write(31, '(a,1x,i10,1x,a)' )'object 3 class array type float rank 0 items',nr1*nr2*nr3,' data follows'
             do i3=1,nr3
               do i2=1,nr2
                 do i1=1,nr1
                   write(31,'(es24.14)') grid(i1,i2,i3)
                 end do
               end do
             end do
             write(31, '(a)' )'attribute "dep" string "positions"'
             write(32, '(i6,/)' ) natom
             do iatom=1,natom
               do ii=1,3
                 xcart2(ii)=xcart(ii,iatom)
                 if (xred(ii,iatom)<-1e-4) then
                   xcart2(ii)=xcart2(ii)+rprimd(ii,ii)
                 else if (xred(ii,iatom)>rprimd(ii,ii)+1e-4) then
                   xcart2(ii)=xcart2(ii)-rprimd(ii,ii)
                 end if
               end do
               write(32,'(i8,3(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),Bohr_Ang*xcart2(1:3)
             end do
             write(31, '(a)' )'object "density" class field'
             write(31, '(a)' )'component "positions" value 1'
             write(31, '(a)' )'component "connections" value 2'
             write(31, '(a)' )'component "data" value 3'
             close(31)
             close(32)
             exit

           case(11)

             call hirsh(grid,natom,nr1,nr2,nr3,ntypat,rprimd,xcart,hdr%typat,hdr%zionpsp,hdr%znucltypat)
             exit

           case(12)
#if defined HAVE_TRIO_NETCDF
             ABI_ALLOCATE(nuclz,(natom))

             filnam = trim(filnam)//'.nc'

!            Creating the NetCDF file

             ncerr = nf90_create(filnam, nf90_clobber, ncid)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Creating file")

!            Ask for a title

             write(std_out,*) 'Do you want a title in your NetCDF file? (0 = NO, 1 = YES)'

             read(*,*) titlechoice
             do
               if (titlechoice ==0 .or. titlechoice ==1) exit
               write(std_out,*) 'The answer is not correct, you must enter an integer between 0 and 1 (0 = NO, 1 = YES)'
               read(*,*) titlechoice
             end do
             if (titlechoice ==1) then
               write(std_out,*) 'Enter your file''s title'
               read(*,'(A)') filetitle
             else
               write(std_out,*) 'No title will be added in your NetCDF file'
             end if

             originatt(1:3,1:3)=0

             do igrid = 1,3
               gridwavefun1(igrid,1)=0
               gridwavefun2(igrid,1)=0
               gridwavefun3(igrid,1)=0
               gridwavefun1(igrid,2)=Bohr_Ang*rprimd(igrid,3)/nr3
               gridwavefun2(igrid,2)=Bohr_Ang*rprimd(igrid,2)/nr2
               gridwavefun3(igrid,2)=Bohr_Ang*rprimd(igrid,1)/nr1
             end do

             do igrid =1,natom
               nuclz(igrid)=hdr%znucltypat(hdr%typat(igrid))
             end do

!            Defining dimensions

             ncerr = nf90_def_dim(ncid,"gridsize1",nr1, gridsize1DimID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

             ncerr = nf90_def_dim(ncid,"gridsize2",nr2, gridsize2DimID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

             ncerr = nf90_def_dim(ncid,"gridsize3",nr3, gridsize3DimID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

             ncerr = nf90_def_dim(ncid, "lat",3, latDimID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

             ncerr = nf90_def_dim(ncid, "pos",2, posDimID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

             ncerr = nf90_def_dim(ncid, "nbatom",natom, nbatomDimID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining dimensions")

!            Defining variables

             ncerr = nf90_def_var(ncid, "latticevec",nf90_float, (/ latDimID, latDimID /), latticevecVarID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             ncerr = nf90_def_var(ncid, "origin",nf90_float, (/ latDimID, latDimID /), originVarID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             ncerr = nf90_def_var(ncid, "atomposi",nf90_float, (/ latDimID, nbatomDimID /), atomposiVarID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             ncerr = nf90_def_var(ncid, "atomicnum",nf90_float, nbatomDimID, atomicnumVarID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             ncerr = nf90_def_var(ncid, "grid1",nf90_float, (/latDimID,posDimID/), grid1VarID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             ncerr = nf90_def_var(ncid, "grid2",nf90_float, (/latDimID,posDimID/), grid2VarID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             ncerr = nf90_def_var(ncid, "grid3",nf90_float, (/latDimID,posDimID/), grid3VarID)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             if(fform0==51 .or. fform0==52)then   ! Density case

               ncerr = nf90_def_var(ncid, "density",nf90_float, (/ gridsize1DimID, gridsize2DimID, gridsize3DimID /), dataVarID)
               if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             else if(fform0==101 .or. fform0==102)then    ! Potential case

               ncerr = nf90_def_var(ncid, "potential",nf90_float, (/ gridsize1DimID, gridsize2DimID, gridsize3DimID /), dataVarID)
               if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining variables")

             end if

!            Defining attibutes

             ncerr = nf90_put_att(ncid,latticevecVarID , "field", "latticevec,vector")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,originVarID , "field", "origin,vector")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,atomposiVarID , "field","atomposi,vector")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,latticevecVarID , "positions","origin" )
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,atomicnumVarID , "positions", "atomposi")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,dataVarID , "positions", "grid1,product,compact;grid2,product,compact;grid3,product,compact")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,dataVarID , "connections", (/nr3,nr2,nr1/))
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,latticevecVarID , "long_name", "Lattice Vectors")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,atomposiVarID , "long_name", "Atomic Positions")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid, atomicnumVarID, "long_name", "Atomic Numbers")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,latticevecVarID, "units", "Angstroms")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             ncerr = nf90_put_att(ncid,atomposiVarID, "units", "Angstroms")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             if(fform0==51 .or. fform0==52)then   ! Density case

               ncerr = nf90_put_att(ncid, dataVarID, "long_name", "Density")
               if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             else if(fform0==101 .or. fform0==102)then    ! Potential case

               ncerr = nf90_put_att(ncid, dataVarID, "long_name", "Potential")
               if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             end if

             ncerr = nf90_put_att(ncid, originVarID, "long_name", "Origin of the Lattice")
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

             if (titlechoice ==1) then
               ncerr = nf90_put_att(ncid, nf90_global, "title", filetitle)
               if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")
             end if

!            Add the creation date
             call date_and_time(strdat,strtime,strzone,values)
             yyyy=values(1)
             mm=values(2)
             dd=values(3)
             write(stridate(1:2),'(I2)') dd
             stridate(3:3)=" "
             stridate(4:7)=monnam(mm)
             write(stridate(8:11),'(I4)') yyyy

             ncerr = nf90_put_att(ncid, nf90_global,"date", stridate)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Defining attributes")

!            Ending the define mode and entering data mode

             ncerr = nf90_enddef(ncid)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Entering data mode")

!            Putting the data

             ncerr = nf90_put_var(ncid,latticevecVarID,Bohr_Ang*rprimd)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

             ncerr = nf90_put_var(ncid,originVarID,originatt)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

             ncerr = nf90_put_var(ncid,atomposiVarID, Bohr_Ang*xcart)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

             ncerr = nf90_put_var(ncid,atomicnumVarID,nuclz)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

             ncerr = nf90_put_var(ncid,grid1VarID,gridwavefun1)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

             ncerr = nf90_put_var(ncid,grid2VarID,gridwavefun2)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

             ncerr = nf90_put_var(ncid,grid3VarID,gridwavefun3)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

             ncerr = nf90_put_var(ncid,dataVarID,grid)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Putting data")

!            Closing the file

             ncerr = nf90_close(ncid)
             if (ncerr /= nf90_noerr) call handle_ncerr(ncerr, "Closing file")

             write(std_out,*) 'The NetCDF file is done'
             ABI_DEALLOCATE(nuclz)
             write(std_out,*)
             exit
#else
             write(std_out,*) 'NetCDF is not defined. You must choose another option'
             exit
#endif
           case(14)            ! CUBE file format from GAUSSIAN

             write(std_out,*)
             write(std_out,*) 'Output a cube file of 3D volumetric data'
             write(std_out,*)

!            EXAMPLE FROM THE WEB
!            CPMD CUBE FILE.
!            OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
!            3    0.000000    0.000000    0.000000
!            40    0.283459    0.000000    0.000000
!            40    0.000000    0.283459    0.000000
!            40    0.000000    0.000000    0.283459
!            8    0.000000    5.570575    5.669178    5.593517
!            1    0.000000    5.562867    5.669178    7.428055
!            1    0.000000    7.340606    5.669178    5.111259
!            -0.25568E-04  0.59213E-05  0.81068E-05  0.10868E-04  0.11313E-04  0.35999E-05

             open(unit=31,file=trim(filnam),status='unknown',form='formatted')
             write(31,'(a)') 'ABINIT generated cube file'
             write(31,'(a)') 'from cut3d tool'

             write(31,'(i9,3(1x,f12.6))') natom,0.,0.,0.
             write(31,'(i9,3(1x,f12.6))') nr1,(rprimd(ir2,1)/nr1, ir2=1,3)
             write(31,'(i9,3(1x,f12.6))') nr2,(rprimd(ir2,2)/nr2, ir2=1,3)
             write(31,'(i9,3(1x,f12.6))') nr3,(rprimd(ir2,3)/nr3, ir2=1,3)

             do iatom = 1,natom
               write(31,'(i9,4(3X,ES17.10))') nint(hdr%znucltypat(hdr%typat(iatom))),0.d0, &
&               xcart(1,iatom), &
&               xcart(2,iatom), &
&               xcart(3,iatom)
             end do

!            C ordering of the indexes
             do i1=1,nr1
               do i2=1,nr2
                 do i3=1,nr3
                   write(31,'(6(f12.6,2x))') grid(i1,i2,i3)
!                  write(31,'(4(es22.12))') gridtt(i1,i2,i3), gridmx(i1,i2,i3), gridmy(i1,i2,i3), gridmz(i1,i2,i3)
                 end do
               end do
             end do

             close(31)

             exit

           case(0)
             write(std_out,*)' Exit requested by user'
             exit

             case default
             cycle

         end select

       end do

       write(std_out,*) ' Task ',itask,' has been done !'
       write(std_out,*)
       write(std_out,'(a)') ' More analysis of the 3D file ? ( 0=no ; 1=default=yes ; 2= treat another file - restricted usage)'
       read(*,*) iprompt
       if(iprompt/=1) then
         if(densfileformat==1 .or. densfileformat == 2)then
           call hdr_clean(hdr)
         end if
         exit
       else
         cycle
       end if
     end do

   end if ! WF file or DEN/POT file

!  A maximum number of files had been previously specified, but set the actual number of files to 1 if one does not read at least one other.
   if(ifiles==1)then
     nfiles=1
     if(iprompt==2)nfiles=mfiles

!    A data structure for storing the important information should be created ...
!    Here, one supposes that the files are compatible ...
     if(isdenpot(ifiles)==1)then
       ABI_ALLOCATE(grid_full_stored,(nr1,nr2,nr3,nspden,nfiles))
       nr1_stored=nr1
       nr2_stored=nr2
       nr3_stored=nr3
       nspden_stored=nspden
     else if(isdenpot(ifiles)/=1 .and. iprompt==2)then
       write(std_out,*)' Error : in case of storage mode, the first file must be a density/potential file.'
       stop
     end if
   end if

   if(isdenpot(ifiles)==1) grid_full_stored(:,:,:,:,ifiles)=grid_full(:,:,:,:)
   if(isdenpot(ifiles)==1) filrho_stored(ifiles)=filrho

   if(allocated(xcart))ABI_DEALLOCATE(xcart)
   if(allocated(xred))ABI_DEALLOCATE(xred)
   if(allocated(grid))ABI_DEALLOCATE(grid)
   if(allocated(grid_full))ABI_DEALLOCATE(grid_full)
   if(allocated(gridtt))ABI_DEALLOCATE(gridtt)
   if(allocated(gridmx))ABI_DEALLOCATE(gridmx)
   if(allocated(gridmy))ABI_DEALLOCATE(gridmy)
   if(allocated(gridmz))ABI_DEALLOCATE(gridmz)
   if(allocated(rhomacu))ABI_DEALLOCATE(rhomacu)
   if(allocated(tau2))ABI_DEALLOCATE(tau2)

   if(iprompt/=2) then
     exit
   end if

 end do ! End big loop on files

!Will provide different information on the density and potential files
 do ifiles=1,nfiles
   if(isdenpot(ifiles)==1)then
     write(std_out,*)
     write(std_out,*) ' Provide some global information about the density and/or potential file(s)'
     exit
   end if
 end do
 do ifiles=1,nfiles
   if(isdenpot(ifiles)==1)then
     write(std_out,*)
     write(std_out, '(a,i5,3a)' ) '  File number ',ifiles,', with name "',trim(filrho_stored(ifiles)),'"'
     write(std_out, '(a,i12,a,es14.6)' ) '  Number of grid points =',nr1*nr2*nr3,' ; Volume of real space cell (Bohr^3)=',ucvol
     do ispden=1,nspden
       sumdenpot=sum(grid_full_stored(:,:,:,ispden,ifiles))
       write(std_out, '(a,i5,3a)' ) '   Spin-component number ',ispden
       write(std_out, '(a,3es16.6)' ) '      Sum of values, mean, mean times cell volume=',&
&       sumdenpot,sumdenpot/real(nr1*nr2*nr3),sumdenpot*ucvol/real(nr1*nr2*nr3)
     end do
   end if
 end do

 if(nspden==1)then
!  At present, only nspden=1 is correctly implemented, due to specificities of the treatment of the spin-density
   do ifiles=1,nfiles
     if(isdenpot(ifiles)==1)then
       write(std_out,*)
       write(std_out,'(a)') ' Provide some global joint information about the stored density and potential file(s)'
       exit
     end if
   end do
   do ifiles=1,nfiles
     if(isdenpot(ifiles)==1)then
       do jfiles=ifiles,nfiles
         if(isdenpot(jfiles)==1)then
           write(std_out,*)
           write(std_out, '(a,2i5)' )'  File numbers :',ifiles,jfiles
           do ispden=1,nspden
             dotdenpot=zero
             do ir1=1,nr1
               do ir2=1,nr2
                 do ir3=1,nr3
                   dotdenpot=dotdenpot+grid_full_stored(ir1,ir2,ir3,ispden,ifiles)*grid_full_stored(ir1,ir2,ir3,ispden,jfiles)
                 end do
               end do
             end do
             write(std_out, '(a,i5,3a)' ) '   Spin-component number ',ispden
             write(std_out, '(a,3es16.6)' ) '      Dot product of values, mean, mean times cell volume=',&
!            write(std_out, '(a,3es20.10)' ) '      Dot product of values, mean, mean times cell volume=',&
&             dotdenpot,dotdenpot/real(nr1*nr2*nr3),dotdenpot*ucvol/real(nr1*nr2*nr3)
           end do
         end if
       end do
     end if
   end do
 end if

 if(allocated(grid_full_stored))ABI_DEALLOCATE(grid_full_stored)
 if(allocated(isdenpot))ABI_DEALLOCATE(isdenpot)

 write(std_out,*)
 write(std_out,*) ' Thank you for using me'
 write(std_out,*)

!DEBUG
!write(std_out,*)' cut3d : before hdr_clean, densfileformat=',densfileformat
!write(std_out,*)' allocated(xcart)=',allocated(xcart)
!ENDDEBUG

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program cut3d
!!***
