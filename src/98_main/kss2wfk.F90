!{\src2tex{textfont=tt}}
!!****p* ABINIT/kss2wfk
!! NAME
!! kss2wfk
!!
!! FUNCTION
!!  Utility for converting a KSS file into a WFK file.
!!  KSS files do not support istwf_k>2. If explicitly asked by the user,
!!  the tool can convert the Fourier components in the Abinit WFK format
!!  taking advantage of time-reversal symmetry at high-symmetry k-points
!!  Note that the KSS file must have been generated with a value of kptopt
!!  that is compatible with time reversal. No check is done here (except for
!!  nspinor) as kptopt is not reported in the abinit header.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  Output is written on file.
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,destroy_mpi_enreg,get_kg,hdr_clean,hdr_io,hdr_io_etsf
!!      initmpi_seq,metric,mpi_comm_rank,mpi_comm_size,nullify_mpi_enreg,prompt
!!      read_kss_header,rwwf,unpack_eneocc,wffclose,wffkg,wffopen,wrtout
!!      xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program kss2wfk

 use defs_basis
 use defs_abitypes
 use m_header
 use m_xmpi
 use m_wffile
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools,   only : prompt, get_unit
 use m_gsphere,    only : get_kg
 use m_ebands,     only : unpack_eneocc
 use m_io_kss,     only : read_kss_header

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kss2wfk'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------

!Local variables-------------------------------
!scalars
 integer,parameter :: accesswff=IO_MODE_FORTRAN,prtvol=0,master=0,headform0=0
 integer,parameter :: formeig0=0,optkg1=1,icg0=0,tim_rwwf0=0,option2=2,mkmem1=1
 integer :: ik_ibz,spin,ig_k,ig_kss,istwf_k,npw_k,il,itypat,ntfound
 integer :: kss_unt,kss_nsym,kss_nbnds,kss_npw,kss_mpsang,ierr,wfk_unt
 integer :: spaceComm,my_rank,rdwr,fform,band,nkibz,nsppol,mband,bantot
 integer :: nband_k,nspinor,ispinor,iatom,natom,use_istwkf
 integer :: nband_disk,mcg,spad_kss,spad_wfn
 real(dp) :: kss_nelect,kss_ecut_eff,ucvol
 logical :: found,wfk_exists
 character(len=500) :: msg
 character(len=fnlen) :: kss_fname,wfk_fname
 type(MPI_type) :: MPI_enreg
 type(Hdr_type) :: KSS_Hdr
 type(Wffile_type) :: Wff
!arrays
 integer :: gk_search(3)
 integer,pointer :: kss_gvec(:,:),kg_k(:,:)
 integer,allocatable :: k2gamma(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),kpoint(3)
 real(dp),pointer :: occ_k(:)
 real(dp),target,allocatable :: kss_occ(:,:,:)
 real(dp),allocatable :: kss_enek(:),cg_k(:,:)
 complex(dpc),allocatable :: kss_ugd(:)

! *************************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

 call xmpi_init()

 call nullify_mpi_enreg(MPI_enreg)
 call initmpi_seq(MPI_enreg)

 spaceComm = xmpi_self
#ifdef HAVE_MPI
 MPI_enreg%world_comm=xmpi_world
 MPI_enreg%world_group=MPI_GROUP_NULL
 call MPI_COMM_RANK(xmpi_world,MPI_enreg%me,ierr)
 call MPI_COMM_SIZE(xmpi_world,MPI_enreg%nproc,ierr)
 MPI_enreg%paral_compil=1
 MPI_enreg%paral_compil_respfn=0
 MPI_enreg%paral_level=2
 spaceComm=MPI_enreg%world_comm
#endif

 if (xcomm_size(spaceComm)>1) then
   MSG_ERROR("kss2wfk not programmed for parallel execution, Run it in sequential")
 end if

 my_rank=xcomm_rank(spaceComm)

 call prompt("Enter name of the KSS file to be converted: ",kss_fname)
 call prompt("Enter name of the WFK file to be produced: ",wfk_fname)

 inquire(file=wfk_fname,exist=wfk_exists)
 if (wfk_exists) then
   write(msg,'(3a)')' File ',TRIM(wfk_fname),' already exists. Cannot overwrite!'
   MSG_ERROR(msg)
 end if

!1) Read header with basic dimensions from the KSS file.
 call read_kss_header(kss_unt,kss_fname,accesswff,prtvol,kss_nsym,kss_nbnds,kss_npw,kss_mpsang,kss_nelect,kss_gvec,KSS_Hdr)

 call metric(gmet,gprimd,-1,rmet,KSS_Hdr%rprimd,ucvol)

 nspinor      = KSS_Hdr%nspinor
 natom        = KSS_Hdr%natom
 nsppol       = KSS_Hdr%nsppol
 mband        = MAXVAL(KSS_Hdr%nband)
 bantot       = SUM(KSS_Hdr%nband)
 nkibz        = KSS_Hdr%nkpt
 kss_ecut_eff = KSS_Hdr%ecut_eff

!Time reversal storage mode.
 use_istwkf=0
 if (nspinor==1) then ! FIXME: kptopt cannot be checked.
   call prompt("Do you want to take advange of time-reversal (istwfk) (0-->no, 1-->yes): ",use_istwkf)
   if ( ALL(use_istwkf /= (/0,1/)) ) then
     MSG_ERROR("Please use 0 or 1")
   end if
 end if

 ABI_ALLOCATE(kss_occ,(mband,nkibz,nsppol))

 call unpack_eneocc(nkibz,nsppol,mband,KSS_Hdr%nband,bantot,KSS_Hdr%occ,kss_occ)

!Have to change npwarr(:) using the values of the k-centered basis set.
 do ik_ibz=1,nkibz
   kpoint  = KSS_Hdr%kptns(:,ik_ibz)

   if (use_istwkf==0) then
     istwf_k = KSS_Hdr%istwfk(ik_ibz)
   else ! Use time-reversal for special k-points.
     istwf_k = set_istwfk(kpoint)
   end if

!  K-centered basis set and mapping.
   call get_kg(kpoint,istwf_k,kss_ecut_eff,gmet,npw_k,kg_k)
   KSS_Hdr%npwarr(ik_ibz) = npw_k
   ABI_DEALLOCATE(kg_k)
 end do

 wfk_unt = get_unit()
 call wrtout(std_out," about to write "//TRIM(wfk_fname),"COLL")
!
!2) Init Wff structure for accessing WFK
!
 call WffOpen(accesswff,spaceComm,wfk_fname,ierr,Wff,master,my_rank,wfk_unt)
!
!3) Copy KSS header to the new WFK file
 rdwr=2; fform=2
 if (ANY(Wff%accesswff == (/IO_MODE_FORTRAN, IO_MODE_FORTRAN_MASTER, IO_MODE_MPI/) )) then
   call hdr_io(fform,KSS_Hdr,rdwr,Wff)
   call WffKg(Wff,1) ! Will write Kg vectors.
 else if (Wff%accesswff==IO_MODE_ETSF .and. my_rank==Wff%master) then
   call hdr_io_etsf(fform,KSS_Hdr,rdwr,Wff%unwff)
 end if
!
!4) For each spin and k-point, do:
!a) read waves from KSS
!b) convert from G- to k-centered basis set
!c) Write G vectors, energies, occ and u(G) on file.
!
 ABI_ALLOCATE(kss_enek,(kss_nbnds))
 ABI_ALLOCATE(kss_ugd,(kss_npw*nspinor))

 do spin=1,nsppol
   do ik_ibz=1,nkibz

     kpoint  = KSS_Hdr%kptns(:,ik_ibz)

     if (use_istwkf==0) then
       istwf_k = KSS_Hdr%istwfk(ik_ibz)
     else ! Use time-reversal for special k-points.
       istwf_k = set_istwfk(kpoint)
     end if

     nband_k = KSS_Hdr%nband(ik_ibz+(spin-1)*nkibz)

!    WFK stuff.
     nband_disk= nband_k
     mband     = nband_k

!    K-centered basis set and mapping.
     call get_kg(kpoint,istwf_k,kss_ecut_eff,gmet,npw_k,kg_k)

     mcg = npw_k*nspinor*mband
     ABI_ALLOCATE(cg_k,(2,mcg))

     ABI_ALLOCATE(k2gamma,(npw_k))
     k2gamma=kss_npw+1
     do ig_k=1,npw_k
       gk_search(:) = kg_k(:,ig_k)
       ig_kss=0; found=.FALSE.
       do while (ig_kss<kss_npw.and..not.found)
         ig_kss=ig_kss+1
         if (ALL(kss_gvec(:,ig_kss)==gk_search)) then
           found=.TRUE.
           k2gamma(ig_k)=ig_kss
         end if
       end do
     end do

     ntfound = (COUNT(k2gamma==kss_npw+1))
     if (ntfound/=0) then
       write(msg,"(a,i0)")" kG-vectors not found: ",ntfound
       MSG_BUG(msg)
     end if

     if (KSS_Hdr%usepaw==0) then  ! Skip Kleynmann-Bylander form factor and derivatives for this k.
       do itypat=1,KSS_Hdr%ntypat
         do il=1,kss_mpsang
           read(kss_unt) !vkbdb(:,itypat,il)
           read(kss_unt) !vkbdd(:,itypat,il)
         end do
       end do
     end if

     read(kss_unt) kss_enek(1:kss_nbnds)

     do band=1,kss_nbnds
!      
!      Read and store KSS ug_{bks}.
       read(kss_unt) kss_ugd(1:kss_npw*nspinor)

!      Convert from gamma- to k-centered G-spheres and store result in cg_k
       do ispinor=1,nspinor
         spad_kss=(ispinor-1)*kss_npw
         spad_wfn=(ispinor-1)*npw_k + (band-1)*npw_k*nspinor
         do ig_k=1,npw_k
           ig_kss=k2gamma(ig_k)
!          write(std_out,*)"ig_k",ig_k+spad_wfn,mcg
!          write(std_out,*)"ig_kss",ig_kss+spad_kss,kss_npw*nspinor
           if (ig_kss/=kss_npw+1) then
             cg_k(1,ig_k+spad_wfn)=DBLE( kss_ugd(ig_kss+spad_kss))
             cg_k(2,ig_k+spad_wfn)=AIMAG(kss_ugd(ig_kss+spad_kss))
           else
             cg_k(:,ig_k+spad_wfn)=zero
           end if
         end do
       end do
!      
!      Skip cprj if PAW.
       if (KSS_Hdr%usepaw==1) then
         do ispinor=1,nspinor
           do iatom=1,natom
             read(kss_unt)  !Cprj(iatom,ispinor)%cp
           end do
         end do
       end if
     end do ! band

     occ_k => kss_occ(1:nband_k,ik_ibz,spin)

     call rwwf(cg_k,kss_enek,formeig0,headform0,icg0,ik_ibz,spin,kg_k,mband,mcg,MPI_enreg,nband_k,&
&     nband_disk,npw_k,nspinor,occ_k,option2,optkg1,tim_rwwf0,Wff)

     ABI_DEALLOCATE(k2gamma)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(cg_k)

   end do !ik_ibz
 end do !spin
!
!* Close the wavefunction file (and do NOT delete it !)
 call WffClose(Wff,ierr)
 close(kss_unt)
!
!Free local memory (not needed but it's a good programming habit).
 ABI_DEALLOCATE(kss_gvec)
 ABI_DEALLOCATE(kss_enek)
 ABI_DEALLOCATE(kss_occ)
 ABI_DEALLOCATE(kss_ugd)

 call hdr_clean(KSS_Hdr)
 call destroy_mpi_enreg(mpi_enreg)

 call xmpi_end()

 end program kss2wfk
!!***
