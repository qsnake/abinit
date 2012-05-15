!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_comm
!! NAME
!! hdr_comm
!!
!! FUNCTION
!! This subroutine transmit the header structured datatype
!! initialized on one processor (or a group of processor),
!! to the other processors. It also allocate the needed
!! part of the header.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (XG, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  master = id of the master process
!!  me = id of the current process
!!  spaceComm = id of the space communicator handler
!!
!! OUTPUT
!!  (no output)
!!
!! SIDE EFFECTS
!!  hdr <type(hdr_type)>=the header. For the master, it is already
!!   initialized entirely, while for the other procs, everything has
!!   to be transmitted.
!!
!! NOTES
!! This routine is called only in the case of MPI version of the code.
!!
!! PARENTS
!!      hdr_io,hdr_io_netcdf,m_io_kss,m_io_screening
!!
!! CHILDREN
!!      mpi_barrier,mpi_bcast,rhoij_alloc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine hdr_comm(hdr,master,me,spaceComm)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_xmpi
#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_comm'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in) :: master,me,spaceComm
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: bantot,cplex,iatom,ierr,index,index2,ipsp,ispden,list_size,list_size2,natom,nkpt
 integer :: npsp,nsel,nspden,nsppol,nsym,nrhoij,ntypat
 integer,allocatable :: list_int(:)
 real(dp),allocatable :: list_dpr(:)
 character(len=fnlen) :: list_tmp
 character(len=fnlen),allocatable :: list_char(:)
#endif

! *************************************************************************

!DEBUG
!write(std_out,*)' hdr_comm : enter '
!ENDDEBUG

#if defined HAVE_MPI
!Transmit the integer scalars
 list_size=20
 ABI_ALLOCATE(list_int,(list_size))
 if (master==me)then
   list_int(1)=hdr%bantot
   list_int(2)=hdr%date
   list_int(3)=hdr%headform
   list_int(4)=hdr%intxc
   list_int(5)=hdr%ixc
   list_int(6)=hdr%natom
   list_int(7)=hdr%nkpt
   list_int(8)=hdr%npsp
   list_int(9)=hdr%nspden
   list_int(10)=hdr%nspinor
   list_int(11)=hdr%nsppol
   list_int(12)=hdr%nsym
   list_int(13)=hdr%ntypat
   list_int(14)=hdr%occopt
   list_int(15)=hdr%pertcase
   list_int(16)=hdr%usepaw
   list_int(17:19)=hdr%ngfft(1:3)
   list_int(20)=hdr%usewvl
 end if
 call MPI_BARRIER(spaceComm,ierr)
 call MPI_BCAST(list_int,list_size,MPI_INTEGER,0,spaceComm,ierr)
 if(master/=me)then
   hdr%bantot  =list_int(1)
   hdr%date    =list_int(2)
   hdr%headform=list_int(3)
   hdr%intxc   =list_int(4)
   hdr%ixc     =list_int(5)
   hdr%natom   =list_int(6)
   hdr%nkpt    =list_int(7)
   hdr%npsp    =list_int(8)
   hdr%nspden  =list_int(9)
   hdr%nspinor =list_int(10)
   hdr%nsppol  =list_int(11)
   hdr%nsym    =list_int(12)
   hdr%ntypat  =list_int(13)
   hdr%occopt  =list_int(14)
   hdr%pertcase=list_int(15)
   hdr%usepaw  =list_int(16)
   hdr%ngfft(1:3)=list_int(17:19)
   hdr%usewvl  =list_int(20)
 end if
 ABI_DEALLOCATE(list_int)

 bantot=hdr%bantot
 natom =hdr%natom
 nkpt  =hdr%nkpt
 npsp  =hdr%npsp
 nspden=hdr%nspden
 nsppol=hdr%nsppol
 nsym  =hdr%nsym
 ntypat=hdr%ntypat

 if(master/=me)then
!  Allocate all components of hdr
   ABI_ALLOCATE(hdr%istwfk,(nkpt))
   ABI_ALLOCATE(hdr%nband,(nkpt*nsppol))
   ABI_ALLOCATE(hdr%npwarr,(nkpt))
   ABI_ALLOCATE(hdr%pspcod,(npsp))
   ABI_ALLOCATE(hdr%pspdat,(npsp))
   ABI_ALLOCATE(hdr%pspso,(npsp))
   ABI_ALLOCATE(hdr%pspxc,(npsp))
   ABI_ALLOCATE(hdr%lmn_size,(npsp))
   ABI_ALLOCATE(hdr%so_psp,(npsp))
   ABI_ALLOCATE(hdr%symafm,(nsym))
   ABI_ALLOCATE(hdr%symrel,(3,3,nsym))
   ABI_ALLOCATE(hdr%typat,(natom))
   ABI_ALLOCATE(hdr%kptns,(3,nkpt))
   ABI_ALLOCATE(hdr%occ,(bantot))
   ABI_ALLOCATE(hdr%tnons,(3,nsym))
   ABI_ALLOCATE(hdr%wtk,(nkpt))
   ABI_ALLOCATE(hdr%xred,(3,natom))
   ABI_ALLOCATE(hdr%zionpsp,(npsp))
   ABI_ALLOCATE(hdr%znuclpsp,(npsp))
   ABI_ALLOCATE(hdr%znucltypat,(ntypat))
   ABI_ALLOCATE(hdr%title,(npsp))
 end if

!Transmit the integer arrays
 list_size=nkpt*(2+nsppol)+6*npsp+10*nsym+natom
 ABI_ALLOCATE(list_int,(list_size))
 if (master==me)then
   list_int(1      :nkpt             )=hdr%istwfk ; index=nkpt
   list_int(1+index:nkpt*nsppol+index)=hdr%nband  ; index=index+nkpt*nsppol
   list_int(1+index:nkpt       +index)=hdr%npwarr ; index=index+nkpt
   list_int(1+index:npsp       +index)=hdr%pspcod ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspdat ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspso  ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspxc  ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%lmn_size ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%so_psp ; index=index+npsp
   list_int(1+index:nsym       +index)=hdr%symafm ; index=index+nsym
   list_int(1+index:nsym*3*3   +index)=reshape(hdr%symrel,(/3*3*nsym/)) &
&   ; index=index+nsym*3*3
   list_int(1+index:natom      +index)=hdr%typat   ; index=index+natom
 end if
 call MPI_BARRIER(spaceComm,ierr)
 call MPI_BCAST(list_int,list_size,MPI_INTEGER,0,spaceComm,ierr)
 if(master/=me)then
   hdr%istwfk=list_int(1      :nkpt             ) ; index=nkpt
   hdr%nband =list_int(1+index:nkpt*nsppol+index) ; index=index+nkpt*nsppol
   hdr%npwarr=list_int(1+index:nkpt       +index) ; index=index+nkpt
   hdr%pspcod=list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspdat=list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspso =list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspxc =list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%lmn_size=list_int(1+index:npsp     +index) ; index=index+npsp
   hdr%so_psp =list_int(1+index:npsp   +index) ; index=index+npsp
   hdr%symafm=list_int(1+index:nsym       +index) ; index=index+nsym
   hdr%symrel=reshape(list_int(1+index:nsym*3*3   +index),(/3,3,nsym/))
   index=index+nsym*3*3
   hdr%typat  =list_int(1+index:natom      +index) ; index=index+natom
 end if
 ABI_DEALLOCATE(list_int)

!Transmit the double precision scalars and arrays
 list_size=21+3*nkpt+nkpt+bantot+3*nsym+3*natom+2*npsp+ntypat
 ABI_ALLOCATE(list_dpr,(list_size))
 if (master==me)then
   list_dpr(1)=hdr%ecut_eff
   list_dpr(2)=hdr%etot
   list_dpr(3)=hdr%fermie
   list_dpr(4)=hdr%residm
   list_dpr(5:13)=reshape(hdr%rprimd(1:3,1:3),(/9/))
   list_dpr(14)=hdr%ecut
   list_dpr(15)=hdr%ecutdg
   list_dpr(16)=hdr%ecutsm
   list_dpr(17)=hdr%tphysel
   list_dpr(18)=hdr%tsmear
   list_dpr(19:21)=hdr%qptn(1:3)                                 ; index=21
   list_dpr(1+index:3*nkpt +index)=reshape(hdr%kptns,(/3*nkpt/)) ; index=index+3*nkpt
   list_dpr(1+index:nkpt   +index)=hdr%wtk                       ; index=index+nkpt
   list_dpr(1+index:bantot +index)=hdr%occ                       ; index=index+bantot
   list_dpr(1+index:3*nsym +index)=reshape(hdr%tnons,(/3*nsym/)) ; index=index+3*nsym
   list_dpr(1+index:3*natom+index)=reshape(hdr%xred,(/3*natom/)) ; index=index+3*natom
   list_dpr(1+index:npsp   +index)=hdr%zionpsp                   ; index=index+npsp
   list_dpr(1+index:npsp   +index)=hdr%znuclpsp                  ; index=index+npsp
   list_dpr(1+index:ntypat  +index)=hdr%znucltypat               ; index=index+ntypat
 end if
 call MPI_BARRIER(spaceComm,ierr)
 call MPI_BCAST(list_dpr,list_size,MPI_DOUBLE_PRECISION,0,spaceComm,ierr)
 if(master/=me)then
   hdr%ecut_eff=list_dpr(1)
   hdr%etot    =list_dpr(2)
   hdr%fermie  =list_dpr(3)
   hdr%residm  =list_dpr(4)
   hdr%rprimd  =reshape(list_dpr(5:13),(/3,3/))
   hdr%ecut    =list_dpr(14)
   hdr%ecutdg  =list_dpr(15)
   hdr%ecutsm  =list_dpr(16)
   hdr%tphysel =list_dpr(17)
   hdr%tsmear  =list_dpr(18)
   hdr%qptn(1:3)=list_dpr(19:21)                                    ; index=21
   hdr%kptns   =reshape(list_dpr(1+index:3*nkpt +index),(/3,nkpt/)) ; index=index+3*nkpt
   hdr%wtk     =list_dpr(1+index:nkpt   +index)                     ; index=index+nkpt
   hdr%occ     =list_dpr(1+index:bantot +index)                     ; index=index+bantot
   hdr%tnons   =reshape(list_dpr(1+index:3*nsym +index),(/3,nsym/)) ; index=index+3*nsym
   hdr%xred    =reshape(list_dpr(1+index:3*natom+index),(/3,natom/)); index=index+3*natom
   hdr%zionpsp =list_dpr(1+index:npsp   +index)                     ; index=index+npsp
   hdr%znuclpsp=list_dpr(1+index:npsp   +index)                     ; index=index+npsp
   hdr%znucltypat=list_dpr(1+index:ntypat  +index)                  ; index=index+ntypat
 end if
 ABI_DEALLOCATE(list_dpr)

!Transmit the characters
 list_size=npsp+1
 ABI_ALLOCATE(list_char,(list_size))
 if (master==me)then
   list_char(1)       =hdr%codvsn  ! Only 6 characters are stored in list_char(1)
   list_char(2:npsp+1)=hdr%title
 end if
 call MPI_BARRIER(spaceComm,ierr)
 call MPI_BCAST(list_char,list_size*fnlen,MPI_CHARACTER,0,spaceComm,ierr)
 if(master/=me)then
   list_tmp=list_char(1)
   hdr%codvsn=list_tmp(1:6)
   do ipsp=2,npsp+1
     list_tmp =list_char(ipsp)
     hdr%title =list_tmp(1:fnlen)
   end do
 end if
 ABI_DEALLOCATE(list_char)

!Transmit the structured variables in case of PAW
 if (hdr%usepaw==1) then

   nrhoij=0
   if (master==me)then
     cplex=hdr%pawrhoij(1)%cplex
     nspden=hdr%pawrhoij(1)%nspden
     do iatom=1,natom
       nrhoij=nrhoij+hdr%pawrhoij(iatom)%nrhoijsel
     end do
   end if
   call MPI_BARRIER(spaceComm,ierr)
   call MPI_BCAST(nrhoij,1,MPI_INTEGER,0,spaceComm,ierr)
   call MPI_BCAST(cplex ,1,MPI_INTEGER,0,spaceComm,ierr)
   call MPI_BCAST(nspden,1,MPI_INTEGER,0,spaceComm,ierr)

   list_size=natom+nrhoij;list_size2=nspden*nrhoij*cplex
   ABI_ALLOCATE(list_int,(list_size))
   ABI_ALLOCATE(list_dpr,(list_size2))
   if (master==me)then
     index=0;index2=0
     do iatom=1,natom
       nsel=hdr%pawrhoij(iatom)%nrhoijsel
       list_int(1+index)=nsel
       list_int(2+index:1+nsel+index)=hdr%pawrhoij(iatom)%rhoijselect(1:nsel)
       index=index+1+nsel
       do ispden=1,nspden
         list_dpr(1+index2:nsel*cplex+index2)=hdr%pawrhoij(iatom)%rhoijp(1:nsel*cplex,ispden)
         index2=index2+nsel*cplex
       end do
     end do
   end if
   call MPI_BARRIER(spaceComm,ierr)
   call MPI_BCAST(list_int,list_size,MPI_INTEGER,0,spaceComm,ierr)
   call MPI_BCAST(list_dpr,list_size2,MPI_DOUBLE_PRECISION,0,spaceComm,ierr)
   if(master/=me)then
     index=0;index2=0
     ABI_ALLOCATE(hdr%pawrhoij,(natom))
     call rhoij_alloc(cplex,hdr%lmn_size,nspden,hdr%nspinor,hdr%nsppol,hdr%pawrhoij,hdr%typat)
     do iatom=1,natom
       nsel=list_int(1+index)
       hdr%pawrhoij(iatom)%nrhoijsel=nsel
       hdr%pawrhoij(iatom)%rhoijselect(1:nsel)=list_int(2+index:1+nsel+index)
       index=index+1+nsel
       do ispden=1,nspden
         hdr%pawrhoij(iatom)%rhoijp(1:nsel*cplex,ispden)=list_dpr(1+index2:nsel*cplex+index2)
         index2=index2+nsel*cplex
       end do
     end do
   end if
   ABI_DEALLOCATE(list_int)
   ABI_DEALLOCATE(list_dpr)

 end if

#else
!Dummy arguments
 if(me == -1 .and. master == -1 .and. spacecomm == 0 .and. hdr%bantot == 0) hdr%bantot = 0
#endif

!DEBUG
!write(std_out,*)' hdr_comm : exit '
!ENDDEBUG

end subroutine hdr_comm
!!***
