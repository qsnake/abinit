!{\src2tex{textfont=tt}}
!!****p* ABINIT/ujdet
!! NAME
!! ujdet
!!
!! FUNCTION
!!  Main routine. Determines U from inputfile ujdet.in (reduced inputfile containing mainly the atomic occupancies,
!!  atomic positions, and potential shifts.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2012 ABINIT group (DJA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! INPUTS
!! OUTPUT
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,destroy_mpi_enreg,herald,intagm,isfile,leave_new
!!      mpi_comm_rank,mpi_comm_size,nullify_mpi_enreg,parsefile,pawuj_det
!!      pawuj_free,pawuj_ini,pawuj_nullify,wrtout,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program ujdet

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use defs_parameters
 use m_xmpi
 use m_build_info
#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ujdet'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_42_parser
 use interfaces_51_manage_mpi
 use interfaces_57_iovars
 use interfaces_66_paw
!End of the abilint section

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Local variables-------------------------------
!scalars
!std_in=5, ab_in=5, std_out=6, ab_out=7, tmp_unit=9, tmp_unit2=10, fnlen=264, strlen=2000000
 integer,parameter     :: ndtpawuj=4,nwfchr=6
 integer               :: ierr,ii,jdtset,lenstr,marr,ndtset,tread
 integer               :: nat,nspden,macro_uj,pawujat,pawprtvol
 character(len=24)     :: codename
 character(len=30)     :: token
 character(len=fnlen)  :: filnam(5)
 character(len=strlen) :: string
 character(len=500)    :: message
 type(MPI_type)        :: mpi_enreg
 real(dp)              :: ures
!arrays
 integer,allocatable   :: idumar1(:),intarr(:),jdtset_(:)
 real(dp),allocatable  :: dpdumar(:),dprarr(:)
 type(macro_uj_type),allocatable   :: dtpawuj(:)

! *********************************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world,new_leave_comm=xmpi_world)

!Initialize MPI (not used but necessary)
 call xmpi_init()

!----------------------------------------------------------------------
!one should write a separate routine -init_mpi_enreg- for doing that !!
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
!inside ujdet.F90.

!----------------------------------------------------------------------

!Check that old input file exists (18seqpar/iofn1.F90)
 filnam(1)='ujdet.in'
 call isfile(filnam(1),'old')

!Check that new output file does NOT exist (18seqpar/iofn1.F90)
 filnam(2)='ujdet.out'
 call isfile(filnam(2),'new')

 open(unit=ab_out,file=trim(filnam(2)),form="formatted",status="new")
 rewind(unit=ab_out)

!Print header
 codename='UJDET'//repeat(' ',18)
 call herald(codename,abinit_version,ab_out)
 call herald(codename,abinit_version,std_out)

!Read the file, stringify it and return the number of datasets. (main/abinit.F90)
 call parsefile(filnam(1), lenstr, ndtset, string)

!DEBUG
!write(std_out,*)'ujdet: trim(string):',trim(adjustl(string))
!write(std_out,*)'ujdet: lenstr: ',lenstr
!write(std_out,*)'ujdet: ndtset: ',ndtset
!ENDDEBUG

 ABI_ALLOCATE(dtpawuj,(0:ndtpawuj))

 do jdtset=0,ndtpawuj
   call pawuj_nullify(dtpawuj(jdtset))
 end do

 call pawuj_ini(dtpawuj,ndtset)

 dtpawuj(1:ndtset)%ndtset=ndtset

 marr=ndtset
 ABI_ALLOCATE(idumar1,(marr))
 ABI_ALLOCATE(dpdumar,(marr))
 ABI_ALLOCATE(jdtset_,(ndtset))
 jdtset_=(/ ( ii,ii=1,ndtset )/)

!Read integers (main dimensions)
 do jdtset=1,ndtset
!  Integer
   token = 'pawprtvol'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%pawprtvol=idumar1(1)

   token = 'nat'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%nat=idumar1(1)

   token = 'nspden'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%nspden=idumar1(1)

   token = 'macro_uj'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%macro_uj=idumar1(1)

   token = 'pawujat'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%pawujat=idumar1(1)

   token = 'dmatpuopt'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%dmatpuopt=idumar1(1)

!  DEBUG
!  write(std_out,*)'ujdet: read pawujat ; jdtset, idumar1(1)',jdtset, idumar1(1)
!  write(std_out,*)'ujdet: read pawujat ; dtpawuj(1:ndtset)%pawujat ', dtpawuj(1:ndtset)%pawujat
!  END DEBUG

   token = 'pawujopt'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(1:ndtset)%option=idumar1(1)

!  Real

   token = 'diemix'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(1:ndtset)%diemix=idumar1(1)

   token = 'mdist'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(1:ndtset)%mdist=idumar1(1)

   token = 'pawujga'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(1:ndtset)%pawujga=dpdumar(1)

   token = 'ph0phiint'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(1:ndtset)%ph0phiint=dpdumar(1)

   token = 'pawrad'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
   if(tread==1) dtpawuj(1:ndtset)%pawrad=dpdumar(1)

   token = 'pawujrad'
   call intagm(dpdumar,idumar1,jdtset,marr,1,string(1:lenstr),token,tread,'LEN')
   if(tread==1) dtpawuj(1:ndtset)%pawujrad=dpdumar(1)

 end do !ndtset

 nat=dtpawuj(1)%nat
 nspden=dtpawuj(1)%nspden
 macro_uj=dtpawuj(1)%macro_uj
 pawujat=dtpawuj(1)%pawujat
 pawprtvol=dtpawuj(1)%pawprtvol

!Minimal tests
 if (.not.all(dtpawuj(1:ndtset)%nat==nat).and.all(dtpawuj(1:ndtset)%pawujat==pawujat)&
& .and.all(dtpawuj(1:ndtset)%macro_uj==macro_uj)) then
   write(message,'(a)')'Problems with nat, pawujat or nspden. Values for datasets differ.'
   call wrtout(std_out,message,'COLL')
   write(message,*) ' ujdet: nat ', all(dtpawuj(1:ndtset)%nat==nat)
   call wrtout(std_out,message,'COLL')
   write(message,*) ' ujdet: pawujat ', all(dtpawuj(1:ndtset)%pawujat==pawujat)
   call wrtout(std_out,message,'COLL')
   write(message,*) ' ujdet: macro_uj', all(dtpawuj(1:ndtset)%macro_uj==macro_uj)
   call wrtout(std_out,message,'COLL')
 end if

 if (pawprtvol>1) then
   write(message,fmt='(a,150i3)')' ujdet: dtpawuj(0:ndtset)%macro_uj ',dtpawuj(0:ndtset)%macro_uj
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150i3)')' ujdet: dtpawuj(0:ndtset)%pawujat ',dtpawuj(0:ndtset)%pawujat
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150i3)')' ujdet: dtpawuj(0:ndtset)%nspden ',dtpawuj(0:ndtset)%nspden
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150i3)')' ujdet: nat*nspden ',nat*nspden
   call wrtout(std_out,message,'COLL')
 end if

!Read arrays

 marr=maxval((/nspden*nat*nat, 3*3 ,nat*3/))
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 do jdtset=0,ndtset
!  DEBUG
!  write(std_out,*)'ujdet 2a, jdtset ',jdtset
!  END DEBUG
   ABI_ALLOCATE(dtpawuj(jdtset)%vsh,(nspden,nat))
   ABI_ALLOCATE(dtpawuj(jdtset)%occ,(nspden,nat))
   ABI_ALLOCATE(dtpawuj(jdtset)%xred,(3,nat))

   dtpawuj(jdtset)%iuj=jdtset
   dtpawuj(jdtset)%vsh=zero
   dtpawuj(jdtset)%occ=zero
   dtpawuj(jdtset)%xred=zero
 end do

 do jdtset=1,ndtset
   dprarr=-1_dp
   intarr=-1

!  Integer arrays

!  scdim, wfchr and rprimd allocated in pawuj_ini
   token = 'scdim'
   call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),token,tread,'INT')
   if(tread==1) dtpawuj(jdtset)%scdim(1:3)=intarr(1:3)

   token = 'wfchr'
   call intagm(dprarr,intarr,jdtset,marr,nwfchr,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%wfchr(1:nwfchr)=dprarr(1:nwfchr)

   write(std_out,*)'ujdet: wfchr ',dtpawuj(jdtset)%wfchr

   token = 'wfchr'
   call intagm(dprarr,intarr,jdtset,marr,nwfchr,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%wfchr(1:nwfchr)=dprarr(1:nwfchr)

   write(std_out,*)'ujdet: wfchr ',dtpawuj(jdtset)%wfchr

   token = 'rprimd'
   call intagm(dprarr,intarr,jdtset,marr,3*3,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%rprimd(1:3,1:3)=reshape(dprarr(1:3*3),(/ 3,3/))

   token = 'occ'
   call intagm(dprarr,intarr,jdtset,marr,nspden*nat,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%occ(1:nspden,1:nat)=reshape(dprarr(1:nspden*nat),(/ nspden,nat/))

   token = 'vsh'
   call intagm(dprarr,intarr,jdtset,marr,nat*nspden,string(1:lenstr),token,tread,'ENE')
   if(tread==1) dtpawuj(jdtset)%vsh(1:nspden,1:nat)=reshape(dprarr(1:nspden*nat),(/ nspden,nat/))

   token = 'xred'
   call intagm(dprarr,intarr,jdtset,marr,nat*3,string(1:lenstr),token,tread,'DPR')
   if(tread==1) dtpawuj(jdtset)%xred(1:3,1:nat)=reshape(dprarr(1:3*nat),(/ 3, nat/))

!  DEBUG
!  write(std_out,*)'dtpawuj(',jdtset,')%occ ',dtpawuj(jdtset)%occ,ch10
!  write(std_out,*)'dtpawuj(',jdtset,')%rprimd ',dtpawuj(jdtset)%rprimd,ch10
!  write(std_out,*)'dtpawuj(',jdtset,')%vsh ',dtpawuj(jdtset)%vsh,ch10
!  write(std_out,*)'dtpawuj(',jdtset,')%xred ',dtpawuj(jdtset)%xred,ch10,ch10
!  END DEBUG

 end do ! jdtset

 call pawuj_det(dtpawuj,ndtset, filnam(2)//"_UJDET.nc",ures)

!Close files
 close(ab_out)

!Deallocations
 do jdtset=0,ndtset
   call pawuj_free(dtpawuj(jdtset))
 end do
 ABI_DEALLOCATE(dtpawuj)

 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(idumar1)
 ABI_DEALLOCATE(dpdumar)
 ABI_DEALLOCATE(jdtset_)

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program ujdet
!!***

