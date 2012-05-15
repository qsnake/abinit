!{\src2tex{textfont=tt}}
!! ===================================================
!! This module contains I/O and MPI methods
!! used to manipulate variables of structured
!! datatype cprj_type
!!
!! cprj_type variables are <p_lmn|Cnk> projected
!! quantities where |p_lmn> are non-local projectors
!!                  |Cnk> are wave functions
!! ===================================================

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!!****f* ABINIT/cprj_diskinit_r
!! NAME
!! cprj_diskinit_r
!!
!! FUNCTION
!! Initialize a cprj temporary file for READING
!! Nothing is done if mkmem=0
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atind(natom)=index table for atoms (see iorder below)
!!  dimcp=first dimension of cprj arrays (1 or natom)
!!  iorder=0 if cprj ordering does not change during reading
!!         1 if cprj ordering changes during writing, depending on content of atind array:
!!              - if atind=atindx  (unsorted->type-sorted)
!!              - if atind=atindx1 (type-sorted->unsorted)
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  natom=number of atoms in cell
!!  ncpgr=number of gradients of cprj
!!  nlmn(dimcp)=array of dimensions of cprj datastructure that will contain the read data
!!  nspinor=number of spinorial components of the wavefunctions
!!  uncp=unit number for cprj data (if used)
!!
!! PARENTS
!!      cprj_utils_mpi,datafordmft,dyfnl3,nstpaw3,optics_paw,optics_paw_core
!!      partial_dos_fractions_paw,pawmkaewf,pawmkrhoij,suscep_stat,vtorho3
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

 subroutine cprj_diskinit_r(atind,dimcp,iorder,mkmem,natom,ncpgr,nlmn,nspinor,uncp)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_diskinit_r'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimcp,iorder,mkmem,natom,ncpgr,nspinor,uncp
!arrays
 integer,intent(in) :: atind(natom),nlmn(dimcp)
!Local variables-------------------------------
 integer :: dimcp0,iatm,iatom,ncpgr0,nspinor0
 character(len=500) :: message
 integer,allocatable :: dimlmn(:)

! *************************************************************************

 if (mkmem==0) then

   rewind uncp;read(uncp) dimcp0,ncpgr0,nspinor0
   if (dimcp/=dimcp0.or.ncpgr/=ncpgr0.or.nspinor/=nspinor0) then
     write(message,'(a,a,a,a)')ch10,&
&     ' cprj_diskinit_r : BUG -',ch10,&
&     '  _PAW file was not created with the right options !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   ABI_ALLOCATE(dimlmn,(dimcp))
   read(uncp) dimlmn(1:dimcp)
   do iatom=1,dimcp
     if (iorder==0) then
       iatm=iatom
     else
       iatm=min(atind(iatom),dimcp)
     end if
     if (dimlmn(iatom)/=nlmn(iatm)) then
       write(message,'(a,a,a,a)')ch10,&
&       ' cprj_diskinit_r : BUG -',ch10,&
&       '  _PAW file was not created with the right options !'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
   end do
   ABI_DEALLOCATE(dimlmn)

 end if

end subroutine cprj_diskinit_r
!!***


!!****f* ABINIT/cprj_diskinit_w
!! NAME
!! cprj_diskinit_w
!!
!! FUNCTION
!! Initialize a cprj temporary file for WRITING
!! Nothing is done if mkmem=0
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atind(natom)=index table for atoms (see iorder below)
!!  dimcp=first dimension of cprj arrays (1 or natom)
!!  iorder=0 if cprj ordering does not change during reading
!!         1 if cprj ordering changes during writing, depending on content of atind array:
!!              - if atind=atindx  (type-sorted->unsorted)
!!              - if atind=atindx1 (unsorted->type-sorted)
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  natom=number of atoms in cell
!!  ncpgr=number of gradients of cprj
!!  nlmn(dimcp)=array of dimensions of cprj datastructure that will contain the read data
!!  nspinor=number of spinorial components of the wavefunctions
!!  uncp=unit number for cprj data (if used)
!!
!! PARENTS
!!      ctocprj,vtorho,vtorho3
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

 subroutine cprj_diskinit_w(atind,dimcp,iorder,mkmem,natom,ncpgr,nlmn,nspinor,uncp)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_diskinit_w'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimcp,iorder,mkmem,natom,ncpgr,nspinor,uncp
!arrays
 integer,intent(in) :: atind(natom),nlmn(dimcp)
!Local variables-------------------------------
 integer :: iatm,iatom
 integer,allocatable :: dimlmn(:)

! *************************************************************************

 if (mkmem==0) then

   rewind uncp
   write(uncp) dimcp,ncpgr,nspinor

   if (iorder==0) then
     write(uncp) nlmn(1:dimcp)
   else
     ABI_ALLOCATE(dimlmn,(dimcp))
     do iatom=1,dimcp
       iatm=min(atind(iatom),dimcp)
       dimlmn(iatom)=nlmn(iatm)
     end do
     write(uncp) dimlmn(1:dimcp)
     ABI_DEALLOCATE(dimlmn)
   end if

 end if

end subroutine cprj_diskinit_w
!!***


!!****f* ABINIT/cprj_diskskip
!! NAME
!! cprj_diskskip
!!
!! FUNCTION
!! Skip records in a cprj temporary file for READING
!! Nothing is done if mkmem=0
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  ncpgr=number of gradients of cprj
!!  nres=number of records to be skipped
!!  uncp=unit number for cprj data (if used)
!!
!! PARENTS
!!      vtowfk3
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

 subroutine cprj_diskskip(mkmem,ncpgr,nrec,uncp)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_diskskip'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mkmem,ncpgr,nrec,uncp
!arrays
!Local variables-------------------------------
 integer :: ii

! *************************************************************************

 if (mkmem==0) then

   do ii=1,nrec
     read(uncp)
     if (ncpgr>0) read(uncp)
   end do

 end if

end subroutine cprj_diskskip
!!***


!!****f* ABINIT/cprj_get
!! NAME
!! cprj_get
!!
!! FUNCTION
!! Read the cprj for a given k-point from memory or from a temporary file
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atind(natom)=index table for atoms (see iorder below)
!!  cprj(dimcp,nspinor*mband*mkmem*nsppol)=input cprj (used if mkmem/=0)
!!  dimcp=first dimension of cprj_k,cprj arrays (1 or natom)
!!  iband1=index of first band
!!  ibg=shift if cprj array to locate current k-point
!!  icpgr= (optional argument) if present, only component icpgr of
!!         input cprj gradient is copied into output cprj
!!         Not used if cprj(:,:)%ncpgr<icpgr (mkmem>0)
!!                  or ncpgr(optional)<icpgr (mkmem=0)
!!  ikpt=index of current k-point
!!  iorder=0 if cprj ordering does not change during reading
!!         1 if cprj ordering changes during reading, depending on content of atind array:
!!              - if atind=atindx  (type-sorted=>unsorted)
!!              - if atind=atindx1 (unsorted=>type-sorted)
!!  isppol=index of current spin component
!!  mband=maximum number of bands
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands to import (usually 1 or nband_k)
!!  nband_k=total number of bands for this k-point
!!  ncpgr=(optional argument) second dimension of cprj%dcp(2,ncpgr,nlmn
!!        stored in memory (mkmem>0) or present on disk (mkmem=0))
!!        needed only when optional argument icpgr is present
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  uncp=unit number for cprj data (used if mkmem=0)
!!
!! OUTPUT
!!  cprj_k(dimcp,nspinor*nband) <type(cprj_type)>= output cprj datastructure
!!
!! PARENTS
!!      berry_linemin,berryphase_new,cgwf,cprj_utils_mpi,dyfnl3,extrapwf
!!      mag_loc_k,make_grad_berry,nstpaw3,optics_paw,optics_paw_core
!!      partial_dos_fractions_paw,pawmkaewf,pawmkrhoij,smatrix_pawinit
!!      suscep_stat,update_mmat,vtorho,vtorho3,vtowfk3
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

 subroutine cprj_get(atind,cprj_k,cprj,dimcp,iband1,ibg,ikpt,iorder,isppol,mband,&
&                    mkmem,mpi_enreg,natom,nband,nband_k,nspinor,nsppol,uncp,&
&                    icpgr,ncpgr) ! optionals arguments

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_get'
 use interfaces_51_manage_mpi, except_this_one => cprj_get
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimcp,iband1,ibg,ikpt,iorder,isppol,mband,mkmem,natom,nband,nband_k,nspinor,nsppol,uncp
 integer,intent(in),optional :: icpgr,ncpgr
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atind(natom)
 type(cprj_type),intent(in) :: cprj(dimcp,nspinor*mband*mkmem*nsppol)
 type(cprj_type),intent(inout) :: cprj_k(dimcp,nspinor*nband)

!Local variables-------------------------------
!scalars
 integer :: iatm,iatom,ib,ibsp,icpgr_,isp,ispinor,jband,me,nband0,ncpgr_
 logical :: has_icpgr
 character(len=500) :: message
!arrays
 real(dp),allocatable :: tmp(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 ncpgr_=cprj_k(1,1)%ncpgr;if (present(ncpgr)) ncpgr_=ncpgr
 icpgr_=-1;if(present(icpgr)) icpgr_=icpgr
 has_icpgr=(icpgr_>0.and.icpgr_<=ncpgr_)
 if (present(icpgr).and.(.not.present(ncpgr))) then
   message='  ncpgr must be present when icpgr is present !'
   MSG_BUG(message)
 end if
 if (has_icpgr.and.cprj_k(1,1)%ncpgr<1) then
   message='  cprj_k%ncpgr not consistent with icpgr !'
   MSG_BUG(message)
 end if

 call xme_init(mpi_enreg,me)

 if (mkmem==0) then

   if (iband1==1) then
     read(uncp) nband0
     if (nband_k/=nband0) then
       message='  _PAW file was not created with the right options !'
       MSG_BUG(message)
     end if
   end if

   isp=0;jband=iband1-1
   do ib=1,nband
     jband=jband+1
     if (mpi_enreg%paral_compil_kpt==1) then
       if (abs(mpi_enreg%proc_distrb(ikpt,jband,isppol)-me)/=0) then
         isp=isp+nspinor
         cycle
       end if
     end if
     do ispinor=1,nspinor
       isp=isp+1
       if (iorder==0) then
         if (ncpgr_==0) then
           do iatom=1,dimcp
             read(uncp) cprj_k(iatom,isp)%cp(:,:)
           end do
         else
           if (has_icpgr) then
             do iatom=1,dimcp
               ABI_ALLOCATE(tmp,(2,ncpgr_,cprj_k(iatom,1)%nlmn))
               read(uncp) cprj_k(iatom,isp)%cp(:,:),tmp(:,:,:)
               cprj_k(iatom,isp)%dcp(:,1,:)=tmp(:,icpgr_,:)
               ABI_DEALLOCATE(tmp)
             end do
           else
             do iatom=1,dimcp
               read(uncp) cprj_k(iatom,isp)%cp(:,:),cprj_k(iatom,isp)%dcp(:,:,:)
             end do
           end if
         end if
       else
         if (ncpgr_==0) then
           do iatom=1,dimcp
             iatm=min(atind(iatom),dimcp)
             read(uncp) cprj_k(iatm,isp)%cp(:,:)
           end do
         else
           if (has_icpgr) then
             do iatom=1,dimcp
               iatm=min(atind(iatom),dimcp)
               ABI_ALLOCATE(tmp,(2,ncpgr_,cprj_k(iatm,1)%nlmn))
               read(uncp) cprj_k(iatm,isp)%cp(:,:),tmp(:,:,:)
               cprj_k(iatm,isp)%dcp(:,1,:)=tmp(:,icpgr_,:)
               ABI_DEALLOCATE(tmp)
             end do
           else
             do iatom=1,dimcp
               iatm=min(atind(iatom),dimcp)
               read(uncp) cprj_k(iatm,isp)%cp(:,:),cprj_k(iatm,isp)%dcp(:,:,:)
             end do
           end if
         end if
       end if
     end do
   end do

 else

   isp=0;ibsp=ibg+nspinor*(iband1-1);jband=iband1-1
   do ib=1,nband
     jband=jband+1
     if (mpi_enreg%paral_compil_kpt==1) then
       if (abs(mpi_enreg%proc_distrb(ikpt,jband,isppol)-me)/=0) then
         isp=isp+nspinor;ibsp=ibsp+nspinor
         cycle
       end if
     end if
     do ispinor=1,nspinor
       isp=isp+1;ibsp=ibsp+1
       if (iorder==0) then
         if (ncpgr_==0) then
           do iatom=1,dimcp
             cprj_k(iatom,isp)%cp(:,:)=cprj(iatom,ibsp)%cp(:,:)
           end do
         else
           if (has_icpgr) then
             do iatom=1,dimcp
               cprj_k(iatom,isp)%cp(:,:)   =cprj(iatom,ibsp)%cp(:,:)
               cprj_k(iatom,isp)%dcp(:,1,:)=cprj(iatom,ibsp)%dcp(:,icpgr_,:)
             end do
           else
             do iatom=1,dimcp
               cprj_k(iatom,isp)%cp(:,:)   =cprj(iatom,ibsp)%cp(:,:)
               cprj_k(iatom,isp)%dcp(:,:,:)=cprj(iatom,ibsp)%dcp(:,:,:)
             end do
           end if
         end if
       else
         if (ncpgr_==0) then
           do iatom=1,dimcp
             iatm=min(atind(iatom),dimcp)
             cprj_k(iatm,isp)%cp(:,:)=cprj(iatom,ibsp)%cp(:,:)
           end do
         else
           if (has_icpgr) then
             do iatom=1,dimcp
               iatm=min(atind(iatom),dimcp)
               cprj_k(iatm,isp)%cp(:,:)   =cprj(iatom,ibsp)%cp(:,:)
               cprj_k(iatm,isp)%dcp(:,1,:)=cprj(iatom,ibsp)%dcp(:,icpgr_,:)
             end do
           else
             do iatom=1,dimcp
               iatm=min(atind(iatom),dimcp)
               cprj_k(iatm,isp)%cp(:,:)   =cprj(iatom,ibsp)%cp(:,:)
               cprj_k(iatm,isp)%dcp(:,:,:)=cprj(iatom,ibsp)%dcp(:,:,:)
             end do
           end if
         end if
       end if
     end do
   end do

 end if

 DBG_EXIT("COLL")

end subroutine cprj_get
!!***


!!****f* ABINIT/cprj_put
!! NAME
!! cprj_put
!!
!! FUNCTION
!! Write the cprj for a given set of (n,k) into memory or into a temporary file
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atind(natom)=index table for atoms (see iorder below)
!!  cprj_k(dimcp,nspinor*nband) <type(cprj_type)>= input cprj datastructure
!!  dimcp=first dimension of cprj_k,cprjnk arrays (1 or natom)
!!  iband1=index of first band
!!  ibg=shift if cprjnk array to locate current k-point
!!  ikpt=index of current k-point
!!  iorder=0 if cprj ordering does not change during reading
!!         1 if cprj ordering changes during writing, depending on content of atind array:
!!              - if atind=atindx  (type-sorted->unsorted)
!!              - if atind=atindx1 (unsorted->type-sorted)
!!  isppol=index of current spin component
!!  mband=maximum number of bands
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands to export (usually 1, nband_k or nblockbd)
!!  nband_k=total number of bands for this k-point
!!  nlmn(dimcp)=array of dimensions of cprj_k,cprjnk datastructures
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  spaceComm_band=communicator used for bands in case of band-fft parallelism
!!  to_be_gathered=TRUE if cprj_k arrays have to gathered between procs,
!!                 (band-fft parallelism only) - OPTIONAL
!!  uncp=unit number for cprj data (used if mkmem=0)
!!
!! SIDE EFFECTS
!!  cprj(dimcp,nspinor*mband*mkmem*nsppol)=output cprj (used if mkmem/=0)
!!
!! PARENTS
!!      berry_linemin,berryphase_new,cgwf,cprj_utils_mpi,ctocprj,extrapwf
!!      vtorho,vtowfk,vtowfk3
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

 subroutine cprj_put(atind,cprj_k,cprj,dimcp,iband1,ibg,ikpt,iorder,isppol,mband,&
&           mkmem,mpi_enreg,natom,nband,nband_k,nlmn,nspinor,nsppol,spaceComm_band,uncp,&
&           to_be_gathered) ! Optional argument

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_put'
 use interfaces_51_manage_mpi, except_this_one => cprj_put
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband1,ibg,ikpt,iorder,isppol,dimcp,mband,mkmem
 integer,intent(in) :: natom,nband,nband_k,nspinor,nsppol,uncp,spaceComm_band
 logical,optional,intent(in) :: to_be_gathered
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer :: atind(natom),nlmn(dimcp)
 type(cprj_type),intent(out) :: cprj(dimcp,nspinor*mband*mkmem*nsppol)
 type(cprj_type),intent(in) :: cprj_k(dimcp,nspinor*nband)
!Local variables-------------------------------
 integer :: iatm,iatom,iband,ibsp,icpgr,ierr,ii,ilmn,isp,ispinor,jband,jj,lmndim,me,ncpgr
 logical :: to_be_gathered_
 real(dp),allocatable :: buffer1(:),buffer2(:)

! *************************************************************************

 ncpgr=cprj_k(1,1)%ncpgr

 to_be_gathered_=.false.;if (present(to_be_gathered)) to_be_gathered_=to_be_gathered

 call xme_init(mpi_enreg,me)

 if (mpi_enreg%mode_para/='b'.or.mpi_enreg%nproc_band==1.or.(.not.to_be_gathered_)) then

   if (mkmem==0) then

     if (iband1==1) write(uncp) nband_k

     isp=0;jband=iband1-1
     do iband=1,nband
       jband=jband+1
       if (mpi_enreg%paral_compil_kpt==1) then
         if (abs(mpi_enreg%proc_distrb(ikpt,jband,isppol)-me)/=0) then
           isp=isp+nspinor
           cycle
         end if
       end if
       do ispinor=1,nspinor
         isp=isp+1
         if (iorder==0) then
           do iatom=1,dimcp
             if (ncpgr==0) then
               write(uncp) cprj_k(iatom,isp)%cp(:,:)
             else
               write(uncp) cprj_k(iatom,isp)%cp(:,:),cprj_k(iatom,isp)%dcp(:,:,:)
             end if
           end do
         else
           do iatom=1,dimcp
             iatm=min(atind(iatom),dimcp)
             if (ncpgr==0) then
               write(uncp) cprj_k(iatm,isp)%cp(:,:)
             else
               write(uncp) cprj_k(iatm,isp)%cp(:,:),cprj_k(iatm,isp)%dcp(:,:,:)
             end if
           end do
         end if
       end do
     end do

   else

     isp=0;ibsp=ibg+nspinor*(iband1-1);jband=iband1-1
     do iband=1,nband
       jband=jband+1
       if (mpi_enreg%paral_compil_kpt==1)then
         if (abs(mpi_enreg%proc_distrb(ikpt,jband,isppol)-me)/=0) then
           isp=isp+nspinor;ibsp=ibsp+nspinor
           cycle
         end if
       end if
       do ispinor=1,nspinor
         isp=isp+1;ibsp=ibsp+1
         if (iorder==0) then
           do iatom=1,dimcp
             cprj(iatom,ibsp)%cp(:,:)=cprj_k(iatom,isp)%cp(:,:)
             if (ncpgr>0) cprj(iatom,ibsp)%dcp(:,:,:)=cprj_k(iatom,isp)%dcp(:,:,:)
           end do
         else
           do iatom=1,dimcp
             iatm=min(atind(iatom),dimcp)
             cprj(iatom,ibsp)%cp(:,:)=cprj_k(iatm,isp)%cp(:,:)
             if (ncpgr>0) cprj(iatom,ibsp)%dcp(:,:,:)=cprj_k(iatm,isp)%dcp(:,:,:)
           end do
         end if
       end do
     end do

   end if

 else ! mode_para==b and nband>1

   lmndim=2*sum(nlmn(1:dimcp))*(1+ncpgr)*nspinor
   ABI_ALLOCATE(buffer1,(lmndim))
   ABI_ALLOCATE(buffer2,(lmndim*mpi_enreg%nproc_band))
   isp=0;ibsp=ibg+nspinor*(iband1-1)
   do iband=1,nband  ! must be nblockbd for band-fft parallelism
     jj=1
     do ispinor=1,nspinor
       isp=isp+1
       do iatom=1,dimcp
         if (iorder==0) then
           iatm=iatom
         else
           iatm=min(atind(iatom),dimcp)
         end if
         do ilmn=1,nlmn(iatm)
           buffer1(jj:jj+1)=cprj_k(iatm,isp)%cp(1:2,ilmn)
           jj=jj+2
         end do
         if (ncpgr>0) then
           do ilmn=1,nlmn(iatm)
             do icpgr=1,ncpgr
               buffer1(jj:jj+1)=cprj_k(iatm,isp)%dcp(1:2,icpgr,ilmn)
               jj=jj+2
             end do
           end do
         end if
       end do
     end do
     call xallgather_mpi(buffer1,lmndim,buffer2,spaceComm_band,ierr)
     jj=1
     do ii=1,mpi_enreg%nproc_band
       do ispinor=1,nspinor
         ibsp=ibsp+1
         do iatom=1,dimcp
           if (iorder==0) then
             iatm=iatom
           else
             iatm=min(atind(iatom),dimcp)
           end if
           do ilmn=1,nlmn(iatm)
             cprj(iatom,ibsp)%cp(1:2,ilmn)=buffer2(jj:jj+1)
             jj=jj+2
           end do
         end do
         if (ncpgr>0) then
           do ilmn=1,nlmn(iatm)
             do icpgr=1,ncpgr
               cprj(iatom,ibsp)%dcp(1:2,icpgr,ilmn)=buffer2(jj:jj+1)
               jj=jj+2
             end do
           end do
         end if
       end do
     end do
   end do
   ABI_DEALLOCATE(buffer1)
   ABI_DEALLOCATE(buffer2)

 end if ! mode_para=b, nband

end subroutine cprj_put
!!***

!!****f* ABINIT/cprj_exch
!! NAME
!! cprj_exch
!!
!! FUNCTION
!! Exchange a cprj_type between two processors inside a MPI communicator.
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom=Number of atoms (size of first dimension of Cprj_send and Cprj_recv).
!!  n2dim=Size of the second dimension.
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  Cprj_send= The datatype to be transmitted.
!!  receiver=ID of the receiver in spaceComm.
!!  sender=ID of the sender in spaceComm.
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  ierr=Error status.
!!  Cprj_recv=The datatype copied on proc. receiver.
!!
!! NOTES
!!  If sender==receiver, Cprj_send is copied into Cprj_recv.
!!  It should be easy to avoid this additional copy in the calling routine.
!!
!! PARENTS
!!      outkss
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

subroutine cprj_exch(natom,n2dim,nlmn,ncpgr,Cprj_send,Cprj_recv,sender,receiver,spaceComm,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi

 use m_errors, only : assert

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_exch'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr
 integer,intent(in) :: sender,receiver,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(cprj_type),intent(in) :: Cprj_send(:,:)
 type(cprj_type),intent(inout) :: Cprj_recv(:,:)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,t2dim,tcpgr,n1dim,nn
 integer :: ntotcp,ipck,rank
 character(len=500) :: message
!arrays
 real(dp),allocatable :: buffer_cp(:,:),buffer_cpgr(:,:,:)

! *************************************************************************

 n1dim=0
 t2dim=0
 tcpgr=0
 ierr=0
 if (sender==receiver) then
   call cprj_copy(Cprj_send,Cprj_recv)
   return
 end if

 rank = xcomm_rank(spaceComm)

!#if defined DEBUG_MODE
 nn=size(nlmn,dim=1)
 if (rank==sender) then
   n1dim=size(Cprj_send,dim=1)
   t2dim=size(Cprj_send,dim=2)
   tcpgr=Cprj_send(1,1)%ncpgr
 end if
 if (rank==receiver) then
   n1dim=size(Cprj_recv,dim=1)
   t2dim=size(Cprj_recv,dim=2)
   tcpgr=Cprj_recv(1,1)%ncpgr
 end if
 if (rank/=sender.and.rank/=receiver) then
   write(message, '(a,a,3i5,a,a)' ) ch10,&
&   ' BUG: cprj_exch: rank is not equal to sender or receiver ', &
&   rank, sender, receiver,ch10,&
&   '   Action: contact abinit group'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 call assert(   (nn==n1dim),'cprj_exch: size mismatch in natom!')
 call assert((t2dim==n2dim),'cprj_exch: size mismatch in dim=2!')
 call assert((tcpgr==ncpgr),'cprj_exch: size mismatch in ncpgr!')
!#endif

 ntotcp=n2dim*SUM(nlmn(:))

 ABI_ALLOCATE(buffer_cp,(2,ntotcp))
 if (ncpgr/=0)  then
   ABI_ALLOCATE(buffer_cpgr,(2,ncpgr,ntotcp))
 end if

!=== Pack Cprj_send ===
 if (rank==sender) then
   ipck=0
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       buffer_cp(:,ipck+1:ipck+nn)=Cprj_send(iat,jj)%cp(:,1:nn)
       if (ncpgr/=0) buffer_cpgr(:,:,ipck+1:ipck+nn)=Cprj_send(iat,jj)%dcp(:,:,1:nn)
       ipck=ipck+nn
     end do
   end do
 end if

!=== Transmit data ===
 call xexch_mpi(buffer_cp,2*ntotcp,sender,buffer_cp,receiver,spaceComm,ierr)
 if (ncpgr/=0) then
   call xexch_mpi(buffer_cpgr,2*ncpgr*ntotcp,sender,buffer_cpgr,receiver,spaceComm,ierr)
 end if

!=== UnPack buffers into Cprj_recv ===
 if (rank==receiver) then
   ipck=0
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       Cprj_recv(iat,jj)%cp(:,1:nn)=buffer_cp(:,ipck+1:ipck+nn)
       if (ncpgr/=0) Cprj_recv(iat,jj)%dcp(:,:,1:nn)=buffer_cpgr(:,:,ipck+1:ipck+nn)
       ipck=ipck+nn
     end do
   end do
 end if

 ABI_DEALLOCATE(buffer_cp)
 if (ncpgr/=0)  then
   ABI_DEALLOCATE(buffer_cpgr)
 end if

end subroutine cprj_exch
!!***

!!****f* ABINIT/cprj_mpi_send
!! NAME
!! cprj_mpi_send
!!
!! FUNCTION
!! Send a cprj_type inside a MPI communicator.
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom=Number of atoms (size of first dimension of cprj_out).
!!  n2dim=Size of the second dimension.
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  ncpgr = number of gradients in cprj_out
!!  cprj_out= The datatype to be transmitted.
!!  receiver=ID of the receiver in spaceComm.
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  ierr=Error status.
!!
!! NOTES
!!   perhaps in general it is more efficient to use cprj_exch but it is
!!   convenient for coding to have separate send and recieve routines.
!!
!! PARENTS
!!      berryphase_new
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

subroutine cprj_mpi_send(natom,n2dim,nlmn,ncpgr,cprj_out,receiver,spaceComm,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi

 use m_errors, only : assert

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_mpi_send'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr
 integer,intent(in) :: receiver,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(cprj_type),intent(in) :: cprj_out(:,:)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,t2dim,tcpgr,n1dim,nn
 integer :: ntotcp,ipck,tag
!arrays
 real(dp),allocatable :: buffer_cp(:,:),buffer_cpgr(:,:,:)

! *************************************************************************

 n1dim=0
 t2dim=0
 tcpgr=0
 ierr=0

 nn=size(nlmn,dim=1)
 n1dim=size(cprj_out,dim=1)
 t2dim=size(cprj_out,dim=2)
 tcpgr=cprj_out(1,1)%ncpgr

 call assert(   (nn==n1dim),'cprj_mpi_send: size mismatch in natom!')
 call assert((t2dim==n2dim),'cprj_mpi_send: size mismatch in dim=2!')
 call assert((tcpgr==ncpgr),'cprj_mpi_send: size mismatch in ncpgr!')

 ntotcp=n2dim*SUM(nlmn(:))

 ABI_ALLOCATE(buffer_cp,(2,ntotcp))
 if (ncpgr/=0)  then
   ABI_ALLOCATE(buffer_cpgr,(2,ncpgr,ntotcp))
 end if

!=== Pack cprj_out ====
 ipck=0
 do jj=1,n2dim
   do iat=1,natom
     nn=nlmn(iat)
     buffer_cp(:,ipck+1:ipck+nn)=cprj_out(iat,jj)%cp(:,1:nn)
     if (ncpgr/=0) buffer_cpgr(:,:,ipck+1:ipck+nn)=cprj_out(iat,jj)%dcp(:,:,1:nn)
     ipck=ipck+nn
   end do
 end do

!=== Transmit data ===
 tag = 2*ntotcp
 call xsend_mpi(buffer_cp,receiver,tag,spaceComm,ierr)
 if (ncpgr/=0) then
   tag=tag*ncpgr
   call xsend_mpi(buffer_cpgr,receiver,tag,spaceComm,ierr)
 end if

!=== Clean up ===
 ABI_DEALLOCATE(buffer_cp)
 if (ncpgr/=0)  then
   ABI_DEALLOCATE(buffer_cpgr)
 end if

end subroutine cprj_mpi_send
!!***

!!****f* ABINIT/cprj_mpi_recv
!! NAME
!! cprj_mpi_recv
!!
!! FUNCTION
!! Receive a cprj_type inside a MPI communicator.
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom=Number of atoms (size of first dimension of Cprj_in).
!!  n2dim=Size of the second dimension.
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  ncpgr = number of gradients in cprj_in
!!  sender=ID of the sender in spaceComm.
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  ierr=Error status.
!!  cprj_in=The datatype copied on proc. receiver.
!!
!! NOTES
!!   perhaps in general it is more efficient to use cprj_exch but it is
!!   convenient for coding to have separate send and receive routines.
!!
!! PARENTS
!!      berryphase_new
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

subroutine cprj_mpi_recv(natom,n2dim,nlmn,ncpgr,cprj_in,sender,spaceComm,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi

 use m_errors, only : assert

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_mpi_recv'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr
 integer,intent(in) :: sender,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(cprj_type),intent(inout) :: cprj_in(:,:)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,t2dim,tcpgr,n1dim,nn
 integer :: ntotcp,ipck,tag
!arrays
 real(dp),allocatable :: buffer_cp(:,:),buffer_cpgr(:,:,:)

! *************************************************************************

 n1dim=0
 t2dim=0
 tcpgr=0
 ierr=0

 nn=size(nlmn,dim=1)
 n1dim=size(cprj_in,dim=1)
 t2dim=size(cprj_in,dim=2)
 tcpgr=cprj_in(1,1)%ncpgr

 call assert(   (nn==n1dim),'cprj_recv: size mismatch in natom!')
 call assert((t2dim==n2dim),'cprj_recv: size mismatch in dim=2!')
 call assert((tcpgr==ncpgr),'cprj_recv: size mismatch in ncpgr!')

 ntotcp=n2dim*SUM(nlmn(:))

 ABI_ALLOCATE(buffer_cp,(2,ntotcp))
 if (ncpgr/=0)  then
   ABI_ALLOCATE(buffer_cpgr,(2,ncpgr,ntotcp))
 end if

!=== Receive data ===
 tag = 2*ntotcp
 call xrecv_mpi(buffer_cp,sender,tag,spaceComm,ierr)
 if (ncpgr/=0) then
   tag=tag*ncpgr
   call xrecv_mpi(buffer_cpgr,sender,tag,spaceComm,ierr)
 end if

!=== UnPack buffers into cprj_in ===
 ipck=0
 do jj=1,n2dim
   do iat=1,natom
     nn=nlmn(iat)
     cprj_in(iat,jj)%cp(:,1:nn)=buffer_cp(:,ipck+1:ipck+nn)
     if (ncpgr/=0) cprj_in(iat,jj)%dcp(:,:,1:nn)=buffer_cpgr(:,:,ipck+1:ipck+nn)
     ipck=ipck+nn
   end do
 end do

!=== Clean up ===
 ABI_DEALLOCATE(buffer_cp)
 if (ncpgr/=0)  then
   ABI_DEALLOCATE(buffer_cpgr)
 end if

end subroutine cprj_mpi_recv
!!***

!!****f* ABINIT/cprj_mpi_allgather
!! NAME
!! cprj_mpi_allgather
!!
!! FUNCTION
!! Perform MPI_ALLGATHER on a cprj_type inside a MPI communicator.
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprj_loc= The cprj on the local proc being all-gathered
!!  natom=Number of atoms (size of first dimension of cprj_loc).
!!  n2dim=Size of the second dimension of cprj_loc.
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  ncpgr = number of gradients in cprj_loc
!!  nproc=number of processors being gathered
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  cprj_gat=the gathered cprjs
!!  ierr=Error status.
!!
!! NOTES
!!
!! PARENTS
!!      berryphase_new,cgwf,suscep_stat,vtorho
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

subroutine cprj_mpi_allgather(cprj_loc,cprj_gat,natom,n2dim,nlmn,ncpgr,nproc,spaceComm,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi

 use m_errors, only : assert

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_mpi_allgather'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr
 integer,intent(in) :: nproc,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(cprj_type),intent(in) :: cprj_loc(:,:)
 type(cprj_type),intent(out) :: cprj_gat(:,:)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,t2dim,tcpgr,tg2dim,n1dim,nn
 integer :: ntotcp,ipck,iproc
!arrays
 real(dp),allocatable :: buffer_cpgr(:,:,:),buffer_cpgr_all(:,:,:)

! *************************************************************************

 n1dim=0
 t2dim=0
 tg2dim=0
 tcpgr=0
 ierr=0

 nn=size(nlmn,dim=1)
 n1dim=size(cprj_loc,dim=1)
 t2dim=size(cprj_loc,dim=2)
 tg2dim=size(cprj_gat,dim=2)
 tcpgr=cprj_loc(1,1)%ncpgr

 call assert(   (nn==n1dim),'cprj_mpi_allgather: size mismatch in natom!')
 call assert((t2dim==n2dim),'cprj_mpi_allgather: size mismatch in dim=2!')
 call assert((tg2dim==n2dim*nproc),'cprj_mpi_allgather: size mismatch in dim=2!')
 call assert((tcpgr==ncpgr),'cprj_mpi_allgather: size mismatch in ncpgr!')

 ntotcp=n2dim*SUM(nlmn(:))

 ABI_ALLOCATE(buffer_cpgr,(2,1+ncpgr,ntotcp))
 ABI_ALLOCATE(buffer_cpgr_all,(2,1+ncpgr,nproc*ntotcp))

!=== Pack cprj_loc ====
 ipck=0
 do jj=1,n2dim
   do iat=1,natom
     nn=nlmn(iat)
     buffer_cpgr(:,1,ipck+1:ipck+nn)=cprj_loc(iat,jj)%cp(:,1:nn)
     if (ncpgr/=0) buffer_cpgr(:,2:1+ncpgr,ipck+1:ipck+nn)=cprj_loc(iat,jj)%dcp(:,:,1:nn)
     ipck=ipck+nn
   end do
 end do

!=== allgather data ===
 call xallgather_mpi(buffer_cpgr,2*(ncpgr+1)*ntotcp,buffer_cpgr_all,spaceComm,ierr)

!=== unpack gathered data into cprj(natom,n2dim*nproc)
!=== second dimension is rank-ordered
 ipck=0
 do iproc = 1, nproc
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       cprj_gat(iat,(iproc-1)*n2dim+jj)%cp(:,1:nn)=buffer_cpgr_all(:,1,ipck+1:ipck+nn)
       if (ncpgr/=0) cprj_gat(iat,(iproc-1)*n2dim+jj)%dcp(:,1:ncpgr,1:nn)=&
&       buffer_cpgr_all(:,2:1+ncpgr,ipck+1:ipck+nn)
       ipck=ipck+nn
     end do
   end do
 end do

!=== Clean up ===
 ABI_DEALLOCATE(buffer_cpgr)
 ABI_DEALLOCATE(buffer_cpgr_all)

end subroutine cprj_mpi_allgather
!!***

!!****f* ABINIT/cprj_bcast
!! NAME
!! cprj_bcast
!!
!! FUNCTION
!! Broadcast a cprj_type from master to all nodes inside a MPI communicator.
!!
!! COPYRIGHT
!! Copyright (C) 2007-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom=Number of atoms (size of the first dimension of Cprj).
!!  n2dim=Size of the second dimension of Cprj.
!!  ncpgr=Number of gradients that have to be cast. It is a bit redundant but, it can be used to
!!   broad cast only the %cp"s without caring about the gradients. Just set it to 0 but be careful!
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  master=ID of the sending node in spaceComm.
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  ierr=Error status.
!!  Cprj(natom,n2dim)<cprj_type>=The datatype to be transmitted by master and received by the others nodes.
!!
!! PARENTS
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

subroutine cprj_bcast(natom,n2dim,nlmn,ncpgr,Cprj,master,spaceComm,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_xmpi

 use m_errors, only : assert

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_bcast'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr
 integer,intent(in) :: master,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(cprj_type),intent(inout) :: Cprj(natom,n2dim)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,n1dim,nn
! integer :: tcpgr
 integer :: ntotcp,ipck,rank,nprocs
!arrays
 real(dp),allocatable :: buffer_cp(:,:),buffer_cpgr(:,:,:)

! *************************************************************************

 ierr=0
 nprocs = xcomm_size(spaceComm)
 if (nprocs==1) return

 rank = xcomm_rank(spaceComm)

!#if defined DEBUG_MODE
 nn=size(nlmn,dim=1)
 n1dim=size(Cprj,dim=1)
 call assert((nn==n1dim),'cprj_bcast: size mismatch in natom!')
!!!tcpgr=Cprj(1,1)%ncpgr
!!!call assert((tcpgr==ncpgr),'cprj_bcast: size mismatch in ncpgr!')
!#endif

 ntotcp=n2dim*SUM(nlmn(:))

 ABI_ALLOCATE(buffer_cp,(2,ntotcp))
 if (ncpgr/=0)  then
   ABI_ALLOCATE(buffer_cpgr,(2,ncpgr,ntotcp))
 end if

!=== Master packs Cprj ===
!Write a routine to pack/unpack?
 if (rank==master) then
   ipck=0
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       buffer_cp(:,ipck+1:ipck+nn)=Cprj(iat,jj)%cp(:,1:nn)
       if (ncpgr/=0) buffer_cpgr(:,:,ipck+1:ipck+nn)=Cprj(iat,jj)%dcp(:,:,1:nn)
       ipck=ipck+nn
     end do
   end do
 end if

!=== Transmit data ===
 call xcast_mpi(buffer_cp,master,spaceComm,ierr)
 if (ncpgr/=0) then
   call xcast_mpi(buffer_cpgr,master,spaceComm,ierr)
 end if

!=== UnPack the received buffer ===
 if (rank/=master) then
   ipck=0
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       Cprj(iat,jj)%cp(:,1:nn)=buffer_cp(:,ipck+1:ipck+nn)
       if (ncpgr/=0) Cprj(iat,jj)%dcp(:,:,1:nn)=buffer_cpgr(:,:,ipck+1:ipck+nn)
       ipck=ipck+nn
     end do
   end do
 end if

 ABI_DEALLOCATE(buffer_cp)
 if (ncpgr/=0)  then
   ABI_DEALLOCATE(buffer_cpgr)
 end if

end subroutine cprj_bcast
!!***

!!****f* ABINIT/cprj_transpose
!! NAME
!! cprj_transpose
!!
!! FUNCTION
!! Transpose a cprj datstructure FOR A GIVEN (K,SPIN)
!! in order to change the parallel distribution from atom to band (or the contrary).
!! At input, cprj is distributed over bands (or atoms); at output, it is distributed over atoms (or bands)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprjin(n1indim,n2indim)<cprj_type>=the input cprj datastructure
!!  cprj_bandpp=number of bands to be treated simultaneoulsy by a processor
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands
!!  nspinor=number of spinorial components
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  cprjout(n1outdim,n2outdim)<cprj_type>=the output cprj datastructure with another distribution
!!
!! NOTES
!!  On the dimensions:
!!   To transfer cprj from band distribution to atom distribution, dimensions should be:
!!    n1indim =natom       n2indim =nband/nproc*nspinor
!!    n1outdim=natom/nproc n2outdim=nband*nspinor
!!   To transfer cprj from atom distribution to band distribution, dimensions should be:
!!    n1indim =natom       n2indim =nband/nproc*nspinor
!!    n1outdim=natom/nproc n2outdim=nband*nspinor
!!
!! PARENTS
!!      cprj_utils_mpi
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

 subroutine cprj_transpose(cprjin,cprjout,cprj_bandpp,natom,nband,nspinor,spaceComm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_transpose'
 use interfaces_44_abitypes_defs
!End of the abilint section

 implicit none

!Arguments-------------------------------------
!scalars
 integer :: cprj_bandpp,natom,nband,nspinor,spaceComm
!arrays
 type(cprj_type),intent(in) ::  cprjin(:,:)
 type(cprj_type),intent(out) :: cprjout(:,:)
!Local variables-------------------------------
!scalars
 integer :: bpp,buf_indx
 integer :: iashft,iatom,iatom_max_sd,iatom_max_rc,iatom_1,iatom_2,iatm1_sd,iatm1_rc,iatm2_sd,iatm2_rc
 integer :: ib,iband,iband_1,iband_2,iband_shift,iblock_atom,iblock_band,ibshft
 integer :: ierr,ip,ispinor,me,nba,nbb,nbnp_sd,nbnp_rc,ncpgr,nlmn,np
 integer :: rbufsize,sbufsize,size11,size12,size21,size22,transpose_mode
!arrays
 integer,allocatable :: cprjsz_atom(:),cprjsz_block(:,:)
 integer,allocatable,target :: count_atom(:),count_band(:),displ_atom(:),displ_band(:)
 integer,pointer :: scount(:),sdispl(:),rcount(:),rdispl(:)
 real(dp),allocatable :: rbuf(:),sbuf(:)

! *************************************************************************

!MPI data
 me = xcomm_rank(spaceComm)
 np = xcomm_size(spaceComm)

!Nothing to do if nprocs=1
 if (np==1) then
   call cprj_copy(cprjin,cprjout)
   return
 end if

!Compute bloc sizes
 bpp=cprj_bandpp
 nba=natom/np;if (mod(natom,np)/=0) nba=nba+1
 nbb=nband/(np*bpp)

!Check sizes, select direction of transposition
 transpose_mode=0
 size11=size(cprjin,1);size12=size(cprjin,2)
 size21=size(cprjout,1);size22=size(cprjout,2)
 if (size11==natom.and.size12==nbb*bpp*nspinor.and.&
& size21==nba.and.size22==nband*nspinor) then
   transpose_mode=1
 else if (size11==nba.and.size12==nband*nspinor.and.&
&   size21==natom.and.size22==nbb*bpp*nspinor) then
 else
   MSG_BUG('  Wrong cprjin/cprjout sizes !')
 end if

!Compute size of atom bloc (wr to cprj)
 ABI_ALLOCATE(cprjsz_atom,(natom))
 ABI_ALLOCATE(cprjsz_block,(np,nba))
 cprjsz_atom=0;cprjsz_block=0
 if (transpose_mode==1) then
   do iatom=1,natom
     cprjsz_atom(iatom)=2*cprjin(iatom,1)%nlmn*(1+cprjin(iatom,1)%ncpgr)
   end do
 else
   do iblock_atom=1,nba
     iatom=(iblock_atom-1)*np+1+me
     if (iatom<=natom) cprjsz_atom(iatom)=2*cprjin(iblock_atom,1)%nlmn*(1+cprjin(iblock_atom,1)%ncpgr)
   end do
   call xsum_mpi(cprjsz_atom,spaceComm,ierr)
 end if
 do iblock_atom=1,nba
   iashft=(iblock_atom-1)*np
   iatom_1=iashft+1;iatom_2=iashft+np
   if (iatom_1>natom) cycle
   if (iatom_2>natom) iatom_2=natom
   do iatom=iatom_1,iatom_2
     cprjsz_block(iatom-iashft,iblock_atom)=cprjsz_atom(iatom)+2  ! +2 for nlmn et ncpgr
   end do
 end do
 ABI_DEALLOCATE(cprjsz_atom)

!Allocations for MPI_ALLTOALL
 ABI_ALLOCATE(count_atom,(np))
 ABI_ALLOCATE(displ_atom,(np))
 ABI_ALLOCATE(count_band,(np))
 ABI_ALLOCATE(displ_band,(np))

!Loop on blocks of bands
 do iblock_band=1,nbb !(note: np divides nband)
   ibshft=(iblock_band-1)*np*bpp
   iband_1=ibshft+1;iband_2=ibshft+np*bpp
   if (iband_1>nband.or.iband_2>nband) cycle ! for security

!  Loop on blocks of atoms
   do iblock_atom=1,nba
     iashft=(iblock_atom-1)*np
     iatom_1=iashft+1;iatom_2=iashft+np
     if (iatom_1>natom) cycle
     if (iatom_2>natom) iatom_2=natom

!    Computation of displacements and sizes of blocks when data are band-distributed
     count_band(1)=cprjsz_block(1,iblock_atom)*nspinor*bpp;displ_band(1)=0
     do ip=2,np
       count_band(ip)=cprjsz_block(ip,iblock_atom)*nspinor*bpp
       displ_band(ip)=displ_band(ip-1)+count_band(ip-1)
     end do

!    Computation of displacements and sizes of blocks when data are atom-distributed
     count_atom(1)=cprjsz_block(1+me,iblock_atom)*bpp*nspinor;displ_atom(1)=0
     do ip=2,np
       count_atom(ip)=count_atom(1)
       displ_atom(ip)=displ_atom(ip-1)+count_atom(ip-1)
     end do

!    According to transposition mode, select
!    - displacements and sizes of blocks
!    - shifts in arrays
     if (transpose_mode==1) then
       scount => count_band ; sdispl => displ_band
       rcount => count_atom ; rdispl => displ_atom
       nbnp_sd=bpp;nbnp_rc=np*bpp
       iatm1_sd=iatom_1;iatm2_sd=iatom_2
       iatm1_rc=iblock_atom;iatm2_rc=iatm1_rc
       iatom_max_sd=iatom_2;iatom_max_rc=iashft+1+me
     else
       scount => count_atom ; sdispl => displ_atom
       rcount => count_band ; rdispl => displ_band
       nbnp_sd=np*bpp;nbnp_rc=bpp
       iatm1_sd=iblock_atom;iatm2_sd=iatm1_sd
       iatm1_rc=iatom_1;iatm2_rc=iatom_2
       iatom_max_sd=iashft+1+me;iatom_max_rc=iatom_2
     end if

!    Allocation of buffers
     sbufsize=sdispl(np)+scount(np)
     rbufsize=rdispl(np)+rcount(np)
     ABI_ALLOCATE(sbuf,(sbufsize))
     ABI_ALLOCATE(rbuf,(rbufsize))

!    Coying of input cprj to buffer for sending
     buf_indx=0
     iband_shift=(iblock_band-1)*nbnp_sd-1
     if (iatom_max_sd<=natom) then
       do iatom=iatm1_sd,iatm2_sd
         do ib=1,nbnp_sd
           iband=(iband_shift+ib)*nspinor
           do ispinor=1,nspinor
             iband=iband+1
             nlmn=cprjin(iatom,iband)%nlmn;ncpgr=cprjin(iatom,iband)%ncpgr
             sbuf(buf_indx+1)=dble(nlmn) ;buf_indx=buf_indx+1
             sbuf(buf_indx+1)=dble(ncpgr);buf_indx=buf_indx+1
             sbuf(buf_indx+1:buf_indx+2*nlmn)=reshape(cprjin(iatom,iband)%cp(1:2,1:nlmn),(/2*nlmn/))
             buf_indx=buf_indx+2*nlmn
             if (ncpgr>0) then
               sbuf(buf_indx+1:buf_indx+2*ncpgr*nlmn)=reshape(cprjin(iatom,iband)%dcp(1:2,1:ncpgr,1:nlmn),(/2*ncpgr*nlmn/))
               buf_indx=buf_indx+2*ncpgr*nlmn
             end if
           end do
         end do
       end do
     end if
     if (buf_indx/=sbufsize) then
       MSG_BUG('  wrong buffer size for sending !')
     end if

!    Main call to MPI_ALLTOALL
     call xalltoallv_mpi(sbuf,scount,sdispl,rbuf,rcount,rdispl,spaceComm,ierr)

!    Retrieving of output cprj for received buffer
     buf_indx=0
     iband_shift=(iblock_band-1)*nbnp_rc-1
     if (iatom_max_rc<=natom) then
       do iatom=iatm1_rc,iatm2_rc
         do ib=1,nbnp_rc
           iband=(iband_shift+ib)*nspinor
           do ispinor=1,nspinor
             iband=iband+1
             nlmn =int(rbuf(buf_indx+1));buf_indx=buf_indx+1
             ncpgr=int(rbuf(buf_indx+1));buf_indx=buf_indx+1
             cprjout(iatom,iband)%nlmn=nlmn;cprjout(iatom,iband)%ncpgr=ncpgr
             cprjout(iatom,iband)%cp(1:2,1:nlmn)=reshape(rbuf(buf_indx+1:buf_indx+2*nlmn),(/2,nlmn/))
             buf_indx=buf_indx+2*nlmn
             if (ncpgr>0) then
               cprjout(iatom,iband)%dcp(1:2,1:ncpgr,1:nlmn)=reshape(rbuf(buf_indx+1:buf_indx+2*nlmn*ncpgr),(/2,ncpgr,nlmn/))
               buf_indx=buf_indx+2*nlmn*ncpgr
             end if
           end do
         end do
       end do
     else
       cprjout(iatom,iband)%nlmn=0;cprjout(iatom,iband)%ncpgr=0
       nullify(cprjout(iatom,iband)%cp);nullify(cprjout(iatom,iband)%dcp)
     end if
     if (buf_indx/=rbufsize) then
       MSG_BUG('  Error: wrong buffer size for receiving !')
     end if

!    Deallocation of buffers
     ABI_DEALLOCATE(sbuf)
     ABI_DEALLOCATE(rbuf)

!    End of loops
   end do ! do iblock_atom
 end do ! do iblock_atom

!Free memory
 ABI_DEALLOCATE(count_atom)
 ABI_DEALLOCATE(displ_atom)
 ABI_DEALLOCATE(count_band)
 ABI_DEALLOCATE(displ_band)
 ABI_DEALLOCATE(cprjsz_block)
 nullify(scount,rcount,sdispl,rdispl)

 end subroutine cprj_transpose
!!***


!!****f* ABINIT/cprj_transpose_all
!! NAME
!! cprj_transpose_all
!!
!! FUNCTION
!! Transpose a cprj datstructure FOR ALL (K,SPIN)
!! in order to change the parallel distribution from atom to band (or the contrary).
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprjin(:,:)=the input cprj datastructure
!!  dtfil <type(datafiles_type)>=variables related to files
!!  mband=maximum number of bands
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  nband(nkpt*nsppol)=number of bands at this k point for that spin polarization
!!  nkpt=number of k points
!!  natom=number of atoms in cell.
!!  nspinor=number of spinorial components
!!  nsppol=number of spin components
!!  paral_kgb=1 if (k points,bands,g vectors) parallelism is activated
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  cprjout(:,:)=the output cprj datastructure (bands and atoms transposed)
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE

 subroutine cprj_transpose_all(cprjin,cprjout,dtfil,mband,mkmem,mpi_enreg,natom,nband,&
&                              nkpt,nspinor,nsppol,paral_kgb,spaceComm)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_transpose_all'
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi, except_this_one => cprj_transpose_all
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,natom,nkpt,nspinor,nsppol,paral_kgb,spaceComm
 type(datafiles_type),intent(in) :: dtfil
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer :: nband(nkpt*nsppol)
 type(cprj_type),intent(in) :: cprjin(:,:)
 type(cprj_type),intent(out) :: cprjout(:,:)
!Local variables-------------------------------
!scalars
 integer :: bandpp,iatom,iba,ibg_in,ibg_out,ierr,ikpt,iorder_cprj,isppol,mcprj,me,me_kpt
 integer :: nba,nband_k,nband_in,nband_out,ncprj_in,ncprj_out,np
 integer :: sizei1,sizei2,sizeo1,sizeo2,transpose_mode
 character(len=500) :: msg
!arrays
 integer,allocatable,target :: dimcprj(:),dimcprj_atm(:)
 integer,pointer :: atind_in(:),atind_out(:),dimcprj_in(:),dimcprj_out(:)
 type(cprj_type),allocatable :: cprjin_k(:,:),cprjout_k(:,:)

! *************************************************************************

 DBG_ENTER('COLL')

!Init parallelism
 me = xcomm_rank(spaceComm)
 np = xcomm_size(spaceComm)
 call xme_init(mpi_enreg,me_kpt)
 bandpp=1;if (paral_kgb/=0) bandpp=mpi_enreg%bandpp

!Compatibility tests
 if (mkmem==0.and.np>1) then
   msg=' Not available for mkmem=0 and parallelism over bands!'
   MSG_ERROR(msg)
 end if

!Nothing to do if mkmem==0 (cprj is already written in file)
 if (mkmem==0.and.np==1) return

!Check/retrieve sizes of cprj ; allocate cprjout
 nba=natom/np;if (mod(natom,np)/=0) nba=nba+1
 mcprj=nspinor*mband*mkmem*nsppol
 transpose_mode=0
 ABI_ALLOCATE(dimcprj,(natom))
 ABI_ALLOCATE(dimcprj_atm,(nba))
 dimcprj_atm=0
 sizei1=size(cprjin ,1);sizei2=size(cprjin ,2)
 sizeo1=size(cprjout,1);sizeo2=size(cprjout,2)
 if (sizei1==natom.and.sizei2==mcprj/np.and.sizeo1==nba.and.sizeo2==mcprj) then
   transpose_mode=1
   ncprj_in=natom;ncprj_out=nba
   do iatom=1,natom
     dimcprj(iatom)=cprjin(iatom,1)%nlmn
   end do
   do iba=1,nba
     iatom=(iba-1)*np+1+me
     if (iatom<=natom) dimcprj_atm(iba)=dimcprj(iatom)
   end do
   dimcprj_in => dimcprj
   dimcprj_out => dimcprj_atm
 else if (sizei1==nba.and.sizei2==mcprj.and.sizeo1==natom.and.sizeo2==mcprj/np) then
   transpose_mode=2
   ncprj_in=nba;ncprj_out=natom
   do iba=1,nba
     iatom=(iba-1)*np+1+me
     if (iatom<=natom) then
       dimcprj_atm(iba)=cprjin(iba,1)%nlmn
       dimcprj(iatom)=dimcprj_atm(iba)
     end if
   end do
   call xsum_mpi(dimcprj,spaceComm,ierr)
   dimcprj_in => dimcprj_atm
   dimcprj_out => dimcprj
 else
   MSG_BUG('  Wrong cprjin/cprjout sizes !')
 end if

!Eventually prepare cprj datastructures
 call cprj_free(cprjout)
 call cprj_alloc(cprjout,cprjin(1,1)%ncpgr,dimcprj_out)
 iorder_cprj=0
 atind_in => dimcprj_in ; atind_out => dimcprj_out  ! This are dummy pointers !!!
 call cprj_diskinit_r(atind_in,ncprj_in,iorder_cprj,mkmem,natom,cprjin(1,1)%ncpgr,dimcprj_in,nspinor,dtfil%unpaw)

!Loop over spins
 ibg_in=0;ibg_out=0
 do isppol=1,nsppol

!  Loop over k-points
   do ikpt=1,nkpt

!    Select k point to be treated by this proc
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     if(mpi_enreg%paral_compil_kpt==1) then
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me_kpt))/=0) cycle
     end if

     if (transpose_mode==1) then
       nband_in=nband_k/(np*bandpp);nband_out=nband_k
     else
       nband_in=nband_k;nband_out=nband_k/(np*bandpp)
     end if

!    Allocate cprjin_k and cprjout_k for this (K,SPIN)
!    Extract cprj for this (K,SPIN) according to mkmem
     ABI_ALLOCATE(cprjin_k ,(ncprj_in ,nspinor*nband_in ))
     ABI_ALLOCATE(cprjout_k,(ncprj_out,nspinor*nband_out))
     call cprj_alloc(cprjin_k ,cprjin_k(1,1)%ncpgr,dimcprj_in)
     call cprj_alloc(cprjout_k,cprjin_k(1,1)%ncpgr,dimcprj_out)
     call cprj_get(atind_in,cprjin_k,cprjin,ncprj_in,1,ibg_in,ikpt,iorder_cprj,isppol,&
&     mband,mkmem,mpi_enreg,natom,nband_in,nband_in,nspinor,nsppol,dtfil%unpaw)

!    Transpose cprjin_k into cprjout_k for this (K,SPIN)
     call cprj_transpose(cprjin_k,cprjout_k,bandpp,natom,nband_k,nspinor,spaceComm)

!    Export new cprjout_k to big array cprjout
     call cprj_put(atind_out,cprjout_k,cprjout,ncprj_out,1,ibg_out,ikpt,iorder_cprj,isppol,&
&     mband,mkmem,mpi_enreg,natom,nband_out,nband_out,dimcprj_out,nspinor,nsppol,spaceComm,dtfil%unpaw)

!    Shift array memory
     if (mkmem/=0) then
       ibg_in =ibg_in +nspinor*nband_in
       ibg_out=ibg_out+nspinor*nband_out
     end if

!    Deallocate temporary storage
     call cprj_free(cprjin_k)
     call cprj_free(cprjout_k)
     ABI_DEALLOCATE(cprjin_k)
     ABI_DEALLOCATE(cprjout_k)

!    End loop over k
   end do
!  End loop over spins
 end do

!Deallocate temporary storage
 ABI_DEALLOCATE(dimcprj)
 ABI_DEALLOCATE(dimcprj_atm)
 nullify(atind_in,atind_out,dimcprj_in,dimcprj_out)

 DBG_EXIT('COLL')

 end subroutine cprj_transpose_all
!!***


!!****f* ABINIT/cprj_gather_spin
!! NAME
!! cprj_gather_spin
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MD,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprj(:,:)=the input cprj datastructure
!!  n2size=number of cprj datastructures to be gathered (second dim)
!!  nspinor : number of spinorial component (on current proc)
!!  nspinortot : total number of spinorial component
!!
!! OUTPUT
!!  cprj_gat(:,:) = the cprj containing all nspinor componants
!!
!! NOTES
!! The cprj has been built like the following:
!!   loop on nsppol
!!   loop on k point
!!   loop over band or block of band
!! These quantities were build only if treated by the current proc
!! the inner quantities being nspinor
!!
!! PARENTS
!!      energy,pawmkrhoij
!!
!! CHILDREN
!!      xallgather_mpi
!!
!! SOURCE
 subroutine cprj_gather_spin(cprj,cprj_gat,natom,n2size,nspinor,nspinortot,&
&                            spaceComm_spin,ierr)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cprj_gather_spin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nspinor,nspinortot,n2size
 integer,intent(in) :: spaceComm_spin
 integer,intent(out) :: ierr
!arrays
 type(cprj_type),intent(in) :: cprj(:,:)
 type(cprj_type),intent(out) :: cprj_gat(:,:)
!Local variables-------------------------------
!scalars
 integer :: i1,iatom,ibsp,icpgr,ilmn,isp,ispinor,jj,lmndim,n2dim,n2dim_gat,ncpgr
 character(len=500) :: msg
!arrays
 integer :: nlmn(natom)
 real(dp),allocatable :: buffer1(:),buffer2(:)

! *************************************************************************

 n2dim    =size(cprj,dim=2)
 n2dim_gat=size(cprj_gat,dim=2)
 if (n2dim_gat/=(nspinortot/nspinor)*n2dim) then
   write(msg,'(a)') "  Wrong dims for cprj and cprj_gat !"
   MSG_BUG(msg)
 end if

 do iatom=1,natom
   nlmn(iatom)=size(cprj(iatom,1)%cp(1,:))
 end do
 ncpgr=cprj(1,1)%ncpgr
 lmndim=2*n2size*sum(nlmn(1:natom))*(1+ncpgr)
 ABI_ALLOCATE(buffer1,(lmndim))
 ABI_ALLOCATE(buffer2,(lmndim*nspinortot))

 isp=0;ibsp=0
 jj=1
 do i1=1,n2size
   isp=isp+1
   do iatom=1,natom
     do ilmn=1,nlmn(iatom)
       buffer1(jj:jj+1)=cprj(iatom,isp)%cp(1:2,ilmn)
       jj=jj+2
     end do
     if (ncpgr>0) then
       do ilmn=1,nlmn(iatom)
         do icpgr=1,ncpgr
           buffer1(jj:jj+1)=cprj(iatom,isp)%dcp(1:2,icpgr,ilmn)
           jj=jj+2
         end do
       end do
     end if
   end do
 end do

 call xallgather_mpi(buffer1,lmndim,buffer2,spaceComm_spin,ierr)

 jj=1
 do ispinor=1,nspinortot
   do i1 =1,n2size
     ibsp=(i1-1)*nspinortot + ispinor
     do iatom=1,natom
       do ilmn=1,nlmn(iatom)
         cprj_gat(iatom,ibsp)%cp(1:2,ilmn)=buffer2(jj:jj+1)
         jj=jj+2
       end do
       if (ncpgr>0) then
         do ilmn=1,nlmn(iatom)
           do icpgr=1,ncpgr
             cprj_gat(iatom,ibsp)%dcp(1:2,icpgr,ilmn)=buffer2(jj:jj+1)
             jj=jj+2
           end do
         end do
       end if
     end do
   end do
 end do

 ABI_DEALLOCATE(buffer1)
 ABI_DEALLOCATE(buffer2)

 end subroutine cprj_gather_spin
!!***
