!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmkrhoij
!!
!! NAME
!! pawmkrhoij
!!
!! FUNCTION
!! Calculate the PAW quantities rhoij (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cprj(natom,mcprj)= wave functions projected with non-local projectors:
!!                     cprj_nk(i)=<p_i|Cnk> where p_i is a non-local projector.
!!  dimcprj(natom)=array of dimensions of array cprj (ordered by atom-type)
!!  istwfk(nkpt)=parameter that describes the storage of wfs
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  mband_cprj=maximum number of bands used in the dimensionning of cprj array (usually mband/nproc_band)
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands for all k points
!!  nkpt=number of k points
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k
!!  paral_kgb=Flag related to the kpoint-band-fft parallelism
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  unpaw=unit number for cprj PAW data (if used)
!!  wtk(nkpt)=weight assigned to each k point
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  On input: arrays dimensions
!!  On output:
!!    pawrhoij(:)%rhoij_(lmn2_size,nspden)=
!!          Sum_{n,k} {occ(n,k)*conjugate[cprj_nk(ii)].cprj_nk(jj)} (non symetrized)
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      cprj_alloc,cprj_diskinit_r,cprj_free,cprj_gather_spin,cprj_get
!!      leave_test,pawaccrhoij,print_ij,timab,wrtout,xcomm_init,xme_init
!!      xsum_mpi
!!
!! NOTES
!!  The cprj are distributed over band processors.
!!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projectors
!!  are stored on each proc.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawmkrhoij(atindx1,cprj,dimcprj,istwfk,kptopt,mband,mband_cprj,mcprj,mkmem,mpi_enreg,&
&                      natom,nband,nkpt,nspinor,nsppol,occ,paral_kgb,paw_dmft,&
&                      pawprtvol,pawrhoij,unpaw,wtk)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_xmpi

 use m_paw_dmft, only: paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawmkrhoij'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
 use interfaces_66_paw, except_this_one => pawmkrhoij
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: kptopt,mband,mband_cprj,mcprj,mkmem,natom,nkpt,nspinor,nsppol
 integer,intent(in) :: paral_kgb,pawprtvol,unpaw
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom),istwfk(nkpt),nband(nkpt*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),wtk(nkpt)
 type(cprj_type),target,intent(in) :: cprj(natom,mcprj)
 type(paw_dmft_type),intent(in) :: paw_dmft
 type(pawrhoij_type),intent(inout) :: pawrhoij(natom)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: max_nband_cprj=100
 integer :: bdtot_index,bufdim,cplex,iatom,ib,ib1,iband,iband1,ibc1,ibg,ib_this_proc,ierr,ikpt
 integer :: iorder_cprj,isppol,jb_this_proc,jbg,jdim,me,my_nspinor,natinc,nband_k,nband_k_cprj
 integer :: nbandc1,nband_k_cprj_read,nband_k_cprj_used,nprocband,nsp2
 integer :: option,spaceComm,use_nondiag_occup_dmft
 logical :: locc_test,usetimerev
 real(dp) :: occup,wtk_k
 character(len=500) :: message
!arrays
 integer,allocatable :: dimlmn(:),idum(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 character(len=8),parameter :: dspin(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 type(cprj_type),allocatable :: cprj_tmp(:,:),cwaveprj(:,:),cwaveprjb(:,:)
 type(cprj_type),pointer :: cprj_ptr(:,:)

!************************************************************************

 DBG_ENTER("COLL")

!Init MPI data
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_kpt)
 call xme_init(mpi_enreg,me)

!Check size of cprj
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spin)
 if (mcprj/=my_nspinor*mband_cprj*mkmem*nsppol) then
   MSG_BUG('  wrong size for cprj !')
 end if

!Check if cprj is distributed over bands
 nprocband=(mband/mband_cprj)
 if (paral_kgb==1.and.nprocband/=mpi_enreg%nproc_band) then
   MSG_BUG('  mband/mband_cprj must be equal to mband !')
 end if
 if (paw_dmft%use_dmft/=0.and.nprocband/=1) then
   MSG_BUG('  parallelization over bands not compatible with DMFT !')
 end if

!Initialise and check dmft variables
 if(paw_dmft%use_dmft/=0) then
   nbandc1=(paw_dmft%mbandc-1)*paw_dmft%use_dmft+1
 else
   nbandc1=1
 end if

!Allocate temporary cwaveprj storage
 ABI_ALLOCATE(cwaveprj,(natom,nspinor))
 call cprj_alloc(cwaveprj,0,dimcprj)
 if(paw_dmft%use_dmft/=0) then
   ABI_ALLOCATE(cwaveprjb,(natom,nspinor))
   call cprj_alloc(cwaveprjb,0,dimcprj)
 end if

!Initialize temporary file (if used)
 iorder_cprj=0
 call cprj_diskinit_r(atindx1,natom,iorder_cprj,mkmem,natom,0,dimcprj,my_nspinor,unpaw)

!Initialize output quantities
 do iatom=1,natom
   pawrhoij(iatom)%rhoij_=zero
 end do

!LOOP OVER SPINS
 option=1
 usetimerev=(kptopt>0.and.kptopt<3)
 bdtot_index=0;ibg=0;jbg=0
 do isppol=1,nsppol

!  LOOP OVER k POINTS
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     nband_k_cprj=nband_k/nprocband
     wtk_k=wtk(ikpt)

     if(mpi_enreg%paral_compil_kpt==1)then
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
         bdtot_index=bdtot_index+nband_k
         cycle
       end if
     end if

     cplex=2;if (istwfk(ikpt)>1) cplex=1

!    In case of spinors parallelism, need some extra storage
     if (mpi_enreg%paral_spin==1) then
       nband_k_cprj_used=min(max_nband_cprj,nband_k_cprj)
       ABI_ALLOCATE(cprj_tmp,(natom,my_nspinor*nband_k_cprj_used))
       ABI_ALLOCATE(cprj_ptr,(natom,   nspinor*nband_k_cprj_used))
       call cprj_alloc(cprj_tmp,0,dimcprj)
       call cprj_alloc(cprj_ptr,0,dimcprj)
     else
       cprj_ptr => cprj
     end if

!    LOOP OVER BANDS
     ib_this_proc=0;jb_this_proc=0
     do ib=1,nband_k
       iband=bdtot_index+ib

!      Parallelization: treat only some bands
       if(mpi_enreg%paral_compil_kpt==1)then
         if (paral_kgb==1) then
           if (mod((ib-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band)/=mpi_enreg%me_band) cycle
         else
           if (mpi_enreg%proc_distrb(ikpt,ib,isppol)/=me) cycle
         end if
       end if
       ib_this_proc=ib_this_proc+1

!      In case of spinors parallelism, gather cprj because we need both components together
!      We do that nband_k_cprj_used by nband_k_cprj_used bands
       if (mpi_enreg%paral_spin==1) then
         jb_this_proc=jb_this_proc+1
         if (mod(jb_this_proc,nband_k_cprj_used)==1) then
           ib_this_proc=1
           nband_k_cprj_read=nband_k_cprj_used
           if (nband_k_cprj<jb_this_proc+nband_k_cprj_used-1) nband_k_cprj_read=nband_k_cprj-jb_this_proc+1
           call cprj_get(atindx1,cprj_tmp,cprj,natom,jb_this_proc,jbg,ikpt,iorder_cprj,isppol,&
&           mband_cprj,mkmem,mpi_enreg,natom,nband_k_cprj_read,nband_k_cprj,my_nspinor,nsppol,unpaw)
           call cprj_gather_spin(cprj_tmp,cprj_ptr,natom,nband_k_cprj_read,my_nspinor,nspinor,&
&           mpi_enreg%comm_spin,ierr)
         end if
       end if

!      DMFT: LOOP ON ADDITIONAL BANDS
       do ibc1=1,nbandc1
!        check if dmft and occupations
!        write(std_out,*) 'ib,ibc1          ',ib,ibc1

!        DMFT stuff: extract cprj and occupations for additional band
         if(paw_dmft%use_dmft /= 0) then
           ib1 = paw_dmft%include_bands(ibc1)
!          write(std_out,*) 'use_dmft=1 ib,ib1',ib,ib1
           iband1 = bdtot_index+ib1
!          write(std_out,*) 'ib, ib1          ',paw_dmft%band_in(ib),paw_dmft%band_in(ib1)
           if(paw_dmft%band_in(ib)) then
             if(.not.paw_dmft%band_in(ib1))  stop
             use_nondiag_occup_dmft = 1
             occup = paw_dmft%occnd(ib,ib1,ikpt,isppol)
             locc_test = abs(paw_dmft%occnd(ib,ib1,ikpt,isppol))>tol8
!            write(std_out,*) 'use_dmft=1,band_in(ib)=1, ib,ibc1',ib,ib1,locc_test
             if (locc_test .or. mkmem == 0) then
               call cprj_get(atindx1,cwaveprjb,cprj_ptr,natom,ib1,ibg,ikpt,iorder_cprj,isppol,&
&               mband_cprj,mkmem,mpi_enreg,natom,1,nband_k_cprj,nspinor,nsppol,unpaw)
             end if
           else
             use_nondiag_occup_dmft = 0
             locc_test = (abs(occ(iband))>tol8)
             occup = occ(iband)
             if(ibc1 /= 1 .and. .not.(paw_dmft%band_in(ib))) cycle
           end if
         else  ! nbandc1=1
           use_nondiag_occup_dmft=0
           locc_test = (abs(occ(iband))>tol8)
           occup = occ(iband)
         end if

!        Extract cprj for current band
!        Must read cprj when mkmem=0 (even if unused) to have right pointer inside _PAW file
         if (locc_test.or.mkmem==0) then
           call cprj_get(atindx1,cwaveprj,cprj_ptr,natom,ib_this_proc,ibg,ikpt,iorder_cprj,isppol,&
&           mband_cprj,mkmem,mpi_enreg,natom,1,nband_k_cprj,nspinor,nsppol,unpaw)
         end if

!        Accumulate contribution from (occupied) current band
         if (locc_test) then
           if(use_nondiag_occup_dmft == 1) then
             call pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprjb,0,isppol,natom,&
&             nspinor,occup,option,pawrhoij,usetimerev,wtk_k)
           else
             call pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj ,0,isppol,natom,&
&             nspinor,occup,option,pawrhoij,usetimerev,wtk_k)
           end if
         end if

       end do ! ib1c
     end do ! ib

     if (mpi_enreg%paral_spin==1) then
       call cprj_free(cprj_tmp)
       call cprj_free(cprj_ptr)
       ABI_DEALLOCATE(cprj_tmp)
       ABI_DEALLOCATE(cprj_ptr)
     else
       nullify(cprj_ptr)
     end if

     bdtot_index=bdtot_index+nband_k
     if (mkmem/=0) then
       if (mpi_enreg%paral_spin==0) then
         ibg=ibg+   nspinor*nband_k_cprj
       else
         jbg=jbg+my_nspinor*nband_k_cprj
       end if
     end if

   end do ! ikpt
 end do ! isppol

!deallocate temporary cwaveprj/cprj storage
 call cprj_free(cwaveprj)
 ABI_DEALLOCATE(cwaveprj)
 if(paw_dmft%use_dmft/=0) then
   call cprj_free(cwaveprjb)
   ABI_DEALLOCATE(cwaveprjb)
 end if

!MPI: need to exchange arrays between procs
!==========================================
 if(mpi_enreg%paral_compil_kpt==1)then
   call leave_test()

!  Exchange pawrhoij%rhoij_
   call timab(48,1,tsec)
   ABI_ALLOCATE(dimlmn,(natom))
   dimlmn(1:natom)=pawrhoij(1:natom)%cplex*pawrhoij(1:natom)%lmn2_size
   nsp2=pawrhoij(1)%nsppol;if (pawrhoij(1)%nspden==4) nsp2=4
   bufdim=sum(dimlmn)*nsp2
   ABI_ALLOCATE(buffer1,(bufdim))
   ABI_ALLOCATE(buffer2,(bufdim))
   jdim=0
   do iatom=1,natom
     do isppol=1,nsp2
       buffer1(jdim+1:jdim+dimlmn(iatom))=pawrhoij(iatom)%rhoij_(:,isppol)
       jdim=jdim+dimlmn(iatom)
     end do
   end do
   call xsum_mpi(buffer1,buffer2,bufdim,spaceComm,ierr) !Build sum of everything
   if (paral_kgb==1.and.nprocband>1) then
     call xsum_mpi(buffer2,mpi_enreg%comm_band,ierr) !Build sum over band processors
   end if
   jdim=0
   do iatom=1,natom
     do isppol=1,nsp2
       pawrhoij(iatom)%rhoij_(:,isppol)=buffer2(jdim+1:jdim+dimlmn(iatom))
       jdim=jdim+dimlmn(iatom)
     end do
   end do
   ABI_DEALLOCATE(buffer1)
   ABI_DEALLOCATE(buffer2)
   ABI_DEALLOCATE(dimlmn)
   call timab(48,2,tsec)
 end if ! mpi_enreg%paral_compil_kpt==1

!Print info
 if (abs(pawprtvol)>=1) then
   natinc=1;if(natom>1.and.pawprtvol>=0) natinc=natom-1
   do iatom=1,natom,natinc
     nsp2=pawrhoij(iatom)%nsppol;if (pawrhoij(iatom)%nspden==4) nsp2=4
     write(message, '(4a,i3,a)') ch10," PAW TEST:",ch10,&
&     ' ====== Values of RHOIJ in pawmkrhoij (iatom=',iatom,') ======'
     if (pawrhoij(iatom)%nspden==2.and.pawrhoij(iatom)%nsppol==1) write(message,'(3a)') trim(message),ch10,&
&     '      (antiferromagnetism case: only one spin component)'
     call wrtout(std_out,message,'COLL')
     do isppol=1,nsp2
       if (pawrhoij(iatom)%nspden/=1) then
         write(message, '(3a)') '   Component ',trim(dspin(isppol+2*(pawrhoij(iatom)%nspden/4))),':'
         call wrtout(std_out,message,'COLL')
       end if
       option=2;if (pawrhoij(iatom)%cplex==2.and.pawrhoij(iatom)%nspinor==1) option=1
       call print_ij(pawrhoij(iatom)%rhoij_(:,isppol),pawrhoij(iatom)%lmn2_size,&
&       pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,&
&       -1._dp,1,opt_sym=option)
     end do
   end do
 end if

 DBG_EXIT("COLL")

end subroutine pawmkrhoij
!!***
