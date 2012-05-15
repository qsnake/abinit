!{\src2tex{textfont=tt}}
!!****f* ABINIT/datafordmft
!! NAME
!! datafordmft
!!
!! FUNCTION
!!  Compute psichi (and print some data for check)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2012 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  cryst_struc <type(crystal_structure)>=crystal structure data
!!        -gprimd(3,3)=dimensional reciprocal space primitive translations
!!        -indsym(4,nsym,natom)=indirect indexing array for atom labels
!!        -symrec(3,3,nsym)=symmetry operations in reciprocal space
!!        - nsym= number of symetry operations
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  dimcprj(natom) = dimension for cprj
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fermie= Fermi energy
!!  lda_occup <type(oper_type)> = occupations in the correlated orbitals in LDA
!!  mband=maximum number of bands
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  nkpt=number of k points.
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(mband*nkpt*nsppol) = occupancies of KS states.
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  unpaw = file number for cprj
!!
!! OUTPUT
!!  paw_dmft%psichi(nsppol,nkpt,mband,nspinor,dtset%natom,(2*maxlpawu+1))): projections <Psi|chi>
!!  paw_dmft%eigen(paw_dmft%nsppol,paw_dmft%nkpt,paw_dmft%mband)
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine datafordmft(cryst_struc,cprj,dimcprj,dtset,eigen,fermie,&
& lda_occup,mband,mkmem,mpi_enreg,nkpt,nspinor,nsppol,occ,&
& paw_dmft,paw_ij,pawang,pawtab,psps,unpaw)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_errors

 use m_paw_dmft, only: paw_dmft_type
 use m_matlu, only: matlu_type,init_matlu,sym_matlu,copy_matlu,print_matlu,diff_matlu,destroy_matlu
 use m_crystal, only : crystal_structure
 use m_oper, only : oper_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'datafordmft'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_51_manage_mpi
 use interfaces_68_dmft, except_this_one => datafordmft
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem
 integer,intent(in) :: nkpt,nspinor,nsppol
 integer,intent(in) :: unpaw
 real(dp),intent(in) :: fermie
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(oper_type), intent(out) :: lda_occup
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(crystal_structure),intent(in) :: cryst_struc
!arrays
 integer, intent(in) :: dimcprj(cryst_struc%natom)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
 type(paw_ij_type),intent(in) :: paw_ij(cryst_struc%natom*psps%usepaw)
 type(cprj_type) :: cprj(cryst_struc%natom,nspinor*mband*mkmem*nsppol)
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
!Local variables-------------------------------
!scalars
 integer :: band_index,comm_world,dimpsichi,facpara
 integer :: iat,iatom,iband,ibandc,ibg,ibsp,icat,icount_proj_ilmn,idijeff,ierr,ierrr,ikpt
 integer :: ilmn,im,im1,iorder_cprj,ispinor,ispinor1,isppol,itypat,ilmn1
 integer :: jj1,ldim,lmn_size
 integer :: m1,maxnproju,me_world,natom,nband_k
 integer :: nbandi,nbandf,nnn,nsploop,option,spaceComm
 real(dp) :: ph0phiint_used
 character(len=500) :: message
!arrays
 real(dp) :: chinorm
 complex(dpc), allocatable :: buffer1(:)
 type(matlu_type), allocatable :: loc_occ_check(:)
 type(matlu_type), allocatable :: loc_norm_check(:)
 type(matlu_type), allocatable :: xocc_check(:)
 type(matlu_type), allocatable :: xnorm_check(:)
 type(matlu_type), allocatable :: matlu_temp(:)
 logical :: lprojchi,t2g
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
!************************************************************************

!DBG_ENTER("COLL")

 facpara=1 !mpi_enreg%nproc
 if(abs(dtset%pawprtvol)>=3) then
   write(message,*) " number of k-points used is nkpt=nkpt ", nkpt
   call wrtout(std_out,  message,'COLL')
   write(message,*) " warning: parallelised version        ", nkpt
   call wrtout(std_out,  message,'COLL')
   write(message,*) " weights k-points used is wtk=wtk"
   call wrtout(std_out,  message,'COLL')
 end if

!----------------------------------- MPI-------------------------------------

!Init parallelism
 call xcomm_world(mpi_enreg,comm_world,myrank=me_world)
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_kpt)

!Prepare temporary PAW file if mkmem==0
 iorder_cprj=0
 call cprj_diskinit_r(cryst_struc%atindx1,dtset%natom,iorder_cprj,dtset%mkmem,dtset%natom,0,dimcprj,nspinor,unpaw)
!todo_ab: extract cprj from file unpaw in the following..

!----------------------------------- MPI-------------------------------------

 nbandi=paw_dmft%dmftbandi
 nbandf=paw_dmft%dmftbandf
 lprojchi=.false.
 lprojchi=.true.
 t2g=.false.
 natom=cryst_struc%natom

!if(mpi_enreg%me==0) write(7886,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
!if(mpi_enreg%me==1) write(7887,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
!if(mpi_enreg%me==2) write(7888,*) "in datafordmft", mpi_enreg%me, mpi_enreg%nproc
 write(message,'(2a)') ch10,&
& '  == Prepare data for DMFT calculation  '
 call wrtout(std_out,message,'COLL')
 if(abs(dtset%pawprtvol)>=3) then
   write(message, '(a,a)' ) ch10,&
&   '---------------------------------------------------------------'
   call wrtout(ab_out,message,'COLL');call wrtout(std_out,  message,'COLL')
   write(message, '(a,a,a,a,a,a,a,a,a,a,a,a)' ) ch10,&
&   '  Print useful data (as a check)',ch10,&
&   '  - Overlap of KS wfc with atomic orbital inside sphere',ch10,&
&   '  - Eigenvalues',ch10,&
&   '  - Weights of k-points',ch10,&
&   '  - Number of spins ',ch10,&
&   '  - Number of states'
   call wrtout(ab_out,message,'COLL');call wrtout(std_out,  message,'COLL')
   write(message, '(a,a)' ) ch10,&
&   '---------------------------------------------------------------'
 end if
 if(dtset%nstep==0) then
   write(message,'(a,a,a,a,a,a)')  ch10,&
&   'nstep should be greater than 1',ch10
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!********************* Max Values for U terms.
!maxlpawu=0
 maxnproju=0
 do iatom=1,natom
   if(pawtab(dtset%typat(iatom))%lpawu.ne.-1 .and. pawtab(dtset%typat(iatom))%nproju.gt.maxnproju)&
&   maxnproju=pawtab(dtset%typat(iatom))%nproju
 end do
!*****************   in forlb.eig
 if(me_world.eq.0.and.abs(dtset%pawprtvol)>=3) then
   open(unit=2010,file='forlb.eig',form='formatted',status='unknown')
   rewind(2010)
   write(2010,*) "Number of bands,   spins, and  k-point; and spin-orbit flag"
   write(2010,*) mband,nsppol,nkpt,nspinor,nbandi,nbandf
   write(2010,*) " For each k-point, eigenvalues for each band"
   write(2010,*) (dtset%wtk(ikpt)*facpara,ikpt=1,nkpt)
   band_index=zero
   do isppol=1,nsppol
     write(2010,*) " For spin"
     write(2010,*)  isppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       ibandc=0
       write(2010,*) " For k-point"
       write(2010,*)  ikpt
       do iband=1,mband
         if(paw_dmft%band_in(iband)) then
           ibandc=ibandc+1
           write(2010, '(2i6,4x,f20.15)' ) ibandc,ikpt,eigen(iband+band_index)*2.d0
         end if
       end do
       band_index=band_index+nband_k
     end do
   end do
   close(2010)
 end if ! proc=me

!==   put eigen into eigen_lda
 band_index=zero
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     ibandc=0
     do iband=1,mband
       if(paw_dmft%band_in(iband)) then
         ibandc=ibandc+1
         paw_dmft%eigen_lda(isppol,ikpt,ibandc)=eigen(iband+band_index) ! in Ha
       end if
     end do
     band_index=band_index+nband_k
   end do
 end do

 if(abs(dtset%pawprtvol)>=3) then
   write(message, '(a,a)' ) ch10,&
&   '   datafordmft :  eigenvalues written'
   call wrtout(std_out,  message,'COLL')
 end if

!==========================================================================
!***************** Compute  <Psi|Chi>=\sum_{proja} <Psi|P_a><phi_a|Chi>
!==========================================================================

 ibg=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
     if(mpi_enreg%paral_compil_kpt==1)then
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me_world))/=0) then
!        shift of ibg must not be done here, because, cprj is defined for nkpt/nproc
!        Skip the rest of the k-point loop
         cycle
       end if
     end if
!    write(2011,*) ikpt
     ibsp=ibg
     ibandc=0
     do iband=1,nband_k
       if(paw_dmft%band_in(iband)) ibandc=ibandc+1
       do ispinor=1,nspinor
         ibsp=ibsp+1
         icat=1
!        write(std_out,*) isppol,ikpt,iband,ispinor
         do itypat=1,dtset%ntypat
           lmn_size=pawtab(itypat)%lmn_size
           iat=0 ! to declare
!          write(std_out,*) isppol,ikpt,iband,ispinor
           do iatom=icat,icat+cryst_struc%nattyp(itypat)-1
             iat=iat+1
             jj1=0
             if(pawtab(dtset%typat(iatom))%lpawu.ne.-1) then
!              chinorm=(pawtab(itypat)%phiphjint(1))
               chinorm=1.d0
!              write(std_out,*) isppol,ikpt,iband,ispinor,iat
               do ilmn=1,lmn_size
!                write(std_out,*) ilmn
!                ------------ Select l=lpawu.
                 if (psps%indlmn(1,ilmn,itypat)==pawtab(dtset%typat(iatom))%lpawu) then
!                  ------------ Check that the band is choosen (within nbandi and nbandf)
                   if(paw_dmft%band_in(iband)) then
!                    write(std_out,*) "inside paw_dmft%band_in",iband
                     jj1=jj1+1
                     if(jj1>pawtab(dtset%typat(iatom))%nproju*(2*pawtab(dtset%typat(iatom))%lpawu+1)) then
                       write(message,'(a,a,a,a)')  ch10,&
&                       'BUG: jj1 is not correct in datafordmft',ch10,&
&                       'Action: CONTACT Abinit group'
                       call wrtout(std_out,  message,'COLL')
                       call leave_new('COLL')
                     end if ! jj1
                     icount_proj_ilmn=psps%indlmn(3,ilmn,itypat)  ! n: nb of projector
!                    write(std_out,*) "icount_proj_ilmn",icount_proj_ilmn
                     m1=psps%indlmn(2,ilmn,itypat)+psps%indlmn(1,ilmn,itypat)+1
!                    ---- if lprochi=true do the sum over every projectors
!                    ---- if lprochi=false: . do the sum over only ground state projector
!                    ----                   . and take ph0phiint(icount_proj_ilmn)=1
                     if(lprojchi.or.icount_proj_ilmn==1) then
                       if(.not.lprojchi) ph0phiint_used=one
                       if(lprojchi) ph0phiint_used=pawtab(itypat)%ph0phiint(icount_proj_ilmn)
                       paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)=&
&                       paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)+&
&                       cmplx(cprj(iatom,ibsp)%cp(1,ilmn)*ph0phiint_used,&
&                       cprj(iatom,ibsp)%cp(2,ilmn)*ph0phiint_used,kind=dp)
                     end if
                   end if ! paw_dmft%band_in
                 end if
               end do !ilmn
             end if ! lpawu.ne.-1
           end do ! iatom
           icat=icat+cryst_struc%nattyp(itypat)
         end do ! itypat
       end do ! ispinor
     end do !iband
     ibg=ibg+nband_k*nspinor
   end do !ikpt
 end do ! isppol
 if(abs(dtset%pawprtvol)>=3) then
   write(message,*) "chinorm used here =",chinorm
   call wrtout(std_out,  message,'COLL')
 end if

!==========================================================================
!********************* Gather information for MPI before printing
!==========================================================================

 dimpsichi=2*nsppol*nkpt*mband*nspinor*natom*(2*paw_dmft%maxlpawu+1)
 ABI_ALLOCATE(buffer1,(dimpsichi))
 nnn=0
!write(176,*) "beg",psichi
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do ibandc=1,paw_dmft%mbandc
       do ispinor=1,nspinor
         do iat=1,natom
           do m1=1,2*paw_dmft%maxlpawu+1
!            do m=1,2
             nnn=nnn+1
             buffer1(nnn)=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)
!            enddo
           end do
         end do
       end do
     end do
   end do
 end do
 call xbarrier_mpi(spaceComm)
 call xsum_mpi(buffer1,spaceComm ,ierr)
 call xbarrier_mpi(spaceComm)
 nnn=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do ibandc=1,paw_dmft%mbandc
       do ispinor=1,nspinor
         do iat=1,natom
           do m1=1,2*paw_dmft%maxlpawu+1
!            do m=1,2
             nnn=nnn+1
             paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1)=buffer1(nnn)
!            enddo
           end do
         end do
       end do
     end do
   end do
 end do
 ABI_DEALLOCATE(buffer1)
!write(177,*) "end",psichi

 call xbarrier_mpi(spaceComm)
!if(mpi_enreg%me.eq.0) write(177,*) "end",psichi
!if(mpi_enreg%me.eq.1) write(178,*) "end",psichi
!if(mpi_enreg%me.eq.2) write(179,*) "end",psichi
 call xbarrier_mpi(spaceComm)
!==========================================================================
!********* WRITE psichi in file for reference
!==========================================================================
 if(me_world.eq.0) then
   call psichi_print(dtset,maxnproju,cryst_struc%nattyp,cryst_struc%ntypat,nkpt,nspinor,&
&   nsppol,paw_dmft,pawtab,psps,t2g)
 end if ! proc=me
!==========================================================================
!********************* Check normalization  and occupations ***************
!==========================================================================
 ABI_ALLOCATE(xocc_check,(natom))
 ABI_ALLOCATE(xnorm_check,(natom))
 call init_matlu(natom,nspinor,nsppol,pawtab(cryst_struc%typat(1:natom))%lpawu,xocc_check)
 call init_matlu(natom,nspinor,nsppol,pawtab(cryst_struc%typat(1:natom))%lpawu,xnorm_check)
 call psichi_check(dtset,cryst_struc%nattyp,nkpt,nspinor,&
& nsppol,cryst_struc%ntypat,paw_dmft,pawtab,psps,xocc_check,xnorm_check)
!==========================================================================
!***************  write checks  *******************************************
!==========================================================================
 if(abs(dtset%pawprtvol)>=3) then
   write(message,*) "normalization computed"
   call wrtout(std_out,  message,'COLL')
 end if

 ABI_ALLOCATE(loc_occ_check,(natom))
 ABI_ALLOCATE(loc_norm_check,(natom))
 call init_matlu(natom,nspinor,nsppol,pawtab(cryst_struc%typat(1:natom))%lpawu,loc_occ_check)
 call init_matlu(natom,nspinor,nsppol,pawtab(cryst_struc%typat(1:natom))%lpawu,loc_norm_check)
 call copy_matlu(xocc_check,loc_occ_check,natom)
 call copy_matlu(xnorm_check,loc_norm_check,natom)

 write(message,'(2a,i4)')  ch10," == Check: Occupations and Norm from psichi are"
 call wrtout(std_out,  message,'COLL')

 if(paw_dmft%dmftcheck>=1) then
!  print occupations
   write(message,'(2a,i4)')  ch10,'  ------ Unsymetrised Occupation'
   call wrtout(std_out,  message,'COLL')

   call print_matlu(xocc_check,natom,dtset%pawprtvol)

!  print norms
   write(message,'(2a,i4)')  ch10,'  ------ Unsymetrised Norm'
   call wrtout(std_out,  message,'COLL')

   call print_matlu(xnorm_check,natom,dtset%pawprtvol)
 end if

!symetrise and print occupations
 call sym_matlu(cryst_struc,loc_occ_check,pawang)

 write(message,'(2a,i4)')  ch10,'  ------ Symetrised Occupation'
 call wrtout(std_out,  message,'COLL')

 call print_matlu(loc_occ_check,natom,dtset%pawprtvol)

!symetrise and print norms
 call sym_matlu(cryst_struc,loc_norm_check,pawang)

 write(message,'(2a,i4)')  ch10,'  ------ Symetrised Norm'
 call wrtout(std_out,  message,'COLL')

 call print_matlu(loc_norm_check,natom,dtset%pawprtvol)

!deallocations
 do iatom=1,natom
   lda_occup%matlu(iatom)%mat=loc_occ_check(iatom)%mat
 end do

!Tests density matrix LDA+U and density matrix computed here.
 if(paw_dmft%dmftcheck==2.or.(paw_dmft%dmftbandi==1)) then
   ABI_ALLOCATE(matlu_temp,(natom))
   call init_matlu(natom,paw_dmft%nspinor,paw_dmft%nsppol,pawtab(cryst_struc%typat(1:natom))%lpawu,matlu_temp)
   do iatom=1,natom
     if(pawtab(dtset%typat(iatom))%lpawu.ne.-1) then
       ldim=2*pawtab(dtset%typat(iatom))%lpawu+1
       nsploop=max(paw_dmft%nsppol,paw_dmft%nspinor**2)
       do idijeff=1,nsploop
         ispinor=0
         ispinor1=0
         if(nsploop==2) then
           isppol=spinor_idxs(1,idijeff)
           ispinor=1
           ispinor1=1
         else if(nsploop==4) then
           isppol=1
           ispinor=spinor_idxs(1,idijeff)
           ispinor1=spinor_idxs(2,idijeff)
         else
           write(message,'(2a)') " BUG in ldau_self: nsploop should be equal to 2 or 4"
           call wrtout(std_out,message,'COLL')
         end if
         do im1 = 1 , ldim
           do im = 1 ,  ldim
             if(nspinor==2) matlu_temp(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=&
&             cmplx(paw_ij(iatom)%noccmmp(1,im,im1,idijeff),paw_ij(iatom)%noccmmp(2,im,im1,idijeff),kind=dp)
             if(nspinor==1) matlu_temp(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=&
&             cmplx(paw_ij(iatom)%noccmmp(1,im,im1,idijeff),zero,kind=dp)
           end do
         end do
       end do
     end if
   end do
   if(paw_dmft%dmftcheck==2) option=1
   if(paw_dmft%dmftcheck<=1) option=0
   call diff_matlu("LDA+U density matrix",&
&   "Direct calculation of density matrix with psichi ",&
&   matlu_temp,lda_occup%matlu,natom,option,tol4,ierrr)
   if(ierrr==-1) then
     write(message,'(2a)') ch10,&
&     '    -> To correct this, check that input wavefunctions come from the same Hamiltonien (e.g LDA/GGA) and dmatpuopt==1'
     call wrtout(std_out,message,'COLL')
   end if
!  write(message,'(2a)') ch10,&
!  &   '  ***** => Calculations of density matrices with projections and in LDA+U are coherent****'
!  call wrtout(std_out,message,'COLL')

   call destroy_matlu(matlu_temp,natom)
   ABI_DEALLOCATE(matlu_temp)
 else
   write(message,'(2a)') ch10,&
&   '  Warning: Consistency of density matrices computed from projection has not been checked: use dmftcheck>=2 '
   call wrtout(std_out,message,'COLL')
 end if

 call destroy_matlu(loc_norm_check,natom)
 ABI_DEALLOCATE(loc_norm_check)
 call destroy_matlu(loc_occ_check,natom)
 ABI_DEALLOCATE(loc_occ_check)

 call destroy_matlu(xocc_check,natom)
 call destroy_matlu(xnorm_check,natom)
 ABI_DEALLOCATE(xocc_check)
 ABI_DEALLOCATE(xnorm_check)

!*********************
!deallocate(psichi)
!DEBUG
!write(std_out,*)' datafordmft : exit'
!stop
!ENDDEBUG
 CONTAINS
!===========================================================
!!***

!!****f* datafordmft/psichi_print
!! NAME
!!  psichi_print
!!
!! FUNCTION
!!  Print psichi for reference
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  maxnproju = maximum number of projector for LDA+U species
!!  nattyp(ntypat)= # atoms of each type
!!  mband= number of bands
!!  nkpt=number of k points
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psichi(2,nsppol,nkpt,mband,nspinor,dtset%natom,(2*paw_dmft%maxlpawu+1))) projections for DMFT
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!!
!! SIDE EFFECTS
!!  print psichi in forlb.ovlp
!!
!! PARENTS
!!      datafordmft
!!
!! CHILDREN
!!
!! SOURCE

subroutine psichi_print(dtset,maxnproju,nattyp,ntypat,nkpt,nspinor,&
&nsppol,paw_dmft,pawtab,psps,t2g)

 use m_profiling

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psichi_print'
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: maxnproju,nkpt,nspinor,nsppol,ntypat
!arrays
 integer, intent(in) :: nattyp(ntypat)
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 logical t2g
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables ------------------------------------
 integer :: ibg,isppol,ikpt,iband,ibandc,ispinor,icat,itypat,lmn_size
 integer :: iat,iatom,jj1,ilmn,icount_proj_ilmn,m1,nband_k
 real(dp) :: chinorm

! *********************************************************************

   open(unit=2011,file='forlb.ovlp',form='formatted',status='unknown')
   rewind(2011)
   ibg=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       rewind(1023)
       write(2011,'(a6,2x,i6)') "ikpt =",ikpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       ibandc=0
       do iband=1,nband_k
         if(paw_dmft%band_in(iband)) then
           ibandc=ibandc+1
           write(2011,'(a8,2x,i6)') " iband =",iband
         end if
         do ispinor=1,nspinor
           icat=1
!          write(std_out,*) isppol,ikpt,iband,ispinor
           iat=0 ! to declare
           do itypat=1,dtset%ntypat
             lmn_size=pawtab(itypat)%lmn_size
!            write(std_out,*) isppol,ikpt,iband,ispinor
             do iatom=icat,icat+nattyp(itypat)-1
               iat=iat+1
               jj1=0
               if(pawtab(dtset%typat(iatom))%lpawu.ne.-1) then
!                chinorm=(pawtab(itypat)%phiphjint(1))
                 chinorm=1.d0
!                write(std_out,*) isppol,ikpt,iband,ispinor,iat
                 do ilmn=1,lmn_size
!                  write(std_out,*) ilmn
!                  ------------ Select l=lpawu.  ---------------------------------------
                   if (psps%indlmn(1,ilmn,itypat)==pawtab(dtset%typat(iatom))%lpawu) then
!                    ------------ Check that the band is choosen (within nbandi and nbandf)
                     if(paw_dmft%band_in(iband)) then
!                      write(std_out,*) "inside paw_dmft%band_in",iband
                       jj1=jj1+1
                       if(jj1>maxnproju*(2*paw_dmft%maxlpawu+1)) then
!                        write(std_out,*) "Error ?"
!                        write(std_out,*) jj1
                         stop
                       end if ! jj1
                       icount_proj_ilmn=psps%indlmn(3,ilmn,itypat)  ! n: nb of projector
!                      write(std_out,*) "icount_proj_ilmn",icount_proj_ilmn
                       m1=psps%indlmn(2,ilmn,itypat)+psps%indlmn(1,ilmn,itypat)+1
!                      ----- Print only when the sum over projectors is done
!                      write(std_out,*) ilmn,m1
                       if(psps%indlmn(3,ilmn,itypat)==pawtab(dtset%typat(iatom))%nproju) then
                         if(t2g) then
                           if(m1==1.or.m1==2.or.m1==4) then
                             write(2011,'(3i6,3x,2f23.15)') isppol, iat, m1,&
&                             real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1))/chinorm,&
&                             imag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1))/chinorm
                           end if
                         else
                           write(2011,'(3i6,3x,2f23.15)') isppol, iat, m1,&
&                           real(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1))/chinorm,&
&                           imag(paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m1))/chinorm
                         end if
                       end if
!                      if natom=1 then jj1 maximal value should be 2*lpawu+1
                     end if ! paw_dmft%band_in
                   end if
                 end do !ilmn
               end if ! lpawu.ne.-1
             end do ! iatom
             icat=icat+nattyp(itypat)
           end do ! itypat
         end do ! ispinor
       end do !iband
       ibg=ibg+nband_k*nspinor
     end do !ikpt
   end do ! isppol
   write(2011,*) "Fermi level (in Ryd)="
   write(2011,*) fermie*two
   close(2011)
 end subroutine psichi_print
!!***

!!****f* datafordmft/psichi_check
!! NAME
!!  psichi_check
!!
!! FUNCTION
!!  Check psichi: compute occupations
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  maxnproju = maximum number of projector for LDA+U species
!!  nattyp(ntypat)= # atoms of each type
!!  mband= number of bands
!!  nkpt=number of k points
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat= number of species
!!  paw_dmft <type(paw_dmft)>=paw data for the self-consistency
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!!  OUTPUTS:
!!  xocc_check(nsppol,nspinor,nspinor,natom,2*maxlpawu+1,2*maxlpawu+1): density matrix
!!  xnorm_check(nsppol,nspinor,nspinor,natom,2*maxlpawu+1,2*maxlpawu+1): matrix of norms
!!
!! SIDE EFFECTS
!!  check psichi: compute norm and occupations
!!
!! PARENTS
!!      datafordmft
!!
!! CHILDREN
!!
!! SOURCE

subroutine psichi_check(dtset,nattyp,nkpt,nspinor,&
& nsppol,ntypat,paw_dmft,pawtab,psps,xocc_check,xnorm_check)

 use m_profiling

 use m_matlu, only: matlu_type,init_matlu,sym_matlu

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psichi_check'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nspinor,nsppol,ntypat
!arrays
 integer, intent(in) :: nattyp(ntypat)
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(matlu_type), intent(out) :: xocc_check(dtset%natom)
 type(matlu_type), intent(out) :: xnorm_check(dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
!Local variables ------------------------------------
 integer :: band_index,ibg,isppol,ikpt,iband,ibandc,ispinor,icat,itypat
 integer :: iat,iatom,ilmn,lmn_size,m,m1,nband_k
 complex(dpc) :: psichic,psichic1
 real(dp) :: chinorm

! *********************************************************************

   ibg=0
   band_index=zero
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       ibandc=0
       do iband=1,nband_k
         if(paw_dmft%band_in(iband)) ibandc=ibandc+1
         do ispinor=1,nspinor
           icat=1
           iat=0
           do itypat=1,dtset%ntypat
             lmn_size=pawtab(itypat)%lmn_size
             do iatom=icat,icat+nattyp(itypat)-1
               iat=iat+1
!              ------------ Select correlated atoms
               if(pawtab(dtset%typat(iatom))%lpawu.ne.-1) then
                 chinorm=1.d0
                 do ilmn=1,lmn_size
!                  ------------ Select l=lpawu.
                   if (psps%indlmn(1,ilmn,itypat)==pawtab(dtset%typat(iatom))%lpawu.and.&
&                   psps%indlmn(3,ilmn,itypat)==1) then
                     do ilmn1=1,lmn_size
!                      ------------ Select l=lpawu and do not sum over projectors
!                      (this is already done in paw_dmft%psichi)
                       if (psps%indlmn(1,ilmn1,itypat)==pawtab(dtset%typat(iatom))%lpawu.and.&
&                       psps%indlmn(3,ilmn1,itypat)==1) then
!                        ------------ Check that the band is choosen (within nbandi and nbandf)
                         if(paw_dmft%band_in(iband)) then
                           m=psps%indlmn(2,ilmn,itypat)+psps%indlmn(1,ilmn,itypat)+1
                           m1=psps%indlmn(2,ilmn1,itypat)+psps%indlmn(1,ilmn,itypat)+1
                           if(psps%indlmn(3,ilmn,itypat)==1) then
                             do ispinor1=1,nspinor
                               psichic=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor,iat,m)
                               psichic1=paw_dmft%psichi(isppol,ikpt,ibandc,ispinor1,iat,m1)
!                              ------------ Compute occupation matrix
                               xocc_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)=&
&                               xocc_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)&
!                              &               +conjg(psichic)*psichic1*dtset%wtk(ikpt)*facpara*occ(iband+band_index)
&                               +conjg(psichic1)*psichic*dtset%wtk(ikpt)*facpara*occ(iband+band_index)
!                              ------------ Compute norm (should be equal to noccmmp
!                              (dmatpuopt=1) if all bands are taken into account)
                               xnorm_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)=&
&                               xnorm_check(iatom)%mat(m,m1,isppol,ispinor,ispinor1)&
!                              &               +conjg(psichic)*psichic1*dtset%wtk(ikpt)*facpara
&                               +conjg(psichic1)*psichic*dtset%wtk(ikpt)*facpara
                             end do ! ispinor1
                           end if
                         end if ! paw_dmft%band_in
                       end if
                     end do !ilmn1
                   end if
                 end do !ilmn
               end if ! lpawu.ne.-1
             end do ! iatom
             icat=icat+nattyp(itypat)
           end do ! itypat
         end do ! ispinor
       end do !iband
       band_index=band_index+nband_k
       ibg=ibg+nband_k*nspinor
     end do !ikpt
   end do ! isppol

 end subroutine psichi_check
!DBG_EXIT("COLL")

end subroutine datafordmft
!!***
