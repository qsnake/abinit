!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_kpgio
!! NAME
!! prep_kpgio
!!
!! FUNCTION
!! Used only in BandFFT parallelization.
!! Do initialization of kg information.
!! Transposition of all quantities which depend on kg_k
!! In case of istwfk==2, this routine completes also
!! the kg_k_gather vector by the opposite values.
!! The values are distributed on the processors in function of
!! the value of modulo(-kg_k_gather(2,i),nproc_fft)
!!
!! COPYRIGHT
!!
!! INPUTS
!!  accesswff          = defines the format of the output
!!  ecut_eff           = kinetic energy planewave cutoff (hartree)
!!  exchn2n3d          = if 1, n2 and n3 are exchanged
!!  gmet(3,3)          = reciprocal space metric (bohr^-2)
!!  istwfk(nkpt)       = input option parameter that describes the storage of wfs
!!  kptns(3,nkpt)      = reduced coords of k points
!!  fnametmp_kg        = name of unkg file
!!  mgfft              = maximum single fft dimension (IN)
!!  mkmem              = number of k points which can fit in memory; set to 0 if use disk
!!  mode_paral         = either 'COLL' or 'PERS', tells whether
!!                       the loop over k points must be done by all processors or not,
!!                       in case of parallel execution.
!!  mpi_enreg          = informations about MPI parallelization
!!  mpw                = maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol) = number of bands at each k point
!!  nkpt               = number of k points
!!  nsppol             = 1 for unpolarized, 2 for polarized
!!  unkg               = unit number for storage of basis sphere data: stores indirect
!!                       indexing array and integer coordinates for all planewaves in basis
!!                       sphere for each k point being considered
!!
!! OUTPUT
!!  kg(3,mpw*mkmem)    = dimensionless coords of G vecs in basis sphere at k point
!!  npwarr(nkpt)       = array holding npw for each k point, taking into account
!!                       the effect of istwfk, and the spreading over processors
!!  npwtot(nkpt)       = array holding the total number of plane waves for each k point,
!!---------------------------------------------------------------------
!!  within the mpi_enreg%bandfft_kpt data_type : Initialize and compute
!!---------------------------------------------------------------------
!!  gbound             = sphere boundary info
!!  idatarecv0         = position of the planewave coordinates (0,0,0)
!!  istwf_k            = input option parameter that describes the storage of wfs
!!  kg_k_gather        = planewave coordinates
!!                       (of the processor + sended by other processors band)
!!  kg_k_gather_sym    = planewave coordinates
!!                       (kg_k_gather + opposited planewave coordinates sended by the processors fft)
!!  ndatarecv          = total number of values received by the processor and sended
!!                       by the other processors band
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited
!!                       planewave coordinates
!!  ndatarecv_tot      = total number of received values by the processor
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  recvcounts         = number of values received by the  processor from each processor band
!!  recvcounts_sym     = number of values received by the  processor from each processor fft
!!  recvcounts_sym_tot = number of values received by each processor from the  other processors fft
!!  rdispls            = positions of values received by the processor from each processor band
!!  rdispls_sym        = positions of values received by the processor from each processor fft
!!  sendcounts         = number of values sended   by the  processor to   each processor band
!!  sendcounts_sym     = number of values sended   by the  processor to   each processor fft
!!  sendcounts_sym_all = number of values sended   by each processor to the other processors fft
!!  sdispls            = postions of values sended by the processor to each processor band
!!  sdispls_sym        = postions of values sended by the processor to each processor fft
!!  tab_proc           = positions of opposited planewave coordinates in the list of the
!!                       processors fft
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Note that in case of band parallelism, the number of spin-up
!! and spin-down bands must be equal at each k points
!!
!! PARENTS
!!      gstate,inwffil
!!
!! CHILDREN
!!      flush_unit,kpgio,leave_new,sphereboundary,wrtout,xallgather_mpi
!!      xallgatherv_mpi,xalltoallv_mpi,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_kpgio(accesswff,ecut_eff,exchn2n3d,gmet,istwfk,kg,kptns,fnametmp_kg,mgfft,mkmem,mode_paral,&
&                     mpi_enreg,mpw,nband,nkpt,npwarr,npwtot,nsppol,unkg)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_io_tools, only : flush_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_kpgio'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_53_ffts
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,exchn2n3d,mgfft,mkmem,mpw,nkpt,nsppol,unkg
 real(dp),intent(in) :: ecut_eff
 character(len=4),intent(in) :: mode_paral
 character(len=fnlen),intent(in) :: fnametmp_kg
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol)
 integer,intent(out) :: kg(3,mpw*mkmem),npwarr(nkpt),npwtot(nkpt)
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: comm_band,comm_fft,idatarecv,idatarecv0,ierr,ikg,ikpt,ikpt_this_proc,iproc,isppol,istwf_k,jsendloc,me_band
 integer :: me_fft,me_kpt,ndatarecv,ndatarecv_tot,ndatasend_sym,nproc_band,nproc_fft,npw_k,npw_tot
 character(len=500)  :: message
!arrays
 integer,allocatable :: gbound(:,:)
 integer,allocatable :: kg_k(:,:),kg_k_gather(:,:),kg_k_gather_all(:,:),kg_k_gather_send(:,:),kg_k_gather_sym(:,:)
 integer,allocatable :: npw_per_proc(:)
 integer,allocatable :: rdispls(:),rdispls_all(:),rdispls_sym(:),rdispls_sym_loc(:)
 integer,allocatable :: recvcounts(:),recvcounts_sym(:),recvcounts_sym_loc(:),recvcounts_sym_tot(:)
 integer,allocatable :: sdispls(:),sdispls_sym_loc(:),sdispls_sym(:)
 integer,allocatable :: sendcounts(:),sendcounts_sym(:),sendcounts_sym_all(:),sendcounts_sym_loc(:)
 integer,allocatable :: sum_kg(:),tab_proc(:)



! *********************************************************************

!DEBUG
!write(std_out,*)' prep_kpgio : enter '
!ENDDEBUG

!---------------------------------------------
!Initialisation
!---------------------------------------------
 nproc_fft    = mpi_enreg%nproc_fft
 nproc_band   = mpi_enreg%nproc_band

 me_band      = mpi_enreg%me_band
 me_fft       = mpi_enreg%me_fft
 me_kpt       = mpi_enreg%me_kpt

 comm_band    = mpi_enreg%comm_band
 comm_fft     = mpi_enreg%comm_fft


!=============================================================================
!Initialize all various tabs within the mpi_enreg%bandfft_kpt(ikpt) data_struc
!These tabs are distributed over the kpt processors
!=============================================================================
 ABI_ALLOCATE(mpi_enreg%bandfft_kpt,(mkmem))
 ABI_ALLOCATE(mpi_enreg%tab_kpt_distrib,(nkpt))
 mpi_enreg%bandfft_kpt(:)%flag1_is_allocated=0
 mpi_enreg%bandfft_kpt(:)%flag2_is_allocated=0
 mpi_enreg%bandfft_kpt(:)%flag3_is_allocated=0
 mpi_enreg%tab_kpt_distrib(:)=0
 do isppol=1,nsppol
   ikpt_this_proc=0
   do ikpt=1,nkpt
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,:,isppol)-me_kpt))/=0) then
       cycle
     end if
     ikpt_this_proc=ikpt_this_proc+1
!    This test should be done when dataset are read and slipt of work do between processor
!    If this test is not good for one proc then other procs fall in deadlock->so PERS and MPI_ABORT
     if (ikpt_this_proc > mkmem) then
       write(message, '(a,a,a,a)' ) ch10,' gstate :  BUG -',ch10,&
&       ' this bandfft tab cannot be allocated !'
       call wrtout(std_out,message,'PERS')
       call flush_unit(std_out)
       call leave_new('PERS')
     end if
     mpi_enreg%tab_kpt_distrib(ikpt)=ikpt_this_proc
     ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%gbound    ,(2*mgfft+8,2))
     ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts,(nproc_band))
     ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts,(nproc_band))
     ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls   ,(nproc_band))
     ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls   ,(nproc_band))
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag1_is_allocated=1
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather)

!    Initialize various quantities which will be computed in vtorho
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ffnl_gather)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kinpw_gather)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kpg_k_gather)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ph3d_gather)

!    Initialize various quantities which will be computed below in case of istwf_k=2
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls_sym)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym_tot)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls_sym)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym_all)
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%tab_proc)

!    In case of MPI_IO, initialize the following tab
     nullify(mpi_enreg%bandfft_kpt(ikpt_this_proc)%ind_kg_mpi_to_seq)
   end do
 end do
 if((mpi_enreg%paral_compil_mpio==1).and.accesswff==IO_MODE_MPI) then
   mpi_enreg%flag_ind_kg_mpi_to_seq = 1
 end if

!=============================================================================
!End of the initialization of the mpi_enreg%bandfft_kpt(ikpt) data_struc
!=============================================================================

!Compute kg, npwarr and npwtot
 call kpgio(ecut_eff,exchn2n3d,gmet,istwfk,kg,fnametmp_kg, &
& kptns,mkmem,nband,nkpt,mode_paral,mpi_enreg,&
& mpw,npwarr,npwtot,nsppol,unkg)

!=============================================================================
!Compute and store various tabs in mpi_enreg%bandfft_kpt(ikpt) data_struc
!These ones will be used in following subroutines:
!vtorho, mkrho, prep_nonlop, prep_fourwf, prep_getghc...
!=============================================================================
 ABI_ALLOCATE(sdispls       ,(nproc_band))
 ABI_ALLOCATE(sendcounts    ,(nproc_band))
 ABI_ALLOCATE(rdispls       ,(nproc_band))
 ABI_ALLOCATE(recvcounts    ,(nproc_band))

 do isppol=1,nsppol
   ikg=0
   do ikpt=1,nkpt
     npw_k=npwarr(ikpt)
     istwf_k=istwfk(ikpt)
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,:,isppol)-me_kpt))/=0) then
       cycle
     end if
     ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)
     if ((ikpt_this_proc > mkmem).or.(ikpt_this_proc==0)) then
       write(message, '(a,a,a,a)' ) ch10,' gstate :  BUG -',ch10,&
&       ' this bandfft tab is not allocated !'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if

     call xallgather_mpi(npw_k,recvcounts,comm_band,ierr)
     rdispls(1)=0
     do iproc=2,nproc_band
       rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
     end do
     ndatarecv=rdispls(nproc_band)+recvcounts(nproc_band)

     ABI_ALLOCATE(kg_k_gather,(3,ndatarecv))
     ABI_ALLOCATE(kg_k,(3,mpw))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     call xallgatherv_mpi(kg_k,3*npw_k,kg_k_gather,3*recvcounts(:),3*rdispls(:),comm_band,ierr)

     sendcounts(:)=npw_k*mpi_enreg%bandpp
     do iproc=1,nproc_band
       sdispls(iproc)=(iproc-1)*npw_k*mpi_enreg%bandpp
     end do

!    ============================================================================
!    Here we compute gbound, as well for istwf_k=1 as for istwf_k=2 and store it
!    ============================================================================
     ABI_ALLOCATE(npw_per_proc,(nproc_fft))
     ABI_ALLOCATE(rdispls_all,(nproc_fft))
     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     if (mgfft>0) gbound(:,:)=0
     if (istwf_k==1) then
       call xallgather_mpi(ndatarecv,npw_per_proc,mpi_enreg%comm_fft,ierr)
       rdispls_all(1)=0
       do iproc=2,nproc_fft
         rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
       end do
       npw_tot=rdispls_all(nproc_fft)+npw_per_proc(nproc_fft)
       ABI_ALLOCATE(kg_k_gather_all,(3,npw_tot))
       call xallgatherv_mpi(kg_k_gather,&
&       3*ndatarecv,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,mpi_enreg%comm_fft,ierr)
       if (mgfft>0) then
         call sphereboundary(gbound,istwf_k,kg_k_gather_all,mgfft,npw_tot)
       end if

     else if (istwf_k==2) then

!      ============================================================================
!      In this case, we have to add the opposite values in the kg_k_gather tab
!      before computing gbound
!      ============================================================================

!      Allocation
       ABI_ALLOCATE(tab_proc          ,(ndatarecv))
       ABI_ALLOCATE(sendcounts_sym    ,(nproc_fft))
       ABI_ALLOCATE(sendcounts_sym_all,(nproc_fft*nproc_fft))
       ABI_ALLOCATE(sdispls_sym       ,(nproc_fft))
       ABI_ALLOCATE(recvcounts_sym    ,(nproc_fft))
       ABI_ALLOCATE(recvcounts_sym_tot,(nproc_fft))
       ABI_ALLOCATE(rdispls_sym       ,(nproc_fft))
       ABI_ALLOCATE(sendcounts_sym_loc,(nproc_fft))
       ABI_ALLOCATE(sdispls_sym_loc   ,(nproc_fft))
       ABI_ALLOCATE(recvcounts_sym_loc,(nproc_fft))
       ABI_ALLOCATE(rdispls_sym_loc   ,(nproc_fft))

!      Initialisation
       tab_proc(:)            = 0
       sendcounts_sym(:)      = 0
       sendcounts_sym_all(:)  = 0
       sdispls_sym(:)         = 0
       recvcounts_sym(:)      = 0
       recvcounts_sym_tot(:)  = 0

!      Localisation of kg_k==[0 0 0]
       ABI_ALLOCATE(sum_kg,(ndatarecv))
       idatarecv0    = -1
       ndatasend_sym = ndatarecv
       sum_kg=sum(abs(kg_k_gather),1)
       if (count(sum_kg==0)/=0) then
         do idatarecv=1,ndatarecv
           if (sum_kg(idatarecv)==0) idatarecv0=idatarecv
         end do
         ndatasend_sym = ndatarecv-1
       end if

!      Localisation of the processor where the vector -k2 is
       do idatarecv=1,ndatarecv
         if (idatarecv/=idatarecv0) then
           tab_proc(idatarecv)   = modulo(-kg_k_gather(2,idatarecv),nproc_fft)
         else
           tab_proc(idatarecv) = -1
         end if
       end do

!      Number of values send by the processor to the others
       do iproc=1,nproc_fft
         sendcounts_sym(iproc) = count(tab_proc(:)==(iproc-1))
       end do

!      Save sendcounts_sym for each processor in sendcounts_sym_all
!      knowed by all processors of comm_fft
       rdispls_sym(1)=0
       do iproc=2,nproc_fft
         rdispls_sym(iproc)= nproc_fft*(iproc-1)
       end do
       recvcounts_sym(:)=nproc_fft
       call xallgatherv_mpi(sendcounts_sym(:),nproc_fft,&
       sendcounts_sym_all(:),recvcounts_sym,rdispls_sym,comm_fft,ierr)

!      Calculation of the dimension of kg_k_gather_sym for each processor
!      recvcounts_sym_tot is knowed by all processors of comm_fft
       call xsum_mpi(sendcounts_sym,recvcounts_sym_tot,nproc_fft,comm_fft,ierr)

!      Dimension of kg_k_gather_sym
       ndatarecv_tot = ndatarecv+recvcounts_sym_tot(me_fft+1)

!      Intialize kg_k_gather_sym
       ABI_ALLOCATE(kg_k_gather_sym,(3,ndatarecv_tot))
       kg_k_gather_sym(:,:)=0
       kg_k_gather_sym(:,1:ndatarecv) = kg_k_gather(:,:)

!      Allocation and initialisation
       ABI_ALLOCATE(kg_k_gather_send,(3,ndatasend_sym))
       kg_k_gather_send(:,:)=0

!      The values are sorted in blocks
       jsendloc=0
       do iproc=1,nproc_fft

!        Position of the beginning of the block
         sdispls_sym(iproc)=jsendloc

!        Creation of the blocks
         do idatarecv=1,ndatarecv
           if (tab_proc(idatarecv)==(iproc-1)) then
             jsendloc=jsendloc+1
             kg_k_gather_send(:,jsendloc)  = -kg_k_gather(:,idatarecv)
           end if
         end do
       end do

!      Position of received data
       rdispls_sym(1)= ndatarecv
       recvcounts_sym(1)= sendcounts_sym_all((me_fft+1))
       do iproc=2,nproc_fft
         rdispls_sym(iproc)    = rdispls_sym(iproc-1) + &
         sendcounts_sym_all((me_fft+1)+(iproc-2)*nproc_fft)
         recvcounts_sym(iproc) = sendcounts_sym_all((me_fft+1)+(iproc-1)*nproc_fft)
       end do

!      Exchange of kg_k
       sendcounts_sym_loc = sendcounts_sym*3
       sdispls_sym_loc    = sdispls_sym   *3
       recvcounts_sym_loc = recvcounts_sym*3
       rdispls_sym_loc    = rdispls_sym   *3
       call xalltoallv_mpi(kg_k_gather_send(:,:),sendcounts_sym_loc,sdispls_sym_loc,&
       kg_k_gather_sym(:,:) ,recvcounts_sym_loc,rdispls_sym_loc,comm_fft,ierr)

!      Store the following data in the mpi_enreg%bandfft_kpt data_struc
       ikpt_this_proc=mpi_enreg%tab_kpt_distrib(ikpt)
       ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym,(3,ndatarecv_tot))
       ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls_sym,(nproc_fft))
       ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym,(nproc_fft))
       ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym_tot,(nproc_fft))
       ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls_sym,(nproc_fft))
       ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym,(nproc_fft))
       ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym_all,(nproc_fft*nproc_fft))
       ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%tab_proc,(ndatarecv))
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%flag3_is_allocated=1

       mpi_enreg%bandfft_kpt(ikpt_this_proc)%idatarecv0           =idatarecv0
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%ndatarecv_tot        =ndatarecv_tot
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%ndatasend_sym        =ndatasend_sym
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym(:,:) =kg_k_gather_sym(:,:)
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls_sym(:)       =rdispls_sym(:)
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym(:)    =recvcounts_sym(:)
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts_sym_tot(:)=recvcounts_sym_tot(:)
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls_sym(:)       =sdispls_sym(:)
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym(:)    =sendcounts_sym(:)
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts_sym_all(:)=sendcounts_sym_all(:)
       mpi_enreg%bandfft_kpt(ikpt_this_proc)%tab_proc(:)          =tab_proc(:)

       ABI_DEALLOCATE(tab_proc)
       ABI_DEALLOCATE(sendcounts_sym)
       ABI_DEALLOCATE(sendcounts_sym_all)
       ABI_DEALLOCATE(sdispls_sym)
       ABI_DEALLOCATE(recvcounts_sym)
       ABI_DEALLOCATE(recvcounts_sym_tot)
       ABI_DEALLOCATE(rdispls_sym)
       ABI_DEALLOCATE(kg_k_gather_sym)
       ABI_DEALLOCATE(sendcounts_sym_loc)
       ABI_DEALLOCATE(recvcounts_sym_loc)
       ABI_DEALLOCATE(sdispls_sym_loc)
       ABI_DEALLOCATE(rdispls_sym_loc)
       ABI_DEALLOCATE(kg_k_gather_send)
       ABI_DEALLOCATE(sum_kg)

!      Then compute gbound
       call xallgather_mpi(ndatarecv_tot,npw_per_proc,mpi_enreg%comm_fft,ierr)
       rdispls_all(1)=0
       do iproc=2,nproc_fft
         rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
       end do
       npw_tot=rdispls_all(nproc_fft)+npw_per_proc(nproc_fft)
       ABI_ALLOCATE(kg_k_gather_all,(3,npw_tot))
       call xallgatherv_mpi(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym,&
&       3*ndatarecv_tot,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,mpi_enreg%comm_fft,ierr)
       if (mgfft>0) then
         call sphereboundary(gbound,istwf_k,kg_k_gather_all,mgfft,npw_tot)
       end if

!      Only calculations with istwfk=1 or 2
     else
       write(message, '(a,a,a,a,i2,a)' ) ch10,' gstate :  BUG -',ch10,&
&       ' the value istwfk=',istwf_k,' is not allowed in case of bandfft parallelization!'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     end if
     ABI_DEALLOCATE(kg_k_gather_all)
     ABI_DEALLOCATE(npw_per_proc)
     ABI_DEALLOCATE(rdispls_all)
!    ============================================================================
!    End of gbound
!    ============================================================================

!    Tabs which are common to istwf_k=1 and 2
     ABI_ALLOCATE(mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather,(3,ndatarecv))
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%recvcounts(:)   =recvcounts(:)
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%sendcounts(:)   =sendcounts(:)
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%rdispls(:)      =rdispls(:)
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%sdispls(:)      =sdispls(:)
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%gbound(:,:)     =gbound(:,:)
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%kg_k_gather(:,:)=kg_k_gather(:,:)
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%ndatarecv       =ndatarecv
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%istwf_k         =istwf_k
     mpi_enreg%bandfft_kpt(ikpt_this_proc)%npw_tot         =npw_tot
     ABI_DEALLOCATE(kg_k_gather)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(gbound)

     ikg=ikg+npw_k
   end do
 end do
 ABI_DEALLOCATE(recvcounts)
 ABI_DEALLOCATE(sendcounts)
 ABI_DEALLOCATE(rdispls)
 ABI_DEALLOCATE(sdispls)
!=============================================================================
!End of computation and storage of the mpi_enreg%bandfft_kpt(ikpt) data_struc
!=============================================================================

end subroutine prep_kpgio
!!***
