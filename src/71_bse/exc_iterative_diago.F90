!{\src2tex{textfont=tt}}
!!****f* ABINIT/exc_iterative_diago
!!
!! NAME
!!  exc_iterative_diago
!!
!! FUNCTION
!!  Calculates eigenvalues and eigenvectors of the Hermitian excitonic Hamiltonian (coupling is neglected).
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  Bsp
!!    %nreh=Rank of the resonant block of the Hamiltonian.
!!    %nstates=Number of eigenstates required.
!!    %nline=Max number of line minimizations.
!!    %tolwfr=Tolerance on the residuals.
!!    %nbdbuf
!!    %nstep
!!  BS_files<excfiles>=Datatype storing names and files used in the Bethe-Salpeter code.
!!    %exh=Name of the file storing the excitonin resonant part.
!!    %out_eig_out=Name of the file where final results are store.
!!    %in_eig=Name of the file used to initialize the calculation.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Eigenvalues and eigenvectors are written on file %out_eig
!!
!! NOTES
!!  Concernig the approach followed to parallelize this routine: the most important
!!  bottleneck is represent by the storage of the excitonic Hamiltonian since a
!!  large number of k-points is needed to obtain converged exciton energies.
!!  The number of eigenstates is usually much smaller than the rank of the full matrix,
!!  this is especially true if we are only interested in the binding energy of the
!!  exciton or in the excitonic states close to the single-particle gap.
!!  Therefore good scaling and good performance should be obtained by distributing
!!  the row of the excitonic Hamiltonian among the nodes while the required
!!  eigenvectors are duplicated on each node. The memory needed to stores the eigenvalues
!!  scales like nreh*nstates where nreh is the rank of the Hamiltonian and this might
!!  render the calculation unfeasible when nstates is large.
!!  On the other hand, having the complex set of trial eigenvectors on each node permits to parallelize
!!  tasks such as the application of the Hamiltonian as well as the orthogonalization or the sub-space rotation
!!  the later two algorithms represent the most CPU demanding part in standard KS calculations
!!  as they scale with the third power of the number of atoms.
!!  The conjugate direction and the gradient as well as Hphi are not distributed as the line minimization
!!  requires the evalauation of <cg_dir_|H_exc|cg_dir>.
!!  Note that this routine has been written having in mind an homogeneous network of machines.
!!  A network made of different CPU will lead to unpredictable results as each node has
!!  to check for the converge of the calculation.
!!
!! PARENTS
!!      exc_diago_driver
!!
!! CHILDREN
!!      wrtout,zhegv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

!#define DEV_MG_DEBUG_THIS

subroutine exc_iterative_diago(BSp,BS_files,Hdr_bse,prtvol,comm)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_bs_defs
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools,   only : get_unit, flush_unit
 use m_abilasi,    only : xheev, xhpev
 use m_header,     only : hdr_clean, hdr_mpio_skip
 use m_bse_io,     only : exc_read_rcblock

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_iterative_diago'
 use interfaces_14_hidewrite
 use interfaces_71_bse, except_this_one => exc_iterative_diago
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,prtvol
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) ::  BS_files
 type(Hdr_type),intent(in) :: Hdr_bse

!Local variables ------------------------------
!scalars
 integer,parameter :: STRICT=2,MEDIUM=1,WORST=0
 integer(i8b) :: bsize_hmat,bsize_phi_block
 integer :: hexc_size,exc_nst,nline,nbdbuf,nstep
 integer :: cg_nsteps,my_t1,my_t2,my_nt
 integer :: ii,jj,state,line
 integer :: max_nline,nsppol
 integer :: istat,cg_step,nbdbuf_
 integer :: ierr,nproc,my_rank,master
 real(dp) :: exc_gap,exc_maxene,norm,etrial,etrial_old,deltae
 real(dp) :: tolwfr
! real(dp) :: deold
 real(dp) :: dhd,dhc,den
 real(dp) :: fac,poly,xx
 real(dp) :: root,swap,tan2th,diff
 real(dp) :: tolwfr_
 !complex(dpc) :: cg_gamma,dotgg,old_dotgg
 real(dp) :: cg_gamma,dotgg,old_dotgg
 real(dp) :: max_resid,costh,sinth
 complex(dpc) :: zz,kprc
 logical :: use_mpio,is_resonant,diago_is_real
 character(len=500) :: msg
 character(len=fnlen) :: hexc_fname,ihexc_fname,oeig_fname
!arrays
 integer :: nline_for(Bsp%nstates),convergence_of(Bsp%nstates)
 real(dp) :: resid(Bsp%nstates),exc_energy(Bsp%nstates),rbuf2(2)
! real(dp),allocatable :: gsc(:,:)
! real(dp),allocatable :: cg(:,:)
 !complex,allocatable :: hexc(:,:)
 complex(dpc),allocatable :: hexc(:,:),hji(:),vec_tmp(:)
 complex(dpc),pointer :: my_phi(:)
 real(dp),allocatable :: hexc_diagonal(:)
 complex(dpc),target,allocatable :: phi_block(:,:)
 complex(dpc),allocatable :: hphi(:) !,buffer_dpc(:)
 complex(dpc),allocatable :: cg_dir(:),grad_dir(:),prc_dir(:)
 complex(dpc),allocatable :: old_cg_dir(:)
 !type(MPI_type) :: MPI_enreg_seq

!************************************************************************

 DBG_ENTER("COLL")

!* Fake MPI_type for the sequential part.
! call initmpi_seq(MPI_enreg_seq)

 nsppol = Hdr_bse%nsppol
 MSG_WARNING("nsppol==2 with cg method is still under development")

 if (Bsp%use_coupling>0) then
   MSG_ERROR("CG Method does not support couplng")
 end if

 nproc   = xcomm_size(comm)
 my_rank = xcomm_rank(comm)
 master=0

 use_mpio=.FALSE.
#ifdef HAVE_MPI_IO
 use_mpio = (nproc > 1)
#endif
 !use_mpio=.FALSE.
 use_mpio = .TRUE.

 hexc_size    = SUM(Bsp%nreh)
 ABI_CHECK(hexc_size>=nproc,"hexc_size<nproc!")

 exc_nst= Bsp%nstates
 nline  = Bsp%nline
 nbdbuf = Bsp%nbdbuf
 nstep  = Bsp%niter
 tolwfr = Bsp%cg_tolwfr

 ! Divide the columns of the Hamiltonian among the nodes.
 call xmpi_split_work(hexc_size,comm,my_t1,my_t2,msg,ierr)
 if (ierr/=0) then
   MSG_WARNING(msg)
 end if

 my_nt = my_t2-my_t1+1
 write(msg,'(a,i0,a)')" Will handle ",my_nt," columns of the excitonic Hamiltoninan. "
 call wrtout(std_out,msg,"PERS")

 tolwfr_ = tolwfr
 if (tolwfr < 10**(-30)) then
   tolwfr_ = tol12
   write(msg,'(2(a,es12.4))')" Input tolwfr= ",tolwfr," Using tolwfr= ",tolwfr_
   MSG_WARNING(msg)
 end if

 cg_nsteps = nstep
 if (cg_nsteps<=0) then
   cg_nsteps = 30
   write(msg,'(2(a,es12.4))')" Input nstep= ",nstep," Using cg_nsteps= ",cg_nsteps
   MSG_WARNING(msg)
 end if

 nbdbuf_ = nbdbuf
 if (nbdbuf<=0) then
   nbdbuf_ = 4
   write(msg,'(2(a,i0))')" Input nbdbuf= ",nbdbuf," Using nbdbuf= ",nbdbuf_
   MSG_WARNING(msg)
 end if

 write(std_out,*)" cg_nsteps ",cg_nsteps
 write(std_out,*)" nstates   ",exc_nst
 write(std_out,*)" nline     ",nline
 write(std_out,*)" nbdbuf_   ",nbdbuf_
 write(std_out,*)" tolwfr_   ",tolwfr_
 !
 bsize_hmat = 2*dpc*hexc_size*my_nt
 write(msg,'(a,f8.1,a)')' Allocating excitonic Hamiltonian. Memory requested: ',bsize_hmat*b2Mb,' Mb.'
 call wrtout(std_out,msg,"COLL")

 ABI_ALLOCATE(hexc_diagonal,(my_t1:my_t2))
 ABI_ALLOCATE(hexc,(hexc_size,my_t1:my_t2))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out of memory: excitonic hamiltonian')
 !
 ! Read and construct full excitonic Hamiltonian using Hermiticity.
 if (BS_files%in_hreso /= BSE_NOFILE) then
   hexc_fname = BS_files%in_hreso
 else
   hexc_fname = BS_files%out_hreso
 end if
 !
 ! Read the resonant block from file.
 is_resonant=.TRUE.
 diago_is_real=(.not.BSp%have_complex_ene)
 call exc_read_rcblock(hexc_fname,Bsp,is_resonant,diago_is_real,nsppol,Bsp%nreh,hexc_size,my_t1,my_t2,hexc,use_mpio,comm)
 !
 ! Save diagonal part for preconditioning.
 do jj=my_t1,my_t2
   hexc_diagonal(jj) = REAL(hexc(jj,jj),kind=dp)
 end do
 !
 ! === Initialisation of the excitonic wavefunctions ===
 ! Two cases are possible.
 !   1) Fill trial eigenvectors with random numbers
 !      One needs to initialize wfs in such a way to avoid symmetry traps,
 !      and to avoid linear dependencies between wavefunctions
 !   2) Read eigenstates generated by a previous calculation.

 bsize_phi_block = 2*spc*my_nt*exc_nst
 write(msg,'(a,f8.1,a)')' Allocating BSE eigenvectors. Memory requested: ',bsize_phi_block*b2Mb,' Mb.'
 call wrtout(std_out,msg,"COLL")
 call flush_unit(std_out)

 ABI_ALLOCATE(phi_block,(my_t1:my_t2,exc_nst))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0," out-of-memory phi_block")

 ihexc_fname = ""
 if ( BS_files%in_eig /= BSE_NOFILE ) ihexc_fname = BS_files%in_eig

 call exc_init_phi_block(ihexc_fname,use_mpio,comm)
 !
 ! =========================
 ! === Orthogonalization ===
 ! =========================
 call exc_cholesky_ortho()
 call exc_check_phi_block("After ortho")

 ! * Sub-space rotation.
 call exc_subspace_rotation()
 call exc_check_phi_block("After subspace_rotation")
 !
 ! ===========================
 ! ==== Conjugate gradient ===
 ! ===========================
 ABI_ALLOCATE(hphi,(hexc_size))
 ABI_ALLOCATE(cg_dir,(hexc_size))
 ABI_ALLOCATE(old_cg_dir,(hexc_size))
 ABI_ALLOCATE(grad_dir,(hexc_size))
 ABI_ALLOCATE(prc_dir,(hexc_size))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out-of-memory in cg vectors")

 max_nline=nline
 resid(:)=HUGE(one); nline_for(1:exc_nst)=max_nline; convergence_of(1:exc_nst)=WORST

 do cg_step=1,cg_nsteps
   !
   do state=1,exc_nst
     !
     if (prtvol>=10) then ! Tell us what is going on:
       write(msg,'(a,i6,2x,a,i3,a)')' --- exc_iterative_diago is called for state ',state,'for',nline_for(state),' lines'
       call wrtout(std_out,msg,'PERS')
     end if

     my_phi => phi_block(my_t1:my_t2,state)  ! Extraction of the vector that is iteratively updated.

     do line=1,nline_for(state)
       hphi = czero
       hphi = MATMUL(hexc, my_phi)  ! Compute etrial=<phi|H|phi> and the residual [H-etrial]|phi>.
       call xsum_mpi(hphi,comm,ierr)

       etrial = DOT_PRODUCT(my_phi, hphi(my_t1:my_t2))
       call xsum_mpi(etrial,comm,ierr)
       exc_energy(state) =  etrial

       grad_dir(my_t1:my_t2) = hphi(my_t1:my_t2) - etrial*my_phi ! Compute residual (squared) norm.
       resid(state) =  DOT_PRODUCT(grad_dir(my_t1:my_t2), grad_dir(my_t1:my_t2))
       call xsum_mpi(resid(state),comm,ierr)
       convergence_of(state) = convergence_degree(resid(state))
       !
       ! Check that etrial is decreasing on succeeding lines:
       if (line>1 .and. (etrial > etrial_old+tol12)) then
         write(msg,'(a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&          ' New trial exc_energy at line ',line,' = ',etrial,ch10,&
&          ' is higher than former:',etrial_old,ch10
         MSG_WARNING(msg)
       end if
       etrial_old = etrial
       !
       ! If residual sufficiently small stop line minimization.
       if (convergence_of(state)==STRICT) then
         if (prtvol>=10) then
           write(msg,'(a,i4,a,i2,a,es12.4)')&
&            ' exc_iterative_diago: state ',state,' converged after ',line,&
&            ' line minimizations : resid =',resid(state)
           call wrtout(std_out,msg,'PERS')
         end if
         EXIT !line
       end if

       ! === PROJECT THE STEEPEST DESCENT DIRECTION OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ===
       ! The following projection over the subspace orthogonal to occupied bands
       ! is optional. It is a bit more accurate, but doubles the number of N^3 ops.
       ! It is done only if ortalg>=0.

       ! Project the steepest descent direction: direc(2,npw)=<G|H|Cnk> - \sum_{(i<=n)} <G|H|Cik> , normalized.

       ! Grad_dir is already orthogonal to this band
       ABI_ALLOCATE(hji,(exc_nst))
       hji=czero

       ! MG TODO Don't know why here we sum over i=<=n!!!!!!!!
       do jj=1,exc_nst
         if (jj/=state) hji(jj) = DOT_PRODUCT(phi_block(:,jj), hphi(my_t1:my_t2) )
       end do
       call xsum_mpi(hji,comm,ierr)

       do jj=1,exc_nst
         if (jj/=state) grad_dir(my_t1:my_t2) = grad_dir(my_t1:my_t2) - hji(jj)*phi_block(:,jj)
       end do
       ABI_DEALLOCATE(hji)
       !
       ! === PRECONDITION THE STEEPEST DESCENT DIRECTION ===
       den = DOT_PRODUCT(grad_dir(my_t1:my_t2), hexc_diagonal(my_t1:my_t2)*grad_dir(my_t1:my_t2) )
       call xsum_mpi(den,comm,ierr)

       do ii=my_t1,my_t2
         xx = hexc_diagonal(ii)/den ! Teter polynomial ratio, modified according to Kresse, Furthmuller, PRB 54, 11169 (1996)
         poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
         fac=poly/(poly+16._dp*xx**4)
         kprc = fac*four/(three*den)
         prc_dir(ii) = kprc * grad_dir(ii)
       end do
       !
       ! * PROJECT THE PRECOND. STEEPEST DESCENT DIRECTION OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS.
       ABI_ALLOCATE(hji,(exc_nst))
       hji=czero
       do jj=1,exc_nst
         hji(jj) = DOT_PRODUCT(phi_block(:,jj), prc_dir(my_t1:my_t2) )
       end do
       call xsum_mpi(hji,comm,ierr)

       do jj=1,exc_nst
         prc_dir(my_t1:my_t2) = prc_dir(my_t1:my_t2) - hji(jj)*phi_block(:,jj)
       end do
       ABI_DEALLOCATE(hji)
       !
       ! === COMPUTE THE CONJUGATE-GRADIENT ===
       dotgg = DOT_PRODUCT(prc_dir(my_t1:my_t2),grad_dir(my_t1:my_t2))
       call xsum_mpi(dotgg,comm,ierr)

       if (line==1) then ! At first iteration, cg_gamma is set to zero
         cg_gamma=zero
         old_dotgg=dotgg
         cg_dir = prc_dir
         old_cg_dir = cg_dir
       else
         cg_gamma=dotgg/old_dotgg
         old_dotgg=dotgg
         !write(std_out,*)"cg_gamma= ",cg_gamma
         !cg_dir = prc_dir + cg_gamma*cg_dir
         cg_dir = prc_dir + cg_gamma*old_cg_dir !TODO check this, anyhow it is much faster.
         old_cg_dir =cg_dir  ! old_cg_dir is used to store the previsou CG direction, cg_dir will be orthonormalized to the band
       end if
       !
       ! === PROJECTION OF THE CONJUGATED GRADIENT ===
       zz = DOT_PRODUCT(my_phi, cg_dir(my_t1:my_t2))
       call xsum_mpi(zz,comm,ierr)
       cg_dir(my_t1:my_t2) = cg_dir(my_t1:my_t2) -zz*my_phi(:)

       norm = DOT_PRODUCT(cg_dir(my_t1:my_t2), cg_dir(my_t1:my_t2) )
       call xsum_mpi(norm,comm,ierr)
       norm = SQRT(norm)
       cg_dir = cg_dir/norm ! Have to normalize it.

       ! Line minimization of the Raileigh functional.
       ABI_ALLOCATE(vec_tmp,(hexc_size))
       vec_tmp=czero
       vec_tmp = MATMUL(hexc, cg_dir(my_t1:my_t2))
       call xsum_mpi(vec_tmp,comm,ierr)

#ifdef DEV_MG_DEBUG_THIS
       if (my_rank==master) then
         write(777,*)"cg_step, state, line",cg_step, state, line
         write(777,*)vec_tmp
       end if
#endif

       dhd = DOT_PRODUCT( cg_dir(my_t1:my_t2), vec_tmp(my_t1:my_t2))  ! is this always real?
       dhc = REAL( DOT_PRODUCT( my_phi, vec_tmp(my_t1:my_t2) ))
       ABI_DEALLOCATE(vec_tmp)

       rbuf2 = (/dhd,dhc/)
       call xsum_mpi(rbuf2,comm,ierr)
       dhd = rbuf2(1)
       dhc = rbuf2(2)

#ifdef DEV_MG_DEBUG_THIS
       write(201*(my_rank+1),*)"cg_step, state, line dotgg dhd dhc",cg_step,state,line,dotgg,dhd,dhc
#endif

! Taken from cgwf
       ! Compute tan(2 theta),sin(theta) and cos(theta)
       tan2th=2.0_dp*dhc/(etrial-dhd)

       if (abs(tan2th)<1.d-05) then
         costh=1.0_dp-0.125_dp*tan2th**2
         sinth=0.5_dp*tan2th*(1.0_dp-0.375_dp*tan2th**2)
         ! Check that result is above machine precision
         ! FIXME  This part is not safe on clusters made of different machines or different compiation options.
         if (abs(sinth)<epsilon(0._dp)) then
           write(msg, '(a,es16.4)' ) ' exc_iterative_diago: converged with tan2th= ',tan2th
           call wrtout(std_out,msg,'PERS')
           EXIT !Exit from the loop on line
         end if

       else
         root =sqrt(1.0_dp+tan2th**2)
         costh=sqrt(0.5_dp+0.5_dp/root)
         sinth=sign(sqrt(0.5_dp-0.5_dp/root),tan2th)
       end if
       !
       ! Check for lower of two possible roots (same sign as curvature at theta where slope is zero)
       diff=(etrial-dhd)
       if (diff>zero) then !   Swap c and d if value of diff is positive
         swap=costh
         costh=-sinth
         sinth=swap
         if (prtvol<0 .or. prtvol>=10) then
           write(msg,'(a,2i4)')' exc_iterative_diago: swap roots, line,diff= ',line,diff
           call wrtout(std_out,msg,'PERS')
         end if
       end if
       !
       ! === GENERATE NEW |wf>, H|wf>  =============
       my_phi = costh*my_phi + sinth*cg_dir(my_t1:my_t2)
#ifdef DEV_MG_DEBUG_THIS
       write(100*(my_rank+1),*)"cg_step state, line costh sinth etrial",cg_step,state,line,costh,sinth,etrial
#endif
!end  taken from cgwf

       !norm = SQRT( DOT_PRODUCT(my_phi,my_phi) )
       !my_phi = my_phi /norm
       !write(std_out,*)norm
       !write(std_out,*)DOT_PRODUCT(hphi,my_phi),cos(theta_min)

       ! ======================================================================
       ! =========== CHECK CONVERGENCE AGAINST TRIAL ENERGY ===================
       ! ======================================================================
       ! Compute delta(E)
       !deltae=chc*(costh**2-1._dp)+dhd*sinth**2+2._dp*costh*sinth*dhc
       deltae=etrial*(costh**2-1._dp)+dhd*sinth**2+2._dp*costh*sinth*dhc

       ! Check convergence and eventually exit
       ! if (line==1) then
       !   deold=deltae
       ! else if (abs(deltae)<0.005_dp*abs(deold) .and. line/=nline_for(state))then
       !  if (prtvol>=10)then
       !   write(msg, '(a,i4,1x,a,1p,e12.4,a,e12.4,a)' ) &
       !&   ' cgwf: line',line,&
       !&   ' deltae=',deltae,' < 0.005*',deold,' =>skip lines'
       !   call wrtout(std_out,msg,'PERS')
       !  end if
       !  exc_energy(state) = exc_energy(state) + deltae
       !  EXIT
       ! end if
     end do ! LOOP FOR A GIVEN BAND. Note that there are three "exit" instructions inside
     ! Modify nline_for(state) according to converge degree.
     !if (convergence_of(state) == STRICT) nline_for(state) = MAX(max_nline-2,2)
     !if (convergence_of(state) == MEDIUM) nline_for(state) = MAX(max_nline-1,2)
     !if (convergence_of(state) == WORST ) nline_for(state) = max_nline
   end do !state

   if (prtvol>2) then
     do ii=0,(exc_nst-1)/8
       write(msg,'(a,8es10.2)')' res:',(resid(state),state=1+ii*8,MIN(exc_nst,8+ii*8))
       call wrtout(std_out,msg,'COLL')
     end do
     do ii=0,(exc_nst-1)/8
       write(msg,'(a,8es10.2)')' ene:',(exc_energy(state),state=1+ii*8,MIN(exc_nst,8+ii*8))
       call wrtout(std_out,msg,'COLL')
     end do
   end if

#ifdef DEV_MG_DEBUG_THIS
   write(msg,'(a,i0)')"After cg_step: ",cg_step
   call exc_check_phi_block(msg)
#endif

   ! Find largest residual over bands and Print residuals
   max_resid=MAXVAL( resid(:MAX(1,exc_nst-nbdbuf_)) )

   if (max_resid < tolwfr_) then
     write(msg,'(a,i0,2(a,es10.2),a,i0,a)')&
&      " After ",cg_step," iterations, max_resid= ",max_resid," < tolwfr= ",tolwfr_," ( Excluding nbdbuf= ",nbdbuf_,")"
     call wrtout(std_out,msg,'COLL')
     EXIT ! cg_step
   end if

   if (cg_step==1.or.MOD(cg_step,1)==0) then
     call wrtout(std_out," Subspace rotation + exc_cholesky_ortho ","COLL")
     call exc_subspace_rotation()

     call exc_cholesky_ortho()

     !mcg=hexc_size; mgsc=hexc_size; useoverlap=0
     !allocate(cg(2,mcg),gsc(2,mgsc*useoverlap))
     !do ii=1,exc_nst
     ! cg(1,:) = REAL (phi_block(:,ii))
     ! cg(2,:) = AIMAG(phi_block(:,ii))
     ! call fxphas(cg,gsc,0,0,1,mcg,mgsc,MPI_enreg_seq,1,hexc_size,useoverlap)
     ! phi_block(:,ii)=CMPLX(cg(1,:),cg(2,:))
     !end do
     !deallocate(cg,gsc)
   end if

 end do !cg_step

 ! Release some memory before entering RMM-DIIS
 ABI_DEALLOCATE(hphi)
 ABI_DEALLOCATE(cg_dir)
 ABI_DEALLOCATE(old_cg_dir)
 ABI_DEALLOCATE(grad_dir)
 ABI_DEALLOCATE(prc_dir)

 do ii=0,(exc_nst-1)/8
   write(msg,'(a,8es10.2)')' res:',(resid(state),state=1+ii*8,min(exc_nst,8+ii*8))
   call wrtout(std_out,msg,'COLL')
 end do
 do ii=0,(exc_nst-1)/8
   write(msg,'(a,8es10.2)')' ene:',(exc_energy(state),state=1+ii*8,min(exc_nst,8+ii*8))
   call wrtout(std_out,msg,'COLL')
 end do

 if (max_resid > tolwfr_) then
   write(msg,'(2a,i5,2a,2(a,es10.2),a,i3,a)')ch10,&
&    " WARNING: conjugate-gradient not converged after ",cg_step," iterations.",ch10,&
&    " max_resid= ",max_resid," > tolwfr= ",tolwfr_," ( Excluding nbdbuf= ",nbdbuf_,")"
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 exc_gap    = MINVAL(exc_energy)
 exc_maxene = MAXVAL(exc_energy)

 write(msg,'(a,2(a,f7.2,2a))')ch10,&
&  " First excitonic eigenvalue= ",exc_gap*Ha_eV,   " [eV]. ",ch10,&
&  " Last  excitonic eigenvalue= ",exc_maxene*Ha_eV," [eV]. ",ch10
 call wrtout(std_out,msg,"COLL")
 call wrtout(ab_out,msg,"COLL")

 call exc_check_phi_block("END OF CONJUGATE-GRADIENT")

 ! RMM-DIIS Algorithm.
 !do state=1,exc_nst
 ! call rmm_diis_for(state)
 !end do
 !call exc_cholesky_ortho()
 !
 ABI_DEALLOCATE(hexc)
 ABI_DEALLOCATE(hexc_diagonal)
 !
 ! =====================================
 ! ==== Write final results on file ====
 ! =====================================
 oeig_fname = BS_files%out_eig
 if (oeig_fname== BSE_NOFILE) then
   MSG_WARNING(" oeig_fname was set to "//TRIM(BSE_NOFILE))
   oeig_fname = TRIM(BS_files%out_basename)//"_BSEIG"
   MSG_WARNING("using oeig_fname : "//TRIM(oeig_fname))
 end if

 call exc_write_phi_block(oeig_fname,use_mpio)

 ABI_DEALLOCATE(phi_block)

 call xbarrier_mpi(comm)

 DBG_EXIT("COLL")

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* exc_iterative_diago/exc_init_phi_block
!! NAME
!! exc_init_phi_block
!!
!! FUNCTION
!!  Initialize the eigenstates either from file or fill them with random number
!!  if restart file is not available
!!
!! INPUTS
!!  ihexc_fname=
!!    Name of the file from which the eigenvectors will be read.
!!    Empty string to initialize trial eigenvectors with random numbers.
!!
!! SIDE EFFECTS
!!   phi_block(my_t1:my_t2,exc_nst)=Contains the trial eigenstates.
!!
!! PARENTS
!!      exc_iterative_diago
!!
!! CHILDREN
!!      wrtout,zhegv
!!
!! SOURCE

subroutine exc_init_phi_block(ihexc_fname,use_mpio,comm)

 use m_profiling

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_init_phi_block'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 logical,intent(in) :: use_mpio

!Local variables ------------------------------
 integer :: eig_unt,hexc_size_restart,ii,state,seed
 integer ::  fold1,fold2,foldim,foldre
 character(len=*),intent(in) :: ihexc_fname
!arrays
 complex(dpc),allocatable :: buffer_dpc(:)
#ifdef HAVE_MPI_IO
 integer:: amode,mpi_fh,mpi_err,old_type,etype,eig_type,my_nel,ierr,my_nrows
 integer(XMPI_OFFSET_KIND) :: ehdr_offset,my_offset,my_offpad,fmarker
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecord(:)
 integer :: array_of_sizes(2),array_of_subsizes(2),array_of_starts(2)
#endif

!************************************************************************

 if (LEN_TRIM(ihexc_fname) == 0) then
   call wrtout(std_out," Initializing eigenvectors with random numbers","COLL")
   !
   ! Use random number generator. For portability, use only integer numbers
   ! The series of couples (fold1,fold2) is periodic with a period of
   ! 3x5x7x11x13x17x19x23x29x31, that is, larger than 2**32, the largest integer*4
   ! fold1 is between 0 and 34, fold2 is between 0 and 114. As sums of five
   ! uniform random variables, their distribution is close to a gaussian
   ! the gaussian distributions are folded, in order to be back to a uniform distribution
   ! foldre is between 0 and 20, foldim is between 0 and 18.
   !
   do state=1,exc_nst
     do ii=my_t1,my_t2
       seed=ii+(state-1)*hexc_size ! Different seed for different transitions and bands
       fold1 =mod(seed,3)+mod(seed,5)+mod(seed,7)+mod(seed,11)+mod(seed,13)
       fold2 =mod(seed,17)+mod(seed,19)+mod(seed,23)+mod(seed,29)+mod(seed,31)
       foldre=mod(fold1+fold2,21)
       foldim=mod(3*fold1+2*fold2,19)

       phi_block(ii,state) = DCMPLX(foldre,foldim)
     end do
   end do

 else
   if (.not.use_mpio) then
     call wrtout(std_out," Initializing eigenvectors from file: "//TRIM(ihexc_fname)//" using Fortran IO.","COLL")

     eig_unt=get_unit()
     open(eig_unt,file=ihexc_fname,form='unformatted',status="old")

     read(eig_unt) hexc_size_restart
     ABI_CHECK(hexc_size_restart==hexc_size,"hexc_size_restart /= hexc_size")
     read(eig_unt) !skip DCMPLX(exevl(1:hexc_size))

     ABI_ALLOCATE(buffer_dpc,(hexc_size))
     do ii=1,exc_nst
       read(eig_unt) buffer_dpc
       phi_block(my_t1:my_t2,ii) = buffer_dpc(my_t1:my_t2)
     end do
     ABI_DEALLOCATE(buffer_dpc)

     close(eig_unt)
   else
     call wrtout(std_out," Initializing eigenvectors from file: "//TRIM(ihexc_fname)//" using MPI-IO.","COLL")
#ifdef HAVE_MPI_IO
     !
     ! Open the file with MPI-IO
     amode=MPI_MODE_RDONLY

     call MPI_FILE_OPEN(comm, ihexc_fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
     msg = " MPI_IO error opening file: "//TRIM(ihexc_fname)
     ABI_CHECK_MPI(mpi_err,msg)

     ! Move the file pointer to skip the first two records.
     ehdr_offset=0
     call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,mpi_err)
     write(std_out,*)"fmarker first record ",fmarker
     call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,mpi_err)
     write(std_out,*)"fmarker first record ",fmarker
     !$call hdr_mpio_skip(mpi_fh,fform,ehdr_offset)
     !$ehdr_offset = 4*xmpio_bsize_frm + xmpio_bsize_int + exc_nst*xmpio_bsize_dpc

     etype=MPI_BYTE; old_type=MPI_DOUBLE_COMPLEX

     my_nrows=my_t2-my_t1+1; old_type=MPI_DOUBLE_COMPLEX
     array_of_sizes    = (/hexc_size,exc_nst/)
     array_of_subsizes = (/my_nrows,exc_nst/)
     array_of_starts   = (/my_t1,1/)
     call xmpio_create_fsubarray_2D(array_of_sizes,array_of_subsizes,array_of_starts,old_type,eig_type,my_offpad,mpi_err)
     ABI_CHECK_MPI(mpi_err,"fsubarray_2D")
     !
     ! Each node uses a different offset to skip the header and the blocks written by the other CPUs.
     my_offset = ehdr_offset + my_offpad

     call MPI_FILE_SET_VIEW(mpi_fh, my_offset, etype, eig_type, 'native', MPI_INFO_NULL, mpi_err)
     ABI_CHECK_MPI(mpi_err,"SET_VIEW")

     call MPI_TYPE_FREE(eig_type,mpi_err)
     ABI_CHECK_MPI(mpi_err,"MPI_TYPE_FREE")

     my_nel = my_nrows*exc_nst
     call MPI_FILE_READ_ALL(mpi_fh, phi_block, my_nel, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpi_err)
     ABI_CHECK_MPI(mpi_err,"FILE_WRITE")

     ! It seems that personal calls make the code stuck
     ! check the fortran markers.
     ABI_ALLOCATE(bsize_frecord,(exc_nst))
     bsize_frecord = hexc_size * xmpi_bsize_dpc
     ! ehdr_offset points to the end of the header.
     call xmpio_check_frmarkers(mpi_fh,ehdr_offset,xmpio_at_all,exc_nst,bsize_frecord,ierr)
     ABI_CHECK(ierr==0,"Error in Fortran markers")
     ABI_DEALLOCATE(bsize_frecord)
     !
     ! Close the file.
     call MPI_FILE_CLOSE(mpi_fh, mpi_err)
     ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")
#else
     MSG_ERROR("You should not be here")
#endif
   end if
   !
 end if

end subroutine exc_init_phi_block
!!***

!----------------------------------------------------------------------

!!****f* exc_iterative_diago/exc_write_phi_block
!! NAME
!! exc_write_phi_block
!!
!! FUNCTION
!!  Write phi_block on the Fortran file oeig_fname.
!!
!! PARENTS
!!      exc_iterative_diago
!!
!! CHILDREN
!!      wrtout,zhegv
!!
!! SOURCE

subroutine exc_write_phi_block(oeig_fname,use_mpio)

 use m_profiling

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_write_phi_block'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: oeig_fname
 logical,intent(in) :: use_mpio

!Local variables ------------------------------
 integer :: eig_unt,ios,state,mpi_err !,fform
 character(len=500) :: msg
! type(Hdr_type) :: hexc_Hdr
!!arrays
 complex(dpc),allocatable :: buffer_dpc(:)
#ifdef HAVE_MPI_IO
 integer:: amode,mpi_fh,old_type,etype,eig_type,my_nel,ierr,my_nrows
 integer(XMPI_OFFSET_KIND) :: ehdr_offset,my_offset,my_offpad,fmarker
 integer(XMPI_OFFSET_KIND),allocatable :: bsize_frecord(:)
 integer :: array_of_sizes(2),array_of_subsizes(2),array_of_starts(2)
#endif

!************************************************************************

 if (.not.use_mpio) then
   !
   ! * Master writes the header.
   if (my_rank==master) then
     call wrtout(std_out," Writing eigenstates on file "//TRIM(oeig_fname)//" via Fortran-IO","COLL")
     eig_unt=get_unit()
     open(eig_unt,file=oeig_fname,form='unformatted',iostat=ios)
     msg = "Opening file: "//TRIM(oeig_fname)
     ABI_CHECK(ios==0,msg)
     write(eig_unt) exc_nst
     write(eig_unt) CMPLX(exc_energy(1:exc_nst),kind=dpc)
   end if

   ! Wavefunctions are gathered on the master node band-by-band.
   ! TODO bands should be treated in blocks to minimize the number of MPI calls.
   ABI_ALLOCATE(buffer_dpc,(hexc_size))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out of memory buffer_dpc")
   !
   do state=1,exc_nst
     buffer_dpc=czero
     buffer_dpc(my_t1:my_t2) = phi_block(:,state)
     call xsum_master(buffer_dpc,master,comm,mpi_err)
     if (my_rank==master) write(eig_unt) buffer_dpc(1:hexc_size)
   end do
   ABI_DEALLOCATE(buffer_dpc)

   if (my_rank==master) close(eig_unt)

 else

#ifdef HAVE_MPI_IO
   call wrtout(std_out," Writing eigenstates on file "//TRIM(oeig_fname)//" with MPI-IO","COLL")
   !
   ! Write the header.
   if (my_rank==master) then ! Write header using Fortran primitives.
     eig_unt = get_unit()
     open(unit=eig_unt,file=oeig_fname,form='unformatted')
     write(eig_unt) exc_nst
     write(eig_unt) CMPLX(exc_energy(1:exc_nst),kind=dpc)
     !!  fform = 1002 ! TODO: change setup_bse so that Hdr_bse reflects the parameters of the run.
     !!  call hdr_io_int(fform,Hdr_bse,rdwr2,eig_unt)
     close(eig_unt)
   end if

   call xbarrier_mpi(comm)
   !
   ! Open the file with MPI-IO
   amode=MPI_MODE_RDWR

   call MPI_FILE_OPEN(comm, oeig_fname, amode, MPI_INFO_NULL, mpi_fh, mpi_err)
   msg = " MPI_IO error opening file: "//TRIM(oeig_fname)
   ABI_CHECK_MPI(mpi_err,msg)

   ! Move the file pointer to skip the first two records.
   ehdr_offset=0
   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,mpi_err)
   !write(std_out,*)"fmarker first record ",fmarker
   call xmpio_read_frm(mpi_fh,ehdr_offset,xmpio_at_all,fmarker,mpi_err)
   !write(std_out,*)"fmarker first record ",fmarker
   !$call hdr_mpio_skip(mpi_fh,fform,ehdr_offset)
   !$ehdr_offset = 4*xmpio_bsize_frm + xmpio_bsize_int + exc_nst*xmpio_bsize_dpc

   etype=MPI_BYTE; old_type=MPI_DOUBLE_COMPLEX

   my_nrows=my_t2-my_t1+1; old_type=MPI_DOUBLE_COMPLEX
   array_of_sizes    = (/hexc_size,exc_nst/)
   array_of_subsizes = (/my_nrows,exc_nst/)
   array_of_starts   = (/my_t1,1/)
   call xmpio_create_fsubarray_2D(array_of_sizes,array_of_subsizes,array_of_starts,old_type,eig_type,my_offpad,mpi_err)
   ABI_CHECK_MPI(mpi_err,"fsubarray_2D")
   !
   ! Each node uses a different offset to skip the header and the blocks written by the other CPUs.
   my_offset = ehdr_offset + my_offpad

   call MPI_FILE_SET_VIEW(mpi_fh, my_offset, etype, eig_type, 'native', MPI_INFO_NULL, mpi_err)
   ABI_CHECK_MPI(mpi_err,"SET_VIEW")

   call MPI_TYPE_FREE(eig_type,mpi_err)
   ABI_CHECK_MPI(mpi_err,"MPI_TYPE_FREE")

   my_nel = my_nrows*exc_nst
   call MPI_FILE_WRITE_ALL(mpi_fh, phi_block, my_nel, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpi_err)
   ABI_CHECK_MPI(mpi_err,"FILE_WRITE")

   ! It seems that personal calls make the code stuck
   ABI_ALLOCATE(bsize_frecord,(exc_nst))
   bsize_frecord = hexc_size * xmpi_bsize_dpc
   ! ehdr_offset points to the end of the header.
   call xmpio_write_frmarkers(mpi_fh,ehdr_offset,xmpio_at_all,exc_nst,bsize_frecord,ierr)
   ABI_CHECK(ierr==0,"Error while writing Fortran markers")
   ABI_DEALLOCATE(bsize_frecord)
   !
   ! Close the file.
   call MPI_FILE_CLOSE(mpi_fh, mpi_err)
   ABI_CHECK_MPI(mpi_err,"FILE_CLOSE")
#else
   MSG_ERROR("MPI-IO support not enabled")
#endif
 end if

end subroutine exc_write_phi_block
!!***

!----------------------------------------------------------------------

!!****f* exc_iterative_diago/exc_subspace_rotation
!! NAME
!! exc_subspace_rotation
!!
!! FUNCTION
!!  This routine performs the subspace rotation.
!!
!! SIDE EFFECTS
!!   phi_block(my_t1:my_t2,exc_nst)=Contains the trial eigenstates.
!!
!! PARENTS
!!      exc_iterative_diago
!!
!! CHILDREN
!!      wrtout,zhegv
!!
!! SOURCE

subroutine exc_subspace_rotation()

 use m_profiling

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_subspace_rotation'
!End of the abilint section

 implicit none

!Local variables ------------------------------
 integer :: ii,jj,ipack,ierr
!arrays
 real(dp),allocatable :: sub_ene(:)
! real(dp),allocatable :: evec(:,:)
 complex(dpc),allocatable :: sub_ham(:,:),sub_pham(:),hphi_tot(:)
! complex(dpc),allocatable :: phi_tmp(:,:)

!************************************************************************

 ! * Sub-space rotation. Calculate <phi_i|H|phi_j> in packed form.
 ! TODO: this part can be rewritten using BLAS3 routines.

 ABI_ALLOCATE(hphi_tot,(hexc_size))
 ABI_ALLOCATE(sub_pham,(exc_nst*(exc_nst+1)/2))
 sub_pham=czero; ipack=0

 do jj=1,exc_nst
   hphi_tot = czero
   hphi_tot(:) = MATMUL(hexc, phi_block(:,jj))
   call xsum_mpi(hphi_tot,comm,ierr)

   do ii=1,jj
     ipack=ipack+1
     sub_pham(ipack) = DOT_PRODUCT(phi_block(my_t1:my_t2,ii), hphi_tot(my_t1:my_t2) )
     if (ii==jj) sub_pham(ipack) = REAL(sub_pham(ipack),kind=dp)
   end do
 end do
 call xsum_mpi(sub_pham,comm,ierr)

 ABI_ALLOCATE(sub_ham,(exc_nst,exc_nst))
 sub_ham=czero
 ABI_ALLOCATE(sub_ene,(exc_nst))

 call xhpev("Vectors","Upper",exc_nst,sub_pham,sub_ene,sub_ham,exc_nst) !,comm)

 ABI_DEALLOCATE(hphi_tot)
 ABI_DEALLOCATE(sub_pham)
 ABI_DEALLOCATE(sub_ene)

 !do ii=1,exc_nst
 ! norm = DOT_PRODUCT(sub_ham(:,ii),sub_ham(:,ii))
 ! write(std_out,*)"norm subspac",norm
 ! sub_ham(:,ii) = sub_ham(:,ii)/norm
 !end do

 !allocate(evec(2*exc_nst,exc_nst))

 !do ii=1,exc_nst
 ! do jj=1,exc_nst
 ! evec(jj,  ii) = REAL (sub_ham(jj,ii))
 ! evec(jj+1,ii) = AIMAG(sub_ham(jj,ii))
 ! end do
 !end do

 !call normev(evec,exc_nst,exc_nst)

 !do ii=1,exc_nst
 ! do jj=1,exc_nst
 !  sub_ham(jj,ii) = CMPLX( evec(jj,ii),evec(jj+1,ii) )
 ! end do
 !end do
 !deallocate(evec)

#if 0
 ABI_ALLOCATE(phi_tmp,(my_nt,exc_nst))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out-of-memory in phi_tmp")
 phi_tmp = phi_block

 call ZGEMM('N','N',my_nt,exc_nst,exc_nst,cone,phi_tmp,my_nt,sub_ham,exc_nst,czero,phi_block,my_nt)

 ABI_DEALLOCATE(phi_tmp)
#else
 phi_block = MATMUL(phi_block,sub_ham)
#endif

 ABI_DEALLOCATE(sub_ham)

end subroutine exc_subspace_rotation
!!***

!----------------------------------------------------------------------

!!****f* exc_iterative_diago/exc_cholesky_ortho
!! NAME
!! exc_cholesky_ortho
!!
!! FUNCTION
!!  This routine performs the orthogonalization of the trial eigenstates using the
!!  Cholesky Algorithm.
!!
!! SIDE EFFECTS
!!   phi_block(my_t1:my_t2,exc_nst)=Contains the trial eigenstates.
!!
!! PARENTS
!!      exc_iterative_diago
!!
!! CHILDREN
!!      wrtout,zhegv
!!
!! SOURCE

subroutine exc_cholesky_ortho()

 use m_profiling

 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_cholesky_ortho'
!End of the abilint section

 implicit none

!Local variables ------------------------------
 integer :: my_info,ii,jj,ipack,ierr
!arrays
 complex(dpc),allocatable :: overlap(:,:),povlp(:)

!************************************************************************

 ! 1) overlap_ij =  <phi_i|phi_j>
 ABI_ALLOCATE(overlap,(exc_nst,exc_nst))

#if defined HAVE_BSE_UNPACKED
 overlap = czero

 call ZGEMM('C','N',exc_nst,exc_nst,my_nt,cone,phi_block,my_nt,phi_block,my_nt,czero,overlap,exc_nst)
 call xsum_mpi(overlap,comm,ierr)

 do ii=1,exc_nst
   overlap(ii,ii)=REAL(overlap(ii,ii),kind=dp)
 end do

 ! 2) Cholesky factorization: overlap = U^H U with U upper triangle matrix.
 call ZPOTRF('U',exc_nst,overlap,exc_nst,my_info)
 if (my_info/=0)  then
   write(msg,'(a,i3)')' ZPOTRF returned info= ',my_info
   MSG_ERROR(msg)
 end if

#else

 ! 1) Calculate overlap_ij =  <phi_i|phi_j> in packed form.
 ABI_ALLOCATE(povlp,(exc_nst*(exc_nst+1)/2))
 povlp = czero; ipack=0
 do jj=1,exc_nst
   do ii=1,jj
     ipack=ipack+1
     povlp(ipack) = DOT_PRODUCT( phi_block(my_t1:my_t2,ii), phi_block(my_t1:my_t2,jj) )
     if (ii==jj) povlp(ipack) = REAL(povlp(ipack),kind=dp)
   end do
 end do
 call xsum_mpi(povlp,comm,ierr)

 ! 2) Cholesky factorization: overlap = U^H U with U upper triangle matrix.
 call ZPPTRF("U",exc_nst,povlp,my_info)
 if (my_info/=0)  then
   write(msg,'(a,i3)')' ZPPTRF returned info= ',my_info
   MSG_ERROR(msg)
 end if
 !call xsum_mpi(povlp,comm,ierr)
 !povlp=povlp/nproc

 !unpack povlp to prepare call to ZTRSM.
 ipack=0
 do jj=1,exc_nst
   do ii=1,jj
     ipack=ipack+1
     if (ii/=jj) then
       overlap(ii,jj)=      povlp(ipack)
       overlap(jj,ii)=CONJG(povlp(ipack))
     else
       overlap(ii,ii)=REAL(povlp(ipack),kind=dp)
     end if
   end do
 end do
 ABI_DEALLOCATE(povlp)
#endif

 ! Check if this can be done with Scalapack. Direct PZTRSM is not provided

 ! 3) Solve X U = phi_block, on exit the phi_block treated by this node is orthonormalized.
 !call ZTRSM('R','U','N','N',hexc_size,exc_nst,cone,overlap,exc_nst,phi_block,hexc_size)
 call ZTRSM('Right','Upper','Normal','Normal',my_nt,exc_nst,cone,overlap,exc_nst,phi_block,my_nt)
 ABI_DEALLOCATE(overlap)

end subroutine exc_cholesky_ortho
!!***

!----------------------------------------------------------------------

!!****f* exc_iterative_diago/convergence_degree
!! NAME
!! convergence_degree
!!
!! FUNCTION
!!  Return the degree of convergence from the input residual.
!!
!! INPUTS
!!  resid=Residual.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function convergence_degree(resid)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'convergence_degree'
!End of the abilint section

 implicit none

!Arguments
 integer :: convergence_degree
 real(dp),intent(in) :: resid

!************************************************************************

 if (resid<tolwfr_) then
   convergence_degree = STRICT
 else
   convergence_degree = WORST
   if (resid<tolwfr_*10**5) convergence_degree = MEDIUM
 end if

end function convergence_degree
!!***

!----------------------------------------------------------------------

!!****f* exc_iterative_diago/exc_check_phi_block
!! NAME
!! exc_check_phi_block
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      exc_iterative_diago
!!
!! CHILDREN
!!      wrtout,zhegv
!!
!! SOURCE

subroutine exc_check_phi_block(string)

 use m_profiling

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_check_phi_block'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: string

!Local variables ------------------------------
!scalars
 integer :: ii,jj,ierr
 real(dp) :: err,rdum
!arrays
 complex(dpc),allocatable :: lbuff(:,:)

!************************************************************************

#if 0
 ABI_ALLOCATE(lbuff,(hexc_size,exc_nst))
 err = -one
 do irank=1,nproc-1
   call xexch_mpi(phi_block,hexc_size*exc_nst,irank,lbuff,master,comm,ierr)
   if (my_rank==master) then
     lbuff = lbuff-phi_block
     err = MAX(err,MAXVAL(MAXVAL(ABS(lbuff),DIM=1)))
   end if
   call xbarrier_mpi(comm)
 end do
 ABI_DEALLOCATE(lbuff)
#else

 ABI_ALLOCATE(lbuff,(exc_nst,exc_nst))
 lbuff=czero
 do jj=1,exc_nst
   do ii=1,jj
     lbuff(ii,jj) = DOT_PRODUCT( phi_block(my_t1:my_t2,ii), phi_block(my_t1:my_t2,jj) )
   end do
 end do
 call xsum_mpi(lbuff,comm,ierr)

 err = -one
 do jj=1,exc_nst
   do ii=1,jj
     if (ii==jj) then
       rdum =  ABS(lbuff(ii,jj)-one)
     else
       rdum =  ABS(lbuff(ii,jj))
     end if
     err = MAX(err,rdum)
   end do
 end do
 ABI_DEALLOCATE(lbuff)
#endif

 if (my_rank==master) then
   write(std_out,*)" After ",TRIM(string),", MAX inconsistency error in phi_block= ",err
 end if

 !write(std_out,*)"master casts its own data"
 !call xcast_mpi(phi_block,master,comm,ierr)

end subroutine exc_check_phi_block
!!***

!----------------------------------------------------------------------

!!****f* exc_iterative_diago/rmm_diis_for
!! NAME
!! rmm_diis_for
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout,zhegv
!!
!! SOURCE

subroutine rmm_diis_for(state)

 use m_profiling

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rmm_diis_for'
 use interfaces_14_hidewrite
 use interfaces_71_bse
!End of the abilint section

 implicit none

!Arguments
 integer,intent(in) :: state

!Local variables
!scalars
 integer,parameter :: DIIS_DIM=10
 integer :: idiis,lwork,info,ii
 real(dp) :: etrial_old,etrial,norm,den
 real(dp) :: fac,poly,xx,lambda
 complex(dpc) :: kprc
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: diis_ene(:),rwork(:)
 complex(dpc),allocatable :: work(:)
 complex(dpc),allocatable :: hphi(:)
 complex(dpc),allocatable,target :: phi_diis(:,:),res_diis(:,:)
 complex(dpc),pointer :: res(:),phi(:),hexc_diagonal(:)
 complex(dpc),allocatable :: diis_mat1(:,:),diis_mat2(:,:)

!************************************************************************

 write(msg,'(a,i0)')" Entering rmm-diis for band: ",state
 call wrtout(std_out,msg,"COLL")

 ABI_ALLOCATE(hphi,(hexc_size))

 ABI_ALLOCATE(phi_diis,(hexc_size,DIIS_DIM))
 ABI_ALLOCATE(res_diis,(hexc_size,DIIS_DIM))

 ! phi_block is assumed to contains a good set of orthonormal vectors to be used as starting points.
 ! Switch to RMM-DIIS.

 phi_diis(:,1)   = phi_block(:,state) ! Exctract the state to be optimized
 phi => phi_diis(:,1)
 res => res_diis(:,1)

 hphi = MATMUL(hexc, phi)
 norm = SQRT( DOT_PRODUCT(phi, phi) )
 etrial = DOT_PRODUCT(phi, hphi) / norm
 res = hphi - etrial*phi
 etrial_old = etrial

 ! TODO find optimal value for lambda.
 lambda = 0.5

 do idiis=2,DIIS_DIM

   phi => phi_diis(:,idiis)
   res => res_diis(:,idiis)

   ! Update trial vector.
   phi = phi_diis(:,idiis-1) + lambda*res_diis(:,idiis-1)

   hphi = MATMUL(hexc, phi)
   norm = SQRT( DOT_PRODUCT(phi, phi) )
   etrial = DOT_PRODUCT(phi, hphi) / norm
   res = hphi - etrial*phi

   resid(state) =  DOT_PRODUCT(res,res)
   convergence_of(state) = convergence_degree(resid(state))

   ! Check that etrial is decreasing on succeeding lines:
   if (idiis>1 .and. (etrial > etrial_old+tol12)) then
     write(msg,'(a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&      ' DIIS: New trial exc_energy at idiis',idiis,' = ',etrial,ch10,&
&      ' is higher than former:',etrial_old,ch10
     MSG_WARNING(msg)
   end if
   etrial_old = etrial

   ! If residual sufficiently small, stop line minimization.
   if (convergence_of(state)==STRICT) then
     if (prtvol>=10) then
       write(msg,'(a,i4,a,i2,a,es12.4)')&
&        ' rmm-diis: band ',state,' converged after ',idiis,' RMM-DIIS iterations : resid= ',resid(state)
       call wrtout(std_out,msg,'PERS')
     end if
     EXIT !line
   end if

   ! preconditioning
   den = DOT_PRODUCT(res, hexc_diagonal(:)*res )
   !call xsum_mpi(den,comm,ierr)

   do ii=1,hexc_size
     xx = hexc_diagonal(ii)/den ! Teter polynomial ratio, modified according to Kresse, Furthmuller, PRB 54, 11169 (1996)
     poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
     fac=poly/(poly+16._dp*xx**4)
     kprc = fac*four/(three*den)
     res(ii) = kprc * res(ii)
   end do

   ! Direct inversion in the iterative subspace.
   ABI_ALLOCATE(diis_mat1,(idiis,idiis))
   ABI_ALLOCATE(diis_mat2,(idiis,idiis))

   do jj=1,idiis
     do ii=jj,idiis
       diis_mat1(ii,jj) = DOT_PRODUCT( res_diis(:,ii), res_diis(:,jj) )
       diis_mat2(ii,jj) = DOT_PRODUCT( phi_diis(:,ii), phi_diis(:,jj) )
     end do
   end do

   ! HereI can use zhegvx that however is not shipped with abinit.
   ! On exit, if JOBZ = 'V', then if INFO = 0, A contains the
   ! matrix Z of eigenvectors.  The eigenvectors are normalizedas follows:
   ! if ITYPE = 1 or 2, Z**H*B*Z = I;
   ABI_ALLOCATE(diis_ene,(idiis))

   lwork = MAX(1,2*idiis-1)
   ABI_ALLOCATE(work,(lwork))
   ABI_ALLOCATE(rwork,(MAX(1,3*idiis-2)))

   call ZHEGV(1,"Vectors","Upper",idiis,diis_mat1,idiis,diis_mat2,idiis,diis_ene,work,lwork,rwork,info)

   if (info /=0) then
     write(msg,'(a,i4)')" ZHEGV returned info :",info
     MSG_ERROR(msg)
   end if

   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(rwork)
   ABI_DEALLOCATE(diis_ene)
   ABI_DEALLOCATE(diis_mat2)

   ! Linear combination in the iterative subspace.
   !call ZGEMM('N','N',hexc_size,exc_nst,exc_nst,cone,phi_block,hexc_size,sub_ham,exc_nst,czero,phi_block,hexc_size)
   hphi = MATMUL( phi_diis(:,1:idiis), diis_mat1(:,1) )
   phi_diis(:,idiis) = hphi

   ABI_DEALLOCATE(diis_mat1)

   !converged = .FALSE.
   !if (converged) ! Save result in phi_block(:,state)
   ! phi_block(:,state) = phi_diis(:,idiis); EXIT
   !end if
 end do

 ABI_DEALLOCATE(hphi)
 ABI_DEALLOCATE(phi_diis)
 ABI_DEALLOCATE(res_diis)

 ! Now phi_block contains the converged eigenvectors.
 ! To be orthogonalized when all states have been optimized.
 write(std_out,*)" RMM-DISS resid(state): ",state,resid(state)

end subroutine rmm_diis_for

!----------------------------------------------------------------------

end subroutine exc_iterative_diago
!!***
