!!****f* ABINIT/haydock_psherm
!! NAME
!! haydock_psherm
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file and construct the Lanczos set of vectors 
!!  by iterative matrix-vector multiplications.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!! hize
!! my_t1,my_t2
!! hreso(hsize,my_t1:my_t2)
!! hcoup(hsize,my_t1:my_t2)
!! nkets
!! kets(hsize,nkets)
!! comm=MPI communicator.
!!
!! OUTPUT
!!  green(BSp%nomega)=The imaginary part of the macroscopic dielectric function.
!!
!! PARENTS
!!      haydock
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine haydock_psherm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hreso,hcoup,nkets,kets,green,comm)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use defs_datatypes
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools,       only : get_unit, file_exist, delete_file, flush_unit
 use m_numeric_tools,  only : continued_fract
 use m_blas,           only : xdotc, xgemv
 use defs_abitypes,    only : Hdr_type
 use m_crystal,        only : crystal_structure

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_psherm'
 use interfaces_14_hidewrite
 use interfaces_71_bse, except_this_one => haydock_psherm
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize,my_t1,my_t2,nkets,comm
 type(crystal_structure),intent(in) :: Cryst
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
!arrays
 complex(dp),intent(out) :: green(BSp%nomega,BSp%nq)
 complex(dpc),intent(in) :: hreso(hsize,my_t1:my_t2) 
 complex(dpc),intent(in) :: kets(hsize,nkets)
 complex(dpc),intent(in) :: hcoup(hsize,my_t1:my_t2)

!Local variables ------------------------------
!scalars
 integer :: inn,it,out_unt,ios,nproc,my_rank,master,ierr
 integer :: niter_file,niter_max,niter_done,nsppol,iq,my_nt,term_type
 real(dp) :: ket0_hbar_norm,nfact
 logical :: can_restart,is_converged
 complex(dpc) :: factor
 character(len=fnlen),parameter :: tag_file="_HAYDC_SAVE"
 character(len=500) :: msg
 character(len=fnlen) :: restart_file,out_file
!arrays
 real(dp),pointer :: bb_file(:)
 real(dp),allocatable :: bb(:)
 complex(dpc),allocatable :: aa(:),cc(:),phi_np1(:),phi_n(:),phi_nm1(:),cbuff(:)
 complex(dpc),pointer :: aa_file(:),phi_n_file(:),phi_np1_file(:),cc_file(:)
 complex(dpc),allocatable :: ket0(:)
 logical :: check(2)

!************************************************************************

 MSG_WARNING("Haydock + coupling is still under development")

 nproc  = xcomm_size(comm)
 my_rank= xcomm_rank(comm)
 master = 0 
 nsppol = Hdr_bse%nsppol

 my_nt = my_t2-my_t1+1
 ABI_CHECK(my_nt>0,"One of the processors has zero columns")

 ! Multiplicative factor (k-point sampling and unit cell volume)  
 ! TODO be careful with the spin here
 ! TODO four_pi comes from the coulomb term 1/|q| is already included in the 
 ! oscillators hence the present approach wont work if a cutoff interaction is used.
 nfact = four_pi/(Cryst%ucvol*BSp%nkbz)
 if (nsppol==1) nfact=two*nfact 

 write(msg,'(a,i0)')' Haydock algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")
 !
 ! Check for presence of the restart file.
 can_restart=.FALSE.

 if ( BS_files%in_haydock_basename /= BSE_NOFILE) then
   restart_file = TRIM(BS_files%in_haydock_basename)//TRIM(tag_file)
   if (file_exist(restart_file) ) then
     can_restart=.TRUE.
     msg = " Restarting Haydock calculation from file: "//TRIM(restart_file)
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")
     MSG_ERROR("Restart is not tested")
   else 
     can_restart=.FALSE.
     call wrtout(ab_out," WARNING: cannot find restart file: "//TRIM(restart_file),"COLL")
   end if
 end if
 !
 ! Open the file and writes basic dimensions and info.
 if (my_rank==master) then 
   out_unt = get_unit()
   out_file = TRIM(BS_files%out_basename)//TRIM(tag_file)
   open(unit=out_unt,file=out_file,form="unformatted",iostat=ios)
   ABI_CHECK(ios==0," Opening file: "//TRIM(out_file))
   ! write header TODO: standardize this part.
   write(out_unt)hsize,Bsp%use_coupling,BSE_HAYD_IMEPS,nkets,Bsp%broad
 end if
 !
 ! Select the terminator for the continued fraction.
 term_type=0 !; if (Bsp%hayd_term>0) term_type=2
 write(msg,'(a,i0)')" Using terminator type: ",term_type
 call wrtout(std_out,msg,"COLL")
 !
 ! Calculate green(w) for the different starting kets.
 green=czero
 do iq=1,nkets
   ABI_ALLOCATE(ket0,(my_nt))
   ket0 = kets(my_t1:my_t2,iq)
   !
   niter_file=0
   nullify(aa_file)
   nullify(bb_file)
   nullify(cc_file)
   nullify(phi_np1_file)
   nullify(phi_n_file)

   if (can_restart) then
     call haydock_restart(BSp,restart_file,BSE_HAYD_IMEPS,iq,hsize,&
&      niter_file,aa_file,bb_file,phi_np1_file,phi_n_file,comm)
   end if 
   !
   ABI_ALLOCATE(phi_nm1,(my_nt))
   ABI_ALLOCATE(phi_n,(my_nt))
   ABI_ALLOCATE(phi_np1,(my_nt))
   !
   ! TODO: Note the different convention used for the coefficients
   ! Should use the same convention in the Hermitian case.
   niter_max = niter_file + Bsp%niter
   ABI_ALLOCATE(aa,(niter_max))
   ABI_ALLOCATE(bb,(niter_max+1))
   ABI_ALLOCATE(cc,(niter_max+1))
   aa=czero; bb=czero; cc=czero

   if (niter_file==0) then ! Calculation from scratch.
     phi_n   = ket0
     phi_np1 = MATMUL(hreso,ket0) - MATMUL(hcoup,CONJG(ket0)) 
     ket0_hbar_norm = SQRT(two*DBLE(DOT_PRODUCT(phi_n,phi_np1)))  
     phi_n   = phi_n  /ket0_hbar_norm
     phi_np1 = phi_np1/ket0_hbar_norm
     !ket0    = ket0/ket0_hbar_norm
     cc(1)=zero ! <P|F|P>
     !cc(1) =  DOT_PRODUCT(ket0,phi_np1)
     write(std_out,*)" cc(1), ket0_hbar_norm =",cc(1),ket0_hbar_norm  

     phi_nm1 = czero
     niter_done=0  ! TODO Be careful here

   else ! Use the previously calculates a and b.
     niter_done=niter_file
     MSG_ERROR("Restart not coded")
     !aa(1:niter_done) = aa_file
     !bb(1:niter_done) = bb_file
     !phi_np1=phi_np1_file(my_t1:my_t2)   ! Select the slice treated by this node.
     !phi_n  =phi_n_file  (my_t1:my_t2)   
   end if

   if (associated(aa_file     ))  then
     ABI_DEALLOCATE(aa_file)
   end if
   if (associated(bb_file     ))  then
     ABI_DEALLOCATE(bb_file)
   end if
   if (associated(cc_file     ))  then
     ABI_DEALLOCATE(cc_file)
   end if
   if (associated(phi_np1_file))  then
     ABI_DEALLOCATE(phi_np1_file)
   end if
   if (associated(phi_n_file  ))  then
     ABI_DEALLOCATE(phi_n_file)
   end if

   ! This factor gives the correct results
   factor = -nfact*ket0_hbar_norm / SQRT(two)

   ! Which quantity should be checked for convergence?
   check = (/.TRUE.,.TRUE./) 
   if (ABS(Bsp%haydock_tol(2)-one)<tol6) check = (/.TRUE. ,.FALSE./) 
   if (ABS(Bsp%haydock_tol(2)-two)<tol6) check = (/.FALSE.,.TRUE./) 

   call haydock_psherm_optalgo(niter_done,niter_max,BSp%nomega,BSp%omega,BSp%haydock_tol(1),check,hsize,&
&    my_t1,my_t2,hreso,hcoup,factor,term_type,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,green(:,iq),inn,is_converged,comm)
   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed 
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !    |n+1>
   !
   if (my_rank==master) then ! Open the file and writes basic dimensions and info.
     write(out_unt)Bsp%q(:,iq)
     write(out_unt)MIN(inn,niter_max)  ! NB: if the previous loop completed inn=niter_max+1
     do it=1,MIN(inn,niter_max)        !     if we exited then inn is not incremented by one.
       write(out_unt)it,aa(it),bb(it)
     end do
   end if
   !
   ! cbuff is used as workspace to gather |n-1>, |n> and |n+1>.
   ABI_ALLOCATE(cbuff,(hsize))
   cbuff=czero; cbuff(my_t1:my_t2) = phi_nm1
   call xsum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n-1>

   cbuff=czero; cbuff(my_t1:my_t2) = phi_n
   call xsum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n>

   cbuff=czero; cbuff(my_t1:my_t2) = phi_np1
   call xsum_master(cbuff,master,comm,ierr)
   if (my_rank==master) write(out_unt) cbuff ! |n+1>

   ABI_DEALLOCATE(phi_nm1)
   ABI_DEALLOCATE(phi_n)
   ABI_DEALLOCATE(phi_np1)
   ABI_DEALLOCATE(cbuff)
   ABI_DEALLOCATE(aa)
   ABI_DEALLOCATE(bb)
   ABI_DEALLOCATE(cc)
   ABI_DEALLOCATE(ket0)
 end do ! iq

 if (my_rank==master) close(out_unt)

 call xbarrier_mpi(comm)

end subroutine haydock_psherm
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/haydock_psherm_optalgo
!! NAME
!! haydock_psherm_optalgo
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  niter_done=Number of iterations already performed (0 if the run starts from scratch).
!!  niter_tot=Max number of iterations. Always > niter_done
!!  nomega=Number of Frequency points for the evaluation of the matrix element.
!!  omega(nomega)=Frequency set (imaginary part is already included).
!!  tol_iter=Tollerance used to stop the the algorithm.
!!  check(2)=Logical flags to specify where both the real and the imaginary part of the 
!!    matrix elements of the Green functions have to be checked for convergence. 
!!  hsize=Size of the blocks.
!!  my_t1,my_t2=Indeces of the first and last column stored treated by this done.
!!  term_type=0 if no terminator is used, 1 otherwise.
!!  hreso(hsize,my_t1:my_t2)=The columns of the resonant block.
!!  hcoup(hsize,my_t1:my_t2)=The columns of the coupling block.
!!  factor
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  green(nomega)=Output matrix elements.
!!  inn=Last iteration performed.
!!  is_converged=.TRUE. of the algorithm converged.
!!
!! SIDE EFFECTS
!!  phi_nm1(my_t2-my_t1+1), phi_n(my_t2-my_t1+1)
!!    input: vectors used to initialize the iteration
!!    output: the vectors obtained in the last iteration 
!!  aa(niter_tot) and bb(niter_tot+1)
!!    if niter_done>0: aa(1:niter_done), bb(1:niter_done) store the coefficients of the previous run.
!!    when the routine returns aa(1:inn) and bb(1:inn) contain the matrix elements of the tridiagonal form.
!!  cc(niter_tot+1)
!!
!! PARENTS
!!      haydock_psherm
!!
!! CHILDREN
!!      continued_fract,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine haydock_psherm_optalgo(niter_done,niter_tot,nomega,omega,tol_iter,check,hsize,my_t1,my_t2,hreso,hcoup,&
&  factor,term_type,aa,bb,cc,ket0,ket0_hbar_norm,phi_nm1,phi_n,phi_np1,green,inn,is_converged,comm)

 use m_profiling

 use defs_basis
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_numeric_tools,  only : continued_fract, print_arr
 use m_blas,           only : xdotc, xgemv, xgemm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_psherm_optalgo'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter_tot,niter_done,nomega,comm,hsize,my_t1,my_t2,term_type
 integer,intent(out) :: inn
 logical,intent(out) :: is_converged
 real(dp),intent(in) :: tol_iter,ket0_hbar_norm
 complex(dpc),intent(in) :: factor
!arrays
 real(dp),intent(inout) :: bb(niter_tot+1)
 complex(dpc),intent(out) :: green(nomega)
 complex(dpc),intent(in) :: omega(nomega) 
 complex(dpc),intent(inout) :: aa(niter_tot),cc(niter_tot+1)
 complex(dpc),intent(in) :: hreso(hsize,my_t1:my_t2)
 complex(dpc),intent(in) :: hcoup(hsize,my_t1:my_t2)
 complex(dpc),intent(in) :: ket0(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_n  (my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_np1(my_t2-my_t1+1)
 logical,intent(in) :: check(2)

!Local variables ------------------------------
!scalars
 integer :: my_nt,niter_min,nconv,parity,ii,jj,tdim,istat !ierr
 integer :: row_max,col_max,nlev
 character(len=500) :: msg
 real(dp) :: max_err,mean_err,mean_err2,std_dev,err
 logical :: keep_vectors=.TRUE.
!arrays
 real(dp) :: abs_err(nomega,2) !,ww_err(nomega,2)
 complex(dpc) :: gn0(nomega,niter_tot)
 complex(dpc),allocatable :: oldg(:),newg(:) 
 complex(dpc),allocatable :: hphi_n(:),save_phi(:,:)
 complex(dpc),allocatable ::  alpha(:,:),beta(:,:),ovlp(:,:)
 complex(dpc),allocatable :: phi_test(:),phi_test2(:),g00(:)
 logical :: test(2)

!************************************************************************

 ABI_UNUSED(ket0_hbar_norm)

 my_nt = my_t2-my_t1+1

 ABI_ALLOCATE(oldg,(nomega))
 oldg=czero 
 ABI_ALLOCATE(newg,(nomega))
 newg=czero 
 ABI_ALLOCATE(g00,(nomega))
 g00=czero 
 nconv=0

 keep_vectors = (keep_vectors.and.xcomm_size(comm)==1)
 if (keep_vectors) then 
   ABI_ALLOCATE(save_phi,(my_t2-my_t1+1,niter_tot))
   istat = ABI_ALLOC_STAT
   
   ABI_CHECK(istat==0,"out of memory in save_phi")        
   save_phi=czero
 end if

 ABI_ALLOCATE(hphi_n,(hsize))

 do inn=niter_done+1,niter_tot
   !
   ! a(n) = <Vn+1|F|Vn+1> = <Vn|HFH|Vn>) = 0 by symmetry.
   aa(inn)=zero

   ! |n+1> = |n+1> - a(n)|Vn> - a(n)|n-1>
   phi_np1 = phi_np1 - bb(inn)*phi_nm1
   !
   ! |n-1> = |n> 
   ! |n>   = |n+1> 
   phi_nm1 = phi_n
   phi_n   = phi_np1
   !
   !|n+1> = H |n> using resonant eh components.
   parity = (-1)**(inn+1)
   phi_np1 = MATMUL(hreso,phi_n) + parity * MATMUL(hcoup,CONJG(phi_n))
   !call xsum_mpi(hphi_np1,comm,ierr)
   !
   ! B(n+1)= <n|F|n+1>^(1/2) = <n|FH|n>^(1/2))= (2*Re(<n|V+1>))^(1/2) 
   ! by symmetry, where the dot_product is done in the resonant eh sub-space. 
   !
   bb(inn+1)=SQRT(two*DBLE(DOT_PRODUCT(phi_n,phi_np1)))
   !bb(inn+1)=two*DBLE(DOT_PRODUCT(phi_n,phi_np1))
   !call xsum_mpi(bb(inn+1),comm,ierr)
   !bb(inn+1)=SQRT(bb(inn+1)
   !
   !|n+1> =|n+1>/B(n+1)
   phi_n   = phi_n  /bb(inn+1)
   phi_np1 = phi_np1/bb(inn+1)

   if (keep_vectors) save_phi(:,inn) = phi_n

   parity = (-1)**(inn+1) 
   !if (parity==-1) then 
   !  cc(inn+1)=czero
   !else 
     cc(inn+1)=DOT_PRODUCT(ket0,phi_n) + parity * DOT_PRODUCT(phi_n,ket0)
   !end if
   !call xsum_mpi(cc(inn+1),comm,ierr)

   write(msg,'(a,i0,a,3es12.4)')' Iteration number ',inn,', b_i RE(c_i+1) IM(c_i+1) ',bb(inn),REAL(cc(inn+1)),AIMAG(cc(inn+1)) 
   call wrtout(std_out,msg,"COLL")

   call continued_fract(inn,term_type,aa,bb(2:),nomega,omega,g00)
   gn0(:,1) = g00

   if (.FALSE.) then
     gn0(:,2) = (one - omega(:)*g00(:))/bb(2)
     do ii=3,inn
       gn0(:,ii) = -(-bb(ii)*gn0(:,ii-2) -omega(:)*gn0(:,ii-1))/bb(ii+1)
     end do
   else 
     do ii=2,inn
       nlev = inn-ii
       call continued_fract(nlev,term_type,aa,bb(ii+1:),nomega,omega,g00)
       gn0(:,ii) = +bb(ii+1) * g00 * gn0(:,ii-1)
     end do
   end if

   newg=czero
   do ii=1,inn
     newg(:) = newg + cc(ii)* gn0(:,ii)
   end do
   newg = factor*newg
   !
   ! Avoid spurious convergence.
   niter_min=4; if (niter_done>1) niter_min=niter_done+1
   if (inn>niter_min) then
     test=.TRUE.
     abs_err(:,1) = ABS(DBLE (newg-oldg))
     abs_err(:,2) = ABS(AIMAG(newg-oldg))
     !
     if (tol_iter>zero) then 
       ! Test on the L1 norm.
       if (check(1)) test(1) = SUM(abs_err(:,1)) < tol_iter*SUM(ABS(DBLE (newg)))
       if (check(2)) test(2) = SUM(abs_err(:,2)) < tol_iter*SUM(ABS(AIMAG(newg)))
     else
       ! Stringent test for each point.
       if (check(1)) test(1) = ALL( abs_err(:,1) < -tol_iter*ABS(DBLE (newg))) 
       if (check(2)) test(2) = ALL( abs_err(:,2) < -tol_iter*ABS(AIMAG(newg)))
     end if
     !
     if (ALL(test)) then 
       nconv = nconv+1
     else 
       nconv = 0
     end if
     if (nconv==2) then 
       write(msg,'(a,es10.2,a,i0,a)')&
&        " >>> Haydock algorithm converged twice within haydock_tol= ",tol_iter," after ",inn," iterations." 
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       EXIT
     end if
   end if
   !
   oldg = newg
 end do ! inn

 green = newg
 if (nconv/=2) then
   write(msg,'(a,es10.2,a,i0,a)')&
&    " WARNING: Haydock algorithm did not converge within ",tol_iter," after ",niter_tot," iterations."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 is_converged = (nconv==2)

 ABI_DEALLOCATE(oldg)
 ABI_DEALLOCATE(newg)
 ABI_DEALLOCATE(g00)
 ABI_DEALLOCATE(hphi_n)

 if (keep_vectors) then
   tdim = MIN(inn,niter_tot)
   ABI_ALLOCATE(ovlp,(tdim,tdim))

   ABI_ALLOCATE(phi_test,(hsize))
   ABI_ALLOCATE(phi_test2,(hsize))

   max_err=smallest_real; mean_err=zero; mean_err2=zero; row_max=-1
   do ii=1,tdim
     parity = (-1)**(ii+1) 
     phi_test  = save_phi(:,ii)
     phi_test2 = MATMUL(hreso,phi_test) + parity * MATMUL(hcoup,CONJG(phi_test))
     ovlp(ii,ii) = DOT_PRODUCT(phi_test,phi_test2) + DOT_PRODUCT(phi_test2,phi_test) 
     err = ABS(ovlp(ii,ii)-cone)
     mean_err  = mean_err + err
     mean_err2 = mean_err2 + err**2
     if (err > max_err) then
       max_err = err 
       row_max = ii
     end if
   end do
   mean_err = mean_err/tdim
   std_dev = mean_err2/tdim -mean_err**2
   write(std_out,'(a,i0,1x,3es14.6)')&
&   " Error in normalization (ii, max_err,mean,std_dev): ",row_max,max_err,mean_err,std_dev

   ABI_DEALLOCATE(phi_test)
   ABI_DEALLOCATE(phi_test2)
                                             
   ABI_ALLOCATE(alpha,(hsize,tdim))
   alpha = MATMUL(hreso,save_phi(:,1:tdim))

   do ii=1,tdim
     parity = (-1)**(ii+1) 
     alpha(:,ii) =  alpha(:,ii) + parity*MATMUL(hcoup,CONJG(save_phi(:,ii))) 
   end do

   ovlp = MATMUL(TRANSPOSE(CONJG(save_phi(:,1:tdim))),alpha)

   ABI_ALLOCATE(beta,(hsize,tdim))
   do ii=1,tdim
     parity = (-1)**(ii+1) 
     beta(:,ii)  =  parity*save_phi(:,ii)
     alpha(:,ii) = -parity*alpha(:,ii)
   end do

   ovlp = ovlp - MATMUL(TRANSPOSE(CONJG(beta)),alpha)

   max_err=smallest_real; row_max=-1; col_max=-1
   mean_err=zero; mean_err2=zero
   do jj=1,tdim
     do ii=1,jj
       err = ABS(ovlp(ii,jj))
       if (ii==jj) err = ABS(err - one)
       mean_err  = mean_err + err
       mean_err2 = mean_err2 + err**2
       if (err > max_err) then
         max_err = err
         row_max=ii
         col_max=jj
       end if
     end do
   end do

   mean_err = mean_err/(tdim*(tdim+1)/2)
   std_dev = mean_err2/(tdim*(tdim+1)/2) - mean_err**2
   write(std_out,'(a,2(i0,1x),3es14.6)')&
&     " Error in Hbar-ortho (i,j), max_err, mean, std_dev ",row_max,col_max,max_err,mean_err,std_dev
   !call print_arr(ovlp,max_r=185,max_c=10,unit=std_out)

   ABI_DEALLOCATE(alpha)
   ABI_DEALLOCATE(beta)
   ABI_DEALLOCATE(ovlp)
   ABI_DEALLOCATE(save_phi)
 end if

end subroutine haydock_psherm_optalgo
!!***

!----------------------------------------------------------------------
