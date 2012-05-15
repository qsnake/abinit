!{\src2tex{textfont=tt}}
!!****f* ABINIT/exc_haydock_driver
!! NAME
!! exc_haydock_driver
!!
!! FUNCTION
!!  Calculate the imaginary part of the macroscopic dielectric function via the Haydock recursive method.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi, Y. Gillet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! BSp<type(excparam)=The parameter for the Bethe-Salpeter run.
!!  inclvkb=If different from zero, the commutator [Vnl,r] is included in the calculation of 
!!    the matrix element of the velocity operator. Meaningless for PAW.
!! Kmesh<type(bz_mesh_type)>=The list of k-points in the BZ, IBZ and symmetry tables.
!! Cryst<type(crystal_structure)>=Info on the crystalline structure.
!! KS_BSt=The KS energies.
!! QP_BSt=The QP energies.
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials.
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data.
!! Hur(Cryst%natom*usepaw)<type(HUr_commutator)>=Only for PAW and LDA+U, quantities used to evaluate the commutator [H_u,r].
!! Wfd<wfs_descriptor>=Handler for the wavefunctions.
!!   %nsppol=Number of independent spin polarizations.
!!   %nspinor=Number of spinorial components.
!!   %usepaw=1 for PAW, 0 otherwise.
!!   %comm=MPI communicator.
!!
!! OUTPUT
!!  The imaginary part of the macroscopic dielectric function is written on the external file _EXC_MDF
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      matrginv,zgesv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine exc_haydock_driver(BSp,BS_files,Cryst,Kmesh,Hdr_bse,KS_BSt,QP_Bst,Wfd,Psps,Pawtab,Hur)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_bs_defs
 use m_xmpi
 use m_errors

 use defs_abitypes,       only : Hdr_type
 use m_blas,              only : xdotc
 use m_numeric_tools,     only : print_arr, symmetrize, hermitianize
 use m_crystal,           only : crystal_structure 
 use m_bz_mesh,           only : bz_mesh_type
 use m_commutator_vkbr,   only : kb_potential
 use m_paw_commutator,    only : HUr_commutator
 use m_wfs,               only : wfs_descriptor
 use m_bse_io,            only : exc_read_rcblock

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_haydock_driver'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_71_bse, except_this_one => exc_haydock_driver
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(bz_mesh_type),intent(in) :: Kmesh
 type(crystal_structure),intent(in) :: Cryst
 type(Hdr_type),intent(in) :: Hdr_bse
 type(wfs_descriptor),intent(inout) :: Wfd
 type(pseudopotential_type),intent(in) :: Psps
 type(Bandstructure_type),intent(in) :: KS_BSt,QP_Bst
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(HUr_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: io,my_rank,master,iq,istat,it,ierr
 integer :: hsize,comm,my_t1,my_t2,nsppol,nkets,nproc
 integer :: spin,spad,ik_bz,iv,ic,trans_idx,minb,maxb
 integer :: max_r,max_c
 real(dp) :: omegaev,rand_phi !,norm
 complex(dpc) :: ks_avg,gw_avg,exc_avg
 logical :: is_resonant,use_mpio,diago_is_real,prtdos
 character(len=500) :: msg
 character(len=fnlen) :: hreso_fname,hcoup_fname
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: dos(:),dos_gw(:),dos_ks(:)
 complex(dpc),allocatable :: green(:,:),hreso(:,:),hcoup(:,:),test(:,:)
 complex(dpc),allocatable :: opt_cvk(:,:,:,:,:),kets(:,:)
 complex(dpc),allocatable :: eps_rpanlf(:,:),eps_gwnlf(:,:)
 complex(dpc),allocatable :: tensor_cart(:,:),tensor_cart_rpanlf(:,:),tensor_cart_gwnlf(:,:)
 complex(dpc),allocatable :: tensor_red(:,:),tensor_red_rpanlf(:,:),tensor_red_gwnlf(:,:)

!************************************************************************

 call timab(680,1,tsec) ! exc_haydock_driver
 call timab(681,1,tsec) ! exc_haydock_driver(read)
 if (BSp%have_complex_ene) then
   MSG_ERROR("Complex energies are not supported yet")
 end if

 my_rank = Wfd%my_rank
 master  = Wfd%master
 comm    = Wfd%comm
 nsppol  = Wfd%nsppol
 nproc   = Wfd%nproc

 use_mpio=.FALSE.
#ifdef HAVE_MPI_IO
 use_mpio = (nproc > 1)
 !use_mpio = .TRUE. 
#endif
 use_mpio=.FALSE.
 !use_mpio = .TRUE. 

 ! Hsize refers to the size of the individual blocks (resonant and coupling). 
 ! Thanks to the symmetry property of the starting vector, the Haydock method 
 ! can be reformulated in terms of matrix-vector multiplication involving the 
 ! blocks thus avoiding to allocation of the full matrix ( R   C )
 !                                                        -C* -R*)
 hsize=SUM(BSp%nreh)
 !
 ! Divide the columns of the Hamiltonian among the nodes.
 call xmpi_split_work(hsize,comm,my_t1,my_t2,msg,ierr)
 if (ierr/=0) then
   MSG_WARNING(msg)
 end if

 ABI_CHECK(my_t2-my_t1+1>0,"found processor with 0 rows")
                                                             
 ABI_ALLOCATE(hreso,(hsize,my_t1:my_t2))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in hreso")
 !
 ! Read the resonant block from file.
 if (BS_files%in_hreso /= BSE_NOFILE) then
   hreso_fname = BS_files%in_hreso
 else 
   hreso_fname = BS_files%out_hreso
 end if

 is_resonant=.TRUE.; diago_is_real=(.not.BSp%have_complex_ene)
 call exc_read_rcblock(hreso_fname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,my_t1,my_t2,hreso,use_mpio,comm)

 !call hermitianize(hreso,"All")

!BEGIN DEBUG
 if (use_mpio) then
   MSG_WARNING("Testing MPI-IO routines")
   ABI_ALLOCATE(test,(hsize,my_t1:my_t2))
   istat = ABI_ALLOC_STAT
   ABI_CHECK(istat==0,"out of memory in hreso")
   diago_is_real=(.not.BSp%have_complex_ene)
   call exc_read_rcblock(hreso_fname,Bsp,is_resonant,diago_is_real,nsppol,Bsp%nreh,hsize,my_t1,my_t2,test,.FALSE.,comm)
   test = test-hreso
   write(std_out,*)"DEBUG: Diff MPI-IO - Fortran ",MAXVAL(ABS(test))
   max_r=20; max_c=10
   write(std_out,*)" **** Testing resonant block **** "
   call print_arr(test,max_r=max_r,max_c=max_c,unit=std_out)
   if (nsppol==2) then
     write(std_out,*)" **** D down down ****"
     call print_arr(test(hsize/2+1:,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
     write(std_out,*)" **** V up down ****"
     call print_arr(test(1:hsize/2,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
     write(std_out,*)" **** V down up ****"
     call print_arr(test(hsize/2+1:,1:hsize/2),max_r=max_r,max_c=max_c,unit=std_out)
   end if
   ABI_DEALLOCATE(test)
 end if
!END DEBUG
 !
 ! Read coupling block.
 if (BSp%use_coupling>0) then 
   if (BS_files%in_hcoup /= BSE_NOFILE) then
     hcoup_fname = BS_files%in_hcoup
   else 
     hcoup_fname = BS_files%out_hcoup
   end if

   ABI_ALLOCATE(hcoup,(hsize,my_t1:my_t2))
   istat = ABI_ALLOC_STAT
   is_resonant=.FALSE.; diago_is_real=.FALSE.
   call exc_read_rcblock(hcoup_fname,Bsp,is_resonant,diago_is_real,nsppol,BSp%nreh,hsize,my_t1,my_t2,hcoup,use_mpio,comm)

   !call symmetrize(hcoup,"ALL")

   if (use_mpio) then
     MSG_WARNING("Testing MPI-IO routines")
     ABI_ALLOCATE(test,(hsize,my_t1:my_t2))
     istat = ABI_ALLOC_STAT
     ABI_CHECK(istat==0,"out of memory in hreso")
     diago_is_real=.FALSE.
     call exc_read_rcblock(hcoup_fname,Bsp,is_resonant,diago_is_real,nsppol,Bsp%nreh,hsize,my_t1,my_t2,test,.FALSE.,comm)
     test = test-hcoup
     write(std_out,*)"DEBUG: Diff MPI-IO - Fortran ",MAXVAL(ABS(test))
     max_r=20; max_c=10
     write(std_out,*)" **** Testing coupling block **** "
     call print_arr(test,max_r=max_r,max_c=max_c,unit=std_out)
     if (nsppol==2) then
       write(std_out,*)" **** D down down ****"
       call print_arr(test(hsize/2+1:,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
       write(std_out,*)" **** V up down ****"
       call print_arr(test(1:hsize/2,hsize/2+1:),max_r=max_r,max_c=max_c,unit=std_out)
       write(std_out,*)" **** V down up ****"
       call print_arr(test(hsize/2+1:,1:hsize/2),max_r=max_r,max_c=max_c,unit=std_out)
     end if
     ABI_DEALLOCATE(test)
   end if
 end if

 call timab(681,2,tsec) ! exc_haydock_driver(read)
 call timab(682,1,tsec) ! exc_haydock_driver(prep)

 !
 ! Prepare the starting vectors for the Lanczos chain.
 nkets=Bsp%nq

 prtdos=.FALSE.
 !prtdos=.TRUE.
 if (prtdos) then
   nkets=nkets+1
   if (Bsp%use_coupling>0) then 
     MSG_ERROR("DOS with coupling not coded")
     nkets=nkets+1
   end if
 end if

 ABI_ALLOCATE(kets,(hsize,nkets))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in kets")
 kets=czero
 !
 ! Prepare the kets for the macroscopic dielectric function.
 minb=Bsp%lomo; maxb=Bsp%nbnds
 ABI_ALLOCATE(opt_cvk,(minb:maxb,minb:maxb,BSp%nkbz,Wfd%nsppol,BSp%nq))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in opt_cvk")

 do iq=1,Bsp%nq
   !
   ! KS_BSt is used here to calculate the commutator.
   call calc_optical_mels(Wfd,Kmesh,KS_BSt,Cryst,Psps,Pawtab,Hur,BSp%inclvkb,minb,maxb,BSp%nkbz,BSp%q(:,iq),opt_cvk(:,:,:,:,iq))
   !
   ! Fill ket0 using the same ordering for the indeces as the one used for the excitonic Hamiltonian.
   ! Note that only the resonant part is used here.
   do spin=1,nsppol
     spad=(spin-1)*BSp%nreh(1)
     do ik_bz=1,BSp%nkbz
       do iv=BSp%lomo,BSp%homo
         do ic=BSp%lumo,BSp%nbnds
           trans_idx = BSp%vcks2t(iv,ic,ik_bz,spin)
          if (trans_idx>0) kets(trans_idx+spad,iq)=opt_cvk(ic,iv,ik_bz,spin,iq)
         end do
       end do
     end do
   end do
 end do

 call timab(682,2,tsec) ! exc_haydock_driver(prep)
 call timab(683,1,tsec) ! exc_haydock_driver(wo lf    - that is, without local field
 !
 ! =======================================================
 ! === Make EPS RPA and GW without local-field effects ===
 ! =======================================================

 call wrtout(std_out," Calculating RPA NLF and QP NLF epsilon","COLL")

 ABI_ALLOCATE(eps_rpanlf,(BSp%nomega,BSp%nq))
 ABI_ALLOCATE(dos_ks,(BSp%nomega))
 call exc_eps_rpa(BSp%nbnds,BSp%lomo,BSp%homo,Kmesh,KS_BSt,BSp%nq,nsppol,opt_cvk,Cryst%ucvol,BSp%broad,BSp%nomega,BSp%omega,&
&  eps_rpanlf,dos_ks)

 ABI_ALLOCATE(eps_gwnlf ,(BSp%nomega,BSp%nq))
 ABI_ALLOCATE(dos_gw,(BSp%nomega))
 call exc_eps_rpa(BSp%nbnds,BSp%lomo,BSp%homo,Kmesh,QP_BSt,BSp%nq,nsppol,opt_cvk,Cryst%ucvol,Bsp%broad,BSp%nomega,BSp%omega,&
&  eps_gwnlf,dos_gw)

 if (my_rank==master) then ! Only master works.
   !
   ! Master node writes final results on file.
   call exc_write_data(BSp,BS_files,"RPA_NLF_MDF",eps_rpanlf,dos=dos_ks)

   call exc_write_data(BSp,BS_files,"GW_NLF_MDF",eps_gwnlf,dos=dos_gw)

   ! Computing and writing tensor in files
   ! It works only with 6 q-points
   if (BSp%nq .eq. 6) then

     ! RPA_NLF
     ABI_ALLOCATE(tensor_cart_rpanlf,(BSp%nomega,6))
     ABI_ALLOCATE(tensor_red_rpanlf,(BSp%nomega,6))
     ! Computing tensor
     call haydock_mdf_to_tensor(BSp,Cryst,eps_rpanlf,tensor_cart_rpanlf, tensor_red_rpanlf)

     ! Writing tensor
     call exc_write_tensor(BSp,BS_files,"RPA_NLF_TSR_CART",tensor_cart_rpanlf)
     call exc_write_tensor(BSp,BS_files,"RPA_NLF_TSR_RED",tensor_red_rpanlf)

     ABI_DEALLOCATE(tensor_cart_rpanlf)
     ABI_DEALLOCATE(tensor_red_rpanlf)

     ! GW_NLF
     ABI_ALLOCATE(tensor_cart_gwnlf,(BSp%nomega,6))
     ABI_ALLOCATE(tensor_red_gwnlf,(BSp%nomega,6))
     ! Computing tensor
     call haydock_mdf_to_tensor(BSp,Cryst,eps_gwnlf,tensor_cart_gwnlf, tensor_red_gwnlf)

     ! Writing tensor
     call exc_write_tensor(BSp,BS_files,"GW_NLF_TSR_CART",tensor_cart_gwnlf)
     call exc_write_tensor(BSp,BS_files,"GW_NLF_TSR_RED",tensor_red_gwnlf)

     ABI_DEALLOCATE(tensor_cart_gwnlf)
     ABI_DEALLOCATE(tensor_red_gwnlf)
   
   end if
 
   !call wrtout(std_out," Checking Kramers Kronig on Excitonic Macroscopic Epsilon","COLL")
   !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_exc(:,1))

   !call wrtout(std_out," Checking Kramers Kronig on RPA NLF Macroscopic Epsilon","COLL")
   !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_rpanlf(:,1))

   !call wrtout(std_out," Checking Kramers Kronig on GW NLF Macroscopic Epsilon","COLL")
   !call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_gwnlf(:,1))

   !call wrtout(std_out," Checking f-sum rule on Excitonic Macroscopic Epsilon","COLL")

   !if (BSp%exchange_term>0) then 
   !  MSG_COMMENT(' f-sum rule should be checked without LF')
   !end if
   !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_exc(:,1)),drude_plsmf)

   !call wrtout(std_out," Checking f-sum rule on RPA NLF Macroscopic Epsilon","COLL")
   !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_rpanlf(:,1)),drude_plsmf)

   !call wrtout(std_out," Checking f-sum rule on GW NLF Macroscopic Epsilon","COLL")
   !call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_gwnlf(:,1)),drude_plsmf)
 end if ! my_rank==master

 ABI_DEALLOCATE(opt_cvk)
 !call xbarrier_mpi(comm)
 !
 ! The ket for the approximated DOS.
 if (prtdos) then 
   MSG_WARNING("Calculating DOS with Haydock method")
   ABI_CHECK(BSp%use_coupling==0,"DOS with coupling not coded")
   iq = BSp%nq + 1
   if (my_rank==master) then
     !call random_seed()
     do it=1,SUM(Bsp%nreh)
       call RANDOM_NUMBER(rand_phi)
       rand_phi = two_pi*rand_phi
       kets(it,iq) = CMPLX( COS(rand_phi), SIN(rand_phi) )
     end do
     ! Normalize the vector.
     !norm = SQRT( DOT_PRODUCT(kets(:,iq), kets(:,iq)) ) 
     !kets(:,iq) = kets(:,iq)/norm
   end if
   call xcast_mpi(kets(:,iq),master,comm,ierr)
 end if

 call timab(683,2,tsec) ! exc_haydock_driver(wo lf    - that is, without local field
 call timab(684,1,tsec) ! exc_haydock_driver(apply

 ABI_ALLOCATE(green,(BSp%nomega,nkets))

 if (BSp%use_coupling==0) then 
   call haydock_herm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hreso,nkets,kets,green,comm)
 else
   call haydock_psherm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hreso,hcoup,nkets,kets,green,comm)
 end if
 !
 ! Add 1 to have the real part right.
 green = one + green

 ABI_DEALLOCATE(kets)

 call timab(684,2,tsec) ! exc_haydock_driver(apply
 call timab(685,1,tsec) ! exc_haydock_driver(end)

 if (my_rank==master) then ! Master writes the final results.
   !
   if (prtdos) then
     ABI_ALLOCATE(dos,(BSp%nomega))
     dos = -AIMAG(green(:,BSp%nq+1))
     call exc_write_data(BSp,BS_files,"EXC_MDF",green,dos=dos)
     ABI_DEALLOCATE(dos)
   else 
     call exc_write_data(BSp,BS_files,"EXC_MDF",green)
   end if
   !
   ! =========================
   ! === Write out Epsilon ===
   ! =========================

   ! It works only with 6 q-points
   if (BSp%nq .eq. 6) then

     ABI_ALLOCATE(tensor_cart,(BSp%nomega,6))
     ABI_ALLOCATE(tensor_red,(BSp%nomega,6))
     ! Computing tensor
     call haydock_mdf_to_tensor(BSp,Cryst,green,tensor_cart, tensor_red)

     ! Writing tensor
     call exc_write_tensor(BSp,BS_files,"EXC_TSR_CART",tensor_cart)
     call exc_write_tensor(BSp,BS_files,"EXC_TSR_RED",tensor_red)

     ABI_DEALLOCATE(tensor_cart)
     ABI_DEALLOCATE(tensor_red)

   end if
   !
   ! This part will be removed when fldiff will be able to compare two mdf files.
   write(ab_out,*)" "
   write(ab_out,*)"Macroscopic dielectric function:"
   write(ab_out,*)"omega [eV] <KS_RPA_nlf>  <GW_RPA_nlf>  <BSE> "
   do io=1,MIN(BSp%nomega,10)
     omegaev = REAL(BSp%omega(io))*Ha_eV
     ks_avg  = SUM( eps_rpanlf(io,:)) / Bsp%nq
     gw_avg  = SUM( eps_gwnlf (io,:)) / Bsp%nq
     exc_avg = SUM( green     (io,:)) / BSp%nq
     write(ab_out,'(7f9.4)')omegaev,ks_avg,gw_avg,exc_avg
   end do
   write(ab_out,*)" "
 end if 

 ABI_DEALLOCATE(green)

 ABI_DEALLOCATE(eps_rpanlf)
 ABI_DEALLOCATE(eps_gwnlf)
 ABI_DEALLOCATE(dos_ks)
 ABI_DEALLOCATE(dos_gw)

 ABI_DEALLOCATE(hreso)
 if (allocated(hcoup))  then
   ABI_DEALLOCATE(hcoup)
 end if

 call timab(685,2,tsec) ! exc_haydock_driver(end)
 call timab(680,2,tsec) ! exc_haydock_driver

end subroutine exc_haydock_driver
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/haydock_herm
!! NAME
!! haydock_herm
!!
!! FUNCTION
!!  Reads the excitonic Hamiltonian from file and construct the Lanczos set of vectors 
!!  by iterative matrix-vector multiplications.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi, Y. Gillet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! BSp<excparam>=Parameters for the Bethe-Salpeter calculation.
!! BS_files<excparam>=Files associated to the bethe_salpeter code.
!! Cryst<crystal_structure>=Info on the crystalline structure.
!! hize=Size of the excitonic matrix.
!! my_t1,my_t2=First and last columns treated by this node.
!! hmat(hsize,my_t1:my_t2)=Excitonic matrix.
!! nkets=Number of starting vectors for Haydock method.
!! kets(hsize,nkets)=The kets in the eh representation.
!! comm=MPI communicator.
!!
!! OUTPUT
!!  green(BSp%nomega,nkets)=
!!
!! PARENTS
!!      haydock
!!
!! CHILDREN
!!      matrginv,zgesv
!!
!! SOURCE

subroutine haydock_herm(BSp,BS_files,Cryst,Hdr_bse,hsize,my_t1,my_t2,hmat,nkets,kets,green,comm)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use defs_datatypes
 use m_xmpi
 use m_errors

 use m_io_tools,       only : get_unit, file_exist, delete_file, flush_unit
 use m_numeric_tools,  only : continued_fract
 use m_blas,           only : xdotc, xgemv
 use defs_abitypes,    only : Hdr_type
 use m_crystal,        only : crystal_structure
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_herm'
 use interfaces_14_hidewrite
 use interfaces_71_bse, except_this_one => haydock_herm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: hsize,my_t1,my_t2,nkets,comm
 type(crystal_structure),intent(in) :: Cryst
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(Hdr_type),intent(in) :: Hdr_bse
!arrays
 complex(dp),intent(out) :: green(BSp%nomega,nkets)
 complex(dpc),intent(in) :: hmat(hsize,my_t1:my_t2),kets(hsize,nkets)

!Local variables ------------------------------
!scalars
 integer :: inn,it,out_unt,ios,nproc,my_rank,master,ierr
 integer :: niter_file,niter_max,niter_done,nsppol,iq,my_nt,term_type
 integer :: n_all_omegas
 real(dp) :: norm,nfact
 logical :: can_restart,is_converged
 complex(dpc) :: factor
 character(len=500) :: msg
 character(len=fnlen),parameter :: tag_file="_HAYDR_SAVE"
 character(len=fnlen) :: restart_file,out_file 
!arrays
 real(dp),pointer :: bb_file(:)
 real(dp),allocatable :: bb(:)
 complex(dpc),allocatable :: aa(:),phi_nm1(:),phi_n(:),hphi_n(:)
 complex(dpc),pointer :: aa_file(:),phi_n_file(:),phi_nm1_file(:)
 complex(dpc),allocatable :: ket0(:)
 complex(dpc),allocatable :: all_omegas(:)
 complex(dpc),allocatable :: green_temp(:,:)
 logical :: check(2)
 
!************************************************************************

 nproc  = xcomm_size(comm)
 my_rank= xcomm_rank(comm)
 master = 0 
 nsppol = Hdr_bse%nsppol

 my_nt = my_t2-my_t1+1
 ABI_CHECK(my_nt>0,"One of the processors has zero columns")

 write(msg,'(a,i0)')' Haydock algorithm with MAX number of iterations: ',BSp%niter
 call wrtout(std_out,msg,"COLL")
 !
 ! Select the terminator for the continued fraction.
 term_type=0; if (Bsp%hayd_term>0) term_type=1 
 write(msg,'(a,i0)')" Using terminator type: ",term_type
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
   else 
     can_restart=.FALSE.
     call wrtout(ab_out," WARNING: cannot find restart file: "//TRIM(restart_file),"COLL")
   end if
 end if
 !
 ! Open the file and write basic dimensions and info.
 if (my_rank==master) then 
   out_unt = get_unit()
   out_file = TRIM(BS_files%out_basename)//TRIM(tag_file)
   open(unit=out_unt,file=out_file,form="unformatted",iostat=ios)
   ABI_CHECK(ios==0," Opening file: "//TRIM(out_file))
   ! write header TODO: standardize this part.
   write(out_unt)hsize,Bsp%use_coupling,BSE_HAYD_IMEPS,nkets,Bsp%broad
 end if
 !
 ! Calculate green(w) for the different starting points.
 green=czero
 do iq=1,nkets
   ABI_ALLOCATE(ket0,(hsize))
   ket0=kets(:,iq)
   !
   niter_file=0
   nullify(aa_file)
   nullify(bb_file)
   nullify(phi_nm1_file)
   nullify(phi_n_file)

   if (can_restart) then
     call haydock_restart(BSp,restart_file,BSE_HAYD_IMEPS,iq,hsize,&
&      niter_file,aa_file,bb_file,phi_nm1_file,phi_n_file,comm)
   end if 
   !
   ! For n>1, we have:
   !  1) a_n = <n|H|n>
   !  2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
   !  3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
   !
   ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
   !  a_1 = <1|H|1>
   !  b_1 = || H|1> - a_1|1> ||
   !  |2> = [H|1> - a_1|1>]/b_1
   !
   ABI_ALLOCATE(hphi_n,(hsize))
   ABI_ALLOCATE(phi_nm1,(my_nt))
   ABI_ALLOCATE(phi_n,(my_nt))

   niter_max = niter_file + Bsp%niter
   ABI_ALLOCATE(aa,(niter_max))
   ABI_ALLOCATE(bb,(niter_max))
   aa=czero; bb=zero

   if (niter_file==0) then       ! Calculation from scratch.
     phi_nm1=ket0(my_t1:my_t2)   ! Select the slice treated by this node.
     norm = DZNRM2(hsize,ket0,1) ! Normalization  
     phi_nm1=phi_nm1/norm      
                                                                                
     ! hphi_n = MATMUL(hmat,phi_nm1)
     call xgemv('N',hsize,my_nt,cone,hmat,hsize,phi_nm1,1,czero,hphi_n,1)
     call xsum_mpi(hphi_n,comm,ierr)

     aa(1)=xdotc(my_nt,phi_nm1,1,hphi_n(my_t1:),1)
     call xsum_mpi(aa(1:1),comm,ierr)

     phi_n = hphi_n(my_t1:my_t2) - aa(1)*phi_nm1

     bb(1) = xdotc(my_nt,phi_n,1,phi_n,1)
     call xsum_mpi(bb(1:1),comm,ierr)
     bb(1) = SQRT(bb(1))

     phi_n = phi_n/bb(1)
     niter_done=1

   else ! Use the previous a and b.
     niter_done=niter_file
     aa(1:niter_done) = aa_file
     bb(1:niter_done) = bb_file
     phi_nm1=phi_nm1_file(my_t1:my_t2)   ! Select the slice treated by this node.
     phi_n  =phi_n_file  (my_t1:my_t2)   
   end if

   if (associated(aa_file     ))  then
     ABI_DEALLOCATE(aa_file)
   end if
   if (associated(bb_file     ))  then
     ABI_DEALLOCATE(bb_file)
   end if
   if (associated(phi_nm1_file))  then
     ABI_DEALLOCATE(phi_nm1_file)
   end if
   if (associated(phi_n_file  ))  then
     ABI_DEALLOCATE(phi_n_file)
   end if

   ! Multiplicative factor (k-point sampling and unit cell volume)  
   ! TODO be careful with the spin here
   ! TODO four_pi comes from the coulomb term 1/|q| is already included in the 
   ! oscillators hence the present approach wont work if a cutoff interaction is used.
   nfact = -four_pi/(Cryst%ucvol*BSp%nkbz)
   if (nsppol==1) nfact=two*nfact 

   factor = nfact*(DZNRM2(hsize,ket0,1)**2)

   ! Which quantity should be checked for convergence?
   check = (/.TRUE.,.TRUE./) 
   if (ABS(Bsp%haydock_tol(2)-one)<tol6) check = (/.TRUE. ,.FALSE./) 
   if (ABS(Bsp%haydock_tol(2)-two)<tol6) check = (/.FALSE.,.TRUE./) 

   ! Create new frequencies "mirror" in negative range to add 
   ! their contributions. Can be improved by computing only once
   ! zero frequency, but loosing clearness
   n_all_omegas = 2*BSp%nomega

   ABI_ALLOCATE(all_omegas,(n_all_omegas))
   ! Put all omegas with frequency > 0 in table
   all_omegas(BSp%nomega+1:n_all_omegas) = BSp%omega
   ! Put all omegas with frequency < 0
   ! Warning, the broadening must be kept positive
   all_omegas(1:BSp%nomega) = -DBLE(BSp%omega(BSp%nomega:1:-1)) &
& + j_dpc*AIMAG(BSp%omega(BSp%nomega:1:-1))   

   ABI_ALLOCATE(green_temp,(n_all_omegas,nkets))

   ! Calling haydock_herm_algo with green_temp with full range of frequencies
   call haydock_herm_algo(niter_done,niter_max,n_all_omegas,all_omegas,BSp%haydock_tol(1),check,hsize,&
&    my_t1,my_t2,hmat,factor,term_type,aa,bb,phi_nm1,phi_n,green_temp(:,iq),inn,is_converged,comm)

   ! Computing result from two ranges of frequencies
   ! The real part is added, the imaginary part is substracted
   green(:,iq) = green_temp(BSp%nomega+1:n_all_omegas,iq)+CONJG(green_temp(BSp%nomega:1:-1,iq))

   ABI_DEALLOCATE(all_omegas)
   ABI_DEALLOCATE(green_temp)
   !
   ! Save the a"s and the b"s for possible restarting.
   ! 1) Info on the Q.
   ! 2) Number of iterations performed.
   ! 3) do iter=1,niter_performed 
   !      aa(iter),bb(iter)
   !    end do
   ! 4) |n-1>
   !    |n>
   !
   if (my_rank==master) then ! Open the file and writes basic dimensions and info.
     write(out_unt)Bsp%q(:,iq)
     write(out_unt)MIN(inn,niter_max)  ! NB if the previous loop completed inn=niter_max+1
     do it=1,MIN(inn,niter_max)        ! if we exited then inn is not incremented by one.
       write(out_unt)it,aa(it),bb(it)
     end do
   end if
   !
   ! hphi_n is used as workspace to gather |n> and |n+1>.
   hphi_n = czero
   hphi_n(my_t1:my_t2) = phi_nm1
   call xsum_master(hphi_n,master,comm,ierr)
   if (my_rank==master) write(out_unt)hphi_n ! |n-1>

   hphi_n = czero
   hphi_n(my_t1:my_t2) = phi_n
   call xsum_master(hphi_n,master,comm,ierr)
   if (my_rank==master) write(out_unt)hphi_n ! |n>

   ABI_DEALLOCATE(hphi_n)
   ABI_DEALLOCATE(phi_nm1)
   ABI_DEALLOCATE(phi_n)
   ABI_DEALLOCATE(aa)
   ABI_DEALLOCATE(bb)
   ABI_DEALLOCATE(ket0)
 end do ! iq

 if (my_rank==master) close(out_unt)

 call xbarrier_mpi(comm)

end subroutine haydock_herm
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/haydock_herm_algo
!! NAME
!! haydock_herm_algo
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  niter_done=Number of iterations already performed (0 if the run starts from scratch).
!!  niter_max=Max number of iterations. Always > niter_done
!!  nomega=Number of Frequency points for the evaluation of the matrix element.
!!  omega(nomega)=Frequency set (imaginary part is already included).
!!  tol_iter=Tolerance used to stop the algorithm.
!!  check(2)=Logical flags to specify where both the real and the imaginary part of the 
!!    matrix elements of the Green functions have to be checked for convergence. 
!!  hsize=Size of the blocks.
!!  my_t1,my_t2=Indices of the first and last column stored treated by this done.
!!  term_type=0 if no terminator is used, 1 otherwise.
!!  hmat(hsize,my_t1:my_t2)=The columns of the block.
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
!!  aa(niter_max) and bb(niter_max)
!!    if niter_done>0: aa(1:niter_done), bb(1:niter_done) store the coefficients of the previous run.
!!    when the routine returns aa(1:inn) and bb(1:inn) contain the matrix elements of the tridiagonal form.
!!
!! PARENTS
!!      haydock,m_shexc
!!
!! CHILDREN
!!      matrginv,zgesv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine haydock_herm_algo(niter_done,niter_max,nomega,omega,tol_iter,check,hsize,my_t1,my_t2,hmat,&
&  factor,term_type,aa,bb,phi_nm1,phi_n,green,inn,is_converged,comm)

 use m_profiling

 use defs_basis
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_numeric_tools,  only : continued_fract
 use m_blas,           only : xdotc, xgemv

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_herm_algo'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: niter_max,niter_done,nomega
 integer,intent(in) :: comm,hsize,my_t1,my_t2,term_type
 integer,intent(out) :: inn
 logical,intent(out) :: is_converged
 real(dp),intent(in) :: tol_iter
 complex(dpc),intent(in) :: factor
!arrays
 real(dp),intent(inout) :: bb(niter_max)
 complex(dpc),intent(out) :: green(nomega)
 complex(dpc),intent(in) :: omega(nomega) 
 complex(dpc),intent(inout) :: aa(niter_max)
 complex(dpc),intent(in) :: hmat(hsize,my_t1:my_t2)
 complex(dpc),intent(inout) :: phi_nm1(my_t2-my_t1+1)
 complex(dpc),intent(inout) :: phi_n  (my_t2-my_t1+1)
 logical,intent(in) :: check(2)

!Local variables ------------------------------
!scalars
 integer :: ierr,istat,my_nt,niter_min,nconv
 character(len=500) :: msg
 logical,parameter :: force_real=.TRUE.
!arrays
 real(dp) :: abs_err(nomega,2) !,rel_err(nomega,2)
 complex(dpc),allocatable :: oldg(:),newg(:)
 complex(dpc),allocatable :: phi_np1(:),hphi_n(:),cfact(:)
 logical :: test(2)

!************************************************************************

 ! The sequences starts with |1> normalized to 1 and b_0 =0, therefore:
 !  a_1 = <1|H|1>
 !  b_1 = || H|1> - a_1|1> ||
 !  |2> = [H|1> - a_1|1>]/b_1
 !
 ! For n>1 we have
 !  1) a_n = <n|H|n>
 !  2) b_n = || H|n> - a_n|n> -b_{n-1}|n-1> ||
 !  3) |n+1> = [H|n> -a_n|n> -b_{n-1}|n-1>]/b_n
 !
 my_nt = my_t2-my_t1+1

 ABI_ALLOCATE(hphi_n,(hsize))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out-of-memory hphi_n")
 
 ABI_ALLOCATE(phi_np1,(my_nt))

 ABI_ALLOCATE(oldg,(nomega))
 oldg=czero 
 ABI_ALLOCATE(newg,(nomega))
 newg=czero 
 ABI_ALLOCATE(cfact,(nomega))
 cfact=czero 

 nconv=0
 do inn=niter_done+1,niter_max
   !
   ! hphi_n = MATMUL(hmat,phi_n)
   call xgemv('N',hsize,my_nt,cone,hmat,hsize,phi_n,1,czero,hphi_n,1)
   call xsum_mpi(hphi_n,comm,ierr)

   aa(inn) = xdotc(my_nt,phi_n,1,hphi_n(my_t1:),1)
   call xsum_mpi(aa(inn:inn),comm,ierr)
   if (force_real) aa(inn) = DBLE(aa(inn)) ! Matrix is Hermitian.

   ! |n+1> = H|n> - A(n)|n> - B(n-1)|n-1>
   phi_np1 = hphi_n(my_t1:my_t2) - aa(inn)*phi_n - bb(inn-1)*phi_nm1

   bb(inn) = xdotc(my_nt,phi_np1,1,phi_np1,1)
   call xsum_mpi(bb(inn),comm,ierr)
   bb(inn) = SQRT(bb(inn))

   phi_np1 = phi_np1/bb(inn)
   
   phi_nm1 = phi_n
   phi_n   = phi_np1

   write(msg,'(a,i0,a,3es12.4)')' Iteration number ',inn,', b_i RE(a_i) IM(a_i) ',bb(inn),REAL(aa(inn)),AIMAG(aa(inn)) 
   call wrtout(std_out,msg,"COLL")
   !if (MOD(inn,2)==0) then
   !  write(100,*)inn,bb(inn),REAL(aa(inn)),AIMAG(aa(inn)) 
   !else 
   !  write(101,*)inn,bb(inn),REAL(aa(inn)),AIMAG(aa(inn)) 
   !end if
   call continued_fract(inn,term_type,aa,bb,nomega,omega,cfact)

   newg= factor*cfact
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

   oldg = newg
 end do ! inn

 green = newg
 if (nconv/=2) then
   write(msg,'(a,es10.2,a,i0,a)')&
&    " WARNING: Haydock algorithm did not converge within ",tol_iter," after ",niter_max," iterations."
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 is_converged = (nconv==2)

 ABI_DEALLOCATE(oldg)
 ABI_DEALLOCATE(newg)
 ABI_DEALLOCATE(cfact)
 ABI_DEALLOCATE(hphi_n)
 ABI_DEALLOCATE(phi_np1)

end subroutine haydock_herm_algo
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/haydock_restart
!! NAME
!! haydock_restart
!!
!! FUNCTION
!! Restart the Haydock method from file reading the data produced in a previous run.
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
!!  iq_search=The index of the q-point to be searched.
!!  hsize
!!  comm=MPI communicator.
!!  nsppol
!!  restart_file
!!
!! OUTPUT
!!  niter_file=Number of iterations already performed. 0 to signal that an error occurred during the reading
!!
!! SIDE EFFECTS
!!  bb_file(:)
!!  aa_file(:)
!!  phi_n_file(:)
!!  phi_nm1_file(:)
!!
!! PARENTS
!!      haydock,haydock_psherm
!!
!! CHILDREN
!!      matrginv,zgesv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine haydock_restart(BSp,restart_file,ftype,iq_search,hsize,niter_file,aa_file,bb_file,phi_nm1_file,phi_n_file,comm)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use m_xmpi
 use m_errors
#if defined HAVE_MPI2
 use mpi
#endif

 use m_io_tools,  only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_restart'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,hsize,iq_search,ftype
 integer,intent(out) :: niter_file
 character(len=*),intent(in) :: restart_file
 type(excparam),intent(in) :: BSp
!arrays
 real(dp),pointer :: bb_file(:)
 complex(dpc),pointer :: aa_file(:),phi_n_file(:),phi_nm1_file(:)

!Local variables ------------------------------
!scalars
 integer :: inn,it,restart_unt,ios,nproc,my_rank,master,ierr,op_file
 integer :: hsize_file,nq_file,iq_file,use_coupling_file
 real(dp) :: broad_file
 logical :: found_q
 character(len=500) :: msg
!arrays
 real(dp) :: qfile(3)

!************************************************************************

 nproc  = xcomm_size(comm)
 my_rank= xcomm_rank(comm)
 master = 0
 !
 if (my_rank==master) then
   restart_unt = get_unit()
   open(unit=restart_unt,file=restart_file,form="unformatted",status="old",iostat=ios)
   ABI_CHECK(ios==0," Opening file: "//TRIM(restart_file))

   read(restart_unt)hsize_file,use_coupling_file,op_file,nq_file,broad_file 
   !write(std_out,*)"hsize_file",hsize_file,nq_file,broad_file

   if (op_file/=ftype) then
     write(msg,"(2(a,i0))")" Expecting restart file with filetype: ",ftype," but found ",op_file
     MSG_ERROR(msg)
   end if

   if (hsize_file/=hsize) then
     write(msg,"(2(a,i0))")&
&      " Rank of H_exc read from file: ",hsize_file," differs from the one used in this run: ",hsize
     MSG_ERROR(msg)
   end if

   if (use_coupling_file /= BSp%use_coupling) then
     write(msg,'(2(a,i0))')&
&      " use_coupling_file: ",use_coupling_file," differs from input file value: ",BSp%use_coupling
     MSG_ERROR(msg)
   end if

   found_q=.FALSE.
   do iq_file=1,nq_file
     read(restart_unt) qfile(:)
     read(restart_unt)niter_file
     if ( ALL(ABS(qfile-BSp%q(:,iq_search)) < tol6) ) then
       found_q=.TRUE.; EXIT
     else 
       ! Skip data for this q.
       do it=1,niter_file
         read(restart_unt) ! it,aa(it),bb(it)
       end do
       read(restart_unt)
       read(restart_unt)
     end if
   end do

   if (.not.found_q) then
     niter_file=0
     write(msg,"(a,3f8.4,3a)")&
&      " Could not find q-point: ",BSp%q(:,iq_search)," in file ",TRIM(restart_file),&
&      " Cannot restart Haydock iterations for this q-point"
     MSG_COMMENT(msg)
   else
     write(msg,'(a,i0)')" Number of iterations already performed: ",niter_file
     call wrtout(std_out,msg,"COLL")
     call wrtout(ab_out,msg,"COLL")

     if ( ABS(broad_file - BSp%broad) > tol6) then
       write(msg,'(2a,2(a,f8.4),a)')&
&        " Restart file has been produced with a different Lorentzian broadening: ",ch10,&
&        " broad_file: ",broad_file," input broadening: ",BSp%broad," Continuing anyway. "
       MSG_WARNING(msg)
     end if

     ABI_ALLOCATE(aa_file,(niter_file))
     ABI_ALLOCATE(bb_file,(niter_file))
     do inn=1,niter_file
       read(restart_unt)it,aa_file(inn),bb_file(inn)
       if (inn/=it) then 
         write(msg,'(2(a,i0))')" Found it_file: ",it," while it should be: ",inn
         MSG_ERROR(msg)
       end if
     end do
     ABI_ALLOCATE(phi_nm1_file,(hsize))
     ABI_ALLOCATE(phi_n_file,(hsize))
     read(restart_unt)phi_nm1_file
     read(restart_unt)phi_n_file
   end if
   close(restart_unt)
 end if
 !
 ! Master broadcasts the data.
 call xcast_mpi(niter_file,master,comm,ierr)

 if (my_rank/=master) then 
   ABI_ALLOCATE(aa_file,(niter_file))
   ABI_ALLOCATE(bb_file,(niter_file))
   ABI_ALLOCATE(phi_nm1_file,(hsize))
   ABI_ALLOCATE(phi_n_file,(hsize))
 end if

 call xcast_mpi(aa_file,master,comm,ierr)
 call xcast_mpi(bb_file,master,comm,ierr)
 call xcast_mpi(phi_nm1_file,master,comm,ierr)
 call xcast_mpi(phi_n_file,master,comm,ierr)

end subroutine haydock_restart
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/haydock_mdf_to_tensor
!! NAME
!! haydock_mdf_to_tensor
!!
!! FUNCTION
!! Transform macroscopic dielectric function from green function to each components of the tensor in red and cart coord.
!! WARNING : Works only with 6 q-points in optical limit
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (YG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  BSp<type(excparam)>=Parameters defining the Bethe-Salpeter calculation.
!!    omega(BSp%nomega)=Frequency mesh for the macroscopic dielectric function (broadening is already included).
!!  Cryst=Parameters of the crystal
!!  eps(BSp%nomega,BSp%nq) = Macroscopic dielectric function to be written.
!!
!! OUTPUT
!!  tensor_cart(BSp%nomega,6) = dielectric tensor for each frequency, order (11,22,33,12,13,23) in cart. coord.
!!  tensor_red(BSp%nomega, 6) = idem in reduced coordinated
!!
!! PARENTS
!!      haydock
!!
!! CHILDREN
!!      matrginv,zgesv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine haydock_mdf_to_tensor(BSp,Cryst,eps,tensor_cart,tensor_red)

 use m_profiling

 use defs_basis
 use m_bs_defs
 use m_errors

 use m_crystal,       only : crystal_structure
 use m_geometry,      only : normv

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'haydock_mdf_to_tensor'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(excparam),intent(in) :: BSp
 type(crystal_structure),intent(in) :: Cryst
!arrays
 complex(dpc),intent(in) :: eps(BSp%nomega,BSp%nq)
 complex(dpc),intent(out) :: tensor_cart(BSp%nomega,6), tensor_red(BSp%nomega,6)

!Local variables ------------------------------
!scalars
 integer :: iq,info
 real(dp) :: normqcart, normqred
!arrays
 integer,allocatable :: ipiv(:)
 real(dp) :: qcart(3), qtmet(3)
 real(dp) :: qred2cart(3,3),qcart2red(3,3)
 complex(dpc) :: qqcart(BSp%nq,6), qqred(BSp%nq,6)
 complex(dpc) :: b(6,BSP%nomega)

!************************************************************************

 ! Transformation matrices from reduced coordinates to cartesian coordinates
 qred2cart = two_pi*Cryst%gprimd
 qcart2red = qred2cart
 call matrginv(qcart2red,3,3)
 do iq = 1, BSp%nq

   ! Computing cartesian q-vector
   qcart = MATMUL(qred2cart, BSp%q(:,iq))

   ! Computing product 'metric - qred' to form quadratic form
   qtmet = (two_pi**2)*MATMUL(Cryst%gmet, BSp%q(:,iq))
 
   ! squared norms
   normqcart = qcart(1)**2+qcart(2)**2+qcart(3)**2
   normqred = (normv(BSp%q(:,iq),Cryst%gmet,"G"))**2

   ! Compute line 'iq' for matrix in cartesian coord
   qqcart(iq,1) = (qcart(1))**2
   qqcart(iq,2) = (qcart(2))**2
   qqcart(iq,3) = (qcart(3))**2
   qqcart(iq,4) = 2*(qcart(1)*qcart(2))
   qqcart(iq,5) = 2*(qcart(1)*qcart(3))
   qqcart(iq,6) = 2*(qcart(2)*qcart(3))

   ! Compute line 'iq' for matrix in reduced coord
   qqred(iq,1) = (qtmet(1))**2
   qqred(iq,2) = (qtmet(2))**2
   qqred(iq,3) = (qtmet(3))**2
   qqred(iq,4) = 2*(qtmet(1)*qtmet(2))
   qqred(iq,5) = 2*(qtmet(1)*qtmet(3))
   qqred(iq,6) = 2*(qtmet(2)*qtmet(3))

   ! Renormalize line
   qqcart(iq,:) = qqcart(iq,:)/normqcart
   qqred(iq,:) = qqred(iq,:)/normqred

 end do

 ABI_ALLOCATE(ipiv,(6))

 ! Solving linear system
 b = TRANSPOSE(eps)
 call ZGESV(6,BSp%nomega,qqcart,6,ipiv, b, 6, info)
 tensor_cart = TRANSPOSE(b)

 b = TRANSPOSE(eps)
 call ZGESV(6,BSp%nomega,qqred,6,ipiv,b,6,info)
 tensor_red = TRANSPOSE(b)

 ABI_DEALLOCATE(ipiv)

end subroutine haydock_mdf_to_tensor
!!***
