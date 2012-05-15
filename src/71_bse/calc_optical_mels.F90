!!****f* ABINIT/calc_optical_mels
!! NAME
!!  calc_optical_mels
!!
!! FUNCTION
!!  Calculate all optical matrix elements in the BZ.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! minb,maxb=minimum and max band index to be calculated.
!! nkbz=Number of points in the full Brillouin zone.
!! inclvkb=if different from 0, [Vnl,r] is included in the calculation of the matrix element of the velocity operator
!!   No meaning for PAW (except for LDA+U)
!! qpt(3)
!! Kmesh<bz_mesh_type>=Info on the k-point sampling for wave functions. 
!! Cryst<crystal_structure>=Structure defining the crystalline structure.
!! KS_Bst<bandstructure_type>
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data
!! Psps <pseudopotential_type>=variables related to pseudopotentials.
!! Hur(Cryst%natom*usepaw)<HUr_commutator>=Only for PAW and LDA+U, quantities used to evaluate the commutator [H_u,r].
!! Wfd<wfs_descriptor>=Handler for the wavefunctions.
!!   nkibz=Number of irreducible k-points.
!!   nsppol=Number of independent spin polarizations.
!!   usepaw=1 for PAW, 0 otherwise.
!!   nspinor=Number of spinorial components.
!!   spaceComm=MPI Communicator.
!!
!! OUTPUT
!! opt_cvk(minb:maxb,minb:maxb,nkbz,nsppol)=Matrix elements <c k|e^{+iqr}|v k> for the different qs.
!!
!! PARENTS
!!      exc_spectra,haydock
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,destroy_kb_potential,get_bz_item,init_kb_potential
!!      matrginv,nullify_kb_potential,wfd_distribute_bbp,wfd_get_cprj,wrtout
!!      xbarrier_mpi,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_optical_mels(Wfd,Kmesh,KS_Bst,Cryst,Psps,Pawtab,Hur,inclvkb,minb,maxb,nkbz,qpoint,opt_cvk)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_bs_defs
 use m_xmpi
 use m_errors

 use m_bz_mesh,           only : bz_mesh_type, get_BZ_item
 use m_crystal,           only : crystal_structure
 use m_paw_commutator,    only : HUr_commutator, paw_ihr_comm
 use m_commutator_vkbr,   only : kb_potential, nullify_kb_potential, destroy_kb_potential, init_kb_potential, &
&                                nc_ihr_comm
 use m_wfs,               only : wfs_descriptor, wfd_get_cprj, wfd_distribute_bbp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_optical_mels'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkbz,inclvkb,minb,maxb
 type(bz_mesh_type),intent(in) :: Kmesh
 type(crystal_structure),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(Bandstructure_type),intent(in) :: KS_Bst
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 real(dp),intent(in) :: qpoint(3)
 complex(dpc),intent(out) :: opt_cvk(minb:maxb,minb:maxb,nkbz,Wfd%nsppol)
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(HUr_commutator),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: nsppol,usepaw,nspinor,spaceComm
 integer :: ik_bz,ik_ibz,itim_k,isym_k,ib_c,ib_v,ierr,my_rank
 integer :: spin,npw_k,istwf_k,my_nbbp
 real(dp) :: ediff
 complex(dpc) :: emcvk
 type(kb_potential) :: KBgrad_k
!arrays
 integer :: bbp_distrb(Wfd%mband,Wfd%mband),got(Wfd%nproc)
 integer,pointer :: kg_k(:,:)
 real(dp) :: mat_dp(3,3),qrot(3),b1(3),b2(3),b3(3),kbz(3)
 complex(dpc),allocatable :: ir_kibz(:,:,:,:,:)
 complex(gwpc),pointer :: ug_c(:),ug_v(:)
 complex(gwpc) :: ihrc(3,Wfd%nspinor**2)
 logical :: bbp_mask(Wfd%mband,Wfd%mband)
 type(Cprj_type),allocatable :: Cp_v(:,:),Cp_c(:,:)

!************************************************************************

 call wrtout(std_out," Calculating optical matrix elements in the IBZ","COLL")
 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")

 spaceComm = Wfd%comm
 my_rank   = Wfd%my_rank

 nsppol  = Wfd%nsppol   
 nspinor = Wfd%nspinor
 usepaw  = Wfd%usepaw

 if (usepaw==1) then
   ABI_ALLOCATE(Cp_v,(Wfd%natom,nspinor))
   call cprj_alloc(Cp_v,0,Wfd%nlmn_atm)
   ABI_ALLOCATE(Cp_c,(Wfd%natom,nspinor))
   call cprj_alloc(Cp_c,0,Wfd%nlmn_atm)
 end if

 if (inclvkb==1.and.usepaw==0) then
   MSG_ERROR("inclvkb==1 not coded,using inclvkb==2")
 end if
 !
 ! Calculate the matrix elements of ir in the IBZ.
 call nullify_kb_potential(KBgrad_k)

 ABI_ALLOCATE(ir_kibz,(3,minb:maxb,minb:maxb,Wfd%nkibz,nsppol))
 ir_kibz=czero

 got=0; 
 bbp_mask=.FALSE.; bbp_mask(minb:maxb,minb:maxb)=.TRUE.
 do spin=1,nsppol
   do ik_ibz=1,Wfd%nkibz
    !
    ! Distribute the (b,b') entries.
    call wfd_distribute_bbp(Wfd,ik_ibz,spin,"All",my_nbbp,bbp_distrb,got=got,bbp_mask=bbp_mask) 
    if ( ALL(bbp_distrb/=my_rank) ) CYCLE

    istwf_k = Wfd%istwfk(ik_ibz)
    ABI_CHECK(istwf_k==1,"istwf_k/=1 not coded") ! KB stuff is missing.
    npw_k = Wfd%npwarr(ik_ibz)
    kg_k  => Wfd%Kdata(ik_ibz)%kg_k

    if (inclvkb/=0.and.usepaw==0) then ! Prepare term i <n,k|[Vnl,r]|n"k>
      call init_kb_potential(KBgrad_k,Cryst,Psps,inclvkb,istwf_k,npw_k,Kmesh%ibz(:,ik_ibz),kg_k)
    end if

    ! Note: spinorial case is not coded therefore we work with ihrc(:,1).
    ! TODO: The lower triangle can be Reconstructed by symmetry.
    do ib_v=minb,maxb ! Loop over bands
      if ( ALL(bbp_distrb(ib_v,:)/=my_rank) ) CYCLE
      ug_v => Wfd%Wave(ib_v,ik_ibz,spin)%ug
      if (usepaw==1) call wfd_get_cprj(Wfd,ib_v,ik_ibz,spin,Cryst,Cp_v,sorted=.FALSE.)

      do ib_c=minb,maxb
       if (bbp_distrb(ib_v,ib_c)/=my_rank) CYCLE
       ug_c => Wfd%Wave(ib_c,ik_ibz,spin)%ug

       if (usepaw==0) then  ! Calculate matrix elements of i[H,r] for NC pseudopotentials.        
         ihrc = nc_ihr_comm(nspinor,npw_k,istwf_k,inclvkb,Kmesh%ibz(:,ik_ibz),KBgrad_k,ug_c,ug_v,kg_k) 

       else ! Matrix elements of i[H,r] for PAW.
         call wfd_get_cprj(Wfd,ib_c,ik_ibz,spin,Cryst,Cp_c,sorted=.FALSE.)

         ihrc = paw_ihr_comm(spin,nspinor,npw_k,istwf_k,Kmesh%ibz(:,ik_ibz),Cryst,Pawtab,ug_c,ug_v,kg_k,Cp_c,Cp_v,HUr)
       end if
       !
       ! Save matrix elements of i*r in the IBZ
       ediff = KS_Bst%eig(ib_c,ik_ibz,spin) - KS_BSt%eig(ib_v,ik_ibz,spin)
       if (ABS(ediff)<tol16) ediff=tol6  ! Treat a possible degeneracy between v and c.
       ir_kibz(:,ib_c,ib_v,ik_ibz,spin) = ihrc(:,1)/ediff

      end do !ib_c
    end do !ib_v

    call destroy_kb_potential(KBgrad_k)

   end do !spin
 end do !ik_ibz

 ! Collect results on each node.
 call xsum_mpi(ir_kibz,spaceComm,ierr)

 if (usepaw==1) then
   call cprj_free(Cp_v)
   ABI_DEALLOCATE(Cp_v)
   call cprj_free(Cp_c)
   ABI_DEALLOCATE(Cp_c)
 end if
 !
 ! ======================================================
 ! ==== Calculate Fcv(kBZ) in the full Brilouin zone ====
 ! ======================================================
 !
 ! Symmetrization of the matrix elements of the position operator.
 ! <Sk b|r|Sk b'> = R <k b|r|k b'> + \tau \delta_{bb'}
 !   where S is one of the symrec operations in reciprocal space, R is the 
 !   corresponding operation in real space, \tau being the associated fractional translations.
 !
 ! q.Mcv( Sk) =  S^{-1}q. Mcv(k) 
 ! q.Mcv(-Sk) = -S^{-1}q. CONJG(Mcv(k)) if time-reversal is used.

 b1=Cryst%gprimd(:,1)*two_pi
 b2=Cryst%gprimd(:,2)*two_pi
 b3=Cryst%gprimd(:,3)*two_pi

 opt_cvk = czero
 do spin=1,nsppol
   !
   do ik_bz=1,nkbz
    !
    ! * Get ik_ibz, and symmetries index from ik_bz.
    call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)
    
    mat_dp = DBLE(Cryst%symrec(:,:,isym_k))
    call matrginv(mat_dp,3,3) ! Invert
    qrot = (3-2*itim_k) * MATMUL(mat_dp,qpoint)

    do ib_v=minb,maxb !  Loops over the bands C and V start
      do ib_c=minb,maxb
        !if (ib_c==ib_v) CYCLE 
        emcvk = pdtqrc(qrot,ir_kibz(:,ib_c,ib_v,ik_ibz,spin),b1,b2,b3)
        if (itim_k==2) emcvk = CONJG(emcvk)
        opt_cvk(ib_c,ib_v,ik_bz,spin) = emcvk
      end do !ib_c
    end do !ib_v

   end do !ik_bz
 end do !spin

 ABI_DEALLOCATE(ir_kibz)

 call xbarrier_mpi(spaceComm)

end subroutine calc_optical_mels
!!***
