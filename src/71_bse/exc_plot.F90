!!****f* ABINIT/exc_plot
!! NAME
!!  exc_plot
!!
!! FUNCTION
!!  Plots the excitonic wavefunction in real space.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Bsp<excparam>=Container storing the input variables of the run.
!! BS_files<excfiles>=filenames used in the Bethe-Salpeter part.
!! Kmesh<bz_mesh_type>=Info on the k-point sampling for wave functions. 
!! Cryst<crystal_structure>=Structure defining the crystalline structure.
!! KS_Bst<bandstructure_type>
!! Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data.
!! Psps <pseudopotential_type>=variables related to pseudopotentials.
!! Wfd<wfs_descriptor>=Handler for the wavefunctions.
!! ngfftf(18)=Info on the dense FFT mesh used for plotting the excitonic wavefunctions.
!! nrcell(3)=Number of cell replicas (the code will enforce odd number so that the e or the 
!!  h can be centered in the bix box.
!! eh_rcoord(3)=Reduced coordinated of the (electron|hole)
!! which_fixed= 1 to fix the electron at eh_red, 2 for the hole.
!! spin_opt=1 for resolving the up-component, 2 for the down component, 3 if spin resolution is not wanted.
!! paw_add_onsite=.TRUE. if the onsite contribution is taken into account.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,destroy_paw_pwaves_lmn,exc_read_eigen
!!      init_paw_pwaves_lmn,nhatgrid,pawfgrtab_free,pawfgrtab_init
!!      pawfgrtab_print,printxsf,wfd_change_ngfft,wfd_sym_ur,wrap2_zero_one
!!      wrtout,xcast_mpi,xsum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine exc_plot(Bsp,Bs_files,Wfd,Kmesh,Cryst,Psps,Pawtab,Pawrad,paw_add_onsite,spin_opt,which_fixed,eh_rcoord,nrcell,ngfftf)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_bs_defs
 use m_xmpi
 use m_errors

 use m_io_tools,          only : get_unit
 use m_numeric_tools,     only : iseven
 use m_bz_mesh,           only : bz_mesh_type, get_BZ_item
 use m_crystal,           only : crystal_structure
 use m_wfs,               only : wfs_descriptor, wfd_change_ngfft, wfd_get_cprj, wfd_sym_ur
 use m_paw_toolbox,       only : pawfgrtab_init, pawfgrtab_free, pawfgrtab_print,&
&                                init_paw_pwaves_lmn, destroy_paw_pwaves_lmn, paw_pwaves_lmn_t
 use m_bse_io,            only : exc_read_eigen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exc_plot'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_44_abitypes_defs
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: which_fixed,spin_opt
 logical,intent(in) :: paw_add_onsite
 type(excparam),intent(in) :: Bsp
 type(excfiles),intent(in) :: BS_files
 type(bz_mesh_type),intent(in) :: Kmesh
 type(crystal_structure),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18),nrcell(3)
 real(dp),intent(in) :: eh_rcoord(3)
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(in) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: cplex1=1
 integer :: nsppol,usepaw,nspinor,comm,cond,val,istat,ios
 integer :: optcut,optgr0,optgr1,optgr2,optrad
 integer :: ik_bz,ierr,my_rank !ik_ibz, istwf_k, isym_k,itim_k, my_nbbp, npw_k,
 integer :: spin,spin_start,spin_stop,reh 
 integer :: rt_idx,art_idx,ii,iatom,nr_tot
 integer :: irc,ir1,ir2,ir3,wp1,wp2,wp3,wp_idx,eh_fft_idx,eh_rr,rr
 integer :: hsize,xsf_unt,ncells,nvec
 integer :: sc_natom,master
 real(dp) :: ene_rt,k_dot_r12
 complex(dpc) :: res_coeff,antres_coeff,eikr12
 character(len=500) :: msg
 character(len=fnlen) :: eig_fname,out_fname       
!arrays
 integer :: nrcl(3),bbox(3)
 !integer :: bbp_distrb(Wfd%mband,Wfd%mband),got(Wfd%nproc)
 integer,allocatable :: rcl2fft(:),sc_typat(:),l_size_atm(:),vec_idx(:)
 real(dp),parameter :: origin0(3)=(/zero,zero,zero/)
 real(dp) :: k_bz(3),eh_red(3),eh_shift(3),r12(3),sc_rprimd(3,3),scart_shift(3)
 real(dp),allocatable :: rclred(:,:),exc_phi2(:),sc_xcart(:,:) !,sc_znucl(:)
 complex(gwpc),allocatable :: exc_phi(:)
 complex(gwpc),allocatable :: ur_v(:),ur_c(:)
 complex(dpc),allocatable :: vec_list(:,:)
 !logical :: bbp_mask(Wfd%mband,Wfd%mband)
 type(Cprj_type),allocatable :: Cp_v(:,:),Cp_c(:,:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

!************************************************************************

 MSG_WARNING("Exc plot is still under development")

 call wrtout(std_out," Calculating excitonic wavefunctions in real space","COLL")
 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")
 ABI_CHECK(Wfd%usepaw==0,"PAW not coded")
 ABI_CHECK(Wfd%nproc==1,"nproc>1 not supported")

 ! If needed, prepare FFT tables to have u(r) on the ngfftf mesh.
 if ( ANY(ngfftf(1:3) /= Wfd%ngfft(1:3)) ) call wfd_change_ngfft(Wfd,Cryst,Psps,ngfftf) 

 comm    = Wfd%comm
 my_rank = Wfd%my_rank
 master  = Wfd%master

 nsppol  = Wfd%nsppol   
 nspinor = Wfd%nspinor
 usepaw  = Wfd%usepaw
 !
 ! TODO recheck this part.
 ! Prepare tables needed for splining the onsite term on the FFT box.
 if (Wfd%usepaw==1.and.paw_add_onsite) then
   MSG_WARNING("Testing the calculation of AE PAW wavefunctions.")
   ! Use a local pawfgrtab to make sure we use the correction in the paw spheres
   ! the usual pawfgrtab uses r_shape which may not be the same as r_paw.
   ABI_ALLOCATE(Pawfgrtab,(Cryst%natom))
   ABI_ALLOCATE(l_size_atm,(Cryst%natom))
   do iatom=1,Cryst%natom
     l_size_atm(iatom) = Pawtab(Cryst%typat(iatom))%l_size
   end do

   call pawfgrtab_init(Pawfgrtab,cplex1,l_size_atm,Wfd%nspden)
   ABI_DEALLOCATE(l_size_atm)

   optcut=1                     ! use rpaw to construct Pawfgrtab.
   optgr0=0; optgr1=0; optgr2=0 ! dont need gY terms locally.
   optrad=1                     ! do store r-R.

   call nhatgrid(Cryst%atindx1,Cryst%gmet,Wfd%MPI_enreg,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
&    optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%ucvol,Cryst%xred)
   !Pawfgrtab is ready to use

   if (Wfd%pawprtvol>0) then
     call pawfgrtab_print(Pawfgrtab,unit=std_out,prtvol=Wfd%pawprtvol,mode_paral="COLL")
   end if

   ABI_ALLOCATE(Paw_onsite,(Cryst%natom))
   call init_paw_pwaves_lmn(Paw_onsite,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%rprimd,&
&    Psps,Pawtab,Pawrad,Pawfgrtab)
 end if 
 !
 ! Number of cells in the big box. Odd numbers are needed to place the e or the h in the middle.
 do ii=1,3
   nrcl(ii) = nrcell(ii)
   if (iseven(nrcell(ii))) then
     nrcl(ii) = nrcell(ii)+1
     write(msg,'(2(a,i0))')" Enforcing odd number of cell replicas ",nrcell(ii)," --> ",nrcl(ii)
     MSG_WARNING(msg)
   end if
 end do

 ncells = PRODUCT(nrcl) 
 nr_tot = Wfd%nfftot * ncells  ! Total number of points in the big box.
 bbox = nrcl*ngfftf(1:3)
 !
 ! rcl2fft: The image of the point in the small box.
 ! rcl2red: The reduced coordinates of the point in the big box in terms of rprimd.
 ABI_ALLOCATE(rcl2fft,(nr_tot))
 ABI_ALLOCATE(rclred,(3,nr_tot))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0," out of memory in rcl tables")

 irc = 0
 do ir3=0,bbox(3)-1 ! Loop over the points in the big box.
   do ir2=0,bbox(2)-1
     do ir1=0,bbox(1)-1  
       irc = 1+irc
       !
       wp1=MODULO(ir1,ngfftf(1)) ! The FFT index of the point wrapped into the original cell.
       wp2=MODULO(ir2,ngfftf(2))
       wp3=MODULO(ir3,ngfftf(3))
       wp_idx = 1 + wp1 + wp2*ngfftf(1) + wp3*ngfftf(1)*ngfftf(2)
       rcl2fft(irc)  = wp_idx 
       rclred(1,irc) = DBLE(ir1)/ngfftf(1) ! Reduced coordinates in terms of the origina cell.
       rclred(2,irc) = DBLE(ir2)/ngfftf(2)
       rclred(3,irc) = DBLE(ir3)/ngfftf(3)
     end do
   end do
 end do
 !
 ! Wrap the point to [0,1[ where 1 is not included (tol12)
 call wrap2_zero_one(eh_rcoord(1),eh_red(1),eh_shift(1))
 call wrap2_zero_one(eh_rcoord(2),eh_red(2),eh_shift(2))
 call wrap2_zero_one(eh_rcoord(3),eh_red(3),eh_shift(3))
 write(std_out,*)"Initial Position of (e|h) in reduced coordinates:",eh_red," shift: ",eh_shift
 !
 ! Initial position on the FFT grid (the closest one)
 eh_red(:)=NINT(eh_red(:)*ngfftf(1:3))
 !
 ! Not translate the (electron|hole) in the center of the big box.
 do ii=1,3
   if (nrcl(ii)/=1) eh_red(ii) = eh_red(ii) + (nrcl(ii)/2) * ngfftf(ii)
 end do
 write(std_out,*)"Position of (e|h) in the center of the big box in FFT coordinates:",eh_red(:)

 eh_rr = 1 + eh_red(1) + eh_red(2)*bbox(1) + eh_red(3)*bbox(1)*bbox(2)
 eh_fft_idx = rcl2fft(eh_rr)

 write(std_out,*)" Reduced coordinates of (e|h) in the center of the big box ",rclred(:,eh_rr)
 !
 ! Allocate wavefunctions.
 ABI_ALLOCATE(ur_v,(Wfd%nfft*nspinor))
 ABI_ALLOCATE(ur_c,(Wfd%nfft*nspinor))

 if (usepaw==1) then
   ABI_ALLOCATE(Cp_v,(Wfd%natom,nspinor))
   call cprj_alloc(Cp_v,0,Wfd%nlmn_atm)
   ABI_ALLOCATE(Cp_c,(Wfd%natom,nspinor))
   call cprj_alloc(Cp_c,0,Wfd%nlmn_atm)
 end if
 !
 ! Read excitonic eigenvector.
 hsize = SUM(BSp%nreh); if (BSp%use_coupling>0) hsize=2*hsize
 nvec=1

 ABI_ALLOCATE(vec_list,(hsize,nvec))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in vec_list")

 ABI_ALLOCATE(vec_idx,(nvec))
 vec_idx = (/(ii, ii=1,nvec)/)

 if (BS_files%in_eig /= BSE_NOFILE) then
   eig_fname = BS_files%in_eig
 else 
   eig_fname = BS_files%out_eig
 end if

 if (my_rank==master) then
   call exc_read_eigen(eig_fname,hsize,nvec,vec_idx,vec_list,Bsp=Bsp)
 end if

 call xcast_mpi(vec_list,master,comm,ierr)

 ABI_DEALLOCATE(vec_idx)
 !
 ! Allocate the excitonic wavefunction on the big box.
 ABI_ALLOCATE(exc_phi,(nr_tot))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0," out of memory in exc_phi")
 exc_phi=czero
 !
 !got=0; 
 !bbp_mask=.FALSE.; bbp_mask(minb:maxb,minb:maxb)=.TRUE.
 spin_start=1; spin_stop=Wfd%nsppol
 if (spin_opt==1) spin_stop =1
 if (spin_opt==2) spin_start=2

 ! TODO 
 ! Distribute (b,b',k,s) states where k is in the *full* BZ.
 !call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,got,bmask)

 do spin=spin_start,spin_stop
   do reh=1,Bsp%nreh(spin)
     ene_rt = Bsp%Trans(reh,spin)%en
     ik_bz  = Bsp%Trans(reh,spin)%k
     val    = Bsp%Trans(reh,spin)%v
     cond   = Bsp%Trans(reh,spin)%c
     k_bz   = Kmesh%bz(:,ik_bz)
     !
     ! The resonant component.
     rt_idx = reh + (spin-1)*Bsp%nreh(1) 
     res_coeff = vec_list(rt_idx,1)
     !
     ! The anti-resonant component.
     antres_coeff = czero
     if (Bsp%use_coupling>0) then 
       art_idx  = reh + (spin-1)*Bsp%nreh(1) + SUM(Bsp%nreh) 
       antres_coeff = vec_list(art_idx,1)
     end if

     call wfd_sym_ur(Wfd,Cryst,Kmesh,val ,ik_bz,spin,ur_v) ! TODO recheck this one.
     call wfd_sym_ur(Wfd,Cryst,Kmesh,cond,ik_bz,spin,ur_c)

     !call wfd_paw_get_aeur(Wfd,band,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae,ur_ae_onsite,ur_ps_onsite)

     if (which_fixed==1) then ! electron
       do rr=1,nr_tot
         wp_idx = rcl2fft(rr)
         r12 = eh_red - rclred(:,rr)
         k_dot_r12 = two_pi * DOT_PRODUCT(k_bz,r12)
         eikr12 = DCMPLX(COS(k_dot_r12), SIN(k_dot_r12))
         exc_phi(rr) = exc_phi(rr) + eikr12 * ur_v(wp_idx) * CONJG(ur_c(wp_idx))
       end do
     else ! hole
       MSG_ERROR("Not coded")
     end if
     !
  end do
 end do

 ABI_DEALLOCATE(vec_list)
 ABI_DEALLOCATE(rcl2fft)
 ABI_DEALLOCATE(rclred)
 ABI_DEALLOCATE(ur_v)
 ABI_DEALLOCATE(ur_c)
 !
 if (Wfd%usepaw==1) then
   call cprj_free(Cp_v)
   ABI_DEALLOCATE(Cp_v)
   call cprj_free(Cp_c)
   ABI_DEALLOCATE(Cp_c)
   if (paw_add_onsite) then
     call pawfgrtab_free(Pawfgrtab)
     ABI_DEALLOCATE(Pawfgrtab)
     call destroy_paw_pwaves_lmn(Paw_onsite)
     ABI_DEALLOCATE(Paw_onsite)
   end if
 end if
 !
 ! Collect results on master node.
 call xsum_master(exc_phi,master,comm,ierr)
 !
 ! The exciton density probability.
 ABI_ALLOCATE(exc_phi2,(nr_tot))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0," out of memory in exc_phi2")

 do rr=1,nr_tot
   exc_phi2(rr) = ABS(exc_phi(rr))**2
 end do
 !
 ! Construct the supercell.
 sc_natom = ncells*Cryst%natom
 ABI_ALLOCATE(sc_xcart,(3,sc_natom))
 ABI_ALLOCATE(sc_typat,(sc_natom))

 ii=0
 do ir3=0,nrcl(3)-1  ! Loop over the replicaed cells.
   do ir2=0,nrcl(2)-1
     do ir1=0,nrcl(1)-1  
      !
      sc_rprimd(:,1) = ir1 * Cryst%rprimd(:,1) ! Workspace.
      sc_rprimd(:,2) = ir2 * Cryst%rprimd(:,2)
      sc_rprimd(:,3) = ir3 * Cryst%rprimd(:,3)
      scart_shift = SUM(sc_rprimd,DIM=2)
      do iatom=1,Cryst%natom
        ii=ii+1
        sc_xcart(:,ii) = Cryst%xcart(:,iatom) + scart_shift
        sc_typat(ii)   = Cryst%typat(iatom)
      end do
      !
    end do
   end do
 end do
 !
 ! Lattice vectors of the big box.
 do ii=1,3
   sc_rprimd(:,ii) = nrcl(ii) * Cryst%rprimd(:,ii)
 end do
 !
 ! Master writes the data.
 if (my_rank==master) then
   out_fname = TRIM(BS_files%out_basename)//"_EXC_WF"

   xsf_unt = get_unit()
   open(xsf_unt,file=out_fname,form='formatted',iostat=ios)
   ABI_CHECK(ios==0," Opening file: "//TRIM(out_fname))

   call printxsf(bbox(1),bbox(2),bbox(3),exc_phi2,sc_rprimd,origin0,sc_natom,Cryst%ntypat,sc_typat,&
&    sc_xcart,Cryst%znucl,xsf_unt,0)

   close(xsf_unt)
 end if

 ABI_DEALLOCATE(sc_xcart)
 ABI_DEALLOCATE(sc_typat)

 ABI_DEALLOCATE(exc_phi)
 ABI_DEALLOCATE(exc_phi2)
 !call xbarrier_mpi(comm)

end subroutine exc_plot
!!***
