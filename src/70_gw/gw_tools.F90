!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_transitions
!! NAME
!! make_transitions
!!
!! FUNCTION
!!  Calculate transition energies entering the espression for the irreducible polarizability.
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nsspol=1 for spin unpolarized, 2 for spin polarized calculations
!!  nbnds=total number of bands
!!  kmesh<bz_mesh_type>=datatype gathering info on the k-mesh:
!!   | %nbz=number of k-points in the full BZ
!!   | %nibz=number of k-points in the IBZ
!!   | %tab(nkbz)=table giving for each k-point in the BZ, the corresponding irreducible point in the IBZ array
!!   | %bz(3,nkbz)=reduced coordinated of k-points
!!  TOL_DELTA_OCC=tolerance on the difference of the occupation numbers
!!  gw_energy(nbnds,kmesh%nkibz,nsppol)=quasi-particle energies energies 
!!  occ(nbnds,kmesh%nkibz,nsppol)=occupation numbers
!!  chi0alg=integer defining the method used to calculate chi0
!!   0 ==> calculate chi0 using the Adler-Wiser expression
!!   1 ==> use spectral method 
!!  timrev=if 2, time-reversal symmetry is considered; 1 otherwise
!!
!! OUTPUT
!! my_max_rest,my_min_rest=Maximum and minimum resonant (posite) transition energy.
!! max_rest,min_rest=Maximun and minimum resonant (posite) transition energy treated by this node.
!!
!! PARENTS
!!      cchi0,cchi0q0
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine make_transitions(Wfd,chi0alg,nbnds,nbvw,nsppol,symchi,timrev,TOL_DELTA_OCC,&
& max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,gw_energy,occ,qpoint,bbp_ks_distrb)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_header,  only : hdr_clean
 use m_bz_mesh, only : bz_mesh_type, has_BZ_item, little_group
 use m_wfs,     only : wfs_descriptor

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'make_transitions'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: chi0alg,nbnds,nbvw,nsppol,symchi,timrev
 real(dp),intent(in) :: TOL_DELTA_OCC
 real(dp),intent(out) :: max_rest,min_rest
 real(dp),intent(out) :: my_max_rest,my_min_rest
 type(bz_mesh_type),intent(in) :: Kmesh
 type(little_group),intent(in) :: Ltg_q
 type(wfs_descriptor),intent(in) :: Wfd
!arrays
 real(dp),intent(in) :: gw_energy(nbnds,Kmesh%nibz,nsppol)
 real(dp),intent(in) :: occ(nbnds,Kmesh%nibz,nsppol),qpoint(3)
 integer,intent(in) :: bbp_ks_distrb(Wfd%mband,Wfd%mband,Kmesh%nbz,Wfd%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2,ii,ik_bz,ik_ibz,ikmq_bz,ikmq_ibz,is,nt,ntrans,my_ntrans,iloop
 real(dp) :: delta_ene,delta_occ,spin_fact
 character(len=500) :: msg
!arrays
 integer :: G0(3)
 real(dp) :: kmq(3)

!************************************************************************

 DBG_ENTER("COLL")

 if (chi0alg<0 .or. chi0alg>=2) then 
   write(msg,'(a,i3,a)')' chi0alg = ',chi0alg,' not allowed '
   MSG_BUG(msg)
 end if 
 if (timrev/=1 .and. timrev/=2) then 
   write(msg,'(a,i3,a)')' timrev = ',timrev,' not allowed'
   MSG_BUG(msg)
 end if 

 ABI_UNUSED(nbvw)
 !
 ! In the first loop calculate total number of transitions for this q-point
 ! as well min and max transition without taking into account distribution of bands. 
 ! In the second iteration calculate min and Max transition for this processor.
 !
 spin_fact=half; if (nsppol==2) spin_fact=one
 my_max_rest=smallest_real; my_min_rest=greatest_real
    max_rest=smallest_real;    min_rest=greatest_real

 do iloop=1,2
   nt=0
   do ik_bz=1,Kmesh%nbz
     ik_ibz=Kmesh%tab(ik_bz)
     kmq(:)=Kmesh%bz(:,ik_bz)-qpoint(:)

     if (symchi==1) then  
       if (Ltg_q%ibzq(ik_bz)/=1) cycle ! This point does not belong to the IBZ defined by the little group
     end if 
     !
     ! Find kp=k-q-G0 and also G0 where kp is in the first BZ
     if (.not.has_BZ_item(Kmesh,kmq,ikmq_bz,g0)) then ! Stop as the weight 1.0/nkbz is wrong.
       write(msg,'(4a,2(2a,3f12.6),2a)')ch10,&
&        ' make_transitions : ERROR - ',ch10,&
&        ' kp  = k-q-G0 not found in the BZ mesh',ch10,&
&        ' k   = ',(Kmesh%bz(ii,ik_bz),ii=1,3),ch10,&
&        ' k-q = ',(kmq(ii),ii=1,3),ch10,&
&        ' weight in cchi0/cchi0q is wrong ' 
       MSG_ERROR(msg)
     end if 

     ikmq_ibz=Kmesh%tab(ikmq_bz)
     do is=1,nsppol
       do ib1=1,nbnds           
         do ib2=1,nbnds

           if (iloop==2) then
             if (bbp_ks_distrb(ib1,ib2,ik_bz,is)/=Wfd%my_rank) cycle
           end if

           if (timrev==2 .and. ib1<ib2) cycle ! Thanks to time-reversal we gain a factor ~2.

           delta_occ=spin_fact*(occ(ib1,ikmq_ibz,is)-occ(ib2,ik_ibz,is))
           delta_ene=gw_energy(ib1,ikmq_ibz,is)-gw_energy(ib2,ik_ibz,is)

           if (chi0alg==0)  then ! Adler-Wiser expression. Skip only if factor due to occupation number is smaller than TOL_DELTA_OCC
             if (abs(delta_occ) < abs(TOL_DELTA_OCC)) cycle
           else if (chi0alg==1) then
             ! Spectral method with time-reversal, only resonant transitions 
             ! This has to changed to include spectral method without time-reversal
             if (delta_ene < -abs(TOL_DELTA_OCC) .or. abs(delta_occ) < abs(TOL_DELTA_OCC)) cycle
           end if
           !
           ! We have a new transition
           nt=nt+1

           if (iloop==1) then 
             max_rest=MAX(max_rest,zero,delta_ene)
             if (delta_ene>=-tol6) min_rest=MIN(min_rest,delta_ene)
           end if
           if (iloop==2) then 
             my_max_rest=MAX(my_max_rest,zero,delta_ene)
             if (delta_ene>=-tol6) my_min_rest=MIN(my_min_rest,delta_ene)
           end if

         end do 
       end do 
     end do
   end do
   if (iloop==1) ntrans=nt
   if (iloop==2) my_ntrans=nt
 end do !iloop

 write(msg,'(2a,i9,2a,f8.3,3a,f8.3,a)')ch10,&
&  ' Total number of transitions = ',ntrans,ch10,&
&  ' min resonant     = ',min_rest*Ha_eV,' [eV] ',ch10,&
&  ' Max resonant     = ',max_rest*Ha_eV,' [eV] '
 call wrtout(std_out,msg,'COLL')

 if (Wfd%nproc/=1) then 
   write(msg,'(2a,i9,2a,f8.3,3a,f8.3,a)')ch10,&
&    ' Total number of transitions for this processor= ',my_ntrans,ch10,&
&    ' min resonant     = ',my_min_rest*Ha_eV,' [eV] ',ch10,&
&    ' Max resonant     = ',my_max_rest*Ha_eV,' [eV] '
   call wrtout(std_out,msg,'PERS')
 end if

 DBG_EXIT("COLL")

end subroutine make_transitions 
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/get_rhor
!! NAME
!! get_rhor
!!
!! FUNCTION
!!  Read the DEN file with name fname reporting the density on the real FFT mesh 
!!  specified through the input variable ngfft_asked. If the FFT mesh asked and that found 
!!  on file differ, perform a FFT interpolation and renormalize the density so that it 
!!  integrates to the correct number of electrons. If the two FFT meshes coincides 
!!  just report the array stored on file.
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! fname=Name of the density file
!! accesswff=File format (Fortran, NETCDF-ETSF)
!! ngfft_asked(18)=All info on the FFT mesh required.
!! paral_kgb=Flag related to the kpoint-band-fft parallelism (FIXME not implemented)
!! MPI_enreg<MPI_type>=Information about MPI parallelization
!!  
!! OUTPUT
!! rhor_out(nfft_asked,nspden)=The interpolated density.
!!
!! NOTES
!!  1) The renormalization of the charge is not done in case of PAW since the onsite
!!     contributions have to be added. This is left to the caller.
!!  2) No check is done on sensible dimensions (nsppol, nspden). 
!!     TODO Maybe I can pass the header?
!!  3) Wavelet case not treated.
!!  4) Complex densities not treated.
!!
!! PARENTS
!!      mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_rhor(fname,accesswff,nspden,nfft_asked,ngfft_asked,paral_kgb,MPI_enreg,rhor_out,get_pawden)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_errors

 use m_io_tools, only : get_unit
 use m_header,   only : hdr_get_nelect_byocc, hdr_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_rhor'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_66_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: accesswff,paral_kgb,nfft_asked,nspden
 character(len=fnlen),intent(in) :: fname
 logical,intent(in), optional :: get_pawden
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 real(dp),intent(out) :: rhor_out(nfft_asked,nspden)
 integer,intent(in) :: ngfft_asked(18)

!Local variables-------------------------------
!scalars
 integer :: unt,ios,rdwr,fform,rdwrpaw,nfft_found,accessfil
 integer :: cplex,optin,optout
 real(dp) :: etotal,nelect,ratio,ucvol 
 logical :: do_interpolation,ltest
 type(Hdr_type) :: Hdr
 type(Dataset_type) :: Dtset
 type(wvl_internal_type) :: wvl
 character(len=500) :: msg
!arrays
 integer :: ngfft_found(18)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp) :: rhogdum(1)
 real(dp),allocatable :: rhor_found(:,:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)

! *************************************************************************

 accessfil=0
 if (accesswff==IO_MODE_NETCDF) accessfil=1
 if (accesswff==IO_MODE_ETSF)   accessfil=3
 if (accesswff==IO_MODE_MPI)    accessfil=4

 ltest=(nfft_asked==PRODUCT(ngfft_asked(1:3)))
 ABI_CHECK(ltest,'FFT parallelism not implemented')

 unt=get_unit()
 open(unit=unt,file=fname,form='unformatted',status='old',iostat=ios)
 ABI_CHECK(ios==0,'Opening file: '//TRIM(fname)//' as old')
 rdwr=1
 call hdr_io(fform,Hdr,rdwr,unt)
 close(unt)
 
 if (fform/=52) then
  write(msg,'(3a,i4)')' File ',TRIM(fname),' is not a density file: fform= ',fform
  MSG_ERROR(msg)
 end if

 ngfft_found(1:3) =Hdr%ngfft(1:3)  
 ngfft_found(4)   =2*(ngfft_found(1)/2)+1 ! 4:18 are not used in ioarr 
 ngfft_found(5)   =2*(ngfft_found(2)/2)+1 ! but they are used in fourdp
 ngfft_found(6)   =ngfft_found(3)
 ngfft_found(7:18)=ngfft_asked(7:18)
 nfft_found=PRODUCT(Hdr%ngfft(1:3))    ! TODO complex densities are not treated.

 do_interpolation=ANY(ngfft_found(1:3)/=ngfft_asked(1:3))

 rdwr=1; rdwrpaw=Hdr%usepaw !; if(dtfil%ireadwf/=0) rdwrpaw=0
 if (present(get_pawden)) then
   if (get_pawden) rdwrpaw=0
 end if
 ABI_ALLOCATE(Pawrhoij,(Hdr%natom*Hdr%usepaw*rdwrpaw))

 !We hack a bit Dtset, wavelets are not treated.
 Dtset%nspden=Hdr%nspden 
 Dtset%usewvl=0

 ABI_CHECK((Hdr%nspden==nspden),'Mismatch in nspden')

 if (.not.do_interpolation) then
   call ioarr(accessfil,rhor_out,Dtset,etotal,fform,fname,Hdr,MPI_enreg,nfft_found,Pawrhoij,rdwr,rdwrpaw,wvl)
 else 
   ABI_ALLOCATE(rhor_found,(nfft_found,nspden))
   call ioarr(accessfil,rhor_found,Dtset,etotal,fform,fname,Hdr,MPI_enreg,nfft_found,Pawrhoij,rdwr,rdwrpaw,wvl)

   ! * Get number of electrons, used to renormalize rhor after the interpolation.
   ! * Cannot use znucl because we might have additional charge or alchemy.
   nelect = hdr_get_nelect_byocc(Hdr)

   call wrtout(std_out," get_rhor: FFT meshes differ, performing Fourier interpolation.","COLL")

   cplex =1 
   optin =0 ! Input is taken from rhor
   optout=0 ! Output is only in real space but we might add an option to return rhog_out

   call fourier_interpol(cplex,nspden,optin,optout,nfft_found,ngfft_found,nfft_asked,ngfft_asked,&
&   paral_kgb,MPI_enreg,rhor_found,rhor_out,rhogdum,rhogdum)

   ABI_DEALLOCATE(rhor_found)

   ! === Renormalize charge to avoid errors due to the interpolation ===
   ! * Do this only for NC since for PAW we should add the onsite contribution.
   if (Hdr%usepaw==0) then
     call metric(gmet,gprimd,-1,rmet,Hdr%rprimd,ucvol)
     ratio=nelect/( SUM(rhor_out(:,1))*ucvol/PRODUCT(ngfft_asked(1:3)) )
     rhor_out(:,:)=rhor_out(:,:)*ratio
     write(msg,'(a,f8.2,a,f8.4)')' Expected nelect: ',nelect,' renormalization ratio: ',ratio
     call wrtout(std_out,msg,'COLL')
   end if
 end if ! Did ngfft_asked agree with ngfft_found?

 ! === Free memory ===
 call hdr_clean(Hdr)

 ! TODO here possible memory leak or crash as we do not know if rhoij is Initialized
 !if (rdwrpaw==1.and.restartpaw/=0) call rhoij_free(Pawrhoij) 
 if (rdwrpaw==1) call rhoij_free(Pawrhoij)
 ABI_DEALLOCATE(Pawrhoij)

end subroutine get_rhor
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/sigma_distribution
!! NAME
!!  sigma_distribution
!!
!! FUNCTION
!!  Distribute the loop over (b,k,s) used to calculate the self-energy matrix elements 
!!  taking into account the MPI distribution of the wavefunctions and the use of 
!!  symmetries to reduce the BZ sum to an appropriate irreducible wedge.
!!
!! INPUTS
!! nsppol
!! can_symmetrize(nsppol)=.TRUE if symmetries can be used to reduce the number of k-points to be summed.
!! Kmesh<bz_mesh_type>
!! Qmesh<bz_mesh_type>
!! Ltg_kgw<little_group>
!! Wfd(wfs_descriptor),intent(inout) :: 
!! mg0(3)
!! kptgw(3)
!! [bks_mask(Wfd%mband,Kmesh%nbz,nsppol)]
!! [got(Wfd%nproc)]=The number of tasks already assigned to the nodes.
!! [global]=If true, an MPI global communication is performed such that each node will have the same table. Useful
!!   if for implementing algorithms in which each node needs to know the global distribution of the tasks, not only
!!   the task it has to complete. Defaults to .FALSE.
!!
!! OUTPUT
!!  my_nbks
!!  proc_distrb(Wfd%mband,Kmesh%nbz,nsppol)
!!
!! SIDE EFFECTS
!!  Wfd%bks_tab
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cohsex_me
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine sigma_distribution(Wfd,Kmesh,Ltg_kgw,Qmesh,nsppol,can_symmetrize,kptgw,mg0,my_nbks,proc_distrb,got,bks_mask,global)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

 use m_bz_mesh,  only : bz_mesh_type, little_group, findqg0
 use m_wfs,      only : wfs_descriptor, wfd_distribute_bands, wfd_update_bkstab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma_distribution'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsppol
 integer,intent(out) :: my_nbks
 logical,optional,intent(in) :: global
 type(bz_mesh_type),intent(in) :: Kmesh,Qmesh
 type(little_group),intent(in) :: Ltg_kgw
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: mg0(3)
 integer,optional,intent(inout) :: got(Wfd%nproc)
 integer,intent(out) :: proc_distrb(Wfd%mband,Kmesh%nbz,nsppol)
 real(dp),intent(in) :: kptgw(3)
 logical,intent(in) :: can_symmetrize(Wfd%nsppol)
 logical,optional,intent(in) :: bks_mask(Wfd%mband,Kmesh%nbz,nsppol)

!Local variables-------------------------------
!scalars
 integer :: ierr,ik_bz,ik_ibz,spin,iq_bz,my_nband
 !character(len=500) :: msg
!arrays
 integer :: g0(3)
 real(dp) :: kgwmk(3)
 integer :: get_more(Wfd%nproc),my_band_list(Wfd%mband)
 !integer :: test(Wfd%mband,Kmesh%nbz,nsppol)
 logical :: bmask(Wfd%mband)

!************************************************************************

 call wfd_update_bkstab(Wfd)

 get_more=0; if (PRESENT(got)) get_more=got

 ! Different distribution of tasks depending whether symmetries can be used or not.
 proc_distrb= xmpi_undefined_rank

 do spin=1,Wfd%nsppol

   if (can_symmetrize(spin)) then 
     do ik_bz=1,Kmesh%nbz
       ik_ibz = Kmesh%tab(ik_bz)
       kgwmk= kptgw-Kmesh%bz(:,ik_bz) ! kptgw must be inside the BZ
       call findqg0(iq_bz,g0,kgwmk,Qmesh%nbz,Qmesh%bz,mG0) ! Identify q_bz and G0 where q_bz+G0=k_gw-k_bz
       if (Ltg_kgw%ibzq(iq_bz)==1) then
         bmask=.FALSE.; bmask(1:Wfd%nband(ik_ibz,spin))=.TRUE.
         if (PRESENT(bks_mask)) bmask = bks_mask(:,ik_bz,spin)
         call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,got=get_more,bmask=bmask)
         if (my_nband>0) proc_distrb(my_band_list(1:my_nband),ik_bz,spin)=Wfd%my_rank
       end if
     end do

   else ! No symmetries for this spin. Divide the full BZ among procs.
     do ik_bz=1,Kmesh%nbz
       ik_ibz = Kmesh%tab(ik_bz)
       bmask=.FALSE.; bmask(1:Wfd%nband(ik_ibz,spin))=.TRUE.
       if (PRESENT(bks_mask)) bmask = bks_mask(:,ik_bz,spin)
       call wfd_distribute_bands(Wfd,ik_ibz,spin,my_nband,my_band_list,got=get_more,bmask=bmask)
       if (my_nband>0) proc_distrb(my_band_list(1:my_nband),ik_bz,spin)=Wfd%my_rank
     end do
   end if
 end do ! spin

 if (PRESENT(global)) then 
   if (global) then ! Each node will have the same table so that it will know how the tasks are distributed.
     proc_distrb = proc_distrb + 1
     where (proc_distrb == xmpi_undefined_rank+1)
       proc_distrb = 0
     end where
     call xsum_mpi(proc_distrb,Wfd%comm,ierr)
     where (proc_distrb == 0)
       proc_distrb = xmpi_undefined_rank
     elsewhere
       proc_distrb = proc_distrb -1
     end where
     !where (proc_distrb /= xmpi_undefined_rank) 
     !  ltest = (ANY(proc_distrb == (/(ii,ii=0,Wfd%nproc-1)/)))
     !end where
     !if (.not.ltest) then
     !  write(std_out,*)proc_distrb
     !  MSG_BUG("Bug in the generation of proc_distrb table")
     !end if
   end if
 end if

 my_nbks = COUNT(proc_distrb==Wfd%my_rank)
 if (PRESENT(got)) got=get_more

end subroutine sigma_distribution
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/chi0_bbp_mask
!! NAME
!!  chi0_bbp_mask
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      cchi0,cchi0q0,check_completeness
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chi0_bbp_mask(Ep,use_tr,QP_BSt,mband,ikmq_ibz,ik_ibz,spin,spin_fact,bbp_mask)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

 use m_gwdefs,        only : GW_TOL_DOCC, epsilonm1_parameters

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chi0_bbp_mask'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin,ik_ibz,ikmq_ibz,mband
 real(dp),intent(in) :: spin_fact
 logical,intent(in) :: use_tr
 type(epsilonm1_parameters),intent(in) :: Ep
 type(Bandstructure_type),intent(in) :: QP_BSt
!arrays
 logical,intent(out) :: bbp_mask(mband,mband)

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2
 real(dp) :: deltaeGW_b1kmq_b2k,deltaf_b1kmq_b2k,e_b1_kmq,f_b1_kmq 
 character(len=500) :: msg

!arrays
 real(dp),pointer :: qp_energy(:,:,:),qp_occ(:,:,:)

!************************************************************************

 qp_energy => QP_BSt%eig(:,:,:)
 qp_occ    => QP_BSt%occ(:,:,:)

 bbp_mask=.FALSE.
 !use_tr = (Ep%awtr==1)

 SELECT CASE (Ep%gwcomp) 

 CASE (0)
 
 do ib1=1,Ep%nbnds ! Loop over "conduction" states.
   e_b1_kmq = qp_energy(ib1,ikmq_ibz,spin)
   f_b1_kmq =    qp_occ(ib1,ikmq_ibz,spin)

   do ib2=1,Ep%nbnds ! Loop over "valence" states.
     deltaf_b1kmq_b2k   = spin_fact*(f_b1_kmq-qp_occ(ib2,ik_ibz,spin))
     deltaeGW_b1kmq_b2k = e_b1_kmq-qp_energy(ib2,ik_ibz,spin)

     SELECT CASE (Ep%spmeth)
     CASE (0) ! Standard Adler-Wiser expression.

       if (ABS(deltaf_b1kmq_b2k) >= GW_TOL_DOCC) then 
         bbp_mask(ib1,ib2)=.TRUE.
         if (use_tr .and. ib1<ib2) bbp_mask(ib1,ib2)=.FALSE. ! GAIN a factor ~2 thanks to time-reversal.
       end if

     CASE (1,2) ! Spectral method, WARNING time-reversal here is always assumed!
       if (ABS(deltaf_b1kmq_b2k) >= GW_TOL_DOCC) then 
         bbp_mask(ib1,ib2)=.TRUE.
         if (deltaeGW_b1kmq_b2k<zero) bbp_mask(ib1,ib2)=.FALSE.  ! Only positive frequencies are needed for the Hilbert transform.
         !$if (use_tr .and. ib1<ib2) bbp_mask(ib1,ib2)=.FALSE. ! GAIN a factor ~2 thanks to time-reversal.
       end if

     CASE  DEFAULT
       write(msg,'(a,i0)')" Wrong value for spmeth: ",Ep%spmeth
       MSG_ERROR(msg)
     END SELECT

   end do !ib2
 end do !ib1

 CASE (1) ! Extrapolar technique
   ABI_CHECK(Ep%spmeth==0,"Hilbert transform and extrapolar method are not compatible")

   do ib1=1,Ep%nbnds ! Loop over "conduction" states.
     e_b1_kmq=qp_energy(ib1,ikmq_ibz,spin)
     f_b1_kmq=   qp_occ(ib1,ikmq_ibz,spin)

     do ib2=1,Ep%nbnds ! Loop over "valence" states.

       deltaf_b1kmq_b2k  =spin_fact*(f_b1_kmq-qp_occ(ib2,ik_ibz,spin))
       deltaeGW_b1kmq_b2k=e_b1_kmq-qp_energy(ib2,ik_ibz,spin)

       ! * When the completeness correction is used,
       !   we need to also consider transitions with vanishing deltaf
       !
       ! Rangel: This is to compute chi in metals correctly with the extrapolar method.
       bbp_mask(ib1,ib2)=.TRUE.
       !if (qp_occ(ib2,ik_ibz,is) < GW_TOL_DOCC) CYCLE
       if (qp_occ(ib2,ik_ibz,spin) < GW_TOL_DOCC .and. (ABS(deltaf_b1kmq_b2k) < GW_TOL_DOCC .or. ib1<ib2)) then
         bbp_mask(ib1,ib2)=.FALSE.
       end if

     end do 
   end do

  CASE DEFAULT
    write(msg,'(a,i0)')" Wrong value of gwcomp: ",Ep%gwcomp 
    MSG_ERROR(msg)
  END SELECT

end subroutine chi0_bbp_mask
!!***
