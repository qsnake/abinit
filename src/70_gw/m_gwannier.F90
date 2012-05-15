!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gwannier
!! NAME
!! m_gwannier
!!
!! FUNCTION
!!  This module contains procedures used to interpolate quasi-particle corrections
!!  using the Wannier representations previously obtained with Wannier90.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gwannier

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_wannier2abinit

 use m_io_tools,       only : prompt
 use m_numeric_tools,  only : set2unit
 use m_header,         only : hdr_clean
 use m_crystal,        only : destroy_crystal, crystal_structure
 use m_bz_mesh,        only : bz_mesh_type, init_kmesh, has_BZ_item, destroy_BZ_mesh_type, isamek, &
&                             nullify_BZ_mesh, make_mesh, make_path
 use m_ebands,         only : SelectBands, ExpandBands, print_bandstructure, bst_plot_bands, get_dos,&
&                             copy_bandstructure, bstruct_init, bstruct_clean
 use m_sigma_results,  only : abi_etsf_get_QP, print_Sigma_perturbative, destroy_Sigma_results, sigma_results

 implicit none

 private

 public :: my_GWannier

CONTAINS  !=========================================================================================================================
!!***

!!****f* m_gwannier/my_GWannier
!! NAME
!! my_GWannier
!!
!! FUNCTION
!!  Main routine for Wannier interpolation of GW results.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!      abi_etsf_get_qp,bst_plot_bands,bstruct_clean,bstruct_init
!!      copy_bandstructure,destroy_bz_mesh_type,destroy_crystal
!!      destroy_sigma_results,destroywandata,get_dos,hdr_check,hdr_clean
!!      init_kmesh,make_mesh,make_path,makewannierhr,metric,nullify_bz_mesh
!!      print_bandstructure,print_sigma_perturbative,printwandata,prompt
!!      readwandata,set2unit,wanmatinterpol,wannierinterpol
!!
!! SOURCE

subroutine my_GWannier

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'my_GWannier'
 use interfaces_42_geometry
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter :: prtvol=0
 integer :: accesswff,ndiv_small,mwan,ib,ii,nkbounds,nkintp,jj,info
 integer :: ikibz,restartpaw,restart,ikcalc,itask,ierr
 integer :: bantot,isppol,iband,spad,nband_k,method,ikpt,kptopt
 real(dp) :: ucvol,broad,dosdeltae
 !real(dp) :: fixmom
 logical :: use_afm,use_tr,ishermitian,witheader,found
 character(len=50) :: task
 character(len=500) :: msg
 character(len=fnlen) :: fname_wan,fname_bands,fname_qps,fname_dos
 type(BZ_mesh_type) :: Kmesh,Kmesh4dos,Kpath
 type(Crystal_structure) :: Cryst
 type(WannierData) :: WData
 type(Bandstructure_type) :: QP_BSt,KS_BSt,QPtmp,QP4Wan,KS_intp,QP_intp
 type(Sigma_results) :: Sr
 type(Hdr_type) :: Hdr
!arrays
 integer :: kptrlatt(3,3),mp_dense(3)
 integer :: G0(3),kptrlatt4dos(3,3)
 integer,allocatable :: ndiv(:),kcalc2bz(:)
 integer,allocatable :: nband(:) 
 integer,pointer :: dummy(:),npwarr(:),istwfk(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),shift(3)
 real(dp),allocatable :: kbounds(:,:),wtk(:)
 real(dp),pointer :: kintp(:,:)
 real(dp),allocatable :: matrix_out(:,:,:,:),matrix_in(:,:,:,:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:)
 logical,allocatable :: bandselect(:,:,:)

! *************************************************************************

 ! === Read WAN file containing the unitary transformation ===
 call prompt(' Enter name of the WAN file: ',fname_wan)

 accesswff=IO_MODE_FORTRAN
 call ReadWanData(WData,fname_wan,accesswff)

 call PrintWanData(WData,prtvol=1)

 call metric(gmet,gprimd,-1,rmet,WData%Hdr%rprimd,ucvol)
 
 ! === Construct the Hamiltonian in the Wannier representation ===
 call prompt(' Enter kptrlatt used during the Wannierization: ',kptrlatt)

 call MakeWannierHR(WData,kptrlatt)

 ! === Read GW results in NETCDF format ===
 call prompt(' Enter name of the QPS file: ',fname_qps)

 call abi_etsf_get_QP(Sr,KS_BSt,Hdr,Cryst,fname_qps)

 ! This is just to do a check, the file format is wrong!
 call hdr_check(1002,1002,Hdr,WData%Hdr,'COLL',restart,restartpaw)

 ! === Copy the KS bands to QP_Bst thus initializing the object ===
 ! * Apply the GW corrections.
 ! FIXME this works only if full k-mesh.
 call copy_bandstructure(KS_Bst,QP_BSt)
 do isppol=1,QP_BSt%nsppol
  do ikpt=1,QP_BSt%nkpt
   do iband=Sr%minbnd(ikpt,isppol),Sr%maxbnd(ikpt,isppol)
    QP_BSt%eig(iband,ikpt,isppol)=QP_BSt%eig(iband,ikpt,isppol)+REAL(Sr%degw(iband,ikpt,isppol))
   end do
  end do
 end do

 ! === Initialize the K-mesh (the same as that used of wavefunctions ====
 kptopt=1
 call init_kmesh(Kmesh,Cryst,Hdr%nkpt,Hdr%kptns,kptopt)

 ! === Find index of GW points in full BZ ===
 ! TODO the indexing has to be changed, everything should be defined in the IBZ or BZ
 ABI_ALLOCATE(kcalc2bz,(Sr%nkptgw))

 do ikcalc=1,Sr%nkptgw                                                                   
  found = has_BZ_item(Kmesh,Sr%kptgw(:,ikcalc),kcalc2bz(ikcalc),G0)
  !found = has_IBZ_item(Kmesh,Sigp%kptgw(:,ikcalc),kcalc2bz,G0)
  if (.not.found) then 
   write(msg,'(a,3(f6.3,1x),a)')' k-point ',Sr%kptgw(:,ikcalc),' not in the set of kbz'
   MSG_ERROR(msg)
  end if
 end do

 do isppol=1,Sr%nsppol
  do ikcalc=1,Sr%nkptgw
   ikibz=Kmesh%tab(kcalc2bz(ikcalc)) ! Irred k-point for GW
   do iband=Sr%minbnd(ikcalc,isppol),Sr%maxbnd(ikcalc,isppol) 
    witheader=.FALSE. ; if (iband==Sr%minbnd(ikcalc,isppol)) witheader=.TRUE.
    call print_Sigma_perturbative(Sr,ikibz,iband,isppol,witheader=witheader)
   end do
  end do
 end do

 ! === Extract GW bands for Wannier ===
 ! Here I am assuming that nwan is equal to the GW bands
 ! better integration will be done in the following.
 if (ANY(Sr%minbnd/=Sr%minbnd(1,1))) STOP 'GW bands are not constant'
 if (ANY(Sr%maxbnd/=Sr%maxbnd(1,1))) STOP 'GW bands are not constant'
 if (ANY(Sr%minbnd(:,1)/=1)) STOP 'GW bands should start at 1'
 !if (ANY(Sr%maxbnd(:)/=WData%mwan)) STOP 'GW bands should be equal to mwan'

 ii=Sr%minbnd(1,1) ; jj=Sr%maxbnd(1,1)
 if ((jj-ii+1)/=WData%mwan) STOP 'GW and Wannier bands do not agree'

 ABI_ALLOCATE(bandselect,(QP_BSt%mband,QP_BSt%nkpt,QP_BSt%nsppol))
 bandselect=.FALSE. ; bandselect(ii:jj,:,:)=.TRUE.

 QPtmp = SelectBands(QP_BSt,bandselect=bandselect)
 ABI_DEALLOCATE(bandselect)

 call print_bandstructure(QPtmp)

 use_tr=(Cryst%timrev==2) 
 use_afm=Cryst%use_antiferro

 QP4Wan = ExpandBands(QPtmp,WData%nkpt,WData%Hdr%kptns,use_tr,use_afm,Cryst%nsym,Cryst%symrec,Cryst%symafm,info)

 if (info/=0) then
  msg=' GW and Wannier mesh do not agree, check k-meshes. '
  MSG_ERROR(msg)
 end if
 call bstruct_clean(QPtmp)

 ! === Start Wannier Interpolation ===
 call prompt(' Enter task, 1 for Bands, 2 for DOS: ',itask)
 task='BANDS' ; if (itask/=1) task='DOS'

 call nullify_BZ_mesh(Kpath)
 call nullify_BZ_mesh(Kmesh4dos)

 SELECT CASE (task)

 CASE ('BANDS')
  call prompt(' Enter number of boundaries for path    : ',nkbounds)
  ABI_ALLOCATE(kbounds,(3,nkbounds))
  call prompt(' Enter boundaries [r.l.u]               : ',kbounds)
  call prompt(' Enter the smallest number of divisions : ',ndiv_small)

  ! === Interpolate GW corrections ===
  ! * Make path in reciprocal space
  ABI_ALLOCATE(ndiv,(nkbounds-1))
  call make_path(nkbounds,kbounds,gmet,'G',ndiv_small,ndiv,nkintp,kintp)
  ABI_DEALLOCATE(kbounds)

  ABI_ALLOCATE(wtk,(nkintp))
  wtk=one/nkintp 

  Kpath%nibz=nkintp
  ABI_ALLOCATE(Kpath%ibz,(3,Kpath%nibz))
  ABI_ALLOCATE(Kpath%wt,(Kpath%nibz))
  Kpath%ibz=kintp
  Kpath%wt=one/Kpath%nibz

 ! === Interpolate KS bands first ===
 ! * All the data we need are already stored in WData.
 call WannierInterpol(WData,Kpath,KS_intp)
 fname_bands='KS_bands'
 call bst_plot_bands(KS_intp,gmet,fname_bands,ierr)

 CASE ('DOS')
  ! === Construct mesh in BZ, fold it back to IBZ and find weights ===
  call prompt(' Enter MP divisions for the dense k-mesh: ',mp_dense)
  call set2unit(kptrlatt4dos)
  kptrlatt4dos(1,1)=mp_dense(1)
  kptrlatt4dos(2,2)=mp_dense(2)
  kptrlatt4dos(3,3)=mp_dense(3)
  call prompt(' Enter the shift for the mesh (multiple shifts not supported) ',shift)

  kptopt=1
  call make_mesh(Kmesh4dos,Cryst,kptopt,kptrlatt4dos,1,shift)
  nkintp=Kmesh4dos%nibz 
  ABI_ALLOCATE(kintp,(3,nkintp))
  kintp=Kmesh4dos%ibz 

  ABI_ALLOCATE(wtk,(nkintp))
  wtk=Kmesh4dos%wt 

  ! === Interpolate KS bands on Kmesh4dos === 
  call WannierInterpol(WData,Kmesh4dos,KS_intp)
  fname_dos='KS_tetra'//'_DOS' ; method=2
  broad=0.01/Ha_eV ; dosdeltae=0.001/Ha_eV

  ! === Calculate DOS ===
  call get_dos(KS_intp,Kmesh4dos,method,fname_dos,broad,dosdeltae)
  
 CASE DEFAULT
  MSG_BUG('Wrong task.')
 END SELECT

 mwan=WData%mwan
 if (mwan/=QP4Wan%mband) stop 'mwan/=QP4Wan%mband'
 ABI_ALLOCATE(matrix_in ,(mwan,mwan,WData%nkpt,WData%nsppol))
 ABI_ALLOCATE(matrix_out,(mwan,mwan,nkintp    ,WData%nsppol))
 matrix_in=czero
 do ib=1,WData%mwan
  matrix_in(ib,ib,:,:)=QP4Wan%eig(ib,:,:)
 end do

 ishermitian=.TRUE.
 call WanMatInterpol(WData,QP4Wan%nkpt,nkintp,kintp,ishermitian,QP4Wan%nsppol,mwan,matrix_in,matrix_out)

 ! Now matrix_out has the interpolated band structure.
 !do ii=1,nkintp
 ! write(54,'(8f8.4)')(matrix_out(ib,ib,ii,1),ib=1,mwan) 
 ! write(55,'(8f8.4)')(KS_intp%eig(ib,ii,1),ib=1,mwan) 
 !end do

 ! === Initialize new object containing the interpolated GW bands ===
 bantot=QP4Wan%nsppol*nkintp*mwan
 ABI_ALLOCATE(doccde,(bantot))
 ABI_ALLOCATE(eigen,(bantot))
 ABI_ALLOCATE(occfact,(bantot))
 doccde(:)=zero ; eigen(:)=zero ; occfact(:)=zero 

 ABI_ALLOCATE(nband,(nkintp*QP4Wan%nsppol))
 do isppol=1,QP4Wan%nsppol
  spad=(isppol-1)*nkintp
  nband(spad+1:spad+nkintp)=WData%nwan(isppol)
 end do
 ABI_ALLOCATE(dummy,(nkintp))
 dummy=1
 istwfk => dummy
 npwarr => dummy

 call bstruct_init(bantot,QP_intp,QP4Wan%nelect,doccde,eigen,istwfk,kintp,nband,nkintp,npwarr,&
& QP4Wan%nsppol,QP4Wan%nspinor,QP4Wan%tphysel,QP4Wan%tsmear,QP4Wan%occopt,occfact,wtk)

 nullify(istwfk,npwarr)
 ABI_DEALLOCATE(dummy)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(occfact)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(wtk)

 ! === Put interpolated GW energies ===
 do isppol=1,QP_intp%nsppol
  do ikpt=1,QP_intp%nkpt
   nband_k=QP_intp%nband(ikpt+(isppol-1)*QP_intp%nkpt)
   do ib=1,nband_k
    QP_intp%eig(ib,ikpt,isppol)=matrix_out(ib,ib,ikpt,isppol)
   end do
  end do
 end do
 ABI_DEALLOCATE(matrix_out)
 ABI_DEALLOCATE(matrix_in)

 !this is dangerous
 !fixmom=99.99_dp
 !call update_occ(QP_intp,fixmom)

 if (task=='BANDS') then
  fname_bands='gw_bands'
  call bst_plot_bands(QP_intp,gmet,fname_bands,ierr)
 else if (task=='DOS') then
  broad=0.01/Ha_eV ; dosdeltae=0.001/Ha_eV
  fname_dos='GW_gauss'//'_DOS' ; method=1 
  !call get_dos(QP_intp,Kmesh4dos,method,fname_dos,broad,dosdeltae)
  fname_dos='GW_tetra'//'_DOS' ; method=2
  call get_dos(QP_intp,Kmesh4dos,method,fname_dos,broad,dosdeltae)
 end if

 ! this does not work yet
 !fixmom=99.99_dp
 !call update_occ(KS_intp,fixmom)
 !call ReportGap(KS_intp,'Ks Gaps',unit=std_out)

 !call update_occ(QP_intp,fixmom)
 !call ReportGap(QP_intp,'QP Gaps',unit=std_out)

 ! === Free memory ===
 ABI_DEALLOCATE(kcalc2bz)
 ABI_DEALLOCATE(kintp)

 call DestroyWanData(WData)
 call destroy_Sigma_results(Sr)
 call destroy_crystal(Cryst)
 call destroy_BZ_mesh_type(Kmesh)
 call destroy_BZ_mesh_type(Kmesh4dos)
 call destroy_BZ_mesh_type(Kpath)
 call bstruct_clean(KS_BSt) 
 call bstruct_clean(QP_BSt) 
 !$call bstruct_clean(QP_intp) 
 !$call bstruct_clean(KS_intp) 
 call hdr_clean(Hdr)

 STOP 'myGWannier OK'

end subroutine my_GWannier

END MODULE m_gwannier
!!***
