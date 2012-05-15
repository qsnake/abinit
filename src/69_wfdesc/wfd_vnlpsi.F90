!{\src2tex{textfont=tt}}
!!****f* ABINIT/wfd_vnlpsi
!! NAME
!! wfd_vnlpsi
!!
!! FUNCTION
!!  Evaluates Vnl |psi> in reciprocal space.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Wfd<wfs_descriptor>=Datatype gathering info on the wavefunctions
!! band=Band index.
!! ik_ibz=K-point index.
!! spin=Spin index
!! npw_k=Number of PW for this k-point (used to dimension output arrays)
!! Cryst<Crystal_structure>= data type gathering info on symmetries and unit cell
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! GS_hamk<gs_hamiltonian_type)=Information about the Hamiltonian for this k-point.
!!  Note that no check is done for the consistency btw ik_ibz and the content of GS_hamk.
!! [Kext]<Kdata_t>=Datatype storing form factors and tables that depend of the K.
!!   In the standard mode the routine uses the tables stored in Wfd. 
!!   Kext is used to calculate Vnl e^{ik'r}|u_k> where k is defined by ik_ibz whereas k' is 
!!   the k-point that has been used to initialize GS_hamk and Kext. This option is used for the
!!   Bloch-state-based interpolation.
!!
!! OUTPUT 
!!  vnl_psi(2,npw_k*Wfd%nspinor)= <G+k|V_nl e^{ik'r}|u_k> 
!!  opaw_psi(2,npw_k*Wfd%nspinor*Wfd%usepaw)  <G+k|1+S|Cnk> e^{ik'r}|u_k> 
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,mkkpg,nonlop,wfd_get_cprj
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wfd_vnlpsi(Wfd,band,ik_ibz,spin,npw_k,Cryst,Psps,GS_hamk,vnl_psi,opaw_psi,Kext)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

 use m_crystal,        only : crystal_structure
 use m_wfs,            only : wfs_descriptor, kdata_t, wfd_get_cprj

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_vnlpsi'
 use interfaces_44_abitypes_defs
 use interfaces_65_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin,npw_k
 type(Crystal_structure),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
 type(gs_hamiltonian_type),intent(in) :: GS_hamk
 type(Kdata_t),optional,intent(in) :: Kext
 type(wfs_descriptor),intent(inout) :: Wfd
!arrays
 real(dp),intent(out) :: vnl_psi(2,npw_k*Wfd%nspinor)  ! <G|Vnl|Cnk>
 real(dp),intent(out) :: opaw_psi(2,npw_k*Wfd%nspinor*Wfd%usepaw) ! <G|1+S|Cnk>

!Local variables ------------------------------
!scalars
 integer,parameter :: idir0=0,ider0=0,nnlout0=0,tim_nonlop0=0
 integer :: istat,natom,nkpg,dimffnl,cp_dim
 integer :: choice,cpopt,matblk,paw_opt,signs,istwf_k
 integer :: dimenl1,dimenl2,nspinor
 real(dp),parameter :: lambda0=zero
 logical :: use_kptext,ltest
 character(len=500) :: msg
!arrays
 integer :: nloalg(5)
 integer,pointer :: kg_k(:,:)
 real(dp) :: kpoint(3)
 real(dp) :: dum_enlout(0)
 real(dp),pointer :: ffnl(:,:,:,:),kpg_k(:,:),ph3d(:,:,:)
 real(dp),allocatable :: vectin(:,:)
 type(cprj_type),allocatable :: Cprj(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")

 use_kptext = PRESENT(Kext)

 natom   = Cryst%natom
 nloalg  = Wfd%nloalg
 nspinor = Wfd%nspinor

 dimenl1 = GS_hamk%dimekb1
 dimenl2 = GS_hamk%dimekb2
 dimffnl=1+ider0 ! Derivatives are not needed. 

 !matblk=nloalg(4); if (nloalg(1)>0) matblk=natom
 matblk=natom

 !npw_k   =  Wfd%Kdata(ik_ibz)%npw
 kg_k    => Wfd%Kdata(ik_ibz)%kg_k

 if (.not.use_kptext) then
   istwf_k =  Wfd%Kdata(ik_ibz)%istwfk
   kpoint  =  Wfd%kibz(:,ik_ibz)
   ffnl => Wfd%Kdata(ik_ibz)%fnl_dir0der0
   ph3d => Wfd%Kdata(ik_ibz)%ph3d
 else 
   istwf_k =  Kext%istwfk
   kpoint  =  GS_hamk%kpoint
   ffnl    => Kext%fnl_dir0der0
   ph3d    => Kext%ph3d
   ltest = (istwf_k==1 .and. Wfd%istwfk(ik_ibz)==1)
   if (.not.ltest) then
     write(msg,'(2a,2(a,i0))')&
&      " istwfk and Wfd%istwfk(ik_ibz) must be 1 when Kext is present, however, ",ch10,&
&      " istwfk= ",istwf_k," and Wfd%istwfk= ",Wfd%istwfk(ik_ibz)
     MSG_BUG(msg)
   end if                                                                               
   ltest = (ik_ibz==1 .and. ALL(ABS(Wfd%kibz(:,ik_ibz))<tol6))
   if (.not.ltest) then
     write(msg,'(3a,i0,a,3(f8.3,1x))')&
&      " ik_ibz must be 1 and kibz should be 0 when Kext is present, however, ",ch10,&
&      " ik_ibz= ",ik_ibz," and kpoint= ",Wfd%kibz(:,ik_ibz)
     MSG_BUG(msg)
   end if
 end if
 !
 ! Input wavefunction coefficients <G|Cnk>
 ABI_ALLOCATE(vectin,(2,npw_k*nspinor))
 vectin(1,:) = DBLE (Wfd%Wave(band,ik_ibz,spin)%ug) 
 vectin(2,:) = AIMAG(Wfd%Wave(band,ik_ibz,spin)%ug)

 signs  = 2  ! => apply the non-local operator to a function in G-space.
 choice = 1  ! => <G|V_nonlocal|vectin>.
 cpopt  =-1
 paw_opt= 0
 if (Wfd%usepaw==1) then 
   paw_opt=4 ! both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
   cpopt=3   ! <p_lmn|in> are already in memory 
   if (use_kptext) cpopt=-1 ! Since we have to change the k-point <p_lmn|in> are recomputed in nonlocal and not saved 
 end if

 cp_dim = ((cpopt+5)/5)
 ABI_ALLOCATE(Cprj,(natom,nspinor*cp_dim))

 if (cp_dim>0) then
   call cprj_alloc(Cprj,0,Wfd%nlmn_sort)
   call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cprj,sorted=.TRUE.)
 end if
 !
 ! Compute (k+G) vectors 
 ! Not needed here, however they might be stored in Kdata_t to avoid the calculation inside the loop over bands
 !nkpg=3*nloalg(5)  
 nkpg=0 
 ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
 if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

 call nonlop(Cryst%atindx1,choice,cpopt,Cprj,dimenl1,dimenl2,dimffnl,dimffnl,&
&  GS_hamk%ekb,dum_enlout,ffnl,ffnl,Cryst%gmet,Cryst%gprimd,idir0,Psps%indlmn,istwf_k,&
&  kg_k,kg_k,kpg_k,kpg_k,kpoint,kpoint,lambda0,Psps%lmnmax,matblk,Wfd%mgfft,&
&  Wfd%MPI_enreg,Psps%mpsang,Psps%mpssoang,Cryst%natom,Cryst%nattyp,Wfd%ngfft,nkpg,nkpg,nloalg,&
&  nnlout0,npw_k,npw_k,nspinor,nspinor,Cryst%ntypat,0,paw_opt,GS_hamk%phkxred,GS_hamk%phkxred,&
&  Wfd%ph1d,ph3d,ph3d,signs,GS_hamk%sij,&
&  opaw_psi,tim_nonlop0,Cryst%ucvol,Psps%useylm,vectin,vnl_psi)

 ABI_DEALLOCATE(kpg_k)
 istat = ABI_ALLOC_STAT
 ABI_DEALLOCATE(vectin)

 if (cp_dim>0) call cprj_free(Cprj)
 ABI_DEALLOCATE(Cprj)
 istat = ABI_ALLOC_STAT

 DBG_EXIT("COLL")

end subroutine wfd_vnlpsi       
!!***
