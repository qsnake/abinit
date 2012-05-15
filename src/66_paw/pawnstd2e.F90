!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawnstd2e
!! NAME
!! pawnstd2e
!!
!! FUNCTION
!! This routine compute the PAW on-site contributions to non-stationary expression for the
!! second derivative of the total energy, for couple of mixed derivatives.
!!  These contributions are equal to:
!!    Int{ VHxc[n1^(j2);nc^(j2)].delta_n1^(j1) }
!!   -Int{ VHxc[tild_n1^(j2)+hat_n1^(j2);tild_n_c^(j2)].delta_(tild_n1+hat_n1)^(j1) }
!!  See PRB 78, 035105 (2008), Eq.(80)
!!  delta_n^(j) is the first order density only due to change of WF overlap
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ipert1,ipert2=indexes of perturbations
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms on current process, size of PAW arrays
!!  natom_tot=total number of atoms in cell
!!  ntypat=number of types of atoms in unit cell.
!!  nzlmopt1= For the (j1) perturbation:
!!            if -1, compute all LM-moments of 1st-order densities and use non-zero LM-moments
!!            if  0, compute all LM-moments of densities and use all LM-moments
!!            if +1, compute only non-zero LM-moments of 1st-order densities (stored before)
!!  nzlmopt2= For the (j2) perturbation:
!!            if -1, compute all LM-moments of 1st-order densities and use non-zero LM-moments
!!            if  0, compute all LM-moments of densities and use all LM-moments
!!            if +1, compute only non-zero LM-moments of 1st-order densities (stored before)
!!  paw_an0(natom) <type(paw_an_type)>=paw arrays for 0th-order quantities given on angular mesh
!!  paw_an1(natom) <type(paw_an_type)>=paw arrays for 1st-order quantities given on angular mesh
!!                                     This corresponds to (j1) perturbation
!!  paw_ij1(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!                                     This corresponds to (j1) perturbation
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij1(natom) <type(pawrhoij_type)>= paw rhoij 1st-order occupancies for the (j1) perturbation
!!  pawrhoij2(natom) <type(pawrhoij_type)>= paw rhoij 1st-order occupancies for the (j2) perturbation
!!                                          (note: this is delta_rhoij^(j2))
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  epawnst(2)= real and imaginary parts of contributions to non-stationary expression for the
!!              second derivative of the total energy
!!
!! SIDE EFFECTS
!!    ==== if paw_an1(:)%has_vxc<2, compute 1st-order XC potentials
!!      paw_an1(natom)%vxc1(cplex1*mesh_size,:,nspden) =AE 1st-order XC potential Vxc^(j1)
!!      paw_an1(natom)%vxct1(cplex1*mesh_size,:,nspden)=PS 1st-order XC potential tVxc^(j1)
!!    ==== if paw_ij1(:)%has_dijhartree<2, compute 1st-order Dij_hartree
!!      paw_ij1(natom)%dijhartree(cplex1*lmn2_size)=Hartree contribution to Dij^(j1)
!!
!! PARENTS
!!      nstpaw3
!!
!! CHILDREN
!!      pawdensities,pawdijhartree,pawxc3,pawxcm3,timab,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawnstd2e(epawnst,ipert1,ipert2,mpi_enreg,natom,natom_tot,ntypat,nzlmopt1,nzlmopt2,&
&                    paw_an0,paw_an1,paw_ij1,pawang,pawprtvol,pawrad,pawrhoij1,pawrhoij2,&
&                    pawtab,pawxcdev,xclevel)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawnstd2e'
 use interfaces_18_timing
 use interfaces_56_xc
 use interfaces_66_paw, except_this_one => pawnstd2e
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert1,ipert2,natom,natom_tot,ntypat,nzlmopt1,nzlmopt2,pawprtvol,pawxcdev,xclevel
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 real(dp),intent(out) :: epawnst(2)
 type(paw_an_type),intent(in) :: paw_an0(natom)
 type(paw_an_type),intent(inout) :: paw_an1(natom)
 type(paw_ij_type),intent(inout) :: paw_ij1(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij1(natom),pawrhoij2(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: cplex1,cplex2,cplex_dijh1,iatom,iatom_tot,ierr,irhoij,ispden,itypat,jrhoij,klmn
 integer :: lm_size1,lm_size2,mesh_size,nspden,nspdiag,opt_compch,optexc,optvxc
 integer :: usecore,usetcore,usexcnhat
 real(dp) :: compch,eexc,eexc_im
 character(len=500) :: msg
!arrays
 logical,allocatable :: lmselect1(:),lmselect2(:),lmselect_tmp(:)
 real(dp) :: dij(2),epawnst_h(2),epawnst_xc(2),ro(2),tsec(2)
 real(dp),allocatable :: kxc_dum(:,:,:),nhat1(:,:,:),rho1(:,:,:),trho1(:,:,:)


! *************************************************************************

 DBG_ENTER("COLL")

 call timab(567,1,tsec)

 if (.not.(ipert1==natom_tot+1.or.ipert1==natom_tot+5.or.&
& ipert2==natom_tot+1.or.ipert2==natom_tot+5)) then
   if((abs(nzlmopt1)/=1.and.nzlmopt1/=0).or.(abs(nzlmopt2)/=1.and.nzlmopt2/=0)) then
     msg='invalid value for nzlmopt !'
     MSG_BUG(msg)
   end if
   if(paw_ij1(1)%has_dijhartree==0) then
     msg='dijhartree must be allocated !'
     MSG_BUG(msg)
   end if
   if(paw_an1(1)%has_vxc==0) then
     msg='vxc1 and vxct1 must be allocated !'
     MSG_BUG(msg)
   end if
   if(paw_an0(1)%has_kxc==0) then
     msg='kxc1 must be allocated !'
     MSG_BUG(msg)
   end if
   if ((ipert1<=natom_tot.or.ipert1==natom_tot+1).and.paw_an0(1)%has_kxc/=2) then
     msg='XC kernels for ground state must be in memory !'
     MSG_BUG(msg)
   end if
   if (paw_ij1(1)%cplex/=paw_an1(1)%cplex) then
     msg='paw_ij1()%cplex and paw_an1()%cplex must be equal !'
     MSG_BUG(msg)
   end if
   if (pawrhoij1(1)%cplex<paw_an1(1)%cplex.or.pawrhoij2(1)%cplex<paw_an1(1)%cplex) then
     msg='pawrhoij()%cplex must be >=paw_an1()%cplex  !'
     MSG_BUG(msg)
   end if
   if (pawrhoij1(1)%nspden/=pawrhoij2(1)%nspden) then
     msg='pawrhoij1()%nspden must =pawrhoij2()%nspden  !'
     MSG_BUG(msg)
   end if
   if (mpi_enreg%nproc_atom>1) then
     if (natom/=mpi_enreg%natom) then
       msg='natom not equal to mpi_enreg%natom !'
       MSG_BUG(msg)
     end if
   end if
 end if

!Init contribution to 2nd-order energy
 epawnst(1:2)=zero

!For some perturbations, nothing else to do
 if (ipert1==natom_tot+1.or.ipert1==natom_tot+5.or.&
& ipert2==natom_tot+1.or.ipert2==natom_tot+5) then
   return
 end if

!Various inits
 opt_compch=0;optvxc=1;optexc=3
 usecore=0;usetcore=0  ! This is true for phonons and Efield pert.
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)
 epawnst_xc(1:2)=zero;epawnst_h(1:2)=zero
 dij(1:2)=zero;ro(1:2)=zero


!================ Loop on atomic sites =======================
 do iatom=1,natom
   iatom_tot=iatom;if (mpi_enreg%nproc_atom>1) iatom_tot=mpi_enreg%atom_indx(iatom)

   itypat=pawrhoij1(iatom)%itypat
   mesh_size=pawrad(itypat)%mesh_size
   nspden=pawrhoij1(iatom)%nspden
   cplex1=pawrhoij1(iatom)%cplex
   cplex2=pawrhoij2(iatom)%cplex
   cplex_dijh1=paw_ij1(iatom)%cplex
   lm_size1=paw_an1(iatom)%lm_size
   lm_size2=paw_an1(iatom)%lm_size

!  If Vxc potentials are not in memory, compute them
   if (paw_an1(iatom)%has_vxc/=2) then
     ABI_ALLOCATE(rho1 ,(cplex1*mesh_size,lm_size1,nspden))
     ABI_ALLOCATE(trho1,(cplex1*mesh_size,lm_size1,nspden))
     ABI_ALLOCATE(nhat1,(cplex1*mesh_size,lm_size1,nspden*usexcnhat))
     ABI_ALLOCATE(lmselect1,(lm_size1))
     lmselect1(:)=paw_an1(iatom)%lmselect(:)
     ABI_ALLOCATE(lmselect_tmp,(lm_size1))
     lmselect_tmp(:)=.true.
     if (nzlmopt1==1) lmselect_tmp(:)=lmselect1(:)
!    Compute on-site 1st-order densities
     call pawdensities(compch,cplex1,iatom_tot,lmselect_tmp,lmselect1,&
&     lm_size1,nhat1,nspden,nzlmopt1,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&     pawrad(itypat),pawrhoij1(iatom),pawtab(itypat),rho1,trho1)
     ABI_DEALLOCATE(lmselect_tmp)
!    Compute on-site 1st-order xc potentials
     if (pawxcdev/=0) then
       call pawxcm3(pawtab(itypat)%coredens,cplex1,cplex1,eexc,paw_an0(iatom)%kxc1,&
&       lm_size1,lmselect1,nhat1,paw_an0(iatom)%nkxc1,nspden,optvxc,&
&       pawang,pawrad(itypat),pawxcdev,rho1,usecore,0,&
&       paw_an1(iatom)%vxc1,xclevel)
       call pawxcm3(pawtab(itypat)%tcoredens,cplex1,cplex1,eexc,paw_an0(iatom)%kxct1,&
&       lm_size1,lmselect1,nhat1,paw_an0(iatom)%nkxc1,nspden,optvxc,&
&       pawang,pawrad(itypat),pawxcdev,trho1,usetcore,2*usexcnhat,&
&       paw_an1(iatom)%vxct1,xclevel)
     else
       call pawxc3(pawtab(itypat)%coredens,cplex1,cplex1,eexc,paw_an0(iatom)%kxc1,&
&       lm_size1,lmselect1,nhat1,paw_an0(iatom)%nkxc1,nspden,optvxc,&
&       pawang,pawrad(itypat),rho1,usecore,0,&
&       paw_an1(iatom)%vxc1,xclevel)
       call pawxc3(pawtab(itypat)%tcoredens,cplex1,cplex1,eexc,paw_an0(iatom)%kxct1,&
&       lm_size1,lmselect1,nhat1,paw_an0(iatom)%nkxc1,nspden,optvxc,&
&       pawang,pawrad(itypat),trho1,usetcore,2*usexcnhat,&
&       paw_an1(iatom)%vxct1,xclevel)
     end if
     paw_an1(iatom)%has_vxc=2
     ABI_DEALLOCATE(lmselect1)
     ABI_DEALLOCATE(rho1)
     ABI_DEALLOCATE(trho1)
     ABI_DEALLOCATE(nhat1)
   end if ! has_vxc

!  If Dij_hartree are not in memory, compute them
   if (paw_ij1(iatom)%has_dijhartree/=2) then
     call pawdijhartree(cplex_dijh1,iatom,natom,ntypat,paw_ij1,pawrhoij1,pawtab)
   end if

!  Compute contribution to 2nd-order energy from 1st-order XC potential
   ABI_ALLOCATE(rho1 ,(cplex2*mesh_size,lm_size2,nspden))
   ABI_ALLOCATE(trho1,(cplex2*mesh_size,lm_size2,nspden))
   ABI_ALLOCATE(nhat1,(cplex2*mesh_size,lm_size2,nspden*usexcnhat))
   ABI_ALLOCATE(lmselect2,(lm_size2))
   lmselect2(:)=paw_an1(iatom)%lmselect(:)
   ABI_ALLOCATE(lmselect_tmp,(lm_size2))
   lmselect_tmp(:)=.true.
   if (nzlmopt2==1) lmselect_tmp(:)=lmselect2(:)
!  Compute on-site 1st-order densities
   call pawdensities(compch,cplex2,iatom_tot,lmselect_tmp,lmselect2,&
&   lm_size2,nhat1,nspden,nzlmopt2,opt_compch,1-usexcnhat,-1,0,pawang,pawprtvol,&
&   pawrad(itypat),pawrhoij2(iatom),pawtab(itypat),rho1,trho1)
   ABI_DEALLOCATE(lmselect_tmp)
!  Compute contributions to 2nd-order energy
   if (pawxcdev/=0) then
     ABI_ALLOCATE(kxc_dum,(mesh_size,pawang%angl_size,0))
     call pawxcm3(pawtab(itypat)%coredens,cplex2,cplex1,eexc,kxc_dum,&
&     lm_size2,lmselect2,nhat1,0,nspden,optexc,pawang,pawrad(itypat),&
&     pawxcdev,rho1,usecore,0,paw_an1(iatom)%vxc1,xclevel,d2enxc_im=eexc_im)
     epawnst_xc(1)=epawnst_xc(1)+eexc
     epawnst_xc(2)=epawnst_xc(2)+eexc_im
     call pawxcm3(pawtab(itypat)%tcoredens,cplex2,cplex1,eexc,kxc_dum,&
&     lm_size2,lmselect2,nhat1,0,nspden,optexc,pawang,pawrad(itypat),&
&     pawxcdev,trho1,usetcore,2*usexcnhat,paw_an1(iatom)%vxct1,xclevel,&
&     d2enxc_im=eexc_im)
     ABI_DEALLOCATE(kxc_dum)
     epawnst_xc(1)=epawnst_xc(1)-eexc
     epawnst_xc(2)=epawnst_xc(2)-eexc_im
   else
     ABI_ALLOCATE(kxc_dum,(mesh_size,lm_size2,0))
     call pawxc3(pawtab(itypat)%coredens,cplex2,cplex1,eexc,kxc_dum,&
&     lm_size2,lmselect2,nhat1,0,nspden,optexc,pawang,pawrad(itypat),&
&     rho1,usecore,0,paw_an1(iatom)%vxc1,xclevel,d2enxc_im=eexc_im)
     epawnst_xc(1)=epawnst_xc(1)+eexc
     epawnst_xc(2)=epawnst_xc(2)+eexc_im
     call pawxc3(pawtab(itypat)%tcoredens,cplex2,cplex1,eexc,kxc_dum,&
&     lm_size2,lmselect2,nhat1,0,nspden,optexc,pawang,pawrad(itypat),&
&     trho1,usetcore,2*usexcnhat,paw_an1(iatom)%vxct1,xclevel,&
&     d2enxc_im=eexc_im)
     ABI_DEALLOCATE(kxc_dum)
     epawnst_xc(1)=epawnst_xc(1)-eexc
     epawnst_xc(2)=epawnst_xc(2)-eexc_im
   end if
   ABI_DEALLOCATE(lmselect2)
   ABI_DEALLOCATE(rho1)
   ABI_DEALLOCATE(trho1)
   ABI_DEALLOCATE(nhat1)

!  Compute contribution to 2nd-order energy from 1st-order Hartree potential
   nspdiag=1;if (nspden==2) nspdiag=2
   do ispden=1,nspdiag
     if (cplex_dijh1==1) then
       jrhoij=1
       do irhoij=1,pawrhoij2(iatom)%nrhoijsel
         klmn=pawrhoij2(iatom)%rhoijselect(irhoij)
         dij(1)=paw_ij1(iatom)%dijhartree(klmn)
         ro(1)=pawrhoij2(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         epawnst_h(1)=epawnst_h(1)+ro(1)*dij(1)
         if (cplex2==2) then
           ro(2)=pawrhoij2(iatom)%rhoijp(jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
           epawnst_h(2)=epawnst_h(2)+ro(2)*dij(1)
         end if
         jrhoij=jrhoij+cplex2
       end do
     else ! cplex_dijh1==2
       jrhoij=1
       do irhoij=1,pawrhoij2(iatom)%nrhoijsel
         klmn=pawrhoij2(iatom)%rhoijselect(irhoij)
         dij(1:2)=paw_ij1(iatom)%dijhartree(2*klmn-1:2*klmn)
         ro(1)=pawrhoij2(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         epawnst_h(1)=epawnst_h(1)+ro(1)*dij(1)
         epawnst_h(2)=epawnst_h(2)-ro(1)*dij(2)
         if (cplex2==2) then
           ro(2)=pawrhoij2(iatom)%rhoijp(jrhoij+1,ispden)*pawtab(itypat)%dltij(klmn)
           epawnst_h(1)=epawnst_h(1)+ro(2)*dij(2)
           epawnst_h(2)=epawnst_h(2)+ro(2)*dij(1)
         end if
         jrhoij=jrhoij+cplex2
       end do
     end if
   end do

!  ================ End loop oon atomic sites =======================
 end do

!Final building of 2nd-order non-stationnary energy
!See PRB 78, 035105 (2008), Eq.(80) (last two terms)
 epawnst(1:2)=epawnst_xc(1:2)+epawnst_h(1:2)

!Reduction in case of parallelism
 if (mpi_enreg%nproc_atom>1) then
   call timab(48,1,tsec)
   call xsum_mpi(epawnst,mpi_enreg%comm_atom,ierr)
   call timab(48,2,tsec)
 end if

 call timab(567,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawnstd2e
!!***
