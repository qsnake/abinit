!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_rpa_functional
!! NAME
!! calc_rpa_functional
!!
!! FUNCTION
!!  Routine used to calculate the RPA approximation to the correlation energy
!!  from the irreducible polarizability. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2012 ABINIT group (FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  iq=index of the q-point in the array Qmesh%ibz where epsilon^-1 has to be calculated
!!  Ep<Epsilonm1_parameters>=Structure with parameters and dimensions related to the inverse dielectric matrix.
!!  Pvc<vcoul_t>=Structure gathering data on the Coulombian interaction
!!  Qmesh<BZ_mesh_type>=Data type with information on the q-sampling
!!  Dtfil<Datafiles_type)>=variables related to files
!!  spaceComm=MPI communicator.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      coeffs_gausslegint,wrtout,xginv,xheev,xsum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_rpa_functional(gwrpacorr,iqcalc,iq,Ep,Pvc,Qmesh,Dtfil,gmet,chi0,spaceComm,ec_rpa)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

 use m_gwdefs,        only : GW_TOLQ0, epsilonm1_parameters
 use m_io_tools,      only : get_unit
 use m_abilasi,       only : xginv,xheev
 use m_geometry,      only : normv
 use m_bz_mesh,       only : bz_mesh_type
 use m_vcoul,         only : vcoul_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_rpa_functional'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqcalc,iq,gwrpacorr,spaceComm
 type(BZ_mesh_type),intent(in) :: Qmesh
 type(vcoul_t),intent(in) :: Pvc
 type(Datafiles_type),intent(in) :: Dtfil
 type(Epsilonm1_parameters),intent(in) :: Ep
!arrays
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(inout) :: ec_rpa(gwrpacorr)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,ilambda,io,istat,master,rank,nprocs,unt,ierr
 real(dp) :: ecorr
 real(dp) :: lambda
 logical :: qeq0
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: z(:),zl(:),zlw(:),zw(:)
 complex(gwpc),allocatable :: chi0_diag(:),chitmp(:,:),epsm_lf(:),epsm_nlf(:)
 complex(gwpc),pointer :: vc_sqrt(:)
 real(gwpc),allocatable :: eig(:)

! *************************************************************************

 DBG_ENTER("COLL")

! initialize MPI data
 master=0
 rank   = xcomm_rank(spaceComm)
 nprocs = xcomm_size(spaceComm)

 !if (rank==master) then ! presently only master has chi0 in screening

 ABI_ALLOCATE(epsm_lf,(Ep%nomega))
 ABI_ALLOCATE(epsm_nlf,(Ep%nomega))
 ! vc_sqrt contains vc^{1/2}(q,G), complex-valued to allow for a possible cutoff
 vc_sqrt => Pvc%vc_sqrt(:,iq)  ; qeq0=(normv(Qmesh%ibz(:,iq),gmet,'G')<GW_TOLQ0)

 ! Calculate Gauss-Legendre quadrature knots and weights for the omega integration
 ABI_ALLOCATE(zw,(Ep%nomegaei))
 ABI_ALLOCATE(z,(Ep%nomegaei))
 call coeffs_gausslegint(zero,one,z,zw,Ep%nomegaei)

 ! Calculate Gauss-Legendre quadrature knots and weights for the lambda integration
 ABI_ALLOCATE(zlw,(gwrpacorr))
 ABI_ALLOCATE(zl,(gwrpacorr))
 call coeffs_gausslegint(zero,one,zl,zlw,gwrpacorr)


 ABI_ALLOCATE(chi0_diag,(Ep%npwe))
 ABI_ALLOCATE(chitmp,(Ep%npwe,Ep%npwe))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out-of-memory in chitmp")

 do io=2,Ep%nomega 

   if(gwrpacorr==1) then ! exact integration over the coupling constant 

     if(modulo(io-2,nprocs)/=rank) cycle ! distributing the workload

     do ig2=1,Ep%npwe
       do ig1=1,Ep%npwe
         chitmp(ig1,ig2) = vc_sqrt(ig1) * vc_sqrt(ig2) * chi0(ig1,ig2,io)
       end do !ig1
     end do !ig2
     ABI_ALLOCATE(eig,(Ep%npwe))
     call xheev('V','U',Ep%npwe,chitmp,eig)

     do ig1=1,Ep%npwe
       ec_rpa(:) = ec_rpa(:) &
&         - zw(io-1) / ( z(io-1) * z(io-1) ) &
&              * Qmesh%wt(iq) * (-log( 1.0_dp-eig(ig1) )  - eig(ig1) ) / (2.0_dp * pi ) 
     end do
     ABI_DEALLOCATE(eig)

   else ! numerical integration over the coupling constant

     if(modulo( (ilambda-1)+gwrpacorr*(io-2),nprocs)/=rank) cycle ! distributing the workload

     do ilambda=1,gwrpacorr
       lambda=zl(ilambda)
       do ig1=1,Ep%npwe
         chi0_diag(ig1) = vc_sqrt(ig1)**2 * chi0(ig1,ig1,io)
       end do
     
       do ig2=1,Ep%npwe
         do ig1=1,Ep%npwe
           chitmp(ig1,ig2) = - lambda * vc_sqrt(ig1) * vc_sqrt(ig1) * chi0(ig1,ig2,io)
         end do !ig1
         chitmp(ig2,ig2) = chitmp(ig2,ig2) + 1.0_dp
       end do !ig2
       call xginv(chitmp(:,:),Ep%npwe)
       chitmp(:,:) = matmul( chi0(:,:,io) , chitmp(:,:) )
     
       do ig1=1,Ep%npwe
         chi0_diag(ig1) = vc_sqrt(ig1) * vc_sqrt(ig1) * chitmp(ig1,ig1) - chi0_diag(ig1)
       end do
     
       do ig1=1,Ep%npwe
         ec_rpa(ilambda) = ec_rpa(ilambda) &
&           - zw(io-1) / ( z(io-1) * z(io-1) ) * Qmesh%wt(iq) * real(  chi0_diag(ig1) ) / (2.0_dp * pi )
       end do

     end do ! ilambda

   end if ! exact or numerical integration over the coupling constant

 end do ! io

 
 ! Output the correlation energy when the last q-point to be calculated is reached
 ! This would allow for a manual parallelization over q-points
 if(iqcalc==Ep%nqcalc) then 

   call xsum_master(ec_rpa,master,spaceComm,ierr)

   if(rank==master) then
     ecorr = sum( zlw(:)*ec_rpa(:) ) 
     unt=get_unit() ; open(unit=unt,file=dtfil%fnameabo_rpa)
     write(unt,'(a,(2x,f14.8))') '#RPA',ecorr
     write(msg,'(2a,(2x,f14.8))') ch10,' RPA energy [Ha] :',ecorr
     call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
     if(gwrpacorr>1) then
       do ilambda=1,gwrpacorr
         write(unt,'(i6,2x,f10.6,2x,e12.6)') ilambda,zl(ilambda),ec_rpa(ilambda)
         write(msg,'(i6,2x,f10.6,2x,e12.6)') ilambda,zl(ilambda),ec_rpa(ilambda)
         call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
       end do
     end if
   end if

 end if

 ABI_DEALLOCATE(chi0_diag)
 ABI_DEALLOCATE(chitmp)
 ABI_DEALLOCATE(zl)
 ABI_DEALLOCATE(zlw)
 ABI_DEALLOCATE(z)
 ABI_DEALLOCATE(zw)

 DBG_EXIT("COLL")

end subroutine calc_rpa_functional
!!***
