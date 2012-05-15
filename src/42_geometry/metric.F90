!{\src2tex{textfont=tt}}
!!****f* ABINIT/metric
!! NAME metric
!! metric
!!
!! FUNCTION
!! Compute first dimensional primitive translation vectors in reciprocal space
!! gprimd from rprimd, and eventually writes out.
!! Then, computes metrics for real and recip space rmet and gmet using length
!! dimensional primitive translation vectors in columns of rprimd(3,3) and gprimd(3,3).
!!  gprimd is the inverse transpose of rprimd.
!!  i.e. $ rmet_{i,j}= \sum_k ( rprimd_{k,i}*rprimd_{k,j} )  $
!!       $ gmet_{i,j}= \sum_k ( gprimd_{k,i}*gprimd_{k,j} )  $
!! Also computes unit cell volume ucvol in $\textrm{bohr}^3$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  iout=unit number of output file.  If iout<0, do not write output.
!!
!! OUTPUT
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  ucvol=unit cell volume ($\textrm{bohr}^{3}$).
!!
!! PARENTS
!!      afterscfloop,bethe_salpeter,chkinp,clnup1,conducti_nc,conducti_paw
!!      conducti_paw_core,dist2,dyfnl3,eltfrhar3,eltfrkin3,eltfrnl3,eltfrxc3
!!      emispec_paw,energy,extrapwf,fftprof,forces,forstrnps,getkgrid,gw_tools
!!      hirsh,ingeo,initaim,inkpts,inqpt,invacuum,invars2m,irred_perts
!!      ks_ddiago,kss2wfk,ladielmt,lavnl,linear_optics_paw,localorb_S,loper3
!!      m_ab6_symmetry,m_commutator_vkbr,m_crystal,m_fft_prof,m_gwannier
!!      m_hamiltonian,m_paw_pwij,m_screening,m_use_ga,m_vcoul,m_wannier2abinit
!!      mlwfovlp_qp,moddiel,mover,mrgscr,newrho,newsp,newvtr,newvtr3,nres2vres
!!      optic,overlap_wf,partial_dos_fractions,pawgrnl,prcref,prcref_PMA
!!      prctfvw1,prctfvw2,prctfw3,pred_bfgs,pred_delocint,pred_diisrelax
!!      pred_isokinetic,pred_isothermal,pred_langevin,pred_nose,pred_srkna14
!!      pred_steepdesc,pred_verlet,prt_cif,prtimg,rdddb9,rdm,rhohxc,scfcv
!!      scfcv3,screening,setup1,setup_bse,setup_screening,setup_sigma,sigma
!!      smallprim,stress,suscep,testkgrid,tetrahedron,thmeig,wffile,wfread
!!      wrt_moldyn_netcdf,xfpack_vin2x,xfpack_x2vin
!!
!! CHILDREN
!!      matr3inv,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

 use m_profiling

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'metric'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout
 real(dp),intent(out) :: ucvol
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(out) :: gmet(3,3),gprimd(3,3),rmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: nu
 character(len=500) :: message
!arrays
 real(dp) :: angle(3)

! *************************************************************************

!DEBUG
!write(std_out,*)' metric : enter '
!ENDDEBUG

!Compute unit cell volume
 ucvol=rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
& rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
& rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))

!Check that the input primitive translations are not linearly dependent (and none is zero); i.e. ucvol~=0
!Also ask that the mixed product is positive.
 if (abs(ucvol)<tol12) then 
!  write(std_out,*)"rprimd",rprimd,"ucvol",ucvol
   write(message,'(5a)')&
&   '  Input rprim and acell gives vanishing unit cell volume.',ch10,&
&   '  This indicates linear dependency between primitive lattice vectors',ch10,&
&   '  Action : correct either rprim or acell in input file.'
   MSG_ERROR(message)
 end if
 if (ucvol<zero)then
   write(message,'(2a,3(a,3es16.6,a),7a)')&
&   '  Current rprimd gives negative (R1xR2).R3 . ',ch10,&
&   '  Rprimd =',rprimd(:,1),ch10,&
&   '          ',rprimd(:,2),ch10,&
&   '          ',rprimd(:,3),ch10,&
&   '  Action : if the cell size and shape are fixed (optcell==0),',ch10,&
&   '   exchange two of the input rprim vectors;',ch10,&
&   '   if you are optimizing the cell size and shape (optcell/=0),',ch10,&
&   '   maybe the move was too large, and you might try to decrease strprecon.'
   MSG_ERROR(message)
 end if

!Generates gprimd
 call matr3inv(rprimd,gprimd)

!Write out rprimd, gprimd and ucvol
 if (iout>=0) then
   write(message,'(2a)')' Real(R)+Recip(G) ','space primitive vectors, cartesian coordinates (Bohr,Bohr^-1):'
   call wrtout(iout,message,'COLL')
   do nu=1,3
     write(message, '(1x,a,i1,a,3f11.7,2x,a,i1,a,3f11.7)' ) &
&     'R(',nu,')=',rprimd(:,nu)+tol10,&
&     'G(',nu,')=',gprimd(:,nu)+tol10
     call wrtout(iout,message,'COLL')
   end do
   write(message,'(a,1p,e15.7,a)') ' Unit cell volume ucvol=',ucvol+tol10,' bohr^3'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Compute real space metric.
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)   

!Compute reciprocal space metric.
 gmet = MATMUL(TRANSPOSE(gprimd),gprimd)

!Write out the angles
 if (iout>=0) then
   angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0d0
   angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0d0
   angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0d0
   write(message, '(a,3es16.8,a)' )' Angles (23,13,12)=',angle(1:3),' degrees'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

end subroutine metric
!!***
