!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dyson_solver
!! NAME
!!  m_dyson_solver
!!
!! FUNCTION
!!  This module contains procedures to solve the Dyson equation to find QP energies.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_dyson_solver

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_profiling

 use m_gwdefs,        only : czero_gw, sigma_parameters
 use m_numeric_tools, only : linfit, pade, dpade, newrap_step
 use m_abilasi,       only : xheev
 use m_bz_mesh,       only : bz_mesh_type, get_BZ_item
 use m_sigma_results, only : sigma_results

 implicit none

 private

 public :: solve_dyson 

 integer,private,parameter :: NR_MAX_NITER=1000            
  ! Max no of iterations in the Newton-Raphson method. 

 real(dp),private,parameter :: NR_ABS_ROOT_ERR=0.0001/Ha_eV 
  ! Tolerance on the absolute error on the Newton-Raphson root.

CONTAINS  !====================================================================
!!***

!!****f* m_dyson_solver/solve_dyson
!! NAME
!! solve_dyson
!!
!! FUNCTION
!!  Solve the Dyson equation for the QP energies. Two different methods are coded:
!!  The first one is based on the standard perturbative approach in which the self-energy
!!  is linearly expanded around the previous single-particle energy (KS energy if one-shot)
!!  and the derivative is evaluated by finite differences.
!!  In the second method (AC), the values of the self-energy operator on the real axis are obtained 
!!  by means of an analitic continuation based on the Pade extrapolation. 
!!
!! INPUTS
!!  ikcalc=Index of the considered k-point in the Sigp%kptgw2bz array.
!!  nomega_sigc=Number of frequencies used to evaluate the correlation part of Sigma.
!!  Sigp<Sigma_parameters>=Structure gathering parameters on the calculation of Sigma.
!!     %minbnd and %maxbnd= min and Max band index for GW correction (for this k-point)
!!     %gwcalctyp=Type of the GW calculation.
!!     %soenergy=Scissor energy
!!  Sr<Sigma_results>=Structure containing the matrix elements of the self-energy INOUT
!!     %nbnds=Number of bands in G0.
!!     %nsppol=Number of independent spin polarizations.
!!     %nsig_ab=Numner of components in the self-energy operator.
!!     %nomega_r=Number of real frequencies for spectral function.
!!     %nomega4sd=Number of real frequencies used to evalute the derivative of Sigma.
!!     %nomega_i=Number of imaginary frequencies for AC.
!!     %omega_i=Purely imaginary frequencies for AC.
!!  Kmesh<BZ_mesh_type>=Info on the K-mesh for the wavefunctions.
!!     %nkibz=Number of points in the IBZ
!!  sigxme_tmp(ib1:ib2,ib1:ib2,nsppol)=Matrix elements of Sigma_x.
!!  sigcme_tmp=(nomega_sigc,ib1:ib2,ib1:ib2,nsppol)=Matrix elements of Sigma_c.
!!  qp_ene(nbnds,nkibz,nsppol)= KS or QP energies, only used in case of calculation with scissor operator.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Sr<Sigma_results>=Structure containing the matrix elements of the self-energy:
!!     %sigxme(ib1:ib2,jkibz,nsspol)=Diagonal elements of Sigma_x
!!     %sigcmee0(ib1:ib2,jkibz,nsppol)=Matrix elements of Sigma_c at the initial energy E0.
!!     %dsigmee0(jb,ib1:ib2,nsppol)=Derivate of sigma at the energy E0.
!!     %ze0(ib1:ib2,jkibz,is)=Renormalization factor at the energy E0.
!!     %degw(ib1:ib2,jkibz,is)= QP correction  i.e DeltaE_GW=E-E0 
!!     %egw(ib1:ib2,jkibz,is)=QP energy
!!     %sigmee(ib1:ib2,jkibz,is)=Self-energy evaluated at the QP energy.
!!     %sigcme (ib1:ib2,jkibz,io,is)= Sigma_c as a function of frequency.
!!     %sigxcme(ib1:ib2,jkibz,io,is)= Sigma_xc as a function of frequency.
!!     %sigcme4sd (ib1:ib2,jkibz,io,is)= Diagonal matrix elements of \Sigma_c  at frequencies around the KS eigenvalue
!!     %sigxcme4sd(ib1:ib2,jkibz,io,is)= Diagonal matrix elements of \Sigma_xc at frequencies around the KS eigenvalue
!!    where ib1 and ib2 are the band indeces included in the GW calculation for this k-point.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      int2char,wrtout
!!
!! SOURCE

subroutine solve_dyson(ikcalc,minbnd,maxbnd,nomega_sigc,Sigp,Kmesh,sigxme_tmp,sigcme_tmp,qp_ene,Sr,prtvol,Dtfil,comm)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'solve_dyson'
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,nomega_sigc,prtvol,minbnd,maxbnd,comm
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Datafiles_type),intent(in) :: Dtfil
 type(Sigma_parameters),intent(in) :: Sigp
 type(Sigma_results),intent(inout) :: Sr
!arrays
 real(dp),intent(in) :: qp_ene(Sr%nbnds,Sr%nkibz,Sr%nsppol)
 complex(dpc),intent(in) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd,minbnd:maxbnd,Sigp%nsppol*Sigp%nsig_ab)
 complex(dpc),intent(in) :: sigxme_tmp(minbnd:maxbnd,minbnd:maxbnd,Sigp%nsppol*Sigp%nsig_ab)

!Local variables-------------------------------
!scalars
 integer :: iab,ib1,ib2,ikbz_gw,io,ioe0j,is,is_idx,isym,iter,itim,jb
 integer :: jkibz,kb,ld_matrix,mod10,nsploop,my_rank,master
 real(dp) :: alpha,beta,smrt
 complex(dpc) :: ctdpc,dct,dsigc,sigc,zz,phase
 logical :: converged,ltest
 character(len=500) :: msg
!arrays
 real(dp) :: kbz_gw(3),tsec(2)
 real(dp),allocatable :: e0pde(:),eig(:),scme(:) 
 complex(dpc),allocatable :: hdp(:,:),tmpcdp(:) 
 complex(dpc),allocatable :: hhartree(:,:,:),htotal(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(490,1,tsec) ! csigme(Dyson)

 master=0
 my_rank = xcomm_rank(comm)

 mod10=MOD(Sigp%gwcalctyp,10)

 ltest=(nomega_sigc==Sr%nomega_r+Sr%nomega4sd)
 if (mod10==1) ltest=(nomega_sigc==Sr%nomega_i)
 ABI_CHECK(ltest,'Wrong number of frequencies')

 ! Index of the KS or QP energy.
 ioe0j=Sr%nomega4sd/2+1

 ! min and Max band index for GW corrections (for this k-point).
 ib1=MINVAL(Sigp%minbnd(ikcalc,:))
 ib2=MAXVAL(Sigp%maxbnd(ikcalc,:))

 ! Index of this k-point in the IBZ array.
 ikbz_gw=Sigp%kptgw2bz(ikcalc) 
 call get_BZ_item(Kmesh,ikbz_gw,kbz_gw,jkibz,isym,itim,phase)
 !$call get_IBZ_item(Kmesh,jkibz,kibz,wtk)

 sigc=czero; dsigc=czero

 ! ===========================================================
 ! ==== Solve the Dyson Equation and store results in Sr% ====
 ! ===========================================================

 ! === Save elements or ab components of Sigma_x (hermitian) ===
 ! TODO It should be hermitian also in the spinorial case re-check
 do is=1,Sr%nsppol
   do jb=ib1,ib2

     do iab=1,Sr%nsig_ab
       is_idx=is ; if (Sr%nsig_ab>1) is_idx=iab
       Sr%sigxme(jb,jkibz,is_idx) = DBLE(sigxme_tmp(jb,jb,is_idx))
     end do
     if (Sr%nsig_ab>1) then 
       write(std_out,'(i3,4f8.3,a,f8.3)')jb,Sr%sigxme(jb,jkibz,:)*Ha_eV,' Tot ',SUM(Sr%sigxme(jb,jkibz,:))*Ha_eV
     end if

   end do 
 end do 

 if (mod10/=1) then 
  ! ===============================
  ! ==== Perturbative approach ==== 
  ! ===============================
  do is=1,Sr%nsppol
   do jb=ib1,ib2

    ! === Get matrix elements of Sigma_c at energy E0 ===
    ! * SigC(w) is linearly interpolated and the slope alpha is assumed as dSigC/dE
    do iab=1,Sr%nsig_ab
     is_idx=is ; if (Sr%nsig_ab>1) is_idx=iab

     Sr%sigcmee0(jb,jkibz,is_idx) = sigcme_tmp(Sr%nomega_r+ioe0j,jb,jb,is_idx)

     ABI_ALLOCATE(scme,(Sr%nomega4sd))
     ABI_ALLOCATE(e0pde,(Sr%nomega4sd))
     e0pde(:) = Sr%omega4sd(jb,jkibz,:,is)
     scme(:)  = REAL(sigcme_tmp(Sr%nomega_r+1:Sr%nomega_r+Sr%nomega4sd,jb,jb,is_idx))

     if (Sr%nomega4sd==1) then
       smrt=zero; alpha=zero
     else
       smrt=linfit(Sr%nomega4sd,e0pde(:),scme(:),alpha,beta)
     end if

     if (smrt>0.1/Ha_eV) then
       write(msg,'(3a,i4,a,i3,2a,2(f22.15,2a))')&
&        ' Values of Re Sig_c are not linear ',ch10,&
&        ' band index      = ',jb,' spin|component = ',is_idx,ch10,& 
&        ' root mean square= ',smrt,ch10,&
&        ' estimated slope = ',alpha,ch10,'Omega [eV] SigC [eV]' 
       MSG_WARNING(msg)
       do io=1,Sr%nomega4sd
         write(msg,'(2f8.4)')e0pde(io)*Ha_eV,scme(io)*Ha_eV
         call wrtout(std_out,msg,"COLL")
       end do
     end if
     ABI_DEALLOCATE(scme)
     ABI_DEALLOCATE(e0pde)
     !   
     ! === Evaluate renormalization factor and QP correction ===
     ! * Z=(1-dSigma/domega(E0))^-1
     ! * DeltaE_GW=E-E0= (Sigma(E0)-V_xc)/(1-dSigma/domega)
     ! * If nspinor==2, this part is done at the end.
     !
     Sr%dsigmee0(jb,jkibz,is_idx)=CMPLX(alpha,zero)

     if (Sr%nsig_ab==1) then
       Sr%ze0(jb,jkibz,is)=one/(one-Sr%dsigmee0(jb,jkibz,is))

       if (ABS(Sigp%soenergy) < tol6) then
         Sr%degw(jb,jkibz,is) = &
&          (Sr%sigxme(jb,jkibz,is)+Sr%sigcmee0(jb,jkibz,is)-Sr%e0(jb,jkibz,is)+Sr%hhartree(jb,jb,jkibz,is))*Sr%ze0(jb,jkibz,is) 

         Sr%egw(jb,jkibz,is)=Sr%e0(jb,jkibz,is)+Sr%degw(jb,jkibz,is)

         ! Estimate Sigma at the QP-energy: Sigma(E_qp)=Sigma(E0)+(E_qp-E0)*dSigma/dE
         Sr%sigmee(jb,jkibz,is)= &
&          Sr%sigxme(jb,jkibz,is)+Sr%sigcmee0(jb,jkibz,is)+Sr%degw(jb,jkibz,is)*Sr%dsigmee0(jb,jkibz,is)

       else ! If GW+scissor: e0 is replaced by qp_ene which contains the updated energy eigenvalue
         Sr%degw(jb,jkibz,is)= &
&          (Sr%sigxme(jb,jkibz,is)+Sr%sigcmee0(jb,jkibz,is)-qp_ene(jb,jkibz,is)+Sr%hhartree(jb,jb,jkibz,is))*Sr%ze0(jb,jkibz,is) 

         Sr%egw(jb,jkibz,is)=qp_ene(jb,jkibz,is)+Sr%degw(jb,jkibz,is)

         ! Estimate Sigma at the QP-energy: Sigma(E_qp)=Sigma(E0)+(E_qp-E0)*dSigma/dE
         Sr%sigmee(jb,jkibz,is)= &
&          Sr%sigxme(jb,jkibz,is)+Sr%sigcmee0(jb,jkibz,is)+Sr%degw(jb,jkibz,is)*Sr%dsigmee0(jb,jkibz,is)

         ! RS: In the output, the gw corr with respect to e0 without soenergy is reported.
         Sr%degw(jb,jkibz,is)=Sr%egw(jb,jkibz,is)-Sr%e0(jb,jkibz,is)
       end if
     end if !Sigp%nsig_ab==1

      ! Spectrum of Sigma
      do io=1,Sr%nomega_r
        Sr%sigcme (jb,jkibz,io,is_idx)= sigcme_tmp(io,jb,jb,is_idx)
        Sr%sigxcme(jb,jkibz,io,is_idx)= Sr%sigxme(jb,jkibz,is_idx)+Sr%sigcme(jb,jkibz,io,is_idx)
      end do
      do io=1,Sr%nomega4sd
        Sr%sigcme4sd (jb,jkibz,io,is_idx)= sigcme_tmp(Sr%nomega_r+io,jb,jb,is_idx)
        Sr%sigxcme4sd(jb,jkibz,io,is_idx)= Sr%sigxme(jb,jkibz,is_idx)+Sr%sigcme4sd(jb,jkibz,io,is_idx)
      end do

    end do !iab

    if (Sr%nsig_ab>1) then
      ABI_CHECK(ABS(Sigp%soenergy)<0.1d-4,'Scissor with spinor not coded')
      !TODO this should be allocated with nsppol, recheck this part

      ! === Evaluate renormalization factor and QP correction ===
      ! * Z=(1-dSigma/domega(E0))^-1
      ! * DeltaE_GW=E-E0= (Sigma(E0)-V_xc)/(1-dSigma/domega)
      write(std_out,'(a,i2,10f8.3)')' Correlation',jb,Sr%sigcmee0(jb,jkibz,:)*Ha_eV,SUM(Sr%sigcmee0(jb,jkibz,:))*Ha_eV

      Sr%ze0 (jb,jkibz,1) = one/(one-SUM(Sr%dsigmee0(jb,jkibz,:)))

      Sr%degw(jb,jkibz,1) = Sr%ze0(jb,jkibz,1) * &
&       (SUM(Sr%sigxme(jb,jkibz,:)+Sr%sigcmee0(jb,jkibz,:)+Sr%hhartree(jb,jb,jkibz,:))-Sr%e0(jb,jkibz,1))

      Sr%egw(jb,jkibz,1)=Sr%e0(jb,jkibz,1)+Sr%degw(jb,jkibz,1)

      ! === Estimate Sigma at the QP-energy ===
      do iab=1,Sr%nsig_ab
       Sr%sigmee(jb,jkibz,iab)= &
&        Sr%sigxme(jb,jkibz,iab)+Sr%sigcmee0(jb,jkibz,iab)+Sr%degw(jb,jkibz,1)*Sr%dsigmee0(jb,jkibz,iab)
      end do
    end if

   end do !jb
  end do !is

 else 
  ! =============================
  ! === Analytic Continuation ===
  ! =============================
  ABI_CHECK(Sr%nsig_ab==1,"AC with spinor not implemented")
  do is=1,Sr%nsppol
   do jb=ib1,ib2

    ABI_ALLOCATE(tmpcdp,(Sr%nomega_i))
    ! * Calculate Sigc(E0), dSigc(E0)
    zz=CMPLX(Sr%e0(jb,jkibz,is),zero)

    if (Sigp%soenergy>0.1d-4) then
     ! RS: e0 is replaced by qp_ene which contains the updated energy eigenvalue
     zz=CMPLX(qp_ene(jb,jkibz,is),zero)   
    end if

    ! === Diagonal elements of sigcme_tmp ===
    ! * if zz in 2 or 3 quadrant, avoid poles in the complex plane using Sigma(-iw)=Sigma(iw)*.
    do iab=1,Sr%nsig_ab
      is_idx=is; if (Sr%nsig_ab>1) is_idx=iab
      if (REAL(zz)>zero) then
        tmpcdp(:)=sigcme_tmp(:,jb,jb,is_idx)
        Sr%sigcmee0(jb,jkibz,is_idx)=  pade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
        Sr%dsigmee0(jb,jkibz,is_idx)= dpade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
      else
        tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,is_idx))
        Sr%sigcmee0(jb,jkibz,is_idx)=  pade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
        Sr%dsigmee0(jb,jkibz,is_idx)= dpade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
      end if
    end do !iab

    ! Z=(1-dSigma/domega(E0))^-1
    if (Sr%nsig_ab==1) then
      Sr%ze0(jb,jkibz,is)=one/(one-Sr%dsigmee0(jb,jkibz,is))
    else
      Sr%ze0(jb,jkibz,1)=one/(one-SUM(Sr%dsigmee0(jb,jkibz,:)))
    end if

    ! Find roots of E^0-V_xc-V_U+Sig_x+Sig_c(z)-z, i.e E^qp.
    ! using Newton-Raphson method and starting point E^0
    zz=CMPLX(Sr%e0(jb,jkibz,is),zero)

    if (Sigp%soenergy>0.1d-4) then ! e0 is replaced by qp_ene which contains the updated energy eigenvalue.
      zz=CMPLX(qp_ene(jb,jkibz,is),0.0)
    end if

    iter=0; converged=.FALSE.; ctdpc=cone
    do while (ABS(ctdpc)>NR_ABS_ROOT_ERR.or.iter<NR_MAX_NITER)
      iter=iter+1
      sigc=czero ; dsigc=czero
      if (REAL(zz)>tol12) then
        tmpcdp(:)=sigcme_tmp(:,jb,jb,is)
        sigc =  pade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
        dsigc= dpade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
      else
        tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,is))
        sigc =  pade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
        dsigc= dpade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
      end if
      ctdpc=Sr%e0(jb,jkibz,is)-Sr%vxcme(jb,jkibz,is)-Sr%vUme(jb,jkibz,is)+Sr%sigxme(jb,jkibz,is)+sigc-zz
      if (ABS(ctdpc)<NR_ABS_ROOT_ERR) then 
       converged=.TRUE.; EXIT
      end if
      dct=dsigc-one
      zz=newrap_step(zz,ctdpc,dct)
    end do

    if (.not.converged) then 
      write(msg,'(a,i7,3a,f8.4,a,f8.4)')&
&       ' Newton-Raphson method not converged after ',NR_MAX_NITER,' iterations. ',ch10,&
&       ' Absolute Error= ',ABS(ctdpc),' > ',NR_ABS_ROOT_ERR
      MSG_WARNING(msg)
    end if
    !   
    ! Store the final result TODO re-shift everything according to efermi
    Sr%egw(jb,jkibz,is)=zz
    Sr%degw(jb,jkibz,is)=Sr%egw(jb,jkibz,is) - Sr%e0(jb,jkibz,is)
    Sr%sigmee(jb,jkibz,is)=Sr%sigxme(jb,jkibz,is) + sigc
    !   
    ! Spectra of Sigma, remember that Sr%nomega_r does not contains the frequencies used to evaluate the derivative
    ! each frequency is obtained using the pade_expression
    do io=1,Sr%nomega_r
      zz=Sr%omega_r(io)
      if (REAL(zz)>zero) then
        tmpcdp(:)=sigcme_tmp(:,jb,jb,is)
        Sr%sigcme(jb,jkibz,io,is) = pade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
      else
        tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,is))
        Sr%sigcme(jb,jkibz,io,is) = pade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
      end if
      Sr%sigxcme(jb,jkibz,io,is)= Sr%sigxme(jb,jkibz,is)+Sr%sigcme(jb,jkibz,io,is)
    end do
    !   
    ! === Save sigma values along the imaginary axis ===
    do iab=1,Sr%nsig_ab
      is_idx=is ; if (Sr%nsig_ab>1) is_idx=iab
      do io=1,Sr%nomega_i
        Sr%sigcmesi (jb,jkibz,io,is_idx)= sigcme_tmp(io,jb,jb,is_idx)
        Sr%sigxcmesi(jb,jkibz,io,is_idx)= Sr%sigxme(jb,jkibz,is_idx)+Sr%sigcmesi(jb,jkibz,io,is_idx)
      end do
    end do
 
    ABI_DEALLOCATE(tmpcdp)

   end do !jb
  end do !is
 end if ! Analytic continuation.
 !
 ! === Diagonalize the QP Hamiltonian (forced to be Hermitian) ===
 ! * Calculate Sr%en_qp_diago and Sr%eigvec_qp to be written in the QPS file.
 ! TODO in case of AC results are wrong.

 ABI_ALLOCATE(hhartree,(ib1:ib2,ib1:ib2,Sr%nsppol*Sr%nsig_ab))
 hhartree=Sr%hhartree(ib1:ib2,ib1:ib2,jkibz,:)

 ! If non self-consistent erase all off-diagonal elements
 if (Sigp%gwcalctyp<20) then
   do jb=ib1,ib2
     do kb=ib1,ib2
      if (jb==kb) CYCLE
      hhartree(jb,kb,:)=czero
     end do
   end do
 end if

 ABI_ALLOCATE(htotal,(ib1:ib2,ib1:ib2,Sr%nsppol*Sr%nsig_ab))
 htotal = hhartree + sigxme_tmp(:,:,:) + sigcme_tmp(Sr%nomega_r+ioe0j,:,:,:)
 !
 ! === Get the Hermitian part of htotal ===
 ! * In the noncollinear case A_{12}^{ab} = A_{21}^{ba}^* if A is Hermitian.
 nsploop=Sr%nsppol ; if (Sr%nsig_ab/=1) nsploop=2
 do is=1,nsploop
   htotal(:,:,is)= half*(htotal(:,:,is)+TRANSPOSE(CONJG(htotal(:,:,is))))
 end do 
 !
 ! Print the different matrix elements of sigma if QPSC and prtvol>9
 if (Sigp%gwcalctyp>=20.and.prtvol>9.and.my_rank==master) then
   call print_sigma_melems(ikcalc,ib1,ib2,Sr%nsppol*Sr%nsig_ab,htotal,hhartree,&
&               sigxme_tmp,sigcme_tmp(Sr%nomega_r+ioe0j,:,:,:),Dtfil%filnam_ds(4)) 
 end if

 if (Sr%nsig_ab==4) then
   htotal(:,:,3)= half*(htotal(:,:,3)+TRANSPOSE(CONJG(htotal(:,:,4))))
   htotal(:,:,4)= TRANSPOSE(CONJG(htotal(:,:,3)))
 end if
 !
 ! === Solve Herm(htotal)*U = E*U ===
 ld_matrix=ib2-ib1+1
 ABI_ALLOCATE(hdp,(ld_matrix,ld_matrix))
 ABI_ALLOCATE(eig,(ld_matrix))

 do is=1,Sr%nsppol
   if (Sr%nsig_ab==1) then
     hdp=htotal(ib1:ib2,ib1:ib2,is)
   else 
     hdp=SUM(htotal(ib1:ib2,ib1:ib2,:),DIM=3)
   end if

   call xheev("Vectors","Upper",ld_matrix,hdp,eig)

   Sr%eigvec_qp(ib1:ib2,ib1:ib2,jkibz,is)=hdp(:,:)
   Sr%en_qp_diago(ib1:ib2,jkibz,is)=eig(:)
 end do 

 ABI_DEALLOCATE(hdp)
 ABI_DEALLOCATE(eig)
 ABI_DEALLOCATE(htotal)
 ABI_DEALLOCATE(hhartree)

 call timab(490,2,tsec)

 DBG_EXIT("COLL")

end subroutine solve_dyson
!!***

!----------------------------------------------------------------------

!!****f* m_dyson_solver/print_sigma_melems
!! NAME
!!  print_sigma_melems
!!
!! FUNCTION
!!  This routine prints the Hermitian and the non-hermitian part of the matrix
!!  elements of Sigma, as well as the individual contributions.
!!  The first 14x14 are printed to screen, and the full matrices are printed
!!  to files: sigma_melems_, sigma_nonH_melems_, sigma_Hart_melems_,
!!            sigma_x_melems, and sigma_c_melems
!!
!! INPUTS
!!  ikcalc  : index of k-point
!!  ib1,ib2 : starting and ending band indices
!!  nsp     : no. of spin elements
!!  htotal  : Hermitianised matrix elements of Sigma
!!  hhartree : Hartree contribution to matrix elements
!!  sigxme  : Sigma_x contribution to matrix elements
!!  sigcme  : Sigma_c contribution to matrix elements
!!  prefil : prefix for output files.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dyson_solver
!!
!! CHILDREN
!!      int2char,wrtout
!!
!! SOURCE

subroutine print_sigma_melems(ikcalc,ib1,ib2,nsp,htotal,hhartree,sigxme,sigcme,prefil)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_sigma_melems'
 use interfaces_14_hidewrite
 use interfaces_27_toolbox_oop
!End of the abilint section

 implicit none

! Arguments ------------------------------------
 !scalars
 integer,intent(in) :: ikcalc,ib1,ib2,nsp
 character(len=*),intent(in) :: prefil
 !arrays
 complex(dpc),intent(in) :: htotal(ib1:ib2,ib1:ib2,nsp),hhartree(ib1:ib2,ib1:ib2,nsp)
 complex(dpc),intent(in) :: sigxme(ib1:ib2,ib1:ib2,nsp),sigcme(ib1:ib2,ib1:ib2,nsp)

! Local variables ------------------------------ 
 integer,parameter :: MAX_NCOLS=14
 integer :: isp,mc,mr,jj,ii
 character(len=10) :: sidx
 character(len=500) :: msg
 character(len=100) :: fmth,fmt1,fmt2,fmthh,kpt_index,fmtfile
 character(len=fnlen) :: filename
! *************************************************************************

 if (nsp==3.or.nsp>4) then
   MSG_ERROR(' nsp has wrong value in print_sigma_melems')   
 end if

 mc = ib2-ib1+1; if (mc>MAX_NCOLS) mc = MAX_NCOLS
 mr = mc

 write(fmthh,*)'(2(a),2(I2,a))'
 write(fmth,*)'(7x,',mc,'(i2,8x))'
 write(fmt1,*)'(3x,i2,',mc,'f10.5)'
 write(fmt2,*)'(5x   ,',mc,'f10.5,a)'

! First print to screen
 do isp=1,nsp
   write(msg,'(a)') ''
   call wrtout(std_out,msg,'COLL')
   write(msg,fmthh) ch10,' Hermitianised matrix elements of Sigma (spin ',isp,' of ',nsp,'):'
   call wrtout(std_out,msg,'COLL')
   write(msg,fmth)(jj,jj=1,mc)
   call wrtout(std_out,msg,'COLL') !header
   do ii=ib1,ib1+mr-1
     write(msg,fmt1)ii-ib1+1,DBLE(htotal(ii,ib1:(ib1+mc-1),isp))
     call wrtout(std_out,msg,'COLL') !real part
     write(msg,fmt2)  AIMAG(htotal(ii,ib1:(ib1+mc-1),isp)),ch10
     call wrtout(std_out,msg,'COLL') !imag part
   end do
 end do !nsp

 write(msg,'(a,i2,a)')" Max. ",MAX_NCOLS," elements printed. Full matrix output in _HTOTAL files"
 call wrtout(std_out,msg,'COLL')

 do isp=1,nsp
   write(msg,fmthh) ch10,' H_Hartree matrix elements (spin ',isp,' of ',nsp,'):'
   call wrtout(std_out,msg,'COLL') 
   write(msg,fmth)(jj,jj=1,mc)
   call wrtout(std_out,msg,'COLL') !header 
   do ii=ib1,ib1+mr-1 
     write(msg,fmt1)ii-ib1+1,DBLE(hhartree(ii,ib1:(ib1+mc-1),isp))
     call wrtout(std_out,msg,'COLL') !real part
     write(msg,fmt2)  AIMAG(hhartree(ii,ib1:(ib1+mc-1),isp)),ch10
     call wrtout(std_out,msg,'COLL') !imag part
   end do 
 end do !nsp

 write(msg,'(a,i2,a)')" Max. ",MAX_NCOLS," elements printed. Full matrix output in _HHARTREE files"
 call wrtout(std_out,msg,'COLL') 

 do isp=1,nsp
   write(msg,fmthh) ch10,' Sigma_x matrix elements (spin ',isp,' of ',nsp,'):'
   call wrtout(std_out,msg,'COLL')
   write(msg,fmth)(jj,jj=1,mc)
   call wrtout(std_out,msg,'COLL') !header 
   do ii=ib1,ib1+mr-1
     write(msg,fmt1)ii-ib1+1,DBLE(sigxme(ii,ib1:(ib1+mc-1),isp))
     call wrtout(std_out,msg,'COLL') !real part
     write(msg,fmt2)  AIMAG(sigxme(ii,ib1:(ib1+mc-1),isp)),ch10
     call wrtout(std_out,msg,'COLL') !imag part
   end do
 end do !nsp

 write(msg,'(a,i2,a)')" Max. ",MAX_NCOLS," elements printed. Full matrix output _SIGX files"
 call wrtout(std_out,msg,'COLL') 

 do isp=1,nsp
   write(msg,fmthh) ch10,' Sigma_c matrix elements (spin ',isp,' of ',nsp,'):'
   call wrtout(std_out,msg,'COLL')
   write(msg,fmth)(jj,jj=1,mc)
   call wrtout(std_out,msg,'COLL') !header 
   do ii=ib1,ib1+mr-1
     write(msg,fmt1)ii-ib1+1,DBLE(sigcme(ii,ib1:(ib1+mc-1),isp))
     call wrtout(std_out,msg,'COLL') !real part
     write(msg,fmt2)  AIMAG(sigcme(ii,ib1:(ib1+mc-1),isp)),ch10
     call wrtout(std_out,msg,'COLL') !imag part
   end do
 end do !nsp

 write(msg,'(a,i2,a)')" Max ",MAX_NCOLS," elements printed. Full matrix output _SIGC files"
 call wrtout(std_out,msg,'COLL') 

 ! Then print to file
 ! Format is: row, column, value; with a blank space for each full
 ! set of columns for easy plotting with the gnuplot splot command

 write(fmtfile,*)'(3X,I6,2X,I6,',nsp,'(2(ES28.16E3,3x)))'
 
 call int2char(ikcalc,sidx)
 kpt_index = "_KPT"//TRIM(sidx)

 filename = TRIM(prefil)//'_HTOTAL'//TRIM(kpt_index)
 open(unit=tmp_unit,file=filename,status='replace',action='write')
 msg = '#   row    col.      Re(htotal(r,c)) Im(htotal(r,c))  for spin11   ... spin22 ... spin12 ... spin13'
 call wrtout(tmp_unit,msg,'COLL')
 do ii=ib1,ib2
   do jj=ib1,ib2
     write(msg,fmtfile) ii,jj,(htotal(jj,ii,isp),isp=1,nsp)
     call wrtout(tmp_unit,msg,'COLL')
   end do
   call wrtout(tmp_unit,"",'COLL')
 end do
 close(tmp_unit)

 filename = TRIM(prefil)//'_HHARTREE'//TRIM(kpt_index)
 open(unit=tmp_unit,file=filename,status='replace',action='write')
 msg = '#   row    col.      Re(hhartree(r,c))  Im(hhartree(r,c)  for spin11   ... spin22 ... spin12 ... spin13'
 call wrtout(tmp_unit,msg,'COLL')
 do ii=ib1,ib2
   do jj=ib1,ib2
     write(msg,fmtfile) ii,jj,(hhartree(jj,ii,isp),isp=1,nsp)
     call wrtout(tmp_unit,msg,'COLL')
   end do
   call wrtout(tmp_unit,"",'COLL')
 end do
 close(tmp_unit)

 filename = TRIM(prefil)//'_SIGX'//TRIM(kpt_index)
 open(unit=tmp_unit,file=filename,status='replace',action='write')
 write(msg,'(a)')'#   row    col.      Re(Sigx(r,c)) Im(Sigx(r,c) for spin11   ... spin22 ... spin12 ... spin13'
 call wrtout(tmp_unit,msg,'COLL')
 do ii=ib1,ib2
   do jj=ib1,ib2
     write(msg,fmtfile) ii,jj,(sigxme(jj,ii,isp),isp=1,nsp)
     call wrtout(tmp_unit,msg,'COLL')
   end do
   call wrtout(tmp_unit,"",'COLL')
 end do
 close(tmp_unit)

 filename = TRIM(prefil)//'_SIGC'//TRIM(kpt_index)
 open(unit=tmp_unit,file=filename,status='replace',action='write')
 write(msg,'(a)')'#   row    col.      Re(Sigc(r,c)) Im(Sigc(r,c) for spin11   ... spin22 ... spin12 ... spin21'
 call wrtout(tmp_unit,msg,'COLL')
 do ii=ib1,ib2
   do jj=ib1,ib2 
     write(msg,fmtfile) ii,jj,(sigcme(jj,ii,isp),isp=1,nsp)
     call wrtout(tmp_unit,msg,'COLL') 
   end do 
   call wrtout(tmp_unit,"",'COLL') 
 end do 
 close(tmp_unit) 

end subroutine print_sigma_melems

!----------------------------------------------------------------------

END MODULE m_dyson_solver
!!***
