!{\src2tex{textfont=tt}}
!!****f* ABINIT/acfd_dyson
!! NAME
!! acfd_dyson
!!
!! FUNCTION
!! In reciprocal space, computes the dynamical interacting susceptibility matrix
!! from the non-interacting one for different Coulomb coupling constants $\lambda$,
!! solving the Dyson equation
!!  $\chi_\lambda = chi_0 + \chi_0 K_\lambda \chi_\lambda$,
!! where $K_\lambda$ is the Coulomb + exchange-correlation kernel at coupling
!! constant $\lambda$.
!! Also computes:
!!  a - The contribution to the exchange-correlation energy according to the
!!      adiabatic-connection fluctuation-dissipation-theorem, doing traces like
!!       $Tr[V_{C}\chi_{\lambda}] = 4\pi\sum_{\vec G}\frac{\chi_{\lambda}(\vec G,\vec G)}{G^2}$,
!!      and integrating them over the coupling constant. For the latter integration the mesh
!!      is generated automatically.
!!  b - Dipole polarizabilities, obtained like
!!       $\alpha = -\lim_{G\rightarrow 0}\frac{\chi(\vec G,\vec G)}{G^2}$ \Omega_{\rm cell}$.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, MF, XG, GMR, LSI, YMN).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  freq = (imaginary) frequency of the calculation (only used for output purposes).
!!  gsq(npwdiel) = the squared norm of the planewaves.
!!  idyson = 1 solve the Dyson equation as linear system.
!!         = 2 solve the Dyson equation as a differential equation.
!!         = 3 solve the Dyson equation iteratively.
!!  ig_tiny(npw_tiny,3) = planewave indices for the npw_tiny shortest G vectors
!!                        along each primitive translation.
!!  igsq_tiny(npw_tiny) = planewave indices for the npw_tiny shortest G vectors.
!!  ikhxc = option for the TDDFT kernel (see defs_parameters.f)
!!  kg_diel(3,npwdiel) = reduced planewave coordinates for the susceptibility matrices.
!!  khxc = a pointer to a (2,npwdiel,nspden,npwdiel,nspden) array containing the
!!         PGG kernel in the PGG approximation.
!!       = an unassociated/unallocated pointer in any other approximation.
!!  ldgapp = 0 Do not calculate the Lein, Dobson and Gross first-order approximation.
!!         > 0 Calculate the Lein, Dobson and Gross first-order approximation (see below).
!!  mpi_enreg=informations about MPI parallelization
!!  ndyson = number of steps in the coupling constant integration.
!!  nfft = number of fft grid points.
!!  ngfft(1:3) = integer fft box dimensions, see getng for ngfft(4:8).
!!  npw_tiny = number of shortest G vectors.
!!  npwdiel = number of planewaves for the susceptibility matrix.
!!  nspden = number of spin-density components.
!!  option = 0 do not calculate and print polarizabilities and
!!             differential correlation energies.
!!         = 1 calculate and print polarizabilities.
!!         > 1 calculate and print polarizabilities as well as
!!             differential correlation energies.
!!  rcut_coulomb = real space cut-off radius for Coulomb interaction in Bohr.
!!  rhor(nfft,nspden) = electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if nspden = 2).
!!  rhocut = cut-off density for the local kernels (ALDA, EOK),
!!           relative to max(rhor(:,:)).
!!  rprimd(3,3) = dimensional primitive translations for real space in Bohr.
!!  suskxcrs=1 if the product chi_0*Kxc is evaluated in real space , 0 otherwise
!!  ucvol = unit cell volume.
!!
!! OUTPUT
!!  susd_data(npwdiel,1) = real part of the coupling-constant integrated diagonal
!!   of $\chi_{\lambda}-chi_0$ (summed over spin-density components if appropriate).
!!  susd_data(npwdiel,2) = dto. in the Lein, Dobson and Gross first-order approximation,
!!   see below (summed over spin-density components if appropriate).
!!  susd_data(npwdiel,3) = real part of the diagonal of $\chi_{\lambda=1}$
!!   (summed over spin-density components if appropriate).
!! Also to the log file (see input variable option above):
!!  dEc*, d2Ec* = various (differential) correlation energies.
!!  alpha* = various polarizabilities.
!!
!! SIDE EFFECTS
!!  susmat(2,npwdiel,nspden,npwdiel,nspden) =
!!   on input:  the Kohn-Sham susceptibility matrix $\chi_{0}$.
!!   on output: the interacting susceptibility matrix $\chi_{\lambda=1}$.
!!
!! OPTIONS
!! The Dyson equation can be solved in different ways:
!!  idsyon = 1, as a set of linear equations:
!!    $(1-K_{\lambda}\chi_0)\ chi_{\lambda} = \chi_0$,
!!   for $\lambda$ = 0, 1, and ndyson Gauss-Legendre abscissas for the coupling constant
!!   integration. A recursive solution is implemented:
!!    $(1-[K_{\lambda_i}-K_{\lambda_{i-1}}])\chi_{\lambda_i} = \chi{_\lambda_{i-1}}$.
!!  idyson = 2, as a set of differential equations:
!!    $d/d\lambda\chi_{\lambda} = \chi_{\lambda} d/d\lambda K_{\lambda} \chi_{\lambda}$.
!!   ndyson steps are used in a leap frog integration scheme on a linear mesh.
!!  idyson = 3, iteratively, by self-consistently computing the linear density change,
!!    $\delta n_{p+1} = \chi_0(\delta v^{\rm ext}+\delta v^{\rm HXC}_{\lambda, p})$,
!!   and projecting $\delta n$ on $\delta v^{\rm ext}$ which yields $\chi_{\lambda}$.
!!   Here $p$ is the iteration index. Used to save memory: only $\chi_0$ is kept
!!   in memory, and only the diagonal of $\chi_{\lambda}$ is computed. Tradeoff:
!!   Computing time increases with the number of self-consistency iterations, will
!!   likely be several times less efficient than idsyon = 6.
!! ndyson > 0 solves the Dyson equation for ndyson $\lambda$ values between 0 and 1.
!!        = 0 does so only for $\lambda$ = 0 and 1.
!! If ldgapp > 0, the Lein, Dobson and Gross first-order approximation is also calculated:
!!   $\chi_{\lambda}-\chi_0 = \chi_0  [d/d\lambda K]_{\lambda=0} \chi_0$.
!!  see Lein, Dobson and Gross, J. Comput. Chem. 20, 12 (1999).
!!
!! NOTES
!! The Coulomb interaction is treated in different ways. For solving the Dyson
!! equation a cutoff Coulomb interaction (RPA kernel) is used, where the cutoff radius
!! rcut_coulomb is input and should always be smaller than the smallest
!! cell dimension. In this case the G = 0 term does not contribute to traces like
!! $Tr[V_C \chi]$. The trace is also computed with full Coulomb interaction, there
!! the $\vec G=0$ term is nonzero and extrapolated the values at the shortest G vectors
!! along each direction in the reciprocal lattice (o.k. for tetragonal unit cells).
!!
!! WARNINGS
!! Current restrictions are:
!!  a - Spin-polarized case not tested (might work in the RPA,
!!      but will likely be too memory intensive).
!!  b - Extrapolations to $\vec G=0$ components o.k. for tetragonal unit cells only.
!!  c - idyson = 3 (iterative solution of the Dyson equation) is preliminary.
!!
!! TODO
!!  a - The fact that khxc(:,ipw,isp1,ipw,isp2) = khxc(:,ipw,isp2,ipw,isp1) is not taken
!!      into account in this subroutine, which could save memory.
!!
!! NOTES
!!
!! PARENTS
!!      xcacfd
!!
!! CHILDREN
!!      dyson_de,dyson_gl,dyson_ls,dyson_sc,get_susd_null,geteexc_cc,geteexc_uc
!!      getlambda,k_rpa,klocal,kxc_alda,kxc_eok,leave_new,timab,wrtout,zheev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine acfd_dyson(dtset,freq,gsq,idyson,ig_tiny,igsq_tiny,ikhxc,kg_diel,khxc,&
&                      ldgapp,mpi_enreg,ndyson,nfft,ngfft,npw_tiny,npwdiel,nspden,&
&                      option,rcut_coulomb,rhor,rhocut,rprimd,susd_data,suskxcrs,susmat,ucvol)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'acfd_dyson'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_77_suscep, except_this_one => acfd_dyson
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: idyson,ikhxc,ldgapp,ndyson,nfft,npw_tiny,npwdiel,nspden
 integer,intent(in) :: option,suskxcrs
 real(dp),intent(in) :: freq,rcut_coulomb,rhocut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ig_tiny(npw_tiny,3),igsq_tiny(npw_tiny)
 integer,intent(in) :: kg_diel(3,npwdiel),ngfft(18)
 real(dp),intent(in) :: gsq(npwdiel),rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(inout) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),intent(out) :: susd_data(npwdiel,3)
 real(dp),pointer :: khxc(:,:,:,:,:)

!Local variables -------------------------------------------------------
!The Perdew-Wang functional (PW92) yields f_c with the expected n^(-2/3)
!behavior when n -> 0.
!ALDA, Perdew-Wang 92.
 character(len = *), parameter :: fmtd = '(a,t13,5(1x,es12.5))'
 character(len = *), parameter :: fmth1 = '(12x,3(1x,i12))'
 character(len = *), parameter :: fmth2 = '(12x,1x,a12,3(1x,i12))'
 character(len = *), parameter :: fmth3 = '(12x,1x,a12,1x,a12,3(1x,i12))'
!scalars
 integer,parameter :: ixc=7
 integer :: ii,ikernel,ikxc,ipw,ipw1,ipw2,isp,isp1,isp2,jj,nkxc,nlambda
! integer :: npwtestmax
 real(dp) :: deccc,energy_nozero,lambda,lambda_step
! real(dp) :: di,dimax,dr,drmax,dummy
 logical :: khxcalloc,kxcpassed,kxcsaved
! logical :: testflag
 character(len=500) :: message
!arrays
 integer :: ispxc(3)
 real(dp) :: trsusmat(2),tsec(2)
 real(dp),allocatable :: deccc_quad(:),decuc(:),decuc_quad(:,:),energy(:)
 real(dp),allocatable :: krpa(:),kxcdsave(:,:,:),kxcg(:,:,:),kxcgdiff(:,:)
 real(dp),allocatable :: kxcgold(:,:,:),lambda_quad(:),rhor_lambda(:,:)
 real(dp),allocatable :: sus_gabs(:),sus_gavg(:),sus_gdir(:,:),susd_isc(:)
 real(dp),allocatable :: susd_nks(:),susd_tmp(:),weight_quad(:)
 real(dp),pointer :: khxcdiag(:),kxcsave(:,:,:,:,:)
!no_abirules

!***********************************************************************

 call timab(96,1,tsec)

!DEBUG
!write(std_out,*)' acfd_dyson : enter '
!call flush(6)
!ENDDEBUG

!Check input parameters.

 if (nspden > 2) then
   write (message,'(4a)') ch10,&
&   ' acfd_dyson: ERROR - ',ch10,&
&   '  acfd_dyson does not work yet for nspden > 2.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!ikernel = 0 if the exchange-correlation kernel is spin-independent and
!diagonal in reciprocal space (case of the RPA).
!ikernel = 1 otherwise.

!kxcpassed = .true. if the kernel (or part of it) has been precomputed
!once for all outside acfd_dyson and is passed through khxc (case of
!the PGG, BPG and EOK1 kernels).

!kxcsaved = .true. if the whole kernel passed trough khxc has to be saved
!in kxcsave (case of the BPG kernel).

!nkxc =  number of spin channels for the local (part of the) kernels
!(ALDA, BPG, EOK2).

 nkxc = 0
 kxcpassed = .false.
 kxcsaved = .false.
 select case (ikhxc)
   case (ikhxc_NULL, ikhxc_RPA)
     ikernel = 0
   case (ikhxc_ALDA, ikhxc_EOK2)
     ikernel = 1+suskxcrs
     nkxc = 2*nspden-1
   case (ikhxc_PGG, ikhxc_EOK1)
     ikernel = 1+suskxcrs
     kxcpassed = .true.
   case (ikhxc_BPG)
     ikernel = 1+suskxcrs
     kxcpassed = .true.
     kxcsaved  = .true.
     nkxc = 1
     case default
     write (message,'(4a,i10,a)') ch10,&
&     ' acfd_dyson: ERROR - ',ch10,&
&     '  ikhxc = ',ikhxc,' is not a valid xc kernel.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
 end select

!NB : Tests on idyson will be performed later.

!Compute the Hartree kernel with a cut-off in real space.

 ABI_ALLOCATE(krpa,(npwdiel))


!DEBUG
!write(std_out,*)' acfd_dyson : call k_rpa '
!call flush(6)
!ENDDEBUG

 call k_rpa(gsq,krpa,npwdiel,1,rcut_coulomb)

!Eventually save the diagonals of the xc kernel and add
!the RPA kernel...

 if (kxcpassed) then
   ABI_ALLOCATE(kxcdsave,(npwdiel,nspden,nspden))
   do isp2 = 1,nspden
     do isp1 = 1,nspden
       do ipw = 1,npwdiel
         kxcdsave(ipw,isp1,isp2) = khxc(1,ipw,isp1,ipw,isp2)
         khxc(1,ipw,isp1,ipw,isp2) = khxc(1,ipw,isp1,ipw,isp2)+krpa(ipw)
       end do
     end do
   end do
 end if

!...then save the xc kernel itself if needed.

 if (kxcsaved) then
   kxcsave => khxc
   nullify(khxc)
 end if

!Allocate memory for the kernel.

 if (ikernel == 0)  then
   ABI_ALLOCATE(khxcdiag,(npwdiel))
 end if

 khxcalloc = (ikernel == 1).and.((.not.kxcpassed).or.kxcsaved)
 if (khxcalloc)  then
   ABI_ALLOCATE(khxc,(2,npwdiel,nspden,npwdiel,nspden))
 end if

 if (nkxc > 0) then
   ABI_ALLOCATE(kxcgold,(2,nfft,nkxc))
   kxcgold(:,:,:) = 0._dp
 end if

!Set-up abscissas and weights for the coupling-strength integration.

 nlambda = max(0,ndyson)+2

 ABI_ALLOCATE(lambda_quad,(nlambda))
 ABI_ALLOCATE(weight_quad,(nlambda))

 call getlambda(idyson,lambda_quad,nlambda,weight_quad)

!Allocate memory.

 ABI_ALLOCATE(susd_nks,(npwdiel))
 ABI_ALLOCATE(susd_isc,(npwdiel))

 if (option > 0) then
   ABI_ALLOCATE(susd_tmp,(npwdiel))
   ABI_ALLOCATE(sus_gabs,(npw_tiny))
   ABI_ALLOCATE(sus_gavg,(npw_tiny))
   ABI_ALLOCATE(sus_gdir,(npw_tiny,3))
 end if

 if (option > 1) then
   ABI_ALLOCATE(energy,(npw_tiny))
   ABI_ALLOCATE(decuc,(npw_tiny))
   ABI_ALLOCATE(decuc_quad,(npw_tiny,nlambda))
   ABI_ALLOCATE(deccc_quad,(nlambda))
 end if

!Save the diagonal of the Kohn-Sham susceptibility matrix.

 susd_nks(:) = 0._dp

 do isp2 = 1,nspden
   do isp1 = 1,nspden
     do ipw = 1,npwdiel
       susd_nks(ipw) = susd_nks(ipw)+susmat(1,ipw,isp1,ipw,isp2)
     end do
   end do
 end do

 susd_data(:,:) = 0._dp

!Evaluate the diagonal of the first-order approximation to the susceptibility
!matrix [see Lein, Dobson and Gross, J. Comput. Chem. 20, 12 (1999)].

 if ((ldgapp > 0).and.(ikhxc /= ikhxc_NULL)) then

!  The first-order approximation is only implemented for linear kernels.

   if ((ikhxc /= ikhxc_RPA).and.(ikhxc /= ikhxc_PGG).and.(ikhxc /= ikhxc_EOK1)) then
     write (message,'(6a)') ch10,&
&     ' acfd_dyson: ERROR - ',ch10,&
&     '  The first-order approximation to the susceptibility matrix is',ch10,&
&     '  only implemented for the RPA, PGG, and EOK1 kernels at the present stage.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   write (message,'(2a)') ch10, &
&   ' acfd_dyson: Evaluating the first-order approximation to the susceptibility matrix...'
   call wrtout(std_out,message,'COLL')

   if (ikhxc == ikhxc_RPA) khxcdiag(:) = krpa(:)

   call dyson_gl(ikernel,khxcdiag,khxc,npwdiel,nspden,susd_isc,susmat)

!  Coupling-constant integration on the diagonal.

   susd_data(:,2) = 0.5_dp*susd_isc(:)

   if (.false..and.(option > 0)) then

     write (message,fmth2) 'frequency',(jj,jj = 0,npw_tiny-1)
     call wrtout(std_out,message,'COLL')

!    Get the (first-order) polarizabilities by constant, linear and parabolic extrapolation.

     susd_tmp(:) = susd_nks(:)+susd_isc(:)
     call get_susd_null(ig_tiny,igsq_tiny,gsq,npwdiel,npw_tiny,sus_gabs,sus_gavg,sus_gdir,susd_tmp)

!    write (message,fmtd) &
!    &  ' alphaLDG_sq',freq,(-ucvol*sus_gabs(jj),jj = 1,npw_tiny)
!    call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alphaLDG_xx',freq,(-ucvol*sus_gdir(jj,1),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alphaLDG_yy',freq,(-ucvol*sus_gdir(jj,2),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alphaLDG_zz',freq,(-ucvol*sus_gdir(jj,3),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alphaLDG_av',freq,(-ucvol*sus_gavg(jj),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

   end if

!  Calculate the contribution to the correlation energy for the present frequency.

   if (.false..and.(option > 1)) then

!    Using the bare Coulomb interaction:

     call geteexc_uc(energy,energy_nozero,gsq,ig_tiny,npwdiel,npw_tiny,susd_isc)
     decuc(:) = -0.5_dp*energy(:)/two_pi !factor -1/(2*pi) from fluctuation-dissipation theorem.

!    Using the cut-off Coulomb interaction:

     call geteexc_cc(energy,energy_nozero,gsq,npwdiel,npw_tiny,rcut_coulomb,susd_isc)
     deccc = -0.5_dp*energy(1)/two_pi

     write (message,fmtd) &
&     ' dEcLDG_bare',freq,decuc(1:npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' dEcLDG_cut',freq,deccc
     call wrtout(std_out,message,'COLL')

   end if

 end if

!Now solve the Dyson equation along the adiabatic connection path.

 write (message,'(2a)') ch10,&
& ' acfd_dyson: Solving the Dyson equation along the adiabatic connection path...'
 call wrtout(std_out,message,'COLL')

!DEBUG
!write(std_out,*)' acfd_dyson : will solve the Dyson equation'
!call flush(6)
!ENDDEBUG


 susd_isc(:) = 0._dp

 do ii = 1,nlambda

   if ((ii > 1).and.(ikhxc /= ikhxc_NULL)) then

!    Set up the kernel.

     lambda = lambda_quad(ii)
     lambda_step = lambda_quad(ii)-lambda_quad(ii-1)

     if(suskxcrs==0)then

       select case (ikhxc)

         case (ikhxc_RPA)

           if (idyson /= idyson_SC) then
             khxcdiag(:) = krpa(:)*lambda_step
           else
             khxcdiag(:) = krpa(:)*lambda
           end if

         case (ikhxc_PGG, ikhxc_EOK1)

           khxc(:,:,:,:,:) = khxc(:,:,:,:,:)*lambda_step

         case (ikhxc_ALDA, ikhxc_BPG, ikhxc_EOK2)

           ABI_ALLOCATE(kxcg,(2,nfft,nkxc))
           ABI_ALLOCATE(kxcgdiff,(2,nfft))
           ABI_ALLOCATE(rhor_lambda,(nfft,nspden))

           rhor_lambda(:,:) = rhor(:,:)/lambda**3

           select case (ikhxc)
             case (ikhxc_ALDA)
               call kxc_alda(dtset,ixc,kxcg,mpi_enreg,nfft,ngfft,nspden,1,rhor_lambda,rhocut,rprimd)
               ispxc = (/1,2,3/)
             case (ikhxc_BPG)
               call kxc_alda(dtset,ixc,kxcg,mpi_enreg,nfft,ngfft,nspden,2,rhor_lambda,rhocut,rprimd)
               ispxc(1) = nspden
             case (ikhxc_EOK2)
               call kxc_eok(2,kxcg,mpi_enreg,nfft,ngfft,nspden,dtset%paral_kgb,rhor_lambda,rhocut)
               ispxc = (/1,2,3/)
           end select

           kxcg(:,:,:) = kxcg(:,:,:)/lambda

           do ikxc = 1,nkxc
             kxcgdiff(:,:) = kxcg(:,:,ikxc)-kxcgold(:,:,ikxc)
             call klocal(ispxc(ikxc),kg_diel,khxc,kxcgdiff,nfft,ngfft,npwdiel,nspden,1)
           end do

           if (ikhxc == ikhxc_BPG) then
             khxc = khxc+kxcsave*lambda_step
           else
             do isp2 = 1,nspden
               do isp1 = 1,nspden
                 do ipw = 1,npwdiel
                   khxc(1,ipw,isp1,ipw,isp2) = khxc(1,ipw,isp1,ipw,isp2)+krpa(ipw)*lambda_step
                 end do
               end do
             end do
           end if

           kxcgold(:,:,:) = kxcg(:,:,:)

           ABI_DEALLOCATE(kxcg)
           ABI_DEALLOCATE(kxcgdiff)
           ABI_DEALLOCATE(rhor_lambda)

       end select

     else if(suskxcrs==1)then

       if(ikhxc_ALDA/=ikhxc)then
         write(std_out,*)' acfd_dyson : suskxcrs==1 presently implemented only for ikhxc=ALDA '
         stop
       end if
       if(nspden/=1)then
         write(std_out,*)' acfd_dyson : suskxcrs==1 presently implemented only for nspden=1 '
         stop
       end if
       khxcdiag(:) = krpa(:)*lambda

!      rhor_lambda(:,:) = rhor(:,:)/lambda**3

       write(std_out,*)' not yet implemented '
       stop

       ABI_DEALLOCATE(rhor_lambda)
     end if

!    Solve the Dyson equation according to idyson.

     select case (idyson)

       case (idyson_LS)

!        Solve the Dyson equation as a linear system.

         call dyson_ls(ikernel,khxcdiag,khxc,npwdiel,nspden,susmat)

       case (idyson_DE)

!        Solve the Dyson equation as a first-order differential equation.
!        WARNING : Only implemented for the RPA at the present stage.

         if ((ikhxc /= ikhxc_RPA).and.(ikhxc /= ikhxc_PGG)) then
           write (message,'(6a)') ch10,&
&           ' acfd_dyson: ERROR - ',ch10,&
&           '  The solution of the Dyson equation as a first-order differential equation',ch10,&
&           '  is only implemented for the RPA and PGG kernels at the present stage.'
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if

         call dyson_de(ikernel,khxcdiag,khxc,npwdiel,nspden,susmat)

       case (idyson_SC)

!        Solve the Dyson equation as a self-consistent problem.
!        WARNING : Only implemented for the RPA at the present stage.

         if (ikhxc /= ikhxc_RPA) then
           write (message,'(6a)') ch10,&
&           ' acfd_dyson: ERROR - ',ch10,&
&           '  The solution of the Dyson equation as a self-consistent problem',ch10,&
&           '  is only implemented for the RPA at the present stage.'
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')
           call leave_new('COLL')
         end if

         call dyson_sc(khxcdiag,npwdiel,nspden,susd_isc,susmat)

         case default

         write (message,'(4a,i10,a)') ch10,&
&         ' acfd_dyson: ERROR - ',ch10,&
&         '  idyson = ',idyson,' is not a valid method to solve the Dyson equation.'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
         call leave_new('COLL')

     end select

!    Restore the xc kernel if needed.

     if ((ikhxc == ikhxc_PGG).or.(ikhxc == ikhxc_EOK1)) &
&     khxc(:,:,:,:,:) = khxc(:,:,:,:,:)/lambda_step

   end if

!  Compute the difference between the diagonals of the interacting
!  and Kohn-Sham susceptibility matrices.

   if (idyson /= idyson_SC) then

     susd_isc(:) = -susd_nks(:)

     do isp2 = 1,nspden
       do isp1 = 1,nspden
         do ipw = 1,npwdiel
           susd_isc(ipw) = susd_isc(ipw)+susmat(1,ipw,isp1,ipw,isp2)
         end do
       end do
     end do

   end if

!  Coupling-constant integration on the diagonal.

   susd_data(:,1) = susd_data(:,1)+weight_quad(ii)*susd_isc(:)

   if (option > 0) then

!    Get the polarizabilities by constant, linear and parabolic extrapolation.

     susd_tmp(:) = susd_isc(:)+susd_nks(:)
     call get_susd_null(ig_tiny,igsq_tiny,gsq,npwdiel,npw_tiny,sus_gabs,sus_gavg,sus_gdir,susd_tmp)

     write (message,fmth3) 'lambda','frequency',(jj,jj = 0,npw_tiny-1)
     call wrtout(std_out,message,'COLL')

!    write (message,fmtd) &
!    &  ' alpha_sq',lambda_quad(ii),freq,(-ucvol*sus_gabs(jj),jj = 1,npw_tiny)
!    call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alpha_xx',lambda_quad(ii),freq,(-ucvol*sus_gdir(jj,1),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alpha_yy',lambda_quad(ii),freq,(-ucvol*sus_gdir(jj,2),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alpha_zz',lambda_quad(ii),freq,(-ucvol*sus_gdir(jj,3),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alpha_av',lambda_quad(ii),freq,(-ucvol*sus_gavg(jj),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

   end if

   if (option > 1) then

!    Calculate the contribution to the correlation energy for the present lambda and frequency.
!    Using the bare Coulomb interaction:

     call geteexc_uc(energy,energy_nozero,gsq,ig_tiny,npwdiel,npw_tiny,susd_isc)
     decuc_quad(:,ii) = -energy(:)/two_pi !factor -1/(2*pi) from fluctuation-dissipation theorem.

!    Using the cut-off Coulomb interaction:

     call geteexc_cc(energy,energy_nozero,gsq,npwdiel,npw_tiny,rcut_coulomb,susd_isc)
     deccc_quad(ii) = -energy(1)/two_pi

     write (message,fmtd) &
&     ' d2Ec_bare',lambda_quad(ii),freq,decuc_quad(1:npw_tiny,ii)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' d2Ec_cut',lambda_quad(ii),freq,deccc_quad(ii)
     call wrtout(std_out,message,'COLL')

   end if

 end do

!Coupling-constant integration for the correlation energy.
!NOTE: This is now done in xcacfd.f.

 if (.false..and.(option > 1)) then

   write (message,'(2a)') ch10,&
&   ' acfd_dyson: Performing the integration over the coupling constant...'
   call wrtout(std_out,message,'COLL')

   decuc(:) = 0._dp
   deccc    = 0._dp

   do ii = 1,nlambda
     decuc(:) = decuc(:)+weight_quad(ii)*decuc_quad(1,ii)
     deccc    = deccc   +weight_quad(ii)*deccc_quad(ii)
   end do

   write (message,fmth2) &
&   'frequency',(jj,jj = 0,npw_tiny-1)
   call wrtout(std_out,message,'COLL')

   write (message,fmtd) &
&   ' dEc_bare',freq,decuc(1:npw_tiny)
   call wrtout(std_out,message,'COLL')

   write (message,fmtd) &
&   ' dEc_cut',freq,deccc
   call wrtout(std_out,message,'COLL')

 end if

!Calculate the trace of the interacting susceptibility matrix
!and store the diagonal of susmat in susd_data(:,3).

 trsusmat(:) = 0._dp

 do isp = 1,nspden
   do ipw = 1,npwdiel
     trsusmat(:) = trsusmat(:)+susmat(:,ipw,isp,ipw,isp)
     susd_data(ipw,3) = susd_data(ipw,3)+susmat(1,ipw,isp,ipw,isp)
   end do
 end do

!write (message,'(a12,1x,es12.5,1x,es12.5,a,es12.5)') &
!& ' Tr(susmat) ',freq,trsusmat(1),'+i',trsusmat(2)
!call wrtout(std_out,message,'COLL')

!Free memory.

 if (ikernel == 0)  then
   ABI_DEALLOCATE(khxcdiag)
 end if
 if (khxcalloc)  then
   ABI_DEALLOCATE(khxc)
 end if

 if (nkxc > 0)  then
   ABI_DEALLOCATE(kxcgold)
 end if

 ABI_DEALLOCATE(krpa)
 ABI_DEALLOCATE(lambda_quad)
 ABI_DEALLOCATE(weight_quad)
 ABI_DEALLOCATE(susd_nks)
 ABI_DEALLOCATE(susd_isc)

 if (option > 0) then
   ABI_DEALLOCATE(susd_tmp)
   ABI_DEALLOCATE(sus_gabs)
   ABI_DEALLOCATE(sus_gavg)
   ABI_DEALLOCATE(sus_gdir)
 end if

 if (option > 1) then
   ABI_DEALLOCATE(energy)
   ABI_DEALLOCATE(decuc)
   ABI_DEALLOCATE(decuc_quad)
   ABI_DEALLOCATE(deccc_quad)
 end if

!Restore the xc kernel itself if needed...

 if (kxcsaved) khxc => kxcsave

!...including its diagonals.

 if (kxcpassed) then
   do isp2 = 1,nspden
     do isp1 = 1,nspden
       do ipw = 1,npwdiel
         khxc(1,ipw,isp1,ipw,isp2) = kxcdsave(ipw,isp1,isp2)
       end do
     end do
   end do
   ABI_DEALLOCATE(kxcdsave)
 end if

!DEBUG
!Check that susmat is hermitian.
!The test is limited to the first npwtestmax planewaves of the up-up spin channel.
!npwtestmax = 100
!write (std_out,'(a)') ' acfd_dyson: Check that susmat is hermitian.'
!testflag = .true.
!drmax = 0._dp; dimax = 0._dp
!do ipw2 = 1,min(npwtestmax,npwdiel)
!do ipw1 = 1,ipw2
!dr = abs(susmat(1,ipw1,1,ipw2,1)-susmat(1,ipw2,1,ipw1,1))
!di = abs(susmat(2,ipw1,1,ipw2,1)+susmat(2,ipw2,1,ipw1,1))
!if ((dr > tol10).or.(di > tol10)) then
!if (.false.) write (std_out,'(2(1x,i4),4(1x,es12.5))') &
!&    ipw1,ipw2,susmat(1,ipw1,1,ipw2,1),dr,susmat(2,ipw1,1,ipw2,1),di
!drmax = max(drmax,dr); dimax = max(dimax,di)
!testflag = .false.
!end if
!end do
!end do
!if (testflag) then
!write (std_out,'(a)') ' - Test passed - '
!else
!write (std_out,'(a)') ' - Test failed - '
!write (message,'(6a,2(1x,es12.5))') ch10,&
!&  ' acfd_dyson: WARNING - ',ch10,&
!&  '  The interacting susceptibility matrix is not hermitian.',ch10,&
!&  '  Max real/imag diffs:',drmax,dimax
!call wrtout(std_out,message,'COLL')
!end if
!write (std_out,'(a)') ' dietcel: Diagonal elements of the susceptibility matrix:'
!do ipw = 1,npwdiel
!write (std_out,'(1x,i4,2(1x,es12.5))') ipw,susmat(1,ipw,1,ipw,1),susmat(2,ipw,1,ipw,1)
!end do
!ENDDEBUG

!This test is enforced at the present time for some kernels but should be removed at last...
!!DEBUG
!Diagonalize susmat and print the eigenvalues.
 if (((ikhxc == ikhxc_ALDA).or.(ikhxc == ikhxc_BPG).or.&
& (ikhxc == ikhxc_EOK1).or.(ikhxc == ikhxc_EOK2)).and.(nspden == 1)) then

   ABI_ALLOCATE(khxcdiag,(npwdiel))
   ABI_ALLOCATE(energy,(npwdiel))
   ABI_ALLOCATE(susd_tmp,(4*npwdiel))
   ABI_ALLOCATE(sus_gavg,(3*npwdiel))
!  Save the diagonal of susmat.
   do ipw = 1,npwdiel
     khxcdiag(ipw) = susmat(1,ipw,1,ipw,1)
   end do
   call zheev('n','l',npwdiel,susmat,npwdiel,energy,susd_tmp,2*npwdiel,sus_gavg,jj)
   if (jj /= 0) then
     write (std_out,'(2a,i10)') ch10,' acfd_dyson: Error diagonalizing susmat, info = ',jj
   else
     write (std_out,'(2a)') ch10,' acfd_dyson: susmat eigenvalues (should all be negative):'
     do ipw = 1,min(npwdiel,25)
       write (std_out,'(1x,i4,a,es12.5)') ipw,' = ',energy(ipw)
     end do
     write (std_out,'(a)') ' [...]'
     do ipw = max(npwdiel-24,1),npwdiel
       write (std_out,'(1x,i4,a,es12.5)') ipw,' = ',energy(ipw)
     end do
   end if
!  Restore the lower triangle (including the diagonal) of susmat.
   do ipw = 1,npwdiel
     susmat(1,ipw,1,ipw,1) = khxcdiag(ipw)
     susmat(2,ipw,1,ipw,1) = 0._dp
   end do
   do ipw2 = 1,npwdiel
     do ipw1 = ipw2+1,npwdiel
       susmat(1,ipw1,1,ipw2,1) =  susmat(1,ipw2,1,ipw1,1)
       susmat(2,ipw1,1,ipw2,1) = -susmat(2,ipw2,1,ipw1,1)
     end do
   end do
   ABI_DEALLOCATE(khxcdiag)
   ABI_DEALLOCATE(energy)
   ABI_DEALLOCATE(susd_tmp)
   ABI_DEALLOCATE(sus_gavg)

 end if
!!ENDDEBUG

 call timab(96,2,tsec)

!DEBUG
!write(std_out,*)' acfd_dyson : exit '
!call flush(6)
!ENDDEBUG

!Replace directly in the routine to avoid a contains (TD)
!contains

!function polarizability(dummy)

!real(dp),intent(in) :: dummy
!real(dp) :: polarizability

!polarizability = -dummy*ucvol

!end function polarizability

end subroutine acfd_dyson

!!***
