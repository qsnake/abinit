!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcacfd
!! NAME
!! xcacfd
!!
!! FUNCTION
!! Compute the exchange-correlation energy from the imaginary-frequency
!! susceptibility matrix (using the adibatic-connection fluctuation-dissipation
!! theorem for the pair-correlations), from input wavefunctions, band occupations,
!! and k point weights. Also computes the exchange energy directly from the
!! the square modulus of the density matrix.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT group (XG,MF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dielar(7)  = input parameters for dielectric matrix and susceptibility:
!!               diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dtfil <type(datafiles_type)> = variables related to files.
!!  dtset <type(dataset_type)> = all input variables in this dataset.
!!  eigen(mband*nkpt*nsppol)  = array for holding eigenvalues (hartree).
!!  freq(nfreqsus) = array for frequencies (hartree).
!!  gbound_diel(2*mgfftdiel+8,2) =  sphere boundary for the dielectric matrix.
!!  gmet(3,3) = reciprocal space metric (bohr**-2).
!!  gprimd(3,3) = dimensional primitive translations for reciprocal space (bohr**-1).
!!  irrzondiel(nfftdiel**(1-1/nsym),2+(nspden/4),(nspden/nsppol)-3*(nspden/4)) = irreducible zone data.
!!  kg(3,mpw*mkmem) = reduced planewave coordinates.
!!  kg_diel(3,npwdiel) = reduced planewave coordinates for the dielectric matrix.
!!  mband = maximum number of bands.
!!  mgfftdiel = maximum size of 1D FFTs, for the computation of the dielectric matrix.
!!  mkmem = maximum number of k points in core memory.
!!  mpw = maximum allowed value for npw.
!!  nfft = number of fft grid points.
!!  nfftdiel = number of fft grid points for the computation of the diel matrix.
!!  nfreqsus = size of frequency grid.
!!  ngfft(1:3), ngfftdiel(1:3) = fft box dimensions,
!!                               see getng for ngfft(1:8) and ngfftdiel(4:8).
!!  nkpt = number of k points.
!!  npwarr(nkpt) = number of planewaves and boundary planewaves at each k,
!!                 for going from the WF sphere to the medium size FFT grid.
!!  npwdiel = third and fifth dimension of the susmat_dyn array.
!!  nspden = number of spin-density components.
!!  nspinor = number of spinorial components of the wavefunctions.
!!  nsppol = 1 for unpolarized, 2 for spin-polarized.
!!  nsym  =number of symmetry elements in group (at least 1 for identity).
!!  occ(mband*nkpt*nsppol) = occupation numbers for each band (usually 2.0)
!!                           at each k point.
!!  phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4)) = nonsymmorphic translation phases.
!!  rhor(nfft,nspden) = density in real space.
!!  rprimd(3,3) = dimensional primitive translations for real space (bohr).
!!  ucvol = unit cell volume (Bohr**3).
!!  wff1= struct info for current wf disk file
!!  wght_freq(nfreqsus) = integration weight for frequency integration.
!!
!! OUTPUT
!!  (print only)
!!
!! SIDE EFFECTS
!!  mpi_enreg = informations about MPI parallelization.
!!
!! WARNINGS
!! Restrictions (MF):
!!  a - Argument occopt >= 3 is not allowed, since the contributions due to
!!      the Fermi level change are not implemented.
!!  b - Direct computation of exchange energy meant for closed shell cases only.
!!  c - PGG and ALDA kernels are entirely experimental.
!!
!! PARENTS
!!      suscep
!!
!! CHILDREN
!!      acfd_dyson,acfd_intexact,get_g_tiny,get_susd_null,geteexc_cc,geteexc_uc
!!      klocal,kxc_eok,kxc_pgg,leave_new,prtsusd,suscep_dyn,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xcacfd(dielar,&
&  dtfil,dtset,eigen,freq,gbound_diel,gmet,&
&  gprimd,irrzondiel,kg,kg_diel,mband,mgfftdiel,&
&  mkmem,mpi_enreg,mpw,nfft,nfftdiel,nfreqsus,ngfft,&
&  ngfftdiel,nkpt,npwarr,npwdiel,nspden,nspinor,nsppol,&
&  nsym,occ,phnonsdiel,rhor,rprimd,ucvol,wff1,wght_freq)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters
 use m_wffile

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcacfd'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_77_suscep, except_this_one => xcacfd
!End of the abilint section

 implicit none

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: mband,mgfftdiel,mkmem,mpw,nfft,nfftdiel,nfreqsus
 integer,intent(in) :: nkpt,npwdiel,nspden,nsppol,nsym
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(wffile_type),intent(inout) :: wff1
!arrays
 integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2+(nspden/4),(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: kg(3,mpw*mkmem),kg_diel(3,npwdiel),ngfft(18)
 integer,intent(in) :: ngfftdiel(18),npwarr(nkpt)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: freq(nfreqsus),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: wght_freq(nfreqsus)
 real(dp),intent(inout) :: dielar(7)

!Local variables -------------------------------------------------------
!Number of frequencies kept in core memory.
!Number of shortest G vectors.
 character(len = *), parameter :: fmtd = '(a,t13,5(1x,es12.5))'
 character(len = *), parameter :: fmth1 = '(12x,3(1x,i12))'
 character(len = *), parameter :: fmth2 = '(12x,1x,a12,3(1x,i12))'
 character(len = *), parameter :: fmth3 = '(12x,1x,a12,1x,a12,3(1x,i12))'
 character(len = *), parameter :: fmtoi = '(a,t21,1x,i15,1x,a)'
 character(len = *), parameter :: fmtod = '(a,t21,1x,es15.8,1x,a)'
 character(len = *), parameter :: fmtos = '(a,t21,1x,a15,1x,a)'
 character(len = *), parameter :: fmtog = '(a,"(",i1,")",t21,1x,es15.8,1x,a,i1,a)'
!scalars
 integer,parameter :: nfreq_mem=1,npw_tiny=3
 integer :: idyson,ifreq,ifreq_start,ifreq_stop,ii,ikhxc,intexact
 integer :: ipw,ipw1,ipw2,isp,isp1,isp2,jfreq,jj,ldgapp,nbandsus
 integer :: ndyson,susopt
 real(dp) :: alpha_int,alpha_int0,alpha_non,alpha_non0,decs_cc,dummy,ecgl_cc
 real(dp) :: ecs_cc,energy_raw,ex_cc,ex_dm_cc,ex_dm_cc_nozerog,ex_dm_uc_nozerog
 real(dp) :: n_test,rcut_coulomb,rhocut
 character(len=500) :: message
!arrays
 integer :: ig_tiny(npw_tiny,3),igsq_tiny(npw_tiny)
 integer,allocatable :: index_g(:)
 real(dp) :: ecgl_uc(npw_tiny),ecs_uc(npw_tiny),energy(npw_tiny)
 real(dp) :: ex_dm_uc(npw_tiny),sus_gabs(npw_tiny),sus_gavg(npw_tiny)
 real(dp) :: sus_gdir(npw_tiny,3)
 real(dp),allocatable :: freq_mem(:),gsq(:),kxcg(:,:,:),susd_aux_dyn(:,:,:,:)
 real(dp),allocatable :: susd_data(:,:),susd_tmp(:),susmat_dyn(:,:,:,:,:,:)
 real(dp),allocatable :: wght_mem(:)
 real(dp),pointer :: khxc(:,:,:,:,:)

!***********************************************************************

!DEBUG
!write (std_out,*) ' xcacfd: enter'
!!call flush(6)
!write (std_out,*) ' nfft    : ',nfft
!write (std_out,*) ' nfreqsus: ',nfreqsus
!write (std_out,*) ' nspden  : ',nspden
!write (std_out,*) ' freq    : ',freq(:)
!!call flush(6)
!ENDDEBUG

!Initialize some variables.

 nullify(khxc)

 idyson = dtset%idyson
 ikhxc = dtset%ikhxc
 intexact = dtset%intexact
 ldgapp = dtset%ldgapp
 nbandsus = dtset%nbandsus
 ndyson = dtset%ndyson

 ifreq_start = 1
 ifreq_stop  = nfreqsus

!This is EXPERIMENTAL !!

!if ((dtset%userre > 0._dp).and.(dtset%userre < 1._dp)) then
!rhocut = dtset%userre
!else
 rhocut = 0.01_dp
!end if

!Check input parameters.

 if (nspden > 1) then
   write (message,'(4a)') ch10,&
&   ' xcacfd: ERROR - ',ch10,&
&   '  xcacfd does not work yet for nspden > 1.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 select case (ikhxc)
   case (ikhxc_NULL)
   case (ikhxc_RPA)
   case (ikhxc_ALDA)
   case (ikhxc_PGG)
   case (ikhxc_BPG)
   case (ikhxc_EOK1)
   case (ikhxc_EOK2)
     case default
     write (message,'(4a,i10,a)') ch10,&
&     ' xcacfd: ERROR - ',ch10,&
&     '  ikhxc = ',ikhxc,' is not a valid xc kernel.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
 end select

 write (message,'(2a,i2,3a)') ch10,&
& ' xcacfd: Enter with ikhxc = ',ikhxc,' (',trim(ikhxc_name(ikhxc)),').'
 call wrtout(std_out,message,'COLL')

!Number of bands to be used in the calculation of the Kohn-Sham susceptibility matrix.

 if (nbandsus > 0) then
   nbandsus = min(nbandsus,minval(dtset%nband(:)))
 else
   nbandsus = minval(dtset%nband(:))
 end if
 write (message,'(a,i5,a)') &
& ' xcacfd: Use nbandsus = ',nbandsus,' bands to calculate the Kohn-Sham susceptibility matrix.'
 call wrtout(std_out,message,'COLL')

!The integration over the coupling constant can be performed analytically in the RPA
!and for PGG (in spin-compensated, two-electron systems). Moreover, the interacting
!susceptibility matrix need not be calculated in that case, which makes the calculation
!much faster.

 if (intexact > 0) then

   if ((ikhxc /= ikhxc_RPA).and.(ikhxc /= ikhxc_PGG)) then
     write (message,'(6a)') ch10,&
&     ' xcacfd: ERROR - ',ch10,&
&     '  The exact integration over the coupling constant can only be',ch10,&
&     '  performed for the RPA and spin-compensated two-electron PGG kernels.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   write (message,'(6a)') ch10,&
&   ' xcacfd: The integration over the coupling constant will be performed analytically,',ch10,&
&   ' using a Coulomb interaction with a cut-off in real space. The interacting',ch10,&
&   ' susceptiblity matrix, which is not needed in that case, will not be calculated.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

 else

!  Solution mode for the Dyson equation.

   select case (idyson)
     case (idyson_LS) !as a linear system.
       if (ndyson < 0) ndyson = 3
     case (idyson_DE) !as a differential equation.
       if (ndyson < 0) ndyson = 9
     case (idyson_SC) !as a self-consistent problem.
       if (ndyson < 0) ndyson = 3
       case default
       write (message,'(4a,i10,a)') ch10,&
&       ' xcacfd: ERROR - ',ch10,&
&       '  idyson = ',idyson,' is not a valid method to solve the Dyson equation.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
   end select

   write (message,'(a,i1,5a,i2,a)') &
&   ' xcacfd: Use idyson = ',idyson,' (solve the Dyson equation as a ',&
&   trim(idyson_name(idyson)),').',ch10,&
&   ' xcacfd: Use ndyson = ',ndyson,' points (excl. 0,1) for coupling constant integration.'
   call wrtout(std_out,message,'COLL')

 end if

!Print a few messages describing the output.

 write (message,'(8a)') ch10,&
& ' xcacfd:',ch10,&
& ' Nomenclature:',ch10,&
& '  alpha[0]_xx, alpha[0]_yy, alpha[0]_zz: diagonal elements of the polarizability matrix (*).',ch10,&
& '  alpha[0]_av: average of alpha[0]_xx, alpha[0]_yy and alpha[0]_zz (*).'
 call wrtout(std_out,message,'COLL')
 write (message,'(5a)') &
& '  d2Ec_bare: contribution to the correlation energy for a given lambda and frequency,',ch10,&
& '             using the bare Coulomb interaction (*).',ch10,&
& '  d2Ec_cut:  idem, using a Coulomb interaction with a cut-off in real space.'
 if (intexact <= 0) call wrtout(std_out,message,'COLL')
 write (message,'(5a)') &
& '  dE[x/c]_bare: contribution to the exchange/correlation energy for a given frequency',ch10,&
& '                (i.e., integrated over lambda), using the bare Coulomb interaction (*).',ch10,&
& '  dE[x/c]_cut:  idem, using a Coulomb interaction with a cut-off in real space.'
 call wrtout(std_out,message,'COLL')
 write (message,'(5a)') &
& '  alphaLDG_*, dEcLDG_*: same as above, but in the Lein, Dobson and Gross approximation.',ch10,&
& ' (*) The different numbers given for these quantities correspond to extrapolations of',ch10,&
& ' increasing order to G = 0.'
 if (intexact <= 0) call wrtout(std_out,message,'COLL')

!Allocate memory.

 ABI_ALLOCATE(gsq,(npwdiel))
 ABI_ALLOCATE(index_g,(npwdiel))
 ABI_ALLOCATE(freq_mem,(nfreq_mem))
 ABI_ALLOCATE(wght_mem,(nfreq_mem))
 ABI_ALLOCATE(susd_tmp,(npwdiel))
 ABI_ALLOCATE(susd_data,(npwdiel,npw_tiny))
 ABI_ALLOCATE(susd_aux_dyn,(2,npwdiel,nspden,nfreq_mem))
 ABI_ALLOCATE(susmat_dyn,(2,npwdiel,nspden,npwdiel,nspden,nfreq_mem))

!Perform initializations.

 ex_cc = 0._dp
 ecs_cc = 0._dp
 ecs_uc(:) = 0._dp
 ecgl_cc = 0._dp
 ecgl_uc(:) = 0._dp
 n_test = 0._dp

!Index the reciprocal space translations and determine the npw_tiny shortest vectors
!along the primitive translations (ig_tiny), and by length (igsq_tiny).

 call get_g_tiny(gmet,gprimd,gsq,ig_tiny,igsq_tiny,index_g,kg_diel,npwdiel,npw_tiny)

!Determine the real space cut-off radius for the Coulomb interaction.

 rcut_coulomb = 1.d20
 do ii = 1,3
   rcut_coulomb = min(rcut_coulomb,sum(rprimd(:,ii)*rprimd(:,ii)))
 end do
 rcut_coulomb = 0.5_dp*sqrt(rcut_coulomb)

!DEBUG
!write (std_out,*) ' xcacfd: rcut_coulomb = ',rcut_coulomb
!ENDDEBUG

!Compute the exchange energy from the density matrix.
!First get the squared modulus of the density matrix.

 if (.true..or.(ikhxc == ikhxc_NULL)) then

   write (message,'(2a)') ch10,&
&   ' xcacfd: Calculating the exchange energy from the density matrix...'
   call wrtout(std_out,message,'COLL')

   susopt = 2 !Density matrix.

   dummy = dielar(6) !dielam.
   dielar(6) = 0._dp

   call suscep_dyn(dielar,dtset,eigen,freq_mem,&
&   gbound_diel,gprimd,irrzondiel,dtset%istwfk,kg,kg_diel,mband,mgfftdiel,mkmem,&
&   mpi_enreg,mpw,dtset%nband,nbandsus,nfftdiel,nfreq_mem,ngfftdiel,nkpt,&
&   npwarr,npwdiel,nspden,nspinor,nsppol,nsym,occ,dtset%occopt,phnonsdiel,&
&   rprimd,susopt,susd_aux_dyn,susmat_dyn,dtset%symafm,dtset%symrel,&
   dtset%tnons,ucvol,dtfil%unkg,wff1,dtset%wtk)

!  DEBUG
!  write (std_out,*) ' xcacfd: suscep_dyn done for density matrix'
!  call flush(6)
!  ENDDEBUG

   dielar(6) = dummy

   susd_tmp(:) = 0._dp
   do isp = 1,nspden
     do ipw = 1,npwdiel
       susd_tmp(ipw) = susd_tmp(ipw)+susmat_dyn(1,ipw,isp,ipw,isp,1)
     end do
   end do

!  Evaluate the exchange integral.
!  Using the cut-off Coulomb interaction:

!  DEBUG
!  write (std_out,*) ' xcacfd: call geteexc_cc'
!  call flush(6)
!  ENDDEBUG

   call geteexc_cc(energy,energy_raw,gsq,npwdiel,npw_tiny,rcut_coulomb,susd_tmp)

!  DEBUG
!  write (std_out,*) ' xcacfd: exit geteexc_cc'
!  call flush(6)
!  ENDDEBUG

   ex_dm_cc = -0.5_dp*energy(1)
   ex_dm_cc_nozerog = -0.5_dp*energy_raw

!  Using the bare Coulomb interaction:

!  DEBUG
!  write (std_out,*) ' xcacfd: call geteexc_uc'
!  call flush(6)
!  ENDDEBUG

   call geteexc_uc(energy,energy_raw,gsq,ig_tiny,npwdiel,npw_tiny,susd_tmp)

!  DEBUG
!  write (std_out,*) ' xcacfd: exit geteexc_uc'
!  call flush(6)
!  ENDDEBUG

   ex_dm_uc_nozerog = -0.5_dp*energy_raw

   susd_tmp(:) = susd_tmp(:)-susd_tmp(1)
   call geteexc_uc(energy,energy_raw,gsq,ig_tiny,npwdiel,npw_tiny,susd_tmp)
   ex_dm_uc(:) = -0.5_dp*energy(:)

   write (message,fmth1) (jj,jj = 0,npw_tiny-1)
   call wrtout(std_out,message,'COLL')

   write (message,fmtd) &
&   ' ExDM_bare',ex_dm_uc(1:npw_tiny)
   call wrtout(std_out,message,'COLL')

   write (message,fmtd) &
&   ' ExDM_cut',ex_dm_cc
   call wrtout(std_out,message,'COLL')

 end if

!DEBUG
!write (std_out,*) ' xcacfd: before PGG'
!call flush(6)
!ENDDEBUG


!Compute the PGG kernel if needed.

 if (((ikhxc == ikhxc_PGG).and.(intexact <= 0)).or.(ikhxc == ikhxc_BPG)) then

   write (message,'(2a)') ch10,&
&   ' xcacfd: Calculating the PGG kernel...'
   call wrtout(std_out,message,'COLL')

   if (nspden > 1) then
     write (message,'(4a)') ch10,&
&     ' xcacfd: ERROR - ',ch10,&
&     '  The PGG and BPG kernels do not work yet for nspden > 1.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

   susopt = 3 !Density-weighted density matrix.

   dummy = dielar(6) !dielam.
   dielar(6) = 0._dp

   call suscep_dyn(dielar,dtset,eigen,freq_mem,&
&   gbound_diel,gprimd,irrzondiel,dtset%istwfk,kg,kg_diel,mband,mgfftdiel,mkmem,&
&   mpi_enreg,mpw,dtset%nband,nbandsus,nfftdiel,nfreq_mem,ngfftdiel,nkpt,&
&   npwarr,npwdiel,nspden,nspinor,nsppol,nsym,occ,dtset%occopt,phnonsdiel,&
&   rprimd,susopt,susd_aux_dyn,susmat_dyn,dtset%symafm,dtset%symrel,&
   dtset%tnons,ucvol,dtfil%unkg,wff1,dtset%wtk)

!  DEBUG
!  write (std_out,*) ' xcacfd: suscep_dyn done for PGG'
!  call flush(6)
!  ENDDEBUG

   dielar(6) = dummy

!  Convolution of the density-weighted density matrix with the cut-off Coulomb interaction.

   ABI_ALLOCATE(khxc,(2,npwdiel,nspden,npwdiel,nspden))
   call kxc_pgg(gmet,kg_diel,khxc,npwdiel,rcut_coulomb,&
&   susmat_dyn(1,1,1,1,1,1),ucvol)

!  There is apparently a factor one-half missing somewhere in the calculation of the
!  spin-restricted PGG kernel...
!  TODO: should be fixed in kxc_pgg or suscep_dyn_pgg...

   khxc(:,:,:,:,:) = 0.5_dp*khxc(:,:,:,:,:)

 end if

!Compute the linear energy optimized kernel (EOK1) if needed.

 if (ikhxc == ikhxc_EOK1) then

   write (message,'(2a)') ch10,&
&   ' xcacfd: Calculating the linear energy optimized kernel (EOK1)...'
   call wrtout(std_out,message,'COLL')

   ABI_ALLOCATE(kxcg,(2,nfft,2*nspden-1))
   ABI_ALLOCATE(khxc,(2,npwdiel,nspden,npwdiel,nspden))

   call kxc_eok(1,kxcg,mpi_enreg,nfft,ngfft,nspden,dtset%paral_kgb,rhor,rhocut)
   do ii = 1,2*nspden-1
     call klocal(ii,kg_diel,khxc,kxcg(:,:,ii),nfft,ngfft,npwdiel,nspden,1)
   end do

   ABI_DEALLOCATE(kxcg)

 end if

!Compute the dynamic susceptibility matrices.
!The following file will be used to store data from prtsusd.

 open (unit = tmp_unit,file = dtfil%fnametmp_sustr,status = 'unknown')

!Big loop over frequencies.

 do ifreq = ifreq_start,ifreq_stop,nfreq_mem

   freq_mem(1:nfreq_mem) = freq(ifreq:ifreq+nfreq_mem-1)
   wght_mem(1:nfreq_mem) = wght_freq(ifreq:ifreq+nfreq_mem-1)

!  Compute the Kohn-Sham susceptibility matrices.

   write (message,'(2a)') ch10,&
&   ' xcacfd: Calculating the Kohn-Sham susceptibility matrices...'
   call wrtout(std_out,message,'COLL')

   susopt = 1

   call suscep_dyn(dielar,dtset,eigen,freq_mem,&
&   gbound_diel,gprimd,irrzondiel,dtset%istwfk,kg,kg_diel,mband,mgfftdiel,mkmem,&
&   mpi_enreg,mpw,dtset%nband,nbandsus,nfftdiel,nfreq_mem,ngfftdiel,nkpt,&
&   npwarr,npwdiel,nspden,nspinor,nsppol,nsym,occ,dtset%occopt,phnonsdiel,&
&   rprimd,susopt,susd_aux_dyn,susmat_dyn,dtset%symafm,dtset%symrel,&
   dtset%tnons,ucvol,dtfil%unkg,wff1,dtset%wtk)

!  DEBUG
!  write(std_out,*) ' xcacfd: suscep_dyn done 3'
!  call flush(6)
!  ENDDEBUG

!  Print the Kohn-Sham susceptibility matrices.

   if (((nfreqsus == 1).or.(ifreq == min(10,nfreqsus))).and.(ikhxc == ikhxc_NULL)) then
     do jfreq = 1,nfreq_mem
       write (message,'(2a,es12.5,a)') ch10,&
       ' --- Kohn-Sham susceptibility matrices for frequency = ',freq_mem(jfreq),'i'
       call wrtout(std_out,message,'COLL')
       do isp1 = 1,nspden
         do isp2 = 1,nspden
           write (message,'(5x,a,2i2)') 'Susceptibility matrix for spins = ',isp1,isp2
           call wrtout(std_out,message,'COLL')
           write (message,'(10x,a,12x,a,13x,a,9x,a)') "g","g'","real","imag"
           call wrtout(std_out,message,'COLL')
           do ipw1 = 1,min(48,npwdiel)
             do ipw2 = ipw1,min(48,npwdiel)
               write (message,'(a,1x,3i4,1x,3i4,2x,es12.5,1x,es12.5)') '+S',&
&               kg_diel(1:3,ipw1),kg_diel(1:3,ipw2),&
&               susmat_dyn(1,ipw1,isp1,ipw2,isp2,jfreq),susmat_dyn(2,ipw1,isp1,ipw2,isp2,jfreq)
               call wrtout(std_out,message,'COLL')
               write (message,'(a,1x,12x,1x,12x,2x,es12.5,1x,es12.5)') '+S',&
&               susmat_dyn(1,ipw2,isp1,ipw1,isp2,jfreq),susmat_dyn(2,ipw2,isp1,ipw1,isp2,jfreq)
               call wrtout(std_out,message,'COLL')
             end do
           end do
         end do
       end do
     end do
   end if

!  Enforce strict zero 1st column and row (particle number conservation).

   do jfreq = 1,nfreq_mem
     do isp = 1,nspden
       do ipw = 1,npwdiel
         susmat_dyn(:,1,isp,ipw,isp,jfreq) = 0._dp
         susmat_dyn(:,ipw,isp,1,isp,jfreq) = 0._dp
       end do
     end do
   end do

!  Calculate the non-interacting polarizabilities and the ACFD contributions to the
!  exchange energy.

   write (message,'(4a)') ch10,&
&   ' xcacfd: Calculating the non-interacting Kohn-Sham polarizabilities',ch10,&
&   '         and the ACFD contributions to the exchange energy...'
   call wrtout(std_out,message,'COLL')

   do jfreq = 1,nfreq_mem

!    Calculate the non-interacting polarizabilities from the zero G limit of
!    the diagonal of the susceptibility matrix ; the limit is performed in
!    different ways.

     susd_tmp(:) = 0._dp
     do isp = 1,nspden
       do ipw = 1,npwdiel
         susd_tmp(ipw) = susd_tmp(ipw)+susmat_dyn(1,ipw,isp,ipw,isp,jfreq)
       end do
     end do

     call get_susd_null(ig_tiny,igsq_tiny,gsq,npwdiel,npw_tiny,sus_gabs,sus_gavg,sus_gdir,susd_tmp)

     write (message,fmth2) 'frequency',(jj,jj = 0,npw_tiny-1)
     call wrtout(std_out,message,'COLL')

!    write (message,fmtd) &
!    &  ' alpha0_sq',freq_mem(jfreq),(-ucvol*sus_gabs(jj),jj = 1,npw_tiny)
!    call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alpha0_xx',freq_mem(jfreq),(-ucvol*sus_gdir(jj,1),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alpha0_yy',freq_mem(jfreq),(-ucvol*sus_gdir(jj,2),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alpha0_zz',freq_mem(jfreq),(-ucvol*sus_gdir(jj,3),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     write (message,fmtd) &
&     ' alpha0_av',freq_mem(jfreq),(-ucvol*sus_gavg(jj),jj = 1,npw_tiny)
     call wrtout(std_out,message,'COLL')

     if (ifreq*jfreq == 1) alpha_non0 = -ucvol*sus_gavg(npw_tiny)
     alpha_non = -ucvol*sus_gavg(npw_tiny)

!    Calculate the ACFD contributions to the exchange energy.
!    Only works in the spin-restricted case.

     susd_tmp(1:npwdiel) = susd_aux_dyn(1,1:npwdiel,1,jfreq)

!    Using the bare Coulomb interaction:

     call geteexc_uc(energy,energy_raw,gsq,ig_tiny,npwdiel,npw_tiny,susd_tmp)
     energy(:) = -energy(:)/two_pi !factor -1/(2*pi) from fluctuation-dissipation theorem.
     energy_raw = -energy_raw/two_pi

     write (message,fmtd) &
&     ' dEx_bare',freq_mem(jfreq),energy(1:npw_tiny)
     call wrtout(std_out,message,'COLL')

!    Using the cut-off Coulomb interaction:

     call geteexc_cc(energy,energy_raw,gsq,npwdiel,npw_tiny,rcut_coulomb,susd_tmp)
     energy(1) = -energy(1)/two_pi
     energy_raw = -energy_raw/two_pi

     write (message,fmtd) &
&     ' dEx_cut',freq_mem(jfreq),energy(1)
     call wrtout(std_out,message,'COLL')

     ex_cc = ex_cc+wght_mem(jfreq)*energy(1)

!    write (message,fmtd) &
!    &  ' dn_test',freq_mem(jfreq),susd_aux_dyn(2,1:1,1,jfreq)*ucvol/pi
!    call wrtout(std_out,message,'COLL')

     n_test = n_test-wght_mem(jfreq)*susd_aux_dyn(2,1,1,jfreq)*ucvol/pi

   end do

!  DEBUG
!  write(std_out,*)' xcacfd : will call acfd_...'
!  call flush(6)
!  ENDDEBUG

!  Compute the interacting susceptibility matrices along the adiabatic connection path.

   do jfreq = 1,nfreq_mem

     if (intexact > 0) then

       call acfd_intexact(decs_cc,freq_mem(jfreq),gsq,ikhxc,mband,dtset%nband,&
&       nkpt,npwdiel,nspden,nsppol,occ,dtset%occopt,rcut_coulomb,&
&       susmat_dyn(1,1,1,1,1,jfreq))

       ecs_cc = ecs_cc+wght_mem(jfreq)*decs_cc

     else

!      DEBUG
!      write(std_out,*)'xcacfd :  will call acfd_dyson '
!      call flush(6)
!      ENDDEBUG

!      HERE


       call acfd_dyson(dtset,freq_mem(jfreq),gsq,idyson,ig_tiny,igsq_tiny,ikhxc,kg_diel,khxc,&
&       ldgapp,mpi_enreg,ndyson,nfft,ngfft,npw_tiny,npwdiel,nspden,2,rcut_coulomb,&
&       rhor,rhocut,rprimd,susd_data,dtset%suskxcrs,susmat_dyn(1,1,1,1,1,jfreq),ucvol)

!      DEBUG
!      write(std_out,*)' xcacfd :  called acfd_dyson '
!      call flush(6)
!      ENDDEBUG

!      susd_data(:,1) = real part of the coupling-constant integrated diagonal
!      of $\chi_{\lambda}-chi_0$ (summed over spin-density components if appropriate).
!      susd_data(:,2) = dto. in the Lein, Dobson and Gross first-order approximation,
!      (summed over spin-density components if appropriate).
!      susd_data(:,3) = real part of the diagonal of $\chi_{\lambda=1}$
!      (summed over spin-density components if appropriate).

       call prtsusd(gsq,ig_tiny,index_g,npw_tiny,npwdiel,0,susd_data)

!      Calculate the non-interacting polarizabilities and the ACFD contributions
!      to the correlation energy.

       write (message,'(4a)') ch10,&
&       ' xcacfd: Calculating the interacting polarizabilities',ch10,&
&       '         and the ACFD contributions to the correlation energy...'
       call wrtout(std_out,message,'COLL')

!      Calculate the interacting polarizabilities

       susd_tmp(1:npwdiel) = susd_data(1:npwdiel,3)
       call get_susd_null(ig_tiny,igsq_tiny,gsq,npwdiel,npw_tiny,sus_gabs,&
&       sus_gavg,sus_gdir,susd_tmp)

       write (message,fmth2) 'frequency',(jj,jj = 0,npw_tiny-1)
       call wrtout(std_out,message,'COLL')

!      write (message,fmtd) &
!      &   ' alpha_sq',freq_mem(jfreq),(-ucvol*sus_gabs(jj),jj = 1,npw_tiny)
!      call wrtout(std_out,message,'COLL')

       write (message,fmtd) &
&       ' alpha1_xx',freq_mem(jfreq),(-ucvol*sus_gdir(jj,1),jj = 1,npw_tiny)
       call wrtout(std_out,message,'COLL')

       write (message,fmtd) &
&       ' alpha1_yy',freq_mem(jfreq),(-ucvol*sus_gdir(jj,2),jj = 1,npw_tiny)
       call wrtout(std_out,message,'COLL')

       write (message,fmtd) &
&       ' alpha1_zz',freq_mem(jfreq),(-ucvol*sus_gdir(jj,3),jj = 1,npw_tiny)
       call wrtout(std_out,message,'COLL')

       write (message,fmtd) &
&       ' alpha1_av',freq_mem(jfreq),(-ucvol*sus_gavg(jj),jj = 1,npw_tiny)
       call wrtout(std_out,message,'COLL')

       if (ifreq*jfreq == 1) alpha_int0 = -ucvol*sus_gavg(npw_tiny)
       alpha_int = -ucvol*sus_gavg(npw_tiny)

!      Calculate the ACFD contributions to the correlation energy.

       susd_tmp(:) = susd_data(:,1)

!      Using the bare Coulomb interaction:

       call geteexc_uc(energy,energy_raw,gsq,ig_tiny,npwdiel,npw_tiny,susd_tmp)
       energy(:) = -energy(:)/two_pi !factor -1/(2*pi) from fluctuation-dissipation theorem.
       energy_raw = -energy_raw/two_pi

       write (message,fmtd) &
&       ' dEc_bare',freq_mem(jfreq),energy(1:npw_tiny)
       call wrtout(std_out,message,'COLL')

       ecs_uc(:) = ecs_uc(:)+wght_mem(jfreq)*energy(:)

!      Using the cut-off Coulomb interaction:

       call geteexc_cc(energy,energy_raw,gsq,npwdiel,npw_tiny,rcut_coulomb,susd_tmp)
       energy(1) = -energy(1)/two_pi
       energy_raw = -energy_raw/two_pi

       write (message,fmtd) &
&       ' dEc_cut',freq_mem(jfreq),energy(1)
       call wrtout(std_out,message,'COLL')

       ecs_cc = ecs_cc+wght_mem(jfreq)*energy(1)

       if (ldgapp > 0) then

!        Calculate the contribution to the correlation energy
!        in the Lein, Dobson and Gross approximation.

         susd_tmp(:) = susd_data(:,2)

!        Using the bare Coulomb interaction:

         call geteexc_uc(energy,energy_raw,gsq,ig_tiny,npwdiel,npw_tiny,susd_tmp)
         energy(:) = -energy(:)/two_pi !factor -1/(2*pi) from fluctuation-dissipation theorem.
         energy_raw = -energy_raw/two_pi

         write (message,fmtd) &
&         ' dEcLDG_bare',freq_mem(jfreq),energy(1:npw_tiny)
         call wrtout(std_out,message,'COLL')

         ecgl_uc(:) = ecgl_uc(:)+wght_mem(jfreq)*energy(:)

!        Using the cut-off Coulomb interaction:

         call geteexc_cc(energy,energy_raw,gsq,npwdiel,npw_tiny,rcut_coulomb,susd_tmp)
         energy(1) = -energy(1)/two_pi
         energy_raw = -energy_raw/two_pi

         write (message,fmtd) &
&         ' dEcLDG_cut',freq_mem(jfreq),energy(1)
         call wrtout(std_out,message,'COLL')

         ecgl_cc = ecgl_cc+wght_mem(jfreq)*energy(1)

       end if

     end if

   end do

!  End big loop over frequencies
 end do

 close(tmp_unit)

!Print frequency integrated quantities and run parameters.

 write (message,'(2a)') ch10,&
& ' Final results of the ACFD calculation (atomic units):'
 call wrtout(std_out,message,'COLL')

 write (message,fmtoi) &
& ' ikhxc',ikhxc,'['//trim(ikhxc_name(ikhxc))//' kernel used to calculate the ACFD energy]'
 call wrtout(std_out,message,'COLL')
 write (message,fmtoi) &
& ' nfreqsus',nfreqsus,'[Number of points used for frequency integration]'
 call wrtout(std_out,message,'COLL')
 write (message,fmtoi) &
& ' npwdiel',npwdiel,'[Number of planewaves in the susceptibility matrices]'
 call wrtout(std_out,message,'COLL')
 write (message,fmtoi) &
& ' nbandsus',nbandsus,'[Number of bands used for the susceptibility matrices]'
 call wrtout(std_out,message,'COLL')
 write (message,fmtod) &
& ' eigen_max',eigen(nbandsus),'[Highest eigenvalue used for the susceptibility matrices]'
 call wrtout(std_out,message,'COLL')
 if (intexact <= 0) then
   write (message,fmtoi) &
&   ' idyson',idyson,'[Solve the Dyson equation as a '//trim(idyson_name(idyson))//']'
   call wrtout(std_out,message,'COLL')
   write (message,fmtoi) &
&   ' ndyson',ndyson,'[Number of points used for coupling-constant integration, excl. 0, 1]'
   call wrtout(std_out,message,'COLL')
 else
   write (message,fmtoi) &
&   ' intexact',intexact,'[Exact integration over the coupling constant]'
   call wrtout(std_out,message,'COLL')
 end if
 write (message,fmtod) &
& ' rcut_coulomb',rcut_coulomb,'[Cut-off radius for the Coulomb interaction]'
 call wrtout(std_out,message,'COLL')
 write (message,fmtod) &
& ' n_test',n_test,'[Should be equal to the number of electrons]'
 call wrtout(std_out,message,'COLL')

 do jj = npw_tiny,1,-1
   write (message,fmtog) &
&   ' ExDM_bare',jj-1,ex_dm_uc(jj),'[Dens. mat. with G = 0 substracted, ',jj-1,'th order extrap. to G = 0]'
   call wrtout(std_out,message,'COLL')
 end do
 write (message,fmtod) &
& ' ExDM_bare_nozerog',ex_dm_uc_nozerog,'[Dens. mat. without G = 0 component]'
 call wrtout(std_out,message,'COLL')
 write (message,fmtod) &
& ' ExDM_cut',ex_dm_cc,'[Dens. mat.]'
 call wrtout(std_out,message,'COLL')
 write (message,fmtod) &
& ' Ex_cut',ex_cc,'[ACFD exchange energy]'
 call wrtout(std_out,message,'COLL')

 if (intexact <= 0) then
   do jj = npw_tiny,1,-1
     write (message,fmtog) &
&     ' Ec_bare',jj-1,ecs_uc(jj),'[',jj-1,'th order extrap. to G = 0]'
     call wrtout(std_out,message,'COLL')
   end do
 end if
 write (message,fmtod) &
& ' Ec_cut',ecs_cc
 call wrtout(std_out,message,'COLL')

 if (ldgapp > 0) then
   do jj = npw_tiny,1,-1
     write (message,fmtog) &
&     ' EcLDG_bare',jj-1,ecgl_uc(jj),'[Lein, Dobson and Gross app., ',jj-1,'th order extrap. to G = 0]'
     call wrtout(std_out,message,'COLL')
   end do
   write (message,fmtod) &
&   ' EcLDG_cut',ecgl_cc,'[Lein, Dobson and Gross app.]'
   call wrtout(std_out,message,'COLL')
 end if

!Give a message with some results.

 write (message,'(2a,3(2a,i15))') ch10,&
& ' xcacfd: check on dynamical susceptibility for imaginary frequencies -----------',ch10,&
& '  Number of planewaves in the susceptibility matrices          : ',npwdiel,ch10,&
& '  Number of bands used for the susceptibility matrices         : ',nbandsus,ch10,&
& '  Number of imaginary frequencies                              : ',nfreqsus
 call wrtout(ab_out,message,'COLL')
!call wrtout(std_out,message,'COLL')
 write (message,'(a,es15.8,2(2a,es15.8))') &
& '  Min. frequency (Ha)                                          : ',freq(1),ch10,&
& '  Max. frequency (Ha)                                          : ',freq(nfreqsus),ch10,&
& '  Coulomb cut-off radius (bohr)                                : ',rcut_coulomb
 call wrtout(ab_out,message,'COLL')
!call wrtout(std_out,message,'COLL')
 if (intexact <= 0) then
   write (message,'(a,i15,4a)') &
&   '  Number of points used for coupling-constant int. (excl. 0, 1): ',ndyson,ch10,&
&   '  Dyson equation solved as a ',trim(idyson_name(idyson)),'.'
   call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,message,'COLL')
 end if
 write (message,'(a,es15.8,2a,es15.8)') &
& '  Average Kohn-Sham polarizability at min. frequency (1/bohr^3): ',alpha_non0,ch10,&
& '  Average Kohn-Sham polarizability at max. frequency (1/bohr^3): ',alpha_non
 call wrtout(ab_out,message,'COLL')
!call wrtout(std_out,message,'COLL')
 if (intexact <= 0) then
   write (message,'(a,a9,a,es15.8,2a,a9,a,es15.8)') &
&   '  Average ',trim(ikhxc_name(ikhxc)),' polarizability at min. frequency (1/bohr^3): ',alpha_int0,ch10,&
&   '  Average ',trim(ikhxc_name(ikhxc)),' polarizability at max. frequency (1/bohr^3): ',alpha_int
   call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,message,'COLL')
 end if
 write (message,'(a,es15.8,2a,es15.8)') &
& '  Kohn-Sham exchange energy from density matrix (Ha)           : ',ex_dm_cc,ch10,&
& '. ACFD Kohn-Sham exchange energy (Ha)                          : ',ex_cc
 call wrtout(ab_out,message,'COLL')
!call wrtout(std_out,message,'COLL')
 write (message,'(3a,t64,a,es15.8)') &
& '. ACFD-',trim(ikhxc_name(ikhxc)),' correlation energy (Ha)',': ',ecs_cc
 call wrtout(ab_out,message,'COLL')
!call wrtout(std_out,message,'COLL')
 if (ldgapp > 0) then
   write (message,'(a,es15.8)') &
&   '. Lein, Dobson and Gross first-order approximation (Ha)        : ',ecgl_cc
   call wrtout(ab_out,message,'COLL')
!  call wrtout(std_out,message,'COLL')
 end if

!Free memory.

 ABI_DEALLOCATE(gsq)
 ABI_DEALLOCATE(index_g)
 ABI_DEALLOCATE(freq_mem)
 ABI_DEALLOCATE(wght_mem)
 ABI_DEALLOCATE(susd_tmp)
 ABI_DEALLOCATE(susd_data)
 ABI_DEALLOCATE(susd_aux_dyn)
 ABI_DEALLOCATE(susmat_dyn)
 if (associated(khxc))  then
   ABI_DEALLOCATE(khxc)
 end if

!Exit instructions.

 write (message,'(2a)') ch10,' xcacfd: done & exit'
 call wrtout(std_out,message,'COLL')

end subroutine xcacfd
!!***
