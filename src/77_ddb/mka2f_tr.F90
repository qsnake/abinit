!{\src2tex{textfont=tt}}
!!****f* ABINIT/mka2f_tr
!!
!! NAME
!! mka2f_tr
!!
!! FUNCTION
!!  calculates the FS averaged Transport alpha^2F_tr alpha^2F_trout alpha^2F_trin functions
!!  calculates and outputs the associated electrical and thermal conductivities
!!  for the first task : copied from mka2F
!!
!!
!! COPYRIGHT
!! Copyright (C) 2004-2012 ABINIT group (JPC, MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYINGS=
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  alter_int_gam = flag to use alternative integration scheme
!!  elph_ds
!!    elph_ds%gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!    elph_ds%nbranch = number of phonon branches = 3*natom
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_fine%nkpt = number of kpts included in the FS integration
!!    elph_ds%k_fine%wtk = integration weights on the FS
!!    delph_ds%n0 = DOS at the Fermi level calculated from the k_fine integration weights
!!    elph_ds%k_fine%kpt = coordinates of all FS kpoints
!!  gprim = reciprocal lattice vectors (maybe dimensioned...)
!!  gprimd = reciprocal lattice vectors (dimensionful)
!!  mustar = coulomb pseudopotential parameter
!!       eventually for 2 spin channels
!!  ucvol = Unit cell volume 
!!  natom = number of atoms
!!  nrpt = number of real-space points for FT interpolation
!!  ntemper = number of temperature points to calculate, from tempermin to tempermin+ntemper*temperinc
!!  tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!  temperinc = interval for temperature grid on which resistivity etc are calculated (in K)
!!  phon_ds = datastructure with interatomic force constants to interpolate
!!     phonons
!!  rpt = coordinates of real-space points for FT interpolation
!!  wghatm = weights for real-space points for FT interpolation
!!
!! OUTPUT
!!  elph_ds
!!    
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgam,gam_mult_displ,inpphon,leave_new,simpson_int,wrtout,zgemm
!!
!! NOTES
!!   copied from ftiaf9.f
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mka2f_tr(alter_int_gam,elph_ds,gprim,gprimd,ucvol,natom,nrpt,&
&   ntemper,tempermin,temperinc,phon_ds,rpt,wghatm,elph_tr_ds)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_io_tools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mka2f_tr'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_77_ddb, except_this_one => mka2f_tr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: alter_int_gam
 integer,intent(in) :: natom,nrpt
 integer,intent(in) :: ntemper
 real(dp),intent(in) :: ucvol
 real(dp),intent(in) :: tempermin,temperinc
 type(elph_tr_type) :: elph_tr_ds
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 real(dp),intent(in) :: gprim(3,3), gprimd(3,3), rpt(3,nrpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)

!Local variables -------------------------
!x =w/(2kbT)
!scalars
 integer :: iFSqpt,ibranch,iomega,iost,isppol,jbranch,nerr
 integer :: unit_a2f_tr, unit_a2f_trout, unit_a2f_trin
 integer :: idir, iatom, k1, kdir
 integer :: unit_lor,unit_rho,unit_tau,unit_therm
 integer :: itemp
 integer :: itrtensor, icomp, jcomp
 real(dp) :: Temp,chgu,chtu,diagerr,firh,firhT,gaussfactor,domega
 real(dp) :: firh_tau,firhT_tau ! added by BX to get Tau
 real(dp) :: a2fprefactor_in, temp_in
 real(dp) :: a2fprefactor_out, temp_out
 real(dp) :: gaussprefactor,gaussval,lambda_tr,lor0,lorentz,maxerr,maxx,omega
 real(dp) :: rho,tau,tolexp,wtherm,xtr,xx
 real(dp) :: lambda_tr_trace
 real(dp) :: omega_min, omega_max
 real(dp) :: qnorm2
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),parameter :: c0(2)=(/0.d0,0.d0/),c1(2)=(/1.d0,0.d0/)
 real(dp) :: eigval_in(elph_ds%nbranch)
 real(dp) :: eigval_out(elph_ds%nbranch)
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: gam_now_in (2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: gam_now_out(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: pheigval(elph_ds%nbranch)
 real(dp) :: tmpa2f_in (elph_ds%na2f)
 real(dp) :: tmpa2f_out(elph_ds%na2f)
 real(dp) :: tmpgam1(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmpgam2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),allocatable :: phfrq(:,:)
 real(dp),allocatable :: displ(:,:,:,:)
 real(dp),allocatable :: pheigvec(:,:)
 real(dp),allocatable :: integrho(:),integtau(:),matrx(:,:),tointegrho(:),tointegtau(:),zhpev1(:,:)
 real(dp),allocatable :: zhpev2(:)
 real(dp),allocatable :: rho_T(:),tau_T(:)

! *********************************************************************
!calculate a2f_tr for frequencies between 0 and omega_max


 write(std_out,*) 'mka2f_tr : enter '
!
!MG: the step should be calculated locally using nomega and the extrema of the spectrum.
!One should not rely on previous calls for the setup of elph_ds%domega
!I will remove elph_ds%domega since mka2f.F90 will become a method of gamma_t
 domega =elph_ds%domega

 ABI_ALLOCATE(elph_tr_ds%a2f_1d_tr,(elph_ds%na2f,9,elph_ds%nsppol))
 ABI_ALLOCATE(elph_tr_ds%a2f_1d_trin,(elph_ds%na2f,9,elph_ds%nsppol))
 ABI_ALLOCATE(elph_tr_ds%a2f_1d_trout,(elph_ds%na2f,9,elph_ds%nsppol))

!! defaults for number of temperature steps and max T (all in Kelvin...)
!ntemper=1000
!tempermin=zero
!temperinc=one
 ABI_ALLOCATE(rho_T,(ntemper))
 ABI_ALLOCATE(tau_T,(ntemper))


!tolerance on gaussian being = 0
 tolexp = 1.d-100
 maxx = sqrt(-log(tolexp))
 lor0=(pi*kb_HaK)**2/3.

!maximum value of frequency (a grid has to be chosen for the representation of alpha^2 F)
!WARNING! supposes this value has been set in mkelph_linwid.

 gaussprefactor = sqrt(piinv) / elph_ds%a2fsmear
 gaussfactor = one / elph_ds%a2fsmear


 ABI_ALLOCATE(matrx,(2,(3*natom*(3*natom+1))/2))
 ABI_ALLOCATE(zhpev1,(2,2*3*natom-1))
 ABI_ALLOCATE(zhpev2,(3*3*natom-2))

!ENDMG

 elph_tr_ds%a2f_1d_tr = zero
 elph_tr_ds%a2f_1d_trin = zero
 elph_tr_ds%a2f_1d_trout = zero

 maxerr=0.
 nerr=0

 ABI_ALLOCATE(phfrq,(elph_ds%nbranch, elph_ds%k_fine%nkpt))
 ABI_ALLOCATE(displ,(2, elph_ds%nbranch, elph_ds%nbranch, elph_ds%k_fine%nkpt))
 ABI_ALLOCATE(pheigvec,(2*elph_ds%nbranch*elph_ds%nbranch, elph_ds%k_fine%nkpt))
 do iFSqpt=1,elph_ds%k_fine%nkpt
   call inpphon(displ(:,:,:,iFSqpt),pheigval,pheigvec(:,iFSqpt),phfrq(:,iFSqpt),phon_ds,elph_ds%k_fine%kpt(:,iFSqpt))
 end do
 omega_min = minval(phfrq)
 omega_max = maxval(phfrq)

 do isppol=1,elph_ds%nsppol

!  loop over qpoint in full kpt grid (presumably dense)
   do iFSqpt=1,elph_ds%k_fine%nkpt
     qnorm2 = sum(elph_ds%k_fine%kpt(:,iFSqpt)**2)
!    if (flag_to_exclude_soft_modes = .false.) qnorm2 = zero
     do itrtensor=1,9

       if (alter_int_gam == 0) then
!        Do FT from real-space gamma grid to 1 qpt.
         call ftgam(wghatm,gam_now_in, elph_tr_ds%gamma_rpt_trin(:,itrtensor,:,isppol,:),gprim,natom,1,nrpt,0,rpt,&
&         elph_ds%k_fine%kpt(:,iFSqpt))
         call ftgam(wghatm,gam_now_out,elph_tr_ds%gamma_rpt_trout(:,itrtensor,:,isppol,:),gprim,natom,1,nrpt,0,rpt,&
&         elph_ds%k_fine%kpt(:,iFSqpt))


       else if (alter_int_gam == 1) then
!        in the alter_int_gam case the gamma_qpt are already interpolated on the k_fine grid:
         gam_now_in(:,:)  = elph_tr_ds%gamma_qpt_trin(:,itrtensor,:,isppol,iFSqpt)
         gam_now_out(:,:) = elph_tr_ds%gamma_qpt_trout(:,itrtensor,:,isppol,iFSqpt)
       end if



!      Diagonalize gamma matrix at this qpoint (complex matrix).

!      if ep_scalprod==0 we have to dot in the displacement vectors here
       if (elph_ds%ep_scalprod==0) then

         displ_red(:,:,:) = zero
         do jbranch=1,elph_ds%nbranch
           do iatom=1,natom
             do idir=1,3
               ibranch=idir+3*(iatom-1)
               do kdir=1,3
                 k1 = kdir+3*(iatom-1)
                 displ_red(1,ibranch,jbranch) = displ_red(1,ibranch,jbranch) + &
&                 gprimd(kdir,idir)*displ(1,k1,jbranch,iFSqpt) 
                 displ_red(2,ibranch,jbranch) = displ_red(2,ibranch,jbranch) + &
&                 gprimd(kdir,idir)*displ(2,k1,jbranch,iFSqpt)
               end do
             end do
           end do
         end do

         tmpgam2 = reshape (gam_now_in, (/2,elph_ds%nbranch,elph_ds%nbranch/))
         call gam_mult_displ(elph_ds%nbranch, displ_red, tmpgam2, tmpgam1)
         do jbranch=1,elph_ds%nbranch
           eigval_in(jbranch)   = tmpgam1(1, jbranch, jbranch)
         end do

         tmpgam2 = reshape (gam_now_out, (/2,elph_ds%nbranch,elph_ds%nbranch/))
         call gam_mult_displ(elph_ds%nbranch, displ_red, tmpgam2, tmpgam1)
         do jbranch=1,elph_ds%nbranch
           eigval_out(jbranch)   = tmpgam1(1, jbranch, jbranch)
         end do
         
       else if (elph_ds%ep_scalprod == 1) then

!        
!        NOTE: in these calls gam_now and pheigvec do not have the right rank, but blas usually does not care
!        
         call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now_in, 3*natom,&
&         pheigvec(:,iFSqpt), 3*natom, c0, tmpgam1, 3*natom)
         call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec(:,iFSqpt), 3*natom,&
&         tmpgam1, 3*natom, c0, tmpgam2, 3*natom)
         diagerr = zero
         do ibranch=1,elph_ds%nbranch
           eigval_in(ibranch) = tmpgam2(1,ibranch,ibranch)
           do jbranch=1,ibranch-1
             diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
           end do
           do jbranch=ibranch+1,elph_ds%nbranch
             diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
           end do
         end do
         if (diagerr > tol12) then
           nerr=nerr+1
           maxerr=max(diagerr, maxerr)
         end if
         
         call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now_out, 3*natom,&
&         pheigvec(:,iFSqpt), 3*natom, c0, tmpgam1, 3*natom)
         call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec(:,iFSqpt), 3*natom,&
&         tmpgam1, 3*natom, c0, tmpgam2, 3*natom)
         diagerr = zero
         do ibranch=1,elph_ds%nbranch
           eigval_out(ibranch) = tmpgam2(1,ibranch,ibranch)
           do jbranch=1,ibranch-1
             diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
           end do
           do jbranch=ibranch+1,elph_ds%nbranch
             diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))
           end do
         end do
         if (diagerr > tol12) then
           nerr=nerr+1
           maxerr=max(diagerr, maxerr)
         end if
       end if
!      end ep_scalprod if

!      Add all contributions from the phonon modes at this qpoint to
!      a2f and the phonon dos.
       do ibranch=1,elph_ds%nbranch
!        if (abs(phfrq(ibranch,iFSqpt)) < tol10) then
         if ( abs(phfrq(ibranch,iFSqpt)) < tol7 .or. &
&         (phfrq(ibranch,iFSqpt) < tol4 .and. qnorm2 > 0.03 )) then !
!          note: this should depend on the velocity of sound, to accept acoustic
!          modes!
           a2fprefactor_in = zero
           a2fprefactor_out= zero
         else
           a2fprefactor_in  = eigval_in (ibranch)/(two_pi*abs(phfrq(ibranch,iFSqpt))*elph_ds%n0(isppol))
           a2fprefactor_out = eigval_out(ibranch)/(two_pi*abs(phfrq(ibranch,iFSqpt))*elph_ds%n0(isppol))
         end if

!        omega = omega_min
!        tmpa2f_in (:) = zero
!        tmpa2f_out(:) = zero
!        do iomega=1,elph_ds%na2f
!        xx = (omega-phfrq(ibranch,iFSqpt))*gaussfactor
!        gaussval = gaussprefactor*exp(-xx*xx)
!        tmpa2f_in (iomega) = tmpa2f_in (iomega) + gaussval*a2fprefactor_in
!        tmpa2f_out(iomega) = tmpa2f_out(iomega) + gaussval*a2fprefactor_out
!        omega = omega+domega
!        end do

         omega = omega_min
         tmpa2f_in (:) = zero
         tmpa2f_out(:) = zero
         do iomega=1,elph_ds%na2f
           xx = (omega-phfrq(ibranch,iFSqpt))*gaussfactor
           gaussval = gaussprefactor*exp(-xx*xx)

           temp_in = gaussval*a2fprefactor_in
           temp_out = gaussval*a2fprefactor_out

           if (dabs(temp_in) < 1.0d-50) temp_in = zero
           if (dabs(temp_out) < 1.0d-50) temp_out = zero
           tmpa2f_in (iomega) = tmpa2f_in (iomega) + temp_in
           tmpa2f_out(iomega) = tmpa2f_out(iomega) + temp_out
           omega = omega+domega
         end do

         elph_tr_ds%a2f_1d_trin (:,itrtensor,isppol) = elph_tr_ds%a2f_1d_trin (:,itrtensor,isppol) + tmpa2f_in(:)
         elph_tr_ds%a2f_1d_trout(:,itrtensor,isppol) = elph_tr_ds%a2f_1d_trout(:,itrtensor,isppol) + tmpa2f_out(:)

       end do ! end ibranch do
     end do ! end itrtensor do
   end do ! end iFSqpt do
 end do ! end isppol

!second 1 / elph_ds%k_fine%nkpt factor for the integration weights
 elph_tr_ds%a2f_1d_trin  = elph_tr_ds%a2f_1d_trin  / elph_ds%k_fine%nkpt
 elph_tr_ds%a2f_1d_trout = elph_tr_ds%a2f_1d_trout / elph_ds%k_fine%nkpt

 if (elph_ds%ep_scalprod == 1) then
   write(std_out,*) 'mka2f_tr: errors in diagonalization of gamma_tr with phon eigenvectors: ', nerr,maxerr
 end if

 elph_tr_ds%a2f_1d_tr = elph_tr_ds%a2f_1d_trout - elph_tr_ds%a2f_1d_trin

 ABI_DEALLOCATE(matrx)
 ABI_DEALLOCATE(zhpev1)
 ABI_DEALLOCATE(zhpev2)

!output the elph_tr_ds%a2f_1d_tr
 unit_a2f_tr = get_unit()
 fname = trim(elph_ds%elph_base_name) // '_A2F_TR'
 open (unit=unit_a2f_tr,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 unit_a2f_trin = get_unit()
 fname = trim(elph_ds%elph_base_name) // '_A2F_TRIN'
 open (unit=unit_a2f_trin,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
 unit_a2f_trout = get_unit()
 fname = trim(elph_ds%elph_base_name) // '_A2F_TROUT'
 open (unit=unit_a2f_trout,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 write (unit_a2f_tr,'(a)')       '#'
 write (unit_a2f_tr,'(a)')       '# ABINIT package : a2f_tr file'
 write (unit_a2f_tr,'(a)')       '#'
 write (unit_a2f_tr,'(a)')       '# a2f_tr function integrated over the FS. omega in a.u.'
 write (unit_a2f_tr,'(a,I10)')   '#     number of kpoints integrated over : ', elph_ds%k_fine%nkpt
 write (unit_a2f_tr,'(a,I10)')   '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f_tr,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', omega_min, &
& ' Ha and omega_max = ', omega_max, ' Ha'
 write (unit_a2f_tr,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_a2f_tr,'(a)')       '#'

 write (unit_a2f_trin,'(a)')       '#'
 write (unit_a2f_trin,'(a)')       '# ABINIT package : a2f_trin file'
 write (unit_a2f_trin,'(a)')       '#'
 write (unit_a2f_trin,'(a)')       '# a2f_trin function integrated over the FS. omega in a.u.'
 write (unit_a2f_trin,'(a,I10)')   '#     number of kpoints integrated over : ', elph_ds%k_fine%nkpt
 write (unit_a2f_trin,'(a,I10)')   '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f_trin,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', omega_min, &
& ' Ha and omega_max = ', omega_max, ' Ha'
 write (unit_a2f_trin,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_a2f_trin,'(a)')       '#'

 write (unit_a2f_trout,'(a)')       '#'
 write (unit_a2f_trout,'(a)')       '# ABINIT package : a2f_trout file'
 write (unit_a2f_trout,'(a)')       '#'
 write (unit_a2f_trout,'(a)')       '# a2f_trout function integrated over the FS. omega in a.u.'
 write (unit_a2f_trout,'(a,I10)')   '#     number of kpoints integrated over : ', elph_ds%k_fine%nkpt
 write (unit_a2f_trout,'(a,I10)')   '#     number of energy points : ',elph_ds%na2f
 write (unit_a2f_trout,'(a,E16.6,a,E16.6,a)') '#       between omega_min = ', omega_min, &
& ' Ha and omega_max = ', omega_max, ' Ha'
 write (unit_a2f_trout,'(a,E16.6)') '#   and the smearing width for gaussians is ', elph_ds%a2fsmear
 write (unit_a2f_trout,'(a)')       '#'

!done with header
 do isppol=1,elph_ds%nsppol
   write (unit_a2f_tr,'(a,E16.6)') '# The DOS at Fermi level is ', elph_ds%n0(isppol)
   write (unit_a2f_trin,'(a,E16.6)') '# The DOS at Fermi level is ', elph_ds%n0(isppol)
   write (unit_a2f_trout,'(a,E16.6)') '# The DOS at Fermi level is ', elph_ds%n0(isppol)
!  omega = zero
   omega = omega_min
   do iomega=1,elph_ds%na2f
     write (unit_a2f_tr,   '(10D16.6)') omega, elph_tr_ds%a2f_1d_tr   (iomega,:,isppol)
     write (unit_a2f_trin, '(10D16.6)') omega, elph_tr_ds%a2f_1d_trin (iomega,:,isppol)
     write (unit_a2f_trout,'(10D16.6)') omega, elph_tr_ds%a2f_1d_trout(iomega,:,isppol)
     omega=omega+domega
   end do
   write (unit_a2f_tr,*)
   write (unit_a2f_trin,*)
   write (unit_a2f_trout,*)
 end do !isppol

 close (unit=unit_a2f_tr)
 close (unit=unit_a2f_trin)
 close (unit=unit_a2f_trout)

!calculation of transport properties
 ABI_ALLOCATE(integrho,(elph_ds%na2f))
 ABI_ALLOCATE(tointegrho,(elph_ds%na2f))
 ABI_ALLOCATE(integtau,(elph_ds%na2f))
 ABI_ALLOCATE(tointegtau,(elph_ds%na2f))

 unit_rho = get_unit()
 fname = trim(elph_ds%elph_base_name) // '_RHO'
 open (unit=unit_rho,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
!print header to resistivity file
 write (unit_rho,*) '# Resistivity as a function of temperature.'
 write (unit_rho,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_rho,*) '#  '
 write (unit_rho,*) '#  Columns are: '
 write (unit_rho,*) '#  temperature[K]   rho[au]   rho [SI]        rho/temp [au]'
 write (unit_rho,*) '#  '

 unit_tau = get_unit()
 fname = trim(elph_ds%elph_base_name) // '_TAU'
 open (unit=unit_tau,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if
!print header to relaxation time file
 write (unit_tau,*) '# Relaxation time as a function of temperature.'
 write (unit_tau,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_tau,*) '#  '
 write (unit_tau,*) '#  Columns are: '
 write (unit_tau,*) '#  temperature[K]   tau[au]   tau [SI]     '
 write (unit_tau,*) '#  '

 unit_therm = get_unit()
 fname = trim(elph_ds%elph_base_name) // '_WTH'
 open (unit=unit_therm,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!print header to thermal conductivity file
 write (unit_therm,'(a)') '# Thermal conductivity/resistivity as a function of temperature.'
 write (unit_therm,'(a)') '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_therm,'(a)') '#  '
 write (unit_therm,'(a)') '#  Columns are: '
 write (unit_therm,'(a)') '#  temperature[K]   thermal rho[au]   thermal cond [au]   thermal rho [SI]   thermal cond [SI]'
 write (unit_therm,'(a)') '#  '

 unit_lor = get_unit()
 fname = trim(elph_ds%elph_base_name) // '_LOR'
 open (unit=unit_lor,file=fname,status='unknown',iostat=iost)
 if (iost /= 0) then
   write (message,'(3a)')' mka2f_tr : ERROR- opening file ',trim(fname),' as new'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!print header to lorentz file
 write (unit_lor,*) '# Lorentz number as a function of temperature.'
 write (unit_lor,*) '#  the formalism is isotropic, so non-cubic crystals may be wrong'
 write (unit_lor,*) '#  '
 write (unit_lor,*) '#  Columns are: '
 write (unit_lor,*) '#  temperature[K]   Lorentz number[au]   Lorentz quantum = (pi*kb_HaK)**2/3'
 write (unit_lor,*) '#  '

 do isppol=1,elph_ds%nsppol
   lambda_tr_trace = zero
   do itrtensor=1,9
     omega = omega_min
     tointegrho = zero
     do iomega=1,elph_ds%na2f
       if(omega<=0) then
         omega=omega+domega
         cycle
       end if
       tointegrho(iomega)=two*elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol)/omega
       omega=omega+domega
     end do

     integrho = zero
     call simpson_int(elph_ds%na2f,domega,tointegrho,integrho)
     lambda_tr=integrho(elph_ds%na2f)
     write (message, '(a,2i3,a,es16.6)' )&
&     ' mka2f_tr: TRANSPORT lambda for isppol itrtensor', isppol, itrtensor, ' =  ', lambda_tr
     call wrtout(std_out,message,'COLL')
     if (itrtensor == 1 .or. itrtensor == 5 .or. itrtensor == 9) lambda_tr_trace = lambda_tr_trace + lambda_tr
   end do !end itrtensor do

   lambda_tr_trace = lambda_tr_trace / three
   write (message, '(a,i3,a,es16.6)' )&
&   ' mka2f_tr: 1/3 trace of TRANSPORT lambda for isppol ', isppol, ' =  ', lambda_tr_trace
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end do !end isppol do

!constant to change units of rho from au to SI
 chgu=2.173969*1.0d-7
 chtu=2.4188843265*1.0d-17

 do isppol=1,elph_ds%nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itrtensor=(icomp-1)*3+jcomp

!      prefactor for resistivity integral
!      firh=6.d0*pi*ucvol*kb_HaK/(elph_ds%n0(isppol)*elph_tr_ds%FSelecveloc_sq(isppol))
!      FIXME: check factor of 2 which is different from Savrasov paper. 6 below for thermal conductivity is correct.
       firh=2.d0*pi*ucvol*kb_HaK/elph_ds%n0(isppol)/&
&       sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))

!      Add by BX to get Tau_elph
       firh_tau = 4.0d0*pi*kb_HaK
!      End Adding

       write(unit_rho,*) '# Rho for isppol, itrten = ', isppol, itrtensor
       write(unit_tau,*) '# Tau for isppol, itrten = ', isppol, itrtensor

       tointegrho(:)=0.
       do itemp=1,ntemper  ! runs over termperature in K
         Temp=tempermin+temperinc*dble(itemp)
         firhT=firh*Temp
         firhT_tau=firh_tau*Temp
         omega = omega_min
         do iomega=1,elph_ds%na2f
           if(omega<=0) then
             omega=omega+domega
             cycle
           end if
           xtr=omega/(2*kb_HaK*Temp)
           tointegrho(iomega)=firhT*omega*elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol)  &
&           /(((2*Temp*kb_HaK)**2)*((exp(xtr)-exp(-xtr))/2)**2)
!          Add by BX to get Tau
           tointegtau(iomega)=firhT_tau*omega*elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol)  &
&           /(((2*Temp*kb_HaK)**2)*((exp(xtr)-exp(-xtr))/2)**2)

           omega=omega+domega
         end do

         call simpson_int(elph_ds%na2f,domega,tointegrho,integrho)
         call simpson_int(elph_ds%na2f,domega,tointegtau,integtau)
         rho=integrho(elph_ds%na2f)
         tau=1.0d0/integtau(elph_ds%na2f)
         write(unit_rho,'(4D20.10)')temp,rho,rho*chgu,rho/temp
         write(unit_tau,'(3D20.10)')temp,tau,tau*chtu
         rho_T(itemp)=rho
         tau_T(itemp)=tau
       end do ! temperature
       write(unit_rho,*)
       write(unit_tau,*)

     end do ! jcomp
   end do ! icomp
 end do ! isppol

!-----------------------------


 do isppol=1,elph_ds%nsppol
   do icomp=1, 3
     do jcomp=1, 3
       itrtensor=(icomp-1)*3+jcomp
!      prefactor for integral of thermal conductivity
!      firh=(18.*ucvol)/(pi*kb_HaK*elph_ds%n0(isppol)*elph_tr_ds%FSelecveloc_sq(isppol))
       firh=(6.d0*ucvol)/(pi*kb_HaK*elph_ds%n0(isppol))/ &
&       sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))


       write(unit_therm,*) '# Thermal resistivity for isppol, itrten= ', isppol
       write(unit_lor,*) '# Lorentz coefficient for isppol, itrten= ', isppol

       tointegrho(:)=0.
       do itemp=1,ntemper

         Temp=tempermin + temperinc*dble(itemp)
         omega = omega_min
         do iomega=1,elph_ds%na2f
           if(omega<=0) then
             omega=omega+domega
             cycle
           end if
           xtr=omega/(2*kb_HaK*Temp)
           tointegrho(iomega)=xtr**2/omega*&
&           ( elph_tr_ds%a2f_1d_tr(iomega,itrtensor,isppol)+&
&           4*xtr**2*elph_tr_ds%a2f_1d_trout(iomega,itrtensor,isppol)/pi**2+   &
&           2*xtr**2*elph_tr_ds%a2f_1d_trin(iomega,itrtensor,isppol)/pi**2)  &
&           /(((exp(xtr)-exp(-xtr))/2)**2)

           omega=omega+domega
         end do

         call simpson_int(elph_ds%na2f,domega,tointegrho,integrho)
         wtherm=integrho(elph_ds%na2f)*firh

         write(unit_therm,'(5D20.10)')temp,wtherm,1./wtherm,wtherm/3.4057d9,1./(wtherm) *3.4057d9

         lorentz=rho_T(itemp)/(wtherm*temp)
         write(unit_lor,*)temp,lorentz,lor0

       end do
       write(unit_therm,*)
       write(unit_lor,*)
     end do ! jcomp
   end do ! icomp
 end do !end isppol do


 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(displ)
 ABI_DEALLOCATE(pheigvec)

 close (unit=unit_lor)
 close (unit=unit_rho)
 close (unit=unit_tau)
 close (unit=unit_therm)

 ABI_DEALLOCATE(integrho)
 ABI_DEALLOCATE(tointegrho)
 ABI_DEALLOCATE(elph_tr_ds%a2f_1d_tr)
 ABI_DEALLOCATE(elph_tr_ds%a2f_1d_trin)
 ABI_DEALLOCATE(elph_tr_ds%a2f_1d_trout)
 
 ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_trin)
 ABI_DEALLOCATE(elph_tr_ds%gamma_qpt_trout)
 if (alter_int_gam == 0) then
   ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_trin)
   ABI_DEALLOCATE(elph_tr_ds%gamma_rpt_trout)
 end if
 write(std_out,*) ' mka2f_tr : end '


end subroutine mka2f_tr
!!***
