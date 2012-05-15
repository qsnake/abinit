!{\src2tex{textfont=tt}}
!!****f* ABINIT/print_psps
!! NAME
!! print_psps
!!
!! FUNCTION
!!  Method to print the content of a pseudopotential_type derived type
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  psps=<type pseudopotential_type>
!!  unit(optional)=unit number for output
!!  prtvol(optional)=verbosity level
!!  mode_paral(optional): either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only writing 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Should add information coming from pspheads
!!
!! PARENTS
!!      bethe_salpeter,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_psps(psps,unit,prtvol,mode_paral)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_psps'
 use interfaces_14_hidewrite
 use interfaces_57_iovars
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
!scalars
 integer :: ierr,ips,ipsp_alch,ityp_alch,itypat,unt,verb
 character(len=4) :: mode
 character(len=500) :: msg
!arrays
 integer :: cond_values(3)
 character(len=9) :: cond_string(3)

! *************************************************************************

 ! Some initialisations
 verb=0      ; if (PRESENT(prtvol))     verb=prtvol
 unt=std_out ; if (PRESENT(unit))     unt=unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral
 ierr=0 ; cond_string(1:3)=' ' ; cond_values(1:3)=(/0,0,0/)
 !
 ! /** General info including spin-orbit **/
 write(msg,'(2a)')ch10,&
& ' ==== Info on pseudopotentials ==== '
 call wrtout(unt,msg,mode)
 SELECT CASE (psps%usepaw) 
 CASE (0)
  write(msg,'(a)')'  Norm-conserving pseudopotentials '
  call wrtout(unt,msg,mode)
  write(msg,'(a,i4)')'  Max number of Kleinman-Bylander energies ',psps%dimekb
  call wrtout(unt,msg,mode)
  !do itypat=1,psps%ntypat 
  ! write(msg,'(a,i4,a,f9.4)')' Type ',itypat,' K-B energies ',(psps%ekb(ikbe,itypat),ikbe=1,psps%dimekb)
  !end do
 CASE (1)
  write(msg,'(a)')'  PAW calculation'
  call wrtout(unt,msg,mode)
  write(std_out,*)'  Max number of D_ij coefficients ',psps%dimekb
 CASE DEFAULT 
  call chkint_eq(0,0,cond_string,cond_values,ierr,'usepaw',psps%usepaw,2,(/0,1/),unt)
 END SELECT

  !integer :: dimekb
  ! Dimension of Ekb
  ! ->Norm conserving : Max. number of Kleinman-Bylander energies
  !                     for each atom type
  !                     dimekb=lnmax (lnmax: see this file)
  ! ->PAW : Max. number of Dij coefficients connecting projectors
  !                     for each atom type
  !                     dimekb=lmnmax*(lmnmax+1)/2 (lmnmax: see this file)

 !real(dp), pointer :: ekb(:,:)
  ! ekb(dimekb,ntypat*(1-usepaw))
  !  ->NORM-CONSERVING PSPS ONLY:
  !    (Real) Kleinman-Bylander energies (hartree)
  !           for number of basis functions (l,n) (lnmax)
  !           and number of atom types (ntypat)
  ! NOTE (MT) : ekb (norm-conserving) is now diagonal (one dimension
  !             lnmax); it would be easy to give it a second
  !             (symmetric) dimension by putting
  !             dimekb=lnmax*(lnmax+1)/2
  !             in the place of dimekb=lmnmax.
 SELECT CASE (psps%positron)
 CASE (0) 
  !write(std_out,*)' Standard Electron Calculation '
 CASE (1,2)
  write(msg,'(a,i3)')'  Positron Calculation with positron .. ',psps%positron 
  call wrtout(unt,msg,mode)
 CASE DEFAULT
  call chkint_eq(0,0,cond_string,cond_values,ierr,'positron',psps%positron,3,(/0,1,2/),unt)
 END SELECT

 write(msg,'(a,i4,2a,i4)')&
& '  Number of pseudopotentials .. ',psps%npsp,ch10,&
& '  Number of types of atoms   .. ',psps%ntypat 
 call wrtout(unt,msg,mode)

 SELECT CASE (psps%mpspso) 
 CASE (1) 
  write(msg,'(a)')'  Calculation without spin-orbit '
  call wrtout(unt,msg,mode)
 CASE (2)
  write(msg,'(3a,i3)')&
&  '  Calculation with spin-orbit coupling ',ch10,&
&  '  Max number of channels (spin-orbit included) ',psps%mpssoang
  call wrtout(unt,msg,mode)
  do itypat=1,psps%ntypat 
   if (psps%pspso(itypat)==2) then 
    write(msg,'(a,i4,a)')'  - Atom type ',itypat,' has spin-orbit characteristics'
    call wrtout(unt,msg,mode)
   end if 
  end do
 CASE DEFAULT
  call chkint_eq(0,0,cond_string,cond_values,ierr,'mpspso',psps%mpspso,2,(/1,2/),unt)
 END SELECT
 !
 ! /** Info on nonlocal part **/
 !
 SELECT CASE (psps%useylm)
 CASE (0)
  write(msg,'(a)')'  Nonlocal part applied using Legendre polynomials '
 CASE (1)
  write(msg,'(a)')'  Nonlocal part applied using real spherical harmonics '
 CASE DEFAULT
  call chkint_eq(0,0,cond_string,cond_values,ierr,'psps%useylm',psps%useylm,2,(/0,1/),unt)
 END SELECT
 call wrtout(unt,msg,mode)

 !FIXME this does not work, it seems it is always 0 , except for HGH
 !write(msg,'(a,i3)')' Max number of non-local projectors over l and type ',psps%mproj 
 !if (psps%mproj==0) then 
 ! write(msg,'(a)')TRIM(msg)//' (All local) '
 !end if
 !call wrtout(unt,msg,mode)
 write(msg,'(a,i3,2a,i3,2a,i3)')&
& '  Highest angular momentum +1 ....... ',psps%mpsang,ch10,&
& '  Max number of (l,n)   components .. ',psps%lnmax, ch10,&
& '  Max number of (l,m,n) components .. ',psps%lmnmax
 call wrtout(unt,msg,mode)
  !integer :: lnmax
  !  Max. number of (l,n) components over all type of psps
  !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
  !  so, it is equal to the max of lnprojso, see pspheader_type
 !integer :: lmnmax
  !  If useylm=0, max number of (l,m,n) comp. over all type of psps (lnproj)
  !  If useylm=1, max number of (l,n)   comp. over all type of psps (lmnproj)
  !  If mpspso is 2, lmnmax takes into account the spin-orbit projectors,
  !  so, it is equal to the max of lmnprojso or lnprojso, see pspheader_type

!$integer, pointer :: indlmn(:,:,:)
! indlmn(6,lmnmax,ntypat)
! For each type of psp,
! array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0)
!                                or i=lmn (if useylm=1)

 !FIXME for paw n1xccc==1
 !
 ! /** Non-linear Core correction **/
 if (psps%n1xccc/=0) then 
  write(msg,'(3a,2(a,i4,a),2a)')ch10,&
&  ' *** Pseudo-Core Charge Info *** ',ch10,&
&  '  Number of radial points for pseudo-core charge .. ',psps%n1xccc,ch10,&
&  '  XC core-correction treatment (optnlxccc) ........ ',psps%optnlxccc,ch10,&
&  '  Radius for pseudo-core charge for each type ..... ',ch10
  call wrtout(unt,msg,mode)
  do itypat=1,psps%ntypat 
   write(msg,'(a,i4,a,f7.4)')'  - Atom type ',itypat,' has pseudo-core radius .. ',psps%xcccrc(itypat)
   call wrtout(unt,msg,mode)
  end do
 end if
 !
 ! /** Alchemical mixing **/
 if (psps%mtypalch/=0) then 
  write(msg,'(3a,3(a,i4,a))')ch10,&
&  ' *** Calculation with alchemical mixing *** ',ch10,&
   '  Number of pure pseudoatoms .... ',psps%ntyppure,ch10,&
   '  Number of pseudos for mixing .. ',psps%npspalch,ch10,&
   '  Alchemical pseudoatoms ........ ',psps%ntypalch,ch10
  call wrtout(unt,msg,mode)
  do ipsp_alch=1,psps%npspalch 
   do ityp_alch=1,psps%ntypalch 
    write(std_out,*)' mixalch ',psps%mixalch(ipsp_alch,ityp_alch)
   end do
  end do
  do ityp_alch=1,psps%ntypalch 
   write(msg,'(a,i4,a,i4)')' For alchemical atom no. ',ityp_alch,' algalch is .. ',psps%algalch(ityp_alch)
   call wrtout(unt,msg,mode)
  end do
 end if
 !integer :: mtypalch
  ! Maximum number of alchemical pseudo atoms. If non-zero,
  ! the mechanism to generate mixing of pseudopotentials is activated
 !integer :: ntypat
  ! Number of types of atoms (might be alchemy wrt pseudopotentials)
 !integer :: ntyppure
  ! Number of types of pure pseudoatoms
 !integer :: ntypalch
  ! Number of types of alchemical pseudoatoms
 !integer :: npspalch
  ! Number of types of pseudopotentials use for alchemical purposes
 !integer, pointer :: algalch(:)   ! algalch(ntypalch)
  ! For each type of pseudo atom, the algorithm to mix the pseudopotentials
 !real(dp), pointer :: mixalch(:,:)
  ! mixalch(npspalch,ntypalch)
  ! Mixing coefficients to generate alchemical pseudo atoms


 !
 ! /** Info in Q-grid for spline of form factors **/
 !
 write(msg,'(3a,a,i6,a,a,i6)')ch10,&
& ' *** Info on the Q-grid used for form factors in spline form *** ',ch10,&
& '  Number of q-points for radial functions ffspl .. ',psps%mqgrid_ff,ch10,&
& '  Number of q-points for vlspl ................... ',psps%mqgrid_vl 
 call wrtout(unt,msg,mode)
 if (psps%vlspl_recipSpace) then 
  write(msg,'(a)')'  vlspl is computed in Reciprocal Space '
 else 
  write(msg,'(a)')'  vlsp is computed in Real Space '
 end if
 call wrtout(unt,msg,mode)
 !TODO additional stuff tbat might be printed

 !real(dp), pointer :: ffspl(:,:,:,:)
  ! ffspl(mqgrid_ff,2,lnmax,ntypat)
  ! Gives, on the radial grid, the different non-local projectors,
  ! in both the norm-conserving case, and the PAW case
 !real(dp), pointer :: qgrid_ff(:)
  ! qgrid_ff(mqgrid_ff)
  ! The coordinates of all the points of the radial grid for the nl form factors
 !real(dp), pointer :: qgrid_vl(:)
   ! qgrid_vl(mqgrid_vl)
   ! The coordinates of all the points of the radial grid for the local part of psp
 !real(dp), pointer :: vlspl(:,:,:)
  ! vlspl(mqgrid_vl,2,ntypat)
  ! Gives, on the radial grid, the local part of each type of psp.
 ! real(dp), pointer :: dvlspl(:,:,:)
  ! dvlspl(mqgrid_vl,2,ntypat)
  ! Gives, on the radial grid, the first derivative of the local
  ! part of each type of psp (computed when the flag 'vlspl_recipSpace'
  ! is true).
 !real(dp), pointer :: xccc1d(:,:,:)
  ! xccc1d(n1xccc*(1-usepaw),6,ntypat)
  ! Norm-conserving psps only
  ! The component xccc1d(n1xccc,1,ntypat) is the pseudo-core charge
  ! for each type of atom, on the radial grid. The components
  ! xccc1d(n1xccc,ideriv,ntypat) give the ideriv-th derivative of the
  ! pseudo-core charge with respect to the radial distance.

  !write(msg,'(2a)')ch10,' Z_ion pseudo Z_at ' 
  !call wrtout(unt,msg,mode)
  !do ips=1,psps%npsp
  ! write(std_out,*)psps%zionpsp(ips),psps%znuclpsp(ips)
  !end do

   !real(dp), pointer :: zionpsp(:)
   ! zionpsp(npsp)
   ! For each pseudopotential, the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)
   !real(dp), pointer :: ziontypat(:)
   ! ziontypat(ntypat)
   !  For each type of atom (might be alchemy wrt psps), the ionic pseudo-charge
   ! (giving raise to a long-range coulomb potential)
   !real(dp), pointer :: znuclpsp(:)
   ! znuclpsp(npsp)
   ! The atomic number of each pseudopotential
   !real(dp), pointer :: znucltypat(:)
   ! znucltypat(ntypat)
   ! The atomic number of each type of atom (might be alchemy wrt psps)

  do itypat=1,psps%ntypat 
   write(msg,'(a,i3,a,i3)')' XC functional for type ',itypat,' is ',psps%pspxc(itypat)
   call wrtout(unt,msg,mode)
   !write(std_out,*)psps%ziontypat(itypat),psps%znucltypat(itypat)
  end do
  !integer, pointer :: pspxc(:)
   ! pspxc(ntypat)
   ! For each type of psp, the XC functional that was used to generate it,
   ! as given by the psp file

  if (verb>=3) then 
   do ips=1,psps%npsp
    write(std_out,*)' Pseudo number   ',ips,' read from ',TRIM(psps%filpsp(ips))
    write(std_out,*)' Format or Code  ',psps%pspcod(ips)
    write(std_out,*)' Generation Date ',psps%pspdat(ips)
    write(std_out,*)' Content of first line ',TRIM(psps%title(ips))
   end do
  end if

  !character(len=fnlen), pointer :: filpsp(:)
   ! filpsp(ntypat)
   ! The filename of the pseudopotential
  !character(len=fnlen), pointer :: title(:)
   ! title(ntypat)
   ! The content of first line read from the psp file
!  integer, pointer :: pspdat(:)
   ! pspdat(ntypat)
   ! For each type of psp, the date of psp generation, as given by the psp file
  !integer, pointer :: pspcod(:)
   ! pspcod(npsp)
   ! For each type of psp, the format -or code- of psp generation,
   !  as given by the psp file


! Types for pseudo-potentials that are based on parameters. Currently, only
! GTH are supported (see pseudopotential_gth_type). To add one, one should
! create an initialisation method and a destruction method in 02psp (see
! psp2params.F90). These methods are called in driver().
!TODO this is still missing
!  type(pseudopotential_gth_type) :: gth_params

 ! If there was a problem, then stop.
 if (ierr/=0) call leave_new('COLL')

end subroutine print_psps
!!***

!!****f* ABINIT/plot_psps
!! NAME
!! plot_psps
!!
!! FUNCTION
!!  Writes on external files some of the arrays defined in the pseudopotential_type.
!!
!! INPUTS
!!
!! OUTPUT
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


subroutine plot_psps(psps,root_filename)

 use m_profiling

 use defs_basis
 use defs_datatypes

 use m_io_tools, only : get_unit

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'plot_psps'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in),optional :: root_filename
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
!scalars
 integer :: ider,iln,iq,ir,ityp,unt
 character(len=100) :: fmt
 character(len=fnlen) :: fname,root

! *************************************************************************

 root='PSPS'
 if (present(root_filename)) root=root_filename
 unt=get_unit()

 !TODO most of the pointer are not nullified, 
 !this part will print a lot of quantities that are not used

 if (associated(psps%vlspl)) then 
  fname=trim(root)//'_VLSPL'
  open(unit=unt,file=fname,status='new',form='formatted')

  write(unt,*)' Local part of each type of atom '  
  write(unt,*)' q-mesh and, for each type, v_loc(q)? and second derivative '
  write(fmt,*)'(',1+2*psps%ntypat,'(es17.9,1x))'
  do iq=1,psps%mqgrid_vl
   write(unt,fmt)psps%qgrid_vl(iq),((psps%vlspl(iq,ider,ityp),ider=1,2),ityp=1,psps%ntypat)
  end do
  close(unt)
 end if

 if (associated(psps%ffspl)) then 
  !TODO write error handler for open
  !here I need the pseudo_header to avoid writing columns made of zero 
  fname=trim(root)//'_FFSPL'
  open(unit=unt,file=fname,status='new',form='formatted')

  write(unt,*)' Form factors for each type of atom '  
  write(unt,*)' q-mesh and, for each type and each (l,n) channel, ffnl(q) and second derivative '
  write(fmt,*)'(',1+2*psps%lnmax*psps%ntypat,'(es17.9,1x))'
  do iq=1,psps%mqgrid_ff
   write(unt,fmt)psps%qgrid_ff(iq),&
&   (((psps%ffspl(iq,ider,iln,ityp),ider=1,2),iln=1,psps%lnmax),ityp=1,psps%ntypat)
  end do
  close(unt)
 end if

 if (associated(psps%xccc1d)) then 
  !TODO write error handler for open

  fname=trim(root)//'_PSCC'
  open(unit=unt,file=fname,status='new',form='formatted')
  write(unt,*)' Pseudo-core charge for each type of atom, on the radial grid. '
  write(unt,*)' radial-mesh and, for each type rho_pscore(r) '
  write(fmt,*)'(',psps%ntypat,'(es17.9,1x))'
  ! TODO Grid is missing
  do ir=1,psps%n1xccc
   !write(unt,fmt)((psps%xccc1d(ir,1,itypat)),itypat=1,psps%ntypat)
  end do
  close(unt)
 end if

end subroutine plot_psps 
!!***
