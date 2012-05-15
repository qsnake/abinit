!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_check
!! NAME
!! hdr_check
!!
!! FUNCTION
!! This subroutine compare the header structured variable (hdr)
!! from input data (mostly dtset and psps) with the one (hdr0) of
!! an input data file (e.g. wf, density, potential).
!! Various values are checked for agreement or near agreement in the
!! case of floating point numbers.  The program will exit or produce
!! warning messages when unexpected values are found.
!! A record of the comparison of the headers is written to stdout.
!!
!! Decisions have been taken about whether a restart is allowed.
!! In the self-consistent case, a restart will always be allowed, but
!! one has to distinguish between a direct restart and a restart with
!! translation of wavefunction (for which the functionalities of newsp will
!! be needed). In the non-self-consistent case, the conditions below
!! must be fulfilled to allow a restart.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (DCA, XG, GMR, FB).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  fform=integer specification of data type (expected)
!!  fform0=integer specification of data type (from disk file)
!!  mode_paral : COLL or PERS, for all leave_new and wrtout calls
!!  hdr <type(hdr_type)>=the header structured variable from dtset and psps
!!  hdr0<type(hdr_type)>=the header structured variable from the disk file
!!
!! OUTPUT
!!  restart=1 if direct restart, =2 if translation is needed, =0 if no
!!              restart is possible.
!!  restartpaw= deals with the additional informations in the PAW method
!!              =1 if direct restart, =0 if no restart from spherical data is possible.
!!              also 0 if no PAW
!!
!! NOTES
!! In the current version of the user interface restarts are allowed from
!! wavefunction files for self-consistent runs and from densities for
!! non-self-consistent runs. The precise conditions under which we will
!! allow a restart in this release are as follows.
!!
!!           self-consistent case : direct restarts
!!           ======================================
!!
!! A direct restart will be allowed provided the following quantities in
!! old and new calculations are the same:
!!
!!   (A) the primitive vectors                             (tprim)
!!   (B) the plane-wave cutoff                             (tecut)
!!   (C) nkpt, kpt(3,nkpt), wtk(nkpt)                      (tkpt)
!!   (D) istwfk(nkpt), the format of wavefunctions         (twfk)
!!   (E) nspinor, the scalar or spinor wf characteristics  (tspinor)
!! For PAW calculations:
!!   (F) the use of PAW method                             (tpaw)
!!   (G) the number of lmn elements for the paw basis      (tlmn)
!!   (H) the energy cutoff for the double (fine) grid      (tdg)
!! For WVL calculations:
!!   (I) the number of wavelets differs                    (twvl)
!!   (J) the space-grid size differs                       (tgrid)
!!
!!            non-self-consistent restarts
!!            ============================
!!
!! A restart will be allowed provided the following quantities in
!! old and new calculation are the same
!!
!!   (A) the primitive vectors                            (tprim)
!!   (B) the number of atoms of each type                 (tatty)
!!   (C) xred(3,natom)                                    (txred)
!!   (D) pseudopotentials (not just pseudocharges)        (tpseu)
!!   (E) the plane-wave cutoff                            (tecut)
!!   (F) ngfft(1:3)                                       (tng)
!! For PAW calculations:
!!   (G) the use of PAW method                            (tpaw)
!!   (H) the number of lmn elements for the paw basis     (tlmn)
!!   (I) the energy cutoff for the double (fine) grid     (tdg)
!!
!! PARENTS
!!      inwffil,ioarr,m_gwannier,m_io_screening,setup_bse,setup_screening
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hdr_check(fform,fform0,hdr,hdr0,mode_paral,restart,restartpaw)

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hdr_check'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_59_io_mpi, except_this_one => hdr_check
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fform,fform0
 integer,intent(out) :: restart,restartpaw
 character(len=4),intent(in) :: mode_paral
 type(hdr_type),intent(in) :: hdr,hdr0

!Local variables-------------------------------
 character(len=1), parameter :: number(0:10)=(/'0','1','2','3','4','5','6','7','8','9',' '/)
 character(len=13), parameter :: filtypes(5)=(/'wf_planewave ','density      ','potential    ','screening    ', 'wf_wavelet   '/)
 character(len=24), save :: bndfmt='(2x, i4,t41,   a,2x, i4)'
 character(len=28), save :: occfmt='(2x, f4.1,t41,   a,2x, f4.1)'
 character(len=28), save :: wtkfmt='(2x, f7.3,t41,   a,2x, f7.3)'
 character(len=28), save :: zatfmt='(2x, f6.2,t41,   a,2x, f6.2)'
!scalars
 integer,parameter :: mwarning=5,nkpt_max=5
 integer :: bantot,bantot_eff,ii,ipsp,isppol,istart,istop,isym,itest,iwarning
 integer :: jj,mu,natom,nelm,nkpt,npsp,nsppol,nsym,ntypat,tatty,tband,tdg
 integer :: tecut,tgrid,tkpt,tlmn,tng,tpaw,tprim,tpsch,tpseu,tspinor,tsym,twfk
 integer :: twvl,txred
 real(dp) :: rms
 logical :: tfform2,tfform52
 character(len=13) :: filtyp,filtyp0
 character(len=26) :: typfmt
 character(len=500) :: message

! *************************************************************************

 DBG_ENTER("COLL")

!We will adopt convention that if things agree between restart
!and current calculation then the tflag is 0. Begin by assuming
!that there is complete agreement between the files

 tatty = 0; tband = 0; tdg = 0 ; tecut = 0; tkpt = 0;
 tlmn = 0; tng = 0; tpaw = 0; tprim = 0; tpsch = 0; tpseu = 0;
 tspinor=0; tsym = 0; twfk = 0 ; txred = 0 ; twvl = 0 ; tgrid = 0

!Write out a header
 write(message,&
& '(a1,80a,2a1,10x,a,3a1,8x,a,25x,a,a1,8x,19a,25x,12a,a1)' )&
& ch10,('=',ii=1,80),ch10,ch10,&
& '- hdr_check: checking restart file header for consistency -',&
& (ch10,ii=1,3),'current calculation',&
& 'restart file',ch10,('-',ii=1,19),('-',ii=1,12),ch10
 call wrtout(std_out,message,mode_paral)

!Check validity of fform, and find filetype
 ii=1+fform /50
 if(ii>5)then
   if(fform==1002)then
     ii=4
   else
     write(message, '(4a,i8)' ) ch10,&
&     ' hdr_check: BUG -',ch10,&
&     '  Incorrect file format, fform=',fform
     call wrtout(std_out,message,mode_paral)
     call leave_new(mode_paral)
   end if
 end if
 filtyp=filtypes(ii)

!Check validity of fform0, and find filetype
 ii=1+fform0/50
 if(ii>5)then
   if(fform0==1002)then
     ii=4
   else
     write(message, '(4a,i8,a)' ) ch10,&
&     ' hdr_check: ERROR -',ch10,&
&     '  Incorrect file format, fform0=',fform0,&
&     '  Action : it seems that the file you try to read is not an appropriate file.'
     call wrtout(std_out,message,mode_paral)
     call leave_new(mode_paral)
   end if
 end if
 filtyp0=filtypes(ii)

 write(message,'(a,a,3x,3a)') &
& '  calculation expects a ',filtyp,'|',&
& '  input file contains a ',filtyp0
 call wrtout(std_out,message,mode_paral)

 write(message,'(a,a,11x,a,a,a)')&
& '. ABINIT  code version ',hdr%codvsn,'|',&
& '  ABINIT  code version ',hdr0%codvsn
 call wrtout(std_out,message,mode_paral)

!Check fform from input, not from header file
 if ( (fform+1)/2 /= (fform0+1)/2  ) then
   write(message, '(a,a,a,a,i10,a,i10,a)' ) ch10,&
&   ' hdr_check: BUG -',ch10,&
&   '  input fform=',fform,' differs from disk file fform=',fform0,'.'
   call wrtout(std_out,message,mode_paral)
   call leave_new(mode_paral)
 end if

 write(message, '(a,i8,a,i4,a,i4,2x,a,a,i8,a,i4,a,i4)' ) &
& '. date ',hdr %date,' bantot ',hdr %bantot,' natom ',hdr %natom,'|',&
& '  date ',hdr0%date,' bantot ',hdr0%bantot,' natom ',hdr0%natom
 call wrtout(std_out,message,mode_paral)

 write(message, '(a,i4,a,i3,3(a,i4),2x,a,a,i4,a,i3,3(a,i4))' )&
& '  nkpt',hdr %nkpt,' nsym',hdr %nsym,&
& ' ngfft',hdr %ngfft(1),',',hdr %ngfft(2),',',hdr %ngfft(3),&
& '|','  nkpt',hdr0%nkpt,' nsym',hdr0%nsym,&
& ' ngfft',hdr0%ngfft(1),',',hdr0%ngfft(2),',',hdr0%ngfft(3)
 call wrtout(std_out,message,mode_paral)

 if (hdr%usewvl == 0) then
!  Note that the header actually contains ecut_eff=ecut*dilatmx**2
   write(message,'(a,i3,a,f12.7,8x,a,a,i3,a,f12.7)')&
&   '  ntypat',hdr %ntypat,' ecut_eff',hdr %ecut_eff,'|',&
&   '  ntypat',hdr0%ntypat,' ecut_eff',hdr0%ecut_eff
   call wrtout(std_out,message,mode_paral)
 else
   write(message,'(a,i3,a,f12.7,8x,a,a,i3,a,f12.7)')&
&   '  ntypat',hdr %ntypat,' hgrid   ', 2. * hdr %rprimd(1,1) / (hdr %ngfft(1) - 31),'|',&
&   '  ntypat',hdr0%ntypat,' hgrid   ', 2. * hdr0%rprimd(1,1) / (hdr0%ngfft(1) - 31)
   call wrtout(std_out,message,mode_paral)
!  Check hgrid and rprimd values.
   if (hdr0%rprimd(1,2) /= zero .or. hdr0%rprimd(1,3) /= zero .or. &
&   hdr0%rprimd(2,1) /= zero .or. hdr0%rprimd(2,3) /= zero .or. &
&   hdr0%rprimd(3,1) /= zero .or. hdr0%rprimd(3,2) /= zero) then
     write(message, '(a,a,a,a)' ) ch10,&
&     ' hdr_check: ERROR -',ch10,&
&     '  disk file rprimd is not parallelepipedic.'
     call wrtout(std_out,message,mode_paral)
     call leave_new(mode_paral)
   end if
   if (abs(hdr0%rprimd(1,1) / hdr0%ngfft(1) - hdr %rprimd(1,1) / hdr %ngfft(1)) > tol8) then
     write(message, '(a,a,a,a,F7.4,a,F7.4)' ) ch10,&
&     ' hdr_check: WARNING -',ch10,&
&     '  input wvl_hgrid=', 2. * hdr%rprimd(1,1) / hdr%ngfft(1), &
&     ' not equal disk file wvl_hgrid=', 2. * hdr0%rprimd(1,1) / hdr0%ngfft(1)
     call wrtout(std_out,message,mode_paral)
     tgrid = 1
   end if
 end if

 write(message, '(a,i3,29x,a,a,i3)' )&
& '  usepaw',hdr %usepaw,'|',&
& '  usepaw',hdr0%usepaw
 call wrtout(std_out,message,mode_paral)

 write(message, '(a,i3,29x,a,a,i3)' )&
& '  usewvl',hdr %usewvl,'|',&
& '  usewvl',hdr0%usewvl
 call wrtout(std_out,message,mode_paral)

 write(message,'(a,31x,a,a,3(a1,2x,3f12.7,2x,a,2x,3f12.7))')&
& '  rprimd:','|','  rprimd:',ch10,&
& hdr%rprimd(:,1),'|',hdr0%rprimd(:,1),ch10,&
& hdr%rprimd(:,2),'|',hdr0%rprimd(:,2),ch10,&
& hdr%rprimd(:,3),'|',hdr0%rprimd(:,3)
 call wrtout(std_out,message,mode_paral)

 if (hdr%bantot/=hdr0%bantot) then
   tband=1
 end if

 if (hdr%intxc/=hdr0%intxc) then
   write(message, '(a,a,a,a,i5,a,i5)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  input intxc=',hdr%intxc,' not equal disk file intxc=',hdr0%intxc
   call wrtout(std_out,message,mode_paral)
 end if

 if (hdr%ixc/=hdr0%ixc) then
   write(message, '(a,a,a,a,i5,a,i5)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  input ixc=',hdr%ixc,' not equal disk file ixc=',hdr0%ixc
 end if

 if (hdr%natom/=hdr0%natom) then
   write(message, '(a,a,a,a,i8,a,i8)') ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  input natom=',hdr%natom,' not equal disk file natom=',hdr0%natom
   call wrtout(std_out,message,mode_paral)
   tatty=1
 end if

 if ( ANY(hdr%ngfft/=hdr0%ngfft) ) then
!  For sensible rho(r) or V(r) data, fft grid must be identical
!  MG TODO one should perform an FFT interpolation when the two ngfft differ!
   if (fform==52.or.fform==102) then
     write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&     ' hdr_check: ERROR -',ch10,&
&     '  fft grids must be the same for restart from a ',filtyp,' file.',ch10,&
&     '  Action : change your fft grid or your restart file.'
     call wrtout(std_out,message,mode_paral)
     call leave_new(mode_paral)
   end if
   tng=1
 end if

 if (hdr%nkpt/=hdr0%nkpt) then
   if (fform==2) then
     write(message,'(a,a,a,a,i8,a,i8)' ) ch10,&
&     ' hdr_check: WARNING -',ch10,&
&     '  input nkpt=',hdr%nkpt,' not equal disk file nkpt=',hdr0%nkpt
     call wrtout(std_out,message,mode_paral)
   end if
   tkpt=1
   twfk=1
 end if

 if (hdr%nspinor/=hdr0%nspinor) then
   if (fform==2) then
     write(message,'(a,a,a,a,i8,a,i8)' ) ch10,&
&     ' hdr_check: WARNING -',ch10,&
&     '  input nspinor=',hdr%nspinor,' not equal disk file nspinor=',hdr0%nspinor
     call wrtout(std_out,message,mode_paral)
   end if
   tspinor=1
 end if

!No check is present for nspden

 if (hdr%nsppol/=hdr0%nsppol) then
   write(message, '(a,a,a,a,i6,a,i6)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  input nsppol=',hdr%nsppol,&
&   ' not equal disk file nsppol=',hdr0%nsppol
   call wrtout(std_out,message,mode_paral)
 end if

 if (hdr%nsym/=hdr0%nsym) then
   write(message, '(a,a,a,a,i6,a,i6)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  input nsym=',hdr%nsym,' not equal disk file nsym=',hdr0%nsym
   call wrtout(std_out,message,mode_paral)
   tsym=1
 end if

 if (hdr%ntypat/=hdr0%ntypat) then
   write(message, '(a,a,a,a,i5,a,i5)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  input ntypat=',hdr%ntypat,' not equal disk file ntypat=',hdr0%ntypat
   call wrtout(std_out,message,mode_paral)
   tatty=1
 end if

 if (hdr%usepaw/=hdr0%usepaw) then
   write(message, '(a,a,a,a,i6,a,i6)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  input usepaw=',hdr%usepaw,' not equal disk file usepaw=',hdr0%usepaw
   call wrtout(std_out,message,mode_paral)
   tpaw=1
 end if

 if (hdr%usewvl/=hdr0%usewvl) then
   write(message, '(a,a,a,a,i6,a,i6,a,a)' ) ch10,&
&   ' hdr_check: ERROR -',ch10,&
&   '  input usewvl=',hdr%usewvl,' not equal disk file usewvl=',hdr0%usewvl, ch10, &
&   '  Action: change usewvl input variable or your restart file.'
   call wrtout(std_out,message,mode_paral)
   call leave_new(mode_paral)
 end if

!Also examine agreement of floating point data

 if (hdr%usewvl == 0 .and. abs(hdr%ecut_eff-hdr0%ecut_eff)>tol8) then
   write(message, '(a,a,a,a,f12.6,a,f12.6,a)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  input ecut_eff=',hdr%ecut_eff,' /= disk file ecut_eff=',hdr0%ecut_eff,'.'
   call wrtout(std_out,message,mode_paral)
   tecut=1
 end if

 do ii=1,3
   do jj=1,3
     if (abs(hdr%rprimd(ii,jj)-hdr0%rprimd(ii,jj))>tol6) then
       write(message, '(a,a,a,a,i1,a,i1,a,1p,e17.9,a,i1,a,i1,a,e17.9)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  input rprimd(',ii,',',jj,')=',hdr%rprimd(ii,jj),&
&       ' /= disk file rprimd(',ii,',',jj,')=',hdr0%rprimd(ii,jj)
       call wrtout(std_out,message,mode_paral)
       tprim=1
     end if
   end do
 end do

!Below this point many comparisons only make sense if
!certain things agree, e.g. nkpt, natom.  Also have to
!accomodate different amounts of data in general.

 if (hdr%usepaw==1 .and. hdr0%usepaw==1) then

!  Compare ecutdg (PAW)
   write(message, '(a,f12.6,15x,a,a,f12.6)' )&
&   '  PAW: ecutdg',hdr %ecutdg,'|',&
&   '  PAW: ecutdg',hdr0%ecutdg
   call wrtout(std_out,message,mode_paral)
   if (hdr%ecutdg/=hdr0%ecutdg) then
     write(message, '(a,a,a,a,f12.6,a,f12.6)' ) ch10,&
&     ' hdr_check: WARNING -',ch10,&
&     '  input ecutdg=',hdr%ecutdg,&
&     '  not equal disk file ecutdg=',hdr0%ecutdg
     call wrtout(std_out,message,mode_paral)
     tdg=1
   end if
 end if

!Compare nband(nkpt*nsppol) (cannot compare if nkpt and nsppol not same)
 if (hdr%nkpt==hdr0%nkpt .and. hdr%nsppol==hdr0%nsppol) then
   nkpt=hdr%nkpt ; nsppol=hdr%nsppol
   write(message,'(a,32x,a,a)') '  nband:','|','  nband:'
   call wrtout(std_out,message,mode_paral)
   do istart = 1,nsppol*nkpt,9
     istop = min(istart + 8,nsppol*nkpt)
     mu = istop - istart + 1
!    generate a format specifier
     bndfmt(5:5) = number(mu)
     bndfmt(21:21) = number(mu)
     write(message,fmt=bndfmt) hdr%nband(istart:istop),'|',&
&     hdr0%nband(istart:istop)
     call wrtout(std_out,message,mode_paral)
   end do

   do isppol=1,nsppol
     do ii=1,nkpt
       if (hdr%nband(ii)/=hdr0%nband(ii)) then
         tband=1
         if (fform == 2) then
           write(message,'(a,a,a,a,i5,a,i8,a,i8)' ) ch10,&
&           ' hdr_check: WARNING -',ch10,&
&           '  kpt num',ii,' input nband=',hdr%nband(ii),&
&           ' not equal disk file nband=',hdr0%nband(ii)
           call wrtout(std_out,message,mode_paral)
         end if
       end if
     end do
   end do
 end if

!Compare the number of wavelets in each resolution.
 if (hdr%usewvl == 1) then
   if (size(hdr%nwvlarr) /= size(hdr0%nwvlarr) .or. size(hdr%nwvlarr) /= 2) then
     write(message, '(a,a,a,a,i6,a,i6,a,a)' ) ch10,&
&     ' hdr_check: ERROR -',ch10,&
&     '  input nwvlres=',size(hdr%nwvlarr), &
&     ' not equal disk file nwvlres=',size(hdr0%nwvlarr), ' or 2',&
&     '  ABINIT is not implemented for wavelet resolutions different from 2.'
     call wrtout(std_out,message,mode_paral)
     call leave_new(mode_paral)
   end if
   write(message,'(a,30x,a,a)') '  nwvlres:','|', '  nwvlres:'
   call wrtout(std_out,message,mode_paral)
   write(message,'(a,2I6,24x,a,a,2I6)') '    ',hdr%nwvlarr, '|','    ', hdr0%nwvlarr
   call wrtout(std_out,message,mode_paral)
   do ii=1,2
     if (hdr%nwvlarr(ii)/=hdr0%nwvlarr(ii)) then
       twvl=1
       write(message,'(a,a,a,a,i5,a,i8,a,i8)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  nwvl resolution',ii,' input =',hdr%nwvlarr(ii),&
&       ' not equal disk file value=',hdr0%nwvlarr(ii)
       call wrtout(std_out,message,mode_paral)
     end if
   end do
 end if

!Compare symmetry arrays (integers) symafm(nsym)
!-- only for same number of symmetries nsym
 itest=0
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   write(message,'(a,31x,a,a)') '  symafm:','|','  symafm:'
   call wrtout(std_out,message,mode_paral)
   do istart = 1,nsym,12
     istop=min(istart+11,nsym)
     nelm = istop - istart + 1
     call mk_hdr_check_fmt(nelm,typfmt)
     write(message,fmt=typfmt) hdr%symafm(istart:istop),&
&     '|',hdr0%symafm(istart:istop)
     call wrtout(std_out,message,mode_paral)
   end do
 end if

 if (itest/=0) then
   write(message,'(a,a,a,a,i6,a)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  For symmetry number',itest,' input symafm not equal disk file symafm'
   call wrtout(std_out,message,mode_paral)
   tsym=1
 end if

!Compare symmetry arrays (integers) symrel(3,3,nsym)
!-- only for same number of symmetries nsym
 itest=0
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   write(message,'(a,31x,a,a)') '  symrel:','|','  symrel:'
   call wrtout(std_out,message,mode_paral)
   do isym=1,nsym
     write(message,'(2x,9i3,11x,a,2x,9i3)') &
&     hdr%symrel(:,:,isym),'|',hdr0%symrel(:,:,isym)
     call wrtout(std_out,message,mode_paral)
     if(sum(abs(hdr%symrel(:,:,isym)-hdr0%symrel(:,:,isym)))/=0)then
       itest=isym
       exit
     end if
   end do
 end if

 if (itest/=0) then
   write(message,'(a,a,a,a,i6,a)' ) ch10,&
&   ' hdr_check: WARNING -',ch10,&
&   '  For symmetry number',itest,' input symrel not equal disk file symrel'
   call wrtout(std_out,message,mode_paral)
   tsym=1
 end if

!Compare typat(natom)
 if (hdr%natom==hdr0%natom) then
   natom=hdr%natom
   write(message,'(a,32x,a,a)') '  typat:','|','  typat:'
   call wrtout(std_out,message,mode_paral)
   do istart = 1,natom,12
     istop=min(istart+11,natom)
     nelm = istop - istart + 1
     call mk_hdr_check_fmt(nelm,typfmt)
     write(message,fmt=typfmt) hdr%typat(istart:istop),'|',hdr0%typat(istart:istop)
     call wrtout(std_out,message,mode_paral)
   end do
   do ii=1,natom
     if (hdr%typat(ii)/=hdr0%typat(ii)) then
       write(message, '(a,a,a,a,i8,a,i3,a,i3)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  For atom number',ii,' input typat=',hdr%typat(ii),&
&       ' not equal disk file typat=',hdr0%typat(ii)
       call wrtout(std_out,message,mode_paral)
       tatty=1
     end if
   end do
 end if


!Compare so_psp(npsp)
 if (hdr%npsp==hdr0%npsp) then
   npsp=hdr%npsp
   write(message,'(a,29x,a,a)') '  so_psp  :','|','  so_psp  :'
   call wrtout(std_out,message,mode_paral)
   do istart = 1,npsp  ,12
     istop=min(istart+11,npsp  )
     nelm = istop - istart + 1
     call mk_hdr_check_fmt(nelm,typfmt)
     write(message,fmt=typfmt) hdr%so_psp  (istart:istop),'|',hdr0%so_psp  (istart:istop)
     call wrtout(std_out,message,mode_paral)
   end do
   do ii=1,npsp
     if (hdr%so_psp  (ii)/=hdr0%so_psp  (ii)) then
       write(message, '(a,a,a,a,i8,a,i3,a,i3)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  For pseudopotential number',ii,' input so_psp  =',hdr%so_psp(ii),&
&       ' not equal disk file so_psp=',hdr0%so_psp(ii)
       call wrtout(std_out,message,mode_paral)
     end if
   end do
 end if


!Compare istwfk(nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt
   write(message,'(a,31x,a,a)') '  istwfk:','|','  istwfk:'
   call wrtout(std_out,message,mode_paral)
   do istart = 1,nkpt,12
     istop=min(istart+11,nkpt)
     nelm = istop - istart + 1
     call mk_hdr_check_fmt(nelm,typfmt)
     write(message,fmt=typfmt) hdr%istwfk(istart:istop),'|',hdr0%istwfk(istart:istop)
     call wrtout(std_out,message,mode_paral)
   end do
   do ii=1,nkpt
     if (hdr%istwfk(ii)/=hdr0%istwfk(ii)) then
       write(message, '(a,a,a,a,i8,a,i3,a,i3)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  For k point number',ii,' input istwfk=',hdr%istwfk(ii),&
&       ' not equal disk file istwfk=',hdr0%istwfk(ii)
       call wrtout(std_out,message,mode_paral)
       twfk=1
     end if
   end do
 end if

!Compare kpt(3,nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt
   write(message,'(a,34x,a,a)') '  kpt:','|','  kpt:'
   call wrtout(std_out,message,mode_paral)
   do ii = 1,min(nkpt,nkpt_max)
     write(message,'(2x,3f12.7,2x,a,2x,3f12.7)')&
&     hdr%kptns(:,ii),'|',hdr0%kptns(:,ii)
     call wrtout(std_out,message,mode_paral)
     if(ii>nkpt_max)then
       write(message,'(a)')'  The number of printed k points is sufficient ... stop writing them.'
       call wrtout(std_out,message,mode_paral)
       exit
     end if
   end do
   iwarning=0
   do ii=1,nkpt
     itest=0
     do mu=1,3
       if(abs( hdr%kptns(mu,ii)-hdr0%kptns(mu,ii) )>tol6)itest=1
     end do
     if (itest==1) then
       write(message, '(a,a,a,a,i5,a,3es17.7,a,a,3es17.7)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  kpt num',ii,', input kpt=',hdr%kptns(:,ii),ch10,&
&       '  not equal  disk file kpt=',hdr0%kptns(:,ii)
       call wrtout(std_out,message,mode_paral)
       tkpt=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         write(message,'(a,a,a,a,a)') ch10,&
&         ' hdr_check: WARNING -',ch10,&
&         '  The number of warning messages is sufficient ... stop writing them.',ch10
         call wrtout(std_out,message,mode_paral)
         exit
       end if
     end if
   end do
 end if

!Compare wtk(nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt

   write(message,'(a,34x,a,a)') '  wtk:','|','  wtk:'
   call wrtout(std_out,message,mode_paral)
   istop = min(nkpt,nkpt_max)
   do ii = 1, istop, 5
     mu = min(5, istop - ii + 1)
     wtkfmt(5:5) = number(mu)
     wtkfmt(23:23) = number(mu)
     write(message, wtkfmt)&
&     hdr%wtk(ii:min(istop, ii + 5 - 1)),'|',hdr0%wtk(ii:min(istop, ii + 5 - 1))
     call wrtout(std_out,message,mode_paral)
   end do
   iwarning=0
   do ii=1,nkpt
     itest=0
     if (abs( hdr%wtk(ii)-hdr0%wtk(ii) )>tol6) then
       write(message, '(a,a,a,a,i5,a,es17.7,a,a,es17.7)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  kpt num',ii,', input weight=',hdr%wtk(ii),ch10,&
&       '  not equal  disk file weight=',hdr0%wtk(ii)
       call wrtout(std_out,message,mode_paral)
       tkpt=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         write(message,'(a,a,a,a,a)') ch10,&
&         ' hdr_check: WARNING -',ch10,&
&         '  The number of warning messages is sufficient ... stop writing them.',ch10
         call wrtout(std_out,message,mode_paral)
         exit
       end if
     end if
   end do
 end if

!Compare occ(bantot)
 if (hdr%nkpt==hdr0%nkpt.and. hdr%bantot==hdr0%bantot) then
   nkpt=hdr%nkpt
   bantot=hdr%bantot

   write(message,'(a,34x,a,a)') '  occ:','|','  occ:'
   call wrtout(std_out,message,mode_paral)
   bantot_eff=min(bantot,9*nkpt_max)
   do istart = 1,bantot_eff,9
     istop = min(istart+8,bantot_eff)
     mu = istop - istart + 1
     occfmt(5:5) = number(mu)
     occfmt(23:23) = number(mu)
     write(message,fmt=occfmt) &
&     hdr%occ(istart:istop),'|', hdr0%occ(istart:istop)
     call wrtout(std_out,message,mode_paral)
     if(istart>9*nkpt_max)then
       write(message,'(a)')&
&       '  The number of printed occupation numbers is sufficient ... stop writing them.'
       call wrtout(std_out,message,mode_paral)
       exit
     end if
   end do
   iwarning=0
   do ii=1,bantot
     if (abs( hdr%occ(ii)-hdr0%occ(ii) )>tol6) then
       write(message, '(a,a,a,a,i10,a,1p,e15.7,a,e15.7)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  band,k',ii,', input occ=',hdr%occ(ii),&
&       ' disk occ=',hdr0%occ(ii)
       call wrtout(std_out,message,mode_paral)
       tband=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         write(message,'(a,a,a,a,a)') ch10,&
&         ' hdr_check: WARNING -',ch10,&
&         '  The number of warning messages is sufficient ... stop writing them.',ch10
         call wrtout(std_out,message,mode_paral)
         exit
       end if
     end if
   end do
 end if

!Compare tnons(3,nsym)
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   itest=0
   write(message,'(a,32x,a,a)') '  tnons:','|','  tnons:'
   call wrtout(std_out,message,mode_paral)
   do isym=1,nsym
     write(message,'(2x,3f12.7,2x,a,2x,3f12.7)') hdr%tnons(:,isym),'|',hdr0%tnons(:,isym)
     call wrtout(std_out,message,mode_paral)
   end do

   do isym=1,nsym
     if( sum(abs(  hdr%tnons(:,isym)-hdr0%tnons(:,isym) )) > tol6) then
       itest=isym
       exit
     end if
   end do
   if (itest/=0) then
     write(message, '(a,a,a,a,i6,a)' ) ch10,&
&     ' hdr_check: WARNING -',ch10,&
&     '  For symmetry number',itest,' input tnons not equal disk file tnons'
     call wrtout(std_out,message,mode_paral)
   end if
 end if

!Compare znucltypat(ntypat)
 if (hdr%ntypat==hdr0%ntypat) then
   ntypat=hdr%ntypat

   write(message,'(a,31x,a,a)') '   znucl:','|','   znucl:'
   call wrtout(std_out,message,mode_paral)
   do istart = 1,ntypat,6
     istop = min(istart+5,ntypat)
     mu = istop-istart+1
     zatfmt(5:5) = number(mu)
     zatfmt(23:23) = number(mu)
     write(message,fmt=zatfmt) hdr%znucltypat(istart:istop),'|',hdr0%znucltypat(istart:istop)
     call wrtout(std_out,message,mode_paral)
   end do

   do ii=1,ntypat
     if (abs(hdr%znucltypat(ii)-hdr0%znucltypat(ii))>tol6) then
       write(message, '(a,a,a,a,i5,a,f12.6,a,f12.6)' ) ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       ' For atom number',ii,' input znucl=',hdr%znucltypat(ii),&
&       ' not equal disk file znucl=',hdr0%znucltypat(ii)
       call wrtout(std_out,message,mode_paral)
     end if
   end do
 end if

!Should perform some checks related to pertcase and qptn,
!that have been introduced in the header in v4.1
!Warning : a GS file might be read, while the hdr corresponds
!to a RF file (to initialize k+q), and vice-versa (in nonlinear).

!DEBUG
!write(std_out,*)' hdr_check : consistency of floating point data checked '
!stop
!ENDDEBUG

!Now check agreement of psp headers too
 if (hdr%npsp==hdr0%npsp) then
   npsp=hdr%npsp
   itest=0

   do ipsp=1,npsp

     write(message,'(a,i3,a,9x,a,a,i3,a)')&
&     '  pseudopotential atom type',ipsp,':','|',&
&     '  pseudopotential atom type',ipsp,':'
     call wrtout(std_out,message,mode_paral)
     if (hdr%usepaw==1 .and. hdr0%usepaw==1) then
       write(message,'(a,i3,a,i3,a,i3,5x,a,a,i3,a,i3,a,i3)')&
&       '  pspso ',hdr %pspso(ipsp),' pspxc ',hdr %pspxc(ipsp),&
&       '  lmn_size ',hdr%lmn_size(ipsp),'|',&
&       '  pspso ',hdr0%pspso(ipsp),' pspxc ',hdr0%pspxc(ipsp),&
&       '  lmn_size ',hdr0%lmn_size(ipsp)
       call wrtout(std_out,message,mode_paral)
       if (hdr%lmn_size(ipsp)/=hdr0%lmn_size(ipsp)) then
         write(message, '(a,a,a,a,i3,a,i3,a,i3)' ) ch10,&
&         ' hdr_check: WARNING -',ch10,&
&         '  For atom type ',ipsp,' input lmn_size=',hdr%lmn_size(ipsp),&
&         ' not equal disk file lmn_size=',hdr0%lmn_size(ipsp)
         call wrtout(std_out,message,mode_paral)
         tlmn=1
       end if
     else
       write(message,'(a,i3,a,i3,19x,a,a,i3,a,i3)')&
&       '  pspso ',hdr %pspso(ipsp),' pspxc ',hdr %pspxc(ipsp),'|',&
&       '  pspso ',hdr0%pspso(ipsp),' pspxc ',hdr0%pspxc(ipsp)
       call wrtout(std_out,message,mode_paral)
     end if
     write(message,'(a,i6,a,i4,a,f5.1,2x,a,a,i6,a,i4,a,f5.1)')&
&     '  pspdat ',hdr %pspdat(ipsp),' pspcod ',hdr %pspcod(ipsp),&
&     ' zion ',hdr %zionpsp(ipsp),'|',&
&     '  pspdat ',hdr0%pspdat(ipsp),' pspcod ',hdr0%pspcod(ipsp),&
&     ' zion ',hdr0%zionpsp(ipsp)
     call wrtout(std_out,message,mode_paral)

!    Second, test
!    NOTE, XG 000719 : should do something about pspso
!    NOTE, XG 020716 : znucl and zion are not written
     if (abs(hdr%znuclpsp(ipsp)-hdr0%znuclpsp(ipsp))>tol6) itest=1
     if (abs(hdr%zionpsp(ipsp)-hdr0%zionpsp(ipsp))>tol6) then
       itest=1
       tpsch=1
     end if
     if (hdr%pspdat(ipsp)/= hdr0%pspdat(ipsp)) itest=1
     if (hdr%pspcod(ipsp)/= hdr0%pspcod(ipsp)) itest=1
     if (hdr%pspxc(ipsp) /= hdr0%pspxc(ipsp) )  itest=1
   end do

   if (itest==1) then
     write(message,'(a,a,a,a)') ch10,&
&     ' hdr_check: WARNING -',ch10,&
&     '  input psp header does not agree perfectly with disk file psp header.'
     call wrtout(std_out,message,mode_paral)
     tpseu=1
   end if
 end if

!DEBUG
!write(std_out,*)' hdr_check : consistency of psp header checked '
!stop
!ENDDEBUG

!Finally, read residm and etotal ("current value" not known), and check xred.

 if (hdr%natom==hdr0%natom) then

   natom=hdr%natom
   write(message,'(a,33x,a,a)') '  xred:','|','  xred:'
   call wrtout(std_out,message,mode_paral)
   do ii=1,natom
     write(message,'(2x,3f12.7,2x,a,2x,3f12.7)') hdr%xred(:,ii),'|',hdr0%xred(:,ii)
     call wrtout(std_out,message,mode_paral)
   end do

!  check atom positions one atom at a time and allow possibility
!  that there is a harmless translation of atoms by a cell vector.
   do ii=1,natom
     rms=0.0_dp
     do jj=1,3
       rms=rms+(             hdr%xred(jj,ii)-hdr0%xred(jj,ii) &
&       - dble(nint((hdr%xred(jj,ii)-hdr0%xred(jj,ii)))) )**2
     end do
     rms=sqrt(rms/3.0_dp)
     if (rms>tol6) txred=1
   end do
 end if

!DEBUG
!write(std_out,*)' hdr_check : resdim and atom position checked '
!stop
!ENDDEBUG

!Run tests here to establish whether this is a valid restart

!tfform2 will be true if there is a problem for the wavefunctions
 tfform2 = (hdr%usewvl == 0 .and. &
& (tprim /= 0 .or. tecut /= 0 .or. tkpt /= 0 .or. &
& twfk /=0 .or. tspinor /= 0)) .or. &
& (hdr%usepaw == 1 .and. &
& (tpaw /= 0 .or. tlmn /= 0 .or. tdg /= 0)) .or. &
& (hdr%usewvl == 1 .and. &
& (tatty /= 0 .or. tband /= 0))
!tfform52 will be true if there is a problem for the format 52
 tfform52=tprim /= 0 .or. tatty /= 0 .or. txred /= 0 .or.&
& tpseu /= 0 .or. tecut /= 0 .or. tng /= 0 .or. &
& (hdr%usepaw == 1 .and. &
& (tpaw /= 0 .or. tlmn /= 0 .or. tdg /= 0))

 restart=1
 restartpaw=hdr%usepaw

!If there is a problem somewhere
 if ( (fform == 2  .and. tfform2  ) .or.  &
& (fform == 52 .and. tfform52 ) .or.  &
& (fform == 200 .and. tfform2 ) ) then

   if(fform==2)then
     restart=2
     write(message,'(a,a,a,a)') ch10,&
&     ' hdr_check: WARNING -',ch10,&
&     '  Restart of self-consistent calculation need translated wavefunctions.'
   else if(fform==52)then
     restart=0
     write(message,'(a,a,a,a)') ch10,&
&     ' hdr_check: WARNING -',ch10,&
&     '  Illegal restart of non-self-consistent calculation'
   end if
   call wrtout(std_out,message,mode_paral)

   write(message,'(a,a1,a)') &
&   '  Indeed, critical differences between current calculation and',&
&   ch10,'  restart file have been detected in:'
   call wrtout(std_out,message,mode_paral)

   if ( (fform==52 .or. fform == 200) .and. tatty /= 0 ) then
     write(message, '(8x,a)' ) '* the number of atoms of each type'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( fform /= 200 .and. tecut /= 0 ) then
     write(message, '(8x,a)' ) '* the plane-wave cutoff'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( fform == 200 .and. tband /= 0 ) then
     write(message, '(8x,a)' ) '* the band and their occupation'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( fform==2 .and. tkpt /= 0 ) then
     write(message, '(8x,a)' ) '* the number, position, or weight of k-points'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( fform==2 .and. twfk /= 0 ) then
     write(message, '(8x,a)' ) '* the format of wavefunctions (istwfk)'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( fform==2 .and. tspinor /= 0 ) then
     write(message, '(8x,a)' ) '* the scalar/spinor character of the wf (nspinor)'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( fform==52 .and. tng /= 0 ) then
     write(message, '(8x,a)' ) '* the Fourier transform box dimensions'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( tprim /= 0 ) then
     write(message, '(8x,a)' ) &
&     '* the vectors defining the unit cell (obtained from rprim and acell)'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( fform==52 .and. tpseu /= 0 ) then
     write(message, '(8x,a)' ) &
&     '* the pseudopotential files'
     call wrtout(std_out,message,mode_paral)
   end if
   if ( fform==52 .and. txred /= 0 ) then
     write(message, '(8x,a)' ) '* the positions of the ions in the basis'
     call wrtout(std_out,message,mode_paral)
   end if

!  Tests for a restart in the framework of the PAW method
   if (hdr%usepaw/=0 .or. hdr0%usepaw/=0) then
     if (tpaw /= 0 .or. tlmn /= 0) restartpaw=0
     if (restartpaw == 0) then
       write(message,'(8x,a)') 'Critical differences for a restart within PAW method:'
       call wrtout(std_out,message,mode_paral)
       if ( tpaw /= 0 ) then
         write(message, '(8x,a)' ) '* the use of the PAW method'
         call wrtout(std_out,message,mode_paral)
       else
         if(tlmn/=0)then
           write(message, '(8x,a)' ) '* the number of lmn elements for the paw basis'
           call wrtout(std_out,message,mode_paral)
         end if
       end if
     else if (tdg/=0) then
       write(message,'(a,a,a,a,a,a)') ch10,&
&       ' hdr_check: WARNING -',ch10,&
&       '  Restart of calculation within PAW may be inconsistent because of:"'
       call wrtout(std_out,message,mode_paral)
       if(tdg/=0)then
         write(message, '(8x,a)' ) &
&         '* the cutoff energy of the paw double (fine) grid'
         call wrtout(std_out,message,mode_paral)
       end if
     end if
   end if

 else

   if(fform==2 .or. fform == 200)then
     write(message,'(a,a)') ' hdr_check: ',&
&     ' Wavefunction file is OK for direct restart of calculation'
     call wrtout(std_out,message,mode_paral)
   else if(fform==52)then
     write(message,'(a,a)') ' hdr_check: ',&
&     ' Density/Potential file is OK for restart of calculation'
     call wrtout(std_out,message,mode_paral)
   end if
!  MG TODO add screening case! 
!  call wrtout(std_out,message,mode_paral)

 end if

 write(message,'(80a)') ('=',ii=1,80)
 call wrtout(std_out,message,mode_paral)

 CONTAINS
!!***

!!****f* hdr_check/mk_hdr_check_fmt
!! NAME
!! mk_hdr_check_fmt
!!
!! FUNCTION
!! make a format needed in hdr_check, for arrays of nint integers each of format i3
!!
!! INPUTS
!!  nelm=number of elements to be printed
!!
!! OUTPUT
!!  character(len=26), typfmt= format needed
!!
!! PARENTS
!!      hdr_check
!!
!! CHILDREN
!!
!! SOURCE

   subroutine mk_hdr_check_fmt(nelm,typfmt)

 use m_profiling

   use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mk_hdr_check_fmt'
!End of the abilint section

   implicit none

!  Arguments ------------------------------------
!  scalars
   integer,intent(in) :: nelm
   character(len=26),intent(out) :: typfmt

!  Local variables-------------------------------
!  scalars
   integer :: ii
   character(len=1), parameter :: number(0:10)=(/'0','1','2','3','4','5','6','7','8','9',' '/)
   character(len=26), parameter :: templatefmt='(2x,  i3,t41   ,a,2x,  i3)'
!  *************************************************************************

!  Initialize the format
   typfmt=templatefmt

!  Generate the type format specifier
   ii=nelm/10
   if ( ii /= 0 ) then
     typfmt(5:5) = number(ii)
     typfmt(22:22) = number(ii)
   else
     typfmt(5:5) = ' '
     typfmt(22:22) = ' '
   end if
   ii = nelm - 10 * (nelm/10)
   typfmt(6:6) = number(ii)
   typfmt(23:23) = number(ii)

!  DEBUG
!  write(std_out,*)' mk_hdr_check_fmt : exit'
!  write(std_out,*)' mk_hdr_check_fmt : typfmt="',typfmt,'"'
!  stop
!  ENDDEBUG

 end subroutine mk_hdr_check_fmt

end subroutine hdr_check
!!***
