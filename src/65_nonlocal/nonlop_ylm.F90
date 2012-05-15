!{\src2tex{textfont=tt}}
!!****f* ABINIT/nonlop_ylm
!! NAME
!! nonlop_ylm
!!
!! FUNCTION
!! * Compute application of a nonlocal operator Vnl in order to get:
!!    - contracted elements (energy, forces, stresses, ...), if signs=1
!!    - a function in reciprocal space (|out> = Vnl|in>), if signs=2
!!   Operator Vnl, as the following general form:
!!    $Vnl=sum_{R,lmn,l''m''n''} {|P_{Rlmn}> Enl^{R}_{lmn,l''m''n''} <P_{Rl''m''n''}|}$
!!   Operator Vnl is -- in the typical case -- the nonlocal potential.
!!   - With norm-conserving pseudopots, $Enl^{R}_{lmn,l''m''n''}$ is the
!!     Kleinmann-Bylander energy $Ekb^{R}_{ln}$.
!!   - In a PAW calculation, $Enl^{R}_{lmn,l''m''n''}$ are the nonlocal
!!     coefficients to connect projectors $D_{ij}$.
!!   - The |P_{Rlmn}> are the projector functions.
!! * Optionnaly, in case of PAW calculation, compute:
!!   - Application of the overlap matrix in reciprocal space
!!     (<in|S|in> or (I+S)|in>).
!!   - Application of (Vnl-lambda.S) in reciprocal space
!!     (<in|Vnl-lambda.S|in> and derivatives or (Vnl-lambda.S)|in>).
!! * This routine uses spherical harmonics Ylm to express Vnl.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  choice: chooses possible output:
!!    choice=0 => do nothing (only compute WF projected with NL projectors)
!!          =1 => a non-local energy contribution
!!          =2 => a gradient with respect to atomic position(s)
!!          =3 => a gradient with respect to strain(s)
!!          =23=> a gradient with respect to atm. pos. and strain(s)
!!          =4 => a gradient and 2nd derivative with respect to atomic pos.
!!          =24=> a gradient and 2nd derivative with respect to atomic pos.
!!          =5 => a gradient with respect to k wavevector, typically
!!                $\sum_{ij}|p_i\rangle D_{ij}\langle dp_j/dk| + |dp_i/dk\rangle D_{ij}\langle p_j|$
!!          =51 => the right derivative with respect to k wavevector, typically
!!                $\sum_{ij}|p_i\rangle D_{ij}\langle dp_j/dk|$
!!          =52 => the left derivative with respect to k wavevector, typically
!!                $\sum_{ij}|dp_i/dk\rangle D_{ij}\langle p_j|$
!!          =53 => the twist derivative with respect to k, typically
!!                $\sum_{ij}|dp_i/dk_(idir+1)\rangle D_{ij}\langle dp_j/dk_(idir-1)| -
!!                 |dp_i/dk_(idir-1)\rangle D_{ij}\langle dp_j/dk_(idir+1)|$
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  cpopt=flag defining the status of cprjin%cp(:)=<Proj_i|Cnk> scalars (see below, side effects)
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  dimffnlin=second dimension of ffnlin (1+number of derivatives)
!!  dimffnlout=second dimension of ffnlout (1+number of derivatives)
!!  enl(dimenl1,dimenl2,nspinortot**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1, 2 or 4 ====
!!                      (Real or complex, hermitian) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  ffnlin(npwin,dimffnlin,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |in> vector
!!  ffnlout(npwout,dimffnlout,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |out> vector
!!  ---- Taken away in beautification because unused MS -------
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  -----------------------------------------------------------
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                          for choice 53, twisted derivative involves idir+1 and idir-1
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kgin(3,npwin)=integer coords of planewaves in basis sphere, for the |in> vector
!!  kgout(3,npwout)=integer coords of planewaves in basis sphere, for the |out> vector
!!  kpgin(npw,npkgin)= (k+G) components and related data, for the |in> vector
!!  kpgout(npw,nkpgout)=(k+G) components and related data, for the |out> vector
!!  kptin(3)=k point in terms of recip. translations, for the |in> vector
!!  kptout(3)=k point in terms of recip. translations, for the |out> vector
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!         Typically lambda is the eigenvalue (or its guess)
!!  lmnmax=max. number of (l,m,n) components over all types of atoms
!!  matblk=dimension of the arrays ph3din and ph3dout
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpgin,nkpgout=second sizes of arrays kpgin/kpgout
!!  nloalg(5)=governs the choice of the algorithm for nonlocal operator
!!  nnlout=dimension of enlout (when signs=1 and choice>0):
!!         ==== if paw_opt=0, 1 or 2 ====
!!         choice=1=>nnlout=1           choice=2=>nnlout=3*natom
!!         choice=3=>nnlout=6           choice=4=>nnlout=6*natom
!!         choice=5=>nnlout=3           choice=6=>nnlout=6*(3*natom+6)
!!         choice=23=>nnlout=6+3*natom  choice=24=>nnlout=9*natom
!!         ==== if paw_opt=3 ====
!!         choice=1 =>nnlout=1
!!         ==== if paw_opt=4 ====
!!         not available
!!  npwin=number of planewaves for given k point, for the |in> vector
!!  npwout=number of planewaves for given k point, for the |out> vector
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nspinortot=total number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  phkxredin(2,natom)=phase factors exp(2 pi kptin.xred)
!!  phkxredout(2,natom)=phase factors exp(2 pi kptout.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
!!  ph3din(2,npwin,matblk)=3D structure factors, for each atom and plane wave (in)
!!  ph3dout(2,npwout,matblk)=3-dim structure factors, for each atom and plane wave (out)
!!  ---- Taken away in beautification because unused MS -------
!!  pspso(ntypat)=spin-orbit characteristic for each atom type
!!  -----------------------------------------------------------
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  sij(dimenl1,ntypat*(paw_opt/3))=overlap matrix components (only if paw_opt=2, 3 or 4)
!!  ucvol=unit cell volume (bohr^3)
!!  vectin(2,npwin*nspinor)=input cmplx wavefunction coefficients <G|Cnk>
!!  [cprjin_left(natom,nspinor)]=The projected input wave function <p_nlm|in_left>
!!    for the left wavefunction. Data are assumed to be in memory, they are NOT recalculated here.
!!    Only signs==1 and choice==1 are supported.
!!
!! OUTPUT
!! ==== if (signs==1) ====
!! --If (paw_opt==0, 1 or 2)
!!    enlout(nnlout)= contribution to the non-local part of the following properties:
!!      if choice=1 : enlout(1)               -> the energy
!!      if choice=2 : enlout(1:3*natom)       -> the forces
!!      if choice=3 : enlout(1:6)             -> the stresses
!!      if choice=23: enlout(1:6+3*natom)     -> the forces and the stresses
!!      if choice=4 : enlout(1:6*natom)       -> the frozen wf part of dyn. mat.
!!      if choice=5 : enlout(3)               -> the derivatives of energy wrt to k
!!      if choice=24: enlout(1:9*natom)       -> the forces and the frozen wf part of dyn. mat.
!!      if choice=6 : enlout(1:6*(3*natom+6)) -> the frozen wf part of elastic tensor
!! --If (paw_opt==3)
!!    if choice=1 : enlout(nnlout)= contribution to <c|S|c>  (nnlout=1)
!! --If (paw_opt==4)
!!    not available
!! ==== if (signs==2) ====
!! --if (paw_opt=0, 1 or 4)
!!    vectout(2,npwout*nspinor)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1) <G|V_nonlocal|vect_start>
!!      if (choice=2) <G|dV_nonlocal/d(atm coord)|vect_start>
!!      if (choice=3) <G|dV_nonlocal/d(strain)|vect_start>
!!      if (choice=5) <G|dV_nonlocal/dk|vect_start>
!!      if (choice=51) <G|d(right)V_nonlocal/dk|vect_start>
!!      if (choice=52) <G|d(left)V_nonlocal/dk|vect_start>
!!      if (choice=53) <G|d(twist)V_nonlocal/dk|vect_start>
!!  if (paw_opt=2)
!!    vectout(2,npwout*nspinor)=final vector in reciprocal space:
!!      if (choice=1) <G|V_nonlocal-lamdba.(I+S)|vect_start>
!!      if (choice=2) <G|d[V_nonlocal-lamdba.(I+S)]/d(atm coord)|vect_start>
!!      if (choice=3) <G|d[V_nonlocal-lamdba.(I+S)]/d(strain)|vect_start>
!!      if (choice=5) <G|d[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
!!      if (choice=51) <G|d(right)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
!!      if (choice=52) <G|d(left)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
!!      if (choice=53) <G|d(twist)[V_nonlocal-lamdba.(I+S)]/dk|vect_start>
!! --if (paw_opt=3 or 4)
!!    svectout(2,npwout*nspinor)=result of the aplication of Sij (overlap matrix)
!!                  or one of its derivatives to the input vect.:
!!      if (choice=1) <G|I+S|vect_start>
!!      if (choice=2) <G|dS/d(atm coord)|vect_start>
!!      if (choice=3) <G|dS/d(strain)|vect_start>
!!      if (choice=5) <G|dS/dk|vect_start>
!!      if (choice=51) <G|d(right)S/dk|vect_start>
!!      if (choice=52) <G|d(left)S/dk|vect_start>
!!      if (choice=53) <G|d(twist)S/dk|vect_start>
!!
!! SIDE EFFECTS
!!  cprjin(natom,nspinor) <type(cprj_type)>=projected input wave function |in> on non-local projectors
!!                                  =<p_lmn|in> and derivatives
!!                    Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  are computed here (and not saved)
!!                     if cpopt= 0, <p_lmn|in> are computed here and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 2  <p_lmn|in> are already in memory;
!!                                  first (and 2nd) derivatives are computed here and not saved
!!                     if cpopt= 3  <p_lmn|in> are already in memory;
!!                                  first derivatives are computed here and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!!                                  other derivatives are not computed
!!                                  This option is not compatible with choice=4,24 or 6
!!
!! NOTES
!! This application of the nonlocal operator is programmed using a direct
!! implementation of spherical harmonics (Ylm). Abinit used historically
!! Legendre polynomials for the application of nonlocal operator; but the
!! implementation of PAW algorithm enforced the use of Ylm.
!!
!! In the case signs=1, the array vectout is not used, nor modified
!! so that the same array as vectin can be used as a dummy argument;
!! the same is true for the pairs npwin-npwout, ffnlin-ffnlout,
!! kgin-kgout, ph3din-ph3dout, phkredin-phkxredout).
!!
!! TODO
!! * Implementation of spin-orbit
!!
!! PARENTS
!!      m_cprj_bspline,m_shirley,nonlop
!!
!! CHILDREN
!!      mkkpg,opernla_ylm,opernlb_ylm,opernlc_ylm,opernld_ylm,ph1d3d,strconv
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine nonlop_ylm(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimffnlin,dimffnlout,&
&                      enl,enlout,ffnlin,ffnlout,gprimd,idir,indlmn,istwf_k,&
&                      kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                      mpi_enreg,natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,nnlout,&
&                      npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,phkxredin,phkxredout,ph1d,&
&                      ph3din,ph3dout,signs,sij,svectout,ucvol,vectin,vectout,cprjin_left)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nonlop_ylm'
 use interfaces_42_geometry
 use interfaces_65_nonlocal, except_this_one => nonlop_ylm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimffnlin,dimffnlout,idir
 integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,natom,nkpgin,nkpgout,nnlout
 integer,intent(in) :: npwin,npwout,nspinor,nspinortot,ntypat,paw_opt,signs
 real(dp),intent(in) :: lambda,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
 integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(5)
! integer,intent(in) :: pspso(ntypat) unused
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2)
 real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat) !,gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin),kpgout(npwout,nkpgout)
 real(dp),intent(in) :: kptin(3),kptout(3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: phkxredin(2,natom),phkxredout(2,natom)
 real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
 real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
 real(dp),intent(inout) :: vectin(2,npwin*nspinor)
 real(dp),intent(out) :: enlout(nnlout)
 real(dp),intent(out) :: svectout(2,npwout*nspinor*(paw_opt/3))
 real(dp),intent(out) :: vectout (2,npwout*nspinor)
 type(cprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5))
 type(cprj_type),optional,intent(in) :: cprjin_left(natom,nspinor)

!Local variables-------------------------------
!scalars
 integer :: choice_,cplex,cplex_enl,cplex_fac,ia,ia1,ia2,ia3,ia4,ia5,iatm,ierr,ilmn,iln
 integer :: ispinor,itypat,mincat,mu,mua,mub,n1,n2,n3,nd2gxdt,ndgxdt,ndgxdtfac,nincat
 integer :: nkpgin_,nkpgout_,nlmn,nua1,nua2,nub1,nub2,optder
 real(dp) :: enlk
 logical :: testnl
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,allocatable :: indlmn_typ(:,:)
! real(dp) :: kptdum(3)
 real(dp),allocatable :: d2gxdt(:,:,:,:,:),dgxdt(:,:,:,:,:),dgxdtfac(:,:,:,:,:),dgxdtfac_sij(:,:,:,:,:)
 real(dp),allocatable :: enlout_tmp(:,:),ffnlin_typ(:,:,:),ffnlout_typ(:,:,:)
 real(dp),allocatable :: fnlk(:),gx(:,:,:,:),gxfac(:,:,:,:),gxfac_sij(:,:,:,:)
 real(dp),allocatable :: kpgin_(:,:),kpgout_(:,:),sij_typ(:),strnlk(:),work(:)
 real(dp),allocatable :: work_out(:),strnlk_out(:)
 real(dp),allocatable :: gx_left(:,:,:,:) !,dgxdt_left(:,:,:,:,:)
!no_abirules
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))

! **********************************************************************

!DEBUG
!write(std_out,*)' nonlop_ylm : enter ',choice,signs
!ENDDEBUG

!Error on bad choice
 if ((choice<0 .or. choice>6)  &
& .and. choice/=23 .and. choice/=24 &
& .and. choice/=51 .and. choice /= 52 &
& .and. choice/=53 ) then
   write(message,'(a,i4,a)') '  Does not presently support this choice=',choice,'.'
   MSG_BUG(message)
 end if
 if ((choice>=51 .and. choice<=52) .and. signs==1) then
   write(message,'(a,i4,a)') '  signs=1 does not presently support this choice=',choice,'.'
   MSG_BUG(message)
 end if
 if (choice==53 .and. (idir<1 .or. idir>3) .and. signs/=1 ) then
   message='  Choice 53, signs=2 requires 1<=idir<=3.'
   MSG_BUG(message)
 end if
 if (choice==53 .and. (dimffnlin/=4 .or. dimffnlout/=4) ) then
   message='  Choice 53 requires dimffnlin and dimffnlout = 4.'
   MSG_BUG(message)
 end if
 if (cpopt<-1.or.cpopt>4) then
   message='  Bad value for cpopt !'
   MSG_BUG(message)
 end if
 if (cpopt==4.and.(choice==4.or.choice==24.or.choice==6)) then
   message='  This value of cpopt is not allowed for second derivatives !'
   MSG_BUG(message)
 end if
 if (PRESENT(cprjin_left)) then
   if (signs/=1) then
     message='signs must be 1 when cprjin_left is present.'
     MSG_BUG(message)
   end if
   if (choice/=1) then
     message='choice must be 1 when cprjin_left is present.'
     MSG_BUG(message)
   end if
 end if

!Spin-orbit not yet allowed
 if (maxval(indlmn(6,:,:))>1) then
   message='  Spin-orbit not yet allowed.'
   MSG_ERROR(message)
 end if

!Test: size of blocks of atoms
 mincat=min(nloalg(4),maxval(nattyp))
 if (nloalg(1)<=0.and.mincat>matblk) then
   write(message, '(a,a,a,i4,a,i4,a)' ) &
&   '  With nloalg<=0, mincat must be less than matblk.',ch10,&
&   '  Their value is ',mincat,' and ',matblk,'.'
   MSG_BUG(message)
 end if

!Test: sizes of kpgin/kpgout
 if (nkpgin>0.and. &
& ( (choice==2.and.nkpgin<3) .or. &
& (choice==4.and.signs==1.and.nkpgin<9) .or. &
& ((choice==6.or.choice==3.or.choice==23).and.signs==1.and.nkpgin<3) )) then
   message='  Incorrect size for nkpgin array !'
   MSG_BUG(message)
 end if
 if (choice==2.and.signs==2.and.nkpgout>0.and.nkpgout<3) then
   message='  Incorrect size for nkpgout array !'
   MSG_BUG(message)
 end if

!Define some useful variables
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 choice_=choice;if (cpopt>=2) choice_=-choice
 cplex=2;if (istwf_k>1) cplex=1 !Take into account TR-symmetry
 cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1))
 cplex_fac=cplex;if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0) cplex_fac=2

!Define dimensions of projected scalars
 ndgxdt=0;ndgxdtfac=0;nd2gxdt=0
 if (choice==2) then
   if (signs==1) ndgxdt=3
   if (signs==2) ndgxdt=1
   if (signs==2) ndgxdtfac=1
 end if
 if (choice==3) then
   if (signs==1) ndgxdt=6
   if (signs==2) ndgxdt=1
   if (signs==2) ndgxdtfac=1
 end if
 if (choice==23) then
   if (signs==1) ndgxdt=9
 end if
 if (choice==4) then
   if(signs==1) ndgxdt=3
   if(signs==1) ndgxdtfac=3
   if(signs==1) nd2gxdt=6
 end if
 if (choice==24) then
   if(signs==1) ndgxdt=3
   if(signs==1) ndgxdtfac=3
   if(signs==1) nd2gxdt=6
 end if
 if (choice==5) then
   if(signs==1) ndgxdt=3
   if(signs==2) ndgxdt=1
   if(signs==2) ndgxdtfac=1
 end if
 if (choice==51) then
   if(signs==1) ndgxdt=3
   if(signs==2) ndgxdt=1
   if(signs==2) ndgxdtfac=1
 end if
 if (choice==52) then
   if(signs==1) ndgxdt=3
   if(signs==2) ndgxdt=1
   if(signs==2) ndgxdtfac=1
 end if
 if (choice==53) then
   if(signs==1) ndgxdt=3
   if(signs==1) ndgxdtfac=3
   if(signs==2) ndgxdt=2
   if(signs==2) ndgxdtfac=2
 end if
 if (choice==6) then
   if(signs==1) ndgxdt=9
   if(signs==1) ndgxdtfac=9
   if(signs==1) nd2gxdt=54
 end if
 optder=0;if (ndgxdtfac>0) optder=1

!Test: gradients of cprjin
 if (cpopt==4) then
   if (ndgxdt>0.and.cprjin(1,1)%ncpgr<=0) then
     message='  cprjin%ncpgr=0 not allowed with cpopt=4 and these (choice,signs) !'
     MSG_BUG(message)
   end if
 end if
 if (cpopt==1.or.cpopt==3) then
   if (cprjin(1,1)%ncpgr<ndgxdt) then
     message='  should have cprjin%ncpgr>=ndgxdt with cpopt=1 or 3 !'
     MSG_BUG(message)
   end if
 end if

!Initialize output arrays
 if (signs==1) then
   ABI_ALLOCATE(fnlk,(3*natom))
   ABI_ALLOCATE(strnlk,(6))
   ABI_ALLOCATE(strnlk_out,(6))
   enlout(:)=zero
   if (choice==3.or.choice==6.or.choice==23) enlk=zero
   if (choice==6) then
     fnlk=zero;strnlk=zero
   end if
 end if
 if (signs==2) then
   if (paw_opt==0.or.paw_opt==1.or.paw_opt==4) vectout(:,:)=zero
   if (paw_opt==2.and.choice==1) vectout(:,:)=-lambda*vectin(:,:)
   if (paw_opt==3.or.paw_opt==4) then
     if (choice==1) svectout(:,:)=vectin(:,:)
     if (choice/=1) svectout(:,:)=zero
   end if
 end if

!Eventually re-compute (k+G) vectors (and related data)
 nkpgin_=0;nkpgout_=0
 if (nkpgin==0) then
   if ((choice==4.or.choice==24).and.signs==1) nkpgin_=9
   if ((choice==2).or.((choice==6.or.choice==3.or.choice==23).and.signs==1)) nkpgin_=3
   if (nkpgin_>0) then
     ABI_ALLOCATE(kpgin_,(npwin,nkpgin_))
     call mkkpg(kgin,kpgin_,kptin,nkpgin_,npwin)
   end if
 end if
 if (nkpgout==0) then
   if (choice==2.and.signs==2) nkpgout_=3
   if (nkpgout_>0) then
     ABI_ALLOCATE(kpgout_,(npwout,nkpgout_))
     call mkkpg(kgout,kpgout_,kptout,nkpgout_,npwout)
   end if
 end if

!Big loop on atom types.
 ia1=1;iatm=0
 do itypat=1,ntypat

!  Get atom loop indices for different types:
   ia2=ia1+nattyp(itypat)-1;ia5=1

!  Select quantities specific to the current type of atom
   nlmn=count(indlmn(3,:,itypat)>0)

!  Temporary test on local part
   testnl=(paw_opt/=0)
   if (paw_opt==0) then
     do ispinor=1,nspinortot
       do ilmn=1,nlmn
         iln=indlmn(5,ilmn,itypat)
         if (abs(enl(iln,itypat,ispinor))>tol10) testnl=.true.
       end do
     end do
   end if

!  Some non-local part is to be applied for that type of atom
   if (testnl) then

!    Store some quantities depending only of the atom type
     ABI_ALLOCATE(indlmn_typ,(6,nlmn))
     ABI_ALLOCATE(ffnlin_typ,(npwin,dimffnlin,nlmn))
     do ilmn=1,nlmn
       indlmn_typ(:,ilmn)=indlmn(:,ilmn,itypat)
       ffnlin_typ(:,:,ilmn)=ffnlin(:,:,ilmn,itypat)
     end do
     if (signs==2) then
       ABI_ALLOCATE(ffnlout_typ,(npwout,dimffnlout,nlmn))
       do ilmn=1,nlmn
         ffnlout_typ(:,:,ilmn)=ffnlout(:,:,ilmn,itypat)
       end do
     end if
     if (paw_opt>=2) then
       ABI_ALLOCATE(sij_typ,(nlmn*(nlmn+1)/2))
       if (cplex_enl==1) then
         do ilmn=1,nlmn*(nlmn+1)/2
           sij_typ(ilmn)=sij(ilmn,itypat)
         end do
       else
         do ilmn=1,nlmn*(nlmn+1)/2
           sij_typ(ilmn)=sij(2*ilmn-1,itypat)
         end do
       end if
     end if

!    Cut the sum on different atoms in blocks, to allow memory saving.
!    Inner summations on atoms will be done from ia3 to ia4.
!    Note: the maximum range from ia3 to ia4 is mincat (max. increment of atoms).
     do ia3=ia1,ia2,mincat
       ia4=min(ia2,ia3+mincat-1)
!      Give the increment of number of atoms in this subset.
       nincat=ia4-ia3+1

!      Prepare the phase factors if they were not already computed
       if (nloalg(1)<=0) call ph1d3d(ia3,ia4,kgin,matblk,natom,npwin,&
&       n1,n2,n3,phkxredin,ph1d,ph3din)

!      Allocate memory for projected scalars
       ABI_ALLOCATE(gx,(cplex,nlmn,nincat,nspinor))
       ABI_ALLOCATE(dgxdt,(cplex,ndgxdt,nlmn,nincat,nspinor))
       ABI_ALLOCATE(d2gxdt,(cplex,nd2gxdt,nlmn,nincat,nspinor))
       ABI_ALLOCATE(dgxdtfac,(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor))
       ABI_ALLOCATE(gxfac,(cplex_fac,nlmn,nincat,nspinor))
       gx(:,:,:,:)=zero;gxfac(:,:,:,:)=zero
       if (ndgxdt>0) dgxdt(:,:,:,:,:)=zero
       if (ndgxdtfac>0) dgxdtfac(:,:,:,:,:)=zero
       if (nd2gxdt>0) d2gxdt(:,:,:,:,:)=zero
       if (paw_opt>=3) then
         ABI_ALLOCATE(gxfac_sij,(cplex,nlmn,nincat,nspinor))
         ABI_ALLOCATE(dgxdtfac_sij,(cplex,ndgxdtfac,nlmn,nincat,nspinor))
         gxfac_sij(:,:,:,:)=zero
         if (ndgxdtfac>0) dgxdtfac_sij(:,:,:,:,:)=zero
       end if

!      Compute projection of current wave function |c> on each
!      non-local projector: <p_lmn|c>
!      ==============================================================

!      Retrieve eventually <p_lmn|c> coeffs
       if (cpopt>=2) then
         do ispinor=1,nspinor
           do ia=1,nincat
             gx(1:cplex,1:nlmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)
           end do
         end do
       end if
       if (cpopt==4.and.ndgxdt>0) then
         do ispinor=1,nspinor
           do ia=1,nincat
             dgxdt(1:cplex,1:ndgxdt,1:nlmn,ia,ispinor)=cprjin(iatm+ia,ispinor)%dcp(1:cplex,1:ndgxdt,1:nlmn)
           end do
         end do
       end if

!      Computation or <p_lmn|c> (and derivatives) for this block of atoms
       if (cpopt<4.and.choice_/=-1) then
         if (nkpgin_>0) then
           call opernla_ylm(choice_,cplex,dimffnlin,d2gxdt,dgxdt,ffnlin_typ,gx,ia3,idir,indlmn_typ,&
&           istwf_k,kpgin_,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpgin_,nlmn,&
&           nloalg,npwin,nspinor,ph3din,signs,ucvol,vectin)
         else
           call opernla_ylm(choice_,cplex,dimffnlin,d2gxdt,dgxdt,ffnlin_typ,gx,ia3,idir,indlmn_typ,&
&           istwf_k,kpgin,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpgin,nlmn,&
&           nloalg,npwin,nspinor,ph3din,signs,ucvol,vectin)
         end if
       end if

!      Transfer result to output variable cprj (if requested)
       if (cpopt==0.or.cpopt==1) then
         do ispinor=1,nspinor
           do ia=1,nincat
             cprjin(iatm+ia,ispinor)%nlmn=nlmn
             cprjin(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)=gx(1:cplex,1:nlmn,ia,ispinor)
             if (cplex==1) cprjin(iatm+ia,ispinor)%cp(2,1:nlmn)=zero
           end do
         end do
       end if
       if ((cpopt==1.or.cpopt==3).and.ndgxdt>0) then
         do ispinor=1,nspinor
           do ia=1,nincat
             cprjin(iatm+ia,ispinor)%dcp(1:cplex,1:ndgxdt,1:nlmn)=dgxdt(1:cplex,1:ndgxdt,1:nlmn,ia,ispinor)
             if (cplex==1) cprjin(iatm+ia,ispinor)%dcp(2,1:ndgxdt,1:nlmn)=zero
           end do
         end do
       end if

!      If choice==0, that's all for these atoms !
       if (choice>0) then

!        Contraction from <p_i|c> to Sum_j[Dij.<p_j|c>] (and derivatives)
         call opernlc_ylm(atindx1,cplex,cplex_enl,cplex_fac,dgxdt,dgxdtfac,dgxdtfac_sij,&
&         dimenl1,dimenl2,enl,gx,gxfac,gxfac_sij,&
&         iatm,indlmn_typ,itypat,lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,&
&         nincat,nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ)


!        Operate with the non-local potential on the projected scalars,
!        in order to get contributions to energy/forces/stress/dyn.mat
!        ==============================================================
         if (signs==1) then
           if (.not.PRESENT(cprjin_left)) then
             call opernld_ylm(choice,cplex,cplex_fac,d2gxdt,dgxdt,dgxdtfac,&
&             enlk,enlout,fnlk,gx,gxfac,gxfac_sij,ia3,natom,nd2gxdt,ndgxdt,ndgxdtfac,&
&             nincat,nlmn,nnlout,nspinor,paw_opt,strnlk,ucvol)
           else
             ABI_ALLOCATE(gx_left,(cplex,nlmn,nincat,nspinor))
!            Retrieve <p_lmn|c> coeffs
!            if (cpopt>=2) then
             do ispinor=1,nspinor
               do ia=1,nincat
                 gx_left(1:cplex,1:nlmn,ia,ispinor)=cprjin_left(iatm+ia,ispinor)%cp(1:cplex,1:nlmn)
               end do
             end do
!            end if
!            TODO
!            if (cpopt==4.and.ndgxdt>0) then
!            do ispinor=1,nspinor
!            do ia=1,nincat
!            dgxdt_left(1:cplex,1:ndgxdt,1:nlmn,ia,ispinor)=cprjin_left(iatm+ia,ispinor)%dcp(1:cplex,1:ndgxdt,1:nlmn)
!            end do
!            end do
!            end if
             call opernld_ylm(choice,cplex,cplex_fac,d2gxdt,dgxdt,dgxdtfac,&
&             enlk,enlout,fnlk,gx_left,gxfac,gxfac_sij,ia3,natom,nd2gxdt,ndgxdt,ndgxdtfac,&
&             nincat,nlmn,nnlout,nspinor,paw_opt,strnlk,ucvol)
             ABI_DEALLOCATE(gx_left)
!            END MG
           end if

         end if

!        Operate with the non-local potential on the projected scalars,
!        in order to get matrix element
!        ==============================================================
         if (signs==2) then
!          Prepare the phase factors if they were not already computed
           if(nloalg(1)<=0) call ph1d3d(ia3,ia4,kgout,matblk,natom,npwout,&
&           n1,n2,n3,phkxredout,ph1d,ph3dout)
           if (nkpgout_>0) then
             call opernlb_ylm(choice,cplex,cplex_fac,dgxdtfac,dgxdtfac_sij,&
&             dimffnlout,ffnlout_typ,gxfac,gxfac_sij,ia3,&
&             idir,indlmn_typ,kpgout_,matblk,ndgxdtfac,nincat,nkpgout_,nlmn,&
&             nloalg,npwout,nspinor,paw_opt,ph3dout,svectout,ucvol,vectout)
           else
             call opernlb_ylm(choice,cplex,cplex_fac,dgxdtfac,dgxdtfac_sij,&
&             dimffnlout,ffnlout_typ,gxfac,gxfac_sij,ia3,&
&             idir,indlmn_typ,kpgout,matblk,ndgxdtfac,nincat,nkpgout,nlmn,&
&             nloalg,npwout,nspinor,paw_opt,ph3dout,svectout,ucvol,vectout)
           end if
         end if

       end if ! choice==0

!      Deallocate temporary projected scalars
       ABI_DEALLOCATE(gx)
       ABI_DEALLOCATE(gxfac)
       ABI_DEALLOCATE(dgxdt)
       ABI_DEALLOCATE(dgxdtfac)
       ABI_DEALLOCATE(d2gxdt)
       if (paw_opt>=3)  then
         ABI_DEALLOCATE(dgxdtfac_sij)
         ABI_DEALLOCATE(gxfac_sij)
       end if

!      End sum on atom subset loop
       iatm=iatm+nincat;ia5=ia5+nincat
     end do
     if (allocated(indlmn_typ))  then
       ABI_DEALLOCATE(indlmn_typ)
     end if
     if (allocated(ffnlin_typ))  then
       ABI_DEALLOCATE(ffnlin_typ)
     end if
     if (signs==2)  then
       ABI_DEALLOCATE(ffnlout_typ)
     end if
     if (paw_opt>=2)  then
       ABI_DEALLOCATE(sij_typ)
     end if

!    End condition of existence of a non-local part
   else
     if (cpopt==0.or.cpopt==1) then
       do ispinor=1,nspinor
         do ia=1,nattyp(itypat)
           cprjin(iatm+ia,ispinor)%cp(:,1:nlmn)=zero
         end do
       end do
     end if
     if ((cpopt==1.or.cpopt==3).and.ndgxdt>0) then
       do ispinor=1,nspinor
         do ia=1,nattyp(itypat)
           cprjin(iatm+ia,ispinor)%dcp(:,1:ndgxdt,1:nlmn)=zero
         end do
       end do
     end if
     iatm=iatm+nattyp(itypat)
   end if

!  End atom type loop
   ia1=ia2+1
 end do

!Reduction in case of parallelization over spinors
 if (signs==1.and.mpi_enreg%paral_spin==1) then
   if (nnlout/=0) then
     call xsum_mpi(enlout,mpi_enreg%comm_spin,ierr)
   end if
   if (choice==3.or.choice==6.or.choice==23) then
     call xsum_mpi( enlk,mpi_enreg%comm_spin,ierr)
   end if
   if (choice==6) then
     call xsum_mpi(fnlk,  mpi_enreg%comm_spin,ierr)
     call xsum_mpi(strnlk,mpi_enreg%comm_spin,ierr)
   end if
 end if

!Convert stress tensor from reduced to cartesian (reciprocal space) coordinates
 if ((choice==3.or.choice==23).and.signs==1.and.(paw_opt==0.or.paw_opt==1.or.paw_opt==2)) then
   call strconv(enlout(1:6),gprimd,enlout(1:6))
   enlout(1:3)=(enlout(1:3)-enlk)/ucvol
   enlout(4:6)= enlout(4:6)/ucvol
 end if

!Convert elastic tensor from reduced to cartesian coordinates
 if (choice==6.and.signs==1.and.(paw_opt==0.or.paw_opt==1.or.paw_opt==2)) then
   ABI_ALLOCATE(enlout_tmp,(6+3*natom,6))
   ABI_ALLOCATE(work,(6))
   ABI_ALLOCATE(work_out,(6))
   enlout_tmp(:,:)=reshape(enlout(:),(/6+3*natom,6/))
   do mu=1,6
     call strconv(enlout_tmp(1:6,mu),gprimd,enlout_tmp(1:6,mu))
   end do
   do mu=1,6+3*natom
     work(1:6)=enlout_tmp(mu,1:6)
     work_out = zero
     call strconv(work,gprimd,work_out)
     work = work_out
     enlout_tmp(mu,1:6)=work(1:6)
   end do
   strnlk_out = zero
   call strconv(strnlk,gprimd,strnlk_out)
   strnlk = strnlk_out
   enlout(:)=reshape(enlout_tmp(:,:),(/6*(6+3*natom)/))
   ABI_DEALLOCATE(enlout_tmp)
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(work_out)
   do mub=1,6
     nub1=alpha(mub);nub2=beta(mub)
     do mua=1,6
       mu=mua+(3*natom+6)*(mub-1)
       nua1=alpha(mua);nua2=beta(mua)
       if (mua<=3.and.mub<=3) enlout(mu)=enlout(mu)+enlk
       if (mua<=3) enlout(mu)=enlout(mu)-strnlk(mub)
       if (mub<=3) enlout(mu)=enlout(mu)-strnlk(mua)
       if (nub1==nua2) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua1,nub2))
       if (nub2==nua2) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua1,nub1))
       if (nub1==nua1) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua2,nub2))
       if (nub2==nua1) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua2,nub1))
     end do
     if (mub<=3) then
       do ia1=1,natom
         ia2=3*(ia1-1);mu=ia2+6+(3*natom+6)*(mub-1)
         enlout(mu+1:mu+3)=enlout(mu+1:mu+3)-fnlk(ia2+1:ia2+3)
       end do
     end if
   end do
 end if

 if (signs==1)  then
   ABI_DEALLOCATE(fnlk)
   ABI_DEALLOCATE(strnlk)
   ABI_DEALLOCATE(strnlk_out)
 end if
 if (nkpgin_>0)  then
   ABI_DEALLOCATE(kpgin_)
 end if
 if (nkpgout_>0)  then
   ABI_DEALLOCATE(kpgout_)
 end if

!DEBUG
!write(std_out,*)' nonlop_ylm : exit '
!ENDDEBUG

end subroutine nonlop_ylm
!!***
