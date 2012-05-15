!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspatm
!! NAME
!! pspatm
!!
!! FUNCTION
!! Open atomic pseudopotential data file for a given atom,
!! read the three first lines, make some checks, then
!! call appropriate subroutine for the reading of
!! V(r) and wf R(r) data for each angular momentum, and subsequent
!! Fourier and Bessel function transforms for local and nonlocal potentials.
!! Close psp file at end.
!!
!! Handles pseudopotential files produced by (pspcod=1 or 4) Teter code,
!! or from the Goedecker-Teter-Hutter paper (pspcod=2),
!! or from the Hartwigsen-Goedecker-Hutter paper (pspcod=3 or 10)
!! or "Phoney pseudopotentials" (Hamman grid in real space) (pspcod=5)
!! or "Troullier-Martins pseudopotentials" from the FHI (pspcod=6)
!! or "XML format" (pspcod=9)
!! or "UPF PWSCF format" (pspcod=11)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, FrD, AF, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dq= spacing of the q-grid
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | ixc=exchange-correlation choice from main routine data file
!!   | pawxcdev=choice of XC development in PAW formalism
!!   | usexcnhat=choice for use of nhat in Vxc in PAW formalism
!!   | xclevel= XC functional level
!!  ipsp=id in the array of the currently read pseudo.
!!
!! OUTPUT
!!  ekb(dimekb)=
!!    ->NORM-CONSERVING PSPS ONLY (pspcod/=7):
!!      (Real) Kleinman-Bylander energies (hartree)
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for number of basis functions (l,n) (dimekb=lnmax)
!!             If any, spin-orbit components begin at l=mpsang+1
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!  xccc1d(n1xccc*(1-usepaw),6)=1D core charge function and five derivatives,
!!                              from psp file (used in NC only)
!!
!! SIDE EFFECTS
!! Input/Output :
!!  psps <type(pseudopotential_type)>=at output, values depending on the read
!!                                    pseudo are set.
!!   | dimekb(IN)=dimension of ekb (see module defs_psp.f)
!!   | filpsp(IN)=name of formatted external file containing atomic psp data.
!!   | lmnmax(IN)=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!   |           =if useylm=0, max number of (l,n)   comp. over all type of psps
!!   | lnmax(IN)=max. number of (l,n) components over all type of psps
!!   |           angular momentum of nonlocal pseudopotential
!!   | mpsang(IN)= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | mpssoang(IN)= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!   | mqgrid_ff(IN)=dimension of q (or G) grid for nl form factors (array ffspl)
!!   | mqgrid_vl(IN)=dimension of q (or G) grid for Vloc (array vlspl)
!!   | n1xccc(IN)=dimension of xccc1d ; 0 if no XC core correction is used
!!   | optnlxccc(IN)=option for nl XC core correction
!!   | positron(IN)=0 if electron GS calculation
!!   |              1 if positron GS calculation
!!   |              2 if electron GS calculation in presence of the positron
!!   | pspso(INOUT)=spin-orbit characteristics, govern the content of ffspl and ekb
!!   |          if =0 : this input requires NO spin-orbit characteristics of the psp
!!   |          if =2 : this input requires HGH characteristics of the psp
!!   |          if =3 : this input requires HFN characteristics of the psp
!!   |          if =1 : this input will be changed at output to 1, 2, 3, according
!!   |                  to the intrinsic characteristics of the psp file
!!   | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!   | qgrid_vl(mqgrid_vl)(IN)=values of q on grid from 0 to qmax (bohr^-1) for Vloc
!!   | usepaw(IN)= 0 for non paw calculation; =1 for paw calculation
!!   | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!   |            1=using Ylm, 0=using Legendre polynomials
!!   | vlspl_recipSpace(IN)=.true. if pseudo are expressed in reciprocal space.
!!   | znuclpsp(IN)=atomic number of atom as specified in input file to main routine
!!
!! NOTES
!!  Format expected for the three first lines of pseudopotentials
!!  (1) title (character) line
!!  (2) znucl,zion,pspdat
!!  (3) pspcod,pspxc,lmax,lloc,mmax,r2well
!!
!!  Dimensions of form factors and Vloc q grids must be the same in Norm-Conserving case
!!
!! PARENTS
!!      pspini
!!
!! CHILDREN
!!      close_xml_t,leave_new,open_xml_file,parse,psp10in,psp1in,psp2in,psp3in
!!      psp5in,psp6in,psp7in,psp8in,psp9in,psxml2ab,upf2abinit,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pspatm(dq,dtset,ekb,epsatm,ffspl,indlmn,ipsp,pawrad,pawtab,&
                & psps,vlspl,dvlspl,xcccrc,xccc1d)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xml_pseudo_types
 use m_xml_pseudo
#if defined HAVE_TRIO_FOX
 use fox_sax
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspatm'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_45_psp_parser
 use interfaces_65_psp, except_this_one => pspatm
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipsp
 real(dp),intent(in) :: dq
 real(dp),intent(out) :: epsatm,xcccrc
 type(dataset_type),intent(in) :: dtset
 type(pawrad_type),intent(out) :: pawrad
 type(pawtab_type),intent(out) :: pawtab
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: indlmn(6,psps%lmnmax)
 real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2)
 real(dp),intent(out) :: ekb(psps%dimekb*(1-psps%usepaw))
 real(dp),intent(out) :: ffspl(psps%mqgrid_ff,2,psps%lnmax)
 real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)
 real(dp),intent(out) :: xccc1d(psps%n1xccc*(1-psps%usepaw),6)

!Local variables ---------------------------------------
!scalars
 integer :: ii,il,ilmn,iln,iln0,lloc,lmax,mmax,pspcod,pspdat,pspxc,usexml
 integer :: useupf
 real(dp) :: qchrg,r2well,zion,znucl
 character(len=500) :: message
 character(len=70) :: testxml
 character(len=fnlen) :: title
!arrays
 integer,allocatable :: nproj(:)
 real(dp),allocatable :: e990(:),e999(:),ekb1(:),ekb2(:),epspsp(:),rcpsp(:)
 real(dp),allocatable :: rms(:)
!no_abirules
#if defined HAVE_TRIO_FOX
!!  usexml= 0 for non xml ps format ; =1 for xml ps format
 integer :: iostat,maxn_pots
 type(xml_t)                :: fxml
 type(pseudo_t), pointer    :: psxml
#endif

! ******************************************************************************

!DEBUG
!write(std_out,*)' pspatm : enter '
!if(.true.)stop
!ENDDEBUG

!Dimensions of form factors and Vloc q grids must be the same in Norm-Conserving case
 if ((psps%usepaw==0).and.(psps%mqgrid_ff/=psps%mqgrid_vl)) then
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' pspatm: BUG -',ch10,&
&   '  Dimension of q-grid for nl form factors (mqgrid_ff)',ch10,&
&   '  is different from dimension of q-grid for Vloc (mqgrid_vl) !',ch10,&
&   '  This is not allowed for norm-consevring psp.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

 write(message, '(a,t38,a)' ) &
& '- pspatm: opening atomic psp file',trim(psps%filpsp(ipsp))
 call wrtout(ab_out,  message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Check if the file pseudopotential file is written in XML
 usexml = 0
 open (unit=tmp_unit,file=psps%filpsp(ipsp),form='formatted',status='old')
 rewind (unit=tmp_unit)
 read(tmp_unit,*) testxml
 if(testxml(1:5)=='<?xml')then
   usexml = 1
 else
   usexml = 0
 end if
 close (unit=tmp_unit)

!Check if pseudopotential file is a Q-espresso UPF file
 useupf = 0
 open (unit=tmp_unit,file=psps%filpsp(ipsp),form='formatted',status='old')
 rewind (unit=tmp_unit)
 read(tmp_unit,*) testxml ! just a string, no relation to xml.
 if(testxml(1:9)=='<PP_INFO>')then
   useupf = 1
 else
   useupf = 0
 end if
 close (unit=tmp_unit)

!----------------------------------------------------------------------------
!allocate nproj here: can be read in now for UPF
 ABI_ALLOCATE(nproj,(psps%mpssoang))
 nproj(:)=0

 if (usexml /= 1 .and. useupf /= 1) then

!  Open the atomic data file, and read the three first lines
!  These three first lines have a similar format in all allowed psp files

!  Open atomic data file (note: formatted input file)
   open (unit=tmp_unit,file=psps%filpsp(ipsp),form='formatted',status='old')
   rewind (unit=tmp_unit)

!  Read and write some description of file from first line (character data)
   read (tmp_unit,'(a)') title
   write(message, '(a,a)' ) ' ',trim(title)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

!  Read and write more data describing psp parameters
   read (tmp_unit,*) znucl,zion,pspdat
   write(message, '(a,f9.5,f10.5,2x,i8,t47,a)' ) &
&   '-',znucl,zion,pspdat,'znucl, zion, pspdat'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*) pspcod,pspxc,lmax,lloc,mmax,r2well
   write(message, '(4i5,i10,f10.5,t47,a)' ) &
&   pspcod,pspxc,lmax,lloc,mmax,r2well,&
&   'pspcod,pspxc,lmax,lloc,mmax,r2well'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 else if (usexml == 1) then

#if defined HAVE_TRIO_FOX
   write(message,'(a,a)') &
&   '- pspatm: Reading pseudopotential header in XML form from ', trim(psps%filpsp(ipsp))
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   call open_xml_file(fxml,psps%filpsp(ipsp),iostat)
   if (iostat /=0) stop "Cannot open file"

   call parse(fxml,pcdata_chunk,startElement_handler=begin_element,endElement_handler=end_element)

   psxml => pseudo

   lloc   = 0
   r2well = 0

   call psxml2ab( psxml, znucl, zion, pspcod, pspxc, lmax, 1 )

   call close_xml_t(fxml)

   maxn_pots=size(psxml%pswf)

   do il = 1, maxn_pots
     if(associated( psxml%pswf(il)%V%data))   ABI_DEALLOCATE (psxml%pswf(il)%V%data)
     if(associated( psxml%pot(il)%V%data))    ABI_DEALLOCATE (psxml%pot(il)%V%data)
   end do   
   if(associated( psxml%core_charge%data))    ABI_DEALLOCATE (psxml%core_charge%data)
   if(associated (psxml%valence_charge%data)) ABI_DEALLOCATE (psxml%valence_charge%data)
#endif

 else if (useupf == 1) then
   if (psps%usepaw /= 0) then
     stop 'pspatm : Error: UPF format not allowed with PAW (USPP part not read yet)'
   end if

   pspcod = 11
   r2well = 0

!  should initialize znucl,zion,pspxc,lmax,lloc,mmax
   call upf2abinit (psps%filpsp(ipsp), znucl, zion, pspxc, lmax, lloc, mmax, &
&   psps, epsatm, xcccrc, indlmn, ekb, ffspl, nproj, &
&   vlspl, xccc1d)

 end if

!------------------------------------------------------------------------------
!Check data for consistency against main routine input

!Does required spin-orbit characteristics agree with format
!(At present, only HGH and phoney pseudopotentials can have spin-orbit)
!write(std_out,*) pspso
 if((pspcod/=3).and.(pspcod/=5).and.(pspcod/=10))then
!  If pspso requires internal characteristics, set it to 1 for non-HGH psps
   if(psps%pspso(ipsp)==1) psps%pspso(ipsp)=0
   if(psps%pspso(ipsp)/=0)then
     write(message, '(a,a,a,a,a,a,i3,a,a,a)' ) ch10,&
&     ' pspatm: ERROR -',ch10,&
&     '  Pseudopotential file cannot give spin-orbit characteristics,',ch10,&
&     '  while pspso(itypat)=',psps%pspso(ipsp),'.',ch10,&
&     '  Action : check your pseudopotential and input files for consistency.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call leave_new('COLL')
   end if
 end if

!Does nuclear charge znuclpsp agree with psp input znucl
 if (abs(psps%znuclpsp(ipsp)-znucl)>tol8) then
   write(message, '(4a,f10.5,2a,f10.5,5a)' ) ch10,&
&   ' pspatm: BUG -',ch10,&
&   '  Pseudopotential file znucl=',znucl,ch10,&
&   '  does not equal input znuclpsp=',psps%znuclpsp(ipsp),' better than 1e-08 .',ch10,&
&   '  znucl is read from the psp file in pspatm, while',ch10,&
&   '  znuclpsp is read in iofn2.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!Is the highest angular momentum within limits?
!Recall mpsang is 1+highest l for nonlocal correction.
!Nonlocal corrections for s, p, d, and f are supported.
 if (lmax+1>psps%mpsang) then
   write(message, '(a,a,a,a,i10,a,i10,a,a)' )ch10,&
&   ' pspatm: BUG -',ch10,&
&   '  input lmax+1=',lmax+1,' exceeds mpsang=',psps%mpsang,ch10,&
&   '  indicates input lmax too large for dimensions.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!Check several choices for ixc against pspxc
!ixc is from ABINIT code; pspxc is from atomic psp file
 if (dtset%ixc==0) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' pspatm: WARNING -',ch10,&
&   '  Note that input ixc=0 => no xc is being used.'
   call wrtout(std_out,  message,'COLL')
 else if(dtset%ixc/=pspxc) then
   write(message, '(a,a,a,a,i8,a,a,a,i8,a,a,a,a,a,a,a,a,a,a)' )ch10,&
&   ' pspatm: WARNING -',ch10,&
&   '  Pseudopotential file pspxc=',pspxc,',',ch10,&
&   '  not equal to input ixc=',dtset%ixc,'.',ch10,&
&   '  These parameters must agree to get the same xc ',ch10,&
&   '  in ABINIT code as in psp construction.',ch10,&
&   '  Action : check psp design or input file.',ch10,&
&   '  Assume experienced user. Execution will continue.',ch10
   call wrtout(std_out,  message,'COLL')
 end if

 if (lloc>lmax .and. pspcod/=4 .and. pspcod/=8 .and. pspcod/=10) then
   write(message, '(a,a,a,a,2i12,a,a,a,a)' )ch10,&
&   ' pspatm: ERROR -',ch10,&
&   '  lloc,lmax=',lloc,lmax,ch10,&
&   '  chosen l of local psp exceeds range from input data.',ch10,&
&   '  Action : check pseudopotential input file.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!Does the pspcod agree with type of calculation (paw or not)?
 if ((pspcod/=7.and.psps%usepaw==1).or.(pspcod==7.and.psps%usepaw==0)) then
   write(message, '(a,a,a,a,i2,a,a,i1,a)' )ch10,&
&   ' pspatm: BUG -',ch10,&
&   '  In reading atomic psp file, finds pspcod=',pspcod,ch10,&
&   '  This is not an allowed value with usepaw=',psps%usepaw,'.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!Check that pspcod is consistent with psps%vlspl_recipSpace
 if (.not.psps%vlspl_recipSpace .and. &
& (pspcod /= 2 .and. pspcod /= 3 .and. pspcod /= 10)) then
!  The following "if" statement can substitute the one just before once libBigDFT
!  has been upgraded to include pspcod 10
!  if (.not.psps%vlspl_recipSpace .and. (pspcod /= 2 .and. pspcod /= 3 .and. pspcod /= 10)) then
   write(message, '(a,a,a,a,i2,a,a)' )ch10,&
&   ' pspatm: BUG -',ch10,&
&   '  In reading atomic psp file, finds pspcod=',pspcod,ch10,&
&   '  This is not an allowed value with real space computation.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!MJV 16/6/2009 added pspcod 11 for upf format
 if( pspcod<1 .or. pspcod>11 ) then
   write(message, '(a,a,a,a,i12,a,a,a,a,a)' )ch10,&
&   ' pspatm: ERROR -',ch10,&
&   '  In reading atomic psp file, finds pspcod=',pspcod,ch10,&
&   '  This is not an allowed value.  Allowed values are',&
&   ' 1 to 11 .',ch10,&
&   '  Action : check pseudopotential input file.'
   call wrtout(std_out,  message,'COLL')
   call leave_new('COLL')
 end if

!-----------------------------------------------------------------------
!Set various terms to 0 in case not defined below
 ABI_ALLOCATE(e990,(psps%mpssoang))
 ABI_ALLOCATE(e999,(psps%mpssoang))
 ABI_ALLOCATE(rcpsp,(psps%mpssoang))
 ABI_ALLOCATE(rms,(psps%mpssoang))
 ABI_ALLOCATE(epspsp,(psps%mpssoang))
 ABI_ALLOCATE(ekb1,(psps%mpssoang))
 ABI_ALLOCATE(ekb2,(psps%mpssoang))
 e990(:)=zero ;e999(:)=zero
 rcpsp(:)=zero;rms(:)=zero
 ekb1(:)=zero ;ekb2(:)=zero
 epspsp(:)=zero
 qchrg=0

!----------------------------------------------------------------------
 if(pspcod==1 .or. pspcod==4)then

!  Teter pseudopotential (pspcod=1 or 4)
   call psp1in(dq,ekb,ekb1,ekb2,epsatm,epspsp,&
&   e990,e999,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,&
&   mmax,psps%mpsang,psps%mqgrid_ff,nproj,psps%n1xccc,pspcod,qchrg,psps%qgrid_ff,&
&   rcpsp,rms,psps%useylm,vlspl,xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

!  DEBUG
!  write(std_out,*)' pspatm : vlspl(1,2)= ',vlspl(1,2)
!  write(std_out,*)' pspatm : vlspl(1,1)= ',vlspl(1,1)
!  ENDDEBUG

 else if (pspcod==2)then

!  GTH pseudopotential
   call psp2in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps,vlspl,dvlspl,zion)
   xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

 else if (pspcod==3)then

!  HGH pseudopotential
   call psp3in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps, psps%pspso(ipsp), &
&   vlspl,zion)
   xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

 else if (pspcod==5)then

!  Old phoney pseudopotentials
   call psp5in(ekb,ekb1,ekb2,epsatm,epspsp,&
&   e990,e999,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,&
&   mmax,psps%mpsang,psps%mpssoang,psps%mqgrid_ff,nproj,psps%n1xccc,psps%pspso(ipsp),qchrg,psps%qgrid_ff,&
&   rcpsp,rms,psps%useylm,vlspl,xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

 else if (pspcod==6)then

!  FHI pseudopotentials
   call psp6in(ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&   psps%mpsang,psps%mqgrid_ff,nproj,psps%n1xccc,psps%optnlxccc,psps%positron,qchrg,psps%qgrid_ff,psps%useylm,vlspl,&
&   xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

 else if (pspcod==7)then

!  PAW "pseudopotentials"
   pawtab%usexcnhat=dtset%usexcnhat
   call psp7in(epsatm,ffspl,indlmn,dtset%ixc,lmax,psps%lmnmax,psps%lnmax,mmax,psps%mqgrid_ff,psps%mqgrid_vl,&
&   pawrad,pawtab,dtset%pawxcdev,psps%pspso(ipsp),psps%qgrid_ff,psps%qgrid_vl,vlspl,xcccrc,&
&   dtset%xclevel,dtset%xc_denpos,zion,psps%znuclpsp(ipsp))

!  Oblige if the user requests a subset of projectors of the PAW dataset
!  in that case several variables and arrays have to be puned and redefined
!  if (ANY(dtset%paw_use_proj_subset/=-1)) then

!  end if

 else if (pspcod==8)then

!  DRH pseudopotentials
   call psp8in(ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&   psps%mpsang,psps%mqgrid_ff,nproj,psps%n1xccc,qchrg,psps%qgrid_ff,psps%useylm,vlspl,&
&   xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))

 else if (pspcod==9)then

#if defined HAVE_TRIO_FOX
   call psp9in(psps%filpsp(ipsp),ekb,epsatm,ffspl,indlmn,lloc,lmax,psps%lmnmax,psps%lnmax,mmax,&
&   psps%mpsang,psps%mqgrid_ff,nproj,psps%n1xccc,qchrg,psps%qgrid_ff,psps%useylm,vlspl,&
&   xcccrc,xccc1d,zion,psps%znuclpsp(ipsp))
#endif

 else if (pspcod==10)then

!  HGH pseudopotential, full h/k matrix read
   call psp10in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps, psps%pspso(ipsp), &
&   vlspl,zion)
   xccc1d(:,:)=0.0d0 ; qchrg=0.0d0 ; xcccrc=0.0d0

!  NB for pspcod 11 the reading has already been done above.
 end if

 close (unit=tmp_unit)

!----------------------------------------------------------------------
 if (pspcod==2 .or. pspcod==3 .or. pspcod==10)then
   write(message, '(a,a,a,a,a,a,a,a,a,a)' )ch10,&
&   ' pspatm : COMMENT -',ch10,&
&   '  the projectors are not normalized,',ch10,&
&   '  so that the KB energies are not consistent with ',ch10,&
&   '  definition in PRB44, 8503 (1991). ',ch10,&
&   '  However, this does not influence the results obtained hereafter.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 if (pspcod/=7) then
   write(message, '(a,f14.8,a,a)' ) ' pspatm: epsatm=',epsatm,ch10,&
&   '         --- l  ekb(1:nproj) -->'

   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   iln0=0
   do ilmn=1,psps%lmnmax
     iln=indlmn(5,ilmn)
     if (iln>iln0) then
       il=indlmn(1,ilmn)
       if (indlmn(6,ilmn)==1) then
         iln0=iln0+nproj(il+1)
         write(message, '(13x,i1,4f12.6)' ) il,&
&         (ekb(iln+ii),ii=0,nproj(il+1)-1)
       else
         iln0=iln0+nproj(il+psps%mpsang)
         write(message, '(2x,a,i1,4f12.6)' ) 'spin-orbit ',il,&
         (ekb(iln+ii),ii=0,nproj(il+psps%mpsang)-1)
       end if
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
     end if
   end do
 end if

 write(message,'(3a)') ' pspatm: atomic psp has been read ',&
& ' and splines computed',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 ABI_DEALLOCATE(e990)
 ABI_DEALLOCATE(e999)
 ABI_DEALLOCATE(rcpsp)
 ABI_DEALLOCATE(rms)
 ABI_DEALLOCATE(ekb1)
 ABI_DEALLOCATE(ekb2)
 ABI_DEALLOCATE(epspsp)
 ABI_DEALLOCATE(nproj)

 if (dtset%prtvol > 9 .and. psps%usepaw==0 .and. psps%lmnmax>3) then
   open (unit=400, file="pspdata.dat")
   write (400,*) '# Pseudopotential data in reciprocal space as used by ABINIT'
   write (400,'(a)', ADVANCE='NO') '# index       vlocal   '
   if (psps%lnmax > 0) &
&   write (400,'(a,I3)', ADVANCE='NO')   '           1st proj(l=', indlmn(1,1)
   if (psps%lnmax > 1) &
&   write (400,'(a,I3)', ADVANCE='NO')   ')            2nd(l=', indlmn(1,2)
   if (psps%lnmax > 2) &
&   write (400,'(a,I3,a)', ADVANCE='NO') ')            3rd(l=', indlmn(1,3), ')'
   write (400,*)

   do ii = 1, psps%mqgrid_vl
     write(400, '(I5,E24.16)', ADVANCE='NO') ii, vlspl(ii,1)
     if (psps%lnmax > 0) write(400, '(E24.16)', ADVANCE='NO') ffspl(ii,1,1)
     if (psps%lnmax > 1) write(400, '(E24.16)', ADVANCE='NO') ffspl(ii,1,2)
     if (psps%lnmax > 2) write(400, '(E24.16)', ADVANCE='NO') ffspl(ii,1,3)
     write(400, *)
   end do
   close(400)

   open (unit=900, file="nlcc_derivs.dat")
   write (900,*) '# Non-linear core corrections'
   write (900,*) '#  r, pseudocharge, 1st, 2nd, 3rd, 4th, 5th derivatives'
   do ii = 1, psps%n1xccc
     write (900,*) xcccrc*(ii-1)/(psps%n1xccc-1), xccc1d(ii,1), xccc1d(ii,2), &
&     xccc1d(ii,3), xccc1d(ii,4), &
&     xccc1d(ii,5), xccc1d(ii,6)
   end do
   close(900)
 end if

!DEBUG
!write(std_out,*)' pspatm : vlspl(1,2)= ',vlspl(1,2)
!write(std_out,*)' pspatm : vlspl(10,2)= ',vlspl(10,2)
!write(std_out,*)' pspatm : vlspl(100,2)= ',vlspl(100,2)
!write(std_out,*)' pspatm : vlspl(1,1)= ',vlspl(1,1)
!write(std_out,*)' pspatm : vlspl(10,1)= ',vlspl(10,1)
!write(std_out,*)' pspatm : vlspl(100,1)= ',vlspl(100,1)
!write(std_out,*)' pspatm : ffspl p',ffspl(10,1,2,1),ffspl(10,1,2,2)
!write(std_out,*)' pspatm : ffspl p_so',ffspl(10,1,4,1),ffspl(10,1,4,2)
!write(std_out,*)' pspatm : ffspl d',ffspl(10,1,3,1),ffspl(10,1,3,2)
!write(std_out,*)' pspatm : ffspl d_so',ffspl(10,1,5,1),ffspl(10,1,5,2)
!ENDDEBUG

!DEBUG
!write(std_out,*)' pspatm : exit '
!if(.true.)stop
!ENDDEBUG

end subroutine pspatm
!!***
