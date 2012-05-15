!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnup1
!! NAME
!! clnup1
!!
!! FUNCTION
!! Perform "cleanup" at end of execution of gstate routine.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  acell(3)=length scales of primitive translations (bohr)
!!  dosdeltae=DOS delta of Energy
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree) for all bands
!!                           at each k point
!!  enunit=choice for units of output eigenvalues: 0=>hartree,
!!   1=> eV, 2=> hartree and eV
!!  fermie=fermi energy (Hartree)
!!  fnameabo_dos=filename of output DOS file
!!  fnameabo_eig=filename of output EIG file
!!  fred(3,natom)=d(E)/d(xred) (hartree)
!!  iatfix(3,natom)=0 if not fixed along specified direction,
!!                  1 if fixed
!!  iscf=parameter controlling scf or non-scf choice
!!  kptopt=option for the generation of k points
!!  kptns(3,nkpt)=k points in terms of recip primitive translations
!!  mband=maximum number of bands
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in unit cell
!!  nband(nkpt*nsppol)=number of bands
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!            see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nstep=desired number of electron iteration steps
!!  occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band
!!                                    and k point
!!  occopt=option for occupancies
!!  prtdos= if == 1, will print the density of states
!!  prtfor= if >0, will print the forces
!!  prtstm= input variable prtstm
!!  prtvol=control print volume and debugging
!!  resid(mband*nkpt*nsppol)=squared residuals for each band and
!!                           k point where
!!                     resid(n,k)=|<C(n,k)|(H-e(n,k))|C(n,k)>|^2
!!  rhor(nfft,nspden)=electron density (electrons/bohr^3)
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  tphysel="physical" electronic temperature with FD occupations
!!  tsmear=smearing energy or temperature (if metal)
!!  vxcavg=average of vxc potential
!!  wtk(nkpt)=real(dp) array of k-point weights
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  (only print and write to disk)
!!
!! NOTES
!!  MG: 40 input arguments for just 227 lines of code, it is a
!!  world record Well, considering only the executable lines we
!!  have 91 lines of code which means an average of 2.3 lines
!!  for each variable. Unbelievable
!!  And the funny thing is that this routine is supposed to do a
!!  cleanup!
!!  GAF: Done, 20 input arguments for 252 lines of code
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      getnel,metric,prteigrs,prtrhomxmn,prtxf,write_eig,wrtout,xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine clnup1(acell,dtset,eigen,fermie,&
  & fnameabo_dos,fnameabo_eig,fred,&
  & mpi_enreg,nfft,ngfft,occ,prtfor,&
  & resid,rhor,rprimd,vxcavg,xred)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnup1'
 use interfaces_14_hidewrite
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_61_ionetcdf
 use interfaces_62_occeig
 use interfaces_67_common, except_this_one => clnup1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 integer,intent(in) :: prtfor
 real(dp),intent(in) :: fermie
 real(dp),intent(in) :: vxcavg
 character(len=fnlen),intent(in) :: fnameabo_dos,fnameabo_eig
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in)  :: ngfft(18)
 real(dp),intent(in) :: acell(3)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: fred(3,dtset%natom)
 real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer :: comm,iatom,ii,iscf_dum,iwfrc,me,nnonsc,option,unitdos
 real(dp) :: entropy,grmax,grsum,maxocc,nelect,tolwf,ucvol
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 character(len=500) :: message
 character(len=fnlen) filename
!arrays
 real(dp),allocatable :: doccde(:)

! ****************************************************************

!DEBUG
!write(std_out,*)' clnup1 : enter,  rhor(1,1),rhor(1,2)',&
!& rhor(1,1),rhor(1,2)
!stop
!ENDDEBUG

 if(dtset%prtstm==0)then
!  Write reduced coordinates xred
   write(message, '(a,i5,a)' ) &
&   ' reduced coordinates (array xred) for',dtset%natom,' atoms'
   call wrtout(ab_out,message,'COLL')
   do iatom=1,dtset%natom
     write(message, '(1x,3f20.12)' ) xred(:,iatom)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

!Write reduced gradients if iscf > 0 and dtset%nstep>0 and
!dtset%prtstm==0
 if (dtset%iscf>0.and.dtset%nstep>0.and.dtset%prtstm==0) then

!  Compute absolute maximum and root mean square value of gradients
   grmax=0.0_dp
   grsum=0.0_dp
   do iatom=1,dtset%natom
     do ii=1,3
!      To be activated in v5.5
!      grmax=max(grmax,abs(fred(ii,iatom)))
       grmax=max(grmax,fred(ii,iatom))
       grsum=grsum+fred(ii,iatom)**2
     end do
   end do
   grsum=sqrt(grsum/dble(3*dtset%natom))

   write(message, '(1x,a,1p,e12.4,a,e12.4,a)' ) &
&   'rms dE/dt=',grsum,'; max dE/dt=',grmax,&
&   '; dE/dt below (all hartree)'
   call wrtout(ab_out,message,'COLL')
   do iatom=1,dtset%natom
     write(message, '(i5,1x,3f20.12)' ) iatom,fred(1:3,iatom)
     call wrtout(ab_out,message,'COLL')
   end do

 end if

 if(dtset%prtstm==0)then

!  Compute and write out dimensional cartesian coords and forces:
   write(message,*)' '
   call wrtout(ab_out,message,'COLL')

!  (only write forces if iscf > 0 and dtset%nstep>0)
   if (dtset%iscf<=0.or.dtset%nstep<=0.or.prtfor==0) then
     iwfrc=0
   else
     iwfrc=1
   end if

   call prtxf(fred,dtset%iatfix,ab_out,iwfrc,dtset%natom,&
&   rprimd,xred)

!  Write length scales
   write(message, '(1x,a,3f16.12,a)' )&
&   'length scales=',acell,' bohr'
   call wrtout(ab_out,message,'COLL')
   write(message, '(14x,a,3f16.12,a)' )&
&   '=',Bohr_Ang*acell(1:3),' angstroms'
   call wrtout(ab_out,message,'COLL')

 end if

 option=1 ; nnonsc=0 ; tolwf=0.0_dp
!write(std_out,*) 'WRITE_EIG option',option

 if(dtset%iscf<=0 .and. dtset%iscf/=-3)option=3
 iscf_dum=dtset%iscf
 if(dtset%nstep==0)iscf_dum=0
 if(dtset%tfkinfunc==0)then
   call prteigrs(eigen,dtset%enunit,fermie,fnameabo_eig,ab_out,&
&   iscf_dum,dtset%kptns,dtset%kptopt,dtset%mband,&
&   dtset%nband,dtset%nkpt,nnonsc,dtset%nsppol,occ,&
&   dtset%occopt,option,dtset%prteig,dtset%prtvol,resid,tolwf,&
&   vxcavg,dtset%wtk)

#if defined HAVE_TRIO_NETCDF
   if (dtset%prteig==1) then
     call xcomm_world(mpi_enreg,comm,myrank=me)
     if (me==0) then
!      write(std_out,*) 'WRITE_EIG'
       filename=trim(fnameabo_eig)//'.nc'
       call write_eig(eigen,&
&       filename,&
&       dtset%kptns,&
&       dtset%mband,&
&       dtset%nband,&
&       dtset%nkpt,&
&       dtset%nsppol)
     end if
   end if
#endif

 end if

!Compute and print location of maximal and minimal density
!DEBUG
!write(std_out,*) ' clnup1 : before prtrhomxmn '
!write(std_out,*) ' dtset%nspden,  rhor(1,1), rhor(1,2) ',dtset%nspden,rhor(1,1),&
!&  rhor(1,2)
!ENDDEBUG
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call prtrhomxmn(ab_out,mpi_enreg,nfft,ngfft,dtset%nspden,2,rhor,ucvol=ucvol)

!If needed, print DOS
 if(dtset%prtdos==1)then
   unitdos=tmp_unit
   option=2
   open (unit=unitdos,file=fnameabo_dos,status='unknown',&
&   form='formatted')
   rewind(unitdos)
   maxocc=two/(dtset%nspinor*dtset%nsppol)  ! Will not work in the fixed
!  moment case
   ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
   call getnel(doccde,dtset%dosdeltae,eigen,entropy,fermie,&
&   maxocc,dtset%mband,dtset%nband,nelect,dtset%nkpt,&
&   dtset%nsppol,occ,dtset%occopt,option,dtset%tphysel,&
&   dtset%tsmear,unitdos,dtset%wtk)
   ABI_DEALLOCATE(doccde)
 end if

!DEBUG
!write(std_out,*)' clnup1 : exit,  fred=',fred(1,1)
!stop
!ENDDEBUG

end subroutine clnup1
!!***
