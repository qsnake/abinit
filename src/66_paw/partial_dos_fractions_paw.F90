!{\src2tex{textfont=tt}}
!!****f* ABINIT/partial_dos_fractions_paw
!! NAME
!! partial_dos_fractions_paw
!!
!! FUNCTION
!!  Calculate PAW contributions to the partial DOS fractions (tetrahedron method)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (SM,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset     structured datatype, from which one uses :
!!   iatsph(nasph)=number of atoms used to project dos
!!   indlmn(6,lmnmax,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!   lmnmax=max. number of (l,m,n) numbers over all types of atom
!!   kpt(3,nkpt)  =irreducible kpoints
!!   mband        =max number of bands per k-point
!!   mkmem        =number of kpoints in memory
!!   natom        =number of atoms in total
!!   natsph       =number of atoms ofor which the spherical decomposition must be done
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nspinor      =number of spinor components
!!   nsppol       =1 or 2 spin polarization channels
!!  fatbands_flag =1 if pawfatbnd=1 or 2
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  mpi_enreg=informations about MPI parallelization
!!  m_dos_flag=option for the m-contributions to the partial DOS
!!  ndosfraction=natsph*mbesslang
!!  paw_dos_flag=option for the PAW contributions to the partial DOS
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  === If paw_dos_flag==1:
!!   dos_fractions_paw1(ikpt,iband,isppol,natom*mbesslang) = contribution to
!!       dos fractions from the PAW partial waves (phi)
!!   dos_fractions_pawt1(ikpt,iband,isppol,natom*mbesslang) = contribution to
!!       dos fractions from the PAW pseudo partial waves (phi_tild)
!!
!! SIDE EFFECTS
!!  dos_fractions(ikpt,iband,isppol,ndosfraction) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol
!!    As input: contains only the pseudo contribution
!!    As output: contains pseudo contribution + PAW corrections
!!  == if m_dos_flag==1
!!  dos_fractions_m(ikpt,iband,isppol,ndosfraction*mbesslang*m_dos_flag) =
!!              m discretization of partial DOS fractions
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      cprj_alloc,cprj_diskinit_r,cprj_free,cprj_get,leave_new,leave_test
!!      simp_gen,timab,wrtout,xcomm_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine partial_dos_fractions_paw(atindx1,cprj,dimcprj,dos_fractions,dos_fractions_m,&
& dos_fractions_paw1,dos_fractions_pawt1,&
& dtfil,dtset,fatbands_flag,indlmn,lmnmax,mbesslang,mcprj,mkmem,&
& mpi_enreg,m_dos_flag,ndosfraction,paw_dos_flag,pawrad,pawtab)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_errors

 use m_radmesh,   only : simp_gen

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'partial_dos_fractions_paw'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_44_abitypes_defs
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fatbands_flag,lmnmax,mbesslang,mcprj,mkmem,m_dos_flag,ndosfraction
 integer,intent(in) :: paw_dos_flag
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: atindx1(dtset%natom),dimcprj(dtset%natom),indlmn(6,lmnmax,dtset%ntypat)
 real(dp),intent(inout) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
 real(dp),intent(inout) :: dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang*max(m_dos_flag,fatbands_flag))
 real(dp),intent(out) :: dos_fractions_paw1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*paw_dos_flag)
 real(dp),intent(out) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*paw_dos_flag)
 type(cprj_type) :: cprj(dtset%natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: basis_size,cplex,iat,iatom,iband,ibg,ibsp,ierr,ikpt,il,ilang,ilmn,iln,im,iorder_cprj,ispinor,isppol
 integer :: itypat,j0lmn,j0ln,jl,jlmn,jln,jm,klmn,kln,lmn_size,me,my_nspinor,nband_k,spaceComm
 real(dp) :: cpij
 character(len=500) :: message
!arrays
 integer ,allocatable :: dimcprj_atsph(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: int1(:,:),int2(:,:),int1m2(:,:)
 type(cprj_type),allocatable :: cprj_k(:,:),cprj_tmp(:,:)

!******************************************************************************************

 DBG_ENTER("COLL")

!m-decompozed DOS not compatible with PAW-decompozed DOS
 if(m_dos_flag==1.and.paw_dos_flag==1) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' partial_dos_fractions_paw : ERROR -',ch10,&
&   '  m-decompozed DOS not compatible with PAW-decompozed DOS !'
   call wrtout(std_out,message,'COLL')
   call leave_new('PERS')
 end if

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spin)
 if (paw_dos_flag==1.or.fatbands_flag==1) then
   dos_fractions_paw1 =zero
   dos_fractions_pawt1=zero
 end if

!Prepare some useful integrals
 basis_size=pawtab(1)%basis_size
 if (dtset%ntypat>1) then
   do itypat=1,dtset%ntypat
     basis_size=max(basis_size,pawtab(itypat)%basis_size)
   end do
 end if
 ABI_ALLOCATE(int1  ,(basis_size*(basis_size+1)/2,dtset%natsph))
 ABI_ALLOCATE(int2,(basis_size*(basis_size+1)/2,dtset%natsph))
 ABI_ALLOCATE(int1m2,(basis_size*(basis_size+1)/2,dtset%natsph))
 int1=zero;int2=zero;int1m2=zero
 do iat=1,dtset%natsph
   iatom=dtset%iatsph(iat)
   itypat= dtset%typat(iatom)
   do jln=1,pawtab(itypat)%basis_size
     j0ln=jln*(jln-1)/2
     do iln=1,jln
       kln=j0ln+iln
       call simp_gen(int1(kln,iat),pawtab(itypat)%phiphj(:,kln),pawrad(itypat))
       if (dtset%pawprtdos<2) then
         call simp_gen(int2(kln,iat),pawtab(itypat)%tphitphj(:,kln),pawrad(itypat))
         int1m2(kln,iat)=int1(kln,iat)-int2(kln,iat)
       else
         int2(kln,iat)=zero;int1m2(kln,iat)=int1(kln,iat)
       end if
     end do !iln
   end do !jln
 end do

!Antiferro case
 if (dtset%nspden==2.and.dtset%nsppol==1.and.dtset%nspinor==1) then
   int1m2(:,:)=half*int1m2(:,:)
   if (paw_dos_flag==1.or.fatbands_flag==1) then
     int1(:,:)=half*int1(:,:);int2(:,:)=half*int2(:,:)
   end if
 end if

!Init parallelism
 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%comm_kpt)
 call xme_init(mpi_enreg,me)

!Prepare temporary PAW file if mkmem==0
 iorder_cprj=0
 call cprj_diskinit_r(atindx1,dtset%natom,iorder_cprj,dtset%mkmem,dtset%natom,0,dimcprj,my_nspinor,dtfil%unpaw)

!LOOPS OVER SPINS,KPTS
 ibg=0
 do isppol =1,dtset%nsppol
   do ikpt=1,dtset%nkpt

     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     if(mpi_enreg%paral_compil_kpt==1)then
       if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) cycle
     end if

     cplex=2;if (dtset%istwfk(ikpt)>1) cplex=1
     ABI_ALLOCATE(cprj_k,(dtset%natsph,my_nspinor*nband_k))
     ABI_ALLOCATE(dimcprj_atsph,(dtset%natsph))
     do iat=1,dtset%natsph
       dimcprj_atsph(iat)=dimcprj(dtset%iatsph(iat))
     end do
     call cprj_alloc(cprj_k,0,dimcprj_atsph)
     ABI_DEALLOCATE(dimcprj_atsph)

!    Extract cprj for this k-point according to mkmem
     if (mkmem==0) then
       ABI_ALLOCATE(cprj_tmp,(dtset%natom,nband_k*my_nspinor))
       call cprj_alloc(cprj_tmp,0,dimcprj)
       call cprj_get(atindx1,cprj_tmp,cprj,dtset%natom,1,ibg,ikpt,iorder_cprj,isppol,dtset%mband,dtset%mkmem,&
&       mpi_enreg,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,dtfil%unpaw)
       ibsp=0
       do iband=1,nband_k
         do ispinor=1,my_nspinor
           ibsp=ibsp+1
           do iat=1,dtset%natsph
             iatom=dtset%iatsph(iat)
             cprj_k(iat,ibsp)%cp(:,:)=cprj_tmp(iatom,ibsp)%cp(:,:)
           end do
         end do
       end do
       call cprj_free(cprj_tmp)
       ABI_DEALLOCATE(cprj_tmp)
     else ! mkmem/=0
       ibsp=0
       do iband=1,nband_k
         do ispinor=1,my_nspinor
           ibsp=ibsp+1
           do iat=1,dtset%natsph
             iatom=dtset%iatsph(iat)
             cprj_k(iat,ibsp)%cp(:,:)=cprj(iatom,ibsp+ibg)%cp(:,:)
           end do
         end do
       end do
     end if

!    LOOP OVER ATOMS
     do iat=1,dtset%natsph
       iatom=dtset%iatsph(iat)
       itypat= dtset%typat(iatom)
       lmn_size=pawtab(itypat)%lmn_size

!      LOOP OVER BANDS
       do iband=1,nband_k
         if(mpi_enreg%paral_compil_kpt==1)then
           if(abs(mpi_enreg%proc_distrb(ikpt,iband,isppol)-me)/=0) cycle
         end if
         ibsp=(iband-1)*my_nspinor
         do ispinor=1,my_nspinor
           ibsp=ibsp+1

           do ilang=1,mbesslang

             do jlmn=1,lmn_size
               jl=indlmn(1,jlmn,itypat)
               jm=indlmn(2,jlmn,itypat)
               j0lmn=jlmn*(jlmn-1)/2
               do ilmn=1,jlmn
                 il=indlmn(1,ilmn,itypat)
                 im=indlmn(2,ilmn,itypat)
                 klmn=j0lmn+ilmn
                 kln=pawtab(itypat)%indklmn(2,klmn)

                 if (il==ilang-1.and.jl==ilang-1.and.im==jm) then

                   cpij=cprj_k(iat,ibsp)%cp(1,ilmn)*cprj_k(iat,ibsp)%cp(1,jlmn)
                   if (cplex==2) cpij=cpij+cprj_k(iat,ibsp)%cp(2,ilmn)*cprj_k(iat,ibsp)%cp(2,jlmn)
                   cpij=pawtab(itypat)%dltij(klmn)*cpij

                   dos_fractions(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                   dos_fractions(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                   cpij*int1m2(kln,iat)
                   if (m_dos_flag==1) then
                     dos_fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im)= &
&                     dos_fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&                     cpij*int1m2(kln,iat)
                   end if
                   if (fatbands_flag==1) then
                     dos_fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im)= &
&                     dos_fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&                     cpij*int1(kln,iat)
                   end if
                   if (paw_dos_flag==1) then
                     dos_fractions_paw1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                     dos_fractions_paw1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                     cpij*int1(kln,iat)
                     dos_fractions_pawt1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                     dos_fractions_pawt1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                     cpij*int2(kln,iat)
                   end if

                 end if

               end do !ilmn
             end do   !jlmn

           end do ! ilang
         end do ! ispinor
       end do ! iband

     end do !iatom

     if (mkmem/=0) ibg = ibg + my_nspinor*nband_k
     call cprj_free(cprj_k)
     ABI_DEALLOCATE(cprj_k)
   end do ! ikpt
 end do ! isppol

 ABI_DEALLOCATE(int1)
 ABI_DEALLOCATE(int2)
 ABI_DEALLOCATE(int1m2)

!Reduce data in case of parallelism
 if(mpi_enreg%paral_compil_kpt==1)then
   call leave_test()
   call timab(48,1,tsec)
   call xsum_mpi(dos_fractions,spaceComm,ierr)
   if (m_dos_flag==1.or.fatbands_flag==1) then
     call xsum_mpi(dos_fractions_m,spaceComm,ierr)
   end if
   if (paw_dos_flag==1) then
     call xsum_mpi(dos_fractions_paw1,spaceComm,ierr)
     call xsum_mpi(dos_fractions_pawt1,spaceComm,ierr)
   end if
   call timab(48,2,tsec)
 end if
 if (mpi_enreg%paral_spin==1) then
   call xsum_mpi(dos_fractions,mpi_enreg%comm_spin,ierr)
   if (m_dos_flag==1.or.fatbands_flag==1) then
     call xsum_mpi(dos_fractions_m,mpi_enreg%comm_spin,ierr)
   end if
   if (paw_dos_flag==1) then
     call xsum_mpi(dos_fractions_paw1, mpi_enreg%comm_spin,ierr)
     call xsum_mpi(dos_fractions_pawt1, mpi_enreg%comm_spin,ierr)
   end if
 end if

!Averaging: A quick hack for m-decomposed LDOS:
!BA: not valid in presence of spin-orbit coupling  !
 if (m_dos_flag==1.and.fatbands_flag==0) then
!  if pawfatbnd is activated, one think in the cubic harmonics basis
!  whereas m_dos_flag is in the spherical harmonics basis.
!  the following trick is done in order to have everything
!  in the complex spherical basis (not useful for pawfatbnd if we want to
!  have e.g t2g and eg d-orbitals).
   do iat=1,dtset%natsph
     do il = 0, mbesslang-1
       do im = 1, il
         dos_fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im) = &
&         (dos_fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&         dos_fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1-im))/2
         dos_fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1-im) = &
&         dos_fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im)
       end do
     end do
   end do !iatom
 end if

 DBG_EXIT("COLL")

end subroutine partial_dos_fractions_paw
!!***
