!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_io_gkk
!! NAME
!!  m_io_gkk
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2012 ABINIT group (MG,MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

MODULE m_io_gkk

 use m_profiling

 use defs_basis
 use defs_abitypes
 use m_errors

 use m_io_tools,    only : get_unit
 use m_header,      only : hdr_clean
 use m_crystal,     only : crystal_structure
 use m_crystal_io,  only : init_crystal_from_hdr

 implicit none

 private 
!!***

!!****t* m_io_gkk/gkkfd_t
!! NAME
!!  gkkfd_t        
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: gkkfd_t

   integer :: funt
   ! Fortran unit number.

   integer :: n1wf
   ! number of perturbations in the gkk file

   integer :: nqibz
   ! Number of q-points stored on disk (supposed to be in the IBZ)

   integer :: natom
   ! Number of atoms.

   integer :: nband_tot
   ! Total number of bands stored on disk.

   !%integer :: nkibz

   integer :: nsppol

   logical :: in_gsblock

   character(len=fnlen) :: fname 
   ! File name 

   integer,pointer :: pos_q(:)   SET2NULL
   ! pos_q(nqibz)
   ! The position inside the GKK file of the first perturbation associated to the q-point.
   ! pos_q is in units of blocks (Hdr,<u0|H1|u0>)

   integer,pointer :: pos_pertq(:,:)
   ! pos_pertq(3*natom,nqibz)
   ! The position inside the GKK file of the perturbation with q-point iq_ibz
   ! pos_pertq is in units of blocks (Hdr,<u0|H1|u0>)
   ! A negative value is used to signal that the perturbations is not available on file
   ! This usually happens if the perturbation can be reconstructed by symmetry.

   integer,pointer :: npert_q(:)        SET2NULL
   ! npert_q(nqibz)
   ! Number of perturbations calculated for each q-point.

   integer,pointer :: pertcase_q(:,:)   SET2NULL
   ! pertcase_q(MAXVAL(npert_q),nqibz)
   ! The list of perturbations for the different q-points.

   real(dp),pointer :: qibz(:,:)  SET2NULL
   ! qibz(3,nqibz)
   ! The list of q-points in the IBZ.

 end type gkkfd_t
!!***
 
 ! Methods:
 public :: gkkfd_nullify
 public :: gkkfd_free
 public :: gkkfd_init
 public :: gkkfd_read_h1me

CONTAINS  !===========================================================
!!***

!!****f* m_io_gkk/gkkfd_nullify
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gkkfd_nullify(Fd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkkfd_nullify'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkkfd_t),intent(inout) :: Fd

! ************************************************************************

 !@gkkfd_t

 ! integer 
 nullify(Fd%pos_q     )
 nullify(Fd%pos_pertq )
 nullify(Fd%npert_q   )
 nullify(Fd%pertcase_q)

 ! real
 nullify(Fd%qibz)

end subroutine gkkfd_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_io_gkk/gkkfd_free
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gkkfd_free(Fd)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkkfd_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gkkfd_t),intent(inout) :: Fd

! ************************************************************************

 !@gkkfd_t

 ! integer 
 if (associated(Fd%pos_q     ))   then
   ABI_DEALLOCATE(Fd%pos_q)
 end if
 if (associated(Fd%pos_pertq ))   then
   ABI_DEALLOCATE(Fd%pos_pertq)
 end if
 if (associated(Fd%npert_q   ))   then
   ABI_DEALLOCATE(Fd%npert_q)
 end if
 if (associated(Fd%pertcase_q))   then
   ABI_DEALLOCATE(Fd%pertcase_q)
 end if

 ! real
 if (associated(Fd%qibz))   then
   ABI_DEALLOCATE(Fd%qibz)
 end if

 close(Fd%funt)

end subroutine gkkfd_free
!!***

!----------------------------------------------------------------------

!!****f* m_io_gkk/gkkfd_init
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gkkfd_init(Fd,fname,GS_Hdr,funt)

 use defs_basis
 use m_header

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkkfd_init'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: funt
 character(len=*),intent(in) :: fname
 type(gkkfd_t),intent(inout) :: Fd
 type(hdr_type),intent(out) :: GS_Hdr

!Local variables-------------------------------
!scalars
 integer,parameter :: buf_sz=1000
 integer :: istat,rdwr,fform,i1wf
 integer :: nband,nkibz,nsppol,ikpt,spin,nqfound,iq_ibz
 integer :: ik_rf,next,npe
 character(len=500) :: msg
 type(hdr_type) :: PH_Hdr

!arrays
 integer :: buf_posq(buf_sz)
 integer,allocatable :: buf_npert(:),buf_pertcase(:,:),buf_pospertq(:,:)
 real(dp) :: buf_q(3,buf_sz)
 real(dp),allocatable :: eigenGS(:,:,:)

! ************************************************************************

 !@gkkfd_t
 call gkkfd_nullify(Fd)

 call hdr_nullify(PH_Hdr)

 !
 ! Open the file.
 if (.not.PRESENT(funt)) then
   Fd%funt  = get_unit()
   open(unit=Fd%funt,file=fname,form='unformatted',status='old',iostat=istat)
   ABI_CHECK(istat==0,"Opening "//TRIM(fname))
 else 
   Fd%funt  = funt
   rewind(Fd%funt)
 end if

 Fd%fname = fname
 Fd%in_gsblock = .TRUE.
 !
 ! Read in GS header of file without rewinding it
 rdwr = 5 
 call hdr_io(fform,GS_Hdr,rdwr,Fd%funt)
 ABI_CHECK(fform/=0,'GKK header was mis-read. fform == 0')

 if ( ANY(GS_Hdr%nband(:)/= GS_Hdr%nband(1)) ) then
     write (msg,'(3a)')&
&     ' Use the same number of bands for all kpts : ',ch10,&
&     ' could have spurious effects if efermi is too close to the last band '
     MSG_ERROR(msg)
 end if
 !
 ! Echo the header to screen
 rdwr = 4 
 call hdr_io(fform,GS_Hdr,rdwr,6)
 !
 ! Copy dimensions
 nband  = GS_Hdr%nband(1)
 nkibz  = GS_Hdr%nkpt
 nsppol = GS_Hdr%nsppol 
 Fd%natom  = GS_Hdr%natom

 Fd%nband_tot = nband
 !%Fd%nkibz     = nkibz
 Fd%nsppol    = nsppol

 !==================================================
 ! Read GS eigenvalues for each irreducible kpt and
 ! number of 1WF files contributing to the GKK file
 !==================================================
 ABI_ALLOCATE(eigenGS,(nband,nkibz,nsppol))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,'out-of-memory in eigenGS')

 do spin=1,nsppol
   do ikpt=1,nkibz
     read(Fd%funt) eigenGS(:,ikpt,spin)
   end do
 end do

 ABI_DEALLOCATE(eigenGS)
 Fd%in_gsblock=.FALSE.

 !read number of 1WF files contributing to the GKK file
 read(Fd%funt) Fd%n1wf
 write(msg,'(a,i0)')' elphon : number of perturbations in the gkk file = ',Fd%n1wf
 call wrtout(std_out,msg,'COLL')
 !
 ! Analyze the blocks finding the number of q-points stored on file.

 ABI_ALLOCATE(buf_npert,(buf_sz))
 buf_npert    = 0
 ABI_ALLOCATE(buf_pertcase,(3*Fd%natom,buf_sz))
 buf_pertcase = -HUGE(1)
 ABI_ALLOCATE(buf_pospertq,(3*Fd%natom,buf_sz))
 buf_pospertq = -HUGE(1)

 nqfound = 0
 do i1wf=1,Fd%n1wf
   write (msg,'(2a,i4,a,i4)')ch10,' gkkfd_init : reading 1WF header # ',i1wf,' /',Fd%n1wf
   call wrtout(std_out,msg,'COLL')
   !
   !  Could check for compatibility of natom, kpt grids, ecut, qpt with DDB grid...
   !  MG: Also this task should be done in mrggkk
   !
   rdwr = 5 !read without rewinding
   call hdr_io(fform,PH_Hdr,rdwr,Fd%funt)
   if (fform == 0) then
     write (msg,'(a,i0,a)')' 1WF header number ',i1wf,' was mis-read. fform == 0'
     MSG_ERROR(msg)
   end if

   write (msg,'(a,i0)')" Number of k-points for this perturbation: ",PH_Hdr%nkpt
   call wrtout(std_out,msg,'COLL')
   !
   ! Add the q-point to buf_q if not found.
   if (i1wf==1) then
     nqfound = 1
     buf_q(:,1)        = PH_Hdr%qptn(:)
     buf_posq(1)       = i1wf
     buf_npert(1)      = 1 
     buf_pertcase(1,1) = PH_Hdr%pertcase
     buf_pospertq(PH_Hdr%pertcase,1) = i1wf 
     !
   else 
     !
     if ( ANY( ABS(PH_Hdr%qptn(:) - buf_q(:,nqfound)) > tol16 )) then 
       ! Found a new q-point.
       nqfound = nqfound + 1
       ABI_CHECK(nqfound<buf_sz," Too many q-points, increase buf_sz")
       write (msg,'(a,3es16.8)')' Found new q-point: ',PH_Hdr%qptn
       call wrtout(std_out,msg,'COLL')

       buf_q(:,nqfound)         = PH_Hdr%qptn(:)
       buf_posq(nqfound)        = i1wf
       buf_npert(nqfound)       = 1 
       buf_pertcase(1,nqfound)  = PH_Hdr%pertcase
       buf_pospertq(PH_Hdr%pertcase,nqfound) = i1wf 
       !
     else 
       ! Previous q-point but with a different perturbation.
       next = buf_npert(nqfound) + 1 
       ABI_CHECK(next<=3*Fd%natom,"next>3*natom!")
       buf_npert(nqfound) = next
       buf_pertcase(next,nqfound) = PH_Hdr%pertcase
       buf_pospertq(PH_Hdr%pertcase,nqfound) = i1wf 
       !
     end if
   end if
   !
   !  Skip matrix elements.
   do spin=1,PH_Hdr%nsppol 
     do ik_rf=1,PH_Hdr%nkpt   
       read(Fd%funt) ! ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)
     end do 
   end do 

   write(msg,'(a,i0)')' read_gkk : Done completing the kpoints for pert ',PH_Hdr%pertcase
   call wrtout(std_out,msg,'COLL')

   call hdr_clean(PH_Hdr)
 end do ! i1wf

 !close(Fd%funt)
 !
 ! Copy info on q-points and the associated perturbations.
 Fd%nqibz = nqfound

 ABI_ALLOCATE(Fd%qibz,(3,Fd%nqibz))
 Fd%qibz      = buf_q(:,1:Fd%nqibz)
 ABI_ALLOCATE(Fd%npert_q,(Fd%nqibz))
 Fd%npert_q   = buf_npert(1:Fd%nqibz)
 !
 ! Add 1 to accound for the GS block.
 ABI_ALLOCATE(Fd%pos_q,(Fd%nqibz))
 Fd%pos_q     = buf_posq(1:Fd%nqibz) + 1
 ABI_ALLOCATE(Fd%pos_pertq,(3*Fd%natom,Fd%nqibz))
 Fd%pos_pertq = buf_pospertq(:,1:Fd%nqibz) + 1

 ABI_ALLOCATE(Fd%pertcase_q,(MAXVAL(Fd%npert_q),Fd%nqibz))
 Fd%pertcase_q = -HUGE(1)

 do iq_ibz=1,Fd%nqibz 
   npe = buf_npert(iq_ibz)
   Fd%pertcase_q(1:npe,iq_ibz) = buf_pertcase(1:npe,iq_ibz)
 end do

 ABI_DEALLOCATE(buf_npert)
 ABI_DEALLOCATE(buf_pertcase)
 ABI_DEALLOCATE(buf_pospertq)
 !
 ! Consistency check: Perturbations with same q must be packed together since it makes life easier.
 istat=0
 do iq_ibz=1,Fd%nqibz-1
   do next=iq_ibz+1,Fd%nqibz
     if ( ALL( ABS(Fd%qibz(:,iq_ibz) - Fd%qibz(:,next)) < tol16 )) then 
       istat=istat+1
       write(msg,'(3a,2i0,6es16.8)')" Found two identical qpoints in Fd%qibz, this is not allowed ",ch10,&
&       " iq_ibz, next, q1, q2 ",iq_ibz,next,Fd%qibz(:,iq_ibz),Fd%qibz(:,next)
       MSG_WARNING(msg)
     end if
   end do
 end do

 if (istat/=0) then
   msg = " Perturbations with same q must be packed together. Solution: generate a new GKK file"
   MSG_ERROR(msg)
 end if

end subroutine gkkfd_init
!!***

!----------------------------------------------------------------------

!!****f* m_io_gkk/gkkfile_skip_nblock
!! NAME
!!  gkkfile_skip_nblock
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_io_gkk
!!
!! CHILDREN
!!      gkkfile_skip_nblock,hdr_clean,hdr_nullify,insy3,mati3inv,symq3
!!      wrap2_pmhalf,wrtout
!!
!! SOURCE

subroutine gkkfile_skip_nblock(Fd,nblock,PH_Hdr,do_rewind)

 use defs_basis
 use m_header

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkkfile_skip_nblock'
 use interfaces_14_hidewrite
 use interfaces_59_io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nblock
 logical,optional,intent(in) :: do_rewind
 type(gkkfd_t),intent(inout) :: Fd
 type(hdr_type),intent(out) :: PH_Hdr

!Local variables-------------------------------
!scalars
 integer,parameter :: rdwr5=5 
 integer :: skip,fform,spin,ikpt,n1wf
 character(len=500) :: msg
 type(hdr_type) :: Tmp_Hdr

! ************************************************************************

 if (PRESENT(do_rewind)) then
   if (do_rewind) then 
     rewind(Fd%funt)
     Fd%in_gsblock = .TRUE.
   end if
 end if

 call hdr_nullify(Tmp_Hdr)

 do skip=1,nblock
   !
   ! Read the header without rewinding it
   call hdr_io(fform,Tmp_Hdr,rdwr5,Fd%funt)
   ABI_CHECK(fform>0,'GKK header was mis-read. fform == 0')
   !
   do spin=1,Tmp_Hdr%nsppol
     do ikpt=1,Tmp_Hdr%nkpt
       read(Fd%funt) ! ((eigen1(:,ii,ib),ii=1,nband),ib=1,nband)
     end do
   end do
   !
   call hdr_clean(Tmp_Hdr)

   if (skip==1.and.Fd%in_gsblock) then ! GS block is followed by n1wf.
     read(Fd%funt) n1wf
     write(msg,'(a,i0)')' gkkfile_skip_nblock : number of perturbations in the gkk file = ',n1wf
     call wrtout(std_out,msg,'COLL')
     ABI_CHECK(n1wf==Fd%n1wf,"n1wf/=Fd%n1wf")
   end if
   !
 end do
 !
 ! Read the PH header without rewinding it
 call hdr_io(fform,PH_Hdr,rdwr5,Fd%funt)
 ABI_CHECK(fform>0,'GKK header was mis-read. fform == 0')

end subroutine gkkfile_skip_nblock
!!***

!!****f* m_io_gkk/gkkfd_read_h1me
!! NAME
!!  gkk_read_h1me
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gkkfd_read_h1me(Fd,iq_ibz,pertcase,nFSband,minFSband,nkwant,kwanted,Cryst,h1_mat_el)

 use defs_basis
 use m_header

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkkfd_read_h1me'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_72_response
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_ibz,pertcase,nFSband,minFSband,nkwant
 type(gkkfd_t),intent(inout) :: Fd
 type(crystal_structure),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: kwanted(3,nkwant)
 real(dp),intent(out) :: h1_mat_el(2,nFSband**2,nkwant,Fd%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: rfmeth2=2,H1_MISS=-1,H1_SYM=+1,H1_READ=3
 integer :: ib1,ib2,ibb,idir,ipert,nblock,spin,ph_ikibz,ph_nband,ikw
 integer :: ph_sym,ph_nsym,ph_tim,qtimrev,dummy_syuse,istat,jkw_src,ik_img !,ii
 real(dp) :: res,ss,timsign
 logical :: found
 character(len=500) :: msg
 type(hdr_type) :: PH_Hdr
!arrays 
 integer,allocatable :: FSirrtok(:,:),h1_flag(:),ph2kwant(:)
 integer :: ph_symafm(Cryst%nsym),symq(4,2,Cryst%nsym)
 integer :: ph_symrec(3,3,Cryst%nsym),ph_symrel(3,3,Cryst%nsym)
 real(dp) :: ph_tnons(3,Cryst%nsym),kpt(3),redkpt(3)
 real(dp),allocatable :: buf_h1me(:,:,:)

! ************************************************************************

 call hdr_nullify(PH_Hdr)

 nblock = Fd%pos_pertq(pertcase,iq_ibz) - 1
 if (nblock<=0) then
   write(msg,'(a,2i0,a)')" (iq_ibz, pertcase) ",iq_ibz,pertcase," not stored on file!"
   MSG_ERROR(msg)
 end if
 write(std_out,*)" pertcase,iq_ibz,nblock ",pertcase,iq_ibz,nblock
 !
 call gkkfile_skip_nblock(Fd,nblock,PH_Hdr,do_rewind=.TRUE.)

 if (PH_Hdr%pertcase /= pertcase) then
   write(msg,'(a,2(i0,1x))')" PH_Hdr%pertcase/=pertcase ",PH_Hdr%pertcase,pertcase
   MSG_BUG(msg)
 end if

 if ( ANY( ABS(PH_Hdr%qptn - Fd%qibz(:,iq_ibz)) > tol16) ) then
   write(msg,'(a,2(3es15.8,2x))')" mismatch in q-points : PH_Hdr%qptn /= Fd%qibz(:,iq_ibz) ",PH_Hdr%qptn,Fd%qibz(:,iq_ibz)
   MSG_BUG(msg)
 end if

 ph_nband = PH_Hdr%nband(1)
 if ( ANY(PH_Hdr%nband /= PH_Hdr%nband(1)) ) then
   MSG_ERROR("Different number of bands not supported")
 end if

 if (minFSband+nFSband>ph_nband) then
   MSG_ERROR("Requiring too many bands!")
 end if
 !
 ! Examine the symmetries of the q wavevector
 ! these will be used to complete the perturbations for other atoms and idir
 !
 call symq3(Cryst%nsym,Fd%qibz(:,iq_ibz),symq,Cryst%symrec,qtimrev,prtvol=0)

 !  Examine the symmetries of the full perturbation. these will be used to complete the kpoints
 !  DOESNT USE TIME REVERSAL IN insy3 except for gamma
 !
 ! pertcase = idir + (ipert-1)*3 where ipert=iatom in the interesting cases
 idir  = mod(pertcase-1,3)+1
 ipert = int(dble(pertcase-idir)/three)+1

 dummy_syuse=0

 call insy3(Cryst%gprimd,idir,Cryst%indsym,ipert,Cryst%natom,Cryst%nsym,ph_nsym,rfmeth2,Cryst%symafm,&
&  ph_symafm,symq,Cryst%symrec,Cryst%symrel,ph_symrel,dummy_syuse,Cryst%tnons,ph_tnons)

 do ph_sym=1,ph_nsym
   call mati3inv(ph_symrel(:,:,ph_sym),ph_symrec(:,:,ph_sym))
 end do
 !
 ! Mapping between the k-points stored on file and the list of k-points required.
 ABI_ALLOCATE(ph2kwant,(PH_Hdr%nkpt))
 ph2kwant=0

 do ph_ikibz=1,PH_Hdr%nkpt ! Loop over irred k-points, WARNING  PH_Hdr%nkpt depends on q-point and symmetries
   !
   ! Check to see if kpoint is in FS set
   ! WARNING! the kpoints in the file (kptns) could be ordered arbitrarily
   do ikw=1,nkwant
     !kpt(:) = PH_Hdr%kptns(:,ph_ikibz)- kwanted(:,ikw)-qptirred_local(:,iqptirred)
     kpt(:) = PH_Hdr%kptns(:,ph_ikibz) - kwanted(:,ikw) - Fd%qibz(:,iq_ibz)
     call wrap2_pmhalf(kpt(1),redkpt(1),res)
     call wrap2_pmhalf(kpt(2),redkpt(2),res)
     call wrap2_pmhalf(kpt(3),redkpt(3),res)
     ss=redkpt(1)**2+redkpt(2)**2+redkpt(3)**2
     !
     if (ss < tol6) then 
       ph2kwant(ph_ikibz) = ikw
       EXIT  ! The ikpt has been identified to ikw, so no need to search the other kwanted
     end if
     !
   end do !ikw
 end do !ph_ikibz
 !
 !  ========================================================
 !  Loop over irred kpts in file, and fill the default gkk
 !  ========================================================
 ! Check to see if kpoint is in FS set
 ! WARNING! the kpoints in the file (kptns) could be ordered arbitrarily
 ABI_ALLOCATE(h1_flag,(nkwant))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory h1_flag")
 h1_flag=H1_MISS

 ABI_ALLOCATE(buf_h1me,(2,ph_nband,ph_nband))

 h1_mat_el = HUGE(one)
 do spin=1,PH_Hdr%nsppol 
   !
   do ph_ikibz=1,PH_Hdr%nkpt ! Loop over irred k-points, WARNING  PH_Hdr%nkpt depends on q-point and symmetries
     !     
     ! this is the main read of the gkk matrix elements from the file (buf_h1me arrays)
     ! it has to be done exactly nsppol*nkpt times, and the kpt_phon are completed
     ! where appropriate in the loop below (normally succeeding only once for each kpt)
     !     
     ikw = ph2kwant(ph_ikibz) 

     if (ikw == 0 ) then
       read(Fd%funt) ! buf_h1me
     else 
       read(Fd%funt) buf_h1me
       !
       ! save this kpoint
       do ib1=1,nFSband
         do ib2=1,nFSband
           ibb = (ib1-1)*nFSband+ib2
           h1_mat_el(:,ibb,ikw,spin) = buf_h1me(:,minFSband-1+ib2,minFSband-1+ib1)
         end do !ib2
       end do !ib1
       h1_flag(ikw) = H1_READ
       !
       ! ===============================================================
       ! we now have contribution to g(k+q,k; \kappa,\alpha) from one
       ! kpoint,and one perturbation,
       ! NB: each perturbation will contribute to all the modes later!
       ! 
       ! SHOULD ONLY DO THIS FOR THE SYMS NEEDED 
       ! TO COMPLETE THE PERTURBATIONS!!!
       ! ================================================================
       !
       ! can we safely exit here? The ikpt has been identified to ikw,
       ! so no need to search the other kwanted
     end if
   end do !ph_ikibz
   !
 end do !spin

 ABI_DEALLOCATE(ph2kwant)
 !
 ! Find correspondence between points in kwanted provided sym conserves pert as well as qpoint
 ! TODO: Should introduce table to loop on shells.
 !
 ABI_ALLOCATE(FSirrtok,(3,nkwant))
 istat = ABI_ALLOC_STAT
 ABI_CHECK(istat==0,"out of memory in FSirrtok")
 FSirrtok=0

 do ikw=1,nkwant
   !
   if (h1_flag(ikw) == H1_MISS) CYCLE
   !
   do ph_sym=1,ph_nsym
     do ph_tim=0,qtimrev
       !
       timsign = one-two*ph_tim
       kpt = timsign* MATMUL( ph_symrec(:,:,ph_sym), kwanted(:,ikw) )
       call wrap2_pmhalf(kpt(1),redkpt(1),res)
       call wrap2_pmhalf(kpt(2),redkpt(2),res)
       call wrap2_pmhalf(kpt(3),redkpt(3),res)
       found=.FALSE.
       !
       do jkw_src=1,nkwant ! FIXME: use rank scheme to avoid this loop
         ss=  ( &
&           redkpt(1)-kwanted(1,jkw_src))**2 &
&         +(redkpt(2)-kwanted(2,jkw_src))**2 &
&         +(redkpt(3)-kwanted(3,jkw_src))**2

         if (ss < tol6) then
           found=.TRUE.
           FSirrtok(1,jkw_src) = ikw
           FSirrtok(2,jkw_src) = ph_sym
           FSirrtok(3,jkw_src) = ph_tim
           !EXIT !exit jkw_src
         end if
         !
       end do !jkw_src
       !
       if (.not.found) then
         write (msg,'(a,3es16.6,a,i5,a,i4,a)')&
&         ' equivalent of kpt ',kwanted(:,ikw),' by sym ',ph_sym,' and itime ',ph_tim,' was not found'
         MSG_ERROR(msg)
       end if
       !
     end do ! ph_tim
   end do ! isim1
   !
 end do !ikw

 !if (elph_ds%tuniformgrid == 1) then  ! check if irred kpoints found do reconstitute the FS kpts
 do ikw=1,nkwant
   if (FSirrtok(1,ikw) == 0) then
     write(msg,'(a,3es16.6,2a)')&
&     ' kwanted = ',kwanted(:,ikw),ch10,&
&     ' is not the symmetric of one of those found in the GKK file'
     MSG_ERROR(msg)
   end if
 end do !ikw
 !end if
 !
 ! ===============================================================
 ! We now have all irred kpoints : complete the others
 ! complete gkk for symmetric ikw with sym1 which conserve
 ! the full perturbation+qpoint
 ! Not tested explicitly, but the results for Pb using reduced kpts look good
 ! should do same RF calculation with nsym=1 and check
 ! ===============================================================
 do ikw=1,nkwant
   if (h1_flag(ikw) == H1_MISS) then  !   if the element has already been filled with another sym op, cycle
     ik_img  = FSirrtok(1,ikw)
     ph_sym  = FSirrtok(2,ikw)
     ph_tim  = FSirrtok(3,ikw)
     timsign = one-two*ph_tim
     if (h1_flag(ik_img) == H1_MISS) then
      MSG_ERROR("Cannot symmetrize")
      !cycle
     end if
     !
     ! copy kpt to symmetric kpoint
     do spin=1,Fd%nsppol
       do ibb=1,nFSband**2
         h1_mat_el(:,ibb,ikw,spin) = h1_mat_el(:,ibb,ik_img,spin)
       end do
     end do
     h1_flag(ikw) = H1_SYM
   end if
 end do !ikw
 !
 ! normally at this point we have used all the gkk for all kpoints on the FS for the given irred perturbation: check
 do ikw=1,nkwant
   if (h1_flag(ikw) == H1_MISS) then
     write (msg,'(2(a,i0),a,3es15.8,2a)')&
&     ' the h1 element with : pertcase = ',pertcase,' ikw = ',ikw,' kpt = ',kwanted(:,ikw),ch10,&
&     ' was not found by symmetry operations on the irreducible k-points given.'
     MSG_ERROR(msg)
   end if
   !if (FSirrtok(1,ikw) == 0) then
   !  write (msg,'(a)')" One of the kpoints was not equivalent to an irred one found in the gkk file "
   !  MSG_ERROR(msg)
   !end if
 end do !ikw

 write(msg,'(a,i0)')' read_gkk : Done readining and completing the k-points for pertcase: ',pertcase
 call wrtout(std_out,msg,'COLL')

 ABI_DEALLOCATE(h1_flag)
 ABI_DEALLOCATE(FSirrtok)
 ABI_DEALLOCATE(buf_h1me)

 call hdr_clean(PH_Hdr)

end subroutine gkkfd_read_h1me
!!***

END MODULE m_io_gkk
!!***

