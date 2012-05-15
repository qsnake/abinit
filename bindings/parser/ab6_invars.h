#ifndef AB6_INVARS_H
#define AB6_INVARS_H

#include <stdlib.h>

#include "ab6_base.h"

/**
 * Ab6InvarsTypes:
 * @_INT_SCALAR: a 32 bits integer.
 * @_DOUBLE_SCALAR: a 64 bits float.
 * @_INT_ARRAY: an array of 32 bits integers.
 * @_DOUBLE_ARRAY: an array of 64 bits floats.
 *
 * The possible types of the attributes of datasets.
 */
typedef enum
  {
    _INT_SCALAR,
    _INT_ARRAY,
    _DOUBLE_SCALAR,
    _DOUBLE_ARRAY,
    _OTHER
  } Ab6InvarsTypes;

/* This file has been automatically generated, do not modify. */
typedef enum
{
  AB6_INVARS_SYMCHI         ,  /* _INT_SCALAR     */
  AB6_INVARS_GWENCOMP       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_SLABWSRAD      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_RFSTRS         ,  /* _INT_SCALAR     */
  AB6_INVARS_EXCHN2N3D      ,  /* _INT_SCALAR     */
  AB6_INVARS_NSTEP          ,  /* _INT_SCALAR     */
  AB6_INVARS_WVL_NPRCCG     ,  /* _INT_SCALAR     */
  AB6_INVARS_GWGAMMA        ,  /* _INT_SCALAR     */
  AB6_INVARS_GW_TOLDFEIG    ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NSHEPS         ,  /* _INT_SCALAR     */
  AB6_INVARS_TOLMXF         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_CHARGE         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PRT1DM         ,  /* _INT_SCALAR     */
  AB6_INVARS_PREPANL        ,  /* _INT_SCALAR     */
  AB6_INVARS_NSPINOR        ,  /* _INT_SCALAR     */
  AB6_INVARS_NDTSET         ,  /* _INT_SCALAR     */
  AB6_INVARS_USEPAWU        ,  /* _INT_SCALAR     */
  AB6_INVARS_MPW            ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWCROSS       ,  /* _INT_SCALAR     */
  AB6_INVARS_OCCOPT         ,  /* _INT_SCALAR     */
  AB6_INVARS_POSTOLDFF      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_POSTOLDFE      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_BOXCENTER      ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_BXCTMINDG      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_QPTDM          ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_TD_MEXCIT      ,  /* _INT_SCALAR     */
  AB6_INVARS_NCTIME         ,  /* _INT_SCALAR     */
  AB6_INVARS_FRZFERMI       ,  /* _INT_SCALAR     */
  AB6_INVARS_RFPHON         ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTDEN         ,  /* _INT_SCALAR     */
  AB6_INVARS_GWPARA         ,  /* _INT_SCALAR     */
  AB6_INVARS_NPFFT          ,  /* _INT_SCALAR     */
  AB6_INVARS_IRDBSRESO      ,  /* _INT_SCALAR     */
  AB6_INVARS_RECPTROTT      ,  /* _INT_SCALAR     */
  AB6_INVARS_TYPAT          ,  /* _INT_ARRAY      */
  AB6_INVARS_JDTSET         ,  /* _INT_SCALAR     */
  AB6_INVARS_BS_HAYDOCK_NITER,  /* _INT_SCALAR     */
  AB6_INVARS_PAWPRTDEN      ,  /* _INT_SCALAR     */
  AB6_INVARS_USE_GPU_CUDA   ,  /* _INT_SCALAR     */
  AB6_INVARS_ALGALCH        ,  /* _INT_ARRAY      */
  AB6_INVARS_PREPSCPHON     ,  /* _INT_SCALAR     */
  AB6_INVARS_KPTRLATT       ,  /* _INT_ARRAY      */
  AB6_INVARS_BS_ALGORITHM   ,  /* _INT_SCALAR     */
  AB6_INVARS_GWMEM          ,  /* _INT_SCALAR     */
  AB6_INVARS_NLOALG         ,  /* _INT_ARRAY      */
  AB6_INVARS_NPKPT          ,  /* _INT_SCALAR     */
  AB6_INVARS_IEXTRAPWF      ,  /* _INT_SCALAR     */
  AB6_INVARS_USEKDEN        ,  /* _INT_SCALAR     */
  AB6_INVARS_XRED_ORIG      ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_KPTNS          ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_NBERRY         ,  /* _INT_SCALAR     */
  AB6_INVARS_CHKSYMBREAK    ,  /* _INT_SCALAR     */
  AB6_INVARS_FREQREMIN      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_USEWVL         ,  /* _INT_SCALAR     */
  AB6_INVARS_USEEXEXCH      ,  /* _INT_SCALAR     */
  AB6_INVARS_EXCHMIX        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_BS_FREQ_MESH   ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_TFKINFUNC      ,  /* _INT_SCALAR     */
  AB6_INVARS_PTCHARGE       ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_KPTNRM         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_IRDBSEIG       ,  /* _INT_SCALAR     */
  AB6_INVARS_SLABZBEG       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_IONMOV         ,  /* _INT_SCALAR     */
  AB6_INVARS_VIS            ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_CD_IMFRQS      ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_USERRE         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NPWWFN         ,  /* _INT_SCALAR     */
  AB6_INVARS_MFFMEM         ,  /* _INT_SCALAR     */
  AB6_INVARS_BS_EXCHANGE_TERM,  /* _INT_SCALAR     */
  AB6_INVARS_MK1MEM         ,  /* _INT_SCALAR     */
  AB6_INVARS_NFREQSUS       ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTKPT         ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTCS          ,  /* _INT_SCALAR     */
  AB6_INVARS_IMGMOV         ,  /* _INT_SCALAR     */
  AB6_INVARS_NPSPINOR       ,  /* _INT_SCALAR     */
  AB6_INVARS_QMASS          ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_GW_QLWL        ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_ECUTSIGX       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NDYSON         ,  /* _INT_SCALAR     */
  AB6_INVARS_WVL_HGRID      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_USEREC         ,  /* _INT_SCALAR     */
  AB6_INVARS_GETXCART       ,  /* _INT_SCALAR     */
  AB6_INVARS_SPMETH         ,  /* _INT_SCALAR     */
  AB6_INVARS_MKMEM          ,  /* _INT_SCALAR     */
  AB6_INVARS_QPRTRB         ,  /* _INT_ARRAY      */
  AB6_INVARS_CHKGWCOMP      ,  /* _INT_SCALAR     */
  AB6_INVARS_GETDEN         ,  /* _INT_SCALAR     */
  AB6_INVARS_RECEFERMI      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PRTGEO         ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTFC          ,  /* _INT_SCALAR     */
  AB6_INVARS_VACWIDTH       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_MIXALCH        ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_USEDMATPU      ,  /* _INT_SCALAR     */
  AB6_INVARS_TD_MAXENE      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PARAL_KGB      ,  /* _INT_SCALAR     */
  AB6_INVARS_GET1WF         ,  /* _INT_SCALAR     */
  AB6_INVARS_IRDDEN         ,  /* _INT_SCALAR     */
  AB6_INVARS_ZCUT           ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_GETBSEIG       ,  /* _INT_SCALAR     */
  AB6_INVARS_GETXRED        ,  /* _INT_SCALAR     */
  AB6_INVARS_BMASS          ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_W90PRTUNK      ,  /* _INT_SCALAR     */
  AB6_INVARS_RATSPH         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_INTEXACT       ,  /* _INT_SCALAR     */
  AB6_INVARS_USEYLM         ,  /* _INT_SCALAR     */
  AB6_INVARS_MACRO_UJ       ,  /* _INT_SCALAR     */
  AB6_INVARS_DMFT_RSLF      ,  /* _INT_SCALAR     */
  AB6_INVARS_FREQSUSLO      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NPIMAGE        ,  /* _INT_SCALAR     */
  AB6_INVARS_GETCELL        ,  /* _INT_SCALAR     */
  AB6_INVARS_MAXNSYM        ,  /* _INT_SCALAR     */
  AB6_INVARS_TOLSYM         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_RECTOLDEN      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_KSSFORM        ,  /* _INT_SCALAR     */
  AB6_INVARS_USERRD         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_USERRA         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_USERRC         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_USERRB         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_BERRYSTEP      ,  /* _INT_SCALAR     */
  AB6_INVARS_DMFT_MXSF      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_WFOPTALG       ,  /* _INT_SCALAR     */
  AB6_INVARS_VCUTGEO        ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_RF1ELFD        ,  /* _INT_SCALAR     */
  AB6_INVARS_TNONS          ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_IRDWFK         ,  /* _INT_SCALAR     */
  AB6_INVARS_IRD1DEN        ,  /* _INT_SCALAR     */
  AB6_INVARS_NTIMIMAGE      ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWLCUTD       ,  /* _INT_SCALAR     */
  AB6_INVARS_NBAND          ,  /* _INT_ARRAY      */
  AB6_INVARS_PRTEIG         ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWUJV         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_ENUNIT         ,  /* _INT_SCALAR     */
  AB6_INVARS_VDW_NFRAG      ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTNEST        ,  /* _INT_SCALAR     */
  AB6_INVARS_IRDWFQ         ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWSTGYLM      ,  /* _INT_SCALAR     */
  AB6_INVARS_IPRCEL         ,  /* _INT_SCALAR     */
  AB6_INVARS_IATFIX         ,  /* _INT_ARRAY      */
  AB6_INVARS_LDGAPP         ,  /* _INT_SCALAR     */
  AB6_INVARS_RESTARTXF      ,  /* _INT_SCALAR     */
  AB6_INVARS_IRDDDK         ,  /* _INT_SCALAR     */
  AB6_INVARS_DYNIMAGE       ,  /* _INT_ARRAY      */
  AB6_INVARS_SPGROUP        ,  /* _INT_SCALAR     */
  AB6_INVARS_NPWSIGX        ,  /* _INT_SCALAR     */
  AB6_INVARS_NOMEGASI       ,  /* _INT_SCALAR     */
  AB6_INVARS_NOMEGASF       ,  /* _INT_SCALAR     */
  AB6_INVARS_IATSPH         ,  /* _INT_ARRAY      */
  AB6_INVARS_NLINE          ,  /* _INT_SCALAR     */
  AB6_INVARS_GW_EET_NBAND   ,  /* _INT_SCALAR     */
  AB6_INVARS_KPT            ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_BDBERRY        ,  /* _INT_ARRAY      */
  AB6_INVARS_RFMETH         ,  /* _INT_SCALAR     */
  AB6_INVARS_USERIA         ,  /* _INT_SCALAR     */
  AB6_INVARS_NFREQSP        ,  /* _INT_SCALAR     */
  AB6_INVARS_USERIC         ,  /* _INT_SCALAR     */
  AB6_INVARS_USERID         ,  /* _INT_SCALAR     */
  AB6_INVARS_USERIE         ,  /* _INT_SCALAR     */
  AB6_INVARS_GW_FREQSP      ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_RANDOM_ATPOS   ,  /* _INT_SCALAR     */
  AB6_INVARS_NTYPALCH       ,  /* _INT_SCALAR     */
  AB6_INVARS_LOCALRDWF      ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTDOSM        ,  /* _INT_SCALAR     */
  AB6_INVARS_NIMAGE         ,  /* _INT_SCALAR     */
  AB6_INVARS_MDTEMP         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_CD_HALFWAY_FREQ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_DIEMIX         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_DENSTY         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_POSOCC         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_IRDKSS         ,  /* _INT_SCALAR     */
  AB6_INVARS_RECTESTEG      ,  /* _INT_SCALAR     */
  AB6_INVARS_NTYPAT         ,  /* _INT_SCALAR     */
  AB6_INVARS_ICOULOMB       ,  /* _INT_SCALAR     */
  AB6_INVARS_CORECS         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_IKHXC          ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTXML         ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWXCDEV       ,  /* _INT_SCALAR     */
  AB6_INVARS_INCLVKB        ,  /* _INT_SCALAR     */
  AB6_INVARS_NDYNIMAGE      ,  /* _INT_SCALAR     */
  AB6_INVARS_ECUTSM         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_FREQREMAX      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NBDBUF         ,  /* _INT_SCALAR     */
  AB6_INVARS_NATPAWU        ,  /* _INT_SCALAR     */
  AB6_INVARS_DIELAM         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_RF1PHON        ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWNTHETA      ,  /* _INT_SCALAR     */
  AB6_INVARS_BDEIGRF        ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWUJAT        ,  /* _INT_SCALAR     */
  AB6_INVARS_IXC            ,  /* _INT_SCALAR     */
  AB6_INVARS_DELAYPERM      ,  /* _INT_SCALAR     */
  AB6_INVARS_NSCFORDER      ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWOVLP        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_ATVSHIFT       ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PRTPMP         ,  /* _INT_SCALAR     */
  AB6_INVARS_CPUS           ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_ECUTEPS        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_RF3ELFD        ,  /* _INT_SCALAR     */
  AB6_INVARS_SUSKXCRS       ,  /* _INT_SCALAR     */
  AB6_INVARS_IPRCFC         ,  /* _INT_SCALAR     */
  AB6_INVARS_USEXCNHAT      ,  /* _INT_SCALAR     */
  AB6_INVARS_IRDBSCOUP      ,  /* _INT_SCALAR     */
  AB6_INVARS_GENAFM         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_SCPHON_TEMP    ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NGFFTDG        ,  /* _INT_ARRAY      */
  AB6_INVARS_RF3ATPOL       ,  /* _INT_ARRAY      */
  AB6_INVARS_DMFT_SOLV      ,  /* _INT_SCALAR     */
  AB6_INVARS_PREPGKK        ,  /* _INT_SCALAR     */
  AB6_INVARS_GW_RECONST_SCR ,  /* _INT_SCALAR     */
  AB6_INVARS_RPRIMD_ORIG    ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PRTVXC         ,  /* _INT_SCALAR     */
  AB6_INVARS_NFFTDG         ,  /* _INT_SCALAR     */
  AB6_INVARS_PARAL_RF       ,  /* _INT_SCALAR     */
  AB6_INVARS_ISTATIMG       ,  /* _INT_SCALAR     */
  AB6_INVARS_VDW_TYPFRAG    ,  /* _INT_ARRAY      */
  AB6_INVARS_OCC_ORIG       ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_IRDSCR         ,  /* _INT_SCALAR     */
  AB6_INVARS_GETGAM_EIG2NKQ ,  /* _INT_SCALAR     */
  AB6_INVARS_USEPAW         ,  /* _INT_SCALAR     */
  AB6_INVARS_USE_SLK        ,  /* _INT_SCALAR     */
  AB6_INVARS_RCUT           ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PAWSUSHAT      ,  /* _INT_SCALAR     */
  AB6_INVARS_DIEMIXMAG      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NSYM           ,  /* _INT_SCALAR     */
  AB6_INVARS_CD_MAX_FREQ    ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_CD_USE_TANGRID ,  /* _INT_SCALAR     */
  AB6_INVARS_IRANDOM        ,  /* _INT_SCALAR     */
  AB6_INVARS_RFUSER         ,  /* _INT_SCALAR     */
  AB6_INVARS_DMATPAWU       ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PRTFSURF       ,  /* _INT_SCALAR     */
  AB6_INVARS_NPBAND         ,  /* _INT_SCALAR     */
  AB6_INVARS_GETBSRESO      ,  /* _INT_SCALAR     */
  AB6_INVARS_NFFT           ,  /* _INT_SCALAR     */
  AB6_INVARS_RFASR          ,  /* _INT_SCALAR     */
  AB6_INVARS_ACCESSWFF      ,  /* _INT_SCALAR     */
  AB6_INVARS_FFTGW          ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTKDEN        ,  /* _INT_SCALAR     */
  AB6_INVARS_RFELFD         ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTVHXC        ,  /* _INT_SCALAR     */
  AB6_INVARS_MGFFT          ,  /* _INT_SCALAR     */
  AB6_INVARS_DMATUDIAG      ,  /* _INT_SCALAR     */
  AB6_INVARS_NBANDSUS       ,  /* _INT_SCALAR     */
  AB6_INVARS_RECGRATIO      ,  /* _INT_SCALAR     */
  AB6_INVARS_IPRCTFVW       ,  /* _INT_SCALAR     */
  AB6_INVARS_NSHSIGX        ,  /* _INT_SCALAR     */
  AB6_INVARS_NSHWFN         ,  /* _INT_SCALAR     */
  AB6_INVARS_MDWALL         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_RF1DIR         ,  /* _INT_ARRAY      */
  AB6_INVARS_SCISS          ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PAWSPNORB      ,  /* _INT_SCALAR     */
  AB6_INVARS_RFATPOL        ,  /* _INT_ARRAY      */
  AB6_INVARS_NPWKSS         ,  /* _INT_SCALAR     */
  AB6_INVARS_FRICTION       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_GW_NQLWL       ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTCML         ,  /* _INT_SCALAR     */
  AB6_INVARS_GWCALCTYP      ,  /* _INT_SCALAR     */
  AB6_INVARS_BS_HAYD_TERM   ,  /* _INT_SCALAR     */
  AB6_INVARS_IRDSUSCEP      ,  /* _INT_SCALAR     */
  AB6_INVARS_ZEEMANFIELD    ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_CHKPRIM        ,  /* _INT_SCALAR     */
  AB6_INVARS_RF1ATPOL       ,  /* _INT_ARRAY      */
  AB6_INVARS_DOSDELTAE      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PRTEFG         ,  /* _INT_SCALAR     */
  AB6_INVARS_KPTRLEN        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_KPTGW          ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_BS_COUPLING    ,  /* _INT_SCALAR     */
  AB6_INVARS_GETVEL         ,  /* _INT_SCALAR     */
  AB6_INVARS_DMFTBANDF      ,  /* _INT_SCALAR     */
  AB6_INVARS_BDGW           ,  /* _INT_ARRAY      */
  AB6_INVARS_GETPAWDEN      ,  /* _INT_SCALAR     */
  AB6_INVARS_ESMEAR         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_DMFTCHECK      ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTWF          ,  /* _INT_SCALAR     */
  AB6_INVARS_WTATCON        ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_NBDBLOCK       ,  /* _INT_SCALAR     */
  AB6_INVARS_GWCOMP         ,  /* _INT_SCALAR     */
  AB6_INVARS_WVL_CRMULT     ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_MQGRIDDG       ,  /* _INT_SCALAR     */
  AB6_INVARS_BS_NSTATES     ,  /* _INT_SCALAR     */
  AB6_INVARS_GW_NSTEP       ,  /* _INT_SCALAR     */
  AB6_INVARS_EFIELD         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_TL_NPRCCG      ,  /* _INT_SCALAR     */
  AB6_INVARS_BFIELD         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PAWPRTWF       ,  /* _INT_SCALAR     */
  AB6_INVARS_STRFACT        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_IRD1WF         ,  /* _INT_SCALAR     */
  AB6_INVARS_RF2PHON        ,  /* _INT_SCALAR     */
  AB6_INVARS_NELECT         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_RFDDK          ,  /* _INT_SCALAR     */
  AB6_INVARS_DIELNG         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PAWPRTDOS      ,  /* _INT_SCALAR     */
  AB6_INVARS_JPAWU          ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_NQPT           ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWLMIX        ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTLDEN        ,  /* _INT_SCALAR     */
  AB6_INVARS_OPTSTRESS      ,  /* _INT_SCALAR     */
  AB6_INVARS_POSITRON       ,  /* _INT_SCALAR     */
  AB6_INVARS_EFFMASS        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PRTATLIST      ,  /* _INT_ARRAY      */
  AB6_INVARS_SMDELTA        ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWUJRAD       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_WTK            ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_GETWFK         ,  /* _INT_SCALAR     */
  AB6_INVARS_FIXMOM         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_DIISMEMORY     ,  /* _INT_SCALAR     */
  AB6_INVARS_GETWFQ         ,  /* _INT_SCALAR     */
  AB6_INVARS_FBAND          ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NBANDKSS       ,  /* _INT_SCALAR     */
  AB6_INVARS_FREQSPMAX      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_GW_SIGXCORE    ,  /* _INT_SCALAR     */
  AB6_INVARS_PITRANSFORM    ,  /* _INT_SCALAR     */
  AB6_INVARS_DIEGAP         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_ZIONTYPAT      ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PAWMIXDG       ,  /* _INT_SCALAR     */
  AB6_INVARS_STRPRECON      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PRTDIPOLE      ,  /* _INT_SCALAR     */
  AB6_INVARS_ISTATSHFT      ,  /* _INT_SCALAR     */
  AB6_INVARS_NPWEPS         ,  /* _INT_SCALAR     */
  AB6_INVARS_RF3DIR         ,  /* _INT_ARRAY      */
  AB6_INVARS_SPNORBSCL      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_BANDPP         ,  /* _INT_SCALAR     */
  AB6_INVARS_OMEGASIMAX     ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_ECUTWFN        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_AWTR           ,  /* _INT_SCALAR     */
  AB6_INVARS_NKPTGW         ,  /* _INT_SCALAR     */
  AB6_INVARS_NPULAYIT       ,  /* _INT_SCALAR     */
  AB6_INVARS_DMFT_NWLI      ,  /* _INT_SCALAR     */
  AB6_INVARS_RDMNB          ,  /* _INT_SCALAR     */
  AB6_INVARS_RF3PHON        ,  /* _INT_SCALAR     */
  AB6_INVARS_NNOS           ,  /* _INT_SCALAR     */
  AB6_INVARS_GWRPACORR      ,  /* _INT_SCALAR     */
  AB6_INVARS_DILATMX        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_DMATPUOPT      ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWECUTDG      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_TOLIMG         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PRTSTM         ,  /* _INT_SCALAR     */
  AB6_INVARS_VDW_SUPERCELL  ,  /* _INT_ARRAY      */
  AB6_INVARS_ACELL_ORIG     ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_DIEMAC         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_SIGNPERM       ,  /* _INT_SCALAR     */
  AB6_INVARS_RF2DIR         ,  /* _INT_ARRAY      */
  AB6_INVARS_IRDPAWDEN      ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTGDEN        ,  /* _INT_SCALAR     */
  AB6_INVARS_GW_CUSTOM_FREQSP,  /* _INT_SCALAR     */
  AB6_INVARS_SOENERGY       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_BERRYOPT       ,  /* _INT_SCALAR     */
  AB6_INVARS_IXCPOSITRON    ,  /* _INT_SCALAR     */
  AB6_INVARS_GETQPS         ,  /* _INT_SCALAR     */
  AB6_INVARS_TIMOPT         ,  /* _INT_SCALAR     */
  AB6_INVARS_RFDIR          ,  /* _INT_ARRAY      */
  AB6_INVARS_SHIFTK         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_CD_CUSTOM_IMFRQS,  /* _INT_SCALAR     */
  AB6_INVARS_KBERRY         ,  /* _INT_ARRAY      */
  AB6_INVARS_BS_LOBAND      ,  /* _INT_SCALAR     */
  AB6_INVARS_MGFFTDG        ,  /* _INT_SCALAR     */
  AB6_INVARS_GETOCC         ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTPOSCAR      ,  /* _INT_SCALAR     */
  AB6_INVARS_NOMEGASRD      ,  /* _INT_SCALAR     */
  AB6_INVARS_IBOXCUT        ,  /* _INT_SCALAR     */
  AB6_INVARS_ISCF           ,  /* _INT_SCALAR     */
  AB6_INVARS_USEDMFT        ,  /* _INT_SCALAR     */
  AB6_INVARS_NSHIFTK        ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTELF         ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTCIF         ,  /* _INT_SCALAR     */
  AB6_INVARS_OMEGASRDMAX    ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NATRD          ,  /* _INT_SCALAR     */
  AB6_INVARS_NOSEINERT      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_INTXC          ,  /* _INT_SCALAR     */
  AB6_INVARS_NGROUP_RF      ,  /* _INT_SCALAR     */
  AB6_INVARS_SPGAXOR        ,  /* _INT_SCALAR     */
  AB6_INVARS_GET1DEN        ,  /* _INT_SCALAR     */
  AB6_INVARS_GW_EET         ,  /* _INT_SCALAR     */
  AB6_INVARS_FREQSUSIN      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_DMFTBANDI      ,  /* _INT_SCALAR     */
  AB6_INVARS_IRDQPS         ,  /* _INT_SCALAR     */
  AB6_INVARS_VEL_ORIG       ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_OPTNLXCCC      ,  /* _INT_SCALAR     */
  AB6_INVARS_DIECUT         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_GW_USE_POLE_SCR,  /* _INT_SCALAR     */
  AB6_INVARS_JELLSLAB       ,  /* _INT_SCALAR     */
  AB6_INVARS_PPMODEL        ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWOPTMIX      ,  /* _INT_SCALAR     */
  AB6_INVARS_OPTCELL        ,  /* _INT_SCALAR     */
  AB6_INVARS_SCPHON_SUPERCELL,  /* _INT_ARRAY      */
  AB6_INVARS_CD_SUBSET_FREQ ,  /* _INT_ARRAY      */
  AB6_INVARS_TOLWFR         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PAWFATBND      ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTVOLIMG      ,  /* _INT_SCALAR     */
  AB6_INVARS_NATVSHIFT      ,  /* _INT_SCALAR     */
  AB6_INVARS_GW_SCTYPE      ,  /* _INT_SCALAR     */
  AB6_INVARS_STRTARGET      ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_KPTOPT         ,  /* _INT_SCALAR     */
  AB6_INVARS_W90INIPRJ      ,  /* _INT_SCALAR     */
  AB6_INVARS_NSPPOL         ,  /* _INT_SCALAR     */
  AB6_INVARS_NFREQRE        ,  /* _INT_SCALAR     */
  AB6_INVARS_SUPERCELL      ,  /* _INT_ARRAY      */
  AB6_INVARS_BS_COULOMB_TERM,  /* _INT_SCALAR     */
  AB6_INVARS_AMU            ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PTGROUPMA      ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWNHATXC      ,  /* _INT_SCALAR     */
  AB6_INVARS_TPHYSEL        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_SYMREL         ,  /* _INT_ARRAY      */
  AB6_INVARS_PAWCPXOCC      ,  /* _INT_SCALAR     */
  AB6_INVARS_OPTFREQSUS     ,  /* _INT_SCALAR     */
  AB6_INVARS_WVL_FRMULT     ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PAWUSECP       ,  /* _INT_SCALAR     */
  AB6_INVARS_USERIB         ,  /* _INT_SCALAR     */
  AB6_INVARS_IDYSON         ,  /* _INT_SCALAR     */
  AB6_INVARS_SO_PSP         ,  /* _INT_ARRAY      */
  AB6_INVARS_PRTPOT         ,  /* _INT_SCALAR     */
  AB6_INVARS_ISTWFK         ,  /* _INT_ARRAY      */
  AB6_INVARS_OPTFORCES      ,  /* _INT_SCALAR     */
  AB6_INVARS_BRVLTT         ,  /* _INT_SCALAR     */
  AB6_INVARS_BOXCUTMIN      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NTIME          ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTVHA         ,  /* _INT_SCALAR     */
  AB6_INVARS_NPSP           ,  /* _INT_SCALAR     */
  AB6_INVARS_SPBROAD        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_SYMSIGMA       ,  /* _INT_SCALAR     */
  AB6_INVARS_GOPRECPRM      ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PRTVOL         ,  /* _INT_SCALAR     */
  AB6_INVARS_NQPTDM         ,  /* _INT_SCALAR     */
  AB6_INVARS_IPRCCH         ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTBLTZTRP     ,  /* _INT_SCALAR     */
  AB6_INVARS_ISECUR         ,  /* _INT_SCALAR     */
  AB6_INVARS_RHOQPMIX       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NATSPH         ,  /* _INT_SCALAR     */
  AB6_INVARS_GETBSCOUP      ,  /* _INT_SCALAR     */
  AB6_INVARS_GETDDK         ,  /* _INT_SCALAR     */
  AB6_INVARS_NWFSHIST       ,  /* _INT_SCALAR     */
  AB6_INVARS_FREQSPMIN      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_ECUT           ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PAWPRT_K       ,  /* _INT_SCALAR     */
  AB6_INVARS_NKPT           ,  /* _INT_SCALAR     */
  AB6_INVARS_SYMAFM         ,  /* _INT_ARRAY      */
  AB6_INVARS_PAWPRT_B       ,  /* _INT_SCALAR     */
  AB6_INVARS_DMFT_DC        ,  /* _INT_SCALAR     */
  AB6_INVARS_RPRIM_ORIG     ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_VDW_XC         ,  /* _INT_SCALAR     */
  AB6_INVARS_IEIG2RF        ,  /* _INT_SCALAR     */
  AB6_INVARS_RF2ATPOL       ,  /* _INT_ARRAY      */
  AB6_INVARS_GW_EET_INCLVKB ,  /* _INT_SCALAR     */
  AB6_INVARS_BS_EH_CUTOFF   ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_GW_EET_SCALE   ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_GETSUSCEP      ,  /* _INT_SCALAR     */
  AB6_INVARS_ICUTCOUL       ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTSPCUR       ,  /* _INT_SCALAR     */
  AB6_INVARS_NATOM          ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTDOS         ,  /* _INT_SCALAR     */
  AB6_INVARS_SPGORIG        ,  /* _INT_SCALAR     */
  AB6_INVARS_RECRCUT        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_NSPDEN         ,  /* _INT_SCALAR     */
  AB6_INVARS_SPINAT         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_SYMMORPHI      ,  /* _INT_SCALAR     */
  AB6_INVARS_FERMIE_NEST    ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_BS_HAYDOCK_TOL ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_NPSPALCH       ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTNABLA       ,  /* _INT_SCALAR     */
  AB6_INVARS_ZNUCL          ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_GETHAYDOCK     ,  /* _INT_SCALAR     */
  AB6_INVARS_PAWNZLM        ,  /* _INT_SCALAR     */
  AB6_INVARS_NFREQIM        ,  /* _INT_SCALAR     */
  AB6_INVARS_DMFT_NWLO      ,  /* _INT_SCALAR     */
  AB6_INVARS_UPAWU          ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PRTWANT        ,  /* _INT_SCALAR     */
  AB6_INVARS_MBAND          ,  /* _INT_SCALAR     */
  AB6_INVARS_GETSCR         ,  /* _INT_SCALAR     */
  AB6_INVARS_OPTDRIVER      ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTBBB         ,  /* _INT_SCALAR     */
  AB6_INVARS_LEXEXCH        ,  /* _INT_ARRAY      */
  AB6_INVARS_DTION          ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_GOPRECON       ,  /* _INT_SCALAR     */
  AB6_INVARS_PRTGKK         ,  /* _INT_SCALAR     */
  AB6_INVARS_QUADMOM        ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_PAWPRTVOL      ,  /* _INT_SCALAR     */
  AB6_INVARS_TOLVRS         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PAWNPHI        ,  /* _INT_SCALAR     */
  AB6_INVARS_RECNREC        ,  /* _INT_SCALAR     */
  AB6_INVARS_TSMEAR         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_TL_RADIUS      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_CD_FULL_GRID   ,  /* _INT_SCALAR     */
  AB6_INVARS_BS_CALCTYPE    ,  /* _INT_SCALAR     */
  AB6_INVARS_POSNSTEP       ,  /* _INT_SCALAR     */
  AB6_INVARS_QPTN           ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_IRDHAYDOCK     ,  /* _INT_SCALAR     */
  AB6_INVARS_FXCARTFACTOR   ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_RECNPATH       ,  /* _INT_SCALAR     */
  AB6_INVARS_ELPH2_IMAGDEN  ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_ESHIFT         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_GW_NPOLES      ,  /* _INT_SCALAR     */
  AB6_INVARS_ISTATR         ,  /* _INT_SCALAR     */
  AB6_INVARS_XC_DENPOS      ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_GETKSS         ,  /* _INT_SCALAR     */
  AB6_INVARS_NNSCLO         ,  /* _INT_SCALAR     */
  AB6_INVARS_LPAWU          ,  /* _INT_ARRAY      */
  AB6_INVARS_XCLEVEL        ,  /* _INT_SCALAR     */
  AB6_INVARS_NTYPPURE       ,  /* _INT_SCALAR     */
  AB6_INVARS_MQGRID         ,  /* _INT_SCALAR     */
  AB6_INVARS_TOLRFF         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PRTDENSPH      ,  /* _INT_SCALAR     */
  AB6_INVARS_SLABZEND       ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_PIMASS         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_DMFT_ITER      ,  /* _INT_SCALAR     */
  AB6_INVARS_NCONEQ         ,  /* _INT_SCALAR     */
  AB6_INVARS_PPMFRQ         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_CHKEXIT        ,  /* _INT_SCALAR     */
  AB6_INVARS_VACNUM         ,  /* _INT_SCALAR     */
  AB6_INVARS_TOLDFF         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_TOLDFE         ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_ORTALG         ,  /* _INT_SCALAR     */
  AB6_INVARS_VPRTRB         ,  /* _DOUBLE_ARRAY   */
  AB6_INVARS_FFT_OPT_LOB    ,  /* _INT_SCALAR     */
  AB6_INVARS_MKQMEM         ,  /* _INT_SCALAR     */
  AB6_INVARS_STMBIAS        ,  /* _DOUBLE_SCALAR  */
  AB6_INVARS_RF2ELFD        ,  /* _INT_SCALAR     */
  AB6_INVARS_NGFFT          ,  /* _INT_ARRAY      */
  AB6_INVARS_N_IDS
} Ab6InvarsIds;



Ab6InvarsTypes ab6_invars_get_type_from_id(Ab6InvarsIds id);
/**
 * AB6_INVARS_TYPE:
 * @A: an #Ab6InvarsIds id.
 *
 * Get the type of a given attribute of Dtset structure.
 *
 * Returns: a #Ab6InvarsTypes id.
 */
#define AB6_INVARS_TYPE(A) ab6_invars_get_type_from_id(A)
/**
 * AB6_INVARS_STR:
 * @A: an #Ab6InvarsIds id.
 *
 * Get a string corresponding to the attribute name.
 *
 * Returns: a string owned by ABINIT, do not free or modify it.
 */
#define AB6_INVARS_STR(A) #A

/**
 * Ab6Invars:
 *
 * An object to handle an array of ABINIT datasets, read from a file.
 */
typedef int Ab6Invars;

/**
 * ab6_invars_new_from_file:
 * @filename: a string, NULL terminated.
 *
 * Parse the given file using ABINIT routines and allocate a
 * dtsets array. This array must be deallocated after use with
 * ab6_invars_free().
 *
 * Returns: an #Ab6Invars object or NULL on failure.
 */
Ab6Invars* ab6_invars_new_from_file(const char *filename);
/**
 * ab6_invars_new_from_file_with_pseudo:
 * @filename: a string, NULL terminated.
 * @pspfiles: an array of strings, NULL terminated. Can be NULL.
 *
 * Parse the given file using ABINIT routines and allocate a
 * dtsets array. This array must be deallocated after use with
 * ab6_invars_free(). If pseudo files are provided with @pspfiles,
 * some further initialisations of dtset are permitted.
 *
 * Returns: an #Ab6Invars object or NULL on failure.
 */
Ab6Invars* ab6_invars_new_from_file_with_pseudo(const char *filename, const char **pspfiles);
/**
 * ab6_invars_new_from_string:
 * @string: a string, NULL terminated.
 *
 * Parse the given string using ABINIT routines and allocate a
 * dtsets array. This array must be deallocated after use with
 * ab6_invars_free().
 *
 * Returns: an #Ab6Invars object or NULL on failure.
 */
Ab6Invars* ab6_invars_new_from_string(const char *string);
/**
 * ab6_invars_free:
 * @ptr: the dataset array to handle.
 *
 * Clean all allocated memory from the data set allocation.
 */
void ab6_invars_free(Ab6Invars *ptr);

/**
 * ab6_invars_get_ndtset:
 * @ptr: the dataset array to handle.
 * @ndtset: a location to store the returned value.
 *
 * An array of datasets may contain more than one. Test it with this
 * routine. @ndtset will contains the number of allocated datasets (in
 * addition to the default one).
 *
 * Returns: #AB6_NO_ERROR if @ptr is valid and correctly parsed.
 */
Ab6Error ab6_invars_get_ndtset(Ab6Invars *ptr, int *ndtset);
/**
 * ab6_invars_get_integer:
 * @ptr: the dataset array to handle.
 * @id: an attribute id, see dtset_c.h.
 * @idtset: the number of the dtset to read, 0 is default value.
 * @value: a location to store the returned value.
 *
 * Use this method to get the value of an integer attribute. @idtset
 * must be in [0;n] where n is the returned value of
 * ab6_invars_get_ndtset(). If @id is unknown, return value is
 * 0. For real attributes, see ab6_invars_get_real().
 *
 * Returns: #AB6_NO_ERROR if values are correctly read.
 */
Ab6Error ab6_invars_get_integer(Ab6Invars *ptr, Ab6InvarsIds id,
                                int idtset, int *value);
/**
 * ab6_invars_get_real:
 * @ptr: the dataset array to handle.
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 * @value: a location to store the returned value.
 *
 * Use this method to get the value of a double attribute. @idtset
 * must be in [0;n] where n is the return value of
 * ab6_invars_get_ndtset(). If @id is unknown, return value is
 * undefined. For integer attributes, see ab6_invars_get_integer().
 *
 * Returns: #AB6_NO_ERROR if values are correctly read.
 */
Ab6Error ab6_invars_get_real(Ab6Invars *ptr, Ab6InvarsIds id,
                             int idtset, double *value);

/**
 * ab6_invars_get_shape:
 * @ptr: the dataset array to handle.
 * @n: a location to store the number of dimensions.
 * @dims: an array with 7 integers ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to poll the size of an array attribute. The
 * shape of the attribute is stored in @dims. Only the @n first values
 * of @dims are relevant.
 *
 * Returns: #AB6_NO_ERROR if values are correctly read.
 */
Ab6Error ab6_invars_get_shape(Ab6Invars *ptr, int *n, int dims[7],
			      Ab6InvarsIds id, int idtset);
/**
 * ab6_invars_get_integer_array:
 * @ptr: the dataset array to handle.
 * @values: an allocated array of @n values ;
 * @n: the size of the given array ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to read the values of an array. The array must
 * already be allocated. To know its size, use ab6_invars_get_shape().
 *
 * Returns: #AB6_NO_ERROR if values are correctly read.
 */
Ab6Error ab6_invars_get_integer_array(Ab6Invars *ptr, int *values, size_t n,
				      Ab6InvarsIds id, int idtset);
/**
 * ab6_invars_get_real_array:
 * @ptr: the dataset array to handle.
 * @values: an allocated array of @n values ;
 * @n: the size of the given array ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to read the values of an array. The array must
 * already be allocated. To know its size, use ab6_invars_get_shape().
 *
 * Returns: #AB6_NO_ERROR if values are correctly read.
 */
Ab6Error ab6_invars_get_real_array(Ab6Invars *ptr, double *values, size_t n,
				   Ab6InvarsIds id, int idtset);

#endif
