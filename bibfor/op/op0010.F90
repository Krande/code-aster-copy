! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine op0010()
!
! person_in_charge: patrick.massin at edf.fr
!
! aslint: disable=W1501
    implicit none
!
! ----------------------------------------------------------------------
!
! OPERATEUR PROPA_XFEM
!
! CALCUL DE LA FISSURE APRES PROPAGATION AU PAS DE TEMPS SUIVANT
!
! ----------------------------------------------------------------------
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/cncinv.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnscno.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infdbg.h"
#include "asterfort/infmaj.h"
#include "asterfort/ismali.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/pre_traitement.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xajuls.h"
#include "asterfort/xbaslo.h"
#include "asterfort/xenrch.h"
#include "asterfort/x_tmp_ligr.h"
#include "asterfort/xlenri.h"
#include "asterfort/xpraju.h"
#include "asterfort/xprdis.h"
#include "asterfort/xprdom.h"
#include "asterfort/xprgeo.h"
#include "asterfort/xprini.h"
#include "asterfort/xprls.h"
#include "asterfort/xprmil.h"
#include "asterfort/xprpls.h"
#include "asterfort/xprtor.h"
#include "asterfort/xprfastmarching.h"
#include "asterfort/xprupw_fmm.h"
#include "asterfort/xprvit.h"
    integer(kind=8) :: ifm, niv, ibid, ndim, iret, jcaraf, clsm, jma, jconx1, jconx2, nbma, i, ima
    integer(kind=8) :: j, nbrinit
    integer(kind=8) :: iadrma
    real(kind=8) :: lcmin
    character(len=8) :: k8bid, noma, nomo, fiss, fispre, method, fisini, ncrack
    character(len=8) :: ma_grill_pre
    character(len=16) :: k16bid, typdis, operation
    character(len=19) :: cnsvt, cnsvn, grlt, grln, cnslt, cnsln, cnsen, cnsenr
    character(len=19) :: noesom, isozro, cnxinv, cnsbl, cnsdis, cnslj
    character(len=19) :: vpoint, delta, mai
    character(len=24) :: lismae, lisnoe, vcn, grlr, vcnt, grlrt
    real(kind=8) :: meserr(3)
    character(len=8) :: test, msgout(2), typma
    aster_logical :: quad
!     MESSAGES
!
!     CRACK ADVANCEMENT
    real(kind=8) :: damax, dttot, vmax, rayon, dafiss, bmax
    character(len=24) :: vvit, vbeta, vgamma
    character(len=19) :: cnsbet, listp
    integer(kind=8) :: crack, jbeta, jgamma, jvit, nbval, nfiss
    integer(kind=8) :: numfis
!
!     LEVELSET AUXILIARY MESH
    character(len=8) :: unoma
    integer(kind=8) :: jlisno
    character(len=19) :: ucnslt, ucnsln, ugrlt, ugrln, ucnxin, disfr, nodtor
    character(len=19) :: eletor, liggrd, ligr_dnoma, ligr_noma
    aster_logical :: grille, locdom
!
!     DUMMY MESH
    character(len=8) :: dnoma
    character(len=19) :: dcnslt, dcnsln, dgrlt, dgrln, dcnxin
!
!     FIELD PROJECTION
    real(kind=8) :: radtor
    character(len=16) :: corres
    character(len=19) :: ndomp, edomg
!
!     TEST_MAIL
    real(kind=8) :: dist, distol
!
!     DOMAINE LOCALISATION
    integer(kind=8) :: nbno, jgltl, jglnl
    character(len=19) :: grltc, grlnc
    aster_logical :: ldpre
    real(kind=8) :: radimp, radlim
!
!     FRONT SUR LA GRILLE
    aster_logical :: goinop
    character(len=19) :: cnseg, cnseng, cnsljg
    character(len=24) :: lismag, lisnog
    real(kind=8), pointer :: gln(:) => null()
    real(kind=8), pointer :: glnp(:) => null()
    real(kind=8), pointer :: glt(:) => null()
    real(kind=8), pointer :: gltp(:) => null()
    character(len=8), pointer :: vfiss(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infdbg('XFEM', ifm, niv)
!
    damax = 0.d0
    dttot = 0.d0
    vmax = 0.d0
    rayon = 0.d0
    dafiss = 0.d0
    bmax = 0.d0
    radtor = 0.d0
    dist = 0.d0
    distol = 0.d0
    radimp = 0.d0
    radlim = 0.d0
!
! --- NOM DU CONCEPT FISSURE
!
    call getres(fiss, k16bid, k16bid)
!
!
! --- RETRIEVE THE NAME OF THE CRACK THAT MUST BE ELABORATED
!
    call getvid(' ', 'FISS_PROP', scal=fispre, nbret=ibid)
!
!   RECUPERATION DE LA METHODE DE REINITIALISATION A EMPLOYER
!
    call getvtx(' ', 'METHODE', scal=method, nbret=ibid)
!
!   VERIFICATION QUE L'ON TRAITE UNE FISSURE ET NON UNE INTERFACE
    call dismoi('TYPE_DISCONTINUITE', fispre, 'FISS_XFEM', repk=typdis)
!
!   OPERATION DEMANDEE
    call getvtx(' ', 'OPERATION', scal=operation, nbret=ibid)
!
    if (typdis .ne. 'FISSURE' .and. typdis .ne. 'COHESIF') then
        call utmess('F', 'XFEM2_1')
    end if
!
!
! --- NOM DU MODELE
!
    call getvid(' ', 'MODELE', scal=nomo, nbret=ibid)
!
!     SEARCH FOR THE CRACK THAT MUST BE PROPAGATED
    if (operation .ne. 'PROPA_COHESIF') then
        call dismoi('NB_FISS_XFEM', nomo, 'MODELE', repi=nfiss)
        if (nfiss .eq. 0) then
            call utmess('F', 'XFEM2_93', sk=nomo)
        end if
!
!     RETRIEVE THE NAME OF THE DATA STRUCTURE CONTAINING EACH CRACK
        call jeveuo(nomo//'.FISS', 'L', vk8=vfiss)
!
        numfis = 0
        do crack = 1, nfiss
!
            ncrack = vfiss(crack)
            if (ncrack .eq. fispre) numfis = crack
!
        end do
!
        if (numfis .eq. 0) then
            msgout(1) = fispre
            msgout(2) = nomo
            call utmess('F', 'XFEM2_89', nk=2, valk=msgout)
        end if
    else
        numfis = 1
        ncrack = fispre
    end if
!
! --- RETRIEVE THE NAME OF THE MESH THAT SHOULD BE USED AS AN AUXILIARY
!     GRID FOR THE EVALUATION OF THE LEVELSETS.
!
    call jeexin(fispre//'.GRI.MAILLA', ibid)
    if (ibid .eq. 0) then
!        NO AUXILIARY GRID USED
        grille = .false.
    else
!        AUXILIARY GRID USED
        grille = .true.
    end if
!
!     WRITE A WARNING IF THE CRACK HAS BEEN DEFINED GIVING DIRECTLY THE
!     TWO LEVEL SET FIELDS
    call jeexin(fispre//'.CHAMPS.LVS', ibid)
    if ((ibid .gt. 0) .and. (.not. grille)) then
        call jeveuo(fispre//'.CHAMPS.LVS', 'L', ibid)
        if (zl(ibid)) then
            call utmess('A', 'XFEM_69')
        end if
    end if

!   CHECK IF THE AUXILIARY GRID USED WHITH THE SIMPLEXE METHODE
    if (method .eq. 'SIMPLEXE' .and. grille) then
        call utmess('F', 'XFEM_27')
    end if

!
!     CHECK IF THE LOCALIZATION OF THE DOMAIN SHOULD BE ACTIVATED
    locdom = .false.
    call getvtx(' ', 'ZONE_MAJ', scal=k8bid, nbret=ibid)
    radimp = 0.d0
    if (k8bid(1:4) .eq. 'TORE') then
!        OK, THE LOCALIZATION MUST BE ACTIVATED
        locdom = .true.
!        CHECK IF THE USER HAS SPECIFIED THE RADIUS OF THE TORUS
        call getvr8(' ', 'RAYON_TORE', scal=radimp, nbret=ibid)
        if (ibid .eq. 0) then
!           THE USER HAS NOT SPECIFIED THE RADIUS OF THE TORUS
            radimp = -1.d0
        else
            radimp = radimp**2
        end if
    else
!        THE WHOLE GRID MUST BE USED
        locdom = .false.
    end if
!
!
! --- NOM DU MAILLAGE ATTACHE AU MODELE
!
    call jeveuo(nomo(1:8)//'.MODELE    .LGRF', 'L', iadrma)
    noma = zk8(iadrma)
!
! --- DIMENSION DU PROBLEME
!
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
    if ((ndim .lt. 2) .or. (ndim .gt. 3)) then
        call utmess('F', 'XFEM_18')
    end if
!
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
    call jeveuo(noma(1:8)//'.CONNEX', 'L', jconx1)
    call jeveuo(jexatr(noma(1:8)//'.CONNEX', 'LONCUM'), 'L', jconx2)
!
!   blindage : interdiction de traiter des mailles quadratiques
!
    mai = noma//'.TYPMAIL'
    call jeveuo(mai, 'L', jma)
    quad = .false.
    do ima = 1, nbma
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(jma-1+ima)), typma)
        if (.not. ismali(typma)) then
            quad = .true.
            exit
        end if
    end do
!
!!    if (quad) call utmess('F', 'XFEM_86')
!
! --- CONNECTIVITE INVERSEE
!
!
    cnxinv = '&&XPRREO.CNCINV'
    call cncinv(noma, [ibid], 0, 'V', cnxinv)
!  COMPATIBILITE ENTRE LA METHODE ET LE TYPE DE FISSURE
!
    if (typdis .eq. 'COHESIF') then
        ASSERT(operation .eq. 'DETECT_COHESIF' .or. operation .eq. 'PROPA_COHESIF')
    end if
    if (operation .eq. 'DETECT_COHESIF') then
        ASSERT(typdis .eq. 'COHESIF')
    end if
!
    if (operation .eq. 'PROPA_COHESIF') then
        ASSERT(typdis .eq. 'COHESIF')
    end if
!
!     RETRIEVE THE MAXIMUM ADVANCEMENT OF THE CRACK FRONT
    if (operation .ne. 'DETECT_COHESIF') then
        call getvr8(' ', 'DA_MAX', scal=damax, nbret=ibid)
    end if
!
!     RETRIEVE THE VALUE FOR THE "TEST_MAIL" PARAMETER
    call getvtx(' ', 'TEST_MAIL', scal=test, nbret=ibid)
!
!     ISSUE AN ALARM FOR THE USER
    if (test(1:3) .eq. 'OUI') then
        if (ndim .eq. 2) then
            call utmess('F', 'XFEM2_87')
        end if
    end if
!
!     RECUPERATION DES VITESSES DE PROPAGATION, DES ANGLES
!     DE BIFURCATION ET DE L'AVANCEE MAXIMALE DE LA FISSURE
!     A PROPAGER
!
    vvit = '&&OP0010.VVIT'
    vbeta = '&&OP0010.VBETA'
    vgamma = '&&OP0010.VGAMMA'
!
    call getvr8(' ', 'VITESSE', nbval=0, nbret=nbval)
!
    call wkvect(vbeta, 'V V R8', -nbval, jbeta)
    call wkvect(vgamma, 'V V R8', -nbval, jgamma)
    call wkvect(vvit, 'V V R8', -nbval, jvit)
!
    call getvr8(' ', 'ANGLE_BETA', nbval=-nbval, vect=zr(jbeta), nbret=ibid)
    call getvr8(' ', 'ANGLE_GAMMA', nbval=-nbval, vect=zr(jgamma), nbret=ibid)
    call getvr8(' ', 'VITESSE', nbval=-nbval, vect=zr(jvit), nbret=ibid)
!
    if (operation .ne. 'DETECT_COHESIF') then
        call getvr8(' ', 'DA_FISS', scal=dafiss, nbret=ibid)
!
!     RECUPERATION DU NOMBRE DE CYCLES
        call getvr8(' ', 'NB_CYCLES', scal=dttot, nbret=ibid)
    end if
!
! --- RECUPERATION DES LEVEL SETS ET GRADIENTS
!
    cnslt = '&&OP0010.CNSLT'
    cnsln = '&&OP0010.CNSLN'
    grlt = '&&OP0010.GRLT'
    grln = '&&OP0010.GRLN'
    call cnocns(fispre//'.LTNO', 'V', cnslt)
    call cnocns(fispre//'.LNNO', 'V', cnsln)
    call cnocns(fispre//'.GRLTNO', 'V', grlt)
    call cnocns(fispre//'.GRLNNO', 'V', grln)
!
! --- DUPLICATION DES GROUP_MA_ENRI ET GROUP_NO_ENRI
!
    lismae = fiss//'.GROUP_MA_ENRI'
    lisnoe = fiss//'.GROUP_NO_ENRI'
    call jedupo(fispre//'.GROUP_MA_ENRI', 'G', lismae, .false._1)
    call jedupo(fispre//'.GROUP_NO_ENRI', 'G', lisnoe, .false._1)
!
! --- DUPLICATION DE INFO ET MAILLAGE (LA NOUVELLE FISSURE RESTE
!     ATTACHE AU MAILLAGE INITIAL
!
    call jedupo(fispre//'.INFO', 'G', fiss//'.INFO', .false._1)
    call jedupo(fispre//'.MAILLAGE', 'G', fiss//'.MAILLAGE', .false._1)
!
!
    if (operation .eq. 'PROPA_COHESIF') then
        call jedupo(fispre//'.FONDFISS', 'G', fiss//'.FONDFISS', .false._1)
        call jedupo(fispre//'.BASEFOND', 'G', fiss//'.BASEFOND', .false._1)
        call jedupo(fispre//'.FONDMULT', 'G', fiss//'.FONDMULT', .false._1)
    end if
    if (operation .eq. 'PROPA_COHESIF' .or. operation .eq. 'DETECT_COHESIF') then
        call jedupo(fispre(1:8)//'.LISEQ     ', 'G', &
                    fiss(1:8)//'.LISEQ     ', .false._1)
    end if
!
! --- RECUPERATION DES CARACTERISTIQUES DU FOND DE FISSURE
!
    if (typdis .ne. 'COHESIF') then
        call jedupo(fispre//'.CARAFOND', 'G', fiss//'.CARAFOND', .false._1)
        call jeveuo(fiss//'.CARAFOND', 'L', jcaraf)
!
!
!   RETRIEVE THE RADIUS THAT MUST BE USED TO ASSESS THE LOCAL RESIDUAL
        call getvr8(' ', 'RAYON', scal=rayon, nbret=ibid)
    end if
!
!     SET THE DEFAULT VALUES FOR THE DOMAIN RESTRICTION FLAG
    ldpre = .false.
!
!-----------------------------------------------------------------------
!     RETRIEVE THE AUXILIARY GRID FOR THE LEVELSETS, IF THIS IS THE CASE
!-----------------------------------------------------------------------
    if (grille) then
!
!        RETRIEVE THE NAME OF THE MAILLAGE
        call jeveuo(fispre//'.GRI.MAILLA', 'L', iadrma)
        unoma = zk8(iadrma)
!
!        RETRIEVE THE LEVELSETS AND THEIR GRADIENTS ON THE AUXILIARY
!        GRID
        ucnslt = '&&OP0010.UCNSLT'
        ucnsln = '&&OP0010.UCNSLN'
        ugrlt = '&&OP0010.UGRLT'
        ugrln = '&&OP0010.UGRLN'
        call cnocns(fispre//'.GRI.LTNO', 'V', ucnslt)
        call cnocns(fispre//'.GRI.LNNO', 'V', ucnsln)
        call cnocns(fispre//'.GRI.GRLTNO', 'V', ugrlt)
        call cnocns(fispre//'.GRI.GRLNNO', 'V', ugrln)
!
!        CREATE THE INVERSE CONNECTIVITY
        ucnxin = '&&OP0010.UCNCINV'
        call cncinv(unoma, [ibid], 0, 'V', ucnxin)
!
!        CREATE A TEMPORARY JEVEUO OBJECT TO STORE THE "CONNECTION"
!        BETWEEN THE PHYSICAL AND AUXILIARY MESH USED IN THE PROJECTION
        corres = '&&OP0010.CORRES'
!
!        WRITE SOME INFORMATIONS ABOUT THE MODELS USED IN THE
!        PROPAGATION
        if (niv .ge. 0) then
            write (ifm, *) 'UNE GRILLE AUXILIAIRE EST UTILISEE POUR LA'// &
                ' PROPAGATION:'
            write (ifm, *) '   MODELE PHYSIQUE  : ', nomo
            write (ifm, *) '   MAILLAGE GRILLE AUXILIAIRE: ', unoma
        end if
!
    else
!
!        NO PROJECTION REQUIRED
        corres = ' '
!
!        WRITE SOME INFORMATIONS ABOUT THE MODELS USED IN THE
!        PROPAGATION
        if (niv .ge. 0) then
            write (ifm, *) 'LA PROPAGATION EST CALCULEE SUR LE MODELE '// &
                'DE LA STRUCTURE.'
            write (ifm, *) 'AUCUNE GRILLE AUXILIAIRE N''EST UTILISEE.'
            write (ifm, *) '   MODELE STRUCTURE: ', nomo
        end if
!
    end if
!
!-----------------------------------------------------------------------
!     CHECK FOR THE COHERENCE OF THE USE OF THE AUXILIARY GRID AND
!     DOMAIN LOCALISATION BETWEEN THE PREVIOUS AND THE ACTUAL STEP
!-----------------------------------------------------------------------
!
!     CHECK THE CONDITION ON THE VALUE OF DAMAX IF THE DOMAIN
!     LOCALISATION HAS BEEN REQUESTED (ONLY FOR 3D MESHES)
    if (locdom .and. (ndim .eq. 3)) then
        call jeveuo(vbeta, 'L', jbeta)
        call jelira(vbeta, 'LONMAX', j)
        bmax = 0.d0
        do i = 1, j
            if (abs(zr(jbeta-1+i)) .gt. bmax) bmax = abs(zr(jbeta-1+i))
        end do
!        THE CHECK IS MADE ONLY IF THE ANGLE IS GREATER THAN 3 DEGREES
!        AND LOWER OF 90 DEGREES
        if ((bmax .lt. 1.57d0) .and. (bmax .gt. 5.2d-2)) then
            if (dafiss .lt. (rayon/cos(bmax))) then
                meserr(1) = damax
                meserr(2) = dafiss
                meserr(3) = rayon/cos(bmax)*damax/dafiss
                call utmess('A', 'XFEM2_94', nr=3, valr=meserr)
            end if
        end if
    end if
!
    call jeexin(fispre//'.PRO.RAYON_TORE', ibid)
    if (ibid .ne. 0) then
        ldpre = .true.
    else
        ldpre = .false.
    end if
!
!     IF THE DOMAIN LOCALISATION HAS BEEN USED PREVIOUSLY, IT MUST BE
!     USED ALSO IN THIS STEP
    if (ldpre .and. (.not. locdom)) then
        call utmess('F', 'XFEM2_97')
    end if
!
!     IF AN AUXILIARY GRID IS USED IN THIS STEP, STORE ITS NAME FOR THE
!     NEW CRACK
    if (grille) then
        call jeveuo(fispre//'.GRI.MAILLA', 'L', ibid)
        ma_grill_pre = zk8(ibid)
        call wkvect(fiss//'.GRI.MAILLA', 'G V K8', 1, ibid)
        zk8(ibid) = ma_grill_pre
    end if
!
!-----------------------------------------------------------------------
!     SET THE CORRECT VALUE OF THE WORKING MODEL IN ORDER TO ASSESS
!     CORRECTLY THE CASE IN WHICH THE USER WANTS TO USE ONLY ONE MODEL
!     AND THE CASE IN WHICH HE OR SHE WANTS TO USE TWO DIFFERENT MODELS.
!     ALL THE FOLLOWING SUBROUTINES REFER TO A DUMMY MODEL AND ALL THE
!     VARIABLES REFERRED TO THIS MODEL BEGIN WITH THE LETTER "D".
!-----------------------------------------------------------------------
!
    if (grille) then
!
        dnoma = unoma
        dcnslt = ucnslt
        dcnsln = ucnsln
        dgrlt = ugrlt
        dgrln = ugrln
        dcnxin = ucnxin
!
    else
!
        dnoma = noma
        dcnslt = cnslt
        dcnsln = cnsln
        dgrlt = grlt
        dgrln = grln
        dcnxin = cnxinv
!
    end if
!
    ligr_dnoma = '&&OP0010.LIGDMA'
    call x_tmp_ligr(dnoma, ligr_dnoma)

!-----------------------------------------------------------------------
!     INITIALISE UPWIND FAST MARCHING
!-----------------------------------------------------------------------
    noesom = '&&OP0010.NOESOM'

    if (.not. grille) then
        vcn = '&&OP0010.VCN'
        grlr = '&&OP0010.GRLR'
    else
        vcn = unoma//'.GRLI'
        grlr = unoma//'.GRLR'
    end if
!
    if (grille) then
!        RETREIVE THE LENGTH OF THE SHORTEST EDGE IN THE GRID FROM THE
!        SD_GRILLE
        call jeveuo(unoma//'.GRLR', 'L', ibid)
        lcmin = zr(ibid)
    end if
!
    call xprini(dnoma, dcnxin, grille, noesom, vcn, grlr, lcmin, ndim)

!-----------------------------------------------------------------------
!     CALCUL DES POINTS DU FOND DE FISSURE SUR LA GRILLE
!     DANS LE CADRE DE L'UTILISATION D'UN FOND VIRTUEL
!-----------------------------------------------------------------------
!
    goinop = .false.
    if ((grille) .and. (ndim .eq. 3) .and. (method .ne. 'GEOMETRI')) then
        lismag = '&&OP0010.LISTE_MA_ENRICH'
        lisnog = '&&OP0010.LISTE_NO_ENRICH'
        goinop = .true.
        call xlenri(dnoma, fispre, goinop, lismag, lisnog)
!
        cnsljg = '&&OP0010.CNSLJG'
        cnseg = '&&OP0010.CNSEG'
        cnseng = '&&OP0010.CNSENG'
        call xenrch(dnoma, dcnslt, dcnsln, cnsljg, &
                    cnseg, cnseng, ndim, fispre, goinop, &
                    lismag, lisnog)
!
        call jedetr(cnsljg)
        call jedetr(cnseg)
        call jedetr(cnseng)
    end if
!
!-----------------------------------------------------------------------
!     CALCUL DES CHAM_NO_S DES VITESSES DE PROPAGATION
!-----------------------------------------------------------------------
!
!   Si locfom = false, initialisation de radtor à 0
    radtor = 0.d0
!
    if (locdom) then
        if (radimp .le. 0.d0) then
            radtor = (rayon+damax)**2
        end if
    end if
!
    if (niv .ge. 0) then
        write (ifm, *)
        write (ifm, *) 'OP0010-1) CALCUL DU CHAMP DE VITESSE AUX NOEUDS'
        write (ifm, 901)
    end if
!
    cnsvt = '&&OP0010.CNSVT'
    cnsvn = '&&OP0010.CNSVN'
    vpoint = '&&OP0010.VPOINT'
    cnsbet = '&&OP0010.CNSBET'
    listp = '&&OP0010.LISTP'
    disfr = '&&OP0010.DISFR'
    cnsbl = '&&OP0010.CNSBL'
    cnsdis = '&&OP0010.CNSDIS'
    delta = '&&OP0010.DELTA'
!
    if ((method .eq. 'GEOMETRI' .or. method .eq. 'SIMPLEXE') &
        .and. (operation .ne. 'PROPA_COHESIF')) then

        call pre_traitement(dnoma, fispre, ndim, vbeta, vgamma)
!
        call xprvit(dnoma, fispre, ndim, vvit, vbeta, &
                    lcmin, cnsvt, cnsvn, vpoint, cnsbl, &
                    cnsdis, disfr, cnsbet, listp, damax, &
                    locdom, radimp, radtor, delta, ucnslt, &
                    ucnsln)
    else
        call xprvit(dnoma, fispre, ndim, vvit, vbeta, &
                    lcmin, cnsvt, cnsvn, vpoint, cnsbl, &
                    cnsdis, disfr, cnsbet, listp, damax, &
                    locdom, radimp, radtor, delta, ucnslt, &
                    ucnsln)
    end if
!
!
!-----------------------------------------------------------------------
!     DOMAINS USED FOR THE RESTRICTION AND FOR THE PROJECTION
!-----------------------------------------------------------------------
!
    if (niv .ge. 0) then
        write (ifm, *)
        write (ifm, *) 'OP0010-2) DOMAINE DE CALCUL'
        write (ifm, 901)
    end if
!
    if (locdom) then
        vcnt = '&&OP0010.VCNT'
        grlrt = '&&OP0010.GRLRT'
    else
        vcnt = vcn
        grlrt = grlr
    end if
!
!     DEFINE THE PROJECTION DOMAINS FOR THE PHYSICAL AND LEVEL SET
!     MESHES (IF THE AUXILIARY GRID IS USED)
    if (grille) then
        ndomp = '&&OP0010.NDOMP'
        edomg = '&&OP0010.EDOMG'
        call xprdom(dnoma, dcnxin, disfr, noma, cnxinv, &
                    fispre, damax, ndomp, edomg, radtor)
    else if (operation .ne. 'DETECT_COHESIF') then
!        IF THE PROJECTION HAS NOT BEEN SELECTED, THE ESTIMATION OF THE
!        RADIUS OF THE TORUS DEFINING THE LOCAL DOMAIN TO BE USED FOR
!        THE LEVEL SET UPDATE CALCULATIONS MUST BE ESTIMATED HERE
        radtor = (rayon+damax)**2
    end if
!
!     RETREIVE THE RADIUS OF THE TORUS TO BE IMPOSED
    if (locdom) then
!        THE USER HAS NOT SPECIFIED THE RADIUS. THE RADIUS USED IN
!        THE PREVIOUS PROPAGATION SHOULD BE RETREIVED, IF ANY
        if ((radimp .lt. 0.d0) .and. ldpre) then
            call jeveuo(fispre//'.PRO.RAYON_TORE', 'L', ibid)
            radimp = zr(ibid)
        end if
    end if
!
!     DEFINE THE DOMAIN USED FOR THE LEVEL SET COMPUTATION (ONLY IF
!     THE LOCALISATION HAS BEEN SELECTED)
    nodtor = '&&OP0010.NODTOR'
    eletor = '&&OP0010.ELETOR'
    liggrd = '&&OP0010.LIGGRD'
!
    call xprtor(method, dnoma, dcnxin, fispre, &
                fiss, vcn, grlr, dcnsln, dgrln, &
                dcnslt, dgrlt, locdom, radtor, radimp, &
                cnsdis, disfr, cnsbl, nodtor, eletor, &
                liggrd, vcnt, grlrt)
!
!     CHECK IF THE RADIUS OF THE TORUS IS GREATER THAN THE CRITICAL
!     VALUE
    if (ldpre) then
        call jeveuo(fispre//'.PRO.RAYON_TORE', 'L', ibid)
!
!        CALCULATE THE CRITICAL VALUE OF THE RADIUS
        radlim = damax**2+zr(ibid)**2
!
        if (radlim .lt. radtor) then
            meserr(1) = sqrt(radtor)
            meserr(2) = sqrt(radlim)
            call utmess('A', 'XFEM2_88', nr=2, valr=meserr)
        end if
!
    end if
!
    if (locdom .and. (niv .ge. 0)) then
        write (ifm, *) '   LE DOMAINE DE CALCUL EST LOCALISE AUTOUR DU'//&
     &               ' FOND DE LA FISSURE:'
        write (ifm, *) '      RAYON DU TORE DE LOCALISATION = ',&
     &                sqrt(radtor)
!
        call jelira(nodtor, 'LONMAX', i)
        write (ifm, *) '      NOMBRE DE NOEUDS DU DOMAINE   = ', i
    end if
!
    if ((.not. locdom) .and. (niv .ge. 0)) then
        if (.not. grille) then
            write (ifm, *) '   LE DOMAINE DE CALCUL COINCIDE AVEC LE'//&
     &                  ' MODELE PHYSIQUE ', nomo
        else
            write (ifm, *) '   LE DOMAINE DE CALCUL COINCIDE AVEC LA'//&
     &                  ' GRILLE AUXILIAIRE ', unoma
        end if
    end if
!
!     MAKE SOME CHECKS
!
!     THE VALUE OF RAYON MUST BE GREATER THAN THE SHORTEST EDGE IN THE
!     MESH IN ORDER TO BE ABLE TO CALCULATE THE LOCAL RESIDUAL
    if (typdis .ne. 'COHESIF') then
        if (rayon .lt. lcmin) then
            meserr(1) = rayon
            meserr(2) = lcmin
            call utmess('F', 'XFEM2_64', nr=2, valr=meserr)
        end if
!
!       THE VALUE OF DAMAX SHOULD BE GREATER THAN THE SHORTEST EDGE IN THE
!       MESH. IF THIS IS NOT TRUE, THE MESH COULD FAIL TO CORRECTLY
!       REPRESENT THE LEVEL SETS. THIS IS NOT A FATAL ERROR AND A WARNING
!       MESSAGE IS ISSUED.
        if (lcmin .gt. damax) then
            meserr(1) = damax
            meserr(2) = lcmin
            call utmess('A', 'XFEM2_63', nr=2, valr=meserr)
        end if
    end if
!
!-----------------------------------------------------------------------
!     AJUSTEMENT DE VT
!-----------------------------------------------------------------------
!
    if (method .ne. 'GEOMETRI' .and. method .ne. 'SIMPLEXE') then
!
        if (niv .ge. 0) then
            write (ifm, *)
            write (ifm, *) 'OP0010-3) AJUSTEMENT DU CHAMP DES VITESSES VN'
            write (ifm, 903)
        end if
!
        call xpraju(dnoma, fiss, dcnslt, cnsvt, cnsvn, &
                    dttot, vmax)
!
    end if
!
!-----------------------------------------------------------------------
!     PROPAGATION DES LEVEL SETS
!-----------------------------------------------------------------------
    if (niv .ge. 0) then
        write (ifm, *)
        if (method .eq. 'GEOMETRI') then
            write (ifm, *) 'OP0010-3) MISE A JOUR DES LEVEL SETS'
        else
            write (ifm, *) 'OP0010-4) PROPAGATION DES LEVEL SETS'
        end if
        write (ifm, 904)
!
!        WRITE SOME INFORMATIONS
        write (ifm, *) '   AVANCEE MAXIMALE DU FOND DE FISSURE    = '&
     &               , dafiss
        write (ifm, *) '   NOMBRE DE CYCLES DE FATIGUE            = '&
     &               , dttot
    end if
!
    if (method .eq. 'GEOMETRI') then
        write (ifm, *) '   '
        write (ifm, *) '   UTILISATION DE LA METHODE GEOMETRIQUE.'
        call xprgeo(dnoma, dcnsln, dcnslt, dgrln, dgrlt, &
                    vpoint, cnsbl, dttot, nodtor, liggrd, &
                    cnsbet, listp, operation)
        goto 100
    end if
!
    if (method .eq. 'SIMPLEXE') then
        write (ifm, *) '   '
        write (ifm, *) '   UTILISATION DE LA METHODE SIMPLEXE.'
        call xprgeo(dnoma, dcnsln, dcnslt, dgrln, dgrlt, &
                    vpoint, cnsbl, dttot, nodtor, liggrd, &
                    cnsbet, listp, operation)
    else
        call xprls(dnoma, dcnsln, dcnslt, dgrln, dgrlt, &
                   cnsvn, cnsvt, cnsbl, dttot, nodtor, &
                   eletor, liggrd, delta)
    end if
!
    call jedetr(cnsvt)
    call jedetr(cnsvn)

!
!-----------------------------------------------------------------------
!     REINITIALISATION DE LSN
!-----------------------------------------------------------------------
!
    if (niv .ge. 0) then
        write (ifm, *)
        write (ifm, *) 'OP0010-5) REINITIALISATION DE LSN'
        write (ifm, 905)
    end if
!
    isozro = '&&OP0010.ISOZRO'

    if (method .eq. 'UPWIND') then
        nbrinit = 1
        call xprupw_fmm('REINITLN', dnoma, vcnt, grlrt, &
                        noesom, lcmin, dcnsln, dgrln, dcnslt, &
                        dgrlt, isozro, nodtor, eletor, liggrd, &
                        vpoint, cnsbl, dttot, cnsbet, listp, nbrinit)
    end if

    if (method .eq. 'SIMPLEXE') then
        call xprfastmarching('REINITLN', dnoma, cnxinv, &
                             noesom, lcmin, dcnsln, dgrln, dcnslt, &
                             dgrlt, isozro, nodtor, eletor, liggrd, &
                             vpoint, cnsbl, cnsbet, listp)
    end if

    call jedetr(isozro)
!
!-----------------------------------------------------------------------
!     REINITIALISATION DE  LST
!-----------------------------------------------------------------------
!
    if (niv .ge. 0) then
        write (ifm, *)
        write (ifm, *) 'OP0010-7) REINITIALISATION DE LST'
        write (ifm, 907)
    end if

!   On réinitialise deux fois pour redresser l'iso zéro
    if (method .eq. 'UPWIND') then
        do nbrinit = 1, 2
            call xprupw_fmm('REINITLT', dnoma, vcnt, grlrt, &
                            noesom, lcmin, dcnsln, dgrln, dcnslt, &
                            dgrlt, isozro, nodtor, eletor, liggrd, &
                            vpoint, cnsbl, dttot, cnsbet, listp, nbrinit)
            call jedetr(isozro)
        end do
    end if

    if (method .eq. 'SIMPLEXE') then
        call xprfastmarching('REINITLT', dnoma, cnxinv, &
                             noesom, lcmin, dcnsln, dgrln, dcnslt, &
                             dgrlt, isozro, nodtor, eletor, liggrd, &
                             vpoint, cnsbl, cnsbet, listp)
        call jedetr(isozro)
    end if

    call jedetr(cnsbl)
!------------------------------------------------------------------!

100 continue
    call jedetr(vvit)
    call jedetr(vbeta)
    call jedetr(noesom)
    if (method .eq. 'UPWIND') then
        if (.not. grille) then
            call jedetr(vcn)
            call jedetr(grlr)
        end if
        if (locdom) then
            call jedetr(vcnt)
            call jedetr(grlrt)
        end if
    end if
    call jedetr(cnsdis)
    call jedetr(nodtor)
    call jedetr(eletor)
    call jedetr(vpoint)
    call jedetr(cnsbet)
    call jedetr(listp)
!
!-----------------------------------------------------------------------
!     THE NEW VALUES OF THE LEVELSETS FOR THE AUXILIARY MESH (TWO GRIDS
!     CASE ONLY) ARE STORED. AFTER THAT THESE VALUES MUST BE PROJECTED
!     TO THE PHYSICAL MESH FOR THE FRACTURE MECHANICS COMPUTATIONS.
!-----------------------------------------------------------------------
!
    if (grille) then
!
        ligr_noma = '&&OP0010.LIGMAI'
        call x_tmp_ligr(noma, ligr_noma)
!
!       CREATE THE CHAMP_NO WITH THE NEW VALUES OF THE LEVELSETS AND
!       THEIR GRADIENTS. THE EXISTING CHAMP_NO ARE AUTOMATICALLY
!       DESTROYED BY THE SUBROUTINE "CNSCNO"
        call cnscno(dcnslt, ' ', 'OUI', 'G', fiss//'.GRI.LTNO', &
                    'F', ibid)
        call cnscno(dcnsln, ' ', 'OUI', 'G', fiss//'.GRI.LNNO', &
                    'F', ibid)
        call cnscno(dgrlt, ' ', 'OUI', 'G', fiss//'.GRI.GRLTNO', &
                    'F', ibid)
        call cnscno(dgrln, ' ', 'OUI', 'G', fiss//'.GRI.GRLNNO', &
                    'F', ibid)
!
!       PROJECT THE LEVEL SETS
        call xprpls(ligr_noma, dnoma, dcnsln, dcnslt, noma, &
                    cnsln, cnslt, grln, grlt, corres, &
                    ndim, ndomp, edomg)
!
!       STORE THE LIST OF THE NODES OF THE STRUCTURAL MESH WHERE THE
!       PROJECTION HAS BEEN CARRIED OUT
        call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=j)
        call wkvect(fiss//'.PRO.NOEUD_PROJ', 'G V L', j, iret)
!
        call jeveuo(ndomp, 'L', ibid)
        call jelira(ndomp, 'LONMAX', j)
!
        do i = 1, j
            zl(iret-1+zi(ibid-1+i)) = .true.
        end do
!
    end if
!
    call jedetr(disfr)
!
!     NOW I CAN WORK ON THE PHYSICAL MESH
!
!-----------------------------------------------------------------------
!     REAJUSTEMENT DES LEVEL SETS TROP PROCHES DE 0
!-----------------------------------------------------------------------
    if (niv .ge. 0) then
        write (ifm, *)
        if (method .eq. 'GEOMETRI') then
            write (ifm, *) 'OP0010-4) ENRICHISSEMENT DE LA SD FISS_XFEM'
        else
            write (ifm, *) 'OP0010-8) ENRICHISSEMENT DE LA SD FISS_XFEM'
        end if
        write (ifm, 908)
    end if
!
!   ON RAJOUTE UN CRITERE LST PLUS LACHE POUR EVITER DES OSCILLATIONS
!   NUMERIQUE LORS DE LA PROPAGATION / EN THEORIE CE N EST PAS BIEN DE
!   DE RETOUCHER LA GEOMETRIE A LA VOLEE
!
    call xajuls(noma, nbma, cnslt, cnsln, jconx1, &
                jconx2, clsm, typdis, critlst=1.d-3)
!
    if (niv .ge. 0) then
        write (ifm, *) 'NOMBRE DE LEVEL SET REAJUSTEES APRES CONTROLE:', &
            clsm
    end if
!
!-----------------------------------------------------------------------
!     EXTENSION DES LEVEL SETS AUX NOEUDS MILIEUX
!-----------------------------------------------------------------------
!
    call xprmil(noma, cnslt, cnsln)
!
!     IF THE DOMAINE LOCALISATION HAS BEEN USED ON THE PHYSICAL MODEL,
!     THE VALUES OF THE LEVEL SET GRADIENTS OUTSIDE THE DOMAINE MUST
!     BE COPIED FROM THE ORIGINAL VALUES (NOW THEY ARE EQUAL TO ZERO!)
    if ((.not. grille) .and. locdom) then
!        RETRIEVE THE LIST OF THE NODES USED IN THE LAST LOCALIZATION
        call jeveuo(fiss//'.PRO.NOEUD_TORE', 'E', jlisno)
        call jelira(fiss//'.PRO.NOEUD_TORE', 'LONMAX', nbno)
!
        grltc = '&&OP0010.GRLTC'
        grlnc = '&&OP0010.GRLNC'
!
        call cnocns(fispre//'.GRLTNO', 'V', grltc)
        call cnocns(fispre//'.GRLNNO', 'V', grlnc)
!
        call jeveuo(grltc//'.CNSV', 'L', vr=gltp)
        call jeveuo(grlnc//'.CNSV', 'L', vr=glnp)
!
        call jeveuo(grlt//'.CNSV', 'E', vr=glt)
        call jeveuo(grlt//'.CNSL', 'E', jgltl)
        call jeveuo(grln//'.CNSV', 'E', vr=gln)
        call jeveuo(grln//'.CNSL', 'E', jglnl)
!
        do i = 1, nbno
!
            if (.not. zl(jlisno-1+i)) then
                glt(ndim*(i-1)+1) = gltp(ndim*(i-1)+1)
                zl(jgltl-1+ndim*(i-1)+1) = .true.
                glt(ndim*(i-1)+2) = gltp(ndim*(i-1)+2)
                zl(jgltl-1+ndim*(i-1)+2) = .true.
!
                gln(ndim*(i-1)+1) = glnp(ndim*(i-1)+1)
                zl(jglnl-1+ndim*(i-1)+1) = .true.
                gln(ndim*(i-1)+2) = glnp(ndim*(i-1)+2)
                zl(jglnl-1+ndim*(i-1)+2) = .true.
!
                if (ndim .eq. 3) then
                    glt(ndim*(i-1)+3) = gltp(ndim*(i-1)+3)
                    zl(jgltl-1+ndim*(i-1)+3) = .true.
!
                    gln(ndim*(i-1)+3) = glnp(ndim*(i-1)+3)
                    zl(jglnl-1+ndim*(i-1)+3) = .true.
                end if
!
            end if
!
        end do
!
        call jedetr(grltc)
        call jedetr(grlnc)
!
    end if
!
    call cnscno(cnslt, ' ', 'NON', 'G', fiss//'.LTNO', &
                'F', ibid)
    call cnscno(cnsln, ' ', 'NON', 'G', fiss//'.LNNO', &
                'F', ibid)
    call cnscno(grlt, ' ', 'NON', 'G', fiss//'.GRLTNO', &
                'F', ibid)
    call cnscno(grln, ' ', 'NON', 'G', fiss//'.GRLNNO', &
                'F', ibid)
!
!     IF THE DOMAIN LOCALISATION HAS NOT BEEN USED, THE BOOLEAN LIST
!     OF THE NODES IN THE TORE MUST BE DESTROYED
    if (.not. locdom) then
        call jedetr(fiss//'.PRO.NOEUD_TORE')
    end if
!
!     IF THE DOMAIN LOCALISATION HAS BEEN USED, THE RADIUS OF THE TORUS
!     USED IN THE LOCALISATION MUST BE STORED
    if (locdom) then
        call wkvect(fiss//'.PRO.RAYON_TORE', 'G V R', 1, ibid)
!        VALUE OF THE RADIUS USED IN THE ACTUAL PROPAGATION
        zr(ibid) = radtor
    end if
!
    call jedetr(delta)
!----------------------------------------------------------------------+
!                 FIN DE LA PARTIE PROPAGATION :                       |
!                 ----------------------------                         |
!    LA SD FISS_XFEM EST ENRICHIE COMME DANS OP0041 : DEFI_FISS_XFEM   |
!   ( TOUTE MODIF. AFFECTANT OP0041 DOIT ETRE REPERCUTEE PLUS BAS,     |
!     EXCEPTE L'APPEL A SDCONX )                                       |
!----------------------------------------------------------------------+
!
!-----------------------------------------------------------------------
!     CALCUL DE L'ENRICHISSEMENT ET DES POINTS DU FOND DE FISSURE
!-----------------------------------------------------------------------
!
    cnslj = '&&OP0010.CNSLJ'
    cnsen = '&&OP0010.CNSEN'
    cnsenr = '&&OP0010.CNSENR'
    goinop = .false.
    call xenrch(noma, cnslt, cnsln, cnslj, &
                cnsen, cnsenr, ndim, fiss, goinop, &
                lismae, lisnoe, operation_opt=operation)
!
!    call cnscno(cnsenr, ' ', 'NON', 'G', fiss//'.STNOR',&
!                'F', ibid)
    call cnscno(cnsen, ' ', 'NON', 'G', fiss//'.STNO', &
                'F', ibid)
!
!-----------------------------------------------------------------------
!     CALCUL DE LA BASE LOCALE AU FOND DE FISSURE
!-----------------------------------------------------------------------
!
    call xbaslo(noma, fiss, grlt, grln, ndim)
!
!-----------------------------------------------------------------------
!     ELABORATE THE OPTION "TEST_MAIL"
!-----------------------------------------------------------------------
!
!     CHECK THE MESH, IF THIS IS THE CASE
    if (test(1:3) .eq. 'OUI') then
!        RETREIVE THE DISTANCE BETWEEN THE PROPAGATED CRACK AND THE
!        INITIAL FRONT
        call getvr8(' ', 'DISTANCE', scal=dist, nbret=ibid)
!        RETREIVE THE VALUE OF THE TOLERANCE
        call getvr8(' ', 'TOLERANCE', scal=distol, nbret=ibid)
!        RETREIVE THE INITIAL CRACK
        call getvid(' ', 'FISS_INITIALE', scal=fisini, nbret=ibid)
!        CHECK THE CRACK FRONT
        call xprdis(fisini, fiss, dist, distol, lcmin)
!
    end if
!
!-----------------------------------------------------------------------
!     FIN
!-----------------------------------------------------------------------
    call jedetr(cnxinv)
    call detrsd('LIGREL', liggrd)
    call detrsd('LIGREL', ligr_dnoma)
    if (grille) then
        call detrsd('LIGREL', ligr_noma)
    end if
!
901 format(10x, 37('-'))
903 format(10x, 35('-'))
904 format(10x, 26('-'))
905 format(10x, 23('-'))
907 format(10x, 23('-'))
908 format(10x, 33('-'))
!
    call jedema()
end subroutine
