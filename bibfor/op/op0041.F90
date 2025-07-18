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

subroutine op0041()
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
!
! ----------------------------------------------------------------------
!
! OPERATEUR DEFI_FISS_XFEM
!
! INITIALISATION DES CHAMPS NECESSAIRES A XFEM
!        - LEVEL-SETS
!        - GRADIENTS DES LEVEL-SETS
!        - MAILLES ENRICHIES DE LA ZONE FISSURE
!        - POINTS DU FOND DE FISSURE
!
! ----------------------------------------------------------------------
!
! N.B.: TOUTE MODIFICATION EFFECTUE APRES LE CALCUL DES LEVEL SETS&
!        LEURS GRADIENT DOIT ETRE REPERCUTEE DANS OP0010 : PROPA_XFEM
!        (MIS A PART L'APPEL A SDCONX A LA FIN)
!
! CONCEPT SORTANT: FISS DE TYPE FISS_XFEM
!
!     CONDENU DE LA SD FISS_XFEM
!         FISS//'.GROUP_MA_ENRI'
!         FISS//'.GROUP_NO_ENRI'
!         FISS//'.LTNO'
!         FISS//'.LNNO'
!         FISS//'.GRLTNO'
!         FISS//'.GRLNNO'
!         FISS//'.MAILFISS.HEAV'
!         FISS//'.MAILFISS.CTIP'
!         FISS//'.MAILFISS.HECT'
!         FISS//'.MAILFISS.MAFOND'
!         FISS//'.FONDFISS'
!         FISS//'.FONDMULT'
!         FISS//'.BASLOC'
!         FISS//'.BASFOND'
!         FISS//'.INFO'
!         FISS//'.MAILLAGE'
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/cnsprj.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprsd.h"
#include "asterfort/infdbg.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/pj2dco.h"
#include "asterfort/pj3dco.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xbaslo.h"
#include "asterfort/xenrch.h"
#include "asterfort/xgrals.h"
#include "asterfort/xinils.h"
#include "asterfort/xinlsj.h"
#include "asterfort/xlenri.h"
#include "asterfort/xlorie.h"
    integer(kind=8) :: ifm, niv, ibid, mxval, iret
    integer(kind=8) :: me1, me2, me3, me4
    integer(kind=8) :: ndim, jinfo, jmod, jma
    real(kind=8) :: noeud(3), vect1(3), vect2(3), a, b, r, dmax
    character(len=8) :: fiss, nfonf, nfong, mafis, fonfis, noma, meth
    character(len=8) :: maiaux
    character(len=8) :: cote, ncham, chadis, kbid
    character(len=16) :: k16bid, geofis, typdis, corres
    character(len=19) :: cnslt, cnsln, grlt, grln, cnsen, cnsenr, cnslj
    character(len=19) :: cnsltg, cnslng, grltg, grlng
    character(len=19) :: ltno, lnno, grltno, grlnno, stno, info, ltnofa
    character(len=19) :: lnnofa, grltfa, grlnfa
    character(len=24) :: lismae, lisnoe
    aster_logical :: grille, l_dmax, goinop
    character(len=8) :: fisgri
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infdbg('XFEM', ifm, niv)
!
! --- NOM DU CONCEPT FISSURE
!
    call getres(fiss, k16bid, k16bid)
!
! --- NOM DU MAILLAGE
!
    call getvid(' ', 'MAILLAGE', scal=noma, nbret=ibid)
    call wkvect(fiss//'.MAILLAGE', 'G V K8', 1, jma)
    zk8(jma-1+1) = noma
!
! --- DIMENSION DU PROBLEME
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
!
!     POUR DIFFERENCIER OP0041 et OP0010
    goinop = .false.
!
!     CHECK IF THE USER WANTS TO USE AN AUXILIARY GRID
    call getvid(' ', 'MAILLAGE_GRILLE', scal=maiaux, nbret=iret)
    if (iret .gt. 0) then
!        YES
        grille = .true.
        write (ifm, 900) maiaux
!
!        CHECK IF THE MESH IS A SD_GRILLE
        call jeexin(maiaux//'.GRLI', ibid)
        if (ibid .eq. 0) then
            call utmess('F', 'XFEM2_95', sk=maiaux)
        end if
!
!        THE GRID AND MODEL DIMENSIONS MUST BE THE SAME
        call dismoi('DIM_GEOM', maiaux, 'MAILLAGE', repi=ibid)
        if (ibid .ne. ndim) then
            call utmess('F', 'XFEM2_58')
        end if
!
!        STORE THE AUXILIARY GRID ON WHICH THE CRACK WILL BE
!        DEFINED
        call wkvect(fiss(1:8)//'.GRI.MAILLA', 'G V K8', 1, jmod)
        zk8(jmod-1+1) = maiaux
    else
        grille = .false.
    end if
!
!     CHECK IF THE USER HAS GIVEN THE CRACK FROM WHICH THE GRID MUST
!     BE COPIED
    call getvid(' ', 'FISS_GRILLE', scal=fisgri, nbret=iret)
    if (iret .gt. 0) then
!        YES, THE GRID INFOS ARE DUPLICATED FOR THE NEW CRACK.
!        CHECK IF A GRID IS ASSOCIATED TO THE GIVEN CRACK.
        call jeexin(fisgri//'.GRI.MAILLA', ibid)
        if (ibid .eq. 0) then
            call utmess('F', 'XFEM_68')
        end if
!
        call jedupo(fisgri//'.GRI.MAILLA', 'G', fiss(1:8)//'.GRI.MAILLA', .false._1)
        call copisd('CHAMP', 'G', fisgri//'.GRI.LNNO', fiss(1:8)//'.GRI.LNNO')
        call copisd('CHAMP', 'G', fisgri//'.GRI.GRLNNO', fiss(1:8)//'.GRI.GRLNNO')
!
        call jeexin(fisgri//'.GRI.LTNO  .REFE', ibid)
        if (ibid .gt. 0) then
            call copisd('CHAMP', 'G', fisgri//'.GRI.LTNO', fiss(1:8)//'.GRI.LTNO')
            call copisd('CHAMP', 'G', fisgri//'.GRI.GRLTNO', fiss(1:8)//'.GRI.GRLTNO')
        end if
!
        call jeexin(fisgri//'.PRO.RAYON_TORE', ibid)
        if (ibid .gt. 0) then
            call jedupo(fisgri//'.PRO.RAYON_TORE', 'G', fiss(1:8)//'.PRO.RAYON_TORE', .false._1)
            call jedupo(fisgri//'.PRO.NOEUD_TORE', 'G', fiss(1:8)//'.PRO.NOEUD_TORE', .false._1)
        end if
!
        grille = .false.
        write (ifm, *) '  LA GRILLE AUXILIAIRE UTILISEE POUR LA FISSURE ',&
     &                fisgri
        write (ifm, *) '  EST UTILISEE AUSSI POUR LA NOUVELLE FISSURE ',&
     &                fiss
        write (ifm, *) '  ET LES LEVEL SETS DEFINIES SUR CETTE GRILLE ONT'&
     &               //' ETE PRESERVEES.'
    end if
!
! --- OJBET INFORMATIONS : TYPE_DISCONT, CHAM_DISCONT ET TYPE_FOND
    info = fiss//'.INFO'
    call wkvect(info, 'G V K16', 3, jinfo)
!     TYPE DE DISCONTINUITE : FISSURE OU INTERFACE
    call getvtx(' ', 'TYPE_DISCONTINUITE', scal=typdis, nbret=ibid)
!     CHAMP DISCONTINU : DEPLACEMENTS OU CONTRAINTES
    call getvtx(' ', 'CHAM_DISCONTINUITE', scal=chadis, nbret=ibid)
    zk16(jinfo-1+1) = typdis
    zk16(jinfo-1+2) = chadis
    zk16(jinfo-1+3) = '      '
!
! --- MOT-CLEFS DEFINITION FISSURE
!
    call getvid('DEFI_FISS', 'FONC_LT', iocc=1, scal=nfonf, nbret=ibid)
    call getvid('DEFI_FISS', 'FONC_LN', iocc=1, scal=nfong, nbret=me1)
    if (me1 .eq. 1 .and. ibid .eq. 0 .and. typdis .eq. 'FISSURE') then
        call utmess('F', 'XFEM_24', sk='FONC_LT')
    end if
    if (me1 .eq. 1 .and. ibid .eq. 1 .and. typdis .eq. 'INTERFACE') then
        call utmess('A', 'XFEM_25', sk='FONC_LT')
    end if
!
    call getvtx('DEFI_FISS', 'GROUP_MA_FISS', iocc=1, scal=mafis, nbret=me2)
    call getvtx('DEFI_FISS', 'GROUP_MA_FOND', iocc=1, scal=fonfis, nbret=ibid)
    if ((me2 .eq. 1 .or. ibid .eq. 1) .and. typdis .eq. 'COHESIF') then
        ASSERT(.false.)
    end if
    if (me2 .eq. 1 .and. ibid .eq. 0 .and. typdis .eq. 'FISSURE') then
        call utmess('F', 'XFEM_24', sk='GROUP_MA_FOND')
    end if
    if (me2 .eq. 1 .and. ibid .eq. 1 .and. typdis .eq. 'INTERFACE') then
        call utmess('A', 'XFEM_25', sk='GROUP_MA_FOND')
    end if
!
    mxval = 0
    call getvtx('DEFI_FISS', 'FORM_FISS', iocc=1, nbval=mxval, vect=geofis, &
                nbret=me3)
!
    call getvid('DEFI_FISS', 'CHAM_NO_LSN', iocc=1, scal=ncham, nbret=me4)
    call getvid('DEFI_FISS', 'CHAM_NO_LST', iocc=1, scal=ncham, nbret=ibid)
    if ((me4 .eq. 1 .or. ibid .eq. 1) .and. typdis .eq. 'COHESIF') then
        ASSERT(.false.)
    end if
    if (me4 .eq. 1 .and. ibid .eq. 0 .and. typdis .eq. 'FISSURE') then
        call utmess('F', 'XFEM_24', sk='CHAM_NO_LST')
    end if
    if (me4 .eq. 1 .and. ibid .eq. 1 .and. typdis .eq. 'INTERFACE') then
        call utmess('A', 'XFEM_25', sk='CHAM_NO_LST')
    end if
!
    if (me3 .eq. -1) then
        call getvtx('DEFI_FISS', 'FORM_FISS', iocc=1, scal=geofis, nbret=me3)
        call getvr8('DEFI_FISS', 'RAYON_CONGE', iocc=1, scal=r, nbret=ibid)
        if (geofis .eq. 'ELLIPSE' .or. geofis .eq. 'RECTANGLE' .or. geofis .eq. 'CYLINDRE') then
            call getvr8('DEFI_FISS', 'DEMI_GRAND_AXE', iocc=1, scal=a, nbret=ibid)
            call getvr8('DEFI_FISS', 'DEMI_PETIT_AXE', iocc=1, scal=b, nbret=ibid)
            call getvr8('DEFI_FISS', 'CENTRE', iocc=1, nbval=3, vect=noeud, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'VECT_X', iocc=1, nbval=3, vect=vect1, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'VECT_Y', iocc=1, nbval=3, vect=vect2, &
                        nbret=ibid)
            call getvtx('DEFI_FISS', 'COTE_FISS', iocc=1, scal=cote, nbret=ibid)
        else if (geofis .eq. 'DEMI_PLAN') then
            call getvr8('DEFI_FISS', 'PFON', iocc=1, nbval=3, vect=noeud, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'NORMALE', iocc=1, nbval=3, vect=vect1, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'DTAN', iocc=1, nbval=3, vect=vect2, &
                        nbret=ibid)
        else if (geofis .eq. 'SEGMENT') then
            call getvr8('DEFI_FISS', 'PFON_ORIG', iocc=1, nbval=3, vect=vect1, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'PFON_EXTR', iocc=1, nbval=3, vect=vect2, &
                        nbret=ibid)
        else if (geofis .eq. 'DEMI_DROITE') then
            call getvr8('DEFI_FISS', 'PFON', iocc=1, nbval=3, vect=noeud, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'DTAN', iocc=1, nbval=3, vect=vect1, &
                        nbret=ibid)
        else if (geofis .eq. 'DROITE') then
            call getvr8('DEFI_FISS', 'POINT', iocc=1, nbval=3, vect=noeud, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'DTAN', iocc=1, nbval=3, vect=vect1, &
                        nbret=ibid)
        else if (geofis .eq. 'ENTAILLE') then
            call getvr8('DEFI_FISS', 'DEMI_LONGUEUR', iocc=1, scal=a, nbret=ibid)
            call getvr8('DEFI_FISS', 'CENTRE', iocc=1, nbval=3, vect=noeud, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'VECT_X', iocc=1, nbval=3, vect=vect1, &
                        nbret=ibid)
            call getvr8('DEFI_FISS', 'VECT_Y', iocc=1, nbval=3, vect=vect2, &
                        nbret=ibid)
        else
            ASSERT(.false.)
        end if
!
    end if
!
! --- STOCKAGE DES DONNEES ORIENTATION FOND DE FISSURE
!     ON ENRICHI LA SD FISS_XFEM DE
!         FISS//'.CARAFOND'
!
    if (typdis .eq. 'FISSURE') then
        call xlorie(fiss)
    end if
!
! --- RECUPERATION DES GROUP_MA_ENRI ET GROUP_NO_ENRI
!     ON ENRICHI LA SD FISS_XFEM DE
!         FISS//'.GROUP_MA_ENRI'
!         FISS//'.GROUP_NO_ENRI'
!
    lismae = '&&OP0041.LISTE_MA_ENRICH'
    lisnoe = '&&OP0041.LISTE_NO_ENRICH'
    call xlenri(noma, fiss, goinop, lismae, lisnoe)
!
!-----------------------------------------------------------------------
!     CALCUL DES LEVEL-SETS
!
!     ON ENRICHI LA SD FISS_XFEM DE
!         FISS//'.LTNO'
!         FISS//'.LNNO'
!
!-----------------------------------------------------------------------
!
    cnslt = '&&OP0041.CNSLT'
    cnsln = '&&OP0041.CNSLN'
    call cnscre(noma, 'NEUT_R', 1, 'X1', 'V', &
                cnslt)
    call cnscre(noma, 'NEUT_R', 1, 'X1', 'V', &
                cnsln)
    if (me1 .eq. 1) then
        meth = 'FONCTION'
    else if (me2 .eq. 1) then
        meth = 'GROUP_MA'
    else if (me3 .eq. 1) then
        meth = 'GEOMETRI'
    else if (me4 .eq. 1) then
        meth = 'CHAMP'
    else
        ASSERT(.false.)
    end if
!
    call xinils(noma, kbid, .false._1, ndim, meth, &
                nfonf, nfong, geofis, a, b, &
                r, noeud, cote, vect1, vect2, &
                cnslt, cnsln)
!
! --- CREATION DES CHAM_NO DES LEVEL-SETS
!
    ltno = fiss(1:8)//'.LTNO'
    lnno = fiss(1:8)//'.LNNO'
    call cnscno(cnslt, ltno(1:13)//'.NUMEQ', 'NON', 'G', ltno, &
                'F', ibid)
    call cnscno(cnsln, ltno(1:13)//'.NUMEQ', 'NON', 'G', lnno, &
                'F', ibid)
!
    if (niv .ge. 3) then
        call imprsd('CHAMP', ltno, ifm, 'FISSURE.LTNO=')
        call imprsd('CHAMP', lnno, ifm, 'FISSURE.LNNO=')
    end if
!
!-----------------------------------------------------------------------
!     CALCULATE THE LEVEL SETS ON THE AUXILIARY GRID
!-----------------------------------------------------------------------
!
    if (grille) then
!
        cnsltg = '&&OP0041.CNSLTG'
        cnslng = '&&OP0041.CNSLNG'
        call cnscre(maiaux, 'NEUT_R', 1, 'X1', 'V', &
                    cnsltg)
        call cnscre(maiaux, 'NEUT_R', 1, 'X1', 'V', &
                    cnslng)
!
        if (meth(1:5) .ne. 'CHAMP') then
!           THE SAME METHOD "METH" IS USED
            call xinils(noma, maiaux, grille, ndim, meth, &
                        nfonf, nfong, geofis, a, b, &
                        r, noeud, cote, vect1, vect2, &
                        cnsltg, cnslng)
        else
!           IF THE CHAMP_NO_S HAVE BEEN GIVEN, THEY ARE PROJECTED TO THE
!           AUXILIARY GRID. NO OTHER CALCULATIONS ARE POSSIBLE.
            if (typdis .ne. 'INTERFACE') then
                write (ifm, *) '  LES LEVEL SETS DONNEES SONT PROJETEES SUR'&
     &                   //' LA GRILLE AUXILIAIRE.'
            else
                write (ifm, *) '  LA LEVEL SET NORMALE DONNEE EST PROJETEE'&
     &                   //' SUR LA GRILLE AUXILIAIRE.'
            end if
!
            l_dmax = .false.
            dmax = r8maem()
            corres = '&&OP0041.CORRES'
!
!           CREATE THE "CONNECTION" TABLE BETWEEN THE PHYSICAL MESH AND
!           THE AUXILIARY GRID
            if (ndim .eq. 2) then
                call pj2dco('TOUT', noma, maiaux, 0, [0], &
                            0, [0], ' ', ' ', corres, &
                            l_dmax, dmax, 0.d0)
            else
                call pj3dco('TOUT', noma, maiaux, 0, [0], &
                            0, [0], ' ', ' ', corres, &
                            l_dmax, dmax, 0.d0)
            end if
!
!           PROJECT THE NORMAL LEVEL SET
            call cnsprj(cnsln, corres, 'G', cnslng, ibid)
            ASSERT(ibid .eq. 0)
!           PROJECT THE TANGENTIAL LEVEL SET
            call cnsprj(cnslt, corres, 'G', cnsltg, ibid)
            ASSERT(ibid .eq. 0)
!
        end if
!
!
! --- CREATION DES CHAM_NO DES LEVEL-SETS
!
        ltnofa = fiss(1:8)//'.GRI.LTNO'
        lnnofa = fiss(1:8)//'.GRI.LNNO'
        call cnscno(cnslng, ltnofa(1:12)//'L.NUMEQ', 'NON', 'G', lnnofa, &
                    'F', ibid)
        call cnscno(cnsltg, ltnofa(1:12)//'L.NUMEQ', 'NON', 'G', ltnofa, &
                    'F', ibid)
!
        if (niv .ge. 3) then
            call imprsd('CHAMP', lnnofa, ifm, 'FISSURE.GRI.LNNO=')
            call imprsd('CHAMP', ltnofa, ifm, 'FISSURE.GRI.LTNO=')
        end if
!
    end if
!
!
!-----------------------------------------------------------------------
!     CALCUL DES GRADIENTS DES LEVEL-SETS
!
!     ON ENRICHI LA SD FISS_XFEM DE
!         FISS//'.GRLTNO'
!         FISS//'.GRLNNO'
!
!-----------------------------------------------------------------------
!
    grlt = '&&OP0041.GRLT'
    grln = '&&OP0041.GRLN'
!
    call xgrals(noma, lnno, ltno, grlt, grln)
!
! --- CREATION DES CHAM_NO DES GRADIENTS DES LEVEL-SETS
!
    grltno = fiss(1:8)//'.GRLTNO'
    grlnno = fiss(1:8)//'.GRLNNO'
    call cnscno(grlt, grltno(1:13)//'.NUMEQ', 'NON', 'G', grltno, &
                'F', ibid)
    call cnscno(grln, grltno(1:13)//'.NUMEQ', 'NON', 'G', grlnno, &
                'F', ibid)
!
    if (niv .ge. 2) then
        call imprsd('CHAMP', grltno, ifm, 'FISSURE.GRLTNO=')
        call imprsd('CHAMP', grlnno, ifm, 'FISSURE.GRLNNO=')
    end if
!
!-----------------------------------------------------------------------
!     CALCULATE THE GRADIENTS OF THE LEVEL SETS ON THE AUXILIARY GRID
!-----------------------------------------------------------------------
!
    if (grille) then
!
        grltg = '&&OP0041.GRLTG'
        grlng = '&&OP0041.GRLNG'
!
        call xgrals(maiaux, lnnofa, ltnofa, grltg, grlng)
!
!
! --- CREATION DES CHAM_NO DES GRADIENTS DES LEVEL-SETS
!
        grltfa = fiss(1:8)//'.GRI.GRLTNO'
        grlnfa = fiss(1:8)//'.GRI.GRLNNO'
        call cnscno(grltg, grltfa(1:12)//'G.NUMEQ ', 'NON', 'G', grltfa, &
                    'F', ibid)
        call cnscno(grlng, grltfa(1:12)//'G.NUMEQ ', 'NON', 'G', grlnfa, &
                    'F', ibid)
!
        if (niv .ge. 2) then
            call imprsd('CHAMP', grltfa, ifm, 'FISSURE.GRI.GRLTNO=')
            call imprsd('CHAMP', grlnfa, ifm, 'FISSURE.GRI.GRLNNO=')
        end if
!
        call detrsd('CHAM_NO_S', cnsltg)
        call detrsd('CHAM_NO_S', cnslng)
        call detrsd('CHAM_NO_S', grltg)
        call detrsd('CHAM_NO_S', grlng)
!
    end if
!
!
!-----------------------------------------------------------------------
!     CALCUL DE L'ENRICHISSEMENT, DES POINTS DU FOND DE FISSURE
!
!     ON ENRICHI LA SD FISS_XFEM DE
!         FISS//'.MAILFISS.HEAV'
!         FISS//'.MAILFISS.CTIP'
!         FISS//'.MAILFISS.HECT'
!         FISS//'.MAILFISS.MAFOND'
!         FISS//'.FONDFISS'
!         FISS//'.FONDMULT'
!
!-----------------------------------------------------------------------
!
!
! --- RECUPERATION EVENTUELLE DES FISSURES DE JONCTION
!
    cnslj = '&&OP0041.CNSLJ'
    call getvid('JONCTION', 'FISSURE', iocc=1, nbval=0, nbret=me1)
    if (me1 .lt. 0) then
        call xinlsj(noma, ndim, fiss, me1, cnslj)
    end if
    cnsen = '&&OP0041.CNSEN'
    cnsenr = '&&OP0041.CNSENR'
!
    call xenrch(noma, cnslt, cnsln, cnslj, &
                cnsen, cnsenr, ndim, fiss, goinop, &
                lismae, lisnoe)
!
! --- CREATION DU CHAM_NO POUR LE STATUT DES NOEUDS
!
    stno = fiss(1:8)//'.STNO'
    call cnscno(cnsen, stno(1:13)//'.NUMEQ', 'NON', 'G', stno, &
                'F', ibid)
    if (niv .ge. 3) then
        call imprsd('CHAMP', stno, ifm, 'FISSURE.STNO=')
    end if
!
!-----------------------------------------------------------------------
!     CALCUL DE LA BASE LOCALE AU FOND DE FISSURE
!
!     ON ENRICHI LA SD FISS_XFEM DE
!         FISS//'.BASLOC'
!
!-----------------------------------------------------------------------
!
    call xbaslo(noma, fiss, grlt, grln, ndim)
!
! --- MENAGE
!
    call detrsd('CHAM_NO_S', cnslt)
    call detrsd('CHAM_NO_S', cnsln)
    call detrsd('CHAM_NO_S', grlt)
    call detrsd('CHAM_NO_S', grln)
    call detrsd('CHAM_NO_S', cnsen)
    call detrsd('CHAM_NO_S', cnsenr)
    call jedetr(lismae)
    call jedetr(lisnoe)
!
900 format('    LA GRILLE ', a8, ' A ETE ASSOCIEE A LA FISSURE')
!
    call jedema()
end subroutine
