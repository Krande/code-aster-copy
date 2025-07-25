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

subroutine cakg2d(optioz, result, modele, depla, theta, &
                  mate, mateco, lischa, symech, fondf, noeud, &
                  time, iord, nbprup, noprup, &
                  lmoda, puls, compor)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/chpver.h"
#include "asterfort/alchml.h"
#include "asterfort/chpchd.h"
#include "asterfort/xelgano.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcharg.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/impfic.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisnnb.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/mesomm.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajvi.h"
#include "asterfort/tbajvk.h"
#include "asterfort/tbajvr.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
    character(len=8) :: modele, fondf, result, symech
    character(len=8) :: noeud
    character(len=16) :: optioz, noprup(*)
    character(len=24) :: depla, mate, mateco, theta, compor
    character(len=19) :: lischa
    real(kind=8) :: time, puls
    integer(kind=8) :: iord, nbprup
    aster_logical :: lmoda
! ......................................................................
!
!     - FONCTION REALISEE:   CALCUL DES COEFFICIENTS D'INTENSITE DE
!                            CONTRAINTES K1 ET K2 EN 2D
!
! IN   OPTION  --> CALC_K_G   (SI CHARGES REELLES)
!              --> CALC_K_G_F (SI CHARGES FONCTIONS)
!              --> CALC_K_X   (SI FISSURE X-FEM)
! IN   RESULT  --> NOM UTILISATEUR DU RESULTAT ET TABLE
! IN   MODELE  --> NOM DU MODELE
! IN   DEPLA   --> CHAMPS DE DEPLACEMENT
! IN   THETA   --> CHAMP THETA
! IN   MATE    --> CHAMP DE MATERIAUX
! IN   LISCHA  --> LISTE DES CHARGES
! IN   SYMECH  --> SYMETRIE DU CHARGEMENT
! IN   FONDF   --> FOND DE FISSURE
! IN   NOEUD   --> NOM DU NOEUD DU FOND
! IN   TIME    --> INSTANT DE CALCUL
! IN   IORD    --> NUMERO D'ORDRE DE LA SD
! IN   NBPRUP  --> NOMBRE DE PARAMETRES RUPTURE DANS LA TABLE
! IN   NOPRUP  --> NOMS DES PARAMETRES RUPTURE
! IN   LMODA   --> TRUE SI LE TYPE DE LA SD RESULTAT = MODE_MECA
! IN   PULS    --> PULSATION SI LMODA
! IN   COMPOR  --> COMPORTEMENT
! ......................................................................
! CORPS DU PROGRAMME
!
    integer(kind=8) :: nbmxpa
    parameter(nbmxpa=20)
!
    integer(kind=8) :: nbinmx, nboumx
    parameter(nbinmx=50, nboumx=1)
    character(len=8) :: lpain(nbinmx), lpaout(nboumx)
    character(len=24) :: lchin(nbinmx), lchout(nboumx)
    integer(kind=8) :: i, ibid, inorma, nsig, ifm, niv, jnor, jbasfo
    integer(kind=8) :: iadrma, iadrff, icoode, iadrco, iadrno, ino1, ino2, inga
    integer(kind=8) :: lobj2, ndimte, nunoff, ndim, nchin, jfond, numfon
    integer(kind=8) :: iret, livi(nbmxpa), nbchar, pbtype
    real(kind=8) :: fic(5), rcmp(6), livr(nbmxpa), girwin
    integer(kind=8) :: mxstac
    complex(kind=8) :: livc(nbmxpa)
    aster_logical :: lfonc, lxfem
    parameter(mxstac=1000)
    character(len=2) :: codret
    character(len=8) :: noma, fond, licmp(6), is_axi, fiss
    character(len=16) :: option, optio2
    character(len=19) :: ch1d2d, chpres, chrota, chpesa, chvolu, ch2d3d, chepsi
    character(len=19) :: chvref, chvarc
    character(len=19) :: basefo
    character(len=19) :: basloc, pintto, cnseto, heavto, loncha, lnno, ltno, hea_no
    character(len=19) :: pmilto
    character(len=19) :: pinter, ainter, cface, longco, baseco, stano
    character(len=24) :: chgeom, chfond, celmod, sigelno, sigseno
    character(len=24) :: ligrmo, norma
    character(len=24) :: obj1, obj2, coord, coorn, chtime
    character(len=24) :: pavolu, pa1d2d, papres, chpuls, chsigi, livk(nbmxpa)
    real(kind=8), pointer :: valg(:) => null()
!
    data chvarc/'&&CAKG2D.CH_VARC_R'/
    data chvref/'&&CAKG2D.CHVREF'/
!
!
    call jemarq()
    iadrno = 1

!
!     VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!     (VOIR CRS 1404)
!
    call lisnnb(lischa, nbchar)
    ASSERT(nbchar .le. mxstac)
    call infniv(ifm, niv)
    option = optioz
    if (optioz .eq. 'CALC_K_X') option = 'CALC_K_G_XFEM'
!
!   cas FEM ou X-FEM
    call getvid('THETA', 'FISSURE', iocc=1, scal=fiss, nbret=ibid)
    lxfem = .false.
    if (ibid .ne. 0) lxfem = .true.

!   RECUPERATION DU CHAMP GEOMETRIQUE
    call megeom(modele, chgeom)
    noma = chgeom(1:8)

!   Recuperation du LIGREL
    ligrmo = modele//'.MODELE'
!
!   Recuperation de l'etat initial
!   ------------------------------

    chsigi = '&&CAKG2D.CHSIGI'
    celmod = '&&CAKG2D.CELMOD'
    sigelno = '&&CAKG2D.SIGELNO'
    sigseno = '&&CAKG2D.SIGSENO'

!   RECUPERATION DE L'ETAT INITIAL
    call getvid('ETAT_INIT', 'SIGM', iocc=1, scal=chsigi, nbret=nsig)

!   Verification du type de champ + transfo, si necessaire en champ elno
    if (nsig .ne. 0) then

!       chpver renvoit 0 si OK et 1 si PB
        call chpver('C', chsigi(1:19), 'ELNO', 'SIEF_R', ino1)
        call chpver('C', chsigi(1:19), 'NOEU', 'SIEF_R', ino2)
        call chpver('C', chsigi(1:19), 'ELGA', 'SIEF_R', inga)

!       Verification du type de champ
        pbtype = 0
        if (.not. lxfem) then
!         cas FEM : verif que le champ est soit ELNO, soit NOEU, soit ELGA
            if (ino1 .eq. 1 .and. ino2 .eq. 1 .and. inga .eq. 1) pbtype = 1
        elseif (lxfem) then
!         cas X-FEM : verif que le champ est ELGA (seul cas autorise)
            if (inga .eq. 1) pbtype = 1
        end if
        if (pbtype .eq. 1) call utmess('F', 'RUPTURE1_12')

!       transformation si champ ELGA
        if (inga .eq. 0) then

!           traitement du champ pour les elements finis classiques
            call detrsd('CHAMP', celmod)
            call alchml(ligrmo, 'CALC_G_XFEM', 'PSIGINR', 'V', celmod, &
                        iret, ' ')
            call chpchd(chsigi(1:19), 'ELNO', celmod, 'OUI', 'V', &
                        sigelno, modele)
            call chpver('F', sigelno(1:19), 'ELNO', 'SIEF_R', ibid)

!           calcul d'un champ supplementaire aux noeuds des sous-elements si X-FEM
            if (lxfem) call xelgano(modele, chsigi, sigseno)
!            call imprsd('CHAMP',chsigi,6,'chsigi')

        end if
    end if
!
!
!   RECUPERATION (S'ILS EXISTENT) DES CHAMP DE TEMPERATURES (T,TREF)
    call vrcins(modele, mate, 'BIDON', time, chvarc, &
                codret)
    call vrcref(modele, mate(1:8), 'BIDON   ', chvref(1:19))
!
! - TRAITEMENT DES CHARGES
!
    chvolu = '&&CAKG2D.VOLU'
    ch1d2d = '&&CAKG2D.1D2D'
    ch2d3d = '&&CAKG2D.2D3D'
    chpres = '&&CAKG2D.PRES'
    chepsi = '&&CAKG2D.EPSI'
    chpesa = '&&CAKG2D.PESA'
    chrota = '&&CAKG2D.ROTA'
    call gcharg(modele, lischa, chvolu, ch1d2d, ch2d3d, &
                chpres, chepsi, chpesa, chrota, lfonc, &
                time, iord)
    if (lfonc) then
        pavolu = 'PFFVOLU'
        pa1d2d = 'PFF1D2D'
        papres = 'PPRESSF'
        option = 'CALC_K_G_XFEM_F'
        optio2 = 'CALC_K_G_XFEM_F'
    else
        pavolu = 'PFRVOLU'
        pa1d2d = 'PFR1D2D'
        papres = 'PPRESSR'
        optio2 = 'CALC_K_G_XFEM'
    end if
!
! OBJET DECRIVANT LE MAILLAGE
!
    obj1 = modele//'.MODELE    .LGRF'
    call jeveuo(obj1, 'L', iadrma)
    noma = zk8(iadrma)
    coorn = noma//'.COORDO    .VALE'
    coord = noma//'.COORDO    .DESC'
    call jeveuo(coorn, 'L', iadrco)
    call jeveuo(coord, 'L', icoode)
    ndim = -zi(icoode-1+2)
!
!
!   modele AXIS ?
    call dismoi('AXIS', modele, 'MODELE', repk=is_axi)
!
!   OBJET CONTENANT LES NOEUDS DU FOND DE FISSURE
    if (.not. lxfem) then
        fond = fondf(1:8)
        obj2 = fond//'.FOND.NOEU'
        call jelira(obj2, 'LONMAX', lobj2)
        if (lobj2 .ne. 1) then
            call utmess('F', 'RUPTURE1_10')
        end if
        call jeveuo(obj2, 'L', iadrno)
        nunoff = char8_to_int(zk8(iadrno))
!
!       OBJET CONTENANT LA BASE LOCALE AU FOND DE FISSURE
!       SI L'OBJET NORMALE EXISTE, ON LE PREND
!       SINON, ON PREND BASEFOND
        norma = fond//'.NORMALE'
        call jeexin(norma, iret)
        if (iret .ne. 0) then
            call jeveuo(norma, 'L', inorma)
            rcmp(3) = zr(inorma-1+1)
            rcmp(4) = zr(inorma-1+2)
            rcmp(5) = -zr(inorma-1+2)
            rcmp(6) = zr(inorma-1+1)
        else if (iret .eq. 0) then
            basefo = fond//'.BASEFOND'
            call jeveuo(basefo, 'L', jbasfo)
!           ATTENTION, ON NE SE SERT PAS DU VECTEUR NORMAL DE BASEFOND
!           MAIS ON FAIT TOURNER DE 90 DEGRES LE VECTEUR DE PROPA
            rcmp(3) = zr(jbasfo-1+3)
            rcmp(4) = zr(jbasfo-1+4)
            rcmp(5) = -zr(jbasfo-1+4)
            rcmp(6) = zr(jbasfo-1+3)

        end if
    end if
!
!
!   CREATION OBJET CONTENANT COORDONNEES DU NOEUD DE FOND
!   DE FISSURE ET LA NORMALE A LA FISSURE
    chfond = '&&CAKG2D.FOND'
    call wkvect(chfond, 'V V R8', 6, iadrff)
!
    licmp(1) = 'XA'
    licmp(2) = 'YA'
    licmp(3) = 'XTAN'
    licmp(4) = 'YTAN'
    licmp(5) = 'XNORM'
    licmp(6) = 'YNORM'

!Recuperation du numero du fond de fissure pour FEM et XFEM
    call getvis('THETA', 'NUME_FOND', iocc=1, scal=numfon, nbret=ibid)

    if (.not. lxfem) then
!       cas FEM
        rcmp(1) = zr(iadrco+ndim*(nunoff-1))
        rcmp(2) = zr(iadrco+ndim*(nunoff-1)+1)
    else
!       cas X-FEM
        call jeveuo(fiss//'.FONDFISS', 'L', jfond)
        rcmp(1) = zr(jfond-1+4*(numfon-1)+1)
        rcmp(2) = zr(jfond-1+4*(numfon-1)+2)
        call jeveuo(fiss//'.BASEFOND', 'L', jnor)
        rcmp(3) = zr(jnor-1+4*(numfon-1)+3)
        rcmp(4) = zr(jnor-1+4*(numfon-1)+4)
        rcmp(5) = -zr(jnor-1+4*(numfon-1)+4)
        rcmp(6) = zr(jnor-1+4*(numfon-1)+3)
!         rcmp(5) = zr(jnor-1+4*(numfon-1)+1)
!         rcmp(6) = zr(jnor-1+4*(numfon-1)+2)
        write (ifm, *) '   '
        write (ifm, *) '    TRAITEMENT DU FOND DE FISSURE NUMERO ', numfon
        write (ifm, *) '    NOMME ', noeud
        write (ifm, *) '    DE COORDONNEES', rcmp(1), rcmp(2)
        write (ifm, *) '    LA NORMALE A LA FISSURE EN CE POINT EST ',&
     &                                              rcmp(5), rcmp(6)
    end if
    zr(iadrff) = rcmp(1)
    zr(iadrff+1) = rcmp(2)
    zr(iadrff+2) = rcmp(3)
    zr(iadrff+3) = rcmp(4)
    zr(iadrff+4) = rcmp(5)
    zr(iadrff+5) = rcmp(6)
!
    call mecact('V', chfond, 'MAILLA', noma, 'FISS_R', &
                ncmp=6, lnomcmp=licmp, vr=rcmp)
!
!
!   RECUPERATION DES DONNEES XFEM (TOPOSE)
    pintto = modele//'.TOPOSE.PIN'
    cnseto = modele//'.TOPOSE.CNS'
    heavto = modele//'.TOPOSE.HEA'
    loncha = modele//'.TOPOSE.LON'
    pmilto = modele//'.TOPOSE.PMI'
!   ON NE PREND PAS LES LSN ET LST DU MODELE
!   CAR LES CHAMPS DU MODELE SONT DEFINIS QUE AUTOUR DE LA FISSURE
!   OR ON A BESOIN DE LSN ET LST MEME POUR LES ELEMENTS CLASSIQUES
    lnno = fiss//'.LNNO'
    ltno = fiss//'.LTNO'
    basloc = fiss//'.BASLOC'
!
!   RECUPERATION DES DONNEES XFEM (TOPOFAC)
    pinter = modele//'.TOPOFAC.OE'
    ainter = modele//'.TOPOFAC.AI'
    cface = modele//'.TOPOFAC.CF'
    longco = modele//'.TOPOFAC.LO'
    baseco = modele//'.TOPOFAC.BA'
    stano = modele//'.STNO'
!
!   RECUPERATION DES DONNEES XFEM (TOPONO)
    hea_no = modele//'.TOPONO.HNO'
!
    ndimte = 5
    AS_ALLOCATE(vr=valg, size=ndimte)
    lpaout(1) = 'PGTHETA'
    lchout(1) = '&&FICGELE'
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PDEPLAR'
    lchin(2) = depla
    lpain(3) = 'PTHETAR'
    lchin(3) = theta
    lpain(4) = 'PMATERC'
    lchin(4) = mateco
    lpain(5) = 'PVARCPR'
    lchin(5) = chvarc
    lpain(6) = 'PVARCRR'
    lchin(6) = chvref
    lpain(7) = pavolu(1:8)
    lchin(7) = chvolu
    lpain(8) = pa1d2d(1:8)
    lchin(8) = ch1d2d
    lpain(9) = papres(1:8)
    lchin(9) = chpres
    lpain(10) = 'PPESANR'
    lchin(10) = chpesa
    lpain(11) = 'PROTATR'
    lchin(11) = chrota
    lpain(12) = 'PFISSR'
    lchin(12) = chfond
!
    lpain(13) = 'PBASLOR'
    lchin(13) = basloc
    lpain(14) = 'PPINTTO'
    lchin(14) = pintto
    lpain(15) = 'PCNSETO'
    lchin(15) = cnseto
    lpain(16) = 'PHEAVTO'
    lchin(16) = heavto
    lpain(17) = 'PLONCHA'
    lchin(17) = loncha
    lpain(18) = 'PLSN'
    lchin(18) = lnno
    lpain(19) = 'PLST'
    lchin(19) = ltno
!
    lpain(20) = 'PCOMPOR'
    lchin(20) = compor
!
    lpain(21) = 'PPMILTO'
    lchin(21) = pmilto
!
    lpain(22) = 'PPINTER'
    lchin(22) = pinter
    lpain(23) = 'PAINTER'
    lchin(23) = ainter
    lpain(24) = 'PCFACE'
    lchin(24) = cface
    lpain(25) = 'PLONGCO'
    lchin(25) = longco
    lpain(26) = 'PBASECO'
    lchin(26) = baseco
    lpain(27) = 'PHEA_NO'
    lchin(27) = hea_no
    lpain(28) = 'PSTANO'
    lchin(28) = stano
!
    nchin = 28
!
    chtime = '&&CAKG2D.CH_INST_R'
    if (option .eq. 'CALC_K_G_XFEM_F') then
        call mecact('V', chtime, 'MODELE', ligrmo, 'INST_R  ', &
                    ncmp=1, nomcmp='INST   ', sr=time)
        nchin = nchin+1
        lpain(nchin) = 'PINSTR'
        lchin(nchin) = chtime
    end if
!
    if (lmoda) then
        chpuls = '&&CAKG2D.PULS'
        call mecact('V', chpuls, 'MODELE', ligrmo, 'FREQ_R  ', &
                    ncmp=1, nomcmp='FREQ   ', sr=puls)
        nchin = nchin+1
        lpain(nchin) = 'PPULPRO'
        lchin(nchin) = chpuls
    end if

!   CHAMP DE CONTRAINTE INITIALE
    if (nsig .ne. 0) then
        if (inga .eq. 0) then
!       champ de contrainte initiale transforme en ELNO
            lpain(nchin+1) = 'PSIGINR'
            lchin(nchin+1) = sigelno
            nchin = nchin+1

!       si X-FEM : champ de contrainte initiale transforme en SE-ELNO
            if (lxfem) then
                lpain(nchin+1) = 'PSIGISE'
                lchin(nchin+1) = sigseno
                nchin = nchin+1
            end if

        else
!       champ de contrainte initiale donne par l'uutilisateur (NOEUD ou ELNO)
            lpain(nchin+1) = 'PSIGINR'
            lchin(nchin+1) = chsigi
            nchin = nchin+1
        end if
    end if
!
    call calcul('S', optio2, ligrmo, nchin, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')

!  SOMMATION DES FIC ET G ELEMENTAIRES
!
    call mesomm(lchout(1), 5, vr=fic)
!
    do i = 1, 5
        valg(i) = fic(i)
    end do
!
    if (is_axi(1:3) .eq. 'OUI') then
        do i = 1, 5
            valg(i) = valg(i)/rcmp(1)
        end do
    end if
!
    if (symech .eq. 'OUI') then
        valg(1) = 2.d0*valg(1)
        valg(1+1) = 2.d0*valg(1+1)
        valg(1+2) = 0.d0
        valg(1+3) = 2.d0*valg(1+3)
        valg(1+4) = 0.d0
    end if
!
    girwin = valg(1+1)*valg(1+1)+valg(1+2)*valg(1+2)
!
! IMPRESSION DE K1,K2,G ET ECRITURE DANS LA TABLE RESU
!
    if (niv .ge. 2) then
        call impfic(valg, zk8(iadrno), rcmp, ifm, lxfem)
    end if
!
!    if (lxfem .and. (option(1:8).eq.'CALC_K_G') .and. (.not.lmoda)) then
!     call tbajvi(result, nbprup, 'NUME_FOND', numfon, livi)
!    endif
!
    call tbajvi(result, nbprup, 'NUME_FOND', numfon, livi)

    if (.not. lxfem) then
        call tbajvk(result, nbprup, 'NOEUD', zk8(iadrno), livk)
    end if
    if (lmoda) then
        call tbajvi(result, nbprup, 'NUME_MODE', iord, livi)
    else
        call tbajvi(result, nbprup, 'NUME_ORDRE', iord, livi)
        call tbajvr(result, nbprup, 'INST', time, livr)
    end if
    call tbajvr(result, nbprup, 'COOR_X', rcmp(1), livr)
    call tbajvr(result, nbprup, 'COOR_Y', rcmp(2), livr)
    call tbajvr(result, nbprup, 'K1', valg(1+3), livr)
    call tbajvr(result, nbprup, 'K2', valg(1+4), livr)
    call tbajvr(result, nbprup, 'G', valg(1), livr)
    call tbajvr(result, nbprup, 'G_IRWIN', girwin, livr)
    call tbajli(result, nbprup, noprup, livi, livr, &
                livc, livk, 0)
!
    call detrsd('CHAMP_GD', chtime)
    call detrsd('CHAMP_GD', chvolu)
    call detrsd('CHAMP_GD', ch1d2d)
    call detrsd('CHAMP_GD', ch2d3d)
    call detrsd('CHAMP_GD', chpres)
    call detrsd('CHAMP_GD', chepsi)
    call detrsd('CHAMP_GD', chpesa)
    call detrsd('CHAMP_GD', chrota)
    AS_DEALLOCATE(vr=valg)
    call jedetr('&&CAKG2D.FOND')
!
    call jedema()
end subroutine
