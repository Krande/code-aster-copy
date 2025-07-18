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

subroutine mecalg(optioz, result, modele, depla, theta, &
                  mate, mateco, lischa, symech, compor, incr, &
                  time, iord, nbprup, noprup, chvite, &
                  chacce, kcalc, coor, iadnoe)
!
!     - FONCTION REALISEE:   CALCUL DU TAUX DE RESTITUTION D'ENERGIE
!
!     - ARGUMENTS   :
!
! IN/OUT    OPTION       --> CALC_G    (G SI CHARGES REELLES)
!                        --> CALC_G_F  (G SI CHARGES FONCTIONS)
! IN    RESULT       --> NOM UTILISATEUR DU RESULTAT ET TABLE
! IN    MODELE       --> NOM DU MODELE
! IN    DEPLA        --> CHAMP DES DEPLACEMENTS
! IN    THETA        --> CHAMP THETA (DE TYPE CHAM_NO)
! IN    MATE         --> CHAMP DU MATERIAU
! IN    SYMECH       --> SYMETRIE DU CHARGEMENT
! IN    TIME         --> INSTANT DE CALCUL
! IN    IORD         --> NUMERO D'ORDRE DE LA SD
! IN    KCALC        --> = 'NON' : ON RECUPERE LES CHAMPS DE CONTRAINTES
!                                  ET D'ENERGIE DE LA SD RESULTAT
!                        = 'OUI' :ON RECALCULE LES CHAMPS DE CONTRAINTES
!                                  ET D'ENERGIE
! IN    COOR         --> COORDONNEES ET ABSCISSES CURVILIGNES DES NOEUDS
!                   DU FOND DE FISSURE (IADFIS DANS OP0100)
! IN    IADNOE      --> ADRESSE DU FOND DE FISSURE (IADFIS DANS OP0100)
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
!
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/xelgano.h"
#include "asterfort/calcul.h"
#include "asterfort/chpchd.h"
#include "asterfort/chpver.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcharg.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/mesomm.h"
#include "asterfort/rsexch.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajvi.h"
#include "asterfort/tbajvk.h"
#include "asterfort/tbajvr.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
!
    character(len=8) :: modele, result, symech
    character(len=8) :: kcalc
    character(len=19) :: lischa
    character(len=16) :: optioz, noprup(*)
    character(len=24) :: depla, mate, mateco, compor, theta
    character(len=24) :: chvite, chacce
    real(kind=8) :: time
    integer(kind=8) :: iord, nbprup, iadnoe
    aster_logical :: incr
!
!
! DECLARATION VARIABLES LOCALES
!
    character(len=2) :: codret
    character(len=6) :: nompro
    parameter(nompro='MECALG')
!
    integer(kind=8) :: nbmxpa
    parameter(nbmxpa=20)
!
    integer(kind=8) :: ibid, iret, nres, numfon, livi(nbmxpa), coor
    integer(kind=8) :: nchin, nsig, ino1, ino2, inga, pbtype
    real(kind=8) :: g(1), livr(nbmxpa)
    complex(kind=8) :: livc(nbmxpa)
    aster_logical :: lfonc, lxfem
    character(len=8) :: resu, lpain(50), lpaout(2), k8b, resuco
    character(len=8) :: fiss
    character(len=16) :: option
    character(len=19) :: ch1d2d, ch2d3d, chpres, chrota, chpesa, chvolu, chepsi
    character(len=19) :: chvref, chvarc
    character(len=19) :: basloc, pintto, cnseto, heavto, loncha, lnno, ltno, hea_no
    character(len=19) :: pmilto
    character(len=19) :: pinter, ainter, cface, longco, baseco, stano
    character(len=24) :: modelLigrel, chgeom, lchin(50), lchout(2)
    character(len=24) :: chtime, celmod, sigelno, sigseno
    character(len=24) :: pavolu, pa1d2d, pa2d3d, papres, pepsin
    character(len=24) :: chsig, chepsp, chvari, chsigi, livk(nbmxpa)
    parameter(resuco='&&MECALG')
    data chvarc/'&&MECALG.CH_VARC_R'/
    data chvref/'&&MECALG.CHVREF'/
!
    call jemarq()
    chtime = ' '
    option = optioz
!
!     INITIALISATIONS
    g = 0.d0
    nsig = 0
    inga = 0
    ch1d2d = '&&MECALG.1D2D'
    ch2d3d = '&&MECALG.2D3D'
    chepsi = '&&MECALG.EPSI'
    chpesa = '&&MECALG.PESA'
    chpres = '&&MECALG.PRES'
    chrota = '&&MECALG.ROTA'
    chtime = '&&MECALG.CH_INST_R'
    chvolu = '&&MECALG.VOLU'
    chsigi = '&&MECALG.CHSIGI'
    celmod = '&&MECALG.CELMOD'
    sigelno = '&&MECALG.SIGELNO'
    sigseno = '&&MECALG.SIGSENO'
!
!   cas FEM ou X-FEM
    call getvid('THETA', 'FISSURE', iocc=1, scal=fiss, nbret=ibid)
    lxfem = .false.
    if (ibid .ne. 0) lxfem = .true.
!
!   RECUPERATION DU CHAMP GEOMETRIQUE
    call megeom(modele, chgeom)

!   Recuperation du LIGREL
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=modelLigrel)

!
!- RECUPERATION DU COMPORTEMENT
!
    if (incr) then
        call getvid(' ', 'RESULTAT', scal=resu, nbret=nres)
        call rsexch('F', resu, 'SIEF_ELGA', iord, chsig, &
                    iret)
        call rsexch('F', resuco, 'EPSP_ELNO', iord, chepsp, &
                    iret)
        call rsexch('F', resuco, 'VARI_ELNO', iord, chvari, &
                    iret)
    end if
!
!   Recuperation de l'etat initial
!   ------------------------------

    if (incr) then

        call getvid('ETAT_INIT', 'SIGM', iocc=1, scal=chsigi, nbret=nsig)

!       Verification du type de champ + transfo, si necessaire en champ elno
        if (nsig .ne. 0) then

!           chpver renvoit 0 si OK et 1 si PB
            call chpver('C', chsigi(1:19), 'ELNO', 'SIEF_R', ino1)
            call chpver('C', chsigi(1:19), 'NOEU', 'SIEF_R', ino2)
            call chpver('C', chsigi(1:19), 'ELGA', 'SIEF_R', inga)

!           Verification du type de champ
            pbtype = 0
            if (.not. lxfem) then
!             cas FEM : verif que le champ est soit ELNO, soit NOEU, soit ELGA
                if (ino1 .eq. 1 .and. ino2 .eq. 1 .and. inga .eq. 1) pbtype = 1
            elseif (lxfem) then
!             cas X-FEM : verif que le champ est ELGA (seul cas autorise)
                if (inga .eq. 1) pbtype = 1
            end if
            if (pbtype .eq. 1) call utmess('F', 'RUPTURE1_12')

!           transformation si champ ELGA
            if (inga .eq. 0) then

!               traitement du champ pour les elements finis classiques
                call detrsd('CHAMP', celmod)
                call alchml(modelLigrel, 'CALC_G_XFEM', 'PSIGINR', 'V', celmod, &
                            iret, ' ')
                call chpchd(chsigi(1:19), 'ELNO', celmod, 'OUI', 'V', &
                            sigelno, modele)
                call chpver('F', sigelno(1:19), 'ELNO', 'SIEF_R', ibid)

!               calcul d'un champ supplementaire aux noeuds des sous-elements si X-FEM
                if (lxfem) call xelgano(modele, chsigi, sigseno)
!                call imprsd('CHAMP',chsigi,6,'chsigi')

            end if

        end if

    else

        nsig = 0

    end if
!
!- RECUPERATION (S'ILS EXISTENT) DES CHAMP DE TEMPERATURES (T,TREF)
    k8b = '        '
    call vrcins(modele, mate, k8b, time, chvarc, &
                codret)
    call vrcref(modele, mate(1:8), k8b, chvref)
!
!
! - TRAITEMENT DES CHARGES
!
    call gcharg(modele, lischa, chvolu, ch1d2d, ch2d3d, &
                chpres, chepsi, chpesa, chrota, lfonc, &
                time, iord)
!
    if (lfonc) then
        pavolu = 'PFFVOLU'
        pa1d2d = 'PFF1D2D'
        pa2d3d = 'PFF2D3D'
        papres = 'PPRESSF'
        pepsin = 'PEPSINF'
        if (option .eq. 'CALC_DG') then
            option = 'CALC_DG_F'
        else if (option .eq. 'CALC_G') then
            option = 'CALC_G_XFEM_F'
        else if (option .eq. 'CALC_DG_E') then
            option = 'CALC_DG_E_F'
        else if (option .eq. 'CALC_DGG_E') then
            option = 'CALC_DGG_E_F'
        else if (option .eq. 'CALC_DG_FORC') then
            option = 'CALC_DG_FORC_F'
        else if (option .eq. 'CALC_DGG_FORC') then
            option = 'CALC_DGG_FORC_F'
        end if
    else
        pavolu = 'PFRVOLU'
        pa1d2d = 'PFR1D2D'
        pa2d3d = 'PFR2D3D'
        papres = 'PPRESSR'
        pepsin = 'PEPSINR'
        option = 'CALC_G_XFEM'
    end if
!
!
    if (lxfem) then
!       RECUPERATION DES DONNEES XFEM (TOPOSE)
        pintto = modele//'.TOPOSE.PIN'
        cnseto = modele//'.TOPOSE.CNS'
        heavto = modele//'.TOPOSE.HEA'
        loncha = modele//'.TOPOSE.LON'
        pmilto = modele//'.TOPOSE.PMI'
        lnno = fiss//'.LNNO'
        ltno = fiss//'.LTNO'
        basloc = fiss//'.BASLOC'
!
!       RECUPERATION DES DONNEES XFEM (TOPOFAC)
        pinter = modele//'.TOPOFAC.OE'
        ainter = modele//'.TOPOFAC.AI'
        cface = modele//'.TOPOFAC.CF'
        longco = modele//'.TOPOFAC.LO'
        baseco = modele//'.TOPOFAC.BA'
        stano = modele//'.STNO'
!
!       RECUPERATION DES DONNEES XFEM (TOPONO)
        hea_no = modele//'.TOPONO.HNO'
!
    end if
!
    lpaout(1) = 'PGTHETA'
    lchout(1) = '&&'//nompro//'.CH_G'
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
    lpain(9) = pa2d3d(1:8)
    lchin(9) = ch2d3d
    lpain(10) = papres(1:8)
    lchin(10) = chpres
    lpain(11) = 'PPESANR'
    lchin(11) = chpesa
    lpain(12) = 'PROTATR'
    lchin(12) = chrota
    lpain(13) = pepsin(1:8)
    lchin(13) = chepsi
    lpain(14) = 'PCOMPOR'
    lchin(14) = compor
!
    modelLigrel = modele//'.MODELE'
    nchin = 14
!
    if (lxfem) then
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
        lpain(20) = 'PBASLOR'
        lchin(20) = basloc
        lpain(21) = 'PPINTTO'
        lchin(21) = pintto
        lpain(22) = 'PPMILTO'
        lchin(22) = pmilto
        lpain(23) = 'PPINTER'
        lchin(23) = pinter
        lpain(24) = 'PAINTER'
        lchin(24) = ainter
        lpain(25) = 'PCFACE'
        lchin(25) = cface
        lpain(26) = 'PLONGCO'
        lchin(26) = longco
        lpain(27) = 'PBASECO'
        lchin(27) = baseco
        lpain(28) = 'PHEA_NO'
        lchin(28) = hea_no
        lpain(29) = 'PSTANO'
        lchin(29) = stano
!
        nchin = 29
!
    end if
!
    if ((option .eq. 'CALC_G_XFEM_F') .or. (option .eq. 'CALC_DG_F') &
        .or. (option .eq. 'CALC_DG_E_F') &
        .or. (option .eq. 'CALC_DGG_E_F') &
        .or. (option .eq. 'CALC_DGG_FORC_F') &
        .or. (option .eq. 'CALC_DG_FORC_F')) then
        call mecact('V', chtime, 'MODELE', modelLigrel, 'INST_R  ', &
                    ncmp=1, nomcmp='INST   ', sr=time)
        lpain(nchin+1) = 'PINSTR'
        lchin(nchin+1) = chtime
        nchin = nchin+1
    end if
!
    if (incr) then
        lpain(nchin+1) = 'PCONTRR'
        lchin(nchin+1) = chsig
        lpain(nchin+2) = 'PDEFOPL'
        lchin(nchin+2) = chepsp
        lpain(nchin+3) = 'PVARIPR'
        lchin(nchin+3) = chvari
        nchin = nchin+3
!
!       CHAMP DE CONTRAINTE INITIALE
        if (nsig .ne. 0) then
            if (inga .eq. 0) then
!           champ de contrainte initiale transforme en ELNO
                lpain(nchin+1) = 'PSIGINR'
                lchin(nchin+1) = sigelno
                nchin = nchin+1

!           si X-FEM : champ de contrainte initiale transforme en SE-ELNO
                if (lxfem) then
                    lpain(nchin+1) = 'PSIGISE'
                    lchin(nchin+1) = sigseno
                    nchin = nchin+1
                end if

            else
!           champ de contrainte initiale donne par l'utilisateur (NOEUD ou ELNO)
                lpain(nchin+1) = 'PSIGINR'
                lchin(nchin+1) = chsigi
                nchin = nchin+1
            end if
        end if
    end if
!
    if (chvite .ne. ' ') then
        lpain(nchin+1) = 'PVITESS'
        lchin(nchin+1) = chvite
        lpain(nchin+2) = 'PACCELE'
        lchin(nchin+2) = chacce
        nchin = nchin+2
    end if
!
    if (kcalc .eq. 'NON') then
        call getvid(' ', 'RESULTAT', scal=resu, nbret=iret)
        call rsexch(' ', resu, 'SIEF_ELGA', iord, chsig, &
                    iret)
        lpain(nchin+1) = 'PCONTGR'
        lchin(nchin+1) = chsig
        nchin = nchin+1
    end if
!
!
!-  SOMMATION DES G ELEMENTAIRES
    call calcul('S', option, modelLigrel, nchin, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
    call mesomm(lchout(1), 1, vr=g(1))
    if (symech .ne. 'NON') then
        g(1) = 2.d0*g(1)
    end if
!
!- IMPRESSION DE G ET ECRITURE DANS LA TABLE RESULT
!
    call getvis('THETA', 'NUME_FOND', iocc=1, scal=numfon, nbret=ibid)
    call tbajvi(result, nbprup, 'NUME_FOND', numfon, livi)

! NOM DES NOEUDS DU FOND
    if (.not. lxfem) then
        call tbajvk(result, nbprup, 'NOEUD', zk8(iadnoe), livk)
    end if
!
    call tbajvi(result, nbprup, 'NUME_ORDRE', iord, livi)
    call tbajvr(result, nbprup, 'INST', time, livr)
!
    call tbajvr(result, nbprup, 'COOR_X', zr(coor), livr)
    call tbajvr(result, nbprup, 'COOR_Y', zr(coor+1), livr)

    call tbajvr(result, nbprup, 'G', g(1), livr)
    call tbajli(result, nbprup, noprup, livi, livr, &
                livc, livk, 0)
!
    call detrsd('CHAMP_GD', ch1d2d)
    call detrsd('CHAMP_GD', ch2d3d)
    call detrsd('CHAMP_GD', chepsi)
    call detrsd('CHAMP_GD', chpesa)
    call detrsd('CHAMP_GD', chpres)
    call detrsd('CHAMP_GD', chrota)
    call detrsd('CHAMP_GD', chtime)
    call detrsd('CHAMP_GD', chvolu)

    call jedema()
end subroutine
