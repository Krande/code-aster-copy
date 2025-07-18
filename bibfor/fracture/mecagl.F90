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

subroutine mecagl(option, result, modele, depla, thetai, &
                  mate, mateco, compor, lischa, symech, chfond, &
                  nnoff, iord, ndeg, liss, &
                  milieu, ndimte, extim, &
                  time, nbprup, noprup, chvite, chacce, &
                  kcalc, fonoeu, lincr, coor, &
                  norfon, connex)
! aslint: disable=W1504
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alchml.h"
#include "asterfort/calcul.h"
#include "asterfort/chpchd.h"
#include "asterfort/chpver.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcharg.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/gimpgs.h"
#include "asterfort/gmeth1.h"
#include "asterfort/gmeth2.h"
#include "asterfort/gmeth3.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
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
#include "asterfort/wkvect.h"
#include "asterfort/xelgano.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: iord, nbprup, ndimte, coor
!
    real(kind=8) :: time
!
    character(len=19) :: lischa
    character(len=8) :: modele, thetai
    character(len=8) :: result, symech, kcalc
    character(len=16) :: option, noprup(*)
    character(len=24) :: depla, chfond, mate, compor, mateco
    character(len=24) :: chvite, chacce, fonoeu, liss, norfon
!
    aster_logical :: extim, milieu, lincr, connex
! ......................................................................
!
!  - FONCTION REALISEE:   CALCUL DU TAUX DE RESTITUTION LOCAL D'ENERGIE
!
!  IN    OPTION --> CALC_G OU G_LAGR (SI CHARGES REELLES)
!               --> CALC_G_F OU G_LAGR_F (SI CHARGES FONCTIONS)
!  IN    RESULT --> NOM UTILISATEUR DU RESULTAT ET TABLE
!  IN    MODELE --> NOM DU MODELE
!  IN    DEPLA  --> CHAMP DE DEPLACEMENT
!  IN    THETAI --> BASE DE I CHAMPS THETA
!  IN    MATE   --> CHAMP DE MATERIAUX
!  IN    COMPOR --> COMPORTEMENT
!  IN    NCHAR  --> NOMBRE DE CHARGES
!  IN    LCHAR  --> LISTE DES CHARGES
!  IN    SYMECH --> SYMETRIE DU CHARGEMENT
!  IN    CHFOND --> VECTEUR CONTENANT LES ABSCISSES CURVILIGNES DES
!                   NOEUDS DU FOND DE FISSURE
!  IN    NNOFF  --> NOMBRE DE NOEUDS DU FOND DE FISSURE
!  IN    TIME   --> INSTANT DE CALCUL
!  IN    IORD   --> NUMERO D'ORDRE DE LA SD
!  IN    LISS   --> TYPE DE LISSAGE
!  IN    NDEG   --> DEGRE DU POLYNOME DE LEGENDRE
!  IN    KCALC  --> = 'NON' : ON RECUPERE LES CHAMPS DE CONTRAINTES
!                             ET D'ENERGIE DE LA SD RESULTAT
!                   = 'OUI' : ON RECALCULE LES CHAMPS DE CONTRAINTES
!                             ET D'ENERGIE
!  IN    FONOEU --> NOM DES NOEUDS DE FOND DE FISSURE
!  IN    COOR   --> COORDONNEES ET ABSCISSES CURVILIGNES DES NOEUDS
!                   DU FOND DE FISSURE (IADFIS DANS OP0100)
! ......................................................................
!
    integer(kind=8), parameter :: nbmxpa = 20
!
    integer(kind=8) :: i, ibid, iadrg, iret, jresu, nchin
    integer(kind=8) :: nnoff, num, incr, nres, nsig, ino1, ino2, inga, pbtype
    integer(kind=8) :: ndeg, livi(nbmxpa), numfon
    integer(kind=8) :: iadrno, iadgi, iadabs, ifm, niv, ifon
    real(kind=8) :: gthi(1), livr(nbmxpa), xl
    complex(kind=8) :: livc(nbmxpa)
    aster_logical :: fonc, lxfem
    character(len=2) :: codret
    character(len=8) :: resu, fiss
    character(len=8) :: lpain(50), lpaout(1)
    character(len=16) :: opti
    character(len=19) :: chrota, chpesa, cf2d3d, chpres, chvolu, cf1d2d, chepsi
    character(len=19) :: chvarc, chvref
    character(len=19) :: basloc, pintto, cnseto, heavto, loncha, lnno, ltno, stano
    character(len=19) :: pmilto, hea_no
    character(len=19) :: longco, pinter, ainter, cface, baseco
    character(len=24) :: modelLigrel, chgeom, chgthi
    character(len=24) :: chsigi, sigelno, sigseno, celmod
    character(len=24) :: lchin(50), lchout(1), chthet, chtime
    character(len=24) :: objcur, normff, pavolu, papres, pa2d3d
    character(len=24) :: chsig, chepsp, chvari, type, pepsin, livk(nbmxpa)
    real(kind=8), pointer :: valg_s(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call infniv(ifm, niv)
!
    chvarc = '&&MECAGL.VARC'
    chvref = '&&MECAGL.VARC.REF'
    chsigi = '&&MECALG.CHSIGI'
    celmod = '&&MECALG.CELMOD'
    sigelno = '&&MECALG.SIGELNO'
    sigseno = '&&MECALG.SIGSENO'
!- RECUPERATION DU CHAMP GEOMETRIQUE
!
    call megeom(modele, chgeom)
!
!   Recuperation du LIGREL
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=modelLigrel)

    call getvid('THETA', 'FISSURE', iocc=1, scal=fiss, nbret=ibid)
    lxfem = .false.
    if (ibid .ne. 0) lxfem = .true.
!
!- RECUPERATION DU COMPORTEMENT
!
!    call getfac('COMPORTEMENT', incr)
!
    incr = 0
    if (lincr) incr = 1
!
    if (incr .ne. 0) then
        call getvid(' ', 'RESULTAT', scal=resu, nbret=nres)
        call dismoi('TYPE_RESU', resu, 'RESULTAT', repk=type)
        if (type .ne. 'EVOL_NOLI') then
            call utmess('F', 'RUPTURE1_15')
        end if
        call rsexch('F', resu, 'SIEF_ELGA', iord, chsig, &
                    iret)
        call rsexch('F', resu, 'EPSP_ELNO', iord, chepsp, &
                    iret)
        call rsexch('F', resu, 'VARI_ELNO', iord, chvari, &
                    iret)
    end if
!
!   Recuperation de l'etat initial
!   ------------------------------

    if (incr .ne. 0) then

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
    call vrcins(modele, mate, ' ', time, chvarc, &
                codret)
    call vrcref(modele, mate(1:8), '        ', chvref(1:19))
!
! - TRAITEMENT DES CHARGES
!
    chvolu = '&&MECAGL.VOLU'
    cf1d2d = '&&MECAGL.1D2D'
    cf2d3d = '&&MECAGL.2D3D'
    chpres = '&&MECAGL.PRES'
    chepsi = '&&MECAGL.EPSI'
    chpesa = '&&MECAGL.PESA'
    chrota = '&&MECAGL.ROTA'
    call gcharg(modele, lischa, chvolu, cf1d2d, cf2d3d, &
                chpres, chepsi, chpesa, chrota, fonc, &
                time, iord)
    if (fonc) then
        pavolu = 'PFFVOLU'
        pa2d3d = 'PFF2D3D'
        papres = 'PPRESSF'
        pepsin = 'PEPSINF'
        if (option .eq. 'CALC_G') then
            opti = 'CALC_G_XFEM_F'
        else
            opti = 'G_LAGR_F'
        end if
    else
        pavolu = 'PFRVOLU'
        pa2d3d = 'PFR2D3D'
        papres = 'PPRESSR'
        pepsin = 'PEPSINR'
        if (option .eq. 'CALC_G') then
            opti = 'CALC_G_XFEM'
        else
            opti = 'G_LAGR'
        end if
    end if
!
!- CALCUL DES G(THETA_I) AVEC I=1,NDIMTE  NDIMTE = NNOFF  SI TH-LAGRANGE
!                                         NDIMTE = NDEG+1 SI TH-LEGENDRE
    if ((liss .eq. 'LAGRANGE') .or. (liss .eq. 'LAGRANGE_NO_NO') .or. (liss .eq. 'MIXTE')) then
        ndimte = nnoff
    else
        ndimte = ndeg+1
    end if
!
    call wkvect('&&MECAGL.VALG', 'V V R8', ndimte, iadrg)
    call jeveuo(thetai, 'L', jresu)
!
! --- RECUPERATION DES DONNEES X-FEM
    if (lxfem) then
        pintto = modele//'.TOPOSE.PIN'
        cnseto = modele//'.TOPOSE.CNS'
        heavto = modele//'.TOPOSE.HEA'
        hea_no = modele//'.TOPONO.HNO'
        loncha = modele//'.TOPOSE.LON'
        pmilto = modele//'.TOPOSE.PMI'
!       ON NE PREND PAS LES LSN ET LST DU MODELE
!       CAR LES CHAMPS DU MODELE SONT DEFINIS QUE AUTOUR DE LA FISSURE
!       OR ON A BESOIN DE LSN ET LST MEME POUR LES
        lnno = fiss//'.LNNO'
        ltno = fiss//'.LTNO'
        basloc = fiss//'.BASLOC'
        longco = modele//'.TOPOFAC.LO'
        pinter = modele//'.TOPOFAC.OE'
        ainter = modele//'.TOPOFAC.AI'
        cface = modele//'.TOPOFAC.CF'
        baseco = modele//'.TOPOFAC.BA'
        stano = modele//'.STNO'
    end if
!
    do i = 1, ndimte
        chthet = zk24(jresu+i-1)
        call codent(i, 'G', chgthi)
        lpaout(1) = 'PGTHETA'
        lchout(1) = chgthi
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PDEPLAR'
        lchin(2) = depla
        lpain(3) = 'PTHETAR'
        lchin(3) = chthet
        lpain(4) = 'PMATERC'
        lchin(4) = mateco
        lpain(5) = 'PVARCPR'
        lchin(5) = chvarc
        lpain(6) = 'PVARCRR'
        lchin(6) = chvref
        lpain(7) = pavolu(1:8)
        lchin(7) = chvolu
        lpain(8) = pa2d3d(1:8)
        lchin(8) = cf2d3d
        lpain(9) = papres(1:8)
        lchin(9) = chpres
        lpain(10) = 'PPESANR'
        lchin(10) = chpesa
        lpain(11) = 'PROTATR'
        lchin(11) = chrota
        lpain(12) = pepsin(1:8)
        lchin(12) = chepsi
        lpain(13) = 'PCOMPOR'
        lchin(13) = compor
!
        nchin = 13
!
        if (lxfem) then
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
            lpain(20) = 'PBASLOR'
            lchin(20) = basloc
!
            lpain(21) = 'PLONGCO'
            lchin(21) = longco
            lpain(22) = 'PPINTER'
            lchin(22) = pinter
            lpain(23) = 'PAINTER'
            lchin(23) = ainter
            lpain(24) = 'PCFACE'
            lchin(24) = cface
            lpain(25) = 'PPMILTO'
            lchin(25) = pmilto
            lpain(26) = 'PBASECO'
            lchin(26) = baseco
            lpain(27) = 'PHEA_NO'
            lchin(27) = hea_no
            lpain(28) = 'PSTANO'
            lchin(28) = stano
!
            nchin = 28
        end if
!
        if ((opti .eq. 'CALC_G_XFEM_F') .or. (opti .eq. 'G_LAGR_F')) then
            chtime = '&&MECAGL.CH_INST_R'
            call mecact('V', chtime, 'MODELE', modelLigrel, 'INST_R', &
                        ncmp=1, nomcmp='INST', sr=time)
            lpain(nchin+1) = 'PINSTR'
            lchin(nchin+1) = chtime
            nchin = nchin+1
        end if
        if (incr .ne. 0) then
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
                    lpain(nchin+1) = 'PSIGINR'
                    lchin(nchin+1) = sigelno
                    nchin = nchin+1
                    lpain(nchin+1) = 'PSIGING'
                    lchin(nchin+1) = chsigi
                    nchin = nchin+1
!                   si X-FEM : champ de contrainte initiale transforme en SE-ELNO
                    if (lxfem) then
                        lpain(nchin+1) = 'PSIGISE'
                        lchin(nchin+1) = sigseno
                        nchin = nchin+1
                    end if
                else
                    lpain(nchin+1) = 'PSIGINR'
                    lchin(nchin+1) = chsigi
                    nchin = nchin+1
                end if
            end if
        end if
        if (opti .eq. 'CALC_G_XFEM' .or. opti .eq. 'CALC_G_XFEM_F') then
            if (chvite .ne. ' ') then
                lpain(nchin+1) = 'PVITESS'
                lchin(nchin+1) = chvite
                lpain(nchin+2) = 'PACCELE'
                lchin(nchin+2) = chacce
                nchin = nchin+2
            end if
        end if
        if (kcalc .eq. 'NON') then
            call getvid(' ', 'RESULTAT', scal=resu, nbret=iret)
            call rsexch(' ', resu, 'SIEF_ELGA', iord, chsig, &
                        iret)
            lpain(nchin+1) = 'PCONTGR'
            lchin(nchin+1) = chsig
            nchin = nchin+1
        end if
!
        call calcul('S', opti, modelLigrel, nchin, lchin, &
                    lpain, 1, lchout, lpaout, 'V', &
                    'OUI')
        call mesomm(chgthi, 1, vr=gthi(1))
        zr(iadrg+i-1) = gthi(1)
    end do
!
!- CALCUL DE G(S) SUR LE FOND DE FISSURE PAR 4 METHODES
!- PREMIERE METHODE : G_LEGENDRE ET THETA_LEGENDRE
!- DEUXIEME METHODE : G_LEGENDRE ET THETA_LAGRANGE
!- TROISIEME METHODE: G_LAGRANGE ET THETA_LAGRANGE
!    (OU G_LAGRANGE_NO_NO ET THETA_LAGRANGE)
!
    AS_ALLOCATE(vr=valg_s, size=nnoff)
    if ((liss .eq. 'LAGRANGE') .or. (liss .eq. 'LAGRANGE_NO_NO')) then
        call wkvect('&&MECAGL.VALGI', 'V V R8', nnoff, iadgi)
    else
        call wkvect('&&MECAGL.VALGI', 'V V R8', ndeg+1, iadgi)
    end if
! ABSCISSE CURVILIGNE
    call jeveuo(chfond, 'L', ifon)
    objcur = '&&MECAGL.ABSGAMM0'
    call wkvect(objcur, 'V V R', nnoff, iadabs)
    do i = 1, nnoff
        zr(iadabs-1+(i-1)+1) = zr(ifon-1+4*(i-1)+4)
    end do
    xl = zr(iadabs-1+(nnoff-1)+1)
!
! NOM DES NOEUDS DU FOND
    if (.not. lxfem) call jeveuo(fonoeu, 'L', iadrno)
!
    if ((liss .ne. 'LAGRANGE') .and. (liss .ne. 'LAGRANGE_NO_NO') .and. (liss .ne. 'MIXTE')) then
        num = 1
        call gmeth1(nnoff, ndeg, zr(iadrg), valg_s, objcur, &
                    xl, zr(iadgi))
    else if ((liss .eq. 'LAGRANGE') .or. (liss .eq. 'LAGRANGE_NO_NO') .or. (liss .eq. 'MIXTE')) then
        normff = zk24(jresu+nnoff+1-1)
        normff(20:24) = '.VALE'
        if ((liss .ne. 'LAGRANGE') .and. (liss .ne. 'LAGRANGE_NO_NO')) then
            num = 2
            call gmeth2(nnoff, ndeg, zr(iadrg), valg_s, &
                        objcur, xl, zr(iadgi), norfon)
!
        else
            call gmeth3(nnoff, zr(iadrg), milieu, valg_s, &
                        objcur, zr(iadgi), num, connex)
        end if
    end if
!
!- SYMETRIE DU CHARGEMENT ET IMPRESSION DES RESULTATS
!
    if (symech .ne. 'NON') then
        do i = 1, nnoff
            valg_s(i) = 2.d0*valg_s(i)
        end do
    end if
!
!- IMPRESSION ET ECRITURE DANS TABLE(S) DE G(S)
!
    if (niv .ge. 2) then
        call gimpgs(result, nnoff, zr(iadabs), valg_s, num, &
                    zr(iadgi), ndeg, ndimte, zr(iadrg), extim, &
                    time, iord, ifm)
    end if
!
    call getvis('THETA', 'NUME_FOND', iocc=1, scal=numfon, nbret=ibid)
!
    call tbajvi(result, nbprup, 'NUME_FOND', numfon, livi)
!
    call tbajvi(result, nbprup, 'NUME_ORDRE', iord, livi)
    call tbajvr(result, nbprup, 'INST', time, livr)
!
    do i = 1, nnoff
        if (.not. lxfem) then
            call tbajvk(result, nbprup, 'NOEUD', zk8(iadrno+i-1), livk)
        end if
        call tbajvi(result, nbprup, 'NUM_PT', i, livi)
        call tbajvr(result, nbprup, 'ABSC_CURV', zr(coor-1+4*(i-1)+4), livr)
        call tbajvr(result, nbprup, 'COOR_X', zr(coor-1+4*(i-1)+1), livr)
        call tbajvr(result, nbprup, 'COOR_Y', zr(coor-1+4*(i-1)+2), livr)
        call tbajvr(result, nbprup, 'COOR_Z', zr(coor-1+4*(i-1)+3), livr)
        call tbajvr(result, nbprup, 'G', valg_s(i), livr)
        call tbajli(result, nbprup, noprup, livi, livr, &
                    livc, livk, 0)
    end do
!
!- DESTRUCTION D'OBJETS DE TRAVAIL
!
    call jedetr(objcur)
    AS_DEALLOCATE(vr=valg_s)
    call jedetr('&&MECAGL.VALGI')
    call detrsd('CHAMP_GD', chvarc)
    call detrsd('CHAMP_GD', chvref)
    call detrsd('CHAMP_GD', chvolu)
    call detrsd('CHAMP_GD', cf1d2d)
    call detrsd('CHAMP_GD', cf2d3d)
    call detrsd('CHAMP_GD', chpres)
    call detrsd('CHAMP_GD', chepsi)
    call detrsd('CHAMP_GD', chpesa)
    call detrsd('CHAMP_GD', chrota)
    call jedetr('&&MECAGL.VALG')
!
    call jedema()
end subroutine
