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
!
subroutine xgelem(elrefp, ndim, coorse, igeom, jheavt, &
                  ise, nfh, ddlc, ddlm, nfe, &
                  basloc, nnop, idepl, lsn, lst, &
                  igthet, fno, nfiss, jheavn, jstno, incr)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/nmplru.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmelnl.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvarc.h"
#include "asterfort/reeref.h"
#include "asterfort/tecach.h"
#include "asterfort/xcinem.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xkamat.h"
#include "asterfort/xnbddl.h"
!
    character(len=8) :: elrefp
    integer(kind=8) :: igeom, ndim, nfh, ddlc, nfe, nnop, ddlm, jheavn
    integer(kind=8) :: idepl, nfiss, jheavt, ise, jstno
    real(kind=8) :: basloc(3*ndim*nnop), lsn(nnop), lst(nnop)
    real(kind=8) :: fno(ndim*nnop), coorse(*)
!
    aster_logical :: incr
!
!    - FONCTION REALISEE:  CALCUL DU TAUX DE RESTITUTION D'ENERGIE
!                          PAR LA METHODE ENERGETIQUE G-THETA
!                          POUR LES ELEMENTS X-FEM (2D ET 3D)
!
! IN  ELREFP  : ÉLÉMENT DE RÉFÉRENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNÉES DES SOMMETS DU SOUS-ÉLÉMENT
! IN  IGEOM   : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ÉLT
! IN  NFH     : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  NFISS   : NOMBRE DE FISSURES "VUES" PAR L'ÉLÉMENT
! IN  JHEAVN  : POINTEUR VERS LA DEFINITION HEAVISIDE
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE AUX NOEUDS
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  DEPL    : DÉPLACEMENTS
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  FNO     : FORCES VOLUMIQUES AUX NOEUDS DE L'ELEMENT PARENT
! OUT IGTHET  : G
!
!
    integer(kind=8) :: ithet, imate, icomp, igthet, jtab(7), ncomp, idecpg, idebs
    integer(kind=8) :: ipoids, jcoopg, ivf, idfde, jdfd2, jgano, jsigse
    integer(kind=8) :: i, j, k, kpg, n, ino, iret, cpt, ig, in, mxstac, isigi, isigm
    integer(kind=8) :: ndimb, nno, nnos, npgbis, ddld, ddls, matcod, m, irett
    integer(kind=8) :: ncompn, heavn(nnop, 5), hea_se
    integer(kind=8) :: singu, alp
    real(kind=8) :: xg(ndim), he(nfiss)
    real(kind=8) :: xe(ndim), ff(nnop), dfdi(nnop, ndim), f(3, 3)
    real(kind=8) :: eps(6)
    real(kind=8) :: grad(ndim, ndim), dudm(3, 4), poids
    real(kind=8) :: dtdm(3, 4)
    real(kind=8) :: rbid
    real(kind=8) :: tthe, r, rp, ppg
    real(kind=8) :: depla(3), theta(3), tgudm(3), tpn(27), tref
    real(kind=8) :: dfdm(3, 4), dfdx(27), dfdy(27), dfdz(27)
    real(kind=8) :: dtx, dty, dtz
    real(kind=8) :: energi(2), sigl(6), prod, prod2, rac2, sr(3, 3), tcla, divt
    real(kind=8) :: tfor, sigse(6*27)
    real(kind=8) :: fk(27, 3, 3), dkdgl(27, 3, 3, 3), ka, mu2
    character(len=8) :: elrese(6), fami(6), typmod(2)
    character(len=16) :: compor(4)
    aster_logical :: cp, axi, l_temp_noeu
    integer(kind=8) :: irese, ddli, nnoi, indeni, nnops, ifiss
    integer(kind=8) :: iret1, iret2, iret3
    type(Behaviour_Integ) :: BEHinteg
!
    real(kind=8) :: tini, prod1, dsigin(6, 3), sigin(6), epsref(6), epsp(6)
    real(kind=8) :: mu, nu(1), e(1)
    integer(kind=8) ::  icodre(1), ncmp
    character(len=16) :: phenom
!
!
!
    parameter(mxstac=1000)
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!
! - Initialisation of behaviour datastructure
!
    call behaviourInit(BEHinteg)
!
!     VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!     (VOIR CRS 1404)
    ASSERT(ndim .le. mxstac)
    ASSERT(nnop .le. mxstac)
    ASSERT(nfiss .le. mxstac)
!
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
!
    typmod = ' '
    sigl = 0.d0
    cp = .false.
    rac2 = sqrt(2.d0)
    tcla = 0.d0
    tthe = 0.d0
    tfor = 0.d0
    tini = 0.d0
    rbid = 0.d0
!
    if (lteatt('C_PLAN', 'OUI')) then
        typmod(1) = 'C_PLAN'
        cp = .true.
    else if (lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN'
    end if
!
!   NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, ddls, singu)
!
!   NOMBRE DE COMPOSANTES DE PHEAVTO (DANS LE CATALOGUE)
    call tecach('OOO', 'PHEAVTO', 'L', iret, nval=2, &
                itab=jtab)
    ncomp = jtab(2)
!
!   ELEMENT DE REFERENCE PARENT : RECUP DE NNOPS
    call elrefe_info(fami='RIGI', nnos=nnops)
!
    axi = lteatt('AXIS', 'OUI')
!
!   NOMBRE DE COMPOSANTES DES TENSEURS
    ncmp = 2*ndim
!
    call jevech('PTHETAR', 'L', ithet)
    call jevech('PMATERC', 'L', imate)
    matcod = zi(imate)
    call jevech('PCOMPOR', 'L', icomp)
    do i = 1, 4
        compor(i) = zk16(icomp+i-1)
    end do
!
!   SOUS-ELEMENT DE REFERENCE
    call elrefe_info(elrefe=elrese(ndim+irese), &
                     fami=fami(ndim+irese), &
                     ndim=ndimb, &
                     nno=nno, &
                     nnos=nnos, &
                     npg=npgbis, &
                     jpoids=ipoids, &
                     jcoopg=jcoopg, &
                     jvf=ivf, &
                     jdfde=idfde, &
                     jdfd2=jdfd2, &
                     jgano=jgano)

    ASSERT(ndim .eq. ndimb)

    if (incr) then
        call jevech('PCONTRR', 'L', isigm)
    end if

!   Recuperation de la contrainte initiale aux noeuds des sous-elts
    call tecach('ONO', 'PSIGISE', 'L', iret, iad=jsigse)

!   Indicateur de contrainte initiale
    isigi = 0
    if (jsigse .ne. 0) isigi = 1

    if (isigi .ne. 0) then
!       Passage de la contrainte initiale aux noeuds des sous-elts
!       dans un tableau local au sous-elt
        do i = 1, nno
            do j = 1, ncmp
                sigse(ncmp*(i-1)+j) = &
                    zr(jsigse-1+ncmp*nno*(ise-1)+ncmp*(i-1)+j)
            end do
        end do

    end if
!
!   TEMPERATURE DE REF
    call rcvarc(' ', 'TEMP', 'REF', 'XFEM', 1, &
                1, tref, irett)
    if (irett .ne. 0) tref = 0.d0
!
!   TEMPERATURE AUX NOEUDS PARENT
    l_temp_noeu = .false.
    do ino = 1, nnop
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, &
                    1, tpn(ino), iret)
        if (iret .ne. 0) tpn(ino) = 0.d0
    end do
    if (iret .eq. 0) l_temp_noeu = .true.
!
!   FONCTION HEAVYSIDE CSTE SUR LE SS-ELT ET PAR FISSURE
    do ifiss = 1, nfiss
        he(ifiss) = zi(jheavt-1+ncomp*(ifiss-1)+ise)
    end do
!
!   RECUPERATION DE LA DEFINITION DES FONCTIONS HEAVISIDES
    if (nfh .gt. 0) then
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
        ASSERT(ncompn .eq. 5)
        do ino = 1, nnop
            do ifiss = 1, ncompn
                heavn(ino, ifiss) = zi(jheavn-1+ncompn*(ino-1)+ifiss)
            end do
        end do
    end if
!
! CALCUL DE L IDENTIFIANT DU SS ELEMENT
    hea_se = xcalc_code(nfiss, he_real=[he])
!     ------------------------------------------------------------------
!     BOUCLE SUR LES POINTS DE GAUSS DU SOUS-TETRA
!     ------------------------------------------------------------------
!
!   indice de decalage pour les champs ELGA
    idecpg = (ise-1)*npgbis
    idebs = ncmp*idecpg
!
    do kpg = 1, npgbis
!
!       INITIALISATIONS
        dtdm(:, :) = 0.d0
        dudm(:, :) = 0.d0
        do i = 1, 6
            sigin(i) = 0.d0
            epsref(i) = 0.d0
            epsp(i) = 0.d0
            do j = 1, 3
                dsigin(i, j) = 0.d0
            end do
        end do
!
!
!       COORDONNEES DU PT DE GAUSS DANS LE REPERE REEL : XG
        xg(:) = 0.d0
        do i = 1, ndim
            do n = 1, nno
                xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+i)
            end do
        end do
!
!       CALCUL DES FF
        call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                    xe, ff, dfdi=dfdi)
!
!       POUR CALCULER LE JACOBIEN DE LA TRANSFO SS-ELT -> SS-ELT REF
!       AINSI QUE LES DERIVEES DES FONCTIONS DE FORMES DU SS-ELT
!       ON ENVOIE DFDM3D/DFDM2D AVEC LES COORD DU SS-ELT
        if (ndim .eq. 3) call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                                     poids, dfdx, dfdy, dfdz)
        if (ndim .eq. 2) call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                                     poids, dfdx, dfdy)
!
!
! -     CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE)
        if (axi) then
            r = 0.d0
            do ino = 1, nnop
                r = r+ff(ino)*zr(igeom-1+2*(ino-1)+1)
            end do
!
            poids = poids*r
            ASSERT(r .gt. 0d0)
!
        end if
!
!         FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
!
!       FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
        if (singu .gt. 0) then
            if (isigi .eq. 0) then
                call xkamat(matcod, ndim, axi, ka, mu2)
            else
                call rccoma(matcod, 'ELAS', 1, phenom, icodre(1))
                call rcvala(matcod, ' ', phenom, 1, ' ', &
                            [rbid], 1, 'NU', nu(1), icodre(1), 1)
                call rcvala(matcod, ' ', phenom, 1, ' ', &
                            [rbid], 1, 'E', e(1), icodre(1), 1)
                mu2 = e(1)/(2.d0*(1.d0+nu(1)))
                ka = 3.d0-4.d0*nu(1)
                if (lteatt('C_PLAN', 'OUI')) ka = (3.d0-nu(1))/(1.d0+nu(1))
            end if
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(1), &
                              lsn, lst, zr(igeom), ka, mu2, ff, fk, dfdi, dkdgl)
        end if
!
!       ---------------------------------------------
!       2) CALCUL DU DEPLACEMENT ET DE SA DERIVEE (DUDM)
!       ---------------------------------------------
!
        depla(:) = 0.d0
!
!       CALCUL DE L'APPROXIMATION DU DEPLACEMENT
        do in = 1, nnop
            if (in .le. nnops) then
                nnoi = 0
                ddli = ddls
            else if (in .gt. nnops) then
                nnoi = nnops
                ddli = ddlm
            end if
            indeni = ddls*nnoi+ddli*(in-nnoi-1)
!
            cpt = 0
!         DDLS CLASSIQUES
            do i = 1, ndim
                cpt = cpt+1
                depla(i) = depla(i)+ff(in)*zr(idepl-1+indeni+cpt)
            end do
!         DDLS HEAVISIDE
            do ig = 1, nfh
                do i = 1, ndim
                    cpt = cpt+1
                    depla(i) = depla(i)+xcalc_heav(heavn(in, ig), hea_se, heavn(in, 5)) &
                               *ff(in)*zr(idepl-1+indeni+cpt)
                end do
            end do
!         DDL ENRICHIS EN FOND DE FISSURE
            do alp = 1, ndim*singu
                cpt = cpt+1
                do i = 1, ndim
                    depla(i) = depla(i)+fk(in, alp, i)*zr(idepl-1+indeni+cpt)
                end do
            end do
!
        end do
!
!       CALCUL DU GRAD DE U AU POINT DE GAUSS
        call xcinem(axi, igeom, nnop, nnos, idepl, &
                    ndim, he, &
                    nfiss, nfh, singu, ddls, ddlm, &
                    fk, dkdgl, ff, dfdi, f, &
                    eps, grad, heavn)
!
!       ON RECOPIE GRAD DANS DUDM (CAR PB DE DIMENSIONNEMENT SI 2D)
        do i = 1, ndim
            do j = 1, ndim
                dudm(i, j) = grad(i, j)
            end do
        end do
!
!       VALEUR DU DEPLACEMENT DANS LA QUATRIEME COLONNE :
        do i = 1, ndim
            dudm(i, 4) = depla(i)
        end do
!
!       TRAITEMENTS DEPENDANT DE LA MODELISATION
        if (cp) then
            dudm(3, 3) = eps(3)
        end if
        if (axi) then
            dudm(3, 3) = dudm(1, 4)/r
        end if
!
!       ------------------------------------------------
!       3) CALCUL DU CHAMP THETA ET DE SA DERIVEE (DTDM)
!       ------------------------------------------------
!
        do i = 1, ndim
!
            theta(i) = 0.d0
            do ino = 1, nnop
                theta(i) = theta(i)+ff(ino)*zr(ithet-1+ndim*(ino-1)+i)
            end do
!
            do j = 1, ndim
                do ino = 1, nnop
                    dtdm(i, j) = dtdm(i, j)+zr(ithet-1+ndim*(ino-1)+i)*dfdi(ino, j)
                end do
            end do
        end do
!
!       valeur du champ theta dans la quatrieme colonne
        do i = 1, ndim
            dtdm(i, 4) = theta(i)
        end do
!
        if (axi) then
            dtdm(3, 3) = dtdm(1, 4)/r
        end if
!
        divt = 0.d0
        do i = 1, 3
            divt = divt+dtdm(i, i)
        end do
!
!       --------------------------------------------------
!       4) CALCUL DU CHAMP DE TEMPERATURE ET DE SA DERIVEE
!       --------------------------------------------------
!
        tgudm(:) = 0.d0
!
!       cas de la varc TEMP, "continue" et donnee au noeud. Calcul
!       de ses derivees partielles
        if (l_temp_noeu) then
            do i = 1, ndim
                do ino = 1, nnop
                    tgudm(i) = tgudm(i)+dfdi(ino, i)*tpn(ino)
                end do
            end do
        end if
!
!       cas des varc DTX DTY DTZ, derivees partielles de la temperature
!       "discontinue". Ces varc sont donnees aux pg xfem
        call rcvarc(' ', 'DTX', '+', 'XFEM', kpg+idecpg, 1, dtx, iret1)
        if (iret1 .eq. 0) then
!           economisons les appels a rcvarc... si DTX est absent, pas
!           besoin de recuperer les autres composantes
            ASSERT(.not. l_temp_noeu)
            call rcvarc(' ', 'DTY', '+', 'XFEM', kpg+idecpg, 1, dty, iret2)
            ASSERT(iret2 .eq. 0)
            tgudm(1) = dtx
            tgudm(2) = dty
            if (ndim .eq. 3) then
                call rcvarc(' ', 'DTZ', '+', 'XFEM', kpg+idecpg, 1, dtz, iret3)
                ASSERT(iret3 .eq. 0)
                tgudm(3) = dtz
            end if
        end if
!
!       --------------------------------------------------
!       5) CALCUL DE LA CONTRAINTE ET DE L ENERGIE
!       --------------------------------------------------
!
        if (incr) then

!           plasticite (en fait juste elasticite + comp_incr pour l'etat initial)

            ppg = 0.d0
            call nmplru('XFEM', kpg+idecpg, 1, '+', ndim, &
                        typmod, matcod, compor, ppg, eps, &
                        epsp, rp, energi)

            do i = 1, 3
                sigl(i) = zr(isigm+idebs-1+ncmp*(kpg-1)+i)
            end do
            sigl(4) = zr(isigm+idebs-1+ncmp*(kpg-1)+4)*rac2
            if (ndim .eq. 3) then
                sigl(5) = zr(isigm+idebs-1+ncmp*(kpg-1)+5)*rac2
                sigl(6) = zr(isigm+idebs-1+ncmp*(kpg-1)+6)*rac2
            end if

        else
            call nmelnl(BEHinteg, &
                        'XFEM', kpg+idecpg, 1, ndim, &
                        typmod, matcod, compor, &
                        eps, 0.d0, 0.d0, sigl, energi)

        end if

!
!       --------------------------------------------------
!       6)   CORRECTIONS LIEES A LA CONTRAINTE INITIALE
!       --------------------------------------------------
!
        if (isigi .ne. 0) then

!           Calcul de la contrainte initiale (somme sur les noeuds du ss-elt)
            do i = 1, nno
                do j = 1, ncmp
                    sigin(j) = sigin(j)+ &
                               sigse(ncmp*(i-1)+j)*zr(ivf-1+nno*(kpg-1)+i)
                end do
            end do
!
!           Calcul du gradient de sigma initial (somme sur les noeuds du ss-elt)
            do i = 1, nno
                do j = 1, ncmp
                    dsigin(j, 1) = dsigin(j, 1)+sigse(ncmp*(i-1)+j)*dfdx(i)
                    dsigin(j, 2) = dsigin(j, 2)+sigse(ncmp*(i-1)+j)*dfdy(i)
                    if (ndim .eq. 3) dsigin(j, 3) = dsigin(j, 3)+sigse(ncmp*(i-1)+j)*dfdz(i)*dfdz(i)
                end do
            end do
!
!           Traitements particuliers des termes croises
            do i = 4, ncmp
                sigin(i) = sigin(i)*rac2
                do j = 1, ndim
                    dsigin(i, j) = dsigin(i, j)*rac2
                end do
            end do
!
!           Calcul de la deformation de reference
            call rccoma(matcod, 'ELAS', 1, phenom, icodre(1))
            call rcvala(matcod, ' ', phenom, 1, ' ', &
                        [rbid], 1, 'NU', nu(1), icodre(1), 1)
            call rcvala(matcod, ' ', phenom, 1, ' ', &
                        [rbid], 1, 'E', e(1), icodre(1), 1)
!
            mu = e(1)/(2.d0*(1.d0+nu(1)))
!
            epsref(1) = -(1.d0/e(1))*(sigin(1)-(nu(1)*(sigin(2)+sigin(3))))
            epsref(2) = -(1.d0/e(1))*(sigin(2)-(nu(1)*(sigin(3)+sigin(1))))
            epsref(3) = -(1.d0/e(1))*(sigin(3)-(nu(1)*(sigin(1)+sigin(2))))
            epsref(4) = -(1.d0/mu)*sigin(4)
            if (ndim .eq. 3) then
                epsref(5) = -(1.d0/mu)*sigin(5)
                epsref(6) = -(1.d0/mu)*sigin(6)
            end if
!
!           Energie elastique (expression wadier)
            do i = 1, ncmp
                energi(1) = energi(1)+(eps(i)-0.5d0*epsref(i))*sigin(i)
            end do
!
        end if
!
!       -----------------------------------------------------------
!       7) CALCUL DES FORCES VOLUMIQUES ET DE LEURS DERIVEES (DFDM)
!       -----------------------------------------------------------
!
        dfdm(:, :) = 0.d0
        do ino = 1, nnop
            do j = 1, ndim
                do k = 1, ndim
                    dfdm(j, k) = dfdm(j, k)+fno(ndim*(ino-1)+j)*dfdi(ino, k)
                end do
!           VALEUR DE LA FORCE DANS LA QUATRIEME COLONNE :
                dfdm(j, 4) = dfdm(j, 4)+fno(ndim*(ino-1)+j)*ff(ino)
            end do
        end do
!
        if (axi) then
            dfdm(3, 3) = dfdm(1, 4)/r
        end if
!
!
!       --------------------------------------------------
!              TERME THERMOELAS. CLASSIQUE :
!       F.SIG:(GRAD(U).GRAD(THET))-ENER*DIVT
!       --------------------------------------------------
!
        sr(1, 1) = sigl(1)
        sr(2, 2) = sigl(2)
        sr(3, 3) = sigl(3)
        sr(1, 2) = sigl(4)/rac2
        sr(2, 1) = sr(1, 2)
        sr(1, 3) = sigl(5)/rac2
        sr(3, 1) = sr(1, 3)
        sr(2, 3) = sigl(6)/rac2
        sr(3, 2) = sr(2, 3)
!
        prod = 0.d0
        prod2 = 0.d0
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    do m = 1, 3
                        prod = prod+f(i, j)*sr(j, k)*dudm(i, m)*dtdm(m, k)
                    end do
                end do
            end do
        end do
!
        prod2 = poids*(prod-energi(1)*divt)
!
        tcla = tcla+prod2
!
!
!       =======================================================
!       TERME THERMIQUE :   -(D(ENER)/DT)(GRAD(T).THETA)
!       =======================================================
        if (irett .eq. 0) then
            prod = 0.d0
            prod2 = 0.d0
            do i = 1, ndim
                prod = prod+tgudm(i)*theta(i)
            end do
            prod2 = -poids*prod*energi(2)
            tthe = tthe+prod2
        end if
!
!       =======================================================
!       TERME FORCE VOLUMIQUE
!       =======================================================
!
        do i = 1, ndim
            prod = 0.d0
            do j = 1, ndim
                prod = prod+dfdm(i, j)*dtdm(j, 4)
            end do
            tfor = tfor+dudm(i, 4)*(prod+dfdm(i, 4)*divt)*poids
        end do
!
!       =======================================================
!       TERME CONTRAINTE INITIALE SIGIN
!       =======================================================
!
        if (isigi .ne. 0) then
            prod1 = 0.d0
            do i = 1, ncmp
                do j = 1, ndim
                    prod1 = prod1-(eps(i)-epsref(i))*dsigin(i, j)*dtdm(j, 4)
                end do
            end do
            tini = tini+prod1*poids
        end if
    end do
!
!   ------------------------------------------------------------------
!   FIN DE LA BOUCLE SUR LES POINTS DE GAUSS DU SOUS-TETRA
!   ------------------------------------------------------------------
!
    zr(igthet) = zr(igthet)+tcla+tthe+tfor+tini

!
end subroutine
