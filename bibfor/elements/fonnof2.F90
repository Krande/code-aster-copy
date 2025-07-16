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

subroutine fonnof2(resu, noma, typfon, nbnoff, basnof)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cgnop0.h"
#include "asterfort/dfflon.h"
#include "asterfort/dfftan.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/gmgnre.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/oreino.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8)  :: noma, resu, typfon
    character(len=19) :: basnof
    integer(kind=8)           :: nbnoff
! FONCTION REALISEE (OPERATEUR DEFI_FOND_FISS) :
!
!      REMPLISSAGE DES OBJETS .SUPNORM.NOEU ET .INFNORM.NOEU
!
!     ENTREES:
!        RESU       : NOM DU CONCEPT RESULTAT DE L'OPERATEUR
!        NOMA       : NOM DU MAILLAGE
!        TYPFON     : TYPE DE FOND
!                     IL PEUT VALOIR OUVERT/FERME
!        NBNOFF     : NOMBRE DE NOEUDS EN FOND DE FISSURE
!        BASNOF     : BASE LOCALE EN FOND DE FISSURE
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nbma, jlima, im, n1, nbnoe, iret, jsup, jnols, inols
    integer(kind=8) :: idlino, nbnols, jcoors, in, nbnoft, inoff, ifm, niv, jnofo
    integer(kind=8) :: nuno, jints, nbtrls, iera, inols2, inoli2
    integer(kind=8) :: numfin, isym, inoli, jinf, jnoli, nbnoli, jcoori, jinti
    integer(kind=8) :: nbtrli, ino, inos, nbs, numun, jts, jti, nbi, inoi
    integer(kind=8) :: jnofos, irlev, ndim
    integer(kind=8) :: ninfsup_norm, ninfsup_norm2
    real(kind=8) :: x0(3), x1, x2, y1, y2, z1, z2, d, vplan(3), dmin
    real(kind=8) :: dmax, prec, preco, ps, vectan(3), precn
    character(len=6) :: nompro
    character(len=8) :: critn
    character(len=24) :: msup, minf, fonnoe, nomnoe
    real(kind=8), pointer :: basefond(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
! DEB-------------------------------------------------------------------
!
    call jemarq()
!   CHANGER LE NOM DANS LA RESTITUTION FINALE
    nompro = 'FONOF2'
    critn = 'RELATIF'
    call infniv(ifm, niv)
!   INDICATEUR DE COMMANDE POUR OREINO: 3-DEFI_FOND_FISS/PREC_NORM
    iera = 3
    jnofos = 0

    ninfsup_norm = 20
    ninfsup_norm2 = 100
!
!     ------------------------------------------------------------------
!                        LE MAILLAGE, LE FOND
!     ------------------------------------------------------------------
!
    call getvid(' ', 'MAILLAGE', scal=noma, nbret=n1)
!
    call getvr8(' ', 'PREC_NORM', scal=precn, nbret=n1)
!
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoe)
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
    nomnoe = noma//'.NOMNOE'
    fonnoe = resu//'.FOND.NOEU'
    call jeexin(fonnoe, irlev)
    if (irlev .ne. 0) then
        call jeveuo(fonnoe, 'L', jnofo)
    else
        ASSERT(.FALSE.)
!        fonnoe =resu//'.FOND_INF.NOEU'
!        call jeveuo(fonnoe, 'L', jnofo)
!        fonnoe =resu//'.FOND_SUP.NOEU'
!        call jeveuo(fonnoe, 'L', jnofos)
    end if
!
!     BASE LOCALE EN FOND DE FISSURE
    call jeveuo(basnof, 'L', vr=basefond)
!
!     ------------------------------------------------------------------
!                  VECTEUR RESULTAT
!     ------------------------------------------------------------------
!
    call wkvect(resu//'.SUPNORM.NOEU', 'G V K8', nbnoff*ninfsup_norm, inols)
    call wkvect(resu//'.SUPNORM.NOEU2', 'G V K8', nbnoff*ninfsup_norm2, inols2)
!
    call jeexin(resu//'.LEVREINF.MAIL', isym)
    if (isym .ne. 0) then
        call wkvect(resu//'.INFNORM.NOEU', 'G V K8', nbnoff*ninfsup_norm, inoli)
        call wkvect(resu//'.INFNORM.NOEU2', 'G V K8', nbnoff*ninfsup_norm2, inoli2)
    end if
!
!     ------------------------------------------------------------------
!         GROUP_MA LEVRE_SUP --> GROUP_NO LEVRE_SUP
!     ------------------------------------------------------------------
    msup = resu//'.LEVRESUP.MAIL'
    call jeveuo(msup, 'L', jsup)
    call jelira(msup, 'LONMAX', nbma)
!
    call wkvect('&&'//nompro//'_TRAV', 'V V I', nbnoe, idlino)
    call wkvect('&&'//nompro//'_NOEU_NORM_SUP', 'V V I', nbnoe, jnols)
!
    call wkvect('&&'//nompro//'_MAILLE_LEV_SUP', 'V V I', nbma, jlima)
    do im = 1, nbma
        zi(jlima-1+im) = char8_to_int(zk8(jsup-1+im))
    end do
    call gmgnre(noma, nbnoe, zi(idlino), zi(jlima), nbma, &
                zi(jnols), nbnols, 'TOUS')
!
    call wkvect('&&PKFOND_COOR_LEV_SUP', 'V V R', 3*nbnols, jcoors)
    do in = 1, nbnols
        zr(jcoors-1+3*(in-1)+1) = vale(3*(zi(jnols-1+in)-1)+1)
        zr(jcoors-1+3*(in-1)+2) = vale(3*(zi(jnols-1+in)-1)+2)
        zr(jcoors-1+3*(in-1)+3) = vale(3*(zi(jnols-1+in)-1)+3)
    end do
!
    call jedetr('&&'//nompro//'_TRAV')
!
!     ------------------------------------------------------------------
!         GROUP_MA LEVRE_INF --> GROUP_NO LEVRE_INF
!     ------------------------------------------------------------------
    if (isym .ne. 0) then
        minf = resu//'.LEVREINF.MAIL'
        call jeveuo(minf, 'L', jinf)
        call jelira(minf, 'LONMAX', nbma)
!
        call wkvect('&&'//nompro//'_TRAV', 'V V I', nbnoe, idlino)
        call wkvect('&&'//nompro//'_NOEU_NORM_INF', 'V V I', nbnoe, jnoli)
        call wkvect('&&'//nompro//'_MAILLE_LEV_INF', 'V V I', nbma, jlima)
        do im = 1, nbma
            zi(jlima-1+im) = char8_to_int(zk8(jinf-1+im))
        end do
        call gmgnre(noma, nbnoe, zi(idlino), zi(jlima), nbma, &
                    zi(jnoli), nbnoli, 'TOUS')
!
        call wkvect('&&PKFOND_COOR_LEV_INF', 'V V R', 3*nbnoli, jcoori)
        do in = 1, nbnoli
            zr(jcoori-1+3*(in-1)+1) = vale(3*(zi(jnoli-1+in)-1)+1)
            zr(jcoori-1+3*(in-1)+2) = vale(3*(zi(jnoli-1+in)-1)+2)
            zr(jcoori-1+3*(in-1)+3) = vale(3*(zi(jnoli-1+in)-1)+3)
        end do
        call jedetr('&&'//nompro//'_TRAV')
    end if
!
!     ------------------------------------------------------------------
!        BOUCLE SUR LES NOEUDS DU FOND DE FISSURE
!     ------------------------------------------------------------------
!
!
! ------ CAS DU FOND_FERME: LE PREMIER ET LE DERNIER NOEUD SONT
!                           IDENTIQUES
    if (typfon .eq. 'FERME') then
        nbnoft = nbnoff-1
    else
        nbnoft = nbnoff
    end if
!
!
    do inoff = 1, nbnoft
!
!        DETERMINATION DU PLAN ORTHOGONAL AU FOND DE FISSURE
!        ET PASSANT PAR LE NOEUD COURANT
!
        nuno = char8_to_int(zk8(jnofo-1+inoff))
!
        x0(1) = vale(3*(nuno-1)+1)
        x0(2) = vale(3*(nuno-1)+2)
        x0(3) = vale(3*(nuno-1)+3)
!
!        VPLAN = VECTEUR ORTHOGONAL AU PLAN RECHERCHE
        call dfftan(ndim, basefond, inoff, vplan)
!
!        D = LONGUEUR DU SEGMENT DU FOND DE FISSURE
        if (ndim .eq. 3) then
            call dfflon(vale, zk8(jnofo), nomnoe, inoff, nbnoff, &
                        typfon, d)
        else if (ndim .eq. 2) then
!          D N'A PAS DE SENS EN 2D, ON PREND ALORS UNE VALEUR CARACT
            d = sqrt(x0(1)**2+x0(2)**2)
        end if
!
!        PREC = PRECISION POUR LA RECHERCHE DES NOEUDS DES LEVRES
        prec = d*precn
!
! ------ CALCUL DE L'INTERSECTION DU PLAN AVEC LES NOEUDS DES MAILLES
!        DEFINISSANT LA LEVRE INFERIEURE - SAUVEGARDE
!
        if (isym .ne. 0) then
            call wkvect('&&PKFOND_INTERS_INF', 'V V I', nbnoe, jinti)
            call cgnop0(nbnoli, zr(jcoori), x0, vplan, prec, &
                        nbtrli, zi(jinti))
            if (nbtrli .le. 2) then
                call jedetr('&&PKFOND_INTERS_INF')
                goto 200
            end if
            call wkvect('&&PKFOND_TRAV_INF', 'V V I', nbtrli, jti)
!
! ---- ORDRE DES NOEUDS
            numfin = nuno
            dmax = 0.d0
            dmin = 100.d0
            nbi = 1
            inoi = 1
            zi(jti-1+inoi) = nuno
            x1 = vale(3*(nuno-1)+1)
            y1 = vale(3*(nuno-1)+2)
            z1 = vale(3*(nuno-1)+3)
!
! identification noeuds sur bon cote de la levre (cas fond ferme)
! a partir du noeud le plus proche du fond
            do in = 1, nbtrli
                ino = jnoli+zi(jinti-1+in)-1
                if (zi(ino) .eq. nuno) goto 310
                x2 = vale(3*(zi(ino)-1)+1)
                y2 = vale(3*(zi(ino)-1)+2)
                z2 = vale(3*(zi(ino)-1)+3)
                d = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                if (d .lt. dmin) then
                    dmin = d
                    numun = zi(ino)
                end if
310             continue
            end do

!
            vectan(1) = vale(3*(numun-1)+1)-x1
            vectan(2) = vale(3*(numun-1)+2)-y1
            vectan(3) = vale(3*(numun-1)+3)-z1
!
            do in = 1, nbtrli
                ino = jnoli+zi(jinti-1+in)-1
                if (zi(ino) .eq. nuno) goto 320
                x2 = vale(3*(zi(ino)-1)+1)
                y2 = vale(3*(zi(ino)-1)+2)
                z2 = vale(3*(zi(ino)-1)+3)
                d = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                ps = (x2-x1)*vectan(1)+(y2-y1)*vectan(2)+(z2-z1)*vectan( &
                     3)
                if (ps .ge. 0.d0) then
                    nbi = nbi+1
                    inoi = inoi+1
                    zi(jti-1+inoi) = zi(ino)
                    if (d .gt. dmax) then
                        dmax = d
                        numfin = zi(ino)
                    end if
                end if
320             continue
            end do
!
            preco = prec*10
            call oreino(noma, zi(jti), nbi, nuno, numfin, &
                        vale, critn, preco, iera, iret)

            do in = 1, min(nbi, ninfsup_norm)
                zk8(inoli-1+ninfsup_norm*(inoff-1)+in) = int_to_char8(zi(jti-1+in))
            end do
            do in = 1, min(nbi, ninfsup_norm2)
                zk8(inoli2-1+ninfsup_norm2*(inoff-1)+in) = int_to_char8(zi(jti-1+in))
            end do
!
            call jedetr('&&PKFOND_INTERS_INF')
            call jedetr('&&PKFOND_TRAV_INF')
        end if
!
! ------ CALCUL DE L'INTERSECTION DU PLAN AVEC LES NOEUDS DES MAILLES
!        DEFINISSANT LA LEVRE SUPERIEURE - SAUVEGARDE
!
        call wkvect('&&PKFOND_INTERS_SUP', 'V V I', nbnoe, jints)
        call cgnop0(nbnols, zr(jcoors), x0, vplan, prec, &
                    nbtrls, zi(jints))
        if (nbtrls .le. 2) then
            call jedetr('&&PKFOND_INTERS_SUP')
            goto 200
        end if
        call wkvect('&&PKFOND_TRAV_SUP', 'V V I', nbtrls, jts)
!
! ---- ORDRE DES NOEUDS
        if (irlev .eq. 0) then
            nuno = char8_to_int(zk8(jnofos-1+inoff))
        end if
        numfin = nuno
        dmax = 0.d0
        dmin = 100.d0
        nbs = 1
        inos = 1
        zi(jts-1+inos) = nuno
        x1 = vale(3*(nuno-1)+1)
        y1 = vale(3*(nuno-1)+2)
        z1 = vale(3*(nuno-1)+3)
!
! identification noeuds sur bon cote de la levre (cas fond ferme)
! a partir du noeud le plus proche du fond
        do in = 1, nbtrls
            ino = jnols+zi(jints-1+in)-1
            if (zi(ino) .eq. nuno) goto 210
            x2 = vale(3*(zi(ino)-1)+1)
            y2 = vale(3*(zi(ino)-1)+2)
            z2 = vale(3*(zi(ino)-1)+3)
            d = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
            if (d .lt. dmin) then
                dmin = d
                numun = zi(ino)
            end if
210         continue
        end do
!
        vectan(1) = vale(3*(numun-1)+1)-x1
        vectan(2) = vale(3*(numun-1)+2)-y1
        vectan(3) = vale(3*(numun-1)+3)-z1
!
        do in = 1, nbtrls
            ino = jnols+zi(jints-1+in)-1
            if (zi(ino) .eq. nuno) goto 220
            x2 = vale(3*(zi(ino)-1)+1)
            y2 = vale(3*(zi(ino)-1)+2)
            z2 = vale(3*(zi(ino)-1)+3)
            d = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
            ps = (x2-x1)*vectan(1)+(y2-y1)*vectan(2)+(z2-z1)*vectan(3)
            if (ps .ge. 0.d0) then
                nbs = nbs+1
                inos = inos+1
                zi(jts-1+inos) = zi(ino)
                if (d .gt. dmax) then
                    dmax = d
                    numfin = zi(ino)
                end if
            end if
220         continue
        end do
        preco = prec*10
        call oreino(noma, zi(jts), nbs, nuno, numfin, &
                    vale, critn, preco, iera, iret)
!
        do in = 1, min(nbs, ninfsup_norm)
            zk8(inols-1+ninfsup_norm*(inoff-1)+in) = int_to_char8(zi(jts-1+in))
        end do
        do in = 1, min(nbs, ninfsup_norm2)
            zk8(inols2-1+ninfsup_norm2*(inoff-1)+in) = int_to_char8(zi(jts-1+in))
        end do
!
        call jedetr('&&PKFOND_INTERS_SUP')
        call jedetr('&&PKFOND_TRAV_SUP')
!
200     continue
    end do
!
    call jedema()
end subroutine
