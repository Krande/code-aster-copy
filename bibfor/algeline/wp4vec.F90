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
subroutine wp4vec(nbfreq, nbvect, neq, shift, vp, &
                  vecp, mxresf, resufi, resufr, lagr, &
                  vauc, omecor)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/wpordc.h"
#include "asterfort/wptest.h"
    integer(kind=8) :: mxresf, neq, nbfreq, nbvect, lagr(*), resufi(mxresf, *)
    complex(kind=8) :: vecp(neq, *), shift, vauc(2*neq, *), vp(*)
    real(kind=8) :: resufr(mxresf, *), omecor
!     RESTITUTION DES VALEURS PROPRES ET DES MODES DU PB QUADRATIQUE
!     AVEC TRI SUIVANT LES PARTIES IMAGINAIRES CROISANTES
!     ET NORMALISATION A LA PLUS GRANDE CMP NON LAGRANGE
!
!     -----------------------------------------------------------------
! IN  NBFREQ : I : NOMBRE DE MODES DEMANDES
! IN  NBVECT : I : NOMBRE DE VECTEURS DE LANCZOS
! IN  NEQ    : I : TAILLE DES MATRICES DU PB QUADRATIQUE
! IN  SHIFT  : C : VALEUR DU DECALAGE
! IN  LAGR   : I : INDICATEUR DES NON-LAGRANGE
! VAR VP     : C : VALEURS PROPRE DU PB QUADRATIQUE
! OUT VECP   : C : MODES DU PB QUADRATIQUE
! IN  VAUC   : C : MODES DU PB QUADRATIQUE COMPLET
! OUT RESUFR : C : TABLEAU DE POST-TRAITEMENT
! IN  OMECOR : R : "ZERO MODAL", SEUIL EN DECA DUQUEL DEUX MODES SONT
!                  CONSIDERES COMME IDENTIQUES
!     -----------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
    real(kind=8) :: am, om, nmabp, seuilp, seuilr, c1, auxrj, auxij, auxrk
    real(kind=8) :: auxik, seuilc
    integer(kind=8) :: i, j, k, av1, av2, iadind, nbreel, nbcmpp, nbcmpc, nbfrga
    integer(kind=8) :: vali(5), ifm, niv, nbfr, ibid
    aster_logical :: trouve, lconj
    character(len=1) :: kmsg
    character(len=16) :: valk, typres
!
!     -----------------------------------------------------------------
    call jemarq()
    call infniv(ifm, niv)
!
! --- 1. PARTITION ,TEST DE CONJUGAISON, ELIMINATION DES CONJUGUES
!        REMARQUE : SI QR N' A PAS EU DE PB ALORS LES MODES REDUITS :
!                   * COMPLEXES ET 2 A 2 CONJUGUES,LES CONJUGUES
!                     APPARAISSENT LES UNS A LA SUITE DES AUTRES
!                   * REELS
!
!        IMPLEMENTATION : DANS UN TABLEAU D' ENTIER (ZI(IADIND))
!        INVARIANT      : T(I) = -2 VP(I) NON CLASSEE
!                         T(I) =  1 <=> VP(I) COMPLEXE AVEC CONJUGUEE
!                                             SELECTIONNEE
!                         T(I) = -1 <=> VP(I) COMPLEXE AVEC CONJUGUEE
!                                             ELIMINEE
!                         T(I) =  0 <=> VP(I) COMPLEXE SANS CONJUGUEE
!                                             OU REELLE
!                                             ELIMINEE
!     -----------------------------------------------------------------
! --- 1.1. PARTITION (OPERATEUR REEL)
    nbcmpp = 0
    nbcmpc = 0
    nbreel = 0
!*****************************************************************
!     SI IM(VP)<SEUILR: VP EST CONSIDEREE COMME REELLE
    seuilr = 1.d-7
!     SI MODULE(VPK-VPJ) < SEUILP: VPK = CONJUGEE DE VPJ
    seuilp = omecor
!     SEUIL POUR LE COUPLAGE HAUT-BAS DES VECTEURS PROPRES
    seuilc = 1.d-4
    call wkvect('&&WP4VEC.INDIC.PART.VP', 'V V I', nbvect, iadind)
    do j = 1, nbvect
        zi(iadind+j-1) = -2
    end do
    do j = 1, nbvect
        auxrj = dble(vp(j))
        auxij = dimag(vp(j))
        if (zi(iadind+j-1) .eq. -2) then
            if (abs(auxij) .lt. seuilr) then
                zi(iadind+j-1) = -3
                nbreel = nbreel+1
            else
                if (abs(auxrj) .lt. seuilr) auxrj = 0.d0
                k = j+1
                trouve = .false.
3               continue
                if ((.not. trouve) .and. (k .le. nbvect)) then
                    auxrk = dble(vp(k))
                    auxik = dimag(vp(k))
                    if (abs(auxrk) .lt. seuilr) auxrk = 0.d0
                    if (abs(auxik) .lt. seuilr) auxik = 0.d0
                    c1 = sqrt((auxrj-auxrk)**2+(auxij+auxik)**2)
                    if (c1 .lt. seuilp) then
                        lconj = .true.
                    else
                        lconj = .false.
                    end if
!
                    if ((zi(iadind+k-1) .eq. -2) .and. lconj .and. (auxij*auxik .le. 0.d0)) then
                        trouve = .true.
                        nbcmpc = nbcmpc+1
                        if (auxij .gt. 0.d0) then
                            zi(iadind+j-1) = 1
                            zi(iadind+k-1) = -1
                        else
                            zi(iadind+j-1) = -1
                            zi(iadind+k-1) = 1
                        end if
                    else
                        k = k+1
                    end if
                    goto 3
                end if
                if (.not. trouve) then
                    nbcmpp = nbcmpp+1
                    zi(iadind+j-1) = 0
                end if
            end if
        end if
    end do
!
    if (zi(iadind+nbvect-1) .eq. -2) then
        zi(iadind+nbvect-1) = 0
        nbcmpp = nbcmpp+1
    end if
!
    if (nbcmpp .gt. 0) then
        vali(1) = nbreel
        vali(2) = nbcmpc
        vali(3) = nbcmpp
        call utmess('A', 'ALGELINE4_87', ni=3, vali=vali)
    end if
!
    if (nbreel .gt. 0) then
        vali(1) = nbreel
        vali(2) = nbcmpc
        vali(3) = nbcmpp
        call utmess('I', 'ALGELINE4_88', ni=3, vali=vali)
    end if
!
!
! --- 1.2. DETERMINATION DE NB FREQUENCES GARDEES
    nbfrga = nbcmpc
!
! --- 1.3. ELIMINATION DES CONJUGUES (OPERATEUR REEL) -- COMPACTAGE --
    k = 1
    do j = 1, nbvect
        if (zi(iadind+j-1) .gt. 0) then
            if (k .ne. j) then
                vp(k) = vp(j)
                zi(iadind+k-1) = zi(iadind+j-1)
                do i = 1, neq, 1
                    vecp(i, k) = vecp(i, j)
                    vauc(i, k) = vauc(i, j)
                    vauc(i+neq, k) = vauc(i+neq, j)
                end do
            end if
            k = k+1
        end if
    end do
    nbfrga = k-1
! NBRE DE VP RECOMPACTEES
    nbfr = k-1
!
!
!     ---------- FIN DE PARTITION TEST ET ELIMINATION -----------------
!     ----------    AU NIVEAU DE L' OPERATEUR REEL    -----------------
!
! --- 2. CALCUL DES SOLUTIONS PROPRES DU PB QUADRATIQUE ---
    call wkvect('&&WP4VEC.VEC.AUX.C1', 'V V C', neq, av1)
    call wkvect('&&WP4VEC.VEC.AUX.C2', 'V V C', neq, av2)
    do j = 1, nbfr
        if (zi(iadind+j-1) .gt. 0) then
            call wptest(lagr, vauc(1, j), vauc(neq+1, j), vp(j), neq, &
                        nmabp)
            if (nmabp .gt. seuilc) then
                zi(iadind+j-1) = 0
                nbfrga = nbfrga-1
            end if
        end if
    end do
!
!
! --- 1.3. ELIMINATION DES VALEURS FAUSSES -- RECOMPACTAGE --
    k = 1
    do j = 1, nbfr
        if (zi(iadind+j-1) .gt. 0) then
            if (k .ne. j) then
                vp(k) = vp(j)
                zi(iadind+k-1) = zi(iadind+j-1)
!
                do i = 1, neq, 1
                    vecp(i, k) = vecp(i, j)
                    vauc(i, k) = vauc(i, j)
                    vauc(i+neq, k) = vauc(i+neq, j)
                end do
            end if
            k = k+1
        end if
    end do
    nbfrga = k-1
!
! --- 3. SELECTION DES VALEURS PROPRES (PB QUADRATIQUE)
    do j = 1, nbfrga, 1
        if ((zi(iadind+j-1) .eq. 1) .and. (dimag(vp(j)) .lt. 0.d0)) then
            vp(j) = dconjg(vp(j))
            do i = 1, neq
                vecp(i, j) = dconjg(vecp(i, j))
                vauc(i, j) = dconjg(vauc(i, j))
                vauc(i+neq, j) = dconjg(vauc(i+neq, j))
            end do
        end if
    end do
!
! --- 4. PREPARATION DE RESUFR
!
    if (nbfreq .gt. nbfrga) then
        vali(1) = nbfreq
        vali(2) = nbfrga
        nbfreq = nbfrga
        if (nbfreq .eq. 0) then
            kmsg = 'F'
        else
            kmsg = 'A'
        end if
        call getvtx(' ', 'TYPE_RESU', scal=typres, nbret=ibid)
        valk = 'FREQ'
        if (typres .ne. 'DYNAMIQUE') valk = 'CHAR_CRIT'
        call utmess(kmsg//'+', 'ALGELINE5_79', sk=valk, ni=2, vali=vali)
        if (kmsg .eq. 'A') then
            call utmess(kmsg//'+', 'ALGELINE5_80', sk=valk)
        end if
        call utmess(kmsg, 'ALGELINE5_81', sk=valk)
    end if
!
! --- 5. TRI (DANS LE SPECTRE ET DE PRESENTATION) DES VALEURS PROPRES-
!          - TRI DANS LE SPECTRE : SUIVANT ABS(SHIFT - VPQ)
    call wpordc(1, shift, vp, vecp, nbfrga, &
                neq)
!          - TRI DE PRESNTATION  : SUIVANT IM(VPQ) - IM(SHIFT)
    call wpordc(0, shift, vp, vecp, nbfreq, &
                neq)
!
! --- 5. PREPARATION DE RESUFR
    do j = 1, nbfreq
        am = dble(vp(j))*dble(vp(j))
        om = dimag(vp(j))*dimag(vp(j))
        resufi(j, 1) = j
        resufr(j, 2) = om
        resufr(j, 3) = -dble(vp(j))/sqrt(om+am)
    end do
!
! --- 6. DESTRUCTION DES OJB TEMPORAIRES
    call jedetr('&&WP4VEC.VEC.AUX.C1')
    call jedetr('&&WP4VEC.VEC.AUX.C2')
    call jedetr('&&WP4VEC.INDIC.PART.VP')
!
    call jedema()
end subroutine
