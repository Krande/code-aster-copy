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
subroutine wp2vec(appr, opt, nbfreq, nbvect, neq, &
                  shift, yh, yb, vr, nlivr, &
                  vpr, vpi, vecp, mxresf, resufi, &
                  resufr, lagr, omecor)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/wpordo.h"
#include "asterfort/wprest.h"
#include "asterfort/wptest.h"
    character(len=1) :: appr
    character(len=*) :: opt
    integer(kind=8) :: nbfreq, nbvect, neq, lagr(*), mxresf, nlivr, resufi(mxresf, *)
    complex(kind=8) :: vecp(neq, *), shift
    real(kind=8) :: resufr(mxresf, *), yh(neq, *), yb(neq, *), vpr(*), vpi(*)
    real(kind=8) :: vr(nlivr, *), omecor
!     RESTITUTION DES VALEURS PROPRES ET DES MODES DU PB QUADRATIQUE
!     AVEC TRI SUIVANT LES PARTIES IMMAGINAIRES CROISANTES
!     ET NORMALISATION A LA PLUS GRANDE CMP NON LAGRANGE
!
! --> EN PUS ON MET UN TEST DE CONJUGAISON
! --> STOCKAGE D'UN SEUL COUPLE PARMI 2 CONJUGUES
!     -----------------------------------------------------------------
! IN  APPR   : K : INDICATEUR D' APPROCHE 'R' OU 'I'
! IN  OPT    : K : OPTION : 'CENTRE' OU 'PLUS_PETITE'
! IN  NBFREQ : I : NOMBRE DE MODES DEMANDES
! IN  NBVECT : I : NOMBRE DE VECTEURS DE LANCZOS
! IN  NEQ    : I : TAILLE DES MATRICES DU PB QUADRATIQUE
! IN  SHIFT  : R : VALEUR DU DECALAGE
! IN  YH     : R : PARTIE HAUTE DES VECTEURS DE LANCZOS
! IN  YB     : R : PARTIE BASSE DES VECTEURS DE LANCZOS
! IN  LAGR   : I : INDICATEUR DES NON-LAGRANGE
! VAR VPR    : R : IN  : PARTIE REELLE DES VALEURS PROPRE DU PB REDUIT
!            :   : OUT : PARTIE REELLE DES VALEURS PROPRE DU PB QUAD
! VAR VPI    : R : IN  : PARTIE IMMAGI DES VALEURS PROPRE DU PB REDUIT
!            :   : OUT : PARTIE IMMAGI DES VALEURS PROPRE DU PB QUAD
! IN  VR     : R : MODES DU PB REDUIT. C'EST L'EQUIVALENCE D'UNE
!                  MATRICE COMPLEXE DE LONGUEUR NBVECT.
! OUT VECP   : C : MODES DU PB QUADRATIQUE
! OUT RESUFR : C : TABLEAU DE POST-TRAITEMENT
! IN  OMECOR : R : "ZERO MODAL", SEUIL EN DECA DUQUEL DEUX MODES SONT
!                  CONSIDERES COMME IDENTIQUES
!     -----------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
    real(kind=8) :: si, mod2, a, b, nmabp, nmabm, am, om, eps, seuilr, seuilp
    real(kind=8) :: c1, auxrj, seuilc, auxij, auxrk, auxik
    integer(kind=8) :: i, j, k, av1, av2, av, iadind, nbfrga, vali(5), nbcmpp, nbcmpc
    integer(kind=8) :: nbreel, nbfr, ibid
    complex(kind=8) :: des, vpq, mhu, vpp, vpm
    aster_logical :: trouve, lconj
    character(len=1) :: kmsg
    character(len=16) :: valk, typres
!
!     -----------------------------------------------------------------
    call jemarq()
    si = dimag(shift)
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
!
!     SI IM(VP)<SEUILR: VP EST CONSIDEREE COMME REELLE
    seuilr = 1.d-7
!     SI MODULE(VPK-VPJ) < SEUILP: VPK = CONJUGEE DE VPJ
    seuilp = omecor
!     SEUIL POUR LE COUPLAGE HAUT-BAS DES VECTEURS PROPRES
    seuilc = 1.d-4
!
    call wkvect('&&WP2VEC.INDIC.PART.VP', 'V V I', nbvect, iadind)
    do j = 1, nbvect
        zi(iadind+j-1) = -2
    end do
    do j = 1, nbvect
        auxrj = vpr(j)
        auxij = vpi(j)
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
                    auxrk = vpr(k)
                    auxik = vpi(k)
                    if (abs(auxrk) .lt. seuilr) auxrk = 0.d0
                    if (abs(auxik) .lt. seuilr) auxik = 0.d0
                    c1 = sqrt((auxrj-auxrk)**2+(auxij+auxik)**2)
                    if (c1 .lt. seuilp) then
                        lconj = .true.
                    else
                        lconj = .false.
                    end if
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
! --- 1.2. DETERMINATION DE NB FREQUENCES GARDEES
    nbfrga = nbcmpc
!
! --- 1.3. ELIMINATION DES CONJUGUES (OPERATEUR REEL) -- COMPACTAGE --
    k = 1
    do j = 1, nbvect
        if (zi(iadind+j-1) .gt. 0) then
            if (k .ne. j) then
                vpr(k) = vpr(j)
                vpi(k) = vpi(j)
                zi(iadind+k-1) = zi(iadind+j-1)
                do i = 1, nlivr, 1
                    vr(i, k) = vr(i, j)
                end do
            end if
            k = k+1
        end if
    end do
    nbfrga = k-1
! NBRE DE VP RECOMPACTEES
    nbfr = k-1
!
!     ---------- FIN DE PARTITION TEST ET ELIMINATION -----------------
!     ----------    AU NIVEAU DE L' OPERATEUR REEL    -----------------
!
! --- 2. CALCUL DES SOLUTIONS PROPRES DU PB QUADRATIQUE ---
    if (opt .eq. 'CENTRE') then
        call wkvect('&&WP2VEC.VEC.AUX.C1', 'V V C', neq, av1)
        call wkvect('&&WP2VEC.VEC.AUX.C2', 'V V C', neq, av2)
        call wkvect('&&WP2VEC.VEC.AUX.C ', 'V V C', neq, av)
    end if
    do j = 1, nbfr
        if (zi(iadind+j-1) .gt. 0) then
            a = vpr(j)
            b = vpi(j)
            mhu = dcmplx(a, b)
            mod2 = a*a+b*b
            mod2 = 1.d0/mod2
            call wprest(yh, vr(1, j), neq, nbvect, vecp(1, j))
            if (opt .eq. 'PLUS_PETITE') then
                a = a*mod2
                b = -b*mod2
            else if (opt .eq. 'CENTRE') then
                call wprest(yb, vr(1, j), neq, nbvect, zc(av))
                if (appr .eq. 'R') then
                    des = dcmplx(1.d0, 0.d0)-dcmplx(4.d0*si*si, 0.d0)*mhu*mhu
                    des = sqrt(des)
                    vpq = .5d0*(dcmplx(1.d0, 0.d0)-dcmplx(0.d0, 2.d0*si)*mhu+des)/mhu
                    vpp = vpq+shift
                    call wptest(lagr, vecp(1, j), zc(av), vpp, neq, &
                                nmabp)
                    vpq = .5d0*(dcmplx(1.d0, 0.d0)-dcmplx(0.d0, 2.d0*si)*mhu-des)/mhu
                    vpm = vpq+shift
                    call wptest(lagr, vecp(1, j), zc(av), vpm, neq, &
                                nmabm)
                else
                    des = -dcmplx(si*si, 0.d0)*mhu*mhu+dcmplx(si, 0.d0)*mhu
                    des = sqrt(des)
                    vpq = -dcmplx(0.d0, si)+des/mhu
                    vpp = vpq+shift
                    call wptest(lagr, vecp(1, j), zc(av), vpp, neq, &
                                nmabp)
                    vpq = -dcmplx(0.d0, si)-des/mhu
                    vpm = vpq+shift
                    call wptest(lagr, vecp(1, j), zc(av), vpm, neq, &
                                nmabm)
                end if
                if (nmabm .lt. nmabp) then
                    a = dble(vpm)
                    b = dimag(vpm)
                    eps = nmabm
                else
                    a = dble(vpp)
                    b = dimag(vpp)
                    eps = nmabp
                end if
                if (eps .gt. seuilc) then
                    zi(iadind+j-1) = 0
                    nbfrga = nbfrga-1
                end if
            end if
            vpr(j) = a
            vpi(j) = b
        end if
    end do
!
! --- 1.3. ELIMINATION DES VALEURS FAUSSES -- RECOMPACTAGE --
    k = 1
    do j = 1, nbfr
        if (zi(iadind+j-1) .gt. 0) then
            if (k .ne. j) then
                vpr(k) = vpr(j)
                vpi(k) = vpi(j)
                zi(iadind+k-1) = zi(iadind+j-1)
                do i = 1, nlivr, 1
                    vr(i, k) = vr(i, j)
                end do
            end if
            k = k+1
        end if
    end do
    nbfrga = k-1
!
! --- 3. SELECTION DES VALEURS PROPRES (PB QUADRATIQUE)
    do j = 1, nbfrga, 1
        if ((zi(iadind+j-1) .eq. 1) .and. (vpi(j) .lt. 0.d0)) then
            vpi(j) = -vpi(j)
            do i = 1, neq
                vecp(i, j) = dconjg(vecp(i, j))
            end do
        end if
    end do
!
! --- 4. PREPARATION DE RESUFR
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
    call wpordo(1, shift, vpr, vpi, vecp, &
                nbfrga, neq)
    call wpordo(0, shift, vpr, vpi, vecp, &
                nbfreq, neq)
!
    do j = 1, nbfreq
        am = vpr(j)*vpr(j)
        om = vpi(j)*vpi(j)
        resufi(j, 1) = j
        resufr(j, 2) = om
        resufr(j, 3) = -vpr(j)/sqrt(om+am)
    end do
!
! --- 6. DESTRUCTION DES OJB TEMPORAIRES
    if (opt .eq. 'CENTRE') then
        call jedetr('&&WP2VEC.VEC.AUX.C1')
        call jedetr('&&WP2VEC.VEC.AUX.C2')
        call jedetr('&&WP2VEC.VEC.AUX.C ')
    end if
    call jedetr('&&WP2VEC.INDIC.PART.VP')
!
    call jedema()
end subroutine
