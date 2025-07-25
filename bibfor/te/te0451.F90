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
subroutine te0451(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/excent.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvala.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
    character(len=16) :: option, nomte
! person_in_charge: jacques.pellet at edf.fr
! ======================================================================
!  BUT:  CALCUL DE L'OPTION EFGE_ELGA
!        POUR LES ELEMENTS DE COQUE A "SOUS-POINTS"
!        ON PART DE SIEF_ELGA ET ON INTEGRE DANS L'EPAISSEUR
! ......................................................................
!
    integer(kind=8) :: j1, nbcou, npgh, jsigm, idec, jeff, npg, itab(7), iret
    integer(kind=8) :: nbsp, kpg, ibid, nbsig, nbeff, icou, jmate, icodre(1), j2
    real(kind=8) :: nxx, nyy, mxx, myy, nxy, mxy, qx, qy, excen
    real(kind=8) :: r8bid, cb, cm, ch, h, hb, hm, hh
    real(kind=8) :: siyyb, siyym, siyyh, sixxb, sixxm, sixxh, sixyb, sixym
    real(kind=8) :: sixyh
    real(kind=8) :: siyzb, siyzm, siyzh, sixzb, sixzm, sixzh, epcou(100), epi(1)
    character(len=8) :: alias8
    character(len=16) :: nomres
    character(len=3) :: cmod, num
    character(len=2) :: val
    aster_logical :: lcoqmu, lreel
!     ------------------------------------------------------------------
    call teattr('S', 'ALIAS8', alias8, ibid)
    cmod = alias8(3:5)
    if (cmod .eq. 'DKT' .or. cmod .eq. 'DST' .or. cmod .eq. 'Q4G' .or. cmod .eq. 'CQ3') then
        nbsig = 6
        nbeff = 8
    else if (cmod .eq. 'CQA' .or. cmod .eq. 'CQC' .or. cmod .eq. 'CQD') then
        nbsig = 4
        nbeff = 6
    else
        ASSERT(.false.)
    end if
!
!     -- EPAISSEUR :
    call jevech('PCACOQU', 'L', j1)
    h = zr(j1)
!
!     -- NOMBRE DE COUCHES :
    call jevech('PNBSP_I', 'L', j2)
    nbcou = zi(j2)
!
!
!     -- SI LE MATERIAU EST 'ELAS_COQMU', LES COUCHES
!        N'ONT PAS LA MEME EPAISSEUR.
!        ON LES STOCKE DANS EPCOU
!     ------------------------------------------------
    lcoqmu = .false.
    call jevech('PMATERC', 'L', jmate)
    call codent(1, 'G', num)
    call codent(1, 'G', val)
    nomres = 'C'//num//'_V'//val
    r8bid = 0.d0
    call rcvala(zi(jmate), ' ', 'ELAS_COQMU', 0, ' ', &
                [r8bid], 1, nomres, epi(1), icodre(1), &
                0)
    if (icodre(1) .eq. 0) lcoqmu = .true.
    if (lcoqmu) then
        ASSERT(nbcou .le. 100)
        do icou = 1, nbcou
            call codent(icou, 'G', num)
            nomres = 'C'//num//'_V'//val
            call rcvala(zi(jmate), ' ', 'ELAS_COQMU', 0, ' ', &
                        [r8bid], 1, nomres, epi(1), icodre(1), &
                        0)
            ASSERT(icodre(1) .eq. 0)
            ASSERT(epi(1) .ge. 0.d0)
            epcou(icou) = epi(1)
        end do
    end if
!
!
!     -- CONTRAINTES DANS LES COUCHES :
!     ----------------------------------
    call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, &
                itab=itab)
    jsigm = itab(1)
    npg = itab(3)
    nbsp = itab(7)
    npgh = 3
    ASSERT(nbsp .eq. nbcou*npgh)
    ASSERT(itab(2) .eq. nbsig*npg)
!
!
!     -- CALCUL DES EFFORTS PAR INTEGRATION DANS L'EPAISSEUR :
!     --------------------------------------------------------
    call tecach('OOO', 'PEFGER', 'E', iret, nval=7, &
                itab=itab)
    jeff = itab(1)
    ASSERT(itab(2) .eq. nbeff*npg)
!
!     -- BOUCLE SUR LES POINTS DE GAUSS :
    do kpg = 1, npg
        nxx = 0.d0
        nyy = 0.d0
        nxy = 0.d0
        mxx = 0.d0
        myy = 0.d0
        mxy = 0.d0
        qx = 0.d0
        qy = 0.d0
!
!       -- BOUCLE SUR LES COUCHES :
        hb = -h/2
        do icou = 1, nbcou
            idec = ((kpg-1)*nbcou+(icou-1))*npgh*nbsig
!
!         -- HB, HM, HH : "HAUTEUR" DES SOUS-POINTS :
            if (lcoqmu) then
                epi(1) = epcou(icou)
            else
                epi(1) = h/nbcou
            end if
            hm = hb+epi(1)/2.d0
            hh = hm+epi(1)/2.d0
!
!         -- SIXXB, SIYYB, ... : CONTRAINTES AU BAS DE LA COUCHE
            sixxb = zr(jsigm-1+idec+1)
            siyyb = zr(jsigm-1+idec+2)
            sixyb = zr(jsigm-1+idec+4)
            if (nbsig .eq. 6) then
                sixzb = zr(jsigm-1+idec+5)
                siyzb = zr(jsigm-1+idec+6)
            end if
!         -- SIXXM, SIYYM, ... : CONTRAINTES AU MILIEU DE LA COUCHE
            sixxm = zr(jsigm-1+idec+1+nbsig)
            siyym = zr(jsigm-1+idec+2+nbsig)
            sixym = zr(jsigm-1+idec+4+nbsig)
            if (nbsig .eq. 6) then
                sixzm = zr(jsigm-1+idec+5+nbsig)
                siyzm = zr(jsigm-1+idec+6+nbsig)
            end if
!
!         -- SIXXH, SIYYH, ... : CONTRAINTES EN HAUT DE LA COUCHE
            sixxh = zr(jsigm-1+idec+1+2*nbsig)
            siyyh = zr(jsigm-1+idec+2+2*nbsig)
            sixyh = zr(jsigm-1+idec+4+2*nbsig)
            if (nbsig .eq. 6) then
                sixzh = zr(jsigm-1+idec+5+2*nbsig)
                siyzh = zr(jsigm-1+idec+6+2*nbsig)
            end if
!
!         -- ON INTEGRE DANS L'EPAISSEUR DE CHAQUE COUCHE
!            AVEC UNE FORRMULE DE NEWTON-COTES A 3 POINTS
!            LES COEFFICIENTS SONT 1/6, 4/6 ET 1/6
            cb = epi(1)/6
            cm = 4.d0*epi(1)/6
            ch = epi(1)/6
!
!         -- NXX, NYY, NXY = SOMME DE SIXX, SIYY, SIXY :
            nxx = nxx+cb*sixxb+cm*sixxm+ch*sixxh
            nyy = nyy+cb*siyyb+cm*siyym+ch*siyyh
            nxy = nxy+cb*sixyb+cm*sixym+ch*sixyh
!
            if (nbeff .eq. 8) then
!           -- QX, QY = SOMME DE SIXZ, SIYZ
                qx = qx+cb*sixzb+cm*sixzm+ch*sixzh
                qy = qy+cb*siyzb+cm*siyzm+ch*siyzh
            end if
!
!         -- MXX, MYY, MXY = MOMENTS DE SIXX, SIYY, SIXY :
            mxx = mxx+cb*sixxb*hb+cm*sixxm*hm+ch*sixxh*hh
            myy = myy+cb*siyyb*hb+cm*siyym*hm+ch*siyyh*hh
            mxy = mxy+cb*sixyb*hb+cm*sixym*hm+ch*sixyh*hh
!
!         -- MISE A JOUR DE HB POUR LA COUCHE SUIVANTE :
            hb = hb+epi(1)
        end do
!
        zr(jeff-1+(kpg-1)*nbeff+1) = nxx
        zr(jeff-1+(kpg-1)*nbeff+2) = nyy
        zr(jeff-1+(kpg-1)*nbeff+4) = mxx
        zr(jeff-1+(kpg-1)*nbeff+5) = myy
        if (nbeff .eq. 8) then
            zr(jeff-1+(kpg-1)*nbeff+3) = nxy
            zr(jeff-1+(kpg-1)*nbeff+6) = mxy
            zr(jeff-1+(kpg-1)*nbeff+7) = qx
            zr(jeff-1+(kpg-1)*nbeff+8) = qy
        end if
    end do
!
!
!     -- POUR LES COQUES EXCENTREES, LES EFFORTS CALCULES SONT
!        DANS LE PLAN 'MOYEN'. IL FAUT LES CALCULER DANS LE PLAN 'MAIL'
!     -----------------------------------------------------------------
    if (cmod .eq. 'DKT' .or. cmod .eq. 'DST') then
        excen = zr(j1-1+5)
        lreel = .true.
        call excent('MAIL', excen, npg, nbeff, lreel, &
                    zr(jeff), zr(jeff), zc(jeff), zc(jeff))
    end if
!
end subroutine
