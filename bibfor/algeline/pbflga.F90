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
subroutine pbflga(umoy, hmoy, rmoy, long, cf0, &
                  fsvr, icoq, imod, nbm, tcoef, &
                  s1, s2, lambda, kcalcu, condit, &
                  gamma)
    implicit none
! COUPLAGE FLUIDELASTIQUE, CONFIGURATIONS DU TYPE "COQUE_COAX"
! RESOLUTION DU PROBLEME FLUIDE INSTATIONNAIRE : CALCUL DE GAMMA(3)
! DANS LE CAS OU UMOY <> 0
! APPELANT : PBFLUI
!-----------------------------------------------------------------------
!  IN : UMOY   : VITESSE DE L'ECOULEMENT MOYEN
!  IN : HMOY   : JEU ANNULAIRE MOYEN
!  IN : RMOY   : RAYON MOYEN
!  IN : LONG   : LONGUEUR DU DOMAINE DE RECOUVREMENT DES DEUX COQUES
!  IN : CF0    : COEFFICIENT DE FROTTEMENT VISQUEUX
!  IN : FSVR   : OBJET .FSVR DU CONCEPT TYPE_FLUI_STRU
!  IN : ICOQ   : INDICE CARACTERISANT LA COQUE SUR LAQUELLE ON TRAVAILLE
!                ICOQ=1 COQUE INTERNE  ICOQ=2 COQUE EXTERNE
!  IN : IMOD   : INDICE DU MODE CONSIDERE
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!  IN : TCOEF  : TABLEAU DES COEFFICIENTS DES DEFORMEES AXIALES
!  IN : S1     : PARTIE REELLE     DE LA FREQUENCE COMPLEXE
!  IN : S2     : PARTIE IMAGINAIRE DE LA FREQUENCE COMPLEXE
!  IN : LAMBDA : VALEURS PROPRES DE L'OPERATEUR DIFFERENTIEL
!  IN : KCALCU : MATRICE RECTANGULAIRE A COEFFICIENTS CONSTANTS
!                PERMETTANT DE CALCULER UNE SOLUTION PARTICULIERE DU
!                PROBLEME FLUIDE INSTATIONNAIRE, LORSQUE UMOY <> 0
! OUT : CONDIT : COEFFICIENTS DE PRECONDITIONNEMENT
! OUT : GAMMA  : COEFFICIENTS DE LA COMBINAISON LINEAIRE DONNANT LA
!                SOLUTION GENERALE DU PROBLEME FLUIDE INSTATIONNAIRE
!                (DECOMPOSITION SUR UNE FAMILLE D'EXPONENTIELLES)
!                LORSQUE UMOY <> 0
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/defaxe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/lcsolz.h"
#include "asterfort/pbflkz.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    real(kind=8) :: umoy, hmoy, rmoy, long, cf0, fsvr(7)
    integer(kind=8) :: icoq, imod, nbm
    real(kind=8) :: tcoef(10, nbm), s1, s2
    complex(kind=8) :: lambda(3), kcalcu(3, 4)
    real(kind=8) :: condit(3)
    complex(kind=8) :: gamma(3)
!
    real(kind=8) :: ln
    complex(kind=8) :: ei, r, t
!
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, imata, iret, itab, j
    real(kind=8) :: cde, cdep, cds, cdsp, reeli, rhof, u
    real(kind=8) :: v, w, x
!-----------------------------------------------------------------------
    call jemarq()
!
    rhof = fsvr(1)
    cde = fsvr(4)
    cdep = fsvr(5)
    cds = fsvr(6)
    cdsp = fsvr(7)
!
    u = 1.d0+cde
    v = 1.d0-cds
    r = -1.d0*dcmplx(s1/umoy+(cf0/hmoy), s2/umoy)
    w = -0.5d0*cdep*umoy*defaxe(icoq, imod, 0.d0, long, nbm, tcoef)
    x = 0.5d0*cdsp*umoy*defaxe(icoq, imod, long, long, nbm, tcoef)
!
    itab = 0
    if (icoq .eq. 2) itab = 5
    ln = tcoef(1+itab, imod)
    gamma(1) = pbflkz(2, 0.d0, long, ln, kcalcu)/rmoy
    gamma(2) = -1.d0*u*pbflkz(1, 0.d0, long, ln, kcalcu)- &
               pbflkz(3, 0.d0, long, ln, kcalcu)/ &
               (rhof*umoy)+dcmplx(w)
    gamma(3) = -1.d0*v*pbflkz(1, long, long, ln, kcalcu)- &
               pbflkz(3, long, long, ln, kcalcu)/ &
               (rhof*umoy)+dcmplx(x)
    do i = 1, 3
        reeli = dble(lambda(i))
        if (reeli .gt. 0.d0) then
            condit(i) = 1.d0
        else
            condit(i) = 0.d0
        end if
    end do
!
    call wkvect('&&PBFLGA.TEMP.MATA', 'V V C', 3*3, imata)
    do j = 1, 3
        ei = dcmplx(exp(-1.d0*condit(j)*lambda(j)*long))
        zc(imata+3*(j-1)) = lambda(j)*ei
        t = lambda(j)*(r-lambda(j))*rmoy*rmoy
        zc(imata+3*(j-1)+1) = (t+dcmplx(u))*ei
        ei = dcmplx(exp((1.d0-condit(j))*lambda(j)*long))
        zc(imata+3*(j-1)+2) = (t+dcmplx(v))*ei
    end do
!
    call lcsolz(zc(imata), gamma, 3, 3, 1, &
                iret)
    if (iret .ne. 0) then
        call utmess('F', 'ALGELINE3_17')
    end if
!
    call jedetr('&&PBFLGA.TEMP.MATA')
    call jedema()
end subroutine
