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
subroutine folocx(vale, n, x, prolgd, i, &
                  epsi, coli, ier)
    implicit none
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: n, i, ier
    real(kind=8) :: vale(n)
    real(kind=8) :: valr(2)
    character(len=*) :: prolgd
    character(len=1) :: coli
! person_in_charge: mathieu.courtois at edf.fr
!     RECHERCHE DE LA PLACE DE X DANS LE VECTEUR VALE ORDONNE CROISSANT
!     ON VERIFIE SI X EST DANS L'INTERVALLE (V(1),V(N))
!                SINON, SUIVANT PROLGD, ON AGIT...
!     ------------------------------------------------------------------
! IN  : VALE   : VECTEUR DES VALEURS DES ABSCISSES
! IN  : N      : NOMBRE DE POINTS DE VALE
! IN  : X      : VALEUR DE L'ABSCISSE
! IN  : PROLGD : PROLONGEMENTS A DROITE ET A GAUCHE
! VAR : I      : NUMERO DU POINT TEL QUE VALE(I) <= X
! IN  : EPSI   : PRECISION A LAQUELLE ON RECHERCHE LA VALEUR
! OUT : COLI   : = 'C',   X = VALE(I)
!                = 'I',   INTERPOLATION VALE(I) < X < VALE(I+1)
!                = 'E',   EXTRAPOLATION PERMISE
!                = 'T',   FONCTION INTERPRETEE
!                = '?',   IER > 0
! OUT : IER    : CODE RETOUR
!                IER = 10 --->  MOINS DE 1 POINT
!                IER = 20 --->  EXTRAPOLATION INCONNUE
!                IER = 30 --->  ON DEBORDE A GAUCHE
!                IER = 40 --->  ON DEBORDE A DROITE
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: id, ie, ind, j
!
    real(kind=8) :: epsi, tole, x
!-----------------------------------------------------------------------
    ier = 0
    coli = '?'
    if (n .lt. 1) then
        ier = 10
        call utmess('E', 'FONCT0_18')
        goto 999
    else if (n .eq. 1) then
!             ON A X = VALE(1) + EPSILON
!             ON A : X < VALE(1) ET PROL GAUCHE AUTORISE
!             ON A : X > VALE(1) ET PROL DROITE AUTORISE
        if ((x .eq. 0.d0 .and. abs(vale(n)) .le. epsi) .or. &
            (x .lt. vale(n) .and. prolgd(1:1) .ne. 'E') .or. &
            (x .gt. vale(n) .and. prolgd(2:2) .ne. 'E')) then
            i = n
            coli = 'C'
        else if (abs((vale(n)-x)/x) .le. epsi) then
            i = n
            coli = 'C'
        else
            ier = 30
            call utmess('E', 'FONCT0_23')
        end if
        goto 999
    end if
!
!     --- PROLONGEMENT A GAUCHE ---
    if (x .le. vale(1)) then
        i = 1
        if (prolgd(1:1) .eq. 'E') then
            tole = epsi*abs(vale(1)-vale(2))
            if (abs(vale(1)-x) .le. tole) then
                coli = 'C'
                goto 999
            end if
            ier = 30
            valr(1) = x
            valr(2) = vale(1)
            call utmess('E', 'FONCT0_19', nr=2, valr=valr)
            goto 999
        else if (prolgd(1:1) .eq. 'L') then
            coli = 'E'
        else if (prolgd(1:1) .eq. 'C') then
            coli = 'C'
        else if (prolgd(1:1) .eq. 'I') then
            coli = 'T'
        else
            ier = 20
            call utmess('E', 'FONCT0_21', sk=prolgd(1:1))
            goto 999
        end if
!
!     --- PROLONGEMENT A DROITE ---
    else if (x .ge. vale(n)) then
        i = n
        if (prolgd(2:2) .eq. 'E') then
            tole = epsi*abs(vale(n)-vale(n-1))
            if (abs(vale(n)-x) .le. tole) then
                coli = 'C'
                goto 999
            end if
            ier = 40
            valr(1) = x
            valr(2) = vale(n)
            call utmess('E', 'FONCT0_20', nr=2, valr=valr)
            goto 999
        else if (prolgd(2:2) .eq. 'C') then
            coli = 'C'
        else if (prolgd(2:2) .eq. 'I') then
            coli = 'T'
        else if (prolgd(2:2) .eq. 'L') then
            i = n-1
            coli = 'E'
        else
            ier = 20
            call utmess('E', 'FONCT0_21', sk=prolgd(2:2))
            goto 999
        end if
!
!     --- RECHERCHE DE LA VALEUR PAR DICHOTOMIE ---
    else
        if (i .lt. 1 .or. i .gt. n) i = n/2
        if (vale(i) .le. x) then
            id = i
            ie = n
        else
            id = 1
            ie = i
        end if
        do j = 1, n
            if (ie .eq. (id+1)) goto 3
            ind = id+(ie-id)/2
            if (x .ge. vale(ind)) then
                id = ind
            else
                ie = ind
            end if
        end do
3       continue
        i = id
        coli = 'I'
    end if
    ASSERT(i .ge. 1 .and. i .le. n)
!
999 continue
!
end subroutine
