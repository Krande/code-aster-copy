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
subroutine mdfext(tinit, dt, neqgen, nbexci, idescf, &
                  nomfon, coefm, liad, inumor, nbpas, &
                  f)
    implicit none
#include "jeveux.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: neqgen, nbexci, nbpas
    integer(kind=8) :: idescf(*), liad(*), inumor(*)
    real(kind=8) :: tinit, dt, t, coefm(*), f(neqgen, *)
    character(len=8) :: nomfon(*)
!
!
!     CALCULE LES FORCES EXTERIEURES A CHAQUE PAS DE TEMPS
!     ------------------------------------------------------------------
! IN  : TINIT  : PAS DE TEMPS INITIAL POUR LE BLOC DE CALCUL
! IN  : NEQGEN : NOMBRE D'EQUATIONS GENERALISEES
! IN  : NBEXCI : NOMBRE D'EXCITATION (MC EXCIT ET EXCIT_RESU)
! IN  : IDESCF : TYPE D'EXCITATION (VECT_ASSE/NUME_ORDRE,FONC_MULT/
!                COEF_MULT)
! IN  : NOMFON : NOM DES FONC_MULT (QUAND IL Y EN A)
! IN  : COEFM  : VALEUR DU COEF_MULT
! IN  : LIAD   : VALEUR DU VECT_ASSE
! IN  : NUMOR  : NUME_ORDRE DU MODE EXCITE
! IN  : NBPAS  : NOMBRE DE PAS DE CALCUL
! OUT : F      : TABLEAU DES FORCES EXTERIEURES A CHAQUE PAS DE TEMPS
! ----------------------------------------------------------------------
!
!
!
!
    integer(kind=8) :: ier
    character(len=4) :: nompar
!
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    real(kind=8) :: alpha
!-----------------------------------------------------------------------
    call jemarq()
!
    nompar = 'INST'
    do i = 1, nbexci
        t = tinit
        do k = 1, nbpas
            if (idescf(i) .eq. 1) then
                call fointe('F ', nomfon(i), 1, [nompar], [t], &
                            alpha, ier)
                do j = 1, neqgen
                    f(j, k) = f(j, k)+alpha*zr(liad(i)+j-1)
                end do
            else if (idescf(i) .eq. 2) then
                call fointe('F ', nomfon(i), 1, [nompar], [t], &
                            alpha, ier)
                f(inumor(i), k) = f(inumor(i), k)+alpha
            else if (idescf(i) .eq. 3) then
                do j = 1, neqgen
                    f(j, k) = f(j, k)+coefm(i)*zr(liad(i)+j-1)
                end do
            else if (idescf(i) .eq. 4) then
                f(inumor(i), k) = f(inumor(i), k)+coefm(i)
            end if
            t = tinit+(k*dt)
        end do
    end do
!
!
    call jedema()
end subroutine
