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
subroutine xmoffc(lact, nlact, nno, ffe, ffc)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "asterfort/assert.h"
    integer(kind=8) :: lact(*), nlact, nno
    real(kind=8) :: ffe(*), ffc(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT AVEC XFEM
! TRAVAIL EFFECTUE EN COLLABORATION AVEC I.F.P.
!
! POUR LA FORMULATION AUX NOEUDS SOMMET,  SI UN NOEUD N'EST PAS ACTIF,
! ON REPARTI SES FF EQUITABLEMENT SUR LES NOEUDS ACTIF
!
! ----------------------------------------------------------------------
!
! IN  NNO    : NOMBRE DE NOEUD DE L'ELEMENT PARENT
! IN LACT    : LITE DES LAGRANGES ACTIFS
! IN NLACT   : NOMBRE TOTAL DE LAGRANGES ACTIFS
! IN FFE     : FONCTION DE FORMES DE L'ELEMENT ESCLAVE PARENT
! OUT FFC    : FONCTION DE FORMES DU POINT DE CONTACT
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j
!
! ----------------------------------------------------------------------
!
!
    do i = 1, nno
        ffc(i) = ffe(i)
    end do
    if (nlact .lt. nno) then
        do i = 1, nno
            if (lact(i) .eq. 0) then
                do j = 1, nno
                    if (i .ne. j .and. lact(j) .ne. 0) then
                        ffc(j) = ffc(j)+ffc(i)/nlact
                    end if
                end do
                ffc(i) = 0.d0
            end if
        end do
    end if
!
    if (nno .gt. 8) then
        ASSERT(.false.)
    end if
!
end subroutine
