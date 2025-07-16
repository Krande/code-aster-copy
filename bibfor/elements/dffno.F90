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
subroutine dffno(elrefe, ndim, nno, nnos, dff)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfdf.h"
    character(len=*) :: elrefe
    integer(kind=8) :: ndim, nno, nnos
    real(kind=8) :: dff(*)
! BUT:   CALCUL DES DERIVEES DES FONCTIONS DE FORMES
!        AUX NOEUDS D'UN ELREFE
!

    real(kind=8) :: x(3, MT_NNOMAX), tab(3, MT_NNOMAX)
    integer(kind=8) :: ino, ideri, ifonc
! ----------------------------------------------------------------------

    call elrfno(elrefe, nno, nnos, ndim, x)
    do ino = 1, nno
        call elrfdf(elrefe, x(1:ndim, ino), tab)
        do ideri = 1, ndim
            do ifonc = 1, nno
                dff((ino-1)*nno*ndim+(ideri-1)*nno+ifonc) = tab(ideri, ifonc)
            end do
        end do
    end do
!
end subroutine
