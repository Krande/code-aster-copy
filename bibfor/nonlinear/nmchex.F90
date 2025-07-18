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

subroutine nmchex(vachap, tychap, tyvari, nomvar)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterfort/nmchai.h"
    character(len=19) :: vachap(*)
    character(len=6) :: tychap
    character(len=6) :: tyvari
    character(len=19) :: nomvar
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
!
! RECUPERATION NOM DE LA VARIABLE DANS UNE VARIABLE CHAPEAU
!
! ----------------------------------------------------------------------
!
!
! IN  VACHAP : VARIABLE CHAPEAU
! IN  TYCHAP : TYPE DE VARIABLE CHAPEAU
!                MEELEM - NOMS DES MATR_ELEM
!                MEASSE - NOMS DES MATR_ASSE
!                VEELEM - NOMS DES VECT_ELEM
!                VEASSE - NOMS DES VECT_ASSE
!                SOLALG - NOMS DES CHAM_NO SOLUTIONS
!                VALINC - VALEURS SOLUTION INCREMENTALE
! IN  TYVARI : TYPE DE LA VARIABLE
! OUT NOMVAR : NOM DE LA VARIABLE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: index
!
! ----------------------------------------------------------------------
!
    call nmchai(tychap, tyvari, index)
    nomvar = vachap(index)
!
end subroutine
