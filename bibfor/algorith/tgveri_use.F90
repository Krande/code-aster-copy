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
subroutine tgveri_use(option, carcri, compor, iret)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    integer(kind=8), intent(out) :: iret
!
! ----------------------------------------------------------------------
! VAR OPTION NOM DE L'OPTION DE CALCUL
!             IN  : CELLE UTILISEE PAR LE TE
!             OUT : 'RAPH_MECA' SI BOUCLE, 'FULL_MECA' SI FIN DE BOUCLE
! IN  CARCRI  : CARCRI(1) = type de matrice tangente
!               0 : ANALYTIQUE, on ne passe pas ici
!               1 : PERTURBATION, on calcule Ktgte (FULL_MECA)
!               2 : VERIFICATION, on calcule Ktgte (FULL_MECA) + Kpertu
!               CARCRI(7) = valeur de la perturbation
! OUT IRET   SI IRET = 0 -> NON UTILISE, SINON -> UTILISE
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: typeMatr
! ----------------------------------------------------------------------
!
!     Calcul de la matrice TGTE par PERTURBATION
!
    iret = 1
    typeMatr = nint(carcri(TYPE_MATR_T))
    if (typeMatr .eq. 0 .or. typeMatr .eq. 3 .or. typeMatr .eq. 4) then
        iret = 0
    else
! INCOMATIBILITE AVEC LES COMPORTEMENTS QUI UTILISENT PVARIMP
        if (compor(PLANESTRESS) .eq. 'DEBORST') then
            iret = 0
        end if
    end if
    if (option(1:9) .eq. 'RIGI_MECA') then
        iret = 0
    end if
end subroutine
