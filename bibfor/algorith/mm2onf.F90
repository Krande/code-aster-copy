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

subroutine mm2onf(ndim, nno, alias, ksi1, ksi2, &
                  ddff)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterfort/assert.h"
#include "asterfort/elrfd2.h"
    character(len=8) :: alias
    real(kind=8) :: ksi1, ksi2
    real(kind=8) :: ddff(3, 9)
    integer(kind=8) :: nno, ndim
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (TOUTES METHODES - UTILITAIRE)
!
! CALCUL DES DERIVEES SECONDES DES FONCTIONS DE FORME EN UN POINT
! DE L'ELEMENT DE REFERENCE
!
! ----------------------------------------------------------------------
!
! IN  ALIAS  : NOM D'ALIAS DE L'ELEMENT
! IN  NNO    : NOMBRE DE NOEUD DE L'ELEMENT
! IN  NDIM   : DIMENSION DE LA MAILLE (2 OU 3)
! IN  KSI1   : POINT DE CONTACT SUIVANT KSI1 DES
!               FONCTIONS DE FORME ET LEURS DERIVEES
! IN  KSI2   : POINT DE CONTACT SUIVANT KSI2 DES
!               FONCTIONS DE FORME ET LEURS DERIVEES
! OUT DDFF   : DERIVEES SECONDES DES FONCTIONS DE FORME EN XI YI
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ibid1, ibid2
    real(kind=8) :: ksi(2)
    real(kind=8) :: d2ff(3, 3, 9)
!
! ----------------------------------------------------------------------
!
!
! --- INITIALISATIONS
!
    ddff(:, :) = 0.d0
    d2ff(:, :, :) = 0.d0
!
    ksi(1) = ksi1
    ksi(2) = ksi2

    if ((nno .lt. 1) .or. (nno .gt. 9) .or. (ndim .lt. 1) .or. (ndim .gt. 3)) then
        ASSERT(.false.)
    end if
! --- RECUP DERIVEES SECONDES DES FONCTIONS DE FORME
!
    call elrfd2(alias, ksi, nno*ndim*ndim, d2ff, ibid1, &
                ibid2)
!
! --- CONVERSION XI-YI/YI-XI -> KSI1-KSI2
!

    ddff(1, :) = d2ff(1, 1, :)
    ddff(2, :) = d2ff(2, 2, :)
    ddff(3, :) = d2ff(1, 2, :)

!
end subroutine
