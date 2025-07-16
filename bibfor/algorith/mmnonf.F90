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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine mmnonf(ndim, nno, alias, ksi1, ksi2, ff)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrfvf.h"
!
    character(len=8) :: alias
    real(kind=8) :: ksi1, ksi2
    real(kind=8) :: ff(9)
    integer(kind=8) :: nno, ndim
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (TOUTES METHODES - UTILITAIRE)
!
! CALCUL DES FONCTIONS DE FORME EN UN POINT DE L'ELEMENT DE REFERENCE
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
! OUT FF     : FONCTIONS DE FORMES EN XI,YI
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: ksi(2)
!
! ----------------------------------------------------------------------
!
    ff(:) = 0.d0
    ksi(1) = ksi1
    ksi(2) = ksi2
    ASSERT(nno .ge. 1)
    ASSERT(nno .le. 9)
    ASSERT(ndim .ge. 1)
    ASSERT(ndim .le. 3)
!
! --- RECUP FONCTIONS DE FORME
!
    call elrfvf(alias, ksi, ff)
!
end subroutine
