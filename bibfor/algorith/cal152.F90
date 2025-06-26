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

subroutine cal152(option, max, may, maz, model, &
                  phib24, iphi1, iphi2, imade, modmec, &
                  chamno, num, vrai, i, j, &
                  mij, cij, kij)
    implicit none
! AUTEUR : G.ROUSSEAU
! OPERATEUR CALCULANT LA MASSE AJOUTEE, L'AMORTISSEMENT
!  ET LA RIGIDITE AJOUTEE EN THEORIE POTENTIELLE : CALC_MATR_AJOU
!     SUR BASE MODALE DE LA STRUCTURE DANS LE VIDE
!---------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calamr.h"
#include "asterfort/calmaj.h"
#include "asterfort/utmess.h"
    aster_logical :: vrai
    integer(kind=8) :: i, j
    integer(kind=8) :: imade
    integer(kind=8) :: iphi1, iphi2
    real(kind=8) :: mij, cij, kij, kij1, cij1, cij2
    real(kind=8) :: valr(2)
    character(len=2) :: model
    character(len=8) :: modmec
    character(len=9) :: option
    character(len=14) :: num
    character(len=19) :: max, may, maz, chamno
    character(len=24) :: phib24
! -----------------------------------------------------------------
    if (option .eq. 'MASS_AJOU') then
        call calmaj(option, max, may, maz, model, &
                    zk24(iphi1+j-1) (1:19), modmec, chamno, num, vrai, &
                    i, j, mij)
    end if
!
    if (option .eq. 'AMOR_AJOU') then
        call calmaj(option, max, may, maz, model, &
                    zk24(iphi2+j-1) (1:19), modmec, chamno, num, vrai, &
                    i, j, cij1)
        call calamr(phib24, zk24(iphi1+j-1) (1:19), zk24(imade+i-1), num, j, &
                    cij2)
        cij = cij1+cij2
        valr(1) = cij1
        valr(2) = cij2
        call utmess('I', 'ALGORITH14_80', nr=2, valr=valr)
    end if
!
    if (option .eq. 'RIGI_AJOU') then
        call calamr(phib24, zk24(iphi2+j-1) (1:19), zk24(imade+i-1), num, j, &
                    kij1)
        kij = kij1
    end if
end subroutine
