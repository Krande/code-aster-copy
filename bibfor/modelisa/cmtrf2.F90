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
subroutine cmtrf2(codcm1, codtrf, ncm1, lcm1, ntrf, &
                  ltrf, nbma, codint, lint, nint)
    implicit none
#include "asterfort/assert.h"
    integer(kind=8) :: codcm1, codtrf, codint, ncm1, ntrf, nint, nbma
    integer(kind=8) :: lint(nbma), lcm1(ncm1), ltrf(ntrf)
!  BUT :
!  -----
!  ETABLIR LA LISTE DES NUMEROS DE MAILLES (LINT) APPARTENANT
!  AUX 2 LISTES LCM1 ET LTRF
! ----------------------------------------------------------------------
!
    integer(kind=8) :: k
! ----------------------------------------------------------------------
    ASSERT(codcm1 .eq. 1 .or. codcm1 .eq. 3)
    ASSERT(codtrf .eq. 1 .or. codtrf .eq. 3)
    ASSERT(nbma .gt. 0)
    codint = 3
!
    if (codcm1 .eq. 1) then
        if (codtrf .eq. 1) then
            codint = 1
            nint = nbma
!
        else
            do k = 1, ntrf
                lint(k) = ltrf(k)
                nint = ntrf
            end do
        end if
!
    else
        if (codtrf .eq. 1) then
            do k = 1, ncm1
                lint(k) = lcm1(k)
                nint = ncm1
            end do
!
        else
!            -- ON NE PEUT PLUS RECULER, IL FAUT CALCULER
!               L'INTERSECTION :
            do k = 1, nbma
                lint(k) = 0
            end do
            do k = 1, ncm1
                lint(lcm1(k)) = 1
            end do
            do k = 1, ntrf
                lint(ltrf(k)) = lint(ltrf(k))+1
            end do
!          -- LES MAILLES COMMUNES CONTIENNENT 2 (1+1) :
            nint = 0
            do k = 1, nbma
                if (lint(k) .eq. 2) then
                    nint = nint+1
                    lint(nint) = k
                end if
            end do
        end if
    end if
!
!
    ASSERT(codint .eq. 1 .or. codint .eq. 3)
    ASSERT(nint .ge. 0)
end subroutine
