! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
function csevl(x, cs, n)
!
! UTILISEE SOUS NT POUR L'EVALUATION DE LA FONCTION D'ERREUR
! ERFC (PROVIENT DE LA BIBLIOTHEQUE SLATEC)
!
!***BEGIN PROLOGUE  CSEVL
!***PURPOSE  EVALUATE A CHEBYSHEV SERIES.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
!***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  FULLERTON, W., (LANL)
!***DESCRIPTION
!
!  EVALUATE THE N-TERM CHEBYSHEV SERIES CS AT X.  ADAPTED FROM
!  A METHOD PRESENTED IN THE PAPER BY BROUCKE REFERENCED BELOW.
!
!       INPUT ARGUMENTS --
!  X    VALUE AT WHICH THE SERIES IS TO BE EVALUATED.
!  CS   ARRAY OF N TERMS OF A CHEBYSHEV SERIES.  IN EVALUATING
!       CS, ONLY HALF THE FIRST COEFFICIENT IS SUMMED.
!  N    NUMBER OF TERMS IN ARRAY CS.
!
!***REFERENCES  R. BROUCKE, TEN SUBROUTINES FOR THE MANIPULATION OF
!                 CHEBYSHEV SERIES, ALGORITHM 446, COMMUNICATIONS OF
!                 THE A.C.M. 16, (1973) PP. 254-256.
!               L. FOX AND I. B. PARKER, CHEBYSHEV POLYNOMIALS IN
!                 NUMERICAL ANALYSIS, OXFORD UNIVERSITY PRESS, 1968,
!                 PAGE 56.
!***END PROLOGUE  DCSEVL
    implicit none
    real(kind=8) :: csevl
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
    real(kind=8) :: b0, b1, b2, cs(*), onepl, twox, x
    aster_logical :: first
    save first, onepl
!-----------------------------------------------------------------------
    integer :: i, n, ni
!-----------------------------------------------------------------------
    data first /.true./
!***FIRST EXECUTABLE STATEMENT  DCSEVL
    if (first) onepl = 1.0d0 + r8prem()
    first = .false.
    ASSERT(n .ge. 1)
    ASSERT(n .le. 1000)
    ASSERT(abs(x) .le. onepl)
!
    b1 = 0.0d0
    b0 = 0.0d0
    twox = 2.0d0*x
    do i = 1, n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
    end do
!
    csevl = 0.5d0*(b0-b2)
!
end function
