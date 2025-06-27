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
!
#include "asterf_types.h"
!
interface
    subroutine wpsorn(appr, lmasse, lamor, lmatra, nbeq,&
                      nbvect, nfreq, tolsor, vect, resid,&
                      workd, workl, lonwl, selec, dsor,&
                      vpr, vpi, sigma, vaux, workv,&
                      ddlexc, ddllag, neqact, maxitr, ifm,&
                      niv, priram, alpha, nconv, flage,&
                      vaur, vauc, vaul, solveu)
        integer(kind=8) :: nfreq
        integer(kind=8) :: nbeq
        character(len=1) :: appr
        integer(kind=8) :: lmasse
        integer(kind=8) :: lamor
        integer(kind=8) :: lmatra
        integer(kind=8) :: nbvect
        real(kind=8) :: tolsor
        complex(kind=8) :: vect(nbeq, *)
        real(kind=8) :: resid(*)
        real(kind=8) :: workd(*)
        real(kind=8) :: workl(*)
        integer(kind=8) :: lonwl
        aster_logical :: selec(*)
        real(kind=8) :: dsor(nfreq+1, *)
        real(kind=8) :: vpr(*)
        real(kind=8) :: vpi(*)
        complex(kind=8) :: sigma
        real(kind=8) :: vaux(*)
        real(kind=8) :: workv(*)
        integer(kind=8) :: ddlexc(*)
        integer(kind=8) :: ddllag(*)
        integer(kind=8) :: neqact
        integer(kind=8) :: maxitr
        integer(kind=8) :: ifm
        integer(kind=8) :: niv
        integer(kind=8) :: priram(8)
        real(kind=8) :: alpha
        integer(kind=8) :: nconv
        aster_logical :: flage
        real(kind=8) :: vaur(2*nbeq, *)
        complex(kind=8) :: vauc(2*nbeq, *)
        complex(kind=8) :: vaul(2*nbeq, *)
        character(len=19) :: solveu
    end subroutine wpsorn
end interface
