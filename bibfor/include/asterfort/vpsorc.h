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
    subroutine vpsorc(lmasse, ldynfa, nbeq, nbvect, nfreq,&
                      tolsor, vect, resid, workd, workl,&
                      lonwl, selec, dsor, sigma, vaux,&
                      workv, rwork, ddlexc, ddllag, neqact,&
                      maxitr, ifm, niv, priram, alpha,&
                      nconv, flage, solveu)
        integer(kind=8) :: nbvect
        integer(kind=8) :: nbeq
        integer(kind=8) :: lmasse
        integer(kind=8) :: ldynfa
        integer(kind=8) :: nfreq
        real(kind=8) :: tolsor
        complex(kind=8) :: vect(nbeq, *)
        complex(kind=8) :: resid(*)
        complex(kind=8) :: workd(*)
        complex(kind=8) :: workl(*)
        integer(kind=8) :: lonwl
        aster_logical :: selec(nbvect)
        complex(kind=8) :: dsor(*)
        complex(kind=8) :: sigma
        complex(kind=8) :: vaux(*)
        complex(kind=8) :: workv(*)
        real(kind=8) :: rwork(*)
        integer(kind=8) :: ddlexc(nbeq)
        integer(kind=8) :: ddllag(nbeq)
        integer(kind=8) :: neqact
        integer(kind=8) :: maxitr
        integer(kind=8) :: ifm
        integer(kind=8) :: niv
        integer(kind=8) :: priram(8)
        real(kind=8) :: alpha
        integer(kind=8) :: nconv
        aster_logical :: flage
        character(len=19) :: solveu
    end subroutine vpsorc
end interface
