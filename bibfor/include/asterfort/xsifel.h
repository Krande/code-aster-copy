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
    subroutine xsifel(elrefp, ndim, coorse, igeom, jheavt,&
                      ise, nfh, ddlc, ddlm, nfe,&
                      puls, basloc, nnop,&
                      idepl, lsn, lst, idecpg, igthet,&
                      fno, nfiss, jheavn, jstno)
        integer(kind=8) :: nfiss
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        character(len=8) :: elrefp
        real(kind=8) :: coorse(*)
        integer(kind=8) :: igeom
        integer(kind=8) :: jheavt
        integer(kind=8) :: ise
        integer(kind=8) :: nfh
        integer(kind=8) :: ddlc
        integer(kind=8) :: ddlm
        integer(kind=8) :: nfe
        real(kind=8) :: rho
        real(kind=8) :: puls
        aster_logical :: lmoda
        real(kind=8) :: basloc(3*ndim*nnop)
        integer(kind=8) :: idepl
        real(kind=8) :: lsn(nnop)
        real(kind=8) :: lst(nnop)
        integer(kind=8) :: idecpg
        integer(kind=8) :: igthet
        real(kind=8) :: fno(ndim*nnop)
        integer(kind=8) :: jheavn
        integer(kind=8) :: jstno
    end subroutine xsifel
end interface
