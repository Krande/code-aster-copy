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
    subroutine xpeshm(nno, nnop, nnops, ndim, nddls,&
                      nddlm, npg, igeom, jpintt, jpmilt, jheavn,&
                      ivf, ipoids, idfde, ivectu, ipesa,&
                      heavt, lonch, cnset, rho, axi,&
                      yaenrm, nfiss, nfh, jfisno)
        integer(kind=8) :: ndim
        integer(kind=8) :: nnop
        integer(kind=8) :: nno
        integer(kind=8) :: nnops
        integer(kind=8) :: nddls
        integer(kind=8) :: nddlm
        integer(kind=8) :: npg
        integer(kind=8) :: igeom
        integer(kind=8) :: jpintt
        integer(kind=8) :: jpmilt
        integer(kind=8) :: jheavn
        integer(kind=8) :: ivf
        integer(kind=8) :: ipoids
        integer(kind=8) :: idfde
        integer(kind=8) :: ivectu
        integer(kind=8) :: ipesa
        integer(kind=8) :: heavt(*)
        integer(kind=8) :: lonch(10)
        integer(kind=8) :: cnset(*)
        real(kind=8) :: rho
        aster_logical :: axi
        integer(kind=8) :: yaenrm
        integer(kind=8) :: nfiss
        integer(kind=8) :: nfh
        integer(kind=8) :: jfisno
    end subroutine xpeshm
end interface 
