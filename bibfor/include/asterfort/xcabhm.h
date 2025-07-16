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
#include "asterf_types.h"
!
interface 
    subroutine xcabhm(ds_thm,&
                      nddls, nddlm, nnop, nnops, nnopm,&
                      dimuel, ndim, kpi, ff, ff2,&
                      dfdi, dfdi2, b, nmec,&
                      addeme, addep1, np1, axi,&
                      ivf, ipoids, idfde, poids, coorse,&
                      nno, geom, yaenrm, adenme, dimenr,&
                      he, heavn, yaenrh, adenhy, nfiss, nfh)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        integer(kind=8) :: dimenr
        integer(kind=8) :: ndim
        integer(kind=8) :: dimuel
        integer(kind=8) :: nnops
        integer(kind=8) :: nnop
        integer(kind=8) :: nddls
        integer(kind=8) :: nddlm
        integer(kind=8) :: nnopm
        integer(kind=8) :: kpi
        real(kind=8) :: ff(nnop)
        real(kind=8) :: ff2(nnops)
        real(kind=8) :: dfdi(nnop, ndim)
        real(kind=8) :: dfdi2(nnops, ndim)
        real(kind=8) :: b(dimenr, dimuel)
        integer(kind=8) :: nmec
        integer(kind=8) :: addeme
        integer(kind=8) :: addep1
        integer(kind=8) :: np1
        aster_logical :: axi
        integer(kind=8) :: ivf
        integer(kind=8) :: ipoids
        integer(kind=8) :: idfde
        integer(kind=8) :: heavn(nnop,5)
        real(kind=8) :: poids
        real(kind=8) :: coorse(81)
        integer(kind=8) :: nno
        real(kind=8) :: geom(ndim, nnop)
        integer(kind=8) :: yaenrm
        integer(kind=8) :: adenme
        real(kind=8) :: he(nfiss)
        integer(kind=8) :: yaenrh
        integer(kind=8) :: adenhy
        integer(kind=8) :: nfiss
        integer(kind=8) :: nfh
    end subroutine xcabhm
end interface 
