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
    subroutine erhmv2(ds_thm, axi, deltat, dimdep, dimdef,&
                      nmec, np1, np2, n2nd, ndim, nno,&
                      nnos, npg, nddls, nddlm,&
                      dimuel, ipoids, ivf, idfde, ipoid2,&
                      ivf2, idfde2, elem_coor, fovo, deplp,&
                      deplm, sielnp, sielnm, nbcmp, biot,&
                      unsurm, fpx, fpy, frx, fry,&
                      addeme, addep1,&
                      addep2, addete, adde2nd, tm2h1v)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        integer(kind=8) :: dimuel
        integer(kind=8) :: nnos
        integer(kind=8) :: nno
        integer(kind=8) :: ndim
        integer(kind=8) :: dimdef
        integer(kind=8) :: dimdep
        aster_logical :: axi
        real(kind=8) :: deltat
        integer(kind=8) :: nmec
        integer(kind=8) :: np1
        integer(kind=8) :: np2
        integer(kind=8) :: n2nd
        integer(kind=8) :: npg
        integer(kind=8) :: nddls
        integer(kind=8) :: nddlm
        integer(kind=8) :: ipoids
        integer(kind=8) :: ivf
        integer(kind=8) :: idfde
        integer(kind=8) :: ipoid2
        integer(kind=8) :: ivf2
        integer(kind=8) :: idfde2
        real(kind=8) :: elem_coor(ndim, nno)
        real(kind=8) :: fovo(ndim)
        real(kind=8) :: deplp(nno*dimdep)
        real(kind=8) :: deplm(nno*dimdep)
        real(kind=8) :: sielnp(140)
        real(kind=8) :: sielnm(140)
        integer(kind=8) :: nbcmp
        real(kind=8) :: biot
        real(kind=8) :: unsurm
        real(kind=8) :: fpx
        real(kind=8) :: fpy
        real(kind=8) :: frx(9)
        real(kind=8) :: fry(9)
        integer(kind=8) :: addeme
        integer(kind=8) :: addep1
        integer(kind=8) :: addep2
        integer(kind=8) :: addete
        integer(kind=8) :: adde2nd
        real(kind=8) :: tm2h1v(3)
    end subroutine erhmv2
end interface
