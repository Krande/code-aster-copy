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
! aslint: disable=W1306
!
subroutine dmatmc(fami, materCodeJv, time, poum, ipg, &
                  ispg, anglNaut, nbsig, dr_, &
                  l_modi_cp, di_)
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dmat3d.h"
#include "asterfort/dmatcp.h"
#include "asterfort/dmatdp.h"
#include "asterfort/lteatt.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: materCodeJv
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: ipg, ispg
    real(kind=8), intent(in) :: anglNaut(3)
    integer(kind=8), intent(in) :: nbsig
    real(kind=8), optional, intent(out) :: dr_(nbsig, nbsig), di_(nbsig, nbsig)
    aster_logical, optional, intent(in) :: l_modi_cp
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Hooke matrix for iso-parametric elements
!
! --------------------------------------------------------------------------------------------------
!
! In  fami             : Gauss family for integration point rule
! In  materCodeJv      : coded material address
! In  time             : current time
! In  poum             : '-' or '+' for parameters evaluation (previous or current temperature)
! In  ipg              : current point gauss
! In  ispg             : current "sous-point" gauss
! In  anglNaut         : nautical angles for definition of basis for non-isotropic elasticity
! In  nbsig            : number of components for stress
! Out dr               : real Hooke matrix
! Out di               : imaginary Hooke matrix
! In  l_modi_cp        : using plane strain Hooke matrix for plane stress case
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: di(nbsig, nbsig)
    real(kind=8) :: dr(nbsig, nbsig)
!
! --------------------------------------------------------------------------------------------------
!
    if (lteatt('DIM_TOPO_MAILLE', '3')) then
        ASSERT(nbsig .eq. 6)
        if (present(di_)) then
            call dmat3d(fami, materCodeJv, time, poum, ipg, &
                        ispg, anglNaut, di_=di)
        end if
        if (present(dr_)) then
            call dmat3d(fami, materCodeJv, time, poum, ipg, &
                        ispg, anglNaut, dr_=dr)
        end if

    else if (lteatt('FOURIER', 'OUI')) then
        ASSERT(nbsig .eq. 6)
        if (present(di_)) then
            call dmat3d(fami, materCodeJv, time, poum, ipg, &
                        ispg, anglNaut, di_=di)
        end if
        if (present(dr_)) then
            call dmat3d(fami, materCodeJv, time, poum, ipg, &
                        ispg, anglNaut, dr_=dr)
        end if

    else if (lteatt('C_PLAN', 'OUI')) then
        ASSERT(nbsig .eq. 4)
        if (present(l_modi_cp)) then
            ASSERT(l_modi_cp)
            if (present(di_)) then
                call dmatdp(fami, materCodeJv, time, poum, ipg, &
                            ispg, anglNaut, di_=di)
            end if
            if (present(dr_)) then
                call dmatdp(fami, materCodeJv, time, poum, ipg, &
                            ispg, anglNaut, dr_=dr)
            end if
        else
            if (present(di_)) then
                call dmatcp(fami, materCodeJv, time, poum, ipg, &
                            ispg, anglNaut, di_=di)
            end if
            if (present(dr_)) then
                call dmatcp(fami, materCodeJv, time, poum, ipg, &
                            ispg, anglNaut, dr_=dr)
            end if
        end if

    else if (lteatt('D_PLAN', 'OUI') .or. lteatt('AXIS', 'OUI')) then
        ASSERT(nbsig .eq. 4)
        if (present(di_)) then
            call dmatdp(fami, materCodeJv, time, poum, ipg, &
                        ispg, anglNaut, di_=di)
        end if
        if (present(dr_)) then
            call dmatdp(fami, materCodeJv, time, poum, ipg, &
                        ispg, anglNaut, dr_=dr)
        end if

    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (present(di_)) then
        di_ = di
    end if
    if (present(dr_)) then
        dr_ = dr
    end if
!
end subroutine
