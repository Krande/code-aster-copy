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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine te0018(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/evalPressure.h"
#include "asterfort/nmpr3d_vect.h"
#include "asterfort/mb_pres.h"
#include "asterfort/tecach.h"
#include "asterfort/lteatt.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D (skin elements)
!           MEMBRANE
!
! Options: CHAR_MECA_PRES_* (for 3D skin elements only)
!          CHAR_MECA_EFON_*
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: mxnoeu = 9, mxnpg = 27
    aster_logical :: l_func, l_time, l_efff
    integer(kind=8) :: jv_geom, jv_time, jv_pres, jv_effe
    integer(kind=8) :: jv_vect
    real(kind=8) :: time
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: nno, npg, ndim, ndofbynode
    integer(kind=8) :: iret, kpg
    real(kind=8) :: pres, pres_pg(mxnpg), coef_mult
!
! --------------------------------------------------------------------------------------------------
!
    pres_pg = 0.d0
    l_func = (option .eq. 'CHAR_MECA_PRES_F') .or. (option .eq. 'CHAR_MECA_EFON_F')
    l_efff = (option .eq. 'CHAR_MECA_EFON_R') .or. (option .eq. 'CHAR_MECA_EFON_F')
!
! - For membrane in small strain, we allow pressure only if it is null
!
    if (lteatt('TYPMOD', 'MEMBRANE')) then
        call mb_pres()
    end if
!
! - Input fields: for pressure, no node affected -> 0
!
    call jevech('PGEOMER', 'L', jv_geom)
    if (l_func) then
        if (l_efff) then
            call jevecd('PPREFFF', jv_pres, 0.d0)
        else
            call jevecd('PPRESSF', jv_pres, 0.d0)
        end if
    else
        if (l_efff) then
            call jevecd('PPREFFR', jv_pres, 0.d0)
        else
            call jevecd('PPRESSR', jv_pres, 0.d0)
        end if
    end if
    if (l_efff) then
        call jevech('PEFOND', 'L', jv_effe)
    end if
!
! - Get time if present
!
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jv_time)
    l_time = ASTER_FALSE
    time = 0.d0
    if (jv_time .ne. 0) then
        l_time = ASTER_TRUE
        time = zr(jv_time)
    end if
!
! - Output fields
!
    call jevech('PVECTUR', 'E', jv_vect)
!
! - Get element parameters
!
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, ndim=ndim, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .le. mxnoeu)
    ASSERT(npg .le. mxnpg)
!
! - Pressure are on skin elements but DOF are volumic
!
    ASSERT(ndim .eq. 2)
    ndofbynode = ndim+1
!
! - Multiplicative ratio for pressure (EFFE_FOND)
!
    coef_mult = 1.d0
    if (l_efff) then
        coef_mult = zr(jv_effe-1+1)
    end if
!
! - Evaluation of pressure at Gauss points (from nodes)
!
    do kpg = 1, npg
        call evalPressure(l_func, l_time, time, &
                          nno, ndim, kpg, &
                          ivf, jv_geom, jv_pres, &
                          pres)
        pres_pg(kpg) = coef_mult*pres
    end do
!
! - Second member
!
    call nmpr3d_vect(nno, npg, ndofbynode, &
                     zr(ipoids), zr(ivf), zr(idfde), &
                     zr(jv_geom), pres_pg, zr(jv_vect))
!
end subroutine
