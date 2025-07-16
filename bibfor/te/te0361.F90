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
subroutine te0361(option, nomte)
!
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/eiangl.h"
#include "asterfort/eimatb.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ngforc.h"
#include "asterfort/terefe.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_INTERFACE
!           PLAN_INTERFACE, AXIS_INTERFACE
!
! Options: FORC_NODA, REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=8) :: lielrf(10)
    aster_logical :: axi
    integer(kind=8) :: nno1, nno2, npg, ndim, neps, nddl, ntrou, jv_vff2, jv_dff2
    integer(kind=8) :: jv_w, jv_vff1, jv_geom, jv_sief, jv_vectu, jv_angmas
    real(kind=8) :: sigref, depref
    real(kind=8), allocatable:: wg(:, :), ni2ldc(:, :), b(:, :, :), ang(:, :)
    real(kind=8), allocatable:: sref(:)
! --------------------------------------------------------------------------------------------------

    ! Get element parameters
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno1, npg=npg, &
                     jpoids=jv_w, jvf=jv_vff1)
    call elrefe_info(elrefe=lielrf(2), fami='RIGI', ndim=ndim, nno=nno2, npg=npg, &
                     jpoids=jv_w, jvf=jv_vff2, jdfde=jv_dff2)
    ndim = ndim+1
    neps = 2*ndim
    nddl = ndim*(2*nno1+nno2)
    axi = lteatt('AXIS', 'OUI')

    allocate (ang(merge(1, 3, ndim .eq. 2), nno2), b(neps, npg, nddl), wg(neps, npg), &
              ni2ldc(neps, npg))
    allocate (sref(neps))

    ! Parametres communs aux deux options
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PCAMASS', 'L', jv_angmas)
    call jevech('PVECTUR', 'E', jv_vectu)

    ! Repere local nodal (par propagation du repere local constant par maille)
    if (nint(zr(jv_angmas)) .eq. -1) call utmess('F', 'JOINT1_47')
    call eiangl(ndim, nno2, zr(jv_angmas+1), ang)

    ! Matrice cinematique
    call eimatb(nomte, ndim, axi, nno1, nno2, npg, &
                zr(jv_w), zr(jv_vff1), zr(jv_vff2), zr(jv_dff2), zr(jv_geom), &
                ang, b, wg, ni2ldc)

    ! Calcul des forces nodales
    if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', jv_sief)
        call ngforc(wg, b, ni2ldc, zr(jv_sief), zr(jv_vectu))

    elseif (option .eq. 'REFE_FORC_NODA') then
        call terefe('SIGM_REFE', 'MECA_INTERFACE', sigref)
        call terefe('DEPL_REFE', 'MECA_INTERFACE', depref)
        sref(1:ndim) = sigref/ndim
        sref(ndim+1:neps) = depref
        call ngforc(wg, abs(b), ni2ldc, transpose(spread(sref, 1, npg)), zr(jv_vectu))

    else
        ASSERT(ASTER_FALSE)
    end if

    deallocate (ang, b, wg, ni2ldc)
    deallocate (sref)

end subroutine
