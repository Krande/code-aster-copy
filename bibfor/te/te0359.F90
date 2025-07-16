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
subroutine te0359(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/eiangl.h"
#include "asterfort/eimatb.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ngpipe.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
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
! Options: PILO_PRED_ELAS
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
    character(len=16), parameter :: typilo = 'PRED_ELAS       '
! --------------------------------------------------------------------------------------------------
    character(len=8) :: typmod(2), lielrf(10)
    aster_logical :: axi
    integer(kind=8) :: nno1, nno2, npg, lgpg, ndim, iret, ntrou, jtab(7), neps, nddl
    integer(kind=8) :: jv_w, jv_vff1, jv_geom, jv_mate, jv_vff2, jv_dff2, jv_dtau
    integer(kind=8) :: jv_copilo, jv_varim, jv_ddlm, jv_ddld, jv_ddl0, jv_ddl1, jv_angmas
    integer(kind=8) :: jv_borne
    real(kind=8):: etamin, etamax
    real(kind=8), allocatable:: wg(:, :), ni2ldc(:, :), b(:, :, :), ang(:, :)
    real(kind=8), allocatable:: sigm(:, :)
    character(len=16), pointer :: compor(:) => null()
! --------------------------------------------------------------------------------------------------
!
!
! - Get element parameters
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno1, npg=npg, &
                     jpoids=jv_w, jvf=jv_vff1)
    call elrefe_info(elrefe=lielrf(2), fami='RIGI', ndim=ndim, nno=nno2, npg=npg, &
                     jpoids=jv_w, jvf=jv_vff2, jdfde=jv_dff2)

    ndim = ndim+1
    nddl = ndim*(2*nno1+nno2)
    neps = 2*ndim

    allocate (ang(merge(1, 3, ndim .eq. 2), nno2), b(neps, npg, nddl), wg(neps, npg), &
              ni2ldc(neps, npg))
    allocate (sigm(neps, npg))

! - Type of finite element
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('S', 'TYPMOD2', typmod(2))
    axi = lteatt('AXIS', 'OUI')

    ! Parametres
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PCAMASS', 'L', jv_angmas)
    call jevech('PMATERC', 'L', jv_mate)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PDEPLMR', 'L', jv_ddlm)
    call jevech('PDDEPLR', 'L', jv_ddld)
    call jevech('PDEPL0R', 'L', jv_ddl0)
    call jevech('PDEPL1R', 'L', jv_ddl1)
    call jevech('PVARIMR', 'L', jv_varim)
    call jevech('PCDTAU', 'L', jv_dtau)
    call jevech('PCOPILO', 'E', jv_copilo)
    call jevech('PBORNPI', 'L', jv_borne)

!    NOMBRE DE VARIABLES INTERNES
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

    ! Repere local
    if (nint(zr(jv_angmas)) .eq. -1) call utmess('F', 'JOINT1_47')
    call eiangl(ndim, nno2, zr(jv_angmas+1), ang)

    ! Calcul la matrice de la cinematique
    call eimatb(nomte, ndim, axi, nno1, nno2, npg, &
                zr(jv_w), zr(jv_vff1), zr(jv_vff2), zr(jv_dff2), zr(jv_geom), &
                ang, b, wg, ni2ldc)

    ! Bornes
    etamin = zr(jv_borne+1)
    etamax = zr(jv_borne)

    ! Pilotage
    sigm = 0
    call ngpipe(typilo, npg, neps, nddl, b, &
                ni2ldc, typmod, zi(jv_mate), compor, lgpg, &
                zr(jv_ddlm), sigm, zr(jv_varim), zr(jv_ddld), zr(jv_ddl0), &
                zr(jv_ddl1), zr(jv_dtau), etamin, etamax, zr(jv_copilo))

    deallocate (ang, b, wg, ni2ldc, sigm)
end subroutine
