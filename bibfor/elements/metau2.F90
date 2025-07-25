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
subroutine metau2(l_meta)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/metaGetType.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/verift.h"
#include "asterfort/Metallurgy_type.h"
!
    aster_logical, intent(out) :: l_meta
!
! --------------------------------------------------------------------------------------------------
!
! Metallurgy
!
! Compute CHAR_MECA_TEMP_R - 3D case
!
! --------------------------------------------------------------------------------------------------
!
! Out l_meta : .true. if metallurgy exists
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: coef
    real(kind=8) :: young, nu
    real(kind=8) :: epsth
    real(kind=8) :: dfdx(27), dfdy(27), dfdz(27)
    real(kind=8) :: poids
    integer(kind=8) :: nb_node, ispg, kp, npg, i_node, elas_id
    integer(kind=8) :: meta_type, nb_phasis
    integer(kind=8) :: ipoids, idfde
    integer(kind=8) :: j_geom, j_mate, j_mater, j_vect
    character(len=16) :: elas_keyword
!
! --------------------------------------------------------------------------------------------------
!
    l_meta = .true.
    ispg = 1
!
! - Get metallurgy type
!
    call metaGetType(meta_type, nb_phasis)
    if (meta_type .eq. META_NONE) then
        l_meta = .false.
        goto 999
    end if
!
! - Finite element informations
!
    call elrefe_info(fami='RIGI', nno=nb_node, npg=npg, jpoids=ipoids, jdfde=idfde)
!
! - Geometry
!
    call jevech('PGEOMER', 'L', j_geom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', j_mate)
!
! - Coded material address
!
    j_mater = zi(j_mate)
!
! - Output field
!
    call jevech('PVECTUR', 'E', j_vect)
!
    do kp = 1, npg
!
! ----- Shape functions derivatives
!
        call dfdm3d(nb_node, kp, ipoids, idfde, zr(j_geom), &
                    poids, dfdx, dfdy, dfdz)
!
! ----- Compute thermic strain
!
        call verift('RIGI', kp, 1, '+', j_mater, &
                    epsth_meta_=epsth)
!
! ----- Get elastic parameters
!
        call get_elas_id(j_mater, elas_id, elas_keyword)
        call get_elas_para('RIGI', j_mater, '+', kp, ispg, &
                           elas_id, elas_keyword, &
                           e_=young, nu_=nu)
        ASSERT(elas_id .eq. 1)
!
! ----- Compute
!
        coef = young/(1.d0-2.d0*nu)
        poids = poids*coef*epsth
!
        do i_node = 1, nb_node
            zr(j_vect+3*i_node-3) = zr(j_vect+3*i_node-3)+poids*dfdx(i_node)
            zr(j_vect+3*i_node-2) = zr(j_vect+3*i_node-2)+poids*dfdy(i_node)
            zr(j_vect+3*i_node-1) = zr(j_vect+3*i_node-1)+poids*dfdz(i_node)
        end do
    end do
!
999 continue
end subroutine
