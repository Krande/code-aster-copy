! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine metau1(l_meta)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/metaGetType.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/verift.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    aster_logical, intent(out) :: l_meta
!
! --------------------------------------------------------------------------------------------------
!
! Metallurgy
!
! Compute CHAR_MECA_TEMP_R - 2D case
!
! --------------------------------------------------------------------------------------------------
!
! Out l_meta : .true. if metallurgy exists
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    real(kind=8) :: young, nu, epsth
    real(kind=8) :: dfdx(MT_NNOMAX2D), dfdy(MT_NNOMAX2D), poids, r
    integer(kind=8) :: nbNode, kpg, npg, iNode
    integer(kind=8) :: metaType, nbPhases
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: jvGeom, jvMater, jvMaterCode, j_vect
    integer(kind=8) :: elasID
    character(len=16) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!
    l_meta = ASTER_TRUE

! - Get metallurgy type
    call metaGetType(metaType, nbPhases)
    ASSERT(nbPhases .le. 5)
    if (metaType .eq. META_NONE) then
        l_meta = ASTER_FALSE

    else
! ----- Finite element informations
        call elrefe_info(fami='RIGI', nno=nbNode, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde)
        ASSERT(nbNode .le. MT_NNOMAX2D)

! ----- Geometry
        call jevech('PGEOMER', 'L', jvGeom)

! ----- Material parameters
        call jevech('PMATERC', 'L', jvMater)
        jvMaterCode = zi(jvMater)

! ----- Output field
        call jevech('PVECTUR', 'E', j_vect)

        do kpg = 1, npg

! --------- Shape functions derivatives
            call dfdm2d(nbNode, kpg, ipoids, idfde, zr(jvGeom), &
                        poids, dfdx, dfdy)

! --------- Axi-symmetric case
            if (lteatt('AXIS', 'OUI')) then
                r = 0.d0
                do iNode = 1, nbNode
                    r = r+zr(jvGeom+2*(iNode-1))*zr(ivf+(kpg-1)*nbNode+iNode-1)
                end do
                poids = poids*r
                do iNode = 1, nbNode
                    dfdx(iNode) = dfdx(iNode)+zr(ivf+(kpg-1)*nbNode+iNode-1)/r
                end do
            end if

! --------- Compute thermic strain
            call verift('RIGI', kpg, ksp, '+', jvMaterCode, &
                        epsth_meta_=epsth)

! --------- Get elastic parameters
            call get_elas_id(jvMaterCode, elasID, elasKeyword)
            call get_elas_para('RIGI', jvMaterCode, '+', kpg, ksp, &
                               elasID, elasKeyword, &
                               e_=young, nu_=nu)
            ASSERT(elasID .eq. ELAS_ISOT)

! --------- Compute
            poids = poids*(young/(1.d0-2.d0*nu))*epsth
            do iNode = 1, nbNode
                zr(j_vect+2*iNode-2) = zr(j_vect+2*iNode-2)+poids*dfdx(iNode)
                zr(j_vect+2*iNode-1) = zr(j_vect+2*iNode-1)+poids*dfdy(iNode)
            end do
        end do
!
    end if
end subroutine
