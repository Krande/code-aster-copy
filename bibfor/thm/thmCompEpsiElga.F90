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
subroutine thmCompEpsiElga(ds_thm)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/thmGetElemInfo.h"
#include "asterfort/thmGetElemDime.h"
#include "asterfort/epsthm.h"
#include "asterfort/thmGetElemRefe.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/thmGetGene.h"
#include "asterfort/thmGetElemIntegration.h"
!
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM - Compute
!
! EPSI_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elrefe, elref2
    integer(kind=8) :: addeme, addep1, addep2, addete, adde2nd
    integer(kind=8) :: ipg, i_cmp
    integer(kind=8) :: jv_geom, jv_disp, jv_strain
    real(kind=8) :: epsm(6, 27)
    integer(kind=8) :: nno, nnos, nnom
    integer(kind=8) :: npi, npi2, npg
    integer(kind=8) :: jv_poids, jv_poids2
    integer(kind=8) :: jv_func, jv_func2, jv_dfunc, jv_dfunc2, jv_gano
    integer(kind=8) :: nddls, nddlm
    integer(kind=8) :: nddl_meca, nddl_p1, nddl_p2, nddl_2nd
    integer(kind=8) :: dimdep, dimdef, dimcon, dimuel
    aster_logical :: l_axi, l_vf
    character(len=3) :: inte_type
    integer(kind=8) :: ndim
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
!
! --------------------------------------------------------------------------------------------------
!

!
! - Get model of finite element
!
    call thmGetElemModel(ds_thm, l_axi, l_vf, ndim)
!
! - Cannot compute for finite volume
!
    ASSERT(.not. l_vf)
!
! - Get type of integration
!
    call thmGetElemIntegration(l_vf, inte_type)
!
! - Get generalized coordinates
!
    call thmGetGene(ds_thm, l_vf, ndim, &
                    mecani, press1, press2, tempe, second)
!
! - Get reference elements
!
    call thmGetElemRefe(l_vf, elrefe, elref2)
!
! - Get informations about element
!
    call thmGetElemInfo(l_vf, elrefe, elref2, &
                        nno, nnos, nnom, &
                        jv_gano, jv_poids, jv_poids2, &
                        jv_func, jv_func2, jv_dfunc, jv_dfunc2, &
                        inte_type, npi, npi2, npg)
    ASSERT(npi .le. 27)
    ASSERT(nno .le. 20)
!
! - Address in generalized strains vector
!
    addeme = mecani(2)
    addete = tempe(2)
    addep1 = press1(3)
    addep2 = press2(3)
    adde2nd = second(2)
!
! - Get dimensions about element
!
    call thmGetElemDime(ndim, nnos, nnom, &
                        mecani, press1, press2, tempe, second, &
                        nddls, nddlm, &
                        nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                        dimdep, dimdef, dimcon, dimuel)
!
! - Input fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PDEPLAR', 'L', jv_disp)
!
! - Compute strains
!
    call epsthm(ds_thm, l_axi, ndim, &
                addeme, addep1, addep2, addete, adde2nd, &
                nno, nnos, &
                dimuel, dimdef, nddls, nddlm, &
                nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                npi, zr(jv_geom), zr(jv_disp), &
                jv_poids, jv_poids2, &
                jv_func, jv_func2, jv_dfunc, jv_dfunc2, &
                epsm)
!
! - Output field
!
    call jevech('PDEFOPG', 'E', jv_strain)
    do ipg = 1, npg
        do i_cmp = 1, 6
            zr(jv_strain+6*(ipg-1)+i_cmp-1) = epsm(i_cmp, ipg)
        end do
    end do
!
end subroutine
