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
subroutine te0459(option, nomte)
!
    use HHO_type
    use HHO_size_module, only: hhoMecaFaceDofs
    use HHO_quadrature_module
    use HHO_Neumann_module
    use HHO_init_module, only: hhoInfoInitFace
    use HHO_eval_module
    use HHO_utils_module
    use HHO_geometry_module
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/tefrep.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/writeVector.h"
!
    character(len=16), intent(in) :: option, nomte
!
!---------------------------------------------------------------------------------------------------
!
!  HHO METHODS
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          CORRESPONDANT A UN CHARGEMENT SURFACIQUE POUR HHO
!          (LE CHARGEMENT PEUT ETRE DONNE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'CHAR_MECA_PRES_R'
!                    'CHAR_MECA_PRES_F'
!                    'CHAR_MECA_FF2D3D'
!                    'CHAR_MECA_FR2D3D'
!                    'CHAR_MECA_FF1D2D'
!                    'CHAR_MECA_FR1D2D'
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!---------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: maxpara = 4
    real(kind=8) :: valpar(maxpara)
    character(len=8) :: nompar(maxpara)
    type(HHO_Data) :: hhoData
    type(HHO_Face) :: hhoFace
    type(HHO_Quadrature) :: hhoQuadFace
    real(kind=8) :: rhs_forces(MSIZE_FACE_VEC), NeumValuesQP(3, MAX_QP_FACE), PresQP(MAX_QP_FACE)
    real(kind=8) :: normal(3)
    integer(kind=8) :: fbs, celldim, ipg, nbpara, idim
    integer(kind=8) :: j_time, j_pres, j_forc
!
! -- Retrieve HHO informations
!
    call hhoInfoInitFace(hhoFace, hhoData, hhoQuadFace=hhoQuadFace)
!
    ASSERT(hhoQuadFace%nbQuadPoints <= MAX_QP_FACE)
!
    celldim = hhoFace%ndim+1
    PresQP = 0.d0
    NeumValuesQP = 0.d0
    nompar(:) = 'XXXXXXXX'
    valpar(:) = 0.d0
!
! ---- Which option ?
!
    if (option .eq. 'CHAR_MECA_PRES_R') then
!
! ----- Evaluate the function PRES
!
        call jevech('PPRESSR', 'L', j_pres)
        call hhoFuncRScalEvalQp(hhoFace, hhoQuadFace, zr(j_pres), PresQP)
!
! ---- Compute the load at the quadrature points T = -p*normal
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
            normal = hhoNormalFaceQP(hhoFace, hhoQuadFace%points_param(1:2, ipg))
            NeumValuesQP(1:3, ipg) = -PresQP(ipg)*normal
        end do
!
    else if (option .eq. 'CHAR_MECA_PRES_F') then
!
! ---- Get Function Parameters
!
        if (celldim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (celldim == 2) then
            nbpara = 3
            nompar(1:2) = (/'X', 'Y'/)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ---- Time
!
        call jevech('PINSTR', 'L', j_time)
        nompar(nbpara) = 'INST'
        valpar(nbpara) = zr(j_time)
!
! ----- Evaluate the analytical function PRES
!
        call jevech('PPRESSF', 'L', j_pres)
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_pres), nbpara, nompar, valpar, &
                                celldim, PresQp)
!
! ---- Compute the load at the quadrature points T = -p*normal
        do ipg = 1, hhoQuadFace%nbQuadPoints
            normal = hhoNormalFaceQP(hhoFace, hhoQuadFace%points_param(1:2, ipg))
            NeumValuesQP(1:3, ipg) = -PresQP(ipg)*normal
        end do
!
    else if (option .eq. 'CHAR_MECA_FF2D3D' .or. option .eq. 'CHAR_MECA_FF1D2D') then
!
! ---- Get Function Parameters
!
        if (celldim == 3) then
            ASSERT(option .eq. 'CHAR_MECA_FF2D3D')
            call jevech('PFF2D3D', 'L', j_forc)
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (celldim == 2) then
            ASSERT(option .eq. 'CHAR_MECA_FF1D2D')
            call jevech('PFF1D2D', 'L', j_forc)
            nbpara = 3
            nompar(1:2) = (/'X', 'Y'/)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ---- Time
!
        call jevech('PINSTR', 'L', j_time)
        nompar(nbpara) = 'INST'
        valpar(nbpara) = zr(j_time)
!
! ----- Evaluate the analytical function (FX,FY,FZ)
!
        do idim = 1, celldim
            call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_forc-1+idim), nbpara, nompar, valpar, &
                                    celldim, NeumValuesQP(idim, 1:MAX_QP_FACE))
        end do
!
    else if (option .eq. 'CHAR_MECA_FR1D2D' .or. option .eq. 'CHAR_MECA_FR2D3D') then
!
! ---- Get Forces
!
        if (celldim == 3) then
            ASSERT(option .eq. 'CHAR_MECA_FR2D3D')
            call tefrep(option, 'PFR2D3D', j_forc)
        else if (celldim == 2) then
            ASSERT(option .eq. 'CHAR_MECA_FR1D2D')
            call tefrep(option, 'PFR1D2D', j_forc)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ---- Compute the load at the quadrature points
!
        call hhoFuncRVecEvalQp(hhoFace, hhoQuadFace, zr(j_forc), NeumValuesQP)
!
    else
!
        ASSERT(ASTER_FALSE)
    end if
!
! ---- compute surface load
!
    call hhoMecaNeumForces(hhoFace, hhoData, hhoQuadFace, NeumValuesQP, rhs_forces)
!
! ---- number of dofs
!
    call hhoMecaFaceDofs(hhoFace, hhoData, fbs)
!
! ---- save result
!
    call writeVector('PVECTUR', fbs, rhs_forces)
!
end subroutine
