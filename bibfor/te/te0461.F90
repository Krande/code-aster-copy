! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine te0461(option, nomte)
!
    use HHO_type
    use HHO_basis_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_Neumann_module
    use HHO_init_module, only: hhoInfoInitFace
    use HHO_eval_module
    use HHO_utils_module
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/readVector.h"
#include "asterfort/writeVector.h"
#include "blas/dcopy.h"
#include "blas/daxpy.h"
!
    character(len=16), intent(in) :: option, nomte
!
!---------------------------------------------------------------------------------------------------
!
!  HHO METHODS
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN THERMIQUE
!          SUR UNE FACE POUR HHO
!          (LE CHARGEMENT PEUT ETRE DONNE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'CHAR_THER_TEXT_R'
!                    'CHAR_THER_TEXT_F'
!                    'CHAR_THER_FLUN_R'
!                    'CHAR_THER_FLUN_F'
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!---------------------------------------------------------------------------------------------------
!
    integer, parameter :: maxpara = 4
    real(kind=8) :: valpar(maxpara)
    character(len=8) :: nompar(maxpara)
    type(HHO_Data) :: hhoData
    type(HHO_Face) :: hhoFace
    type(HHO_Quadrature) :: hhoQuadFace
    type(HHO_basis_face) :: hhoBasisFace
    real(kind=8), dimension(MSIZE_FACE_SCAL) :: rhs, temp_F_curr
    real(kind=8) :: CoefHQP_curr(MAX_QP_FACE), CoefHQP_prev(MAX_QP_FACE)
    real(kind=8) :: ParaQP_curr(MAX_QP_FACE), ParaQP_prev(MAX_QP_FACE)
    real(kind=8) :: ValQP_curr(MAX_QP_FACE), ValQP_prev(MAX_QP_FACE)
    real(kind=8) :: NeumValuesQP(MAX_QP_FACE)
    real(kind=8) :: theta, time_curr, time_prev, temp_eval_curr
    integer :: fbs, celldim, ipg, nbpara, npg
    integer :: j_time, j_coefh, j_para
!
!
! -- Get number of Gauss points
!
    call elrefe_info(fami='RIGI', npg=npg)
!
! -- Retrieve HHO informations
!
    call hhoInfoInitFace(hhoFace, hhoData, npg, hhoQuadFace)
    call hhoBasisFace%initialize(hhoFace)
!
! ---- number of dofs
!
    call hhoTherFaceDofs(hhoFace, hhoData, fbs)
!
    ASSERT(hhoQuadFace%nbQuadPoints <= MAX_QP_FACE)
!
    celldim = hhoFace%ndim+1
    CoefHQP_curr = 0.d0
    CoefHQP_prev = 0.d0
    ParaQP_curr = 0.d0
    ParaQP_prev = 0.d0
    ValQP_curr = 0.d0
    ValQP_prev = 0.d0
    nompar(:) = 'XXXXXXXX'
    valpar(:) = 0.d0
!
    call jevech('PTEMPSR', 'L', j_time)
    time_curr = zr(j_time)
    time_prev = time_curr-zr(j_time+1)
    theta = zr(j_time+2)
!
! ---- Which option ?
!
    if (option .eq. 'CHAR_THER_TEXT_R') then
!
! ----- Get real value COEF_H
!
        call jevech('PCOEFHR', 'L', j_coefh)
        CoefHQP_prev(1:hhoQuadFace%nbQuadPoints) = zr(j_coefh)
        CoefHQP_curr = CoefHQP_prev
!
! ----- Get real value Text
!
        call jevech('PT_EXTR', 'L', j_para)
        ParaQP_prev(1:hhoQuadFace%nbQuadPoints) = zr(j_para)
        ParaQP_curr = ParaQP_prev
!
    elseif (option .eq. 'CHAR_THER_TEXT_F') then
        call jevech('PCOEFHF', 'L', j_coefh)
        call jevech('PT_EXTF', 'L', j_para)
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
! ---- Time +
!
        nompar(nbpara) = 'INST'
        valpar(nbpara) = time_curr
!
! ----- Evaluate the analytical function at T+
!
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_coefh), nbpara, nompar, valpar, &
                                & celldim, CoefHQP_curr)
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para), nbpara, nompar, valpar, &
                                & celldim, ParaQP_curr)
!
! ---- Time -
!
        nompar(nbpara) = 'INST'
        valpar(nbpara) = time_prev
!
! ----- Evaluate the analytical function at T-
!
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_coefh), nbpara, nompar, valpar, &
                                & celldim, CoefHQP_prev)
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para), nbpara, nompar, valpar, &
                                & celldim, ParaQP_prev)
!
    elseif (option .eq. 'CHAR_THER_FLUN_R') then
!
! ----- Get real value FLUXN
!
        call jevech('PFLUXNR', 'L', j_para)
        ParaQP_prev(1:hhoQuadFace%nbQuadPoints) = zr(j_para)
        ParaQP_curr = ParaQP_prev
!
    elseif (option .eq. 'CHAR_THER_FLUN_F') then
        call jevech('PFLUXNF', 'L', j_para)
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
! ---- Time +
!
        nompar(nbpara) = 'INST'
        valpar(nbpara) = time_curr
!
! ----- Evaluate the analytical function at T+
!
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para), nbpara, nompar, valpar, &
                                & celldim, ParaQP_curr)
!
! ---- Time -
!
        nompar(nbpara) = 'INST'
        valpar(nbpara) = time_prev
!
! ----- Evaluate the analytical function at T-
!
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para), nbpara, nompar, valpar, &
                                & celldim, ParaQP_prev)
!
    else

        ASSERT(ASTER_FALSE)
    end if
!
    if (option(1:15) .eq. 'CHAR_THER_TEXT_') then
        ValQP_curr = CoefHQP_curr*ParaQP_curr
!
        call readVector('PTEMPER', fbs, temp_F_curr)
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
            temp_eval_curr = hhoEvalScalFace(hhoFace, hhoBasisFace, hhoData%face_degree(), &
                                             hhoQuadFace%points(1:3, ipg), temp_F_curr, fbs)

            ValQP_prev(ipg) = CoefHQP_prev(ipg)*(ParaQP_prev(ipg)-temp_eval_curr)
        end do
!
    elseif (option(1:15) .eq. 'CHAR_THER_FLUN_') then
        ValQP_curr = ParaQP_curr
        ValQP_prev = ParaQP_prev
    else
        ASSERT(ASTER_FALSE)
    end if
!
! ---- compute surface load
!
    NeumValuesQP = theta*ValQP_curr+(1.d0-theta)*ValQP_prev
    call hhoTherNeumForces(hhoFace, hhoData, hhoQuadFace, NeumValuesQP, rhs)
!
! ---- save result
!
    call writeVector('PVECTTR', fbs, rhs)
!
end subroutine
