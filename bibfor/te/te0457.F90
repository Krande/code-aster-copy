! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine te0457(option, nomte)
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
#include "asterfort/writeMatrix.h"
#include "blas/dcopy.h"
#include "blas/dsyr.h"
!
    character(len=16), intent(in) :: option, nomte
!
!---------------------------------------------------------------------------------------------------
!
!  HHO METHODS
!     BUT: CALCUL DES MATRICES ELEMENTAIRES EN THERMIQUE
!          CORRESPONDANT A ECHANGE PAROI POUR HHO
!          (LE CHARGEMENT PEUT ETRE DONNE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'RIGI_THER_COEH_R'
!                    'RIGI_THER_COEH_F'
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
    real(kind=8), dimension(MSIZE_FACE_SCAL) :: basisScalEval
    real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL) :: lhs
    real(kind=8) :: CoefHQP(MAX_QP_FACE), coeff, time
    integer :: fbs, celldim, ipg, nbpara, npg
    integer :: j_time, j_coefh
!
!
! -- Get number of Gauss points
!
    call elrefe_info(fami='RIGI', npg=npg)
!
! -- Retrieve HHO informations
!
    call hhoInfoInitFace(hhoFace, hhoData, npg, hhoQuadFace)
!
    ASSERT(hhoQuadFace%nbQuadPoints <= MAX_QP_FACE)
!
    celldim = hhoFace%ndim+1
    CoefHQP = 0.d0
    nompar(:) = 'XXXXXXXX'
    valpar(:) = 0.d0
!
    call jevech('PTEMPSR', 'L', j_time)
    time = zr(j_time)
!
! ---- Which option ?
!
    if (option .eq. 'RIGI_THER_COEH_R') then
!
! ----- Get real value COEF_H
!
        call jevech('PCOEFHR', 'L', j_coefh)
        CoefHQP(1:hhoQuadFace%nbQuadPoints) = zr(j_coefh)
!
    elseif (option .eq. 'RIGI_THER_COEH_F') then
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
        nompar(nbpara) = 'INST'
        valpar(nbpara) = time
!
! ----- Evaluate the analytical function COEF_H
!
        call jevech('PCOEFHF', 'L', j_coefh)
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_coefh), nbpara, nompar, valpar, &
                                & celldim, CoefHQP)
!
    else

        ASSERT(ASTER_FALSE)
    end if
!
! ---- number of dofs
!
    call hhoTherFaceDofs(hhoFace, hhoData, fbs)
!
! ---- Compute mass matrix
!
    lhs = 0.d0
!
    call hhoBasisFace%initialize(hhoFace)
!
! ----- Loop on quadrature point
    do ipg = 1, hhoQuadFace%nbQuadPoints
! --------- Eval basis function at the quadrature point
        call hhoBasisFace%BSEval(hhoFace, hhoQuadFace%points(1:3, ipg), 0, &
                                 hhoData%face_degree(), basisScalEval)
! --------  Eval massMat
        coeff = CoefHQP(ipg)*hhoQuadFace%weights(ipg)
        call dsyr('U', fbs, coeff, basisScalEval, 1, lhs, MSIZE_FACE_SCAL)
    end do
!
! ----- Copy the lower part
!
    call hhoCopySymPartMat('U', lhs(1:fbs, 1:fbs))
!
! ---- save result
!
    call writeMatrix('PMATTTR', fbs, fbs, ASTER_TRUE, lhs)
!
end subroutine
