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
! person_in_charge: mickael.abbas at edf.fr
!
module HHO_L2proj_module
!
    use HHO_rhs_module
    use HHO_massmat_module
    use HHO_quadrature_module
    use HHO_type
    use HHO_size_module
    use HHO_eval_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
#include "blas/dposv.h"
#include "blas/dcopy.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Module to compute L2 prjoection
!
! --------------------------------------------------------------------------------------------------
    public :: hhoL2ProjFaceScal, hhoL2ProjFaceVec, hhoL2ProjScal, hhoL2ProjVec
    public :: hhoL2ProjCellScal, hhoL2ProjCellVec, hhoL2ProjFieldScal
!    private  ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoL2ProjFaceScal(hhoFace, hhoQuad, FuncValuesQP, degree, coeff_L2Proj)
!
        implicit none
!
        type(HHO_Face), intent(in)          :: hhoFace
        type(HHO_Quadrature), intent(in)    :: hhoQuad
        real(kind=8), intent(in)            :: FuncValuesQP(MAX_QP_FACE)
        integer, intent(in)                 :: degree
        real(kind=8), intent(out)           :: coeff_L2Proj(MSIZE_FACE_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the L2-prjoection of a scalar given function on P^k(F,R)
!   In hhoFace      : the current HHO Face
!   In hhoQuad      : Quadrature for the face
!   In FuncValuesQP : Values of the function to project at the quadrature points
!   In degree       : degree of the projection k
!   Out coeff_L2Proj: coefficient after projection
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: faceMass(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL)
        integer :: info, mbs
! --------------------------------------------------------------------------------------------------
!
        info = 0
        if (2*degree > hhoQuad%order) then
            call utmess('F', 'HHO1_12')
        end if
!
! ----- Compute face mass matrix
!
        call hhoMassMatFaceScal(hhoFace, 0, degree, faceMass, mbs)
!
! ---- Compute rhs
!
        call hhoMakeRhsFaceScal(hhoFace, hhoQuad, FuncValuesQP, degree, coeff_L2Proj)
!
! ---- Solve the system
!
        call dposv('U', mbs, 1, faceMass, MSIZE_FACE_SCAL, coeff_L2Proj, MSIZE_FACE_SCAL, info)
!
! ---- Sucess ?
!
        if (info .ne. 0) then
            call utmess('F', 'HHO1_4')
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoL2ProjFaceVec(hhoFace, hhoQuad, FuncValuesQP, degree, coeff_L2Proj)
!
        implicit none
!
        type(HHO_Face), intent(in)          :: hhoFace
        type(HHO_Quadrature), intent(in)    :: hhoQuad
        real(kind=8), intent(in)            :: FuncValuesQP(3, MAX_QP_FACE)
        integer, intent(in)                 :: degree
        real(kind=8), intent(out)           :: coeff_L2Proj(MSIZE_FACE_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the L2-prjoection of a vectorial given function on P^k(F;R^d)
!   In hhoFace      : the current HHO Face
!   In hhoQuad      : Quadrature for the face
!   In FuncValuesQP : Values of the function to project at the quadrature points
!   In degree       : degree of the projection k
!   Out coeff_L2Proj: coefficient after projection
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: faceMass(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL)
        integer :: info, mbs
!
! --------------------------------------------------------------------------------------------------
!
        info = 0
        if (2*degree > hhoQuad%order) then
            call utmess('F', 'HHO1_12')
        end if
!
! ----- Compute face mass matrix
!
        call hhoMassMatFaceScal(hhoFace, 0, degree, faceMass, mbs)
!
! ---- Compute rhs
!
        call hhoMakeRhsFaceVec(hhoFace, hhoQuad, FuncValuesQP, degree, coeff_L2Proj)
!
! ---- Solve the system
!
        call dposv('U', mbs, hhoFace%ndim+1, faceMass, MSIZE_FACE_SCAL, coeff_L2Proj, mbs, info)
!
! ---- Sucess ?
!
        if (info .ne. 0) then
            call utmess('F', 'HHO1_4')
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoL2ProjCellScal(hhoCell, hhoQuad, FuncValuesQP, degree, coeff_L2Proj)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Quadrature), intent(in)    :: hhoQuad
        real(kind=8), intent(in)            :: FuncValuesQP(MAX_QP_CELL)
        integer, intent(in)                 :: degree
        real(kind=8), intent(out)           :: coeff_L2Proj(MSIZE_CELL_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the L2-prjoection of a scalar given function on a cell
!   In hhoCell      : the current HHO Cell
!   In hhoQuad      : Quadrature for the Cell
!   In FuncValuesQP : Values of the function to project at the quadrature points
!   In degree       : degree of the projection k
!   Out coeff_L2Proj: coefficient after projection
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: cellMass(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
        integer :: info, mbs
! --------------------------------------------------------------------------------------------------
!
        info = 0
        if (2*degree > hhoQuad%order) then
            call utmess('F', 'HHO1_12')
        end if
!
! ----- Compute Cell mass matrix
!
        call hhoMassMatCellScal(hhoCell, 0, degree, cellMass, mbs)
!
! ---- Compute rhs
!
        call hhoMakeRhsCellScal(hhoCell, hhoQuad, FuncValuesQP, degree, coeff_L2Proj)
!
! ---- Solve the system
!
        call dposv('U', mbs, 1, cellMass, MSIZE_CELL_SCAL, coeff_L2Proj, MSIZE_CELL_SCAL, info)
!
! ---- Sucess ?
!
        if (info .ne. 0) then
            call utmess('F', 'HHO1_4')
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoL2ProjCellVec(hhoCell, hhoQuad, FuncValuesQP, degree, coeff_L2Proj)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Quadrature), intent(in)    :: hhoQuad
        real(kind=8), intent(in)            :: FuncValuesQP(3, MAX_QP_CELL)
        integer, intent(in)                 :: degree
        real(kind=8), intent(out)           :: coeff_L2Proj(MSIZE_CELL_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the L2-prjoection of a vectorial given function on a cell
!   In hhoCell      : the current HHO Cell
!   In hhoQuad      : Quadrature for the cace
!   In FuncValuesQP : Values of the function to project at the quadrature points
!   In degree       : degree of the projection k
!   Out coeff_L2Proj: coefficient after projection
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: cellMass(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
        integer :: info, mbs
!
! --------------------------------------------------------------------------------------------------
!
        info = 0
        if (2*degree > hhoQuad%order) then
            call utmess('F', 'HHO1_12')
        end if
!
! ----- Compute cell mass matrix
!
        call hhoMassMatCellScal(hhoCell, 0, degree, cellMass, mbs)
!
! ---- Compute rhs
!
        call hhoMakeRhsCellVec(hhoCell, hhoQuad, FuncValuesQP, degree, coeff_L2Proj)
!
! ---- Solve the system
!
        call dposv('U', mbs, hhoCell%ndim, cellMass, MSIZE_CELL_SCAL, coeff_L2Proj, mbs, info)
!
! ---- Sucess ?
!
        if (info .ne. 0) then
            call utmess('F', 'HHO1_4')
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoL2ProjScal(hhoCell, hhoData, func, time, coeff_L2Proj, all)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Data), intent(in)          :: hhoData
        character(len=8), intent(in)        :: func
        real(kind=8), intent(in)            :: time
        real(kind=8), intent(out)           :: coeff_L2Proj(MSIZE_TDOFS_SCAL)
        aster_logical, intent(in), optional          :: all
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the L2-prjoection of a scalar given function on a local space (cell + face)
!   In hhoCell      : the current HHO Cell
!   In FuncValuesQP : Values of the function to project at the quadrature points
!   Out coeff_L2Proj: coefficient after projection
!
! --------------------------------------------------------------------------------------------------
!
        integer, parameter :: maxpara = 4
        real(kind=8) :: valpar(maxpara)
        character(len=8) :: nompar(maxpara)
        type(HHO_Face) :: hhoFace
        type(HHO_Quadrature) :: hhoQuadFace, hhoQuadCell
        integer :: cbs, fbs, total_dofs, iFace, ind, nbpara
        real(kind=8) :: FuncValuesCellQP(MAX_QP_CELL), FuncValuesFaceQP(MAX_QP_FACE)
        aster_logical :: with_faces
! --------------------------------------------------------------------------------------------------
!
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
        coeff_L2Proj = 0
        FuncValuesCellQP = 0.d0
        with_faces = ASTER_TRUE
        if (present(all)) with_faces = all
!
! --- Type of function dor a face
!
        if (hhoCell%ndim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
            nompar(nbpara) = 'INST'
            valpar(nbpara) = time
        else if (hhoCell%ndim == 2) then
            nbpara = 3
            nompar(1:2) = (/'X', 'Y'/)
            nompar(nbpara) = 'INST'
            valpar(nbpara) = time
            nompar(4) = 'XXXXXXXX'
            valpar(4) = 0.d0
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Loop on faces
!
        ind = 1
        do iFace = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iFace)
!
! ----- get quadrature
!
            if (with_faces) then
                call hhoQuadFace%GetQuadFace(hhoface, 2*hhoData%face_degree()+1)
!
! -------------- Value of the function at the quadrature point
!
                call hhoFuncFScalEvalQp(hhoQuadFace, func, nbpara, nompar, &
                                        valpar, hhoCell%ndim, FuncValuesFaceQP)
!
!
! -------------- Compute L2 projection
!
                call hhoL2ProjFaceScal(hhoFace, hhoQuadFace, FuncValuesFaceQP, &
                                       hhoData%face_degree(), coeff_L2Proj(ind))
            end if
            ind = ind+fbs
        end do
!
! --- On cell
!
        call hhoQuadCell%GetQuadCell(hhoCell, 2*hhoData%cell_degree()+1)
!
! -------------- Value of the function at the quadrature point
!
        call hhoFuncFScalEvalQp(hhoQuadCell, func, nbpara, nompar, &
                                valpar, hhoCell%ndim, FuncValuesCellQP)
!
        call hhoL2ProjCellScal(hhoCell, hhoQuadCell, FuncValuesCellQP, hhoData%cell_degree(), &
                               coeff_L2Proj(ind))
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoL2ProjVec(hhoCell, hhoData, func, time, coeff_L2Proj, all)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Data), intent(in)          :: hhoData
        character(len=8), intent(in)        :: func(*)
        real(kind=8), intent(in)            :: time
        real(kind=8), intent(out)           :: coeff_L2Proj(MSIZE_TDOFS_VEC)
        aster_logical, intent(in), optional   :: all

!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the L2-prjoection of a vectorial given function on a local space (cell + face)
!   In hhoCell      : the current HHO Cell
!   In FuncValuesQP : Values of the function to project at the quadrature points
!   Out coeff_L2Proj: coefficient after projection
!
! --------------------------------------------------------------------------------------------------
!
        integer, parameter :: maxpara = 4
        real(kind=8) :: valpar(maxpara)
        character(len=8) :: nompar(maxpara)
        type(HHO_Face) :: hhoFace
        type(HHO_Quadrature) :: hhoQuadFace, hhoQuadCell
        integer :: cbs, fbs, total_dofs, iFace, ind, nbpara, idim
        real(kind=8) :: FuncValuesCellQP(3, MAX_QP_CELL), FuncValuesFaceQP(3, MAX_QP_FACE)
        real(kind=8) :: rhs_face(MSIZE_FACE_VEC), rhs_cell(MSIZE_CELL_VEC)
        aster_logical :: with_faces
! --------------------------------------------------------------------------------------------------
!
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
        coeff_L2Proj = 0
        FuncValuesCellQP = 0.d0
        with_faces = ASTER_TRUE
        if (present(all)) with_faces = all
!
! --- Type of function dor a face
!
        if (hhoCell%ndim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
            nompar(nbpara) = 'INST'
            valpar(nbpara) = time
        else if (hhoCell%ndim == 2) then
            nbpara = 3
            nompar(1:2) = (/'X', 'Y'/)
            nompar(nbpara) = 'INST'
            valpar(nbpara) = time
            nompar(4) = 'XXXXXXXX'
            valpar(4) = 0.d0
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Loop on faces
!
        ind = 1
        do iFace = 1, hhoCell%nbfaces
            if (with_faces) then
                hhoFace = hhoCell%faces(iFace)
!
! ----- get quadrature
!
                call hhoQuadFace%GetQuadFace(hhoface, 2*hhoData%face_degree()+1)
!
! -------------- Value of the function at the quadrature point
!
                do idim = 1, hhoCell%ndim
                    call hhoFuncFScalEvalQp(hhoQuadFace, func(idim), nbpara, nompar, &
                                            valpar, hhoCell%ndim, &
                                            FuncValuesFaceQP(idim, 1:MAX_QP_FACE))
                end do
!
! -------------- Compute L2 projection
!
                call hhoL2ProjFaceVec(hhoFace, hhoQuadFace, FuncValuesFaceQP, &
                                      hhoData%face_degree(), rhs_face)
            end if
            call dcopy(fbs, rhs_face, 1, coeff_L2Proj(ind), 1)
            ind = ind+fbs
        end do
!
! --- On cell
!
        call hhoQuadCell%GetQuadCell(hhoCell, 2*hhoData%cell_degree()+1)
!
! -------------- Value of the function at the quadrature point
!
        do idim = 1, hhoCell%ndim
            call hhoFuncFScalEvalQp(hhoQuadCell, func(idim), nbpara, nompar, &
                                    valpar, hhoCell%ndim, FuncValuesCellQP(idim, 1:MAX_QP_CELL))
        end do
!
        call hhoL2ProjCellVec(hhoCell, hhoQuadCell, FuncValuesCellQP, hhoData%cell_degree(), &
                              rhs_cell)
        call dcopy(cbs, rhs_cell, 1, coeff_L2Proj(ind), 1)

!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoL2ProjFieldScal(hhoCell, hhoData, field, coeff_L2Proj)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Data), intent(in)          :: hhoData
        real(kind=8), intent(in)            :: field(*)
        real(kind=8), intent(out)           :: coeff_L2Proj(MSIZE_TDOFS_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the L2-prjoection of a scalar field on a local space (cell + face)
!   In hhoCell      : the current HHO Cell
!   In FuncValuesQP : Values of the function to project at the quadrature points
!   Out coeff_L2Proj: coefficient after projection
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_Face) :: hhoFace
        type(HHO_Quadrature) :: hhoQuadFace, hhoQuadCell
        integer :: cbs, fbs, total_dofs, iFace, ind, ino
        real(kind=8) :: FuncValuesCellQP(MAX_QP_CELL), FuncValuesFaceQP(MAX_QP_FACE)
        real(kind=8) :: FieldValuesNodes(27)
        real(kind=8) :: rhs_face(MSIZE_FACE_SCAL), rhs_cell(MSIZE_CELL_SCAL)

! --------------------------------------------------------------------------------------------------
!
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
        coeff_L2Proj = 0
        FuncValuesCellQP = 0.d0
!
! --- Loop on faces
!
        ind = 1
        do iFace = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iFace)
!
! ----- get quadrature
!
            call hhoQuadFace%GetQuadFace(hhoface, 2*hhoData%face_degree()+1)
!
! -------------- Value of the field at the quadrature point
!
            FieldValuesNodes = 0.d0
            do ino = 1, hhoFace%nbnodes
                FieldValuesNodes(ino) = field(hhoFace%nodes_loc(ino))
            end do
            call hhoFuncRScalEvalQp(hhoFace, hhoQuadFace, FieldValuesNodes, FuncValuesFaceQP)
!
! -------------- Compute L2 projection
!
            call hhoL2ProjFaceScal(hhoFace, hhoQuadFace, FuncValuesFaceQP, &
                                   hhoData%face_degree(), rhs_face)
            call dcopy(fbs, rhs_face, 1, coeff_L2Proj(ind), 1)
            ind = ind+fbs
        end do
!
! --- On cell
!
        call hhoQuadCell%GetQuadCell(hhoCell, 2*hhoData%cell_degree()+1)
!
! -------------- Value of the field at the quadrature point
!
        call hhoFuncRScalEvalCellQp(hhoCell, hhoQuadCell, field, FuncValuesCellQP)
!
        call hhoL2ProjCellScal(hhoCell, hhoQuadCell, FuncValuesCellQP, hhoData%cell_degree(), &
                               rhs_cell)
        call dcopy(cbs, rhs_cell, 1, coeff_L2Proj(ind), 1)
!
    end subroutine
!
end module
