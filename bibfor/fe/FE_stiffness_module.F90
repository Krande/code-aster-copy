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
!
module FE_stiffness_module
!
    use FE_basis_module
    use FE_quadrature_module
    use HHO_utils_module, only: hhoCopySymPartMat
!
    implicit none
!
    private
#include "asterf_types.h"
#include "FE_module.h"
#include "blas/dgemv.h"
#include "asterfort/assert.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE - Finite Element
!
! Module to compute stiffness terms
!
! --------------------------------------------------------------------------------------------------
!
    public :: FEStiffResiScal, FEStiffJacoScal, FEStiffResiScalAdd, FEStiffJacoScalAdd
    public :: FEStiffResiVectSymAdd, FEStiffJacoVectSymAdd, FEStiffGeomVectSymAdd
!    private  ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffResiScalAdd(FEBasis, BGSEval, weight, ValuesQP, vec)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(3, MAX_BS)  :: BGSEval
        real(kind=8), intent(in)            :: weight
        real(kind=8), intent(inout)         :: vec(MAX_BS)
        real(kind=8), intent(in)            :: ValuesQP(3)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the rigidity vector
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
!
!
        call dgemv('T', FEBasis%ndim, FEBasis%size, weight, BGSEval, 3, &
                   ValuesQP, 1, 1.d0, vec, 1)
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffResiScal(FEQuad, FEBasis, ValuesQP, vec)
!
        implicit none
!
        type(FE_Quadrature), intent(in)     :: FEQuad
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(out)           :: vec(MAX_BS)
        real(kind=8), intent(in)            :: ValuesQP(3, MAX_QP)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity vector
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer :: ipg
        real(kind=8), dimension(3, MAX_BS) :: BSEval
!
        vec = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%grad(FEQuad%points_param(1:3, ipg), FEQuad%jacob(1:3, 1:3, ipg))
!
            call FEStiffResiScalAdd(FEBasis, BSEval, FEQuad%weights(ipg), ValuesQP(1:3, ipg), vec)
        end do
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffJacoScalAdd(FEBasis, BGSEval, weight, ValueQP, mat)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(3, MAX_BS) :: BGSEval
        real(kind=8), intent(in)            :: weight
        real(kind=8), intent(inout)           :: mat(MAX_BS, MAX_BS)
        real(kind=8), intent(in)            :: ValueQP(3, 3)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity matrix
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer :: j
        real(kind=8) :: Kgradj(3)
!
        do j = 1, FEBasis%size
            if (FEBasis%ndim == 3) then
                Kgradj(1) = ValueQP(1, 1)*BGSEval(1, j)+ValueQP(1, 2)*BGSEval(2, j)+ &
                            ValueQP(1, 3)*BGSEval(3, j)
                Kgradj(2) = ValueQP(2, 1)*BGSEval(1, j)+ValueQP(2, 2)*BGSEval(2, j)+ &
                            ValueQP(2, 3)*BGSEval(3, j)
                Kgradj(3) = ValueQP(3, 1)*BGSEval(1, j)+ValueQP(3, 2)*BGSEval(2, j)+ &
                            ValueQP(3, 3)*BGSEval(3, j)
            elseif (FEBasis%ndim == 2) then
                Kgradj(1) = ValueQP(1, 1)*BGSEval(1, j)+ValueQP(1, 2)*BGSEval(2, j)
                Kgradj(2) = ValueQP(2, 1)*BGSEval(1, j)+ValueQP(2, 2)*BGSEval(2, j)
            else
                Kgradj(1) = ValueQP(1, 1)*BGSEval(1, j)
            end if
            call dgemv('T', FEBasis%ndim, j, weight, BGSEval, 3, Kgradj, 1, 1.d0, mat(:, j), 1)
        end do
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffJacoScal(FEQuad, FEBasis, ValuesQP, mat)
!
        implicit none
!
        type(FE_Quadrature), intent(in)     :: FEQuad
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(out)           :: mat(MAX_BS, MAX_BS)
        real(kind=8), intent(in)            :: ValuesQP(3, 3, MAX_QP)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity matrix
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer :: ipg
        real(kind=8), dimension(3, MAX_BS) :: BSEval
!
        mat = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%grad(FEQuad%points_param(1:3, ipg), FEQuad%jacob(1:3, 1:3, ipg))

            call FEStiffJacoScalAdd(FEBasis, BSEval, FEQuad%weights(ipg), ValuesQP(1:3, 1:3, ipg), &
                                    mat)
!
        end do
!
! ----- Copy the lower part
!
        call hhoCopySymPartMat('U', mat(1:FEBasis%size, 1:FEBasis%size))
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffResiVectSymAdd(FEBasis, def, weight, stress, vec)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(6, MAX_BS, 3) :: def
        real(kind=8), intent(in)            :: weight
        real(kind=8), intent(inout)         :: vec(*)
        real(kind=8), intent(in)            :: stress(6)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the rigidity vector (with symetric stess)
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i, ind
!
        select case (FEBasis%ndim)
        case (2)
            ind = 0
            do i = 1, FEBasis%size
                vec(ind+1) = vec(ind+1)+weight*(def(1, i, 1)*stress(1)+def(2, i, 1)*stress(2)+ &
                                                def(3, i, 1)*stress(3)+def(4, i, 1)*stress(4))
                vec(ind+2) = vec(ind+2)+weight*(def(1, i, 2)*stress(1)+def(2, i, 2)*stress(2)+ &
                                                def(3, i, 2)*stress(3)+def(4, i, 2)*stress(4))
                ind = ind+2
            end do
        case (3)
            ind = 0
            do i = 1, FEBasis%size
                vec(ind+1) = vec(ind+1)+weight*(def(1, i, 1)*stress(1)+def(2, i, 1)*stress(2)+ &
                                                def(3, i, 1)*stress(3)+def(4, i, 1)*stress(4)+ &
                                                def(5, i, 1)*stress(5)+def(6, i, 1)*stress(6))
                vec(ind+2) = vec(ind+2)+weight*(def(1, i, 2)*stress(1)+def(2, i, 2)*stress(2)+ &
                                                def(3, i, 2)*stress(3)+def(4, i, 2)*stress(4)+ &
                                                def(5, i, 2)*stress(5)+def(6, i, 2)*stress(6))
                vec(ind+3) = vec(ind+3)+weight*(def(1, i, 3)*stress(1)+def(2, i, 3)*stress(2)+ &
                                                def(3, i, 3)*stress(3)+def(4, i, 3)*stress(4)+ &
                                                def(5, i, 3)*stress(5)+def(6, i, 3)*stress(6))
                ind = ind+3
            end do
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffJacoVectSymAdd(FEBasis, def, weight, dsidep, l_matsym, mat)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(6, MAX_BS, 3) :: def
        real(kind=8), intent(in)            :: weight
        aster_logical, intent(in)           :: l_matsym
        real(kind=8), intent(inout)         :: mat(*)
        real(kind=8), intent(in)            :: dsidep(6, 6)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the rigidity vector (with symetric stess)
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_n, i_d, i_tens, j_d, j_n, j1
        integer :: kk, kkd
        real(kind=8) :: tmp, sig(6)
!
        select case (FEBasis%ndim)
        case (2)
            if (l_matsym) then
                do i_n = 1, FEBasis%size
                    do i_d = 1, 2
                        kkd = (2*(i_n-1)+i_d-1)*(2*(i_n-1)+i_d)/2
                        do i_tens = 1, 4
                            sig(i_tens) = 0.d0
                            sig(i_tens) = sig(i_tens)+def(1, i_n, i_d)*dsidep(1, i_tens)
                            sig(i_tens) = sig(i_tens)+def(2, i_n, i_d)*dsidep(2, i_tens)
                            sig(i_tens) = sig(i_tens)+def(3, i_n, i_d)*dsidep(3, i_tens)
                            sig(i_tens) = sig(i_tens)+def(4, i_n, i_d)*dsidep(4, i_tens)
                        end do
                        do j_d = 1, 2
                            do j_n = 1, i_n
                                if (j_n .eq. i_n) then
                                    j1 = i_d
                                else
                                    j1 = 2
                                end if
                                if (j_d .le. j1) then
                                    tmp = def(1, j_n, j_d)*sig(1)+def(2, j_n, j_d)*sig(2)+ &
                                          def(3, j_n, j_d)*sig(3)+def(4, j_n, j_d)*sig(4)
                                    kk = kkd+2*(j_n-1)+j_d
                                    mat(kk) = mat(kk)+tmp*weight
                                end if
                            end do
                        end do
                    end do
                end do
            else
                do i_n = 1, FEBasis%size
                    do i_d = 1, 2
                        do i_tens = 1, 4
                            sig(i_tens) = 0.d0
                            sig(i_tens) = sig(i_tens)+def(1, i_n, i_d)*dsidep(1, i_tens)
                            sig(i_tens) = sig(i_tens)+def(2, i_n, i_d)*dsidep(2, i_tens)
                            sig(i_tens) = sig(i_tens)+def(3, i_n, i_d)*dsidep(3, i_tens)
                            sig(i_tens) = sig(i_tens)+def(4, i_n, i_d)*dsidep(4, i_tens)
                        end do
                        do j_d = 1, 2
                            do j_n = 1, FEBasis%size
                                tmp = def(1, j_n, j_d)*sig(1)+def(2, j_n, j_d)*sig(2)+ &
                                      def(3, j_n, j_d)*sig(3)+def(4, j_n, j_d)*sig(4)
                                kk = 2*FEBasis%size*(2*(i_n-1)+i_d-1)+2*(j_n-1)+j_d
                                mat(kk) = mat(kk)+tmp*weight
                            end do
                        end do
                    end do
                end do
            end if
        case (3)
            if (l_matsym) then
                do i_n = 1, FEBasis%size
                    do i_d = 1, 3
                        kkd = (3*(i_n-1)+i_d-1)*(3*(i_n-1)+i_d)/2
                        do i_tens = 1, 6
                            sig(i_tens) = 0.d0
                            sig(i_tens) = sig(i_tens)+def(1, i_n, i_d)*dsidep(1, i_tens)
                            sig(i_tens) = sig(i_tens)+def(2, i_n, i_d)*dsidep(2, i_tens)
                            sig(i_tens) = sig(i_tens)+def(3, i_n, i_d)*dsidep(3, i_tens)
                            sig(i_tens) = sig(i_tens)+def(4, i_n, i_d)*dsidep(4, i_tens)
                            sig(i_tens) = sig(i_tens)+def(5, i_n, i_d)*dsidep(5, i_tens)
                            sig(i_tens) = sig(i_tens)+def(6, i_n, i_d)*dsidep(6, i_tens)
                        end do
                        do j_d = 1, 3
                            do j_n = 1, i_n
                                if (j_n .eq. i_n) then
                                    j1 = i_d
                                else
                                    j1 = 3
                                end if
                                if (j_d .le. j1) then
                                    tmp = def(1, j_n, j_d)*sig(1)+def(2, j_n, j_d)*sig(2)+ &
                                          def(3, j_n, j_d)*sig(3)+def(4, j_n, j_d)*sig(4)+ &
                                          def(5, j_n, j_d)*sig(5)+def(6, j_n, j_d)*sig(6)
                                    kk = kkd+3*(j_n-1)+j_d
                                    mat(kk) = mat(kk)+tmp*weight
                                end if
                            end do
                        end do
                    end do
                end do
            else
                do i_n = 1, FEBasis%size
                    do i_d = 1, 3
                        do i_tens = 1, 6
                            sig(i_tens) = 0.d0
                            sig(i_tens) = sig(i_tens)+def(1, i_n, i_d)*dsidep(1, i_tens)
                            sig(i_tens) = sig(i_tens)+def(2, i_n, i_d)*dsidep(2, i_tens)
                            sig(i_tens) = sig(i_tens)+def(3, i_n, i_d)*dsidep(3, i_tens)
                            sig(i_tens) = sig(i_tens)+def(4, i_n, i_d)*dsidep(4, i_tens)
                            sig(i_tens) = sig(i_tens)+def(5, i_n, i_d)*dsidep(5, i_tens)
                            sig(i_tens) = sig(i_tens)+def(6, i_n, i_d)*dsidep(6, i_tens)
                        end do
                        do j_d = 1, 3
                            do j_n = 1, FEBasis%size
                                tmp = def(1, j_n, j_d)*sig(1)+def(2, j_n, j_d)*sig(2)+ &
                                      def(3, j_n, j_d)*sig(3)+def(4, j_n, j_d)*sig(4)+ &
                                      def(5, j_n, j_d)*sig(5)+def(6, j_n, j_d)*sig(6)
                                kk = 3*FEBasis%size*(3*(i_n-1)+i_d-1)+3*(j_n-1)+j_d
                                mat(kk) = mat(kk)+tmp*weight
                            end do
                        end do
                    end do
                end do
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffGeomVectSymAdd(FEBasis, pff, weight, stress, l_matsym, mat)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(6, MAX_BS, MAX_BS) :: pff
        real(kind=8), intent(in)            :: weight
        aster_logical, intent(in)           :: l_matsym
        real(kind=8), intent(inout)         :: mat(*)
        real(kind=8), intent(in)            :: stress(6)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the rigidity vector (with symetric stess)
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
        integer :: i_n, i_d, j_n, j1
        integer :: kk, kkd
        real(kind=8) :: tmp
!
        select case (FEBasis%ndim)
        case (2)
            if (l_matsym) then
                do i_n = 1, FEBasis%size
                    do i_d = 1, 2
                        kkd = (2*(i_n-1)+i_d-1)*(2*(i_n-1)+i_d)/2
                        do j_n = 1, i_n
                            if (j_n .eq. i_n) then
                                j1 = i_d
                            else
                                j1 = 2
                            end if
                            tmp = pff(1, i_n, j_n)*stress(1)+pff(2, i_n, j_n)*stress(2)+ &
                                  pff(3, i_n, j_n)*stress(3)+pff(4, i_n, j_n)*stress(4)
                            if (i_d .le. j1) then
                                kk = kkd+2*(j_n-1)+i_d
                                mat(kk) = mat(kk)+tmp*weight
                            end if
                        end do
                    end do
                end do
            else
                do i_n = 1, FEBasis%size
                    do i_d = 1, 2
                        do j_n = 1, FEBasis%size
                            tmp = pff(1, i_n, j_n)*stress(1)+pff(2, i_n, j_n)*stress(2)+ &
                                  pff(3, i_n, j_n)*stress(3)+pff(4, i_n, j_n)*stress(4)
                            kk = 2*FEBasis%size*(2*(i_n-1)+i_d-1)+2*(j_n-1)+i_d
                            mat(kk) = mat(kk)+tmp*weight
                        end do
                    end do
                end do
            end if
        case (3)
            if (l_matsym) then
                do i_n = 1, FEBasis%size
                    do i_d = 1, 3
                        kkd = (3*(i_n-1)+i_d-1)*(3*(i_n-1)+i_d)/2
                        do j_n = 1, i_n
                            if (j_n .eq. i_n) then
                                j1 = i_d
                            else
                                j1 = 3
                            end if
                            if (i_d .le. j1) then
                                tmp = pff(1, i_n, j_n)*stress(1)+pff(2, i_n, j_n)*stress(2)+ &
                                      pff(3, i_n, j_n)*stress(3)+pff(4, i_n, j_n)*stress(4)+ &
                                      pff(5, i_n, j_n)*stress(5)+pff(6, i_n, j_n)*stress(6)
                                kk = kkd+3*(j_n-1)+i_d
                                mat(kk) = mat(kk)+tmp*weight
                            end if
                        end do
                    end do
                end do
            else
                do i_n = 1, FEBasis%size
                    do i_d = 1, 3
                        do j_n = 1, FEBasis%size
                            tmp = pff(1, i_n, j_n)*stress(1)+pff(2, i_n, j_n)*stress(2)+ &
                                  pff(3, i_n, j_n)*stress(3)+pff(4, i_n, j_n)*stress(4)+ &
                                  pff(5, i_n, j_n)*stress(5)+pff(6, i_n, j_n)*stress(6)
                            kk = 3*FEBasis%size*(3*(i_n-1)+i_d-1)+3*(j_n-1)+i_d
                            mat(kk) = mat(kk)+tmp*weight
                        end do
                    end do
                end do
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
end module
