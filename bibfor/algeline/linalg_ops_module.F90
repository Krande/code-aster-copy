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

module linalg_ops_module

    implicit none

    interface as_matmul
        module procedure &
            f_as_matmul_rank22, f_as_matmul_rank21, &
            f_as_matmul_rank12
    end interface as_matmul

contains

    function f_as_matmul_rank22(A, B) result(C)
        implicit none
        real(kind=8), dimension(:, :), intent(in) :: A, B
        real(kind=8), dimension(size(A, 1), size(B, 2)) :: C
        call r_as_matmul(A, B, size(A, 1), size(A, 2), size(B, 2), C)
    end function f_as_matmul_rank22

    function f_as_matmul_rank21(A, B) result(C)
        implicit none
        real(kind=8), dimension(:, :), intent(in) :: A
        real(kind=8), dimension(size(A, 2)), intent(in) :: B
        real(kind=8), dimension(size(A, 1)) :: C
        call r_as_matmul(A, B, size(A, 1), size(A, 2), 1, C)
    end function f_as_matmul_rank21

    function f_as_matmul_rank12(A, B) result(C)
        implicit none
        real(kind=8), dimension(:, :), intent(in) :: B
        real(kind=8), dimension(size(B, 1)), intent(in) :: A
        real(kind=8), dimension(size(B, 2)) :: C
        call r_as_matmul(A, B, 1, size(A), size(B, 2), C)
    end function f_as_matmul_rank12

    subroutine r_as_matmul(A, B, n1, n2, n3, C)
        implicit none
        integer(kind=8) :: n1, n2, n3, i, j, k
        real(kind=8) :: A(n1, n2), B(n2, n3), C(n1, n3)

        C(:, :) = 0.d0

        do k = 1, n3
            do j = 1, n2
                do i = 1, n1
                    C(i, k) = A(i, j)*B(j, k)+C(i, k)
                end do
            end do
        end do

    end subroutine r_as_matmul

end module linalg_ops_module
