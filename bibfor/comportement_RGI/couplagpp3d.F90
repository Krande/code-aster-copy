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

subroutine couplagpp3d(a, ngf, na, nc, nf0, ig, &
                       dgfa_ds, dpfa_ds, dpfa_dpg, dpg_depspg6, &
                       raideur66p, dpfa_dr, dra_dl, dpg_depsa6)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
!   autres couplages : fluage -> plasticite

!   ************************************************************************
    implicit none

    integer(kind=8), intent(in) :: ngf, na, nc, nf0, ig(nc)
    real(kind=8), intent(inout) :: a(ngf, ngf+1)
    real(kind=8), intent(in) :: dgfa_ds(nc, 6), dpfa_ds(nc, 6), dpfa_dpg(nc)
    real(kind=8), intent(in) :: dpg_depspg6(6), raideur66p(6, 6), dpfa_dr(nc)
    real(kind=8), intent(in) :: dra_dl(nc), dpg_depsa6(6)
!
    integer(kind=8) :: i, j, k, l
    real(kind=8) :: som1, som2
!
    do i = 1, na
        do j = 1, na
            som2 = 0.d0
            do k = 1, 6
!               influence des autres deformations plastiques
                som1 = 0.d0
                do l = 1, 6
                    som1 = som1-raideur66p(k, l)*dgfa_ds(j, l)
                end do
                som2 = som2+dpfa_ds(i, k)*som1
!               influence de la pression intra poreuse
                if ((ig(j) .ge. 7) .and. (ig(j) .le. 9)) then
                    som2 = som2+dpfa_dpg(i)*dpg_depspg6(k)*dgfa_ds(j, k)
                else
                    som2 = som2+dpfa_dpg(i)*dpg_depsa6(k)*dgfa_ds(j, k)
                end if
            end do
            a(nf0+i, nf0+j) = som2
            if (i .eq. j) then
!              prise en compte de l ecrouissage (dra_dl derive de
!              la resistance par rapport au multiplicateur plastique)
                a(nf0+i, nf0+j) = a(nf0+i, nf0+j)+dpfa_dr(i)*dra_dl(i)
            end if
        end do
    end do

end subroutine
