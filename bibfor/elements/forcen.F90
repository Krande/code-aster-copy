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
subroutine forcen(rnormc, intsn, nb1, xi, xr, &
                  rho, epais, vomega, vecl1, xa)
    implicit none
!
    integer(kind=8) :: intsn, nb1, intsx, ie(3, 3, 3)
    real(kind=8) :: wgt, rho
    real(kind=8) :: xi(3, *), xr(*), vomega(3), vecl1(42), xa(3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i2, ib, ik, ip, iq
    integer(kind=8) :: ir, j, jb, k, l1
    real(kind=8) :: epais, rnormc
!-----------------------------------------------------------------------
    wgt = xr(127-1+intsn)
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                ie(i, j, k) = 0
            end do
        end do
    end do
    ie(1, 2, 3) = 1
    ie(1, 3, 2) = -1
    ie(2, 1, 3) = -1
    ie(2, 3, 1) = 1
    ie(3, 1, 2) = 1
    ie(3, 2, 1) = -1
!
    l1 = 135
    intsx = 8*(intsn-1)
!
    i1 = l1+intsx
!
    do ib = 1, nb1
        i2 = 5*(ib-1)
        do iq = 1, 3
            do ip = 1, 3
                do ik = 1, 3
                    do ir = 1, 3
                        do jb = 1, nb1
                            vecl1(i2+1) = vecl1(i2+1)-wgt*rho*epais* &
                                          rnormc*ie(iq, ip, ik)*ie(1, ir, iq)*vomega( &
                                          ir)*vomega(ip)*(xi(ik, jb)-xa(ik))*xr(i1+ &
                                                                                ib)*xr(i1+jb)
!
                            vecl1(i2+2) = vecl1(i2+2)-wgt*rho*epais* &
                                          rnormc*ie(iq, ip, ik)*ie(2, ir, iq)*vomega( &
                                          ir)*vomega(ip)*(xi(ik, jb)-xa(ik))*xr(i1+ &
                                                                                ib)*xr(i1+jb)
!
                            vecl1(i2+3) = vecl1(i2+3)-wgt*rho*epais* &
                                          rnormc*ie(iq, ip, ik)*ie(3, ir, iq)*vomega( &
                                          ir)*vomega(ip)*(xi(ik, jb)-xa(ik))*xr(i1+ &
                                                                                ib)*xr(i1+jb)
                        end do
                    end do
                end do
            end do
        end do
    end do
end subroutine
