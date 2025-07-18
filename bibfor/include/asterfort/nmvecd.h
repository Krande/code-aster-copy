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
!
interface
    subroutine nmvecd(imate, mate, nmat, matcst, loi,&
                      hook, dt, tp, p, np,&
                      beta, nb, ep, rm, dm,&
                      dsgde, dsgdb, dsgdp, drbde, drpde,&
                      rb, rp, drbdb, drbdp, drpdb,&
                      drpdp, etatf, ier)
        integer(kind=8) :: nb
        integer(kind=8) :: np
        integer(kind=8) :: nmat
        integer(kind=8) :: imate
        real(kind=8) :: mate(nmat, 2)
        character(len=3) :: matcst
        character(len=16) :: loi
        real(kind=8) :: hook(6, 6)
        real(kind=8) :: dt
        real(kind=8) :: tp
        real(kind=8) :: p(np)
        real(kind=8) :: beta(nb)
        real(kind=8) :: ep(*)
        real(kind=8) :: rm
        real(kind=8) :: dm
        real(kind=8) :: dsgde(nb, nb)
        real(kind=8) :: dsgdb(nb, nb)
        real(kind=8) :: dsgdp(nb, np)
        real(kind=8) :: drbde(nb, nb)
        real(kind=8) :: drpde(np, nb)
        real(kind=8) :: rb(nb)
        real(kind=8) :: rp(np)
        real(kind=8) :: drbdb(nb, nb)
        real(kind=8) :: drbdp(nb, np)
        real(kind=8) :: drpdb(np, nb)
        real(kind=8) :: drpdp(np, np)
        character(len=7) :: etatf(3)
        integer(kind=8) :: ier
    end subroutine nmvecd
end interface
