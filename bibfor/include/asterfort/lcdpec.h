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
interface
    subroutine lcdpec(BEHinteg, &
                      vind, nbcomm, nmat, ndt, cpmono,&
                      materf, iter, nvi, itmax, toler,&
                      pgl, nfs, nsg, toutms, hsr,&
                      dt, dy, yd, vinf,&
                      sigf, df, nr, mod,&
                      codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        integer(kind=8) :: nsg
        integer(kind=8) :: nfs
        integer(kind=8) :: nmat
        real(kind=8) :: vind(*)
        integer(kind=8) :: nbcomm(nmat, 3)
        integer(kind=8) :: ndt
        character(len=24) :: cpmono(5*nmat+1)
        real(kind=8) :: materf(nmat*2)
        integer(kind=8) :: iter
        integer(kind=8) :: nvi
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        real(kind=8) :: pgl(3, 3)
        real(kind=8) :: toutms(nfs, nsg, 6)
        real(kind=8) :: hsr(nsg, nsg)
        real(kind=8) :: dt
        real(kind=8) :: dy(*)
        real(kind=8) :: yd(*)
        real(kind=8) :: vinf(*)
        real(kind=8) :: sigf(6)
        real(kind=8) :: df(3, 3)
        integer(kind=8) :: nr
        character(len=8) :: mod
        integer(kind=8) :: codret
    end subroutine lcdpec
end interface
