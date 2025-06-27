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
    subroutine pipepe(BEHinteg,&
                      pilo, ndim, nno, npg, ipoids,&
                      ivf, idfde, geom, typmod, mate,&
                      compor, lgpg, deplm, sigm, vim,&
                      ddepl, depl0, depl1, copilo,&
                      iborne, ictau)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: ndim
        character(len=16) :: pilo
        integer(kind=8) :: nno
        integer(kind=8) :: ipoids
        integer(kind=8) :: ivf
        integer(kind=8) :: idfde
        real(kind=8) :: geom(ndim, *)
        character(len=8) :: typmod(*)
        integer(kind=8) :: mate
        character(len=16) :: compor(*)
        real(kind=8) :: deplm(*)
        real(kind=8) :: sigm(2*ndim, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: ddepl(*)
        real(kind=8) :: depl0(*)
        real(kind=8) :: depl1(*)
        real(kind=8) :: copilo(5, npg)
        integer(kind=8) :: iborne
        integer(kind=8) :: ictau
    end subroutine pipepe
end interface
