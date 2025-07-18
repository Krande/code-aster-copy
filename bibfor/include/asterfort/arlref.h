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
    subroutine arlref(elrefe, fami, nomte, ndim, nno,&
                      nnos, npg, jpoids, jcoopg, jvf,&
                      jdfde, jdfd2, jgano)
        character(len=*), intent(in), optional :: elrefe
        character(len=*), intent(in)    :: fami
        character(len=16), intent(in)   :: nomte
        integer(kind=8), intent(out), optional  :: ndim
        integer(kind=8), intent(out), optional  :: nno
        integer(kind=8), intent(out), optional  :: nnos
        integer(kind=8), intent(out), optional  :: npg
        integer(kind=8), intent(out), optional  :: jpoids
        integer(kind=8), intent(out), optional  :: jcoopg
        integer(kind=8), intent(out), optional  :: jvf
        integer(kind=8), intent(out), optional  :: jdfde
        integer(kind=8), intent(out), optional  :: jdfd2
        integer(kind=8), intent(out), optional  :: jgano
    end subroutine arlref
end interface
