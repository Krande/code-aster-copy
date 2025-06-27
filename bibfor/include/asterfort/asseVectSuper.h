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
    subroutine asseVectSuper(model, mesh, vectElem,&
                             vectScalType, vectElemCoef,&
                             nomacr, meshNbNode,&
                             nec, nbecmx, nbCmp,&
                             icodla, icodge,&
                             idprn1, idprn2,&
                             iapsdl, ianueq, jvale, jresl)
        character(len=8), intent(in) :: model, mesh
        character(len=19), intent(in) :: vectElem
        integer(kind=8), intent(in) :: vectScalType
        real(kind=8), intent(in) :: vectElemCoef
        character(len=8), pointer :: nomacr(:)
        integer(kind=8), intent(in) :: meshNbNode, nbCmp, nec, nbecmx
        integer(kind=8), intent(in) :: iapsdl, ianueq, jvale, jresl
        integer(kind=8), intent(in) :: idprn1, idprn2
        integer(kind=8), intent(inout) :: icodla(nbecmx), icodge(nbecmx)
    end subroutine asseVectSuper
end interface
