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
#include "asterf_types.h"
!
interface
    subroutine meamme(modelz,&
                      matez, matecoz, caraElemz,&
                      time, basez,&
                      matrRigiz,  matrMassz,&
                      matrElemz, &
                      variz, comporz, sddyna)
        character(len=*), intent(in) :: modelz
        character(len=*), intent(in) :: matez, matecoz, caraElemz
        real(kind=8), intent(in) :: time
        character(len=*), intent(in) :: basez
        character(len=*), intent(in) :: matrRigiz, matrMassz, matrElemz
        character(len=*), intent(in) :: variz, comporz
        character(len=19), intent(in) :: sddyna
    end subroutine meamme
end interface
