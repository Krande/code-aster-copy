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
#include "asterf_types.h"
!
interface
    subroutine caeihm(ds_thm, nomte, l_axi, mecani, press1,&
                      press2, tempe, dimdef, dimcon, ndim,&
                      nno1, nno2, npi, npg, dimuel,&
                      iw, ivf1, idf1, ivf2, idf2,&
                      jgano1, iu, ip, ipf, iq,&
                      inte_type)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        character(len=16) :: nomte
        aster_logical :: l_axi, l_steady
        integer(kind=8) :: mecani(8)
        integer(kind=8) :: press1(9)
        integer(kind=8) :: press2(9)
        integer(kind=8) :: tempe(5)
        integer(kind=8) :: dimdef
        integer(kind=8) :: dimcon
        integer(kind=8) :: ndim
        integer(kind=8) :: nno1
        integer(kind=8) :: nno2
        integer(kind=8) :: npi
        integer(kind=8) :: npg
        integer(kind=8) :: dimuel
        integer(kind=8) :: iw
        integer(kind=8) :: ivf1
        integer(kind=8) :: idf1
        integer(kind=8) :: ivf2
        integer(kind=8) :: idf2
        integer(kind=8) :: jgano1
        integer(kind=8) :: iu(3, 18)
        integer(kind=8) :: ip(2, 9)
        integer(kind=8) :: ipf(2, 2, 9)
        integer(kind=8) :: iq(2, 2, 9)
        character(len=3), intent(out) :: inte_type
    end subroutine caeihm
end interface
