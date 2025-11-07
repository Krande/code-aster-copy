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
    subroutine ace_masse_repartie(nbocc, infdonn, infcarte, grplmax, grpnbma, zjdlm, nbdisc)
!
        use cara_elem_info_type
        use cara_elem_carte_type
!
        integer(kind=8)         :: nbocc
        type (cara_elem_info)   :: infdonn
        type (cara_elem_carte)  :: infcarte(*)
        character(len=24)       :: grplmax(*)
        integer(kind=8)         :: grpnbma(*)
        integer(kind=8)         :: zjdlm(*)
        integer(kind=8)         :: nbdisc
    end subroutine ace_masse_repartie
end interface
