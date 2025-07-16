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
    subroutine xcaehm(ds_thm, nomte, l_axi, type_elem, inte_type,&
                      mecani, press1, press2, tempe, dimdef,&
                      dimcon, nmec, np1, np2, ndim,&
                      nno, nnos, nnom, npi, npg,&
                      nddls, nddlm, dimuel, ipoids, ivf,&
                      idfde, ddld, ddlm, ddlp, enrmec, nenr,&
                      dimenr, nnop, nnops, nnopm, enrhyd, ddlc, nfh)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        character(len=16) :: nomte
        character(len=3), intent(out) :: inte_type
        integer(kind=8) :: mecani(5)
        integer(kind=8) :: press1(7)
        integer(kind=8) :: press2(7)
        integer(kind=8) :: tempe(5)
        integer(kind=8) :: dimdef
        integer(kind=8) :: dimcon
        integer(kind=8) :: nmec
        integer(kind=8) :: np1
        integer(kind=8) :: np2
        integer(kind=8) :: nno
        integer(kind=8) :: nnos
        integer(kind=8) :: nnom
        integer(kind=8) :: npi
        integer(kind=8) :: npg
        integer(kind=8) :: nddls
        integer(kind=8) :: nddlm
        integer(kind=8) :: dimuel
        integer(kind=8) :: ipoids
        integer(kind=8) :: ivf
        integer(kind=8) :: idfde
        integer(kind=8) :: ddld
        integer(kind=8) :: ddlm
        integer(kind=8) :: ddlp
        integer(kind=8) :: enrmec(3)
        integer(kind=8) :: nenr
        integer(kind=8) :: dimenr
        integer(kind=8) :: nnop
        integer(kind=8) :: nnops
        integer(kind=8) :: nnopm
        integer(kind=8) :: enrhyd(3)
        integer(kind=8) :: ddlc
        integer(kind=8) :: nfh
        aster_logical, intent(out) :: l_axi
        integer(kind=8), intent(out) :: ndim
        character(len=8), intent(out) :: type_elem(2)
    end subroutine xcaehm
end interface
