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
    subroutine xpoajc(nnm, inm, inmtot, nbmac, ise,&
                      npg, jcesd1, jcesd2, jcvid1, jcvid2,&
                      ima, ndim, ndime, iadc, iadv,&
                      jcesv1, jcesl2, jcesv2, jcviv1, jcvil2,&
                      jcviv2)
        integer(kind=8) :: nnm
        integer(kind=8) :: inm
        integer(kind=8) :: inmtot
        integer(kind=8) :: nbmac
        integer(kind=8) :: ise
        integer(kind=8) :: npg
        integer(kind=8) :: jcesd1
        integer(kind=8) :: jcesd2
        integer(kind=8) :: jcvid1
        integer(kind=8) :: jcvid2
        integer(kind=8) :: ima
        integer(kind=8) :: ndim
        integer(kind=8) :: ndime
        integer(kind=8) :: iadc
        integer(kind=8) :: iadv
        integer(kind=8) :: jcesv1
        integer(kind=8) :: jcesl2
        integer(kind=8) :: jcesv2
        integer(kind=8) :: jcviv1
        integer(kind=8) :: jcvil2
        integer(kind=8) :: jcviv2
    end subroutine xpoajc
end interface 
