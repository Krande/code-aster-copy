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
#include "asterf_types.h"
!
interface
    subroutine ircerl(ifi, nbel, ligrel, nbgrel, longr,&
                      ncmpmx, vale, nomcmp, nomel, loc,&
                      celd, connex, point, nomnos, nbcmpt,&
                      nucmpu, nbnot, numnoe, nbmat, nummai,&
                      lsup, borsup, linf, borinf, lmax,&
                      lmin, lcor, ndim, coor, nolili,&
                      formr, ncmpv, nucmp)
        integer(kind=8) :: ifi
        integer(kind=8) :: nbel
        integer(kind=8) :: ligrel(*)
        integer(kind=8) :: nbgrel
        integer(kind=8) :: longr(*)
        integer(kind=8) :: ncmpmx
        real(kind=8) :: vale(*)
        character(len=*) :: nomcmp(*)
        character(len=*) :: nomel(*)
        character(len=*) :: loc
        integer(kind=8) :: celd(*)
        integer(kind=8) :: connex(*)
        integer(kind=8) :: point(*)
        character(len=*) :: nomnos(*)
        integer(kind=8) :: nbcmpt
        integer(kind=8) :: nucmpu(*)
        integer(kind=8) :: nbnot
        integer(kind=8) :: numnoe(*)
        integer(kind=8) :: nbmat
        integer(kind=8) :: nummai(*)
        aster_logical :: lsup
        real(kind=8) :: borsup
        aster_logical :: linf
        real(kind=8) :: borinf
        aster_logical :: lmax
        aster_logical :: lmin
        aster_logical :: lcor
        integer(kind=8) :: ndim
        real(kind=8) :: coor(*)
        character(len=19) :: nolili
        character(len=*) :: formr
        integer(kind=8) :: ncmpv
        integer(kind=8) :: nucmp(*)
    end subroutine ircerl
end interface
