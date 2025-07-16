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
    subroutine fnodil(option, typmod, ds_dil, ndim, nnos, &
                    nnom, npg, nddl, dimdef, iw, vff, &
                    vffb, idff,idffb,geomi, compor, &
                    sief,  fint)
        use dil_type
        character(len=8),intent(in) :: typmod(*)
        character(len=16),intent(in):: option, compor(*)
        type(dil_modelisation)      :: ds_dil
        integer(kind=8),intent(in)          :: ndim,nnos,nnom,npg,nddl,dimdef
        integer(kind=8),intent(in)          :: iw,idff,idffb
        real(kind=8),intent(in)     :: geomi(ndim,nnos+nnom)
        real(kind=8),intent(in)     :: vff(nnos+nnom, npg),vffb(nnos, npg)
        real(kind=8),intent(in)     :: sief(dimdef*npg)
        real(kind=8),intent(out)    :: fint(nddl)
    end subroutine fnodil
end interface
