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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_read_exte(keywf, i_comp, libr_name, subr_name, nb_vari_umat)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
!
    character(len=16), intent(in) :: keywf
    integer(kind=8), intent(in) :: i_comp
    character(len=255), intent(out) :: libr_name
    character(len=255), intent(out) :: subr_name
    integer(kind=8), intent(out) :: nb_vari_umat
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get parameters for external programs (UMAT)
!
! --------------------------------------------------------------------------------------------------
!
! In  keywf            : factor keyword to read (COMPORTEMENT)
! In  i_comp           : factor keyword index
! Out libr_name        : name of library if UMAT or MFront
! Out subr_name        : name of comportement in library if UMAT or MFront
! Out nb_vari_umat     : number of internal variables for UMAT
!
! --------------------------------------------------------------------------------------------------
!
    libr_name = ' '
    subr_name = ' '
    nb_vari_umat = 0
!
! - Get parameters
!
    ASSERT(i_comp .ne. 0)
    call getvtx(keywf, 'LIBRAIRIE', iocc=i_comp, scal=libr_name)
    call getvtx(keywf, 'NOM_ROUTINE', iocc=i_comp, scal=subr_name)
    call getvis(keywf, 'NB_VARI', iocc=i_comp, scal=nb_vari_umat)
!
end subroutine
