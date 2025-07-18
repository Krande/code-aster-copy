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

subroutine as_mficom(nom, hdfok, medok, cret)
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
#include "asterc/hdfopf.h"
#include "asterc/hdfclf.h"
#include "med/mficom.h"
! person_in_charge: nicolas.sellenet at edf.fr
    aster_int :: cret, hdfok, medok
    hid_t :: fid
    character(len=*) :: nom
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_int :: cret4, hdfok4, medok4
#endif
    cret = 0
    ! On verifie par un appel à HDF que le fichier est bien de type hdf avant de vérifier
    ! la compatibilite afin d'eviter les "Erreur à l'ouverture du fichier" dans MED
    fid = hdfopf(nom)
    if (fid .gt. 0) then
        cret = hdfclf(fid)
        ASSERT(cret .eq. 0)
    else
        cret = -1
        hdfok = 0
        medok = 0
    end if
    if (cret .eq. 0) then
#if !ASTER_MED_SAME_INT_IDT
        call mficom(nom, hdfok4, medok4, cret4)
        cret = to_aster_int(cret4)
        hdfok = to_aster_int(hdfok4)
        medok = to_aster_int(medok4)
#else
        call mficom(nom, hdfok, medok, cret)
#endif
    end if
!
#endif
end subroutine
