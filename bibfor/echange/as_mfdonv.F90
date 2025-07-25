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

subroutine as_mfdonv(fid, cha, typent, typgeo, noma, &
                     numdt, numo, pit, nompro, stm, &
                     npr, nomloc, nip, n, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mfdonv.h"
    med_idt :: fid
    aster_int :: typent, typgeo, stm, npr, nip, n, cret, numdt, numo, pit
    character(len=*) :: cha, nompro, nomloc, noma
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: typen4, typge4, stm4, npr4, nip4, n4, cret4
    med_int :: numdt4, numo4, pit4
    fidm = to_med_idt(fid)
    typen4 = typent
    typge4 = typgeo
    numdt4 = numdt
    numo4 = numo
    stm4 = stm
    pit4 = pit
    call mfdonv(fidm, cha, numdt4, numo4, typen4, &
                typge4, noma, pit4, stm4, nompro, &
                npr4, nomloc, nip4, n4, cret4)
    npr = npr4
    nip = nip4
    n = n4
    cret = cret4
#else
    call mfdonv(fid, cha, numdt, numo, typent, &
                typgeo, noma, pit, stm, nompro, &
                npr, nomloc, nip, n, cret)
#endif
!
#endif
end subroutine
