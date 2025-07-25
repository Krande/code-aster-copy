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

subroutine as_mmhmii(fid, indice, maa, dim, type, &
                     desc, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mmhmii.h"
    med_idt :: fid
    aster_int :: dim, cret, indice, type
    character(len=64) :: maa
    character(len=200) :: desc
    character(len=16) :: descdt
    character(len=16) :: nom(3), unit(3)
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: dim4, cret4, indic4, type4, dimb4, typtr4
    med_int :: nbseq4, typre4
    fidm = to_med_idt(fid)
    indic4 = to_med_int(indice)
!   -- initialisation des chaines (prudent car appel C ensuite):
    maa = ' '
    desc = ' '
    descdt = ' '
    nom(:) = ' '
    unit(:) = ' '
    call mmhmii(fidm, indic4, maa, dim4, dimb4, &
                type4, desc, descdt, typtr4, nbseq4, &
                typre4, nom, unit, cret4)
    dim = dim4
    type = type4
    cret = cret4
#else
    aster_int :: dimb, typtri, nbseq, typrep
    call mmhmii(fid, indice, maa, dim, dimb, &
                type, desc, descdt, typtri, nbseq, &
                typrep, nom, unit, cret)
#endif
!
#endif
end subroutine
