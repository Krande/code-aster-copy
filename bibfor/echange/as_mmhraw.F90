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

subroutine as_mmhraw(fid, nomail, typgeo, nomatt, nbrval, &
                     tabval, codret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/utmess.h"
#include "med/mmhraw.h"
    character(len=*) :: nomail, nomatt
    med_idt :: fid
    aster_int :: typgeo, nbrval, codret
    real(kind=8) :: tabval(*)
    aster_int :: numdt, numit
    parameter(numdt=-1)
    parameter(numit=-1)
#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else
!
#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fidm
    med_int :: typge4, nbrva4, codre4, numdt4, numit4
    fidm = to_med_idt(fid)
    typge4 = typgeo
    nbrva4 = nbrval
    numdt4 = numdt
    numit4 = numit
    call mmhraw(fidm, nomail, numdt4, numit4, typge4, &
                nomatt, nbrva4, tabval, codre4)
    codret = codre4
#else
    call mmhraw(fid, nomail, numdt, numit, typgeo, &
                nomatt, nbrval, tabval, codret)
#endif
!
#endif
end subroutine
