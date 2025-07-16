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

subroutine as_msdcrw(fid, lmname, jname, numdt, numit, entlcl, &
                     geolcl, entdst, geodst, ncorr, corrtab, cret)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
    implicit none
#include "asterf_config.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/conv_int.h"
#include "asterfort/utmess.h"
#include "med/msdcrw.h"
    character(len=*) :: lmname, jname
    med_idt :: fid
    aster_int :: numdt, numit, entlcl, geolcl, entdst
    aster_int :: geodst, ncorr, corrtab(*), cret

#ifndef ASTER_HAVE_MED
    call utmess('F', 'FERMETUR_2')
#else

#ifdef ASTER_DEBUG_MED
    write (6, *) '=== as_msdcrw fid=', fid
    write (6, *) '=== as_msdcrw lmname=', lmname
    write (6, *) '=== as_msdcrw jname=', jname
    write (6, *) '=== as_msdcrw numdt=', numdt
    write (6, *) '=== as_msdcrw numit=', numit
    write (6, *) '=== as_msdcrw entlcl=', entlcl
    write (6, *) '=== as_msdcrw geolcl=', geolcl
    write (6, *) '=== as_msdcrw entdst=', entdst
    write (6, *) '=== as_msdcrw geodst=', geodst
    write (6, *) '=== as_msdcrw ncorr=', ncorr
    write (6, *) '=== as_msdcrw corrtab(1:min(3,ncorr))=', corrtab(1:min(3, ncorr))
#endif

#if !ASTER_MED_SAME_INT_IDT
    med_idt :: fid4
    med_int :: numdt4, numit4, entlcl4, geolcl4, entdst4, geodst4, cret4, ncorr4
    med_int, allocatable :: corrtab4(:)
    fid4 = to_med_idt(fid)
    numdt4 = to_med_int(numdt)
    numit4 = to_med_int(numit)
    entlcl4 = to_med_int(entlcl)
    geolcl4 = to_med_int(geolcl)
    entdst4 = to_med_int(entdst)
    geodst4 = to_med_int(geodst)
    ncorr4 = to_med_int(ncorr)
    ASSERT(ncorr .gt. 0)
    allocate (corrtab4(2*ncorr))

    call conv_int('ast->med', 2*ncorr, vi_ast=corrtab, vi_med=corrtab4)
    call msdcrw(fid4, lmname, jname, numdt4, numit4, entlcl4, &
                geolcl4, entdst4, geodst4, ncorr4, corrtab4, cret4)

    deallocate (corrtab4)
    cret = to_aster_int(cret4)

#else
    call msdcrw(fid, lmname, jname, numdt, numit, entlcl, &
                geolcl, entdst, geodst, ncorr, corrtab, cret)
#endif

#ifdef ASTER_DEBUG_MED
    write (6, *) '=== as_msdcrw cret=', cret
#endif

#endif
end subroutine
