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

subroutine utmess(typ, idmess, nk, valk, sk, &
                  ni, vali, si, nr, valr, &
                  sr, num_except, fname)
    use calcul_module, only: calcul_status
! person_in_charge: mathieu.courtois at edf.fr
!
! All messages (informations, warnings, errors) should be printed through this subroutine.
! Only the first two arguments are compulsory.
! Example: call utmess('A', 'SUPERVIS_22')
!
! To pass a single value, just use sk/si/sr=value.
! To pass more values, use nk/ni/nr=<number of values> + valk/vali/valr=<array of values>.
! Example: call utmess('A', 'MECANONLINE_34', nr=2, valr=[a, b])
!
! If 'fname' is provided, it must be valid filename and the message will be written
! into this file instead of standard MESSAGE, RESULTAT, ERREUR files.
!
! See comments in utmess_core for details
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/temess.h"
#include "asterfort/utmess_core.h"
!
    character(len=*), intent(in) :: typ
    character(len=*), intent(in) :: idmess
    integer(kind=8), intent(in), optional :: nk
    character(len=*), intent(in), optional, target :: valk(*)
    character(len=*), intent(in), optional :: sk
    integer(kind=8), intent(in), optional :: ni
    integer(kind=8), intent(in), optional, target :: vali(*)
    integer(kind=8), intent(in), optional :: si
    integer(kind=8), intent(in), optional :: nr
    real(kind=8), intent(in), optional, target :: valr(*)
    real(kind=8), intent(in), optional :: sr
    integer(kind=8), intent(in), optional :: num_except
    character(len=*), optional :: fname
!
!   working variables
    integer(kind=8) :: unk, uni, unr, nexcep
    character(len=256), target :: uvk(1)
!   because it is not supported by older versions of gfortran, we use two different
!   calls to utmess_core
!    character(len=:), pointer :: ptrk(:)
    aster_logical :: use_valk, under_te0000
    integer(kind=8), target :: uvi(1)
    integer(kind=8), pointer :: ptri(:) => null()
    real(kind=8), target :: uvr(1)
    real(kind=8), pointer :: ptrr(:) => null()
    character(len=2) :: typ2
    character(len=256) :: ufname
!-----------------------------------------------------------------------------------

!
    ASSERT(ENSEMBLE2(nk, valk))
    ASSERT(ENSEMBLE2(ni, vali))
    ASSERT(ENSEMBLE2(nr, valr))
    ASSERT(EXCLUS2(valk, sk))
    ASSERT(EXCLUS2(vali, si))
    ASSERT(EXCLUS2(valr, sr))
    ASSERT(absent(num_except) .or. typ == 'Z')
    nexcep = 0
    if (present(num_except)) then
        nexcep = num_except
    end if
!   associate pointers to valk or sk
    unk = 1
    uvk(1) = ' '
    use_valk = .false.
    if (AU_MOINS_UN2(sk, valk)) then
        if (present(nk)) then
            unk = nk
            use_valk = .true.
        else
            unk = 1
            uvk(1) = sk
        end if
    end if
!   associate pointers to vali or si
    uni = 1
    uvi(1) = 0
    ptri => uvi(1:1)
    if (AU_MOINS_UN2(si, vali)) then
        if (present(ni)) then
            uni = ni
            ptri => vali(1:ni)
        else
            uni = 1
            uvi(1) = si
            ptri => uvi(1:1)
        end if
    end if
!   associate pointers to valr or sr
    unr = 1
    uvr(1) = 0.d0
    ptrr => uvr(1:1)
    if (AU_MOINS_UN2(sr, valr)) then
        if (present(nr)) then
            unr = nr
            ptrr => valr(1:nr)
        else
            unr = 1
            uvr(1) = sr
            ptrr => uvr(1:1)
        end if
    end if
!
    ufname = ' '
    if (present(fname)) then
        ufname = fname
    end if

!  1. Faut-il completer le message (si on est dans un calcul elementaire) ?
!  ------------------------------------------------------------------------
    typ2 = typ
    if (calcul_status() .eq. 3 .and. (typ2(1:1) .eq. 'F' .or. typ2(1:1) .eq. 'E')) then
        under_te0000 = .true.
        if (typ2(2:2) .eq. '+') under_te0000 = .false.
    else
        under_te0000 = .false.
    end if
    if (under_te0000) then
        typ2(2:2) = '+'
    end if

!   2. Emission du message demande :
!   --------------------------------
    if (use_valk) then
        call utmess_core(typ2, idmess, unk, valk, uni, &
                         ptri, unr, ptrr, nexcep, ufname)
    else
        call utmess_core(typ2, idmess, unk, uvk, uni, &
                         ptri, unr, ptrr, nexcep, ufname)
    end if

!   3. Complement de message pour un calcul elementaire :
!   -----------------------------------------------------
    if (under_te0000) then
        call temess(typ2(1:1))
    end if
!
end subroutine utmess
