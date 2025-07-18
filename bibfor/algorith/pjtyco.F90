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

subroutine pjtyco(isole, resuin, cham1, lnoeu, lelno, &
                  lelem, lelga)
! person_in_charge: jacques.pellet at edf.fr
!
! COMMANDE:  PROJ_CHAMP
! BUT : DETERMINER LES TYPES DE CHAMP A PROJETER
!
!
    implicit none

! 0.1. ==> ARGUMENTS

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutc4.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
    character(len=8) :: resuin
    character(len=19) :: cham1
    aster_logical :: isole
    aster_logical :: lnoeu, lelno, lelem, lelga

!  LNOEU  : .TRUE.  : IL Y A UN CHAM_NO A PROJETER
!  LELNO  : .TRUE.  : IL Y A UN CHAM_ELEM DE TYPE ELNO A PROJETER
!  LELEM  : .TRUE.  : IL Y A UN CHAM_ELEM DE TYPE ELEM A PROJETER
!  LELGA  : .TRUE.  : IL Y A UN CHAM_ELEM DE TYPE ELGA A PROJETER

! 0.2. ==> COMMUNS
! ----------------------------------------------------------------------

! 0.3. ==> VARIABLES LOCALES

    integer(kind=8) :: i, ie, iret
    integer(kind=8) :: nbordr
    integer(kind=8) :: iordr, isym, nbsym
    aster_logical :: acceno
    real(kind=8) :: prec
    character(len=4) :: tych
    character(len=8) :: crit
    character(len=16) :: nomsym(200)
    integer(kind=8), pointer :: nume_ordre(:) => null()

! DEB ------------------------------------------------------------------
    call jemarq()

    lnoeu = .false.
    lelno = .false.
    lelem = .false.
    lelga = .false.

!   1- CAS CHAMP ISOLE :
!   =====================
    if (isole) then
        call dismoi('TYPE_CHAMP', cham1, 'CHAMP', repk=tych)
        if (tych .eq. 'NOEU') then
            lnoeu = .true.
        else if (tych .eq. 'ELNO') then
            lelno = .true.
        else if (tych .eq. 'ELEM') then
            lelem = .true.
        else if (tych .eq. 'ELGA') then
            lelga = .true.
        else
            ASSERT(.false.)
        end if

!   2- CAS SD_RESULTAT :
!   =====================
    else
        call getvr8(' ', 'PRECISION', scal=prec, nbret=ie)
        call getvtx(' ', 'CRITERE', scal=crit, nbret=ie)
        call rsutnu(resuin, ' ', 0, '&&PJXXCO.NUME_ORDRE', nbordr, &
                    prec, crit, iret)

        if (iret .ne. 0) then
            call utmess('F', 'CALCULEL4_61', sk=resuin)
        end if
        if (nbordr .eq. 0) then
            call utmess('F', 'CALCULEL4_62', sk=resuin)
        end if

        call jeveuo('&&PJXXCO.NUME_ORDRE', 'L', vi=nume_ordre)
        call rsutc4(resuin, ' ', 1, 200, nomsym, &
                    nbsym, acceno)

!       -- DETERMINATION DE LNOEU
        do isym = 1, nbsym
            do i = 1, nbordr
                iordr = nume_ordre(i)
                call rsexch(' ', resuin, nomsym(isym), iordr, cham1, &
                            iret)

                if (iret .eq. 0) then
                    call dismoi('TYPE_CHAMP', cham1, 'CHAMP', repk=tych)
                    if (tych .eq. 'NOEU') then
                        lnoeu = .true.
                        goto 20

                    end if
                end if

            end do
20          continue
        end do

!       -- DETERMINATION DE LELNO
        do isym = 1, nbsym
            do i = 1, nbordr
                iordr = nume_ordre(i)
                call rsexch(' ', resuin, nomsym(isym), iordr, cham1, &
                            iret)

                if (iret .eq. 0) then
                    call dismoi('TYPE_CHAMP', cham1, 'CHAMP', repk=tych)
                    if (tych .eq. 'ELNO') then
                        lelno = .true.
                        goto 40

                    end if
                end if

            end do
40          continue
        end do

!       -- DETERMINATION DE LELEM
        do isym = 1, nbsym
            do i = 1, nbordr
                iordr = nume_ordre(i)
                call rsexch(' ', resuin, nomsym(isym), iordr, cham1, &
                            iret)

                if (iret .eq. 0) then
                    call dismoi('TYPE_CHAMP', cham1, 'CHAMP', repk=tych)
                    if (tych .eq. 'ELEM') then
                        lelem = .true.
                        goto 60

                    end if
                end if

            end do
60          continue
        end do

!       -- DETERMINATION DE LELGA
        do isym = 1, nbsym
            do i = 1, nbordr
                iordr = nume_ordre(i)
                call rsexch(' ', resuin, nomsym(isym), iordr, cham1, &
                            iret)

                if (iret .eq. 0) then
                    call dismoi('TYPE_CHAMP', cham1, 'CHAMP', repk=tych)
                    if (tych .eq. 'ELGA') then
                        lelga = .true.
                        goto 80

                    end if
                end if

            end do
80          continue
        end do
    end if

    if (.not. lelga .and. .not. lelno .and. .not. lelem .and. .not. lnoeu) then
        call utmess('F', 'CALCULEL4_75', sk=resuin)
    end if

    call jedema()
end subroutine
