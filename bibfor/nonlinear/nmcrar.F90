! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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
subroutine nmcrar(result, sddisc, fonact)
!
implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmarex.h"
#include "asterfort/nmarnr.h"
#include "asterfort/nmarpr.h"
#include "asterfort/nmcrpx.h"
#include "asterfort/nmdide.h"
#include "asterfort/wkvect.h"
!
character(len=19) :: sddisc
character(len=8) :: result
integer :: fonact(*)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (STRUCTURES DE DONNES)
!
! CREATION SD ARCHIVAGE
!
! --------------------------------------------------------------------------------------------------
!
! IN  RESULT : NOM DE LA SD RESULTAT
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IN  SDDISC : SD DISCRETISATION
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ifm, niv
    character(len=16), parameter :: keywfact = 'ARCHIVAGE', keywStep = 'PAS_ARCH'
    character(len=1), parameter :: base = 'V'
    integer :: nocc, iocc
    integer :: numder, numrep, numarc, numreo
    character(len=19) :: sdarch
    character(len=24) :: arcinf
    integer, pointer :: storeArcinf(:) => null()
    aster_logical :: lreuse, lDyna
    real(kind=8) :: insder
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I','MECANONLINE13_14')
    endif
!
! - Initializations
!
    iocc   = 1
    numarc = -1
    numreo = -1
    numrep = -1
    call getfac(keywfact, nocc)
    ASSERT(nocc .le. 1)
!
! - Active functionnalites
!
    lreuse = isfonc(fonact, 'REUSE')
    lDyna  = isfonc(fonact, 'DYNAMIQUE')
!
! - Name of datastructures to store
!
    sdarch = sddisc(1:14)//'.ARCH'
    arcinf = sdarch(1:19)//'.AINF'
!
! --- DERNIER NUMERO ARCHIVE DANS L'EVOL  SI REUSE
!
    call nmdide(lreuse, result, numder, insder)
!
! --- LECTURE LISTE INSTANTS D'ARCHIVAGE
!
    call nmcrpx(keywfact, keywStep, iocc, sdarch, base)
!
! --- CONSTRUCTION CHAMPS EXCLUS DE L'ARCHIVAGE
!
    call nmarex(keywfact, sdarch, lDyna)
!
! --- RECUPERATION DU PREMIER NUMERO A ARCHIVER
!
    call nmarpr(result, sddisc, lreuse, numder, insder, numarc)
!
! --- RECUPERATION NUMERO REUSE - TABLE OBSERVATION
!
    call nmarnr(result, 'OBSERVATION', numreo)
!
! --- RECUPERATION NUMERO REUSE - TABLE PARA_CALC
!
    call nmarnr(result, 'PARA_CALC', numrep)
!
! --- NUMERO D'ARCHIVE COURANT ET NUMERO DE REUSE
!
    ASSERT(numarc .ge. 0)
    ASSERT(numreo .ge. 0)
    ASSERT(numrep .ge. 0)
    call wkvect(arcinf, 'V V I', 3, vi = storeArcinf)
    storeArcinf(1) = numarc
    storeArcinf(2) = numreo
    storeArcinf(3) = numrep
!
end subroutine
