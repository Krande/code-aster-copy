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

subroutine op0000()
    implicit none
!   COMMANDE DEBUT OU POURSUITE
!
#include "asterc/getres.h"
#include "asterc/jdcget.h"
#include "asterf_types.h"
#include "asterfort/adjust_memlimit.h"
#include "asterfort/dbg_base.h"
#include "asterfort/foint0.h"
#include "asterfort/fozero.h"
#include "asterfort/gcncon.h"
#include "asterfort/getvtx.h"
#include "asterfort/ibbase.h"
#include "asterfort/ibcata.h"
#include "asterfort/ibdbgs.h"
#include "asterfort/ibtcpu.h"
#include "asterfort/jvinfo.h"
#include "asterfort/uldefi.h"
#include "asterfort/utmess.h"
    character(len=8) :: k8b
    character(len=16) :: nomcmd, k16b
    integer(kind=8) :: icode, ier, n, dummy
    integer(kind=8), save :: ipass = 0
!
    if (ipass .ne. 0) then
        call utmess('F', 'SUPERVIS_2')
    end if
    ipass = 1
!   to be set by 'ExecutionParameter().enable(Options.Debug)' or similar
    dummy = jvinfo(0)
!
    icode = jdcget('TestMode')
    if (icode .ne. 0) then
        call uldefi(15, ' ', 'CODE', 'A', 'A', 'O')
    end if

! --- LECTURE DU MOT CLE FACTEUR DEBUG OU DE GESTION MEMOIRE DEMANDE
    call ibdbgs()
!
! --- LECTURE DU MOT CLEF TEMPS_CPU
    call ibtcpu(ier)
!
! --- LECTURE DU MOT CLE FACTEUR BASE ET ---
! --- ALLOCATION DES BASES DE DONNEES ---
    call ibbase(ier)
    if (ier .eq. 0) then
        call getres(k8b, k16b, nomcmd)
!        -- INITIALISATION DE LA FONCTION NULLE : '&FOZERO'
!           ET DU COMMON FOSAV
        call fozero('&FOZERO')
        call foint0()
    end if
!
! --- POUR EVITER QUE LA CREATION DE '&&_NUM_CONCEPT_UNIQUE'
!        NE SOIT REPROCHE A UNE COMMANDE CREANT UNE SD
!        (DEBUT/DEBUG/SDVERI='OUI')
    call gcncon('.', k8b)
!
! --- LECTURE DU MOT CLE FACTEUR  CATALOGUE ---
    call ibcata(ier)
!
!
! --- DEBUG / VERI_BASE : lu ici et non dans ibdbgs car après ibbase
    call getvtx('DEBUG', 'VERI_BASE', iocc=1, nbret=n)
    if (n .eq. 1) then
        ! message et debug_jeveux forcé dans ibdbgs
        call dbg_base()
    end if
!
    call adjust_memlimit(ASTER_TRUE)
!
end subroutine
