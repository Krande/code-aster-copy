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

subroutine nmcrcv(sdcrit)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
    character(len=19) :: sdcrit
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (STRUCTURES DE DONNES)
!
! CREATION SD POUR ARCHIVAGE DES INFORMATIONS DE CONVERGENCE
!
! ----------------------------------------------------------------------
!
!
! OUT SDCRIT : SD POUR ARCHIVAGE DES INFORMATIONS DE CONVERGENCE
!                (1) NOMBRE ITERATIONS NEWTON
!                (2) NOMBRE ITERATIONS RECHERCHE LINEAIRE
!                (3) RESI_GLOB_RELA
!                (4) RESI_GLOB_MAXI
!                (5) PARAMETRE DE PILOTAGE ETA
!                (6) CHARGEMENT EXTERIEUR
!                (9) RESI_COMP_RELA
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: jcrr, jcrk
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- CREATION
!
    call wkvect(sdcrit(1:19)//'.CRTR', 'V V R8', 9, jcrr)
    call wkvect(sdcrit(1:19)//'.CRDE', 'V V K16', 9, jcrk)
    zk16(jcrk+1-1) = 'ITER_GLOB'
    zk16(jcrk+2-1) = 'ITER_LINE'
    zk16(jcrk+3-1) = 'RESI_GLOB_RELA'
    zk16(jcrk+4-1) = 'RESI_GLOB'
    zk16(jcrk+5-1) = 'ETA_PILOTAGE'
    zk16(jcrk+6-1) = 'CHAR_MINI'
    zk16(jcrk+7-1) = 'RESI_GLOB_MOINS'
    zk16(jcrk+8-1) = 'RESI_REFE'
    zk16(jcrk+9-1) = 'RESI_COMP'
!
    call jedema()
!
end subroutine
