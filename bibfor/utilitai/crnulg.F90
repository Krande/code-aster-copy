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
!
subroutine crnulg(numddl)
    implicit none
#include "asterf_config.h"
#include "asterf.h"
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/crnlgc.h"
#include "asterfort/crnlgn.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jemarq.h"
#include "asterfort/asmpi_info.h"
    character(len=14) :: numddl
!
#ifdef ASTER_HAVE_MPI
!
    integer(kind=8) :: rang, nbproc
    mpi_int :: mrank, msize
!
!----------------------------------------------------------------------
    call jemarq()
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!   Création de la numérotation
    call crnlgn(numddl//'.NUME')

!   Communication de la numérotation
    call crnlgc(numddl//'.NUME')

!   Suppression des objets temporaires
    call jedetc('V', '&&CRNULG', 1)
!
    call jedema()
#else
    character(len=14) :: k14
    k14 = numddl
#endif
!
end subroutine
