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

subroutine pjcovide(noma1, noma2, corres)
    implicit none

! ----------------------------------------------------------------------
!     COMMANDE:  PROJ_CHAMP  METHODE:'ELEM'
! BUT : CREER UNE STRUCTURE DE DONNEE CORRESP_2_MAILLA VIDE
! ----------------------------------------------------------------------
!
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
#include "asterfort/isParallelMesh.h"
    character(len=8) :: noma1, noma2
    character(len=16) :: corres
!
! 0.3. ==> VARIABLES LOCALES
!
!
    aster_logical :: l_parallel_mesh
    integer(kind=8) :: jxxk1, nbno
    integer(kind=8), pointer :: pjef_nb(:) => null()
    integer(kind=8), pointer :: pjef_nu(:) => null()
    real(kind=8), pointer :: pjef_cf(:) => null()
    integer(kind=8), pointer :: pjef_tr(:) => null()
!----------------------------------------------------------------------
    call jemarq()

    ASSERT(isParallelMesh(noma1))

    call wkvect(corres//'.PJXX_K1', 'V V K24', 5, jxxk1)
    zk24(jxxk1-1+1) = noma1
    zk24(jxxk1-1+2) = noma2
    zk24(jxxk1-1+3) = 'COLLOCATION'
    call dismoi('NB_NO_MAILLA', noma1, 'MAILLAGE', nbno)
    call wkvect(corres//'.PJEF_NB', 'V V I', nbno, vi=pjef_nb)
    call wkvect(corres//'.PJEF_NU', 'V V I', 4, vi=pjef_nu)
    call wkvect(corres//'.PJEF_CF', 'V V R', 4, vr=pjef_cf)
    call wkvect(corres//'.PJEF_TR', 'V V I', 1, vi=pjef_tr)
!
    call jedema()
end subroutine
