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
!  Person in charge: mickael.abbas at edf.fr
!
subroutine cagene(load, command, model, mesh, geomDime)
!
    implicit none
!
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: load
    character(len=16), intent(in) :: command
    integer(kind=8), intent(out) :: geomDime
    character(len=8), intent(out) :: mesh
    character(len=8), intent(out) :: model
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Get infos
!
! --------------------------------------------------------------------------------------------------
!
! In  load      : load
! In  command   : command
! Out mesh      : mesh
! Out model     : model
! Out geomDime  : space dimension
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: nomoJv, phenomenon, valk(2)
    character(len=8), pointer :: nomo(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Get model
    call getvid(' ', 'MODELE', scal=model)

! - Mesh
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Check model/loading
    call dismoi('PHENOMENE', model, 'MODELE', repk=phenomenon)
    valk(1) = command
    valk(2) = phenomenon
    if (command .eq. 'AFFE_CHAR_THER' .and. phenomenon .ne. 'THERMIQUE') then
        call utmess('F', 'CHARGES2_64', nk=2, valk=valk)
    else if (command .eq. 'AFFE_CHAR_MECA' .and. phenomenon .ne. 'MECANIQUE') then
        call utmess('F', 'CHARGES2_64', nk=2, valk=valk)
    else if (command .eq. 'DEFI_CONTACT' .and. phenomenon .ne. 'MECANIQUE') then
        call utmess('F', 'CHARGES2_64', nk=2, valk=valk)
    else if (command .eq. 'AFFE_CHAR_ACOU' .and. phenomenon .ne. 'ACOUSTIQUE') then
        call utmess('F', 'CHARGES2_64', nk=2, valk=valk)
    end if

! - Dimension of problem
    call dismoi('DIM_GEOM', model, 'MODELE', repi=geomDime)

! - Create .NOMO
    nomoJv = load(1:8)//'.CH'//phenomenon(1:2)//'.MODEL.NOMO'
    call wkvect(nomoJv, 'G V K8', 1, vk8=nomo)
    nomo(1) = model
!
    call jedema()
end subroutine
