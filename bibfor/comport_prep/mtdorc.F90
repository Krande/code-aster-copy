! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine mtdorc(model, comporMeta)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/comp_init.h"
#include "asterfort/comp_meta_info.h"
#include "asterfort/comp_meta_read.h"
#include "asterfort/comp_meta_save.h"
#include "asterfort/comp_meta_pvar.h"
#include "asterfort/comp_meta_prnt.h"
#include "asterfort/dismoi.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: comporMeta
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Construct map for behaviour in metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  comporMeta       : name of map for behaviour in metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nbCmp
    character(len=8) :: mesh
    character(len=19), parameter :: comporMetaInfo = '&&MTDORC.INFO'
    type(META_PrepPara) :: metaPrepPara
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Create datastructure to prepare comportement
    call comp_meta_info(metaPrepPara)

! - Create COMPOR <CARTE>
    call comp_init(mesh, comporMeta, 'V', nbCmp)

! - Read informations from command file
    call comp_meta_read(metaPrepPara)

! - Save informations in COMPOR <CARTE>
    call comp_meta_save(mesh, comporMeta, nbCmp, metaPrepPara)

! - Prepare informations about internal variables
    call comp_meta_pvar(model, comporMeta, comporMetaInfo)

! - Print informations about internal variables
    call comp_meta_prnt(comporMetaInfo)

! - Clean
    deallocate (metaPrepPara%para)
!
end subroutine
