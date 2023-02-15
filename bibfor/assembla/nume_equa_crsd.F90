! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine nume_equa_crsd(nume_equaz, base, nb_equa, meshz, gran_namez)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jenonu.h"
#include "asterfort/wkvect.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/profchno_crsd.h"
!
!
    character(len=*), intent(in) :: nume_equaz
    character(len=1), intent(in) :: base
    integer, intent(in) :: nb_equa
    character(len=*), intent(in) :: meshz
    character(len=*), intent(in) :: gran_namez
!
! --------------------------------------------------------------------------------------------------
!
! NUME_EQUA
!
! Create object
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_equa    : name of NUME_EQUA
! In  base         : JEVEUX base to create NUME_EQUA
! In  nb_equa      : number of equations
! In  mesh         : name of mesh
! In  gran_name    : name of GRANDEUR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: nume_equa
    integer, pointer :: nequ(:) => null()
    character(len=24), pointer :: refn(:) => null()
    integer, pointer :: delg(:) => null()
    aster_logical :: l_pmesh
!
! --------------------------------------------------------------------------------------------------
!
    nume_equa = nume_equaz
    call detrsd('NUME_EQUA', nume_equa)
    l_pmesh = isParallelMesh(meshz)
    if (.not. l_pmesh) then
        ASSERT(nb_equa > 0)
    end if
!
    call profchno_crsd(nume_equa, base, nb_equa, meshz=meshz, gran_namez=gran_namez)
!
! - Create object NEQU
!
    call wkvect(nume_equa//'.NEQU', base//' V I', max(1, nb_equa), vi=nequ)
    nequ(1) = nb_equa
    nequ(2) = nb_equa
!
! - Create object REFN
!
    call wkvect(nume_equa//'.REFN', base//' V K24', 4, vk24=refn)
    refn(1) = meshz
    refn(2) = gran_namez
!
! - Create object DELG
!
    call wkvect(nume_equa//'.DELG', base//' V I', nb_equa, vi=delg)
    delg(1:nb_equa) = 0
!
end subroutine
