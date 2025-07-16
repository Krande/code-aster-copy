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

subroutine nume_equa_crsd(nume_equaz, base, nb_equa, meshz, gran_namez, l_coll_constz)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jeecra.h"
#include "asterfort/wkvect.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/profchno_crsd.h"
!
!
    character(len=*), intent(in) :: nume_equaz
    character(len=1), intent(in) :: base
    integer(kind=8), intent(in) :: nb_equa
    character(len=*), intent(in) :: meshz
    character(len=*), intent(in) :: gran_namez
    aster_logical, optional, intent(in) :: l_coll_constz
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
    integer(kind=8) :: nb_equa2
    character(len=19) :: nume_equa
    integer(kind=8), pointer :: nequ(:) => null()
    character(len=24), pointer :: refn(:) => null()
    integer(kind=8), pointer :: delg(:) => null()
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
    if (present(l_coll_constz)) then
        call profchno_crsd(nume_equa, base, nb_equa, meshz=meshz, gran_namez=gran_namez, &
                           l_coll_const=l_coll_constz)
    else
        call profchno_crsd(nume_equa, base, nb_equa, meshz=meshz, gran_namez=gran_namez)
    end if
!
! - Create object NEQU
!
    call wkvect(nume_equa//'.NEQU', base//' V I', 2, vi=nequ)
    nequ(1) = nb_equa
    nequ(2) = nb_equa
!
! - Create object REFN
!
    call wkvect(nume_equa//'.REFN', base//' V K24', 5, vk24=refn)
    refn(1) = meshz
    refn(2) = gran_namez
!
! - Create object DELG
!
    nb_equa2 = nb_equa
    if (l_pmesh .and. nb_equa .eq. 0) then
        nb_equa2 = 1
    end if
    call wkvect(nume_equa//'.DELG', base//' V I', nb_equa2, vi=delg)
    call jeecra(nume_equa//'.DELG', 'LONUTI', nb_equa)
    delg(1:nb_equa2) = 0
!
end subroutine
