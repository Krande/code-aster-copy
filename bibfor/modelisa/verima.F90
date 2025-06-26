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
subroutine verima(meshz, list_obj, list_size, typez_objet)
!
    implicit none
!
#include "jeveux.h"
!
#include "asterfort/assert.h"
#include "asterfort/cleanListOfGrpMa.h"
#include "asterfort/cleanListOfGrpNo.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
!
    character(len=*), intent(in) :: meshz, typez_objet
    integer(kind=8), intent(inout) :: list_size
    character(len=*), intent(inout) :: list_obj(list_size)
!
!     VERIFICATION DE L'APPARTENANCE DES OBJETS DE LA LISTE
!     LIMANO AU MAILLAGE MESH
!
! IN       : MESHZ    : NOM DU MAILLAGE
! IO       : LIMANZ   : LISTE DE MAILLES OU DE NOEUDS OU DE GROUP_NO
!                       OU DE GROUP_MA
! IO       : LIST_SIZE: LONGUEUR DE LA LISTE LIMANO
! IN       : TYPZ     : TYPE DES OBJETS DE LA LISTE :
!                       MAILLE OU NOEUD OU GROUP_NO OU GROUP_MA
! ----------------------------------------------------------------------
!
    integer(kind=8) :: iret, ino, ima, nbno, nbma, ier2
    character(len=8) :: mesh, type_obj
    character(len=24) :: noeuma, grnoma, mailma, grmama, object
    character(len=24) :: valk(2)
    aster_logical :: l_parallel_mesh, lnomnoe, lnommai
    aster_logical, parameter :: l_stop = ASTER_TRUE
! ----------------------------------------------------------------------
!
    mesh = meshz
    type_obj = typez_objet
!
    if (list_size .lt. 1) then
        goto 999
    end if
!
    noeuma = mesh//'.NOMNOE'
    grnoma = mesh//'.GROUPENO'
    mailma = mesh//'.NOMMAI'
    grmama = mesh//'.GROUPEMA'
    call jeexin(noeuma, ier2)
    lnomnoe = .false.
    if (ier2 .ne. 0) then
        lnomnoe = .true.
    end if
    call jeexin(mailma, ier2)
    lnommai = .false.
    if (ier2 .ne. 0) then
        lnommai = .true.
    end if
!
    l_parallel_mesh = isParallelMesh(mesh)
!
    if (type_obj .eq. 'GROUP_NO') then
!
!      --VERIFICATION DE L'APPARTENANCE DES GROUP_NO
!        AUX GROUP_NO DU MAILLAGE
!        -------------------------------------------------------
!
        if (.not. l_parallel_mesh) then
            call jeexin(grnoma, iret)
            if ((list_size .ne. 0) .and. (iret .eq. 0)) then
                valk(1) = type_obj
                valk(2) = mesh
                call utmess('F', 'MODELISA7_12', nk=2, valk=valk)
            end if
        end if
!
        call cleanListOfGrpNo(mesh, list_obj, list_size, l_stop, iret)
!
        if (.not. l_parallel_mesh) then
            ASSERT(iret == 1)
        end if
!
    else if (type_obj .eq. 'NOEUD') then
!
!      --VERIFICATION DE L'APPARTENANCE DES NOEUDS
!        AUX NOEUDS DU MAILLAGE
!        -------------------------------------------------------
        if (l_parallel_mesh .and. list_size .ne. 0) then
            call utmess('F', 'MODELISA7_86')
        end if
!
        call jelira(mesh//'.COORDO    .VALE', 'LONMAX', nbno)
        nbno = nbno/3
        do ino = 1, list_size
            object = list_obj(ino)
            if (lnomnoe) then
                call jenonu(jexnom(noeuma, object), iret)
            else
                iret = char8_to_int(object)
            end if
            if ((iret .gt. nbno) .or. (iret .le. 0)) then
                valk(1) = object
                valk(2) = mesh
                call utmess('F', 'MODELISA7_76', nk=2, valk=valk)
            end if
        end do
!
    else if (type_obj .eq. 'GROUP_MA') then
!
!      --VERIFICATION DE L'APPARTENANCE DES GROUP_MA
!        AUX GROUP_MA DU MAILLAGE
!        -------------------------------------------------------
        if (.not. l_parallel_mesh) then
            call jeexin(grmama, iret)
            if ((list_size .ne. 0) .and. (iret .eq. 0)) then
                valk(1) = type_obj
                valk(2) = mesh
                call utmess('F', 'MODELISA7_12', nk=2, valk=valk)
            end if
        end if
!
        call cleanListOfGrpMa(mesh, list_obj, list_size, l_stop, iret)
!
        if (.not. l_parallel_mesh) then
            ASSERT(iret == 1)
        end if
!
    else if (type_obj .eq. 'MAILLE') then
!
!      --VERIFICATION DE L'APPARTENANCE DES MAILLES
!        AUX MAILLES DU MAILLAGE
!        -------------------------------------------------------
        if (l_parallel_mesh .and. list_size .ne. 0) then
            call utmess('F', 'MODELISA7_86')
        end if
!
        call jelira(mesh//'.TYPMAIL', 'LONMAX', nbma)
        do ima = 1, list_size
            object = list_obj(ima)
            if (lnommai) then
                call jenonu(jexnom(mailma, object), iret)
            else
                iret = char8_to_int(object)
            end if
            if ((iret .gt. nbma) .or. (iret .le. 0)) then
                valk(1) = object
                valk(2) = mesh
                call utmess('F', 'MODELISA6_10', nk=2, valk=valk)
            end if
        end do
!
    else
        call utmess('F', 'MODELISA7_79', sk=type_obj)
    end if
999 continue
end subroutine
