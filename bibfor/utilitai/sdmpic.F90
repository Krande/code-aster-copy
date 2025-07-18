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

subroutine sdmpic(typesd, nomsd)
    implicit none
#include "jeveux.h"
#include "asterfort/asmpi_comm_jev.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
    character(len=*) :: nomsd, typesd
! person_in_charge: jacques.pellet at edf.fr
!
! ----------------------------------------------------------------------
!  BUT : "COMPLETER" LE CALCUL D'UNE STRUCTURE DE DONNEES INCOMPLETEMENT
!        CALCULEE  DU FAIT DE L'UTILISATION DE CALCULS PARALLELES (MPI)
!
!  LA ROUTINE ECHANGE LES MORCEAUX CALCULES PARTIELLEMENT SUR LES
!  DIFFERENTS PROCESSEURS (MPI_ALLREDUCE) SAUF SI LE MAILLAGE SUPPORT
!  EST UN PARALLEL_MESH
!
! ----------------------------------------------------------------------
! IN TYPESD (K*) :  TYPE DE LA SD A COMPLETER
! IN NOMSD  (K*) :  NOM DE LA SD A COMPLETER
! ----------------------------------------------------------------------
    character(len=24) :: noms2, types2
    character(len=19) :: k19
    character(len=24) :: k24
    character(len=8) :: kmpic, kbid, mesh
    integer(kind=8) :: ifm, niv, iexi, nbrel, i
    character(len=24), pointer :: noli(:) => null()
    character(len=24), pointer :: refa(:) => null()
    character(len=24), pointer :: celk(:) => null()
    character(len=24), pointer :: relr(:) => null()
    character(len=16), pointer :: valk(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
    noms2 = nomsd
    types2 = typesd
!
!
    k19 = noms2(1:19)
    k24 = noms2(1:24)
    if (types2 .eq. 'CHAM_ELEM') then
!     ----------------------------------
        call dismoi('MPI_COMPLET', k19, 'CHAM_ELEM', repk=kmpic)
        if (kmpic .eq. 'OUI') goto 999
        call dismoi('NOM_MAILLA', k19, 'CHAM_ELEM', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call asmpi_comm_jev('MPI_SUM', k19//'.CELV')
        call jeveuo(k19//'.CELK', 'E', vk24=celk)
        celk(7) = 'MPI_COMPLET'
!
    else if (types2 .eq. 'RESUELEM') then
!     ----------------------------------
        call dismoi('MPI_COMPLET', k19, 'RESUELEM', repk=kmpic)
        if (kmpic .eq. 'OUI') goto 999
        call dismoi('NOM_MAILLA', k19, 'RESUELEM', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call asmpi_comm_jev('MPI_SUM', k19//'.RESL')
        call jeveuo(k19//'.NOLI', 'E', vk24=noli)
        noli(3) = 'MPI_COMPLET'
!
    else if (types2 .eq. 'MATR_ELEM') then
!     ----------------------------------
        call dismoi('NOM_MAILLA', k19, 'MATR_ELEM', repk=mesh)
        if (.not. isParallelMesh(mesh)) then
            call jeveuo(k19//'.RELR', 'L', vk24=relr)
            call jelira(k19//'.RELR', 'LONMAX', nbrel, kbid)
            do i = 1, nbrel
                if (relr(i) .ne. ' ') then
                    k19 = relr(i) (1:19)
                    call dismoi('MPI_COMPLET', k19, 'RESUELEM', repk=kmpic)
                    if (kmpic .eq. 'OUI') cycle
                    call asmpi_comm_jev('MPI_SUM', k19//'.RESL')
                    call jeveuo(k19//'.NOLI', 'E', vk24=noli)
                    noli(3) = 'MPI_COMPLET'
                end if
            end do
        end if
    else if (types2 .eq. 'VECT_ELEM') then
!     ----------------------------------
        call dismoi('NOM_MAILLA', k19, 'VECT_ELEM', repk=mesh)
        if (.not. isParallelMesh(mesh)) then
            call jeveuo(k19//'.RELR', 'L', vk24=relr)
            call jelira(k19//'.RELR', 'LONMAX', nbrel, kbid)
            do i = 1, nbrel
                if (relr(i) .ne. ' ') then
                    k19 = relr(i) (1:19)
                    call dismoi('MPI_COMPLET', k19, 'RESUELEM', repk=kmpic)
                    if (kmpic .eq. 'OUI') cycle
                    call asmpi_comm_jev('MPI_SUM', k19//'.RESL')
                    call jeveuo(k19//'.NOLI', 'E', vk24=noli)
                    noli(3) = 'MPI_COMPLET'
                end if
            end do
        end if
!
    else if (types2 .eq. 'MATR_ASSE') then
!     ----------------------------------
        call dismoi('MPI_COMPLET', k19, 'MATR_ASSE', repk=kmpic)
        if (kmpic .eq. 'OUI') goto 999
        call dismoi('NOM_MAILLA', k19, 'MATR_ASSE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call asmpi_comm_jev('MPI_SUM', k19//'.VALM')
!
        call jeexin(k19//'.CCVA', iexi)
        if (iexi .gt. 0) call asmpi_comm_jev('MPI_SUM', k19//'.CCVA')
!
        call jeveuo(k19//'.REFA', 'E', vk24=refa)
        refa(11) = 'MPI_COMPLET'
!
    else if (types2 .eq. 'SD_APPA') then
!     ----------------------------------
        call jeveuo(k19//'.MPIA', 'E', vk16=valk)
        if (valk(1) .eq. 'MPI_COMPLET') goto 999
        call asmpi_comm_jev('MPI_SUM', k19//'.APPA')
        call asmpi_comm_jev('MPI_SUM', k19//'.DIST')
        call asmpi_comm_jev('MPI_SUM', k19//'.TAU1')
        call asmpi_comm_jev('MPI_SUM', k19//'.TAU2')
        call asmpi_comm_jev('MPI_SUM', k19//'.PROJ')
        valk(1) = 'MPI_COMPLET'
!
    else if (types2 .eq. 'SD_APPA_TGEL') then
!     ----------------------------------
        call jeveuo(k19//'.MPIB', 'E', vk16=valk)
        if (valk(1) .eq. 'MPI_COMPLET') goto 999
        call asmpi_comm_jev('MPI_SUM', k19//'.TGEL')
        valk(1) = 'MPI_COMPLET'
!
    else if (types2 .eq. 'SD_APPA_TGNO') then
!     ----------------------------------
        call jeveuo(k19//'.MPIC', 'E', vk16=valk)
        if (valk(1) .eq. 'MPI_COMPLET') goto 999
        call asmpi_comm_jev('MPI_SUM', k19//'.TGNO')
        valk(1) = 'MPI_COMPLET'
!
    else if (types2 .eq. 'SD_APPA_LAC1') then
!     ----------------------------------
        call jeveuo(k19//'.MPID', 'E', vk16=valk)
        if (valk(1) .eq. 'MPI_COMPLET') goto 999
        call asmpi_comm_jev('MPI_SUM', k19//'.PWT ')
        call asmpi_comm_jev('MPI_SUM', k19//'.NAPP')
        valk(1) = 'MPI_COMPLET'
!
!
!
    else if (types2 .eq. 'SD_APPA_LAC2') then
!     ----------------------------------
        call jeveuo(k19//'.MPIE', 'E', vk16=valk)
        if (valk(1) .eq. 'MPI_COMPLET') goto 999
        call asmpi_comm_jev('MPI_SUM', k19//'.AUX ')
        call asmpi_comm_jev('MPI_SUM', k19//'.AUX2')
        call asmpi_comm_jev('MPI_SUM', k19//'.AUX3')
        call asmpi_comm_jev('MPI_SUM', k19//'.AUX4')
        call asmpi_comm_jev('MPI_SUM', k19//'.AUX5')
        valk(1) = 'MPI_COMPLET'

    else if (types2 .eq. '&&OP29NL') then
!     ----------------------------------
        call jeveuo(k24(1:20)//'.MPI', 'E', vk16=valk)
        if (valk(1) .eq. 'MPI_COMPLET') goto 999
        call asmpi_comm_jev('MPI_SUM', k24)
        valk(1) = 'MPI_COMPLET'
!
!
    else
        ASSERT(.false.)
    end if
!
999 continue
    call jedema()
end subroutine
