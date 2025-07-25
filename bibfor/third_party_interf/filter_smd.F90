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

subroutine filter_smd(nommat, vsmb)
!
#include "asterf_petsc.h"
!
! person_in_charge: natacha.bereux at edf.fr
    use aster_petsc_module
    use petsc_data_module

    implicit none

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=*) :: nommat
    real(kind=8) :: vsmb(*)
!-----------------------------------------------------------------------
! BUT : ON MET A ZERO LES TERMES DU SECOND MEMBRE QUI N'APPARTIENNENT PAS
!       DE FACON EXCLUSIVE (AU SENS PETSC) AU PROCESSEUR COURANT.
!-----------------------------------------------------------------------
!     VARIABLES LOCALES
!-----------------------------------------------------------------------
    integer(kind=8) :: ieql, ieqg, jpddl, neqg, neql
    integer(kind=8) :: iccid, rang
    character(len=14) :: nu
    character(len=19) :: mat
    mpi_int :: mrank, msize
    aster_logical :: is_ddl_cine, iam_sole_owner
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: ccid(:) => null()
    integer(kind=8), pointer :: nulg(:) => null()
    integer(kind=8), pointer :: nequl(:) => null()
    integer(kind=8), pointer :: nequ(:) => null()
!-----------------------------------------------------------------------
!     DEBUT
    call jemarq()
!-----------------------------------------------------------------------
    mat = nommat
!
    call jeveuo(mat//'.REFA', 'L', vk24=refa)
    if (refa(11) .eq. 'MATR_DISTR') then
! Infos du processeur courant
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
! Infos du NUME_DDL
        nu = refa(2) (1:14)
        call jeveuo(nu//'.NUML.NULG', 'L', vi=nulg)
        call jeveuo(nu//'.NUML.PDDL', 'L', jpddl)
        call jeveuo(nu//'.NUML.NEQU', 'L', vi=nequl)
        call jeveuo(nu//'.NUME.NEQU', 'L', vi=nequ)
        neqg = nequ(1)
        neql = nequl(1)
        call jeexin(mat//'.CCID', iccid)
!
        if (iccid .ne. 0) then
            call jeveuo(mat//'.CCID', 'L', vi=ccid)
        end if
!
        do ieql = 1, neql
            ieqg = nulg(ieql)
! Le dl courant est-il fixé par une charge cinématique ?
            if (iccid == 0) then
! Il n'y a pas de charge cinématique sur le modèle
                is_ddl_cine = .false.
            else
! Il existe au moins une charge cinématique. On vérifie si
! le numéro global de dl est concerné
                is_ddl_cine = ccid(ieqg) .eq. 1
            end if
! Suis-je le proriétaire exclusif (PETSc) de ce ddl ?
            iam_sole_owner = zi(jpddl-1+ieql) .eq. rang
            if ((.not. is_ddl_cine) .and. (.not. iam_sole_owner)) then
                vsmb(ieqg) = 0.d0
            end if
        end do
    end if
!
    call jedema()
end subroutine
