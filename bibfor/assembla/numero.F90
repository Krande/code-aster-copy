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

subroutine numero(nume_ddlz, base, &
                  old_nume_ddlz, modelocz, &
                  modelz, list_loadz, &
                  nb_matr_elem, list_matr_elem, &
                  nb_ligrel, list_ligrel, &
                  sd_iden_relaz)
!
    implicit none
!
#include "asterfort/as_deallocate.h"
#include "asterfort/crnulg.h"
#include "asterfort/dismoi.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/numer2.h"
#include "asterfort/numcch.h"
#include "asterfort/numoch.h"
#include "asterfort/uttcpu.h"
#include "asterfort/as_allocate.h"
!
! person_in_charge: jacques.pellet at edf.fr
!
    character(len=*), intent(inout) :: nume_ddlz
    character(len=2), intent(in) :: base
    character(len=*), optional, intent(in) :: modelz
    character(len=*), optional, intent(in) :: list_loadz
    character(len=24), optional, intent(in) :: list_matr_elem(*), list_ligrel(*)
    integer, optional, intent(in) :: nb_matr_elem, nb_ligrel
    character(len=*), optional, intent(in) :: old_nume_ddlz
    character(len=*), optional, intent(in) :: modelocz
    character(len=*), optional, intent(in) :: sd_iden_relaz
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering
!
! --------------------------------------------------------------------------------------------------
!
! IO  nume_ddl       : name of numbering object (NUME_DDL)
! In  base           : JEVEUX base to create objects
!                      base(1:1) => NUME_EQUA objects
!                      base(2:2) => NUME_DDL objects
! In  old_nume_ddl   : name of previous nume_ddl object
! In  modelocz       : local mode for GRANDEUR numbering
! In  model          : name of model
! In  list_load      : list of loads
! In  list_matr_elem : list of elementary matrixes
! In  list_ligrel    : list of ligrel
! In  nb_matr_elem   : number of elementary matrixes
! In  sd_iden_rela   : name of object for identity relations between dof
!
! If old_nume_ddl is present
!   -> try to know if NUME_EQUA in old_nume_ddl can be reuse
!      In this case nume_ddl = old_nume_ddl
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_ligr, igr
    character(len=8) :: nommai
    character(len=14) :: nume_ddl
    character(len=16) :: typsd
    character(len=24) :: modeloc, old_nume_ddl
    character(len=24), pointer :: list_ligr(:) => null()
    character(len=24) :: sd_iden_rela
!
! --------------------------------------------------------------------------------------------------
!
    call uttcpu('CPU.RESO.1', 'DEBUT', ' ')
    call uttcpu('CPU.RESO.2', 'DEBUT', ' ')
!
! - Identity relations between dof
!
    sd_iden_rela = ' '
    if (present(sd_iden_relaz)) then
        sd_iden_rela = sd_iden_relaz
    end if
!
! - Local mode
!
    modeloc = ' '
    if (present(modelocz)) then
        modeloc = modelocz
    end if
    old_nume_ddl = ' '
    if (present(old_nume_ddlz)) then
        old_nume_ddl = old_nume_ddlz
    end if
!
! - Create list of LIGREL for numbering
!
    if (present(list_matr_elem)) then
        call numoch(list_matr_elem, nb_matr_elem, list_ligr, nb_ligr)
    elseif (present(list_ligrel)) then
        ASSERT(present(nb_ligrel))
        AS_ALLOCATE(vk24=list_ligr, size=nb_ligrel)
        nb_ligr = nb_ligrel
        do igr = 1, nb_ligr
            list_ligr(igr) = list_ligrel(igr)
        end do
    else
        call numcch(modelz, list_loadz, list_ligr, nb_ligr)
    end if
!
! - Create numbering
!
    call numer2(nb_ligr, list_ligr, base, nume_ddlz, &
                old_nume_ddl, modeloc, modelz, sd_iden_rela)
!
    if (present(modelz)) then
        call dismoi('NOM_MAILLA', modelz, 'MODELE', repk=nommai)
        call gettco(nommai, typsd)
        if (typsd .eq. 'MAILLAGE_P') then
            nume_ddl = nume_ddlz
            call crnulg(nume_ddl)
        end if
    end if
!
    AS_DEALLOCATE(vk24=list_ligr)
    call uttcpu('CPU.RESO.1', 'FIN', ' ')
    call uttcpu('CPU.RESO.2', 'FIN', ' ')
!
end subroutine
