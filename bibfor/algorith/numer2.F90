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

subroutine numer2(nb_ligr, list_ligr, base, nume_ddlz, &
                  nume_ddl_oldz, modelocz, modele, sd_iden_relaz)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterc/cheksd.h"
#include "asterfort/detrsd.h"
#include "asterfort/idensd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jemarq.h"
#include "asterfort/matdis.h"
#include "asterfort/nueffe.h"
#include "asterfort/nugllo.h"
#include "asterfort/promor.h"
!
! person_in_charge: jacques.pellet at edf.fr
!
    integer, intent(in) :: nb_ligr
    character(len=24), pointer :: list_ligr(:)
    character(len=2), intent(in) :: base
    character(len=*), intent(inout) :: nume_ddlz
    character(len=*), intent(in) :: nume_ddl_oldz
    character(len=*), intent(in) :: modelocz
    character(len=*), intent(in) :: modele
    character(len=*), optional, intent(in) :: sd_iden_relaz
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering - Create objects
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_ligr        : number of LIGREL in list
! In  list_ligr      : pointer to list of LIGREL
! In  base           : JEVEUX base to create objects
!                      base(2:2) => NUME_EQUA objects
!                      base(1:1) => NUME_DDL objects
! IO  nume_ddl       : name of numbering object (NUME_DDL)
! In  modeloc        : local mode for GRANDEUR numbering
! In  nume_ddl_old   : name of previous nume_ddl object
! In  sd_iden_rela   : name of object for identity relations between dof
!
! If nume_ddl_old is present
!   -> try to know if NUME_EQUA in nume_ddl_old can be reuse
!      In this case nume_ddl = nume_ddl_old
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: nume_equa, nume_equa_old
    character(len=14) :: nume_ddl, nume_ddl_old, moloc
    character(len=24) :: sd_iden_rela
    character(len=3) :: matd
    aster_logical :: l_matr_dist, printt
    aster_logical, parameter :: verbose = ASTER_FALSE, debug = ASTER_FALSE
    integer :: iret
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    nume_ddl = nume_ddlz
    moloc = modelocz
    nume_ddl_old = nume_ddl_oldz
    nume_equa = nume_ddl//'.NUME'
    nume_equa_old = nume_ddl_old//'.NUME'
!
    call detrsd('NUME_DDL', nume_ddl)
!
! - Identity relations between dof
!
    sd_iden_rela = ' '
    if (present(sd_iden_relaz)) then
        sd_iden_rela = sd_iden_relaz
    end if
!
    call matdis(matd, verbose)
    ASSERT(matd .eq. 'OUI' .or. matd .eq. 'NON')
    if (matd .eq. 'OUI') then
        l_matr_dist = .true.
    else
        l_matr_dist = .false.
    end if
!
! - Create NUME_EQUA objects
!
    call nueffe(nb_ligr, list_ligr, base, nume_ddl, 'SANS', &
                modele, modelocz=moloc, sd_iden_relaz=sd_iden_rela)
!
! - Create NUML_EQUA objects
!
    if (l_matr_dist) then
        call nugllo(nume_ddlz, base)
    end if
!
! - Trying to reuse old nume_ddl
!
    if (nume_ddl_old .ne. ' ') then
        if (idensd('NUME_EQUA', nume_equa, nume_equa_old)) then
            call detrsd('NUME_DDL', nume_ddl)
            call jedupo(nume_ddl//'     .ADNE', 'V', nume_ddl_old//'     .ADNE', .false._1)
            call jedupo(nume_ddl//'     .ADLI', 'V', nume_ddl_old//'     .ADLI', .false._1)
            call jedetr(nume_ddl//'     .ADLI')
            call jedetr(nume_ddl//'     .ADNE')
            nume_ddl = nume_ddl_old
        end if
    end if
!
! - Create matrix topology
!
    printt = moloc .eq. ' '
    call promor(nume_ddl, base(1:1), printt)
!
! - Cleaning
!
    call jedetr(nume_ddl//'     .ADLI')
    call jedetr(nume_ddl//'     .ADNE')
!
    nume_ddlz = nume_ddl
!
    if (debug) call cheksd(nume_ddlz, 'SD_NUME_DDL', iret)
!
    call jedema()
end subroutine
