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
subroutine nmprof(model, result, listload, numeDof, ds_contact)
!
    use listLoad_module
    use NonLin_Datastructure_type
    use contactAlgo_module
!
    implicit none
!
#include "asterfort/addModelLigrel.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/gnomsd.h"
#include "asterfort/numero.h"
#include "asterfort/rsnume.h"
!
    character(len=24), intent(in) :: model
    character(len=24), intent(out) :: numeDof
    character(len=8), intent(in) :: result
    character(len=19), intent(in) :: listload
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm - Initializations
!
! Create numbering
!
! --------------------------------------------------------------------------------------------------
!
! Out numeDof          : name of numbering object (NUME_DDL)
! In  result           : name of result datastructure (EVOL_NOLI)
! In  model            : name of model datastructure
! In  listload         : list of loads
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=14) :: numeDofOld
    character(len=24) :: noojb, idenRela
    integer(kind=8) :: nbLigr
    character(len=24), pointer :: listLigr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    idenRela = ds_contact%iden_rela

! - Add LIGREL from model
    nbLigr = 0
    call addModelLigrel(model, nbLigr, listLigr)

! - Get list of LIGREL from loads
    call getListLoadLigrel(listLoad, nbLigr, listLigr)

! - Get list of LIGREL from contact
    if (ds_contact%l_contact) then
        call addContactLigrel(ds_contact%sdcont, ds_contact%ligrel_elem_cont, nbLigr, listLigr)
    end if

! - Generate name of numbering object (numeDof)
    numeDof = '12345678.NUMED'
    noojb = '12345678.00000.NUME.PRNO'
    call gnomsd(' ', noojb, 10, 14)
    numeDof = noojb(1:14)

! - Get numbering from displacement
    call rsnume(result, 'DEPL', numeDofOld)

! - Create numbering
    call numero(numeDof, 'VG', &
                nbLigr, listLigr, &
                numeDofOldZ_=numeDofOld, idenRelaZ_=idenRela)
    AS_DEALLOCATE(vk24=listLigr)
!
end subroutine
