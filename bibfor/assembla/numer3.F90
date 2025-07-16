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
subroutine numer3(modelZ, base, listLoadZ, numeDofZ, ds_contact)
!
    use NonLin_Datastructure_type
    use listLoad_module
    use contactAlgo_module
!
    implicit none
!
#include "asterfort/addModelLigrel.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/idenob.h"
#include "asterfort/infbav.h"
#include "asterfort/infmue.h"
#include "asterfort/numero.h"
!
    character(len=*), intent(in) :: modelZ
    character(len=2), intent(in) :: base
    character(len=*), intent(inout) :: numeDofZ
    character(len=*), intent(in) :: listLoadZ
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! (Re)-Numbering - Used for variable element topology (contact)
!
! --------------------------------------------------------------------------------------------------
!
! IO  numeDof          : name of numbering object (NUME_DDL)
! In  base             : JEVEUX base to create objects
!                        base(1:1) => NUME_EQUA objects
!                        base(2:2) => NUME_DDL objects
! In  model            : name of model datastructure
! In  listload         : list of loads
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=14), parameter :: numeDofSave = '&&NUMER3.NUAV'
    integer(kind=8) :: nbLigr
    character(len=24), pointer :: listLigr(:) => null()
    character(len=14) :: numeDofOld
    character(len=24) :: idenRela
!
! --------------------------------------------------------------------------------------------------
!
    call infmue()

! - Management of numeDof
    numeDofOld = numeDofZ
    idenRela = ds_contact%iden_rela
    call copisd('NUME_DDL', 'V', numeDofZ, numeDofSave)
    call detrsd('NUME_DDL', numeDofZ)

! - Add LIGREL from model
    nbLigr = 0
    call addModelLigrel(modelZ, nbLigr, listLigr)

! - Get list of LIGREL from loads
    call getListLoadLigrel(listLoadZ, nbLigr, listLigr)

! - Get list of LIGREL from contact
    if (ds_contact%l_contact) then
        call addContactLigrel(ds_contact%sdcont, ds_contact%ligrel_elem_cont, nbLigr, listLigr)
    end if

! - Numbering
    call numero(numeDofZ, base, &
                nbLigr, listLigr, &
                idenRelaZ_=idenRela)
    AS_DEALLOCATE(vk24=listLigr)

! - Same equations ! The re-numbering works only with MUMPS/MULT_FRONT/PETSc, not with LDLT
    ASSERT(idenob(numeDofOld//'.NUME.DEEQ', numeDofSave//'.NUME.DEEQ'))
!
    call detrsd('NUME_DDL', numeDofSave)
    call infbav()
!
end subroutine
