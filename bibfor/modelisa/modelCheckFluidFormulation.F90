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
subroutine modelCheckFluidFormulation(model)
!
    use model_module, only: getFluidCell
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/teattr.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: model
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_MODELE
!
! Check that we have the same formulation in fluid modelisation
!
! --------------------------------------------------------------------------------------------------
!
! In  model           : name of the model
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv, iret
    character(len=8) :: mesh
    integer(kind=8) :: nbCellFluid
    integer(kind=8), pointer :: cellFluid(:) => null()
    integer(kind=8), pointer :: modelCells(:) => null()
    integer(kind=8) :: iCellFluid, cellTypeNume
    character(len=19) :: ligrel
    character(len=16) :: cellTypeName, FEForm, FEForm2
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
    call jeexin(ligrel//'.TYFE', iret)
    if (iret .ne. 0) then

! ----- Get list of cell in model
        call jeveuo(ligrel//'.TYFE', 'L', vi=modelCells)

! ----- Get list of cells with fluid model (Note that FSI cells are alreary considered
!       as Fluid cell)
        call getFluidCell(model, nbCellFluid, cellFluid)

! ----- Check that all fluid are modeled with the same formulation
        FEForm = ' '
        do iCellFluid = 1, nbCellFluid
            cellTypeNume = modelCells(cellFluid(iCellFluid))
            if (cellTypeNume .ne. 0) then
                call jenuno(jexnum('&CATA.TE.NOMTE', cellTypeNume), cellTypeName)
                call teattr('C', 'FORMULATION', FEForm2, iret, typel=cellTypeName)
                if (iret .ne. 0) then
                    cycle
                end if
                if (FEForm .eq. ' ') then
                    ! Store the type of formulation in the first fluid cell
                    FEForm = FEForm2
                    cycle
                end if
                ! Check that we used the same formulation
                if (FEForm2 .ne. FEForm) then
                    call utmess('F', 'FLUID1_8')
                end if
            end if
        end do
        AS_DEALLOCATE(vi=cellFluid)
    end if
end subroutine
