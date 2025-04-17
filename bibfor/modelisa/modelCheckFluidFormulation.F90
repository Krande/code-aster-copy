! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jeveuo.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/teattr.h"
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
    integer :: ifm, niv, iret, iret2
    character(len=8) :: mesh
    integer :: nbCell, nbCellFluid
    integer, pointer :: cellFluid(:) => null()
    integer, pointer :: modelCells(:)
    integer :: iCellFluid, cellTypeNume
    character(len=16) :: cellTypeName, FEForm, FEForm2
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call jeexin(model//'.MAILLE', iret)
    if (iret .ne. 0) then

! ----- Get list of cell in model
        modelCells => null()
        call jeveuo(model//'.MAILLE', 'L', vi=modelCells)
        
! ----- Get list of cells with fluid model (Note that FSI cells are alreary considered as Fluid cell)
        call getFluidCell(model, nbCellFluid, cellFluid)

! ----- Check that all fluid are modeled with the same formulation 

        if (nbCellFluid .ne. 0) then

! ---------- Get type of formulation in the first fluid cell 
             iCellFluid = 1
             cellTypeNume = modelCells(cellFluid(iCellFluid))
             if (cellTypeNume .ne. 0) then
                 call jenuno(jexnum('&CATA.TE.NOMTE', cellTypeNume), cellTypeName)
             end if
             call teattr('S', 'FORMULATION', FEForm, iret, typel = cellTypeName)
             
! ---------- Get type of formulation in the "iCellFluid" fluid cell 
             do iCellFluid = 2, nbCellFluid
                 cellTypeNume = modelCells(cellFluid(iCellFluid))
                 if (cellTypeNume .ne. 0) then
                     call jenuno(jexnum('&CATA.TE.NOMTE', cellTypeNume), cellTypeName)
                 end if
                 call teattr('S', 'FORMULATION', FEForm2, iret, typel = cellTypeName)

! -------------- Check that we used the same formulation 
                 if (FEForm2 .ne. FEForm) then
                     call utmess('F', 'FLUID1_8')
                 end if 
             end do 
             AS_DEALLOCATE(vi=cellFluid)
             
        end if
        
    end if 
!

end subroutine



