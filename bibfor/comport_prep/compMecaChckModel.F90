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
subroutine compMecaChckModel(iComp, &
                             model, fullElemField, &
                             lAllCellAffe, cellAffe, nbCellAffe, &
                             relaComp, relaCompPY, chmate, typeComp, &
                             lElasByDefault, lNeedDeborst, lIncoUpo)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/lctest.h"
#include "asterfort/asmpi_all.h"
#include "asterfort/asmpi_any.h"
#include "asterfort/cesexi.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/dismoi.h"
#include "asterfort/dismte.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/teattr.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: iComp
    character(len=8), intent(in) :: model
    character(len=19), intent(in) :: fullElemField
    aster_logical, intent(in) :: lAllCellAffe
    character(len=24), intent(in) :: cellAffe
    integer(kind=8), intent(in) :: nbCellAffe
    character(len=16), intent(in) :: relaCompPY, relaComp
    character(len=16), intent(in) :: typeComp
    character(len=8), intent(in) :: chmate
    aster_logical, intent(out) :: lElasByDefault, lNeedDeborst, lIncoUpo
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Checking the consistency of the modelization with the behavior
!
! --------------------------------------------------------------------------------------------------
!
! In  iComp            : occurrence number
! In  model            : name of model
! In  fullElemField    : field for FULL_MECA option
! In  lAllCellAffe     : ASTER_TRUE if affect on all cells where behaviour is defined
! In  nbCellAffe       : number of cells where behaviour is defined
! In  cellAffe         : list of cells where behaviour is defined
! In  relaComp         : comportement RELATION
! In  relaCompPY       : comportement RELATION - Python coding
! In  chmate           : material field (sd_mater)
! Out lElasByDefault   : flag if at least one element use ELAS by default
! Out lNeedDeborst     : flag if at least one element swap to Deborst algorithm
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: elemTypeName, modelType, incoType, isNuFunc, typmod2Type
    character(len=8) :: mesh
    integer(kind=8) :: elemTypeNume, cellNume, nbCmpAffected
    integer(kind=8) :: jvCesd, jvCesl, jvVale
    integer(kind=8) :: modelTypeIret, lctestIret, iCell, incoTypeIret
    integer(kind=8) :: nbCellMesh, nbCell, ibid, ier
    character(len=16), pointer :: cesv(:) => null()
    character(len=19) :: ligrel
    integer(kind=8), pointer :: cellAffectedByModel(:) => null()
    integer(kind=8), pointer :: listCellAffe(:) => null()
    aster_logical :: lAtOneCellAffect, lAllCellAreBound, lPlStressFuncNu, l_kit_thm, l_parallel_mesh
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    lElasByDefault = ASTER_FALSE
    lNeedDeborst = ASTER_FALSE
    lAtOneCellAffect = ASTER_FALSE
    lAllCellAreBound = ASTER_FALSE
    lIncoUpo = ASTER_FALSE
    lPlStressFuncNu = ASTER_FALSE
    call comp_meca_l(relaComp, 'KIT_THM', l_kit_thm)

    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    l_parallel_mesh = isParallelMesh(mesh)

    if (nbCellAffe == 0 .and. l_parallel_mesh) goto 999

!
! - Access to model
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
    call jeveuo(ligrel//'.TYFE', 'L', vi=cellAffectedByModel)
!
! - Access to <CHELEM_S> of FULL_MECA option
!
    call jeveuo(fullElemField//'.CESD', 'L', jvCesd)
    call jeveuo(fullElemField//'.CESL', 'L', jvCesl)
    call jeveuo(fullElemField//'.CESV', 'L', vk16=cesv)
    nbCellMesh = zi(jvCesd-1+1)
!
! - Mesh affectation
!
    if (lAllCellAffe) then
        nbCell = nbCellMesh
    else
        call jeveuo(cellAffe, 'L', vi=listCellAffe)
        nbCell = nbCellAffe
    end if
!
! - Is Poisson coefficient given by a function
!
    call dismoi('NU_FO', chmate, 'CHAM_MATER', repk=isNuFunc)
!
! - Loop on elements
!
    nbCmpAffected = 0
    do iCell = 1, nbCell
! ----- Current cell
        if (lAllCellAffe) then
            cellNume = iCell
        else
            cellNume = listCellAffe(iCell)
        end if

! ----- Number of components affected
        nbCmpAffected = max(nbCmpAffected, zi(jvCesd-1+5+4*(cellNume-1)+3))

! ----- Get adress in field for FULL_MECA option
        call cesexi('C', jvCesd, jvCesl, cellNume, 1, 1, 1, jvVale)

        if (jvVale .gt. 0) then
            lAtOneCellAffect = ASTER_TRUE

! --------- Behaviour on this cell is elastic (default)
            if (cesv(jvVale) .eq. ' ') then
                lElasByDefault = ASTER_TRUE
            end if

! --------- Access to type of finite element
            elemTypeNume = cellAffectedByModel(cellNume)
            call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)

! --------- Type of modelization
            call teattr('C', 'TYPMOD', modelType, modelTypeIret, typel=elemTypeName)
            if (modelTypeIret .eq. 0) then
                if (modelType .eq. 'C_PLAN') then
                    call lctest(relaCompPY, 'MODELISATION', 'C_PLAN', lctestIret)
! ----------------- C_PLAN is not allowed for this behaviour => activation of Deborst algorithm
                    if (lctestIret .eq. 0) then
                        lNeedDeborst = ASTER_TRUE
                    end if

! ----------------- bug: NU provided as a function while native plane stress and COMP_INCR

                    if (lctestIret /= 0) then
                        if (isNuFunc(1:3) == 'OUI' .and. typeComp == 'COMP_INCR') then
                            lPlStressFuncNu = ASTER_TRUE
                        end if
                    end if

                else if (modelType .eq. '1D') then
                    call lctest(relaCompPY, 'MODELISATION', '1D', lctestIret)
                    if (lctestIret .eq. 0) then
                        call utmess('F', 'COMPOR4_32', sk=relaComp)
                    end if

                else if (modelType .eq. '3D') then
                    call lctest(relaCompPY, 'MODELISATION', '3D', lctestIret)

                else
                    call lctest(relaCompPY, 'MODELISATION', modelType, lctestIret)

                end if
            end if
! --------- Check presence of modelization INCO_UPO
            call teattr('C', 'INCO', incoType, incoTypeIret, typel=elemTypeName)
            if (incoTypeIret .eq. 0) then
                if (incoType .eq. 'C2O') then
                    lIncoUpo = ASTER_TRUE
                end if
            end if
! --------- Verification pour KIT_THM
            if (l_kit_thm) then
                call teattr('C', 'TYPMOD2', typmod2Type, modelTypeIret, typel=elemTypeName)
                if (typmod2Type .ne. 'THM' .and. typmod2Type .ne. 'XFEM_HM' &
                    .and. typmod2Type .ne. 'JHMS') then
                    call dismte('MODELISATION', elemTypeName, ibid, modelType, ier)
                    call utmess('F', 'COMPOR1_22', nk=2, valk=[relaComp, modelType])
                end if
            end if
        end if
    end do
!
! - All elements are boundary elements
!
    lAllCellAreBound = nbCmpAffected .eq. 0

999 continue

!
! - Comm for MPI
!
    lAllCellAreBound = asmpi_all(lAllCellAreBound, ASTER_TRUE)
    lAtOneCellAffect = asmpi_any(lAtOneCellAffect, ASTER_TRUE)
    lNeedDeborst = asmpi_any(lNeedDeborst, ASTER_TRUE)
    lElasByDefault = asmpi_any(lElasByDefault, ASTER_TRUE)
    lPlStressFuncNu = asmpi_any(lPlStressFuncNu, ASTER_TRUE)
!
! - Error when nothing is affected by the behavior
!
    if (.not. lAtOneCellAffect .and. .not. l_parallel_mesh) then
        if (lAllCellAreBound) then
            call utmess('F', 'COMPOR1_60', si=iComp)
        else
            call utmess('F', 'COMPOR1_59', si=iComp)
        end if
    end if
!
! - Alarm while nu is a function while native plane stress and COMP_INCR
!
    if (lPlStressFuncNu) then
        call utmess('A', 'COMPOR6_14', si=iComp)
    end if
!
    call jedema()
!
end subroutine
