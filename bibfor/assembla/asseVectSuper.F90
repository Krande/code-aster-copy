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
subroutine asseVectSuper(model, mesh, vectElem, &
                         vectScalType, vectElemCoef, &
                         nomacr, meshNbNode, &
                         nec, nbecmx, nbCmp, &
                         icodla, icodge, &
                         idprn1, idprn2, &
                         iapsdl, ianueq, jvale, jresl)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cordd2.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/ssvalv.h"
!
    character(len=8), intent(in) :: model, mesh
    character(len=19), intent(in) :: vectElem
    integer(kind=8), intent(in) :: vectScalType
    real(kind=8), intent(in) :: vectElemCoef
    character(len=8), pointer :: nomacr(:)
    integer(kind=8), intent(in) :: meshNbNode, nbCmp, nec, nbecmx
    integer(kind=8), intent(in) :: iapsdl, ianueq, jvale, jresl
    integer(kind=8), intent(in) :: idprn1, idprn2
    integer(kind=8), intent(inout) :: icodla(nbecmx), icodge(nbecmx)
!
! --------------------------------------------------------------------------------------------------
!
! Assembly vector from super elements
!
! --------------------------------------------------------------------------------------------------
!
! - Convention: first LIGREL (model) is on mesh
    integer(kind=8), parameter :: ligrelMeshIndx = 1
!
    integer(kind=8) :: iSuperCell, iLoadCase, iNode, iDofSuper
    integer(kind=8) :: nbSuperCell, nbLoadCase, nbNode, nbDofSuper
    integer(kind=8) :: iec, iDof, iret
    integer(kind=8) :: iad1
    integer(kind=8) :: nodeNume, ncmpel, nodeNumeOld
    integer(kind=8), pointer :: sssa(:) => null()
    integer(kind=8), pointer :: relc(:) => null()
    integer(kind=8), pointer :: conx(:) => null()
    integer(kind=8), pointer :: superCellNode(:) => null()
    integer(kind=8), pointer :: prno(:) => null()
    character(len=8) :: loadCaseName, superCellName
    character(len=14) :: superCellNume
    character(len=19) :: ligrel
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    loadCaseName = ' '
    call ssvalv('DEBUT', loadCaseName, model, mesh, 0, jresl, ncmpel)
    call dismoi('NB_SM_MAILLA', model, 'MODELE', repi=nbSuperCell)
    if (nbSuperCell == 0) then
        goto 999
    end if
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
    call jeveuo(ligrel//'.SSSA', 'L', vi=sssa)
    nbLoadCase = 0
    call jeexin(vectElem//'.RELC', iret)
    if (iret .ne. 0) then
        call jelira(vectElem//'.RELC', 'NUTIOC', nbLoadCase)
    end if

! - Loop on load cases
    do iLoadCase = 1, nbLoadCase
        call jenuno(jexnum(vectElem//'.RELC', iLoadCase), loadCaseName)
        call jeveuo(jexnum(vectElem//'.RELC', iLoadCase), 'L', vi=relc)

! ----- Loop on super cells
        do iSuperCell = 1, nbSuperCell
            if (sssa(iSuperCell) .eq. 0) cycle
            if (relc(iSuperCell) .eq. 0) cycle
            call jeveuo(jexnum(mesh//'.SUPMAIL', iSuperCell), 'L', vi=superCellNode)
            call jelira(jexnum(mesh//'.SUPMAIL', iSuperCell), 'LONMAX', nbNode)
            call ssvalv(' ', loadCaseName, model, mesh, iSuperCell, jresl, ncmpel)
            superCellName = nomacr(iSuperCell)
            call dismoi('NOM_NUME_DDL', superCellName, 'MACR_ELEM_STAT', repk=superCellNume)
            call jeveuo(superCellName//'.CONX', 'L', vi=conx)
            call jeveuo(jexnum(superCellNume//'.NUME.PRNO', 1), 'L', vi=prno)
            iDof = 0
            do iNode = 1, nbNode
                nodeNume = superCellNode(iNode)
                if (nodeNume .gt. meshNbNode) then
                    do iec = 1, nbecmx
                        icodge(iec) = icodla(iec)
                    end do
                else
                    nodeNumeOld = conx(3*(iNode-1)+2)
                    do iec = 1, nec
                        icodge(iec) = prno((nec+2)*(nodeNumeOld-1)+2+iec)
                    end do
                end if

! ------------- Get components
                iad1 = zi(idprn1-1+zi(idprn2+ligrelMeshIndx-1)+(nodeNume-1)*(nec+2))
                call cordd2(idprn1, idprn2, ligrelMeshIndx, icodge, nec, &
                            nbCmp, nodeNume, nbDofSuper, zi(iapsdl))

! ------------- Add value in vector (cumul)
                if (vectScalType .eq. 1) then
                    do iDofSuper = 1, nbDofSuper
                        iDof = iDof+1
                        zr(jvale-1+zi(ianueq-1+iad1+zi(iapsdl-1+iDofSuper)-1)) = &
                            zr(jvale-1+zi(ianueq-1+iad1+zi(iapsdl-1+iDofSuper)-1))+ &
                            zr(jresl+iDof-1)*vectElemCoef
                    end do
                else if (vectScalType .eq. 2) then
                    do iDofSuper = 1, nbDofSuper
                        iDof = iDof+1
                        zc(jvale-1+zi(ianueq-1+iad1+zi(iapsdl-1+iDofSuper)-1)) = &
                            zc(jvale-1+zi(ianueq-1+iad1+zi(iapsdl-1+iDofSuper)-1))+ &
                            zc(jresl+iDof-1)*vectElemCoef
                    end do
                else
                    ASSERT(ASTER_FALSE)
                end if
            end do
        end do
    end do
999 continue
!
! - End: clean temporary objects
    call ssvalv('FIN', loadCaseName, model, mesh, 0, jresl, ncmpel)
!
!
    call jedema()
end subroutine
