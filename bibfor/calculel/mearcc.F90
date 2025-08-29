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
subroutine mearcc(option, model, chin, chout)
!
    use postComp_module
    use mesh_module
    use model_module
    implicit none
!
#include "asterc/indik8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getelem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=8), intent(in) :: model
    character(len=16) :: option
    character(len=24) :: chin, chout
!
! --------------------------------------------------------------------------------------------------
!
!     BUT: REDUIRE LE CHAMP DE CONTRAINTES (ELNO) DES ELEMENTS 3D
!          A LEURS FACES POUR L'OPTION SIRO_ELEM
!
! --------------------------------------------------------------------------------------------------
!
!     IN  MO     : NOM DU MODELE
!     IN  OPTION : NOM DE L'OPTION
!     IN  CHIN   : CHAMP DE CONTRAINTE ELNO DES ELEMENTS 3D
!     OUT CHOUT  : CHAMP DE CONTRAINTES ELNO REDUIT AUX MAILLES DE PEAU
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = " "
    integer(kind=8), parameter :: iOcc = 0
    integer(kind=8), parameter :: nbCmpMax = 6
    character(len=8), parameter  :: cmpName(nbCmpMax) = (/'SIXX', 'SIYY', 'SIZZ', &
                                                          'SIXY', 'SIXZ', 'SIYZ'/)
    character(len=19), parameter :: chins = '&&MEARCC.CHIN_S', chous = '&&MEARCC.CHOUT_S'
    integer(kind=8) :: ibid, cellDime, iCellSkin, nbCellSkin
    integer(kind=8) :: jcesd3, jcesk3, jcesl3, jcesd2
    integer(kind=8) :: jcesk2, jcesl2, jcesc2, jlcnx, jcnx, ipt, iCmp, nodeSkinNume, nodeSuppNume
    integer(kind=8) :: jco3, jco2, npt3, npt2, ipt2, ipt3, k, npt
    integer(kind=8) :: iad3, iad2, cmpNume, nbCmpEff
    character(len=8) :: mesh, k8b
    integer(kind=8) :: cellNumeSkin, cellNumeSupp
    character(len=24) :: modelLigrel
    integer(kind=8), pointer :: pt3d(:) => null()
    real(kind=8), pointer :: cesv2(:) => null()
    real(kind=8), pointer :: cesv3(:) => null()
    character(len=8), pointer :: cesc3(:) => null()
    integer(kind=8) :: nbCell
    character(len=24), parameter :: listCellJv = "&&MEDOM2.LISTE_MAILLES"
    integer(kind=8), pointer :: listCell(:) => null()
    integer(kind=8) :: nbCellSiro
    integer(kind=8), pointer :: listCellSiro(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to mesh
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('DIM_GEOM', mesh, 'MAILLAGE', repi=cellDime)
    call jeveuo(mesh//'.CONNEX', 'L', jcnx)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', jlcnx)
    nbCmpEff = nbCmpMax
    if (cellDime .eq. 2) then
        nbCmpEff = 4
    end if

! - Get list of cells
    if (hasCellsDefinedFromCmd(factorKeyword, iOcc)) then
! ----- Get list of cells from user
        call getelem(mesh, factorKeyword, iOcc, 'F', listCellJv, nbCell, model=model)
        ASSERT(nbCell .ge. 1)
        call jeveuo(listCellJv, 'L', vi=listCell)
    else
! ----- Get list of all cells in model
        call getAllCellsAffectedByModel(model, nbCell, listCell)
    end if

! - Get cells only for SIRO_ELEM option
    call getCellSiroElem(model, nbCell, listCell, &
                         nbCellSiro, listCellSiro)
    nbCellSkin = nbCellSiro/2

! - Acces to input field
    call celces(chin, 'V', chins)
    call jeveuo(chins//'.CESV', 'L', vr=cesv3)
    call jeveuo(chins//'.CESD', 'L', jcesd3)
    call jeveuo(chins//'.CESK', 'L', jcesk3)
    call jeveuo(chins//'.CESL', 'L', jcesl3)
    call jeveuo(chins//'.CESC', 'L', vk8=cesc3)

! - Create output field on skin elements
    call cescre('V', chous, 'ELNO', mesh, 'SIEF_R', &
                nbCmpEff, cmpName, [-1], [-1], [-nbCmpEff])
    call jeveuo(chous//'.CESV', 'E', vr=cesv2)
    call jeveuo(chous//'.CESD', 'E', jcesd2)
    call jeveuo(chous//'.CESK', 'E', jcesk2)
    call jeveuo(chous//'.CESL', 'E', jcesl2)
    call jeveuo(chous//'.CESC', 'E', jcesc2)
! - Create pairing: for each skin cell which is the support cell paired to ?
    AS_ALLOCATE(vi=pt3d, size=nbCellSkin*MT_NNOMAX2D)
    do iCellSkin = 1, nbCellSkin
        cellNumeSkin = listCellSiro(iCellSkin)
        cellNumeSupp = listCellSiro(nbCellSkin+iCellSkin)
        if (cellNumeSupp .eq. 0) cycle
        jco3 = jcnx+zi(jlcnx-1+cellNumeSupp)-1
        jco2 = jcnx+zi(jlcnx-1+cellNumeSkin)-1
        npt3 = zi(jcesd3-1+5+4*(cellNumeSupp-1)+1)
        npt2 = zi(jcesd2-1+5+4*(cellNumeSkin-1)+1)
        k = 0
        do ipt2 = 1, npt2
            nodeSkinNume = zi(jco2+ipt2-1)
            do ipt3 = 1, npt3
                nodeSuppNume = zi(jco3+ipt3-1)
                if (nodeSuppNume .eq. nodeSkinNume) then
                    k = k+1
                    pt3d(MT_NNOMAX2D*(iCellSkin-1)+k) = ipt3
                    exit
                end if
            end do
        end do
    end do

!   -- remplissage du champ simple 3d
!   ----------------------------------
    do iCellSkin = 1, nbCellSkin
        cellNumeSkin = listCellSiro(iCellSkin)
        cellNumeSupp = listCellSiro(nbCellSkin+iCellSkin)
        if (cellNumeSupp .eq. 0) cycle
        npt = zi(jcesd2-1+5+4*(cellNumeSkin-1)+1)
        do ipt = 1, npt
            do iCmp = 1, nbCmpEff
                cmpNume = indik8(cesc3, cmpName(iCmp), 1, zi(jcesd3+1))
                call cesexi('C', jcesd3, jcesl3, cellNumeSupp, &
                            pt3d(MT_NNOMAX2D*(iCellSkin-1)+ipt), 1, cmpNume, iad3)
                if (iad3 .eq. 0) then
                    call utmess('F', 'CALCULEL5_52', ni=2, vali=[cellNumeSupp, cellNumeSkin])
                end if
                call cesexi('S', jcesd2, jcesl2, cellNumeSkin, ipt, &
                            1, cmpNume, iad2)
                cesv2(1-iad2-1) = cesv3(iad3)
                zl(jcesl2-iad2-1) = .true.
            end do
        end do
    end do
!
    call cesred(chous, nbCellSkin, listCellSiro, 0, [k8b], 'V', chous)
    call cescel(chous, modelLigrel, option, 'PSIG3D', 'OUI', &
                ibid, 'V', chout, 'F', ibid)
!
    AS_DEALLOCATE(vi=pt3d)
    AS_DEALLOCATE(vi=listCellSiro)
    call detrsd('CHAM_ELEM_S', chous)
    call detrsd('CHAM_ELEM_S', chins)
    if (hasCellsDefinedFromCmd(factorKeyword, iOcc)) then
        call jedetr(listCellJv)
    else
        AS_DEALLOCATE(vi=listCell)
    end if
!
    call jedema()
!
end subroutine
