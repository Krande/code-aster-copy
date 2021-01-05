! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine cm_dclac(meshIn, meshOut)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnmpmc.h"
#include "asterfort/copisd.h"
#include "asterfort/cpifpa.h"
#include "asterfort/cppagn.h"
#include "asterfort/def_list_test.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvtx.h"
#include "asterfort/gtgrma.h"
#include "asterfort/jedetr.h"
#include "asterfort/wkvect.h"
#include "asterfort/infniv.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/utmess.h"
!
character(len=8), intent(in) :: meshIn, meshOut
!
! --------------------------------------------------------------------------------------------------
!
! CREA_MAILLAGE - DECOUPE_LAC
!
! --------------------------------------------------------------------------------------------------
!
! In  meshIn          : name of mesh (input)
! In  meshOut         : name of mesh (output)
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ifm, niv
    character(len=24), parameter :: cninv = '&&CPPAGN.CNINV'
    character(len=24) :: typ_dec_lac
    character(len=16), parameter :: keywfact = 'DECOUPE_LAC'
    character(len=16) :: ligrma
    character(len=8) :: meshAux
    integer :: nbSlavCellGroup, iSlavCellGroup, nbtrav, jvConxInv
    integer :: nbCell, nbCellInit, nb_ma_test
    integer :: typ_dec
    integer, pointer :: li_trav(:) => null()
    integer, pointer :: listCell(:) => null()
    character(len=24), pointer :: slavCellGroup(:) => null()
    character(len=24), pointer :: ptrPatchName(:) => null()
    aster_logical :: same_zone
!
! --------------------------------------------------------------------------------------------------
!
    meshAux = 'MAILAUX'
    ligrma  = '&&OP0167.LIMA'
    nbSlavCellGroup   = 0
    call infniv(ifm, niv)
!
! - Get group
!
    call getvtx(keywfact, 'GROUP_MA_ESCL' , iocc=1, nbval=0, nbret=nbSlavCellGroup)
    nbSlavCellGroup = -nbSlavCellGroup
    call wkvect(ligrma, 'V V K24', nbSlavCellGroup, vk24 = slavCellGroup)
    call getvtx(keywfact, 'GROUP_MA_ESCL' , iocc=1, nbval=nbSlavCellGroup, vect=slavCellGroup)
!
! - Get options
!
    call getvtx(keywfact, 'DECOUPE_HEXA', iocc=1,scal=typ_dec_lac)
    if (typ_dec_lac == "PYRA") then
        typ_dec = 1
    elseif (typ_dec_lac == "HEXA") then
        typ_dec = 0
    endif
! 
    call copisd('MAILLAGE', 'V', meshIn, meshAux)
!
    do iSlavCellGroup = 1,nbSlavCellGroup
! ----- LISTE DE MAILLE DU GROUP_MA
        same_zone  = ASTER_FALSE
        nb_ma_test = 0
        nbtrav     = 0
        call gtgrma(meshIn, meshAux, slavCellGroup(iSlavCellGroup), listCell, nbCell)
        nbCellInit = nbCell
        AS_ALLOCATE(vi=li_trav, size=nbCell)
! ----- Gestion du cas avec des mailles surfaciques possédant des mailles volumiques communes
        do while (nb_ma_test .lt. nbCell)

            call wkvect(cninv,'V V I', nbCell, jvConxInv)
            call cnmpmc(meshAux,nbCell, listCell, zi(jvConxInv))
            call def_list_test(nbCell, jvConxInv, listCell, li_trav, nbtrav)
! --------- CREATION DES PATCHS ET RAFFINEMENT LOCAL
            call cppagn(meshAux, meshOut,  nbtrav, li_trav, iSlavCellGroup, typ_dec, jvConxInv,&
                        same_zone, nb_ma_test)
! --------- COPIE DES DONNEES DANS LE MAILLAGE AUXILIAIRE
            call detrsd('MAILLAGE', meshAux)
            call copisd('MAILLAGE', 'V', meshOut, meshAux)
            call cpifpa(meshOut, meshAux)
            same_zone = ASTER_TRUE
            nbtrav = nb_ma_test
            call jedetr(cninv)
            AS_DEALLOCATE(vi=listCell)
            AS_DEALLOCATE(vi=li_trav)
            call gtgrma(meshIn, meshAux, slavCellGroup(iSlavCellGroup), listCell, nbCell)
            AS_ALLOCATE(vi=li_trav, size=nbCell)
        enddo
        AS_DEALLOCATE(vi=listCell)
        AS_DEALLOCATE(vi=li_trav)
! ----- IMPRESSIONS
        if (niv .ge. 1) then
            call utmess('I', 'MESH1_3', ni = 2, vali = [iSlavCellGroup, nbCellInit])
        endif
! ----- NETTOYAGE
        AS_DEALLOCATE(vi=li_trav)
    enddo
! - CREATION DU POINTEUR VERS LES NOMS 
    call wkvect(meshOut//'.PTRNOMPAT', 'G V K24', nbSlavCellGroup, vk24 = ptrPatchName)
    do iSlavCellGroup = 1, nbSlavCellGroup
       ptrPatchName(iSlavCellGroup) = slavCellGroup(iSlavCellGroup)
    end do
    call jedetr(meshAux//'.PATCH')
    call jedetr(meshAux//'.CONOPA')
    call jedetr(meshAux//'.COMAPA')
    call jedetr(meshAux)
    call jedetr(ligrma)
!
end subroutine
