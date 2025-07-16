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
subroutine imvari(compor_info)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(in) :: compor_info
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Print informations about internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  compor_info      : name of object for information about internal variables and comportement
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iVari, mapZoneNume
    integer(kind=8) :: nbVari, mapNbZone, nbElemZone, nt_vari, nbCellPMF
    character(len=16) :: vari_excl
    character(len=16) :: rela_comp, defo_comp, type_cpla, regu_visc, post_incr
    aster_logical :: l_excl, lPMF
    integer(kind=8), pointer :: comporInfoInfo(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    character(len=16), pointer :: comporInfoVari(:) => null()
    character(len=16), pointer :: comporInfoRela(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lPMF = ASTER_FALSE
    nbCellPMF = 0

! - Access to informations
    call jeveuo(compor_info(1:19)//'.INFO', 'L', vi=comporInfoInfo)
    nt_vari = comporInfoInfo(4)
    if (nt_vari .eq. 0) then
        goto 99
    end if
    call utmess('I', 'COMPOR4_1')
    mapNbZone = comporInfoInfo(2)
    call jeveuo(compor_info(1:19)//'.RELA', 'L', vk16=comporInfoRela)
    call jeveuo(compor_info(1:19)//'.ZONE', 'L', vi=comporInfoZone)

    do mapZoneNume = 1, mapNbZone

        nbElemZone = comporInfoZone(mapZoneNume)
!
        if (nbElemZone .ne. 0) then
! --------- Acces to list of name of internal variables
            call jeveuo(jexnum(compor_info(1:19)//'.VARI', mapZoneNume), 'L', vk16=comporInfoVari)
            call jelira(jexnum(compor_info(1:19)//'.VARI', mapZoneNume), 'LONMAX', nbVari)

! --------- Exceptions ?
            vari_excl = comporInfoVari(1)
            l_excl = vari_excl(1:2) .eq. '&&'

! --------- Get names of relation
            rela_comp = comporInfoRela(5*(mapZoneNume-1)+1)
            defo_comp = comporInfoRela(5*(mapZoneNume-1)+2)
            type_cpla = comporInfoRela(5*(mapZoneNume-1)+3)
            regu_visc = comporInfoRela(5*(mapZoneNume-1)+4)
            post_incr = comporInfoRela(5*(mapZoneNume-1)+5)

! --------- Print name of internal variables
            if (l_excl) then
                if (vari_excl .eq. '&&MULT_COMP') then
                    call utmess('I', 'COMPOR4_4', si=nbElemZone)
                    call utmess('I', 'COMPOR4_11')
                    if (regu_visc .eq. 'VIDE') then
                        call utmess('I', 'COMPOR4_18')
                    else
                        call utmess('I', 'COMPOR4_7', sk=regu_visc)
                    end if
                    if (post_incr .eq. 'VIDE') then
                        call utmess('I', 'COMPOR4_19')
                    else
                        call utmess('I', 'COMPOR4_27', sk=post_incr)
                    end if
                    call utmess('I', 'COMPOR4_9', si=nbVari)
                    call utmess('I', 'COMPOR4_15')
                else if (vari_excl .eq. '&&PROT_COMP') then
                    call utmess('I', 'COMPOR4_4', si=nbElemZone)
                    call utmess('I', 'COMPOR4_10')
                    call utmess('I', 'COMPOR4_6', sk=defo_comp)
                    if (regu_visc .eq. 'VIDE') then
                        call utmess('I', 'COMPOR4_18')
                    else
                        call utmess('I', 'COMPOR4_7', sk=regu_visc)
                    end if
                    if (post_incr .eq. 'VIDE') then
                        call utmess('I', 'COMPOR4_19')
                    else
                        call utmess('I', 'COMPOR4_27', sk=post_incr)
                    end if
                    call utmess('I', 'COMPOR4_9', si=nbVari)
                    call utmess('I', 'COMPOR4_16')
                else if (vari_excl .eq. '&&MULT_PMF') then
                    lPMF = ASTER_TRUE
                    nbCellPMF = nbCellPMF+nbElemZone

                else
                    ASSERT(ASTER_FALSE)
                end if
            else
                call utmess('I', 'COMPOR4_4', si=nbElemZone)
                call utmess('I', 'COMPOR4_5', sk=rela_comp)
                call utmess('I', 'COMPOR4_6', sk=defo_comp)
                if (regu_visc .eq. 'VIDE') then
                    call utmess('I', 'COMPOR4_18')
                else
                    call utmess('I', 'COMPOR4_7', sk=regu_visc)
                end if
                if (post_incr .eq. 'VIDE') then
                    call utmess('I', 'COMPOR4_19')
                else
                    call utmess('I', 'COMPOR4_27', sk=post_incr)
                end if
                if (type_cpla .eq. 'DEBORST') then
                    call utmess('I', 'COMPOR4_8')
                end if
                call utmess('I', 'COMPOR4_9', si=nbVari)
                do iVari = 1, nbVari
                    call utmess('I', 'COMPOR4_20', sk=comporInfoVari(iVari), si=iVari)
                end do
            end if
        end if
    end do
!
99  continue
!
    if (lPMF) then
        call utmess('I', 'COMPOR4_12', si=nbCellPMF)
    end if
!
    call jedema()
!
end subroutine
