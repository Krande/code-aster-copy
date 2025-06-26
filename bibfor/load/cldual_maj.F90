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

subroutine cldual_maj(listLoad, disp)

    implicit none

#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/ischar_iden.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/load_list_info.h"
#include "asterfort/solide_tran_maj.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(in) :: listLoad
    character(len=19), intent(in) :: disp
!
! --------------------------------------------------------------------------------------------------
!
! Loads - Computation
!
! Update dualized relations for non-linear Dirichlet boundary conditions (undead)
!
! --------------------------------------------------------------------------------------------------
!
! In  listLoad         : name of datastructure for list of loads
! In  disp             : displacements
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: dual_type
    character(len=8) :: loadName
    character(len=8), pointer :: dual_prdk(:) => null()
    character(len=8), pointer :: loadType(:) => null()
    aster_logical :: ltran
    aster_logical :: load_empty
    integer(kind=8) :: nb_link, i_link, iexi
    integer(kind=8) :: iLoad, nbLoad
    integer(kind=8) :: i_load_diri
    aster_logical :: ischar_diri
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    character(len=24), pointer :: listLoadName(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    i_load_diri = 0
    ltran = .false._1
    ischar_diri = .false._1

! - Loads
    call load_list_info(load_empty, nbLoad, listLoadName, listLoadInfo, &
                        list_load_=listLoad)
    ASSERT(.not. load_empty)

! - Identify undead Dirichlet load
    do iLoad = 1, nbLoad
        ischar_diri = ischar_iden(listLoadInfo, iLoad, nbLoad, 'DIRI', 'SUIV', listLoadName(iLoad))
        if (ischar_diri) then
            loadName = listLoadName(iLoad) (1:8)

! --------- Some checks
            call jeexin(loadName//'.DUAL.PRDK', iexi)
            ASSERT(iexi .gt. 0)
            call jeveuo(loadName//'.TYPE', 'L', vk8=loadType)
            ASSERT(loadType(1) .eq. 'MECA_RE')

! --------- Datastructure access
            call jeveuo(loadName//'.DUAL.PRDK', 'L', vk8=dual_prdk)
            call jelira(loadName//'.DUAL.PRDK', 'LONUTI', ival=nb_link)
!
! --------- Find type of dual relation
!
            do i_link = 1, nb_link
                dual_type = dual_prdk(i_link)
                if (dual_type .eq. ' ' .or. dual_type .eq. 'LIN') then
! -------- Nothing to do
                else if (dual_type(1:2) .eq. '2D' .or. dual_type(1:2) .eq. '3D') then
                    ltran = .true._1
                else if (dual_type .eq. 'ROTA2D') then
                    call utmess('F', 'CHARGES_35')
                else if (dual_type .eq. 'ROTA3D') then
                    call utmess('F', 'CHARGES_36')
                else if (dual_type .eq. 'NLIN') then
!                  -- Le nom de la charge est malheureusement inexploitable.
!                     C'est une copie temporaire.
                    call utmess('F', 'CHARGES_30')
                else
                    ASSERT(ASTER_FALSE)
                end if
            end do
!
! --------- Update for LIAISON_SOLIDE
!
            if (ltran) then
                call solide_tran_maj(loadName, disp)
            end if
        end if
    end do
!
    call jedema()
end subroutine
