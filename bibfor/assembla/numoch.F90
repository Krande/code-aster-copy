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
subroutine numoch(listMatrElem, nbMatrElem, listLigr, nbLigr)
!
    implicit none
!
#include "asterfort/dismoi.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8), intent(in) :: nbMatrElem
    character(len=24), pointer :: listMatrElem(:)
    character(len=24), pointer :: listLigr(:)
    integer(kind=8), intent(out) :: nbLigr
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Extract list of LIGREL from elementary matrixes
!
! --------------------------------------------------------------------------------------------------
!
! Ptr listMatrElem   : pointer to list of elementary matrixes
! In  nbMatrElem     : number of elementary matrixes
! Ptr listLigr       : pointer to list of LIGREL
! In  nbLigr         : number of LIGREL in list
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: model
    character(len=19) :: matrElem, resuElem
    character(len=19) :: modelLigrel, ligrelName
    integer(kind=8) :: nbResuLigr, iLigr
    integer(kind=8) :: iMatrElem, iResuElem, iret, nb_subs, nbResuElem
    aster_logical :: l_found
    character(len=24), pointer :: resuElemNoli(:) => null()
    character(len=24), pointer :: listResuElem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nbLigr = 0

! - Count number of LIGREL from resuElem
    nbResuLigr = 2
    do iMatrElem = 1, nbMatrElem
        matrElem = listMatrElem(iMatrElem) (1:19)
        call jeexin(matrElem//'.RELR', iret)
        if (iret .ne. 0) then
            call jelira(matrElem//'.RELR', 'LONUTI', nbResuElem)
            nbResuLigr = nbResuLigr+nbResuElem
        end if
    end do

! - Create object (oversized !)
    AS_ALLOCATE(vk24=listLigr, size=nbResuLigr)

! - Set ligrel in object
    do iMatrElem = 1, nbMatrElem
        matrElem = listMatrElem(iMatrElem) (1:19)

! ----- Substructuration matrix
        call dismoi('NB_SS_ACTI', matrElem, 'MATR_ELEM', repi=nb_subs)
        if (nb_subs .gt. 0) then
! --------- Get LIGREL (from model)
            call dismoi("NOM_MODELE", matrElem, 'MATR_ELEM', repk=model)
            call dismoi("NOM_LIGREL", model, 'MODELE', repk=modelLigrel)

! --------- Already in list ?
            l_found = ASTER_FALSE
            do iLigr = 1, nbLigr
                if (modelLigrel .eq. listLigr(iLigr)) then
                    l_found = ASTER_TRUE
                end if
            end do

! --------- No => add it
            if (.not. l_found) then
                nbLigr = nbLigr+1
                listLigr(nbLigr) = modelLigrel
            end if
        end if

! ----- Standard matrix
        call jeexin(matrElem//'.RELR', iret)
        if (iret .ne. 0) then
            call jeveuo(matrElem//'.RELR', 'L', vk24=listResuElem)
            call jelira(matrElem//'.RELR', 'LONUTI', nbResuElem)
            do iResuElem = 1, nbResuElem
! ------------- Get LIGREL (from RESU_ELEM)
                resuElem = listResuElem(iResuElem) (1:19)
                call jeexin(resuElem//'.NOLI', iret)
                if (iret .ne. 0) then
                    call jeveuo(resuElem//'.NOLI', 'L', vk24=resuElemNoli)
                    ligrelName = resuElemNoli(1) (1:19)
! ----------------- Already in list ?
                    l_found = ASTER_FALSE
                    do iLigr = 1, nbLigr
                        if (ligrelName .eq. listLigr(iLigr)) then
                            l_found = ASTER_TRUE
                        end if
                    end do

! ----------------- No => add it
                    if (.not. l_found) then
                        nbLigr = nbLigr+1
                        listLigr(nbLigr) = ligrelName
                    end if
                end if
            end do
        end if
    end do
!
end subroutine
