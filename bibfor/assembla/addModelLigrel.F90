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
subroutine addModelLigrel(modelZ, nbLigr, listLigr)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/dismoi.h"
!
    character(len=*), intent(in) :: modelZ
    integer(kind=8), intent(inout) :: nbLigr
    character(len=24), pointer :: listLigr(:)
!
! --------------------------------------------------------------------------------------------------
!
! Add LIGREL from model
!
! --------------------------------------------------------------------------------------------------
!
! In  model             : name of model datastructure
! IO  nbLigr            : number of LIGREL in list
! Ptr listLigr          : list of LIGREL
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: modelLigrel
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi("NOM_LIGREL", modelZ, "MODELE", repk=modelLigrel)

! - On met toujours le modèle en premier dans la liste
! - Affreuse glute (tant qu'on stockera le modèle dans NUME_EQUA/REFN) (voir numero.F90)
    ASSERT(nbLigr .eq. 0)
    AS_ALLOCATE(vk24=listLigr, size=1)
    nbLigr = 1
    listLigr(nbLigr) = modelLigrel
!
end subroutine
