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
subroutine calmaa(modelInterfaceZ, dir, &
                  lchin, lpain, numeDof, matrAsse)
!
    implicit none
!
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assmam.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/memare.h"
#include "asterfort/numddl.h"
#include "asterfort/promor.h"
#include "asterfort/reajre.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: modelInterfaceZ
    character(len=1), intent(in) :: dir
    character(len=24), intent(in) :: lchin(1)
    character(len=8), intent(in) :: lpain(1)
    character(len=14), intent(out) :: numeDof
    character(len=24), intent(out) :: matrAsse
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE CALCULANT LES MATRICES  N(I)*N(J)*NX OU NY ou NZ SUR LE MODELE (THERMIQUE) D'INTERFACE
!
! IN : modelInterface : MODELE INTERFACE
! IN : MATE : MATERIAU CODE
! IN : DIR : DIRECTION X OU Y
! IN : LIGRMO : LIGREL DU MODELE
! IN: LCHIN,LPAIN,LPAOUT : PARAMETRES DU CALCUL ELEMENTAIRE
! IN : NUM : NUMEROTATION DES DDLS THERMIQUES D INTERFACE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: renumSans = "SANS"
    character(len=8), parameter :: lpaout(1) = (/'PMATTTR'/)
    character(len=16) :: option
    character(len=24) :: matrElem, resuElem, modelLigrel
    character(len=24) :: lchout(1)
    integer(kind=8), parameter :: nbMatrElem = 1
    character(len=24), pointer :: listMatrElem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    numeDof = " "
    matrAsse = " "
    option = 'FLUX_FLUI_'//dir

! - Allocate objects for resu_elem
    matrElem = '&&CA.MA'//dir
    call memare('V', matrElem, modelInterfaceZ, option)
    call reajre(matrElem, ' ', "V")

! - Set output field
    resuElem = matrElem(1:8)//'.ME000'
    call codent(1, 'D0', resuElem(12:14))
    lchout(1) = resuElem

! - Compute elementary matrices
    call dismoi("NOM_LIGREL", modelInterfaceZ, "MODELE", repk=modelLigrel)
    call calcul('S', option, modelLigrel, 1, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
    call reajre(matrElem, resuElem, "V")

! - Create list of elementary matrices
    AS_ALLOCATE(vk24=listMatrElem, size=nbMatrElem)
    listMatrElem(1) = matrElem

! - Create numbering from list of elementary matrices
    numeDof = '&&CALMAA.NUM'//dir
    call numddl(numeDof, renumSans, 'VV', nbMatrElem, listMatrElem)
    AS_DEALLOCATE(vk24=listMatrElem)

! - Create storage for matrices
    call promor(numeDof, 'V')

! - Assemblying matrices
    matrAsse = '&&CA.AA'//dir
    call assmam('V', matrAsse, 1, matrElem, [1.d0], &
                numeDof, 'ZERO', 1)
!
    call jedema()
end subroutine
