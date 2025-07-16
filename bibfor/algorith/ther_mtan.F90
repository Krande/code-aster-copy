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
subroutine ther_mtan(l_stat, &
                     modelZ, caraElemZ, matecoZ, &
                     timePara, varcCurrZ, &
                     comporTherZ, tempIterZ, &
                     resuElemZ, matrElemZ, jvBase)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/multResuElem.h"
#include "asterfort/reajre.h"
!
    aster_logical, intent(in) :: l_stat
    character(len=*), intent(in) :: modelZ, caraElemZ, matecoZ
    real(kind=8), intent(in) :: timePara(2)
    character(len=*), intent(in) :: tempIterZ, comporTherZ, varcCurrZ
    character(len=*), intent(inout) :: resuElemZ
    character(len=*), intent(in) :: matrElemZ
    character(len=1), intent(in) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
! Thermic
!
! Tangent matrix (volumic terms)
!
! --------------------------------------------------------------------------------------------------
!
! In  l_stat           : flag for stationnary computation (no mass term)
! In  model            : name of the model
! In  caraElem         : name of elementary characteristics (field)
! In  mateco           : name of codeing material characteristics (field)
! In  timePara         : timePara(1) = theta
!                        timePara(2) = deltat
! In  varcCurr         : command variable for current time
! In  comporTher       : name of comportment definition (field)
! In  tempIter         : temperature field at current Newton iteration
! IO  resuElem         : name of resu_elem
! In  matrElem         : name of matr_elem result
! In  jvBase           : JEVEUX base for object
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: optionRigi = 'RIGI_THER_TANG', optionMass = 'MASS_THER_TANG'
    integer(kind=8), parameter :: nbIn = 6, nbout = 1
    character(len=8) :: lpain(nbIn), lpaout(nbout)
    character(len=24) :: lchin(nbIn), lchout(nbout)
    character(len=24) :: ligrel_model
    character(len=24) :: chgeom, chcara(18)
    character(len=19) :: resuElem
    real(kind=8) :: theta, deltat
    character(len=8) :: newnom
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=ligrel_model)
    theta = timePara(1)
    deltat = timePara(2)
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "
    resuElem = resuElemZ(1:19)

! - Geometry field
    call megeom(modelZ, chgeom)

! - Elementary characteristics field
    call mecara(caraElemZ, chcara)

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = matecoZ
    lpain(3) = 'PTEMPEI'
    lchin(3) = tempIterZ
    lpain(4) = 'PCOMPOR'
    lchin(4) = comporTherZ
    lpain(5) = 'PVARCPR'
    lchin(5) = varcCurrZ
    lpain(6) = 'PCAMASS'
    lchin(6) = chcara(12)

! - Output fields
    lpaout(1) = 'PMATTSR'
    lchout(1) = resuElemZ

! - Compute rigidity term
    call calcul("S", optionRigi, ligrel_model, nbin, lchin, &
                lpain, nbout, lchout, lpaout, jvBase, &
                'OUI')

! - Multiply values by theta
    call multResuElem(resuElem, theta)

! - Add RESU_ELEM in MATR_ELEM
    call reajre(matrElemZ, resuElem, jvBase)

! - Compute mass term
    if (.not. l_stat) then
! - --- Output fields
        newnom = resuElem(9:16)
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        lpaout(1) = 'PMATTTR'
        lchout(1) = resuElem

! - --- Compute
        call calcul("S", optionMass, ligrel_model, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, jvBase, &
                    'OUI')

! - --- Multiply values by 1/dt
        call multResuElem(resuElem, 1.d0/deltat)

! - --- Add RESU_ELEM in MATR_ELEM
        call reajre(matrElemZ, resuElem, jvBase)

    end if
!
    resuElemZ = resuElem
!
end subroutine
