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
subroutine resi_ther(l_stat, &
                     modelZ, caraElemZ, matecoZ, &
                     timePara, timeMapZ, varcCurrZ, &
                     comporTherZ, tempIterZ, &
                     tempPrevZ, hydrPrevZ, hydrCurrZ, &
                     resuElemZ, vectElemZ, jvBase)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
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
    character(len=*), intent(in) :: tempPrevZ, hydrPrevZ, hydrCurrZ, timeMapZ
    character(len=*), intent(inout) :: resuElemZ
    character(len=*), intent(in) :: vectElemZ
    character(len=1), intent(in) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
! Thermic
!
! Residuals from non-linear laws
!
! --------------------------------------------------------------------------------------------------
!
! In  l_stat           : flag for stationnary computation (no mass term)
! In  model            : name of the model
! In  caraElem         : name of elementary characteristics (field)
! In  mateco           : name of coding material characteristics (field)
! In  timePara         : timePara(1) = theta
!                        timePara(2) = deltat
! In  varcCurr         : command variable for current time
! In  comporTher       : name of comportment definition (field)
! In  tempIter         : temperature field at current Newton iteration
! In  tempPrev         : previous temperature
! In  hydrPrev         : previous hydration
! In  hydrCurr         : current hydration
! IO  resuElem         : name of resu_elem
! In  vectElem         : name of vect_elem result
! In  jvBase           : JEVEUX base for object
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: optionRigi = 'RAPH_THER', optionMass = 'MASS_THER_RESI'
    character(len=16), parameter :: optionHydr = "HYDR_ELGA"
    integer(kind=8), parameter :: nbIn = 10, nbout = 2
    character(len=8) :: lpain(nbIn), lpaout(nbout)
    character(len=24) :: lchin(nbIn), lchout(nbout)
    character(len=24) :: ligrel_model
    character(len=24) :: chgeom, chcara(18)
    character(len=19) :: resuElem
    real(kind=8) :: theta, deltat
    character(len=8) :: newnom
    character(len=3) :: answer
    aster_logical :: l_dry
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=ligrel_model)
    call dismoi('EXI_SECH', modelZ, 'MODELE', repk=answer)
    l_dry = ASTER_FALSE
    if (answer .eq. 'OUI') l_dry = ASTER_TRUE
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
    lpaout(1) = 'PRESIDU'
    lchout(1) = resuElemZ
    lpaout(2) = 'PFLUXPR'
    lchout(2) = "&&RESI_THER.FLUXPR"
    call corich('E', lchout(1), ichin_=-1)

! - Compute rigidity term
    call calcul("S", optionRigi, ligrel_model, &
                nbin, lchin, lpain, &
                nbout, lchout, lpaout, &
                jvBase, 'OUI')

! - Multiply values by theta
    call multResuElem(resuElem, theta)

! - Add RESU_ELEM in VECT_ELEM
    call reajre(vectElemZ, resuElem, jvBase)

! - Compute hydratation
    if (.not. l_stat .and. .not. l_dry) then
! ----- Input fields
        lpain = " "
        lchin = " "
        lpain(1) = 'PMATERC'
        lchin(1) = matecoZ
        lpain(2) = 'PCOMPOR'
        lchin(2) = comporTherZ
        lpain(3) = 'PINSTR'
        lchin(3) = timeMapZ
        lpain(4) = 'PTEMPMR'
        lchin(4) = tempPrevZ
        lpain(5) = 'PTEMPPR'
        lchin(5) = tempIterZ
        lpain(6) = 'PHYDRMR'
        lchin(6) = hydrPrevZ
        lpain(7) = 'PGEOMER'
        lchin(7) = chgeom

! - --- Output fields
        lpaout(1) = 'PHYDRPR'
        lchout(1) = hydrCurrZ

! - --- Compute
        call calcul("S", optionHydr, ligrel_model, &
                    nbin, lchin, lpain, &
                    1, lchout, lpaout, &
                    jvBase, 'OUI')
    end if

! - Compute mass term
    if (.not. l_stat) then
! ----- Input fields
        lpain = " "
        lchin = " "
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
        lpain(6) = 'PHYDRPR'
        lchin(6) = hydrCurrZ

! - --- Output fields
        newnom = resuElem(9:16)
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        lpaout(1) = 'PRESIDU'
        lchout(1) = resuElem
        call corich('E', lchout(1), ichin_=-1)

! - --- Compute
        call calcul("S", optionMass, ligrel_model, &
                    nbin, lchin, lpain, &
                    1, lchout, lpaout, &
                    jvBase, 'OUI')

! - --- Multiply values by 1/dt
        call multResuElem(resuElem, 1.d0/deltat)

! - --- Add RESU_ELEM in VECT_ELEM
        call reajre(vectElemZ, resuElem, jvBase)
    end if
!
    resuElemZ = resuElem
!
end subroutine
