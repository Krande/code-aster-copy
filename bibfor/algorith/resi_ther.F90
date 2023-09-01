! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine resi_ther(model, cara_elem, mate, time, compor, &
                     temp_prev, temp_iter, hydr_prev, hydr_curr, &
                     dry_curr, varc_curr, resu_elem, vect_elem, base, &
                     l_stat, para)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/gcnco2.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/multResuElem.h"
#include "asterfort/reajre.h"
#include "asterfort/inical.h"
!
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: cara_elem
    character(len=24), intent(in) :: time
    character(len=24), intent(in) :: mate
    character(len=24), intent(in) :: temp_prev
    character(len=24), intent(in) :: temp_iter
    character(len=24), intent(in) :: hydr_prev
    character(len=24), intent(in) :: hydr_curr
    character(len=24), intent(in) :: dry_curr
    character(len=24), intent(in) :: compor
    character(len=19), intent(in) :: varc_curr
    character(len=19), intent(inout) :: resu_elem
    character(len=24), intent(in) :: vect_elem
    character(len=1), intent(in) :: base
    aster_logical, intent(in) :: l_stat
    real(kind=8), intent(in) :: para(2)

!
! --------------------------------------------------------------------------------------------------
!
! Thermic
!
! Residuals from non-linear laws
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  cara_elem        : name of elementary characteristics (field)
! In  mate             : name of material characteristics (field)
! In  time             : time (<CARTE>)
! In  temp_prev        : previous temperature
! In  temp_iter        : incrementaltemperature field at current Newton iteration
! In  hydr_prev        : previous hydratation
! In  hydr_curr        : current hydratation
! In  dry_curr         : current drying
! In  compor           : name of comportment definition (field)
! In  varc_curr        : command variable for current time
! In  resu_elem        : name of resu_elem
! In  vect_elem        : name of vect_elem result
! In  base             : JEVEUX base for object
! In  para             : para(1) = theta
!                        para(2) = deltat
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbin = 10
    integer, parameter :: nbout = 2
    character(len=8) :: lpain(nbin), lpaout(nbout), newnom
    character(len=19) :: lchin(nbin), lchout(nbout)
!
    character(len=1) :: stop_calc
    character(len=16) :: option1, option2
    character(len=24) :: ligrel_model
    character(len=24) :: chgeom, chcara(18)
    real(kind=8) :: theta, deltat
!
! --------------------------------------------------------------------------------------------------
!
    stop_calc = 'S'
    option1 = 'RAPH_THER'
    option2 = 'MASS_THER_RESI'
    ligrel_model = model(1:8)//'.MODELE'
    theta = para(1)
    deltat = para(2)
!
! - Init fields
!
    call inical(nbin, lpain, lchin, nbout, lpaout, &
                lchout)
!
! - Geometry field
!
    call megeom(model, chgeom)
!
! - Elementary characteristics field
!
    call mecara(cara_elem, chcara)
!
! - Input fields
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = mate(1:19)
    lpain(3) = 'PTEMPEI'
    lchin(3) = temp_iter(1:19)
    lpain(4) = 'PCOMPOR'
    lchin(4) = compor(1:19)
    lpain(5) = 'PTMPCHF'
    lchin(5) = dry_curr(1:19)
    lpain(6) = 'PVARCPR'
    lchin(6) = varc_curr(1:19)
    lpain(7) = 'PCAMASS'
    lchin(7) = chcara(12) (1:19)
!
! - Rigidity
!

!
! - Output fields
!
    lpaout(1) = 'PRESIDU'
    lchout(1) = resu_elem(1:19)
    lpaout(2) = 'PFLUXPR'
    lchout(2) = "&&RESI_THER.FLUXPR"
!
    call corich('E', lchout(1), ichin_=-1)
!
! - Number of fields
!
    call calcul(stop_calc, option1, ligrel_model, 7, lchin, &
                lpain, 2, lchout, lpaout, base, &
                'OUI')
!
! - Multiply values by theta
!
    call multResuElem(resu_elem, theta)
!
! - Add RESU_ELEM in vect_elem
!
    call reajre(vect_elem, resu_elem, base)

    if (.not. l_stat) then
!
! --- Compute hydratation
!
        lpain(1) = 'PMATERC'
        lchin(1) = mate(1:19)
        lpain(2) = 'PCOMPOR'
        lchin(2) = compor(1:19)
        lpain(3) = 'PTEMPSR'
        lchin(3) = time(1:19)
        lpain(4) = 'PTEMPMR'
        lchin(4) = temp_prev(1:19)
        lpain(5) = 'PTEMPPR'
        lchin(5) = temp_iter(1:19)
        lpain(6) = 'PHYDRMR'
        lchin(6) = hydr_prev(1:19)

        lpaout(1) = 'PHYDRPR'
        lchout(1) = hydr_curr(1:19)

        call calcul(stop_calc, "HYDR_ELGA", ligrel_model, 6, lchin, &
                    lpain, 1, lchout, lpaout, base, 'OUI')
!
! - --- Mass
!
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom(1:19)
        lpain(2) = 'PMATERC'
        lchin(2) = mate(1:19)
        lpain(3) = 'PTEMPEI'
        lchin(3) = temp_iter(1:19)
        lpain(4) = 'PCOMPOR'
        lchin(4) = compor(1:19)
        lpain(5) = 'PVARCPR'
        lchin(5) = varc_curr(1:19)
        lpain(6) = 'PHYDRPR'
        lchin(6) = hydr_curr(1:19)
        lpain(7) = 'PTEMPSR'
        lchin(7) = time(1:19)
!
! - --- Output fields
!
        newnom = resu_elem(9:16)
        call gcnco2(newnom)
        resu_elem(10:16) = newnom(2:8)

        lpaout(1) = 'PRESIDU'
        lchout(1) = resu_elem(1:19)
!
        call corich('E', lchout(1), ichin_=-1)
!
! - --- Number of fields
!
        call calcul(stop_calc, option2, ligrel_model, 7, lchin, &
                    lpain, 1, lchout, lpaout, base, &
                    'OUI')
!
! - --- Multiply values by 1/dt
!
        call multResuElem(resu_elem, 1.d0/deltat)
!
! - --- Add RESU_ELEM in vect_elem
!
        call reajre(vect_elem, resu_elem, base)
    end if
!
end subroutine
