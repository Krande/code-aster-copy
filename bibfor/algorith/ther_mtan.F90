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
subroutine ther_mtan(model, cara_elem, mate, para, varc_curr, &
                     compor, temp_iter, dry_prev, dry_curr, resu_elem, &
                     matr_elem, base, l_stat)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/gcnco2.h"
#include "asterfort/megeom.h"
#include "asterfort/mecara.h"
#include "asterfort/multResuElem.h"
#include "asterfort/reajre.h"
#include "asterfort/inical.h"
!
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: cara_elem
    real(kind=8), intent(in) :: para(2)
    character(len=24), intent(in) :: mate
    character(len=24), intent(in) :: temp_iter
    character(len=24), intent(in) :: dry_prev
    character(len=24), intent(in) :: dry_curr
    character(len=24), intent(in) :: compor
    character(len=19), intent(in) :: varc_curr
    character(len=19), intent(inout) :: resu_elem
    character(len=24), intent(in) :: matr_elem
    character(len=1), intent(in) :: base
    aster_logical, intent(in) :: l_stat
!
! --------------------------------------------------------------------------------------------------
!
! Thermic
!
! Tangent matrix (volumic terms)
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  cara_elem        : name of elementary characteristics (field)
! In  mate             : name of material characteristics (field)
! In  para             : para(1) = theta
!                        para(2) = deltat
! In  varc_curr        : command variable for current time
! In  compor           : name of comportment definition (field)
! In  temp_iter        : temperature field at current Newton iteration
! In  dry_prev         : previous drying
! In  dry_curr         : current drying
! In  resu_elem        : name of resu_elem
! In  matr_elem        : name of matr_elem result
! In  base             : JEVEUX base for object
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_in_maxi, nbout
    parameter(nb_in_maxi=9, nbout=1)
    character(len=8) :: lpain(nb_in_maxi), lpaout(nbout), newnom
    character(len=19) :: lchin(nb_in_maxi), lchout(nbout)
!
    character(len=1) :: stop_calc
    character(len=16) :: option1, option2
    character(len=24) :: ligrel_model
    character(len=24) :: chgeom, chcara(18)
    integer :: nbin
    real(kind=8) :: theta, deltat

!
! --------------------------------------------------------------------------------------------------
!
    stop_calc = 'S'
    option1 = 'RIGI_THER_TANG'
    option2 = 'MASS_THER_TANG'
    ligrel_model = model(1:8)//'.MODELE'
    theta = para(1)
    deltat = para(2)
!
! - Init fields
!
    call inical(nb_in_maxi, lpain, lchin, nbout, lpaout, &
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
    nbin = 7
!
! - Rigidity term
!

!
! - Output fields
!
    lpaout(1) = 'PMATTTR'
    lchout(1) = resu_elem
!
! - Compute
!
    call calcul(stop_calc, option1, ligrel_model, nbin, lchin, &
                lpain, nbout, lchout, lpaout, base, &
                'OUI')
!
! - Multiply values by theta
!
    call multResuElem(resu_elem, theta)
!
! - Add RESU_ELEM in MATR_ELEM
!
    call reajre(matr_elem, resu_elem, base)
!
! - Mass term
!
    if (.not. l_stat) then
!
! - --- Output fields
!
        newnom = resu_elem(9:16)
        call gcnco2(newnom)
        resu_elem(10:16) = newnom(2:8)

        lpaout(1) = 'PMATTTR'
        lchout(1) = resu_elem
!
! - --- Compute
!
        call calcul(stop_calc, option2, ligrel_model, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, &
                    'OUI')
!
! - --- Multiply values by 1/dt
!
        call multResuElem(resu_elem, 1.d0/deltat)

!
! - --- Add RESU_ELEM in MATR_ELEM
!
        call reajre(matr_elem, resu_elem, base)
!
    end if
!
end subroutine
