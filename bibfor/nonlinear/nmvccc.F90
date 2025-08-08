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

subroutine nmvccc(model, nbin, nbout, lpain, lchin, &
                  lpaout, lchout, exis_temp, exis_hydr, exis_ptot, &
                  exis_sech, exis_epsa, calc_meta, base, vect_elem)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/reajre.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: model
    integer(kind=8), intent(in) :: nbout
    integer(kind=8), intent(in) :: nbin
    character(len=8), intent(in) :: lpain(nbin)
    character(len=19), intent(in) :: lchin(nbin)
    character(len=8), intent(in) :: lpaout(nbout)
    character(len=19), intent(inout) :: lchout(nbout)
    aster_logical, intent(in) :: exis_temp
    aster_logical, intent(in) :: exis_hydr
    aster_logical, intent(in) :: exis_ptot
    aster_logical, intent(in) :: exis_sech
    aster_logical, intent(in) :: exis_epsa
    aster_logical, intent(in) :: calc_meta
    character(len=1), intent(in) :: base
    character(len=19), intent(in) :: vect_elem
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Command variables - Compute elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  model          : name of model
! In  nbin           : number of input fields
! In  nbout          : number of output fields
! In  lpain          : list of input field parameters
! In  lchin          : list of input fields
! In  lpaout         : list of output field parameters
! IO  lchout         : list of output fields
! In  exis_temp      : .true. if temperature variable command exists
! In  exis_hydr      : .true. if hydration variable command exists
! In  exis_ptot      : .true. if total pressure (THM) variable command exists
! In  exis_sech      : .true. if drying variable command exists
! In  exis_epsa      : .true. if non-elastic strain variable command exists
! In  calc_meta      : .true. to compute metallurgy variable command
! In  base           : JEVEUX base to create objects
! In  vect_elem      : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6) :: masque
    character(len=16) :: option
    character(len=24) :: ligrmo
    integer(kind=8) :: nbr
!
! --------------------------------------------------------------------------------------------------
!
    nbr = 0
    masque = '.VEXXX'
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)
!
! - Temperature
!
    if (exis_temp) then
        nbr = nbr+1
        call codent(nbr, 'D0', masque(4:6))
        lchout(1) = vect_elem(1:8)//masque
        option = 'CHAR_MECA_TEMP_R'
        call calcul('C', option, ligrmo, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, &
                    'OUI')
        call reajre(vect_elem, lchout(1), base)
    end if
!
! - Hydration
!
    if (exis_hydr) then
        nbr = nbr+1
        call codent(nbr, 'D0', masque(4:6))
        lchout(1) = vect_elem(1:8)//masque
        option = 'CHAR_MECA_HYDR_R'
        call calcul('S', option, ligrmo, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, 'V', &
                    'OUI')
        call reajre(vect_elem, lchout(1), 'V')
    end if
!
! - Total pressure (THM)
!
    if (exis_ptot) then
        nbr = nbr+1
        call codent(nbr, 'D0', masque(4:6))
        lchout(1) = vect_elem(1:8)//masque
        option = 'CHAR_MECA_PTOT_R'
        call calcul('S', option, ligrmo, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, 'V', &
                    'OUI')
        call reajre(vect_elem, lchout(1), 'V')
    end if
!
! - Drying
!
    if (exis_sech) then
        nbr = nbr+1
        call codent(nbr, 'D0', masque(4:6))
        lchout(1) = vect_elem(1:8)//masque
        option = 'CHAR_MECA_SECH_R'
        call calcul('S', option, ligrmo, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, 'V', &
                    'OUI')
        call reajre(vect_elem, lchout(1), 'V')
    end if
!
! - Non-elastic strain
!
    if (exis_epsa) then
        nbr = nbr+1
        call codent(nbr, 'D0', masque(4:6))
        lchout(1) = vect_elem(1:8)//masque
        option = 'CHAR_MECA_EPSA_R'
        call calcul('S', option, ligrmo, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, 'V', &
                    'OUI')
        call reajre(vect_elem, lchout(1), 'V')
    end if
!
! - Metallurgy
!
    if (calc_meta) then
        nbr = 6
        call codent(nbr, 'D0', masque(4:6))
        lchout(1) = vect_elem(1:8)//masque
        option = 'CHAR_MECA_META_Z'
        call calcul('S', option, ligrmo, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, 'V', &
                    'OUI')
        call reajre(vect_elem, lchout(1), 'V')
    end if
!
end subroutine
