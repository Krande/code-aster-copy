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
subroutine vedith(model, loadNameJv, loadInfoJv, time, vect_elemz)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/inical.h"
#include "asterfort/load_list_info.h"
#include "asterfort/detrsd.h"
#include "asterfort/vemare.h"
#include "asterfort/reajre.h"
#include "asterfort/megeom.h"
#include "asterfort/gcnco2.h"
#include "asterfort/corich.h"
#include "asterfort/calcul.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: time, loadNameJv, loadInfoJv
    character(len=24), intent(inout) :: vect_elemz
!
! --------------------------------------------------------------------------------------------------
!
! Compute Dirichlet loads
!
! For Lagrange elements (AFFE_CHAR_THER) - T(given)
!
! --------------------------------------------------------------------------------------------------
!
! In  model          : name of model
! In  loadNameJv     : name of object for list of loads name
! In  loadInfoJv     : name of object for list of loads info
! In  time           : time (<CARTE>)
! IO  vect_elem      : name of vect_elem result
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbin = 3
    integer(kind=8), parameter :: nbout = 1
    character(len=8) :: lpain(nbin), lpaout(nbout)
    character(len=24) :: lchin(nbin), lchout(nbout)
!
    character(len=8) :: load_name, newnom
    character(len=16) :: option
    character(len=19) :: vect_elem
    character(len=24) :: ligrch, chgeom, resu_elem
    integer(kind=8) :: load_nume, nb_load, i_load
    character(len=24), pointer :: v_load_name(:) => null()
    integer(kind=8), pointer :: v_load_info(:) => null()
    aster_logical :: load_empty
    character(len=1) :: base
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    newnom = '.0000000'
    resu_elem = '&&VEDITH.???????'
    base = 'V'
!
! - Init fields
!
    call inical(nbin, lpain, lchin, nbout, lpaout, lchout)
!
! - Result name for vect_elem
!
    vect_elem = vect_elemz(1:19)
    if (vect_elem .eq. ' ') then
        vect_elem = '&&VEDITH'
    end if
!
! - Loads
!
    call load_list_info(load_empty, nb_load, v_load_name, v_load_info, &
                        loadNameJv, loadInfoJv)
!
! - Allocate result
!
    call detrsd('VECT_ELEM', vect_elem)
    call vemare(base, vect_elem, model)
    call reajre(vect_elem, ' ', base)
    if (load_empty) then
        goto 99
    end if
!
! - Geometry field
!
    call megeom(model, chgeom)
!
! - Input fields
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PINSTR'
    lchin(2) = time(1:19)
!
! - Output field
!
    lpaout(1) = 'PVECTTR'
!
! - Computation
!
    do i_load = 1, nb_load
        load_name = v_load_name(i_load) (1:8)
        load_nume = v_load_info(i_load+1)
        if ((load_nume .gt. 0) .and. (load_nume .le. 3)) then
            ligrch = load_name//'.CHTH.LIGRE'
!
! --------- Input field
!
            lchin(3) = load_name//'.CHTH.CIMPO.DESC'
            if (load_nume .eq. 1) then
                option = 'THER_DDLI_R'
                lpain(3) = 'PDDLIMR'
            else if (load_nume .eq. 2) then
                option = 'THER_DDLI_F'
                lpain(3) = 'PDDLIMF'
            else if (load_nume .eq. 3) then
                option = 'THER_DDLI_F'
                lpain(3) = 'PDDLIMF'
            end if
!
! --------- Generate new RESU_ELEM name
!
            call gcnco2(newnom)
            resu_elem(10:16) = newnom(2:8)
            call corich('E', resu_elem, ichin_=i_load)
            lchout(1) = resu_elem
!
! --------- Computation
!
            call calcul('S', option, ligrch, nbin, lchin, &
                        lpain, nbout, lchout, lpaout, base, &
                        'OUI')
!
! --------- Copying output field
!
            call reajre(vect_elem, lchout(1), base)
!
        end if
    end do
!
99  continue
!
    vect_elemz = vect_elem//'.RELR'
!
    call jedema()
end subroutine
