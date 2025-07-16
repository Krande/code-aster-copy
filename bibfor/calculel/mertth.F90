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
subroutine mertth(model, loadNameJv, loadInfoJv, caraElem, mateco, &
                  time, time_move, temp_prev, temp_iter, matr_elem)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/load_list_info.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
!
    character(len=8), intent(in) :: model, caraElem
    character(len=24), intent(in) :: loadNameJv, loadInfoJv
    character(len=24), intent(in) :: mateco
    character(len=24), intent(in) :: time
    character(len=24), intent(in) :: time_move
    character(len=24), intent(in) :: temp_prev
    character(len=24), intent(in) :: temp_iter
    character(len=19), intent(inout) :: matr_elem
!
! --------------------------------------------------------------------------------------------------
!
! Thermic - Matrix
!
! Elementary matrix for transport (volumic and surfacic terms)
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  caraElem         : name of elementary characteristics (field)
! In  loadNameJv       : name of object for list of loads name
! In  loadInfoJv       : name of object for list of loads info
! In  time             : time (<CARTE>)
! In  time_move        : modified time (<CARTE>) for THER_NON_LINE_MO
! In  temp_prev        : previous temperature
! In  temp_iter        : temperature field at current Newton iteration
! IO  matr_elem        : name of matr_elem result
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbchmx
    parameter(nbchmx=3)
    integer(kind=8) :: nbopt(nbchmx), nligr(nbchmx)
    character(len=6) :: nomchp(nbchmx)
    character(len=7) :: nompar(nbchmx), nompaf(nbchmx)
!
    character(len=8) :: lpain(5), lpaout(1), load_name
    character(len=16) :: option, nomopr(nbchmx), nomopf(nbchmx)
    character(len=24) :: ligrel(2), lchin(5), lchout(1)
    character(len=24) :: chgeom, chcara(18)
    integer(kind=8) :: iret, nb_load, i_load, ilires, k, load_nume
    aster_logical :: load_empty
    character(len=24), pointer :: v_load_name(:) => null()
    integer(kind=8), pointer :: v_load_info(:) => null()
    data nomchp/'.FLUNL', '.HECHP', '.COEFH'/
    data nomopr/'                ', 'MTAN_THER_PARO_R', 'RIGI_THER_ECHA_R'/
    data nomopf/'MTAN_THER_FLUXNL', 'MTAN_THER_PARO_F', 'RIGI_THER_ECHA_F'/
    data nompar/'       ', 'PHECHPR', 'PCOEFHR'/
    data nompaf/'PFLUXNL', 'PHECHPF', 'PCOEFHF'/
    data nbopt/4, 5, 3/
    data nligr/1, 2, 1/
!
! --------------------------------------------------------------------------------------------------
!

!
! - Loads
!
    call load_list_info(load_empty, nb_load, v_load_name, v_load_info, &
                        loadNameJv, loadInfoJv)
!
    call megeom(model, chgeom)
    call mecara(caraElem, chcara)
!
    call jeexin(matr_elem(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        matr_elem = '&&METRIG'
        call memare('V', matr_elem, model(1:8), 'RIGI_THER')
    else
        call jedetr(matr_elem(1:19)//'.RELR')
    end if
!
    ligrel(1) = model(1:8)//'.MODELE'
!
    lpaout(1) = 'PMATTTR'
    lchout(1) = matr_elem(1:8)//'.ME001'
    ilires = 0
!
    if (model .ne. '        ') then
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PMATERC'
        lchin(2) = mateco
        lpain(3) = 'PCACOQU'
        lchin(3) = chcara(7)
        lpain(4) = 'PTEMPER'
        lchin(4) = temp_prev
        lpain(5) = 'PTEMPEI'
        lchin(5) = temp_iter
        option = 'RIGI_THER_TRANS'
        ilires = ilires+1
        call codent(ilires, 'D0', lchout(1) (12:14))
        call calcul('S', option, ligrel(1), 5, lchin, &
                    lpain, 1, lchout, lpaout, 'V', &
                    'OUI')
        call reajre(matr_elem, lchout(1), 'V')
    end if
!
    if (nb_load .gt. 0) then
        do i_load = 1, nb_load
            load_name = v_load_name(i_load) (1:8)
            load_nume = v_load_info(nb_load+i_load+1)
            if (load_nume .gt. 0) then
                ligrel(2) = load_name//'.CHTH.LIGRE'
                lpain(1) = 'PGEOMER'
                lchin(1) = chgeom
                lpain(3) = 'PINSTR'
                lchin(3) = time
                lpain(4) = 'PTEMPEI'
                lchin(4) = temp_iter
                lpain(5) = 'PDEPLAR'
                lchin(5) = '&&DEPPLU'
                lpaout(1) = 'PMATTTR'
                lchout(1) = matr_elem(1:8)//'.ME001'
                do k = 1, nbchmx
                    lchin(2) = load_name(1:8)//'.CHTH'//nomchp(k)//'.DESC'
                    call jeexin(lchin(2), iret)
                    if (iret .gt. 0) then
                        if (load_nume .eq. 1) then
                            option = nomopr(k)
                            lpain(2) = nompar(k)
                        else if (load_nume .eq. 2 .or. load_nume .eq. 3) then
                            option = nomopf(k)
                            lpain(2) = nompaf(k)
                        end if
                        if (option(11:14) .eq. 'PARO') then
                            lpain(3) = 'PINSTR'
                            lchin(3) = time_move
                        end if
                        if (k .eq. 2) lchin(4) = temp_iter
                        ilires = ilires+1
                        call codent(ilires, 'D0', lchout(1) (12:14))
                        call calcul('S', option, ligrel(nligr(k)), nbopt(k), lchin, &
                                    lpain, 1, lchout, lpaout, 'V', &
                                    'OUI')
                        call reajre(matr_elem, lchout(1), 'V')
                    end if
                end do
            end if
        end do
    end if
!
end subroutine
