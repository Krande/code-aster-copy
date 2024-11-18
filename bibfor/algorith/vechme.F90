! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine vechme(stop, modelz, loadNameJvZ, loadInfoJvz, inst, &
                  caraElem, mate, mateco, vect_elemz, varc_currz, ligrel_calcz, &
                  nharm, basez)
!
    use loadCompute_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/load_list_info.h"
#include "asterfort/load_neum_prep.h"
#include "asterfort/load_neum_comp.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "LoadTypes_type.h"
!
    character(len=1), intent(in) :: stop
    character(len=*), intent(in) :: modelz
    character(len=*), intent(in) :: loadNameJvZ, loadInfoJvz
    real(kind=8), intent(in) :: inst(3)
    character(len=*), intent(in) :: caraElem
    character(len=*), intent(in) :: mate, mateco
    character(len=*), intent(inout) :: vect_elemz
    character(len=*), optional, intent(in) :: varc_currz
    character(len=*), optional, intent(in) :: ligrel_calcz
    character(len=1), optional, intent(in) :: basez
    integer, optional, intent(in) :: nharm
!
! --------------------------------------------------------------------------------------------------
!
! Compute Neumann loads
!
! Dead and fixed loads
!
! --------------------------------------------------------------------------------------------------
!
! In  stop           : continue or stop computation if no loads on elements
! In  model          : name of model
! In  mate           : name of material characteristics (field)
! In  mateco         : mane of coded material
! In  caraElem       : name of elementary characteristics (field)
! In  loadNameJv     : name of object for list of loads name
! In  loadInfoJv     : name of object for list of loads info
! In  inst           : times informations
! In  ligrel_calc    : LIGREL to compute
! In  varc_curr      : command variable for current time
! IO  vectElem      : name of vectElem result
! In  nharm          : Fourier mode
!
! ATTENTION :
!   LE vectElem (VECELZ) RESULTAT A 1 PARTICULARITE :
!   CERTAINS RESUELEM NE SONT PAS DES RESUELEM MAIS DES CHAM_NO (.VEASS)
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbout = 1
    character(len=8) :: lpain(LOAD_NEUM_NBMAXIN), lpaout(nbout)
    character(len=19) :: lchin(LOAD_NEUM_NBMAXIN), lchout(nbout)
    character(len=4), parameter :: loadApply = "Dead"
    integer :: nbLoad, iLoad
    integer :: loadNume
    integer :: nb_in_prep
    real(kind=8) :: inst_prev, inst_curr, inst_theta
    character(len=8) :: loadName
    character(len=24) :: ligrel_calc, model
    character(len=19) :: vectElem, varc_curr, resuElem
    character(len=24) :: loadNameJv
    character(len=24), pointer :: listLoadName(:) => null()
    character(len=24) :: loadInfoJv
    integer, pointer :: listLoadInfo(:) => null()
    aster_logical :: load_empty
    character(len=1) :: base
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    resuElem = '&&VECHME.0000000'
    model = modelz
    loadNameJv = loadNameJvZ
    loadInfoJv = loadInfoJvz
    varc_curr = ' '
    if (present(varc_currz)) then
        varc_curr = varc_currz
    end if
    ligrel_calc = model(1:8)//'.MODELE'
    if (present(ligrel_calcz)) then
        ligrel_calc = ligrel_calcz
    end if
    inst_prev = inst(1)
    inst_curr = inst(1)+inst(2)
    inst_theta = inst(3)
    base = 'V'
    if (present(basez)) then
        base = basez
    end if
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Result name for vectElem
    vectElem = vect_elemz
    if (vectElem .eq. ' ') then
        vectElem = '&&VECHME'
    end if

! - Loads
    call load_list_info(load_empty, nbLoad, listLoadName, listLoadInfo, &
                        loadNameJv, loadInfoJv)

! - Allocate result
    call detrsd('VECT_ELEM', vectElem)
    call memare(base, vectElem, model, 'CHAR_MECA')
    call reajre(vectElem, ' ', base)
    if (load_empty) then
        goto 99
    end if

! - Preparing input fields
    if (present(nharm)) then
        call load_neum_prep(model, caraElem, mate, mateco, loadApply, inst_prev, &
                            inst_curr, inst_theta, LOAD_NEUM_NBMAXIN, nb_in_prep, lchin, &
                            lpain, varc_curr=varc_curr, nharm=nharm)
    else
        call load_neum_prep(model, caraElem, mate, mateco, 'Dead', inst_prev, &
                            inst_curr, inst_theta, LOAD_NEUM_NBMAXIN, nb_in_prep, lchin, &
                            lpain, varc_curr=varc_curr)
    end if

! - Computation
    do iLoad = 1, nbLoad
        loadName = listLoadName(iLoad) (1:8)
        loadNume = listLoadInfo(nbLoad+iLoad+1)

        if ((loadNume .gt. 0 .and. loadNume .lt. 4) .or. (loadNume .eq. 55)) then
! --------- Standard dead Neumann loads
            call load_neum_comp(stop, iLoad, loadName, loadNume, loadApply, &
                                ligrel_calc, LOAD_NEUM_NBMAXIN, nb_in_prep, lpain, lchin, &
                                base, resuElem, vectElem)

! --------- Composite dead Neumann loads (EVOL_CHAR)
            call compEvolChar(model, caraElem, inst_prev, base, &
                              iLoad, loadName, loadApply, ligrel_calc, &
                              nb_in_prep, lpain, lchin, &
                              resuElem, vectElem)

        end if
    end do
!
99  continue
!
    vect_elemz = vectElem//'.RELR'
!
    call jedema()
end subroutine
