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
subroutine vecgme(model, caraElem, matez, matecoz, loadNameJvZ, loadInfoJvZ, &
                  inst_curr, disp_prevz, disp_cumu_instz, vect_elemz, inst_prev, &
                  compor, ligrel_calcz, vite_currz, acce_currz, strx_prevz)
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
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: caraElem
    character(len=*), intent(in) :: matez, matecoz
    real(kind=8), intent(in) :: inst_curr
    character(len=*), intent(in) :: disp_prevz
    character(len=*), intent(in) :: disp_cumu_instz
    character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
    character(len=*), intent(inout) :: vect_elemz
    real(kind=8), intent(in) :: inst_prev
    character(len=24), intent(in) :: compor
    character(len=*), intent(in) :: ligrel_calcz
    character(len=*), intent(in) :: vite_currz
    character(len=*), intent(in) :: acce_currz
    character(len=*), intent(in) :: strx_prevz
!
! --------------------------------------------------------------------------------------------------
!
! Compute Neumann loads
!
! Undead loads - Depending on geometry or speed - Vector
!
! --------------------------------------------------------------------------------------------------
!
! In  model          : name of model
! In  mate           : name of material characteristics (field)
! In  caraElem       : name of elementary characteristics (field)
! In  loadNameJv     : name of object for list of loads name
! In  loadInfoJv     : name of object for list of loads info
! In  inst_prev      : previous time
! In  inst_curr      : current time
! In  ligrel_calc    : LIGREL to compute
! In  vite_curr      : speed at current time
! In  acce_curr      : acceleration at current time
! In  disp_prev      : displacement at beginning of current time
! In  strx_prev      : fibers information at beginning of current time
! In  disp_cumu_inst : displacement increment from beginning of current time
! In  compor         : name of comportment definition (field)
! IO  vectElem       : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbout = 1
    character(len=8) :: lpain(LOAD_NEUM_NBMAXIN), lpaout(nbout)
    character(len=19) :: lchin(LOAD_NEUM_NBMAXIN), lchout(nbout)
    character(len=1), parameter :: stop = 'S', jvBaseTemporary = "V"
    character(len=4), parameter :: loadApply = "Suiv"
    character(len=8) :: newnom, loadName
    integer :: nbLoad, iLoad
    integer :: loadNume
    integer :: nb_in_prep
    real(kind=8) :: inst_theta
    character(len=24) :: ligrel_calc, mate, mateco
    character(len=19) :: vectElem, resuElem
    character(len=19) :: disp_prev, disp_cumu_inst, vite_curr, acce_curr, strx_prev
    character(len=24) :: loadNameJv
    character(len=24), pointer :: listLoadName(:) => null()
    character(len=24) :: loadInfoJv
    integer, pointer :: listLoadInfo(:) => null()
    aster_logical :: load_empty
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    newnom = '.0000000'
    resuElem = '&&VECGME.0000000'
    mate = matez
    mateco = matecoz
    loadNameJv = loadNameJvZ
    loadInfoJv = loadInfoJvZ
    disp_prev = disp_prevz
    strx_prev = strx_prevz
    disp_cumu_inst = disp_cumu_instz
    vite_curr = vite_currz
    acce_curr = acce_currz
    ligrel_calc = ligrel_calcz
    inst_theta = 0.d0
    if (ligrel_calc .eq. ' ') then
        ligrel_calc = model(1:8)//'.MODELE'
    end if
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Result name for vectElem
    vectElem = vect_elemz
    if (vectElem .eq. ' ') then
        vectElem = '&&VECGME'
    end if

! - Loads
    call load_list_info(load_empty, nbLoad, listLoadName, listLoadInfo, &
                        loadNameJv, loadInfoJv)

! - Allocate result
    call detrsd('VECT_ELEM', vectElem)
    call memare(jvBaseTemporary, vectElem, model, 'CHAR_MECA')
    call reajre(vectElem, ' ', jvBaseTemporary)
    if (load_empty) then
        goto 99
    end if

! - Preparing input fields
    call load_neum_prep(model, caraElem, mate, mateco, loadApply, inst_prev, &
                        inst_curr, inst_theta, LOAD_NEUM_NBMAXIN, nb_in_prep, lchin, &
                        lpain, disp_prev=disp_prev, disp_cumu_inst=disp_cumu_inst, &
                        compor=compor, strx_prev_=strx_prev, vite_curr_=vite_curr, &
                        acce_curr_=acce_curr)

! - Computation
    do iLoad = 1, nbLoad
        loadName = listLoadName(iLoad) (1:8)
        loadNume = listLoadInfo(nbLoad+iLoad+1)

        if (loadNume .eq. 4) then
! --------- Standard undead Neumann loads
            call load_neum_comp(stop, iLoad, loadName, loadNume, loadApply, &
                                ligrel_calc, LOAD_NEUM_NBMAXIN, nb_in_prep, lpain, lchin, &
                                jvBaseTemporary, resuElem, vectElem)

! --------- Composite undead Neumann loads (EVOL_CHAR)
            call compEvolChar(model, caraElem, inst_curr, jvBaseTemporary, &
                              iLoad, loadName, loadApply, ligrel_calc, &
                              nb_in_prep, lpain, lchin, &
                              resuElem, vectElem, &
                              disp_prev, disp_cumu_inst, strx_prev, vite_curr)
        end if
    end do
!
99  continue
!
    vect_elemz = vectElem//'.RELR'
!
    call jedema()
end subroutine
