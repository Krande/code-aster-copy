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
! aslint: disable=W1306
!
subroutine ap_infast_n(mesh, newgeo, pair_tole, dist_ratio, nb_elem_mast, &
                       list_elem_mast, nb_elem_slav, list_elem_slav, elem_slav_flag, &
                       nb_mast_start, elem_mast_start, nb_slav_start, elem_slav_start, &
                       sdappa, list_node_mast, nb_node_mast)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jelira.h"
#include "asterfort/jexatr.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/apcoor.h"
#include "asterfort/aptype.h"
#include "asterfort/prjint_ray.h"
#include "asterfort/gtctma.h"
!#include "asterfort/gtclno.h"
#include "asterfort/gtclno_n.h"
#include "asterfort/gtlmex.h"
#include "asterfort/codent.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: newgeo
    real(kind=8), intent(in) :: pair_tole, dist_ratio
    integer(kind=8), intent(in) :: nb_elem_mast
    integer(kind=8), intent(in) :: list_elem_mast(nb_elem_mast)
    integer(kind=8), intent(in) :: nb_elem_slav
    integer(kind=8), intent(in) :: list_elem_slav(nb_elem_slav)
    integer(kind=8), pointer :: elem_slav_flag(:)
    integer(kind=8), intent(out) :: nb_mast_start
    integer(kind=8), intent(out) :: elem_mast_start(nb_elem_slav)
    integer(kind=8), intent(out) :: nb_slav_start
    integer(kind=8), intent(out) :: elem_slav_start(nb_elem_slav)
    character(len=19), intent(in) :: sdappa
    integer(kind=8), intent(in) :: list_node_mast(*)
    integer(kind=8), intent(in) ::  nb_node_mast
!integer, intent(in) :: i_zone
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Find initial elements for pairing by "fast' PANG method
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  newgeo           : name of field for geometry update from initial coordinates of nodes
! In  pair_tole        : tolerance for pairing
! In  l_not_memory     : .true. is algorithm is the "non-memory" version for standard PANG
! In  nb_elem_mast     : number of master elements on current zone
! In  list_elem_mast   : name of datastructure for list of master elements on current zone
! In  nb_elem_slav     : number of slave elements on current zone
! In  list_elem_slav   : name of datastructure for list of slave elements on current zone
! IO  elem_mast_flag   : flag to mark master elements already tracked
! IO  elem_slav_flag   : flag to mark slave elements already tracked
! Out nb_mast_start    : number of master elements to start algorithm
! Out elem_mast_start  : list of master elements to start algorithm
! Out nb_slav_start    : number of slave elements to start algorithm
! Out elem_slav_start  : list of slave elements to start algorithm
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: elem_type_nume
    integer(kind=8) :: elem_slav_nbnode, elem_slav_dime, elem_slav_nume, elem_slav_indx
    character(len=24) :: conx_inve
    character(len=8) :: elem_slav_type, elem_slav_code, knuzo
    real(kind=8) :: elem_slav_coor(27)
    integer(kind=8) ::  elin_slav_nbnode
    integer(kind=8) :: elem_mast_nbnode, elem_mast_dime, elem_mast_nume, elem_mast_indx
    character(len=8) :: elem_mast_type, elem_mast_code, elem_slav_name, elem_mast_name
    real(kind=8) :: elem_mast_coor(27)
    integer(kind=8) :: elin_mast_nbnode
    character(len=8) :: elin_mast_code, elin_slav_code
    integer(kind=8) :: slav_indx_mini, mast_indx_mini
    integer(kind=8) :: jv_geom, i_zone
    integer(kind=8) :: i_elem_slav, i_elem_mast
    integer(kind=8) :: nb_poin_inte, nume_node_cl, nb_el_ma_ax
    real(kind=8) :: poin_inte_es(32), poin_inte_ma(32), inte_weight, center(3)
    integer(kind=8), pointer :: v_mesh_typmail(:) => null()
    integer(kind=8), pointer :: v_mesh_connex(:) => null()
    integer(kind=8), pointer :: v_connex_lcum(:) => null()
!    integer, pointer :: list_node_mast(:) => null()
    integer(kind=8), pointer :: v_cninv(:) => null()
    integer(kind=8), pointer :: v_cninv_lcum(:) => null()
    integer(kind=8) :: list_el_ma_ax(nb_elem_mast)
    aster_logical :: debug
!
! --------------------------------------------------------------------------------------------------
!
    debug = ASTER_FALSE
    mast_indx_mini = minval(list_elem_mast)
    slav_indx_mini = minval(list_elem_slav)
    nb_mast_start = 0
    nb_slav_start = 0
    elem_slav_start(1:nb_elem_slav) = 0
    elem_mast_start(1:nb_elem_slav) = 0
    i_zone = 0
    ASSERT(i_zone .le. 100)
    call codent(i_zone-1, 'G', knuzo)
    conx_inve = sdappa(1:19)
    call jeveuo(conx_inve, 'L', vi=v_cninv)
    call jeveuo(jexatr(conx_inve, 'LONCUM'), 'L', vi=v_cninv_lcum)
!    call jeveuo(sdappa(1:19)//'.LM'//knuzo(1:2),'L', vi = list_node_mast)
!    call jelira(sdappa(1:19)//'.LM'//knuzo(1:2),'LONMAX',nb_node_mast)
    if (nb_node_mast .eq. 0) then
        go to 100
    end if
!
! - Access to mesh
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=v_mesh_typmail)
    call jeveuo(mesh//'.CONNEX', 'L', vi=v_mesh_connex)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=v_connex_lcum)
!
! - Access to updated geometry
!
    call jeveuo(newgeo(1:19)//'.VALE', 'L', jv_geom)
    if (debug) then
        write (*, *) "Find start elements"
    end if
!
! - Loop on slave elements
!
    do i_elem_slav = 1, nb_elem_slav
!
! ----- Current slave element
!
        elem_slav_nume = list_elem_slav(i_elem_slav)
        elem_slav_indx = elem_slav_nume+1-slav_indx_mini
        elem_type_nume = v_mesh_typmail(elem_slav_nume)
        if (debug) then
            elem_slav_name = int_to_char8(elem_slav_nume)
            write (*, *) "Slave element", i_elem_slav, elem_slav_nume, elem_slav_name
        end if
!
! ----- Already tracked ?
!
        if (elem_slav_flag(elem_slav_indx) .eq. 0) then
            if (debug) then
                write (*, *) "Slave element not yet tracked"
            end if
!
! --------- Get informations about slave element
!
            call jenuno(jexnum('&CATA.TM.NOMTM', elem_type_nume), elem_slav_type)
            call aptype(elem_slav_type, &
                        elem_slav_nbnode, elem_slav_code, elem_slav_dime)
!
! --------- Get coordinates of slave element
!
            call apcoor(v_mesh_connex, v_connex_lcum, jv_geom, &
                        elem_slav_nume, elem_slav_nbnode, elem_slav_dime, &
                        elem_slav_coor)
!
! --------- Cut slave element in linearized sub-elements (SEG2 or TRIA3)
!
            if (elem_slav_code .eq. "TR6") then
                elin_slav_code = "TR3"
                elin_slav_nbnode = 3
            elseif (elem_slav_code .eq. "QU8" .or. elem_slav_code .eq. "QU9") then
                elin_slav_code = "QU4"
                elin_slav_nbnode = 4
            elseif (elem_slav_code .eq. "SE3") then
                elin_slav_code = "SE2"
                elin_slav_nbnode = 2
            else
                elin_slav_code = elem_slav_code
                elin_slav_nbnode = elem_slav_nbnode
            end if
!
! --------- Find the closest master node from center
!
            call gtctma(elem_slav_coor, elem_slav_nbnode, elem_slav_code, elem_slav_dime, center)
            call gtclno_n(jv_geom, list_node_mast, nb_node_mast, center, nume_node_cl)
!
! --------- Loop on master elements next to the closest master node
!
            call gtlmex(v_cninv, v_cninv_lcum, nume_node_cl, nb_elem_mast, list_elem_mast, &
                        list_el_ma_ax, nb_el_ma_ax)
            do i_elem_mast = 1, nb_el_ma_ax
!
! ------------- Current master element
!
                elem_mast_nume = list_el_ma_ax(i_elem_mast)
                elem_mast_indx = elem_mast_nume+1-mast_indx_mini
                elem_type_nume = v_mesh_typmail(elem_mast_nume)
                if (debug) then
                    elem_mast_name = int_to_char8(elem_mast_nume)
                    write (*, *) "Master element", i_elem_mast, elem_mast_nume, elem_mast_name
                end if
!
! ------------- Already tracked ?
!
                if (.true.) then
                    if (debug) then
                        write (*, *) "Master element not yet tracked"
                    end if
!
! ----------------- Get informations about master element
!
                    call jenuno(jexnum('&CATA.TM.NOMTM', elem_type_nume), elem_mast_type)
                    call aptype(elem_mast_type, &
                                elem_mast_nbnode, elem_mast_code, elem_mast_dime)
!
! ----------------- Get coordinates of master element
!
                    call apcoor(v_mesh_connex, v_connex_lcum, jv_geom, &
                                elem_mast_nume, elem_mast_nbnode, elem_mast_dime, &
                                elem_mast_coor)
!
! ----------------- Cut master element in linearized sub-elements (SEG2 or TRIA3)
!
                    if (elem_mast_code .eq. "TR6") then
                        elin_mast_code = "TR3"
                        elin_mast_nbnode = 3
                    elseif (elem_mast_code .eq. "QU8" .or. elem_mast_code .eq. "QU9") then
                        elin_mast_code = "QU4"
                        elin_mast_nbnode = 4
                    elseif (elem_mast_code .eq. "SE3") then
                        elin_mast_code = "SE2"
                        elin_mast_nbnode = 2
                    else
                        elin_mast_code = elem_mast_code
                        elin_mast_nbnode = elem_mast_nbnode
                    end if
!
! ----------------- Projection/intersection of elements in slave parametric space
!
                    call prjint_ray(pair_tole, dist_ratio, elem_mast_dime, &
                                    elin_slav_nbnode, elem_slav_coor, elin_slav_code, &
                                    elin_mast_nbnode, elem_mast_coor, elin_mast_code, &
                                    poin_inte_ma, poin_inte_es, inte_weight, nb_poin_inte)
!
! ----------------- Set start elements
!
                    if (inte_weight .gt. 100*pair_tole) then
                        elem_mast_start(1) = elem_mast_nume
                        nb_mast_start = 1
                        elem_slav_start(1) = elem_slav_nume
                        nb_slav_start = 1
                        elem_slav_flag(elem_slav_indx) = 1
                        if (debug) then
                            elem_mast_name = int_to_char8(elem_mast_nume)
                            elem_slav_name = int_to_char8(elem_slav_nume)
                            write (*, *) "Depart trouvé(M/S): ", elem_mast_name, elem_slav_name
                        end if
                        goto 100
                    end if
                else
                    if (debug) then
                        write (*, *) "Master element not yet tracked"
                    end if
                end if
            end do
        else
            if (debug) then
                write (*, *) "Slave element already tracked"
            end if
        end if
        elem_slav_flag(elem_slav_indx) = 2
    end do
100 continue
!
end subroutine
