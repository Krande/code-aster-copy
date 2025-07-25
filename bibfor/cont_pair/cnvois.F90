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

subroutine cnvois(mesh, list_elem, conx_inve, nb_elem, elem_indx_mini, &
                  elem_indx_maxi, elem_neigh)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jecrec.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeecra.h"
#include "asterfort/jecroc.h"
#include "asterfort/gtvois.h"
#include "asterfort/utlisi.h"
#include "asterfort/assert.h"
!
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: conx_inve
    integer(kind=8), intent(in) :: nb_elem
    integer(kind=8), intent(in) :: list_elem(nb_elem)
    integer(kind=8), intent(in) :: elem_indx_mini
    integer(kind=8), intent(in) :: elem_indx_maxi
    character(len=24), intent(in) :: elem_neigh
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Create object of neigbours of a list of elements
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  conx_inve        : name of object for inverse connectivity
! In  nb_elem          : number of elements
! In  list_elem        : list of elements
! In  elem_indx_mini   : minimum index in list of elements
! In  elem_indx_maxi   : minimum index in list of elements
! In  elem_neigh       : name of object for neigbours of a list of elements
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_elem, i_neigh, aux(1), nb_find, lmail(1), jtypma
    integer(kind=8) :: elem_indx, elem_nume
    integer(kind=8) :: list_neigh(4), nb_neigh
    character(len=8) :: elem_code, elem_type
    integer(kind=8), pointer :: v_elem_neigh(:) => null()
    integer(kind=8), pointer :: v_connex(:) => null()
    integer(kind=8), pointer :: v_connex_lcum(:) => null()
    integer(kind=8), pointer :: v_conx_inv(:) => null()
    integer(kind=8), pointer :: v_inv_lcum(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jeveuo(mesh//'.TYPMAIL', 'L', jtypma)
!
! - Acces to mesh
!
    call jeveuo(mesh//'.CONNEX', 'L', vi=v_connex)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=v_connex_lcum)

!
! - Acces to conx_inve
!
    call jeveuo(conx_inve, 'L', vi=v_conx_inv)
    call jeveuo(jexatr(conx_inve, 'LONCUM'), 'L', vi=v_inv_lcum)

!
! - Create object (collection)
!
    call jecrec(elem_neigh, 'G V I', 'NU', 'CONTIG', 'VARIABLE', elem_indx_maxi+1-elem_indx_mini)
    call jeecra(elem_neigh, 'LONT', nb_elem*4+(elem_indx_maxi+1-elem_indx_mini-nb_elem)*4)
    do i_elem = 1, elem_indx_maxi+1-elem_indx_mini
        elem_nume = i_elem-1+elem_indx_mini
        nb_find = 0
        lmail(1) = elem_nume
        call utlisi('INTER', lmail, 1, list_elem, nb_elem, aux, 1, nb_find)
        if (nb_find .eq. 1) then
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtypma+elem_nume-1)), elem_type)
            select case (elem_type)
            case ('SEG2')
                nb_neigh = 2
            case ('SEG3')
                nb_neigh = 2
            case ('TRIA3')
                nb_neigh = 3
            case ('TRIA6')
                nb_neigh = 3
            case ('QUAD4')
                nb_neigh = 4
            case ('QUAD8')
                nb_neigh = 4
            case ('QUAD9')
                nb_neigh = 4
            case default
                ASSERT(.false.)
            end select
            call jecroc(jexnum(elem_neigh, i_elem))
            call jeecra(jexnum(elem_neigh, i_elem), 'LONMAX', ival=4)
        elseif (nb_find .eq. 0) then
            call jecroc(jexnum(elem_neigh, i_elem))
            call jeecra(jexnum(elem_neigh, i_elem), 'LONMAX', ival=4)
        else
            ASSERT(.false.)
        end if
    end do
!
! - Fill object
!
    do i_elem = 1, nb_elem
        elem_nume = list_elem(i_elem)
        elem_indx = elem_nume+1-elem_indx_mini
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(jtypma+elem_nume-1)), elem_type)
        select case (elem_type)
        case ('SEG2')
            nb_neigh = 2
            elem_code = 'SE2'
        case ('SEG3')
            nb_neigh = 2
            elem_code = 'SE3'
        case ('TRIA3')
            nb_neigh = 3
            elem_code = 'TR3'
        case ('TRIA6')
            nb_neigh = 3
            elem_code = 'TR6'
        case ('QUAD4')
            nb_neigh = 4
            elem_code = 'QU4'
        case ('QUAD8')
            nb_neigh = 4
            elem_code = 'QU8'
        case ('QUAD9')
            nb_neigh = 4
            elem_code = 'QU9'
        case default
            ASSERT(.false.)
        end select
        call jeveuo(jexnum(elem_neigh, elem_indx), 'E', vi=v_elem_neigh)
        call gtvois(v_connex, v_connex_lcum, list_elem, nb_elem, elem_nume, elem_code, &
                    v_conx_inv, v_inv_lcum, nb_neigh, list_neigh)
        do i_neigh = 1, 4
            v_elem_neigh(i_neigh) = list_neigh(i_neigh)
        end do
    end do
end subroutine
