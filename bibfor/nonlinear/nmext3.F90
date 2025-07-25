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

subroutine nmext3(mesh, field, field_type, field_s, &
                  nb_cmp, nb_elem, nb_poin, nb_spoi, &
                  type_extr_elem, type_extr, type_extr_cmp, type_sele_cmp, &
                  list_elem, list_poin, list_spoi, list_cmp, &
                  work_poin, work_elem)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/celces.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmextj.h"
#include "asterfort/sdmpic.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: field
    character(len=24), intent(in) :: field_type
    character(len=24), intent(in) :: field_s
    integer(kind=8), intent(in) :: nb_cmp
    integer(kind=8), intent(in) :: nb_elem
    integer(kind=8), intent(in) :: nb_poin
    integer(kind=8), intent(in) :: nb_spoi
    character(len=8), intent(in) :: type_extr_elem
    character(len=8), intent(in) :: type_extr
    character(len=8), intent(in) :: type_extr_cmp
    character(len=8), intent(in) :: type_sele_cmp
    character(len=24), intent(in) :: list_elem
    character(len=24), intent(in) :: list_poin
    character(len=24), intent(in) :: list_spoi
    character(len=24), intent(in) :: list_cmp
    character(len=19), intent(in) :: work_poin
    character(len=19), intent(in) :: work_elem
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Field extraction datastructure
!
! Extract value(s) at point and store them in working vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  field            : name of field
! In  field_type       : type of field (name in results datastructure)
! In  field_s          : name of reduced field (CHAM_ELEM_S)
! In  nb_cmp           : number of components
! In  nb_elem          : number of elements
! In  nb_poin          : number of points (Gauss)
! In  nb_spoi          : number of subpoints
! In  type_extr_elem   : type of extraction by element
! In  type_extr        : type of extraction
! In  type_extr_cmp    : type of extraction for components
! In  type_sele_cmp    : type of selection for components NOM_CMP or NOM_VARI
! In  list_elem        : name of object contains list of elements
! In  list_poin        : name of object contains list of points (Gauss)
! In  list_spoi        : name of object contains list of subpoints
! In  list_cmp         : name of object contains list of components
! In  work_poin        : working vector to save point (Gauss) values
! In  work_elem        : working vector to save element values
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_para_maxi = 20
    real(kind=8) :: vale_resu(nb_para_maxi)
    integer(kind=8) :: i_elem, i_poin, i_spoi, i_elem_r, i_poin_r, i_spoi_r, elem_nume
    integer(kind=8) :: poin_nume, spoi_nume, iret
    integer(kind=8) :: nb_elem_poin, nb_elem_spoi, nb_poin_r, nb_spoi_r, nb_elem_effe
    integer(kind=8) :: i_vale, nb_vale
    real(kind=8) :: valr, val2r
    integer(kind=8) :: jcesd, jcesl, jcesv, jcesc
    integer(kind=8), pointer :: v_list_elem(:) => null()
    integer(kind=8), pointer :: v_list_poin(:) => null()
    integer(kind=8), pointer :: v_list_spoi(:) => null()
    real(kind=8), pointer :: v_work_poin(:) => null()
    real(kind=8), pointer :: v_work_elem(:) => null()
    real(kind=8), pointer :: v_init_poin(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_cmp .le. nb_para_maxi)
    ASSERT(nb_poin .ge. 1)
    ASSERT(nb_spoi .ge. 1)
!
! - Conversion to reduced field (CHAM_ELEM_S)
!
    call jeexin(field_s, iret)
    if (iret .eq. 0) then
        call sdmpic('CHAM_ELEM', field)
        call celces(field, 'V', field_s)
    end if
!
! - Access to reduced field (CHAM_ELEM_S)
!
    call jeveuo(field_s(1:19)//'.CESD', 'L', jcesd)
    call jeveuo(field_s(1:19)//'.CESL', 'L', jcesl)
    call jeveuo(field_s(1:19)//'.CESV', 'L', jcesv)
    call jeveuo(field_s(1:19)//'.CESC', 'L', jcesc)
!
! - Get access to working vectors
!
    call jeveuo(work_elem, 'E', vr=v_work_elem)
    call jeveuo(work_poin, 'E', vr=v_work_poin)
!
! - Copy initial point values
!
    AS_ALLOCATE(vr=v_init_poin, size=nb_poin*nb_spoi*nb_cmp)
    v_init_poin(:) = v_work_poin(:)
!
! - Get access to lists
!
    call jeveuo(list_elem, 'L', vi=v_list_elem)
    call jeveuo(list_poin, 'L', vi=v_list_poin)
    call jeveuo(list_spoi, 'L', vi=v_list_spoi)
!
! - Effective number of elements (for which there is at least one point and one subpoint)
!
    nb_elem_effe = 0
    do i_elem = 1, nb_elem
        elem_nume = v_list_elem(i_elem)
        nb_elem_poin = zi(jcesd+5+4*(elem_nume-1))
        nb_elem_spoi = zi(jcesd+5+4*(elem_nume-1)+1)
        if (nb_elem_poin*nb_elem_spoi .ge. 1) nb_elem_effe = nb_elem_effe+1
    end do
    ASSERT(nb_elem_effe .ge. 1)
!
! - Loop on elements
!
    do i_elem = 1, nb_elem
!
! ----- Copy initial value
!
        v_work_poin(:) = v_init_poin(:)
!
! ----- Current element
!
        elem_nume = v_list_elem(i_elem)
!
! ----- Number of points/subpoints on current element
!
        nb_elem_poin = zi(jcesd+5+4*(elem_nume-1))
        nb_elem_spoi = zi(jcesd+5+4*(elem_nume-1)+1)
!
! ----- Check
!
        nb_poin_r = min(nb_poin, nb_elem_poin)
        nb_spoi_r = min(nb_spoi, nb_elem_spoi)
        if (nb_poin_r*nb_spoi_r .eq. 0) goto 100
!
! --------- Extract and set point/subpoint value(s) by element
!
        do i_poin = 1, nb_poin_r
            do i_spoi = 1, nb_spoi_r
!
! ----------------- Index on point/subpoint
!
                poin_nume = v_list_poin(i_poin)
                spoi_nume = v_list_spoi(i_spoi)
!
! ----------------- Extract value at Gauss point
!
                call nmextj(field_type, nb_cmp, list_cmp, type_extr_cmp, type_sele_cmp, &
                            poin_nume, spoi_nume, nb_vale, i_elem, elem_nume, &
                            jcesd, jcesv, jcesl, jcesc, vale_resu)
!
! ----------------- Select index in working vectors (point/subpoint)
!
                if (type_extr_elem .eq. 'VALE') then
                    i_poin_r = i_poin
                    i_spoi_r = i_spoi
                else
                    i_poin_r = 1
                    i_spoi_r = 1
                end if
!
! ----------------- Save values in working vector (element)
!
                do i_vale = 1, nb_vale
                    valr = vale_resu(i_vale)
                    val2r = v_work_poin(1+nb_poin*nb_spoi*(i_vale-1)+ &
                                        nb_spoi*(i_poin_r-1)+ &
                                        (i_spoi_r-1))

                    if (type_extr_elem .eq. 'VALE') then
                        v_work_poin(1+nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin_r-1)+ &
                                    (i_spoi_r-1)) = valr
                    else if (type_extr_elem .eq. 'MAX') then
                        v_work_poin(1+nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin_r-1)+ &
                                    (i_spoi_r-1)) = max(valr, val2r)
                    else if (type_extr_elem .eq. 'MIN') then
                        v_work_poin(1+nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin_r-1)+ &
                                    (i_spoi_r-1)) = min(valr, val2r)
                    else if (type_extr_elem .eq. 'MOY') then
                        v_work_poin(1+nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin_r-1)+ &
                                    (i_spoi_r-1)) = val2r+valr/(nb_poin_r*nb_spoi_r)
                    else
                        ASSERT(.false.)
                    end if
                end do
            end do
        end do
!
! --------- Extract and set point/subpoint value(s) by point/subpoint
!
        do i_poin = 1, merge(nb_poin_r, 1, type_extr_elem .eq. 'VALE')
            do i_spoi = 1, merge(nb_spoi_r, 1, type_extr_elem .eq. 'VALE')
!
! ----------------- Save values in working vector (point)
!
                do i_vale = 1, nb_vale
!
! ----------------- Select index in working vector (element)
!
                    if (type_extr .eq. 'VALE') then
                        i_elem_r = i_elem
                    else
                        i_elem_r = 1
                    end if
                    valr = v_work_poin(1+nb_poin*nb_spoi*(i_vale-1)+ &
                                       nb_spoi*(i_poin-1)+ &
                                       (i_spoi-1))

                    val2r = v_work_elem(1+nb_cmp*nb_poin*nb_spoi*(i_elem_r-1)+ &
                                        nb_poin*nb_spoi*(i_vale-1)+ &
                                        nb_spoi*(i_poin-1)+ &
                                        (i_spoi-1))
                    if (type_extr .eq. 'VALE') then
                        v_work_elem(1+nb_cmp*nb_poin*nb_spoi*(i_elem_r-1)+ &
                                    nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin-1)+ &
                                    (i_spoi-1)) = valr
                    else if (type_extr .eq. 'MAX') then
                        v_work_elem(1+nb_cmp*nb_poin*nb_spoi*(i_elem_r-1)+ &
                                    nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin-1)+ &
                                    (i_spoi-1)) = max(valr, val2r)
                    else if (type_extr .eq. 'MIN') then
                        v_work_elem(1+nb_cmp*nb_poin*nb_spoi*(i_elem_r-1)+ &
                                    nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin-1)+ &
                                    (i_spoi-1)) = min(valr, val2r)
                    else if (type_extr .eq. 'MAXI_ABS') then
                        v_work_elem(1+nb_cmp*nb_poin*nb_spoi*(i_elem_r-1)+ &
                                    nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin-1)+ &
                                    (i_spoi-1)) = max(abs(val2r), abs(valr))
                    else if (type_extr .eq. 'MINI_ABS') then
                        v_work_elem(1+nb_cmp*nb_poin*nb_spoi*(i_elem_r-1)+ &
                                    nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin-1)+ &
                                    (i_spoi-1)) = min(abs(val2r), abs(valr))
                    else if (type_extr .eq. 'MOY') then
                        v_work_elem(1+nb_cmp*nb_poin*nb_spoi*(i_elem_r-1)+ &
                                    nb_poin*nb_spoi*(i_vale-1)+ &
                                    nb_spoi*(i_poin-1)+ &
                                    (i_spoi-1)) = val2r+valr/nb_elem_effe
                    else
                        ASSERT(.false.)
                    end if
                end do
            end do
        end do
100     continue
    end do
    AS_DEALLOCATE(vr=v_init_poin)
!
end subroutine
