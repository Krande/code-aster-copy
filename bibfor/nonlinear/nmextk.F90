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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmextk(mesh, model, &
                  keyw_fact, i_keyw_fact, &
                  field, field_type, field_s, field_disc, &
                  list_node, list_elem, list_poin, list_spoi, &
                  nb_node, nb_elem, nb_poin, nb_spoi, &
                  compor, list_cmp, list_vari, nb_cmp, type_sele_cmp)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxliis.h"
#include "asterfort/posddl.h"
#include "asterfort/utmess.h"
#include "asterfort/varinonu.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: mesh, model
    character(len=16), intent(in) :: keyw_fact
    integer(kind=8), intent(in) :: i_keyw_fact
    character(len=19), intent(in) :: field
    character(len=24), intent(in) :: field_type
    character(len=24), intent(in) :: field_s
    character(len=4), intent(in) :: field_disc
    integer(kind=8), intent(in) :: nb_node, nb_elem, nb_poin, nb_spoi
    character(len=24), intent(in) :: list_node, list_elem, list_poin, list_spoi
    character(len=19), optional, intent(in) :: compor
    character(len=24), intent(in) :: list_cmp
    character(len=24), intent(in) :: list_vari
    integer(kind=8), intent(out) :: nb_cmp
    character(len=8), intent(out) :: type_sele_cmp
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Extraction (OBSERVATION/SUIVI_DDL) utilities
!
! Get component(s)
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  keyw_fact        : factor keyword to read extraction parameters
! In  i_keyw_fact      : index of keyword to read extraction parameters
! In  field            : name of field
! In  field_type       : type of field (name in results datastructure)
! In  field_disc       : localization of field (discretization: NOEU, ELGA or ELEM)
! In  field_s          : name of reduced field (CHAM_ELEM_S)
! In  list_node        : name of object contains list of nodes
! In  nb_node          : number of nodes
! In  list_elem        : name of object contains list of elements
! In  nb_elem          : number of elements
! In  list_poin        : name of object contains list of points (Gauss)
! In  nb_poin          : number of points (Gauss)
! In  list_spoi        : name of object contains list of subpoints
! In  nb_spoi          : number of subpoints
! In  compor           : name of <CARTE> COMPOR
! In  list_cmp         : name of object contains list of components (NOM_CMP)
! In  list_vari        : name of object contains list of components (NOM_VARI)
! Out nb_cmp           : number of components
! Out type_sele_cmp    : type of selection for components NOM_CMP or NOM_VARI
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_para_maxi = 20
    integer(kind=8) :: n1
    integer(kind=8) :: iret, iad
    integer(kind=8) :: i_node, i_elem, i_cmp, ipi, ispi, ipar, i_cmp_maxi
    integer(kind=8) :: nb_cmp_maxi, nuno, nuddl
    integer(kind=8) :: nb_elem_poin, nb_elem_spoi, npi, nspi
    integer(kind=8) :: node_nume, elem_nume, num, snum
    character(len=8) :: node_name, elem_name, cmp_name
    character(len=8) :: cmp_vari_name
    integer(kind=8) :: i_vari
    character(len=16) :: valk(2)
    integer(kind=8) :: jcesd, jcesl, jcesv
    integer(kind=8) :: vali(4)
    character(len=8), pointer :: cesc(:) => null()
    character(len=8), pointer :: v_list_cmp(:) => null()
    character(len=16), pointer :: v_list_vari(:) => null()
    integer(kind=8), pointer :: v_list_node(:) => null()
    integer(kind=8), pointer :: v_list_elem(:) => null()
    integer(kind=8), pointer :: v_list_poin(:) => null()
    integer(kind=8), pointer :: v_list_spoi(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_cmp = 0
    type_sele_cmp = ' '
!
! - Get reduced field (CHAM_ELEM_S)
!
    if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
        call jeveuo(field_s(1:19)//'.CESD', 'L', jcesd)
        call jeveuo(field_s(1:19)//'.CESL', 'L', jcesl)
        call jeveuo(field_s(1:19)//'.CESV', 'L', jcesv)
        call jeveuo(field_s(1:19)//'.CESC', 'L', vk8=cesc)
        nb_cmp_maxi = zi(jcesd+4)
    end if
!
! - Number and name of components
!
    call getvtx(keyw_fact, 'NOM_CMP', iocc=i_keyw_fact, nbval=0, nbret=n1)
    if (n1 .lt. 0) then
        nb_cmp = -n1
        if ((nb_cmp .lt. 1) .or. (nb_cmp .gt. nb_para_maxi)) then
            vali(1) = nb_para_maxi
            vali(2) = nb_cmp
            call utmess('F', 'EXTRACTION_12', ni=2, vali=vali)
        end if

        call wkvect(list_cmp, 'V V K8', nb_cmp, vk8=v_list_cmp)
        call getvtx(keyw_fact, 'NOM_CMP', iocc=i_keyw_fact, nbval=nb_cmp, vect=v_list_cmp, &
                    nbret=iret)
        type_sele_cmp = 'NOM_CMP'
    else
        call getvtx(keyw_fact, 'NOM_VARI', iocc=i_keyw_fact, nbval=0, nbret=n1)
        ASSERT(n1 .lt. 0)
        ASSERT(field_type .eq. 'VARI_ELGA')
        nb_cmp = -n1
        call wkvect(list_cmp, 'V V K8', nb_elem*nb_cmp, vk8=v_list_cmp)
        call wkvect(list_vari, 'V V K16', nb_cmp, vk16=v_list_vari)
        call getvtx(keyw_fact, 'NOM_VARI', iocc=i_keyw_fact, nbval=nb_cmp, vect=v_list_vari, &
                    nbret=iret)
        call jeveuo(list_elem, 'L', vi=v_list_elem)
        call varinonu(model, compor, nb_elem, v_list_elem, nb_cmp, v_list_vari, v_list_cmp)
        type_sele_cmp = 'NOM_VARI'
    end if
!
! - Check components
!
    if (field_disc .eq. 'NOEU') then

        call exisd('CHAM_NO', field, iret)
        if (iret .ne. 1) call utmess('F', 'EXTRACTION_1', sk=field_type)
!
! ----- For nodes
!
        if (nb_node > 0) then
            call jeveuo(list_node, 'L', vi=v_list_node)
            do i_node = 1, nb_node
!
! --------- Current node
!
                node_nume = v_list_node(i_node)
                node_name = int_to_char8(node_nume)
                do i_cmp = 1, nb_cmp
                    cmp_name = v_list_cmp(i_cmp)
                    call posddl('CHAM_NO', field, node_name, cmp_name, nuno, &
                                nuddl)
                    if ((nuno .eq. 0) .or. (nuddl .eq. 0)) then
                        valk(1) = node_name
                        valk(2) = cmp_name
                        call utmess('F', 'EXTRACTION_20', nk=2, valk=valk)
                    end if
                end do
            end do
        end if
    else if (field_disc .eq. 'ELGA' .or. field_disc .eq. 'ELEM') then
!
! ----- For elements
!
        if (nb_elem > 0) then
            call jeveuo(list_elem, 'L', vi=v_list_elem)
            call jeveuo(list_poin, 'L', vi=v_list_poin)
            call jeveuo(list_spoi, 'L', vi=v_list_spoi)
            do i_elem = 1, nb_elem
!
! --------- Current element
!
                elem_nume = v_list_elem(i_elem)
                elem_name = int_to_char8(elem_nume)
!
! --------- Number of points/subpoints on current element
!
                nb_elem_poin = zi(jcesd+5+4*(elem_nume-1))
                nb_elem_spoi = zi(jcesd+5+4*(elem_nume-1)+1)
!
! --------- Check
!
                npi = nb_poin
                nspi = nb_spoi
                if (npi .gt. nb_elem_poin) npi = nb_elem_poin
                if (nspi .gt. nb_elem_spoi) nspi = nb_elem_spoi
!
                nb_cmp_maxi = zi(jcesd+4)
                do ipar = 1, nb_cmp
                    if (type_sele_cmp .eq. 'NOM_CMP') then
                        cmp_name = v_list_cmp(ipar)
                    elseif (type_sele_cmp .eq. 'NOM_VARI') then
                        cmp_name = v_list_cmp(nb_cmp*(i_elem-1)+ipar)
                    else
                        ASSERT(.false.)
                    end if
!
! ------------- For VARI_ELGA field
!
                    if (field_type(1:4) .eq. 'VARI') then
                        cmp_vari_name = cmp_name(2:8)//' '
                        call lxliis(cmp_vari_name, i_vari, iret)
                        if (iret .ne. 0) then
                            call utmess('F', 'EXTRACTION_22', sk=cmp_name)
                        end if
                    else
                        i_vari = 0
                    end if

                    if (field_type(1:4) .eq. 'VARI') then
                        i_cmp = i_vari
                    else
                        do i_cmp_maxi = 1, nb_cmp_maxi
                            if (cmp_name .eq. cesc(i_cmp_maxi)) then
                                i_cmp = i_cmp_maxi
                            end if
                        end do
                    end if
                    do ipi = 1, npi
                        num = v_list_poin(ipi)
                        ASSERT(num .ne. 0)
                        do ispi = 1, nspi
                            snum = v_list_spoi(ispi)
                            ASSERT(snum .ne. 0)
                            call cesexi('C', jcesd, jcesl, elem_nume, num, &
                                        snum, i_cmp, iad)
                            if (iad .eq. 0) then
                                valk(1) = elem_name
                                valk(2) = cmp_name
                                vali(1) = num
                                vali(2) = snum
                                call utmess('F', 'EXTRACTION_21', nk=2, valk=valk, ni=2, &
                                            vali=vali)
                            end if
                        end do
                    end do
                end do
            end do
        end if
    else
        ASSERT(.false.)
    end if
!
end subroutine
