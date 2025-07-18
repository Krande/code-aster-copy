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

subroutine nmextp(keyw_fact, i_keyw_fact, field_type, field_disc, field, &
                  field_s, list_poin, list_spoi, nb_poin, nb_spoi, &
                  type_extr_elem)
!
    implicit none
!
#include "asterfort/celces.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/exisd.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sdmpic.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=16), intent(in) :: keyw_fact
    integer(kind=8), intent(in) :: i_keyw_fact
    character(len=19), intent(in) :: field
    character(len=24), intent(in) :: field_type
    character(len=4), intent(in) :: field_disc
    character(len=24), intent(in) :: field_s
    character(len=24), intent(in) :: list_poin
    character(len=24), intent(in) :: list_spoi
    integer(kind=8), intent(out) :: nb_poin
    integer(kind=8), intent(out) :: nb_spoi
    character(len=8), intent(out) :: type_extr_elem
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Extraction (OBSERVATION/SUIVI_DDL) utilities
!
! Get topology (point and subpoints) and type of extraction for element
!
! --------------------------------------------------------------------------------------------------
!
! In  keyw_fact        : factor keyword to read extraction parameters
! In  i_keyw_fact      : index of keyword to read extraction parameters
! In  field            : name of field
! In  field_disc       : type of discretization (ELGA or ELEM)
! In  field_s          : name of reduced field (CHAM_ELEM_S)
! In  field_type       : type of field (name in results datastructure)
! In  list_poin        : name of object contains list of points (Gauss)
! Out nb_poin          : number of points (Gauss)
! In  list_spoi        : name of object contains list of subpoints
! Out nb_spoi          : number of subpoints
! Out type_extr_elem   : type of extraction by element
!                'MIN'  VALEUR MINI SUR TOUS LES POINTS DE GAUSS
!                'MAX'  VALEUR MAXI SUR TOUS LES POINTS DE GAUSS
!                'MOY'  VALEUR MOYENNE SUR TOUS LES POINTS DE GAUSS
!                'VALE' VALEUR SUR POINT/SOUS_POINTS DONNES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_poin, i_spoi
    integer(kind=8) :: n1, n2, n3, iret
    integer(kind=8) :: nb_poin_maxi, nb_spoi_maxi
    integer(kind=8), pointer :: cesd(:) => null()
    integer(kind=8), pointer :: v_list_poin(:) => null()
    integer(kind=8), pointer :: v_list_spoi(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nb_poin = 0
    nb_spoi = 0
!
! - Conversion to reduced field (CHAM_ELEM_S)
!
    call exisd('CHAM_ELEM', field, iret)
    if (iret .ne. 1) call utmess('F', 'EXTRACTION_1', sk=field_type)
    call exisd('CHAM_ELEM_S', field_s, iret)
    if (iret .eq. 0) then
        call sdmpic('CHAM_ELEM', field)
        call celces(field, 'V', field_s)
    end if

    call jeveuo(field_s(1:19)//'.CESD', 'L', vi=cesd)
!
! - Type of extraction on element
!
    if (field_disc .eq. 'ELGA') then
        call getvtx(keyw_fact, 'EVAL_ELGA', iocc=i_keyw_fact, scal=type_extr_elem, nbret=n1)
        if (n1 .eq. 0) then
            type_extr_elem = 'VALE'
            call utmess('A', 'EXTRACTION_6', sk=field_type)
        end if
    else
        type_extr_elem = 'VALE'
    end if
!
! - Max number of points/subpoint for this field
!
    nb_poin_maxi = cesd(3)
    nb_spoi_maxi = cesd(4)
!
! - Number of points/subpoint
!
    if (field_disc .eq. 'ELGA') then
        if (type_extr_elem .eq. 'VALE') then
            call getvis(keyw_fact, 'POINT', iocc=i_keyw_fact, nbval=0, nbret=n2)
            call getvis(keyw_fact, 'SOUS_POINT', iocc=i_keyw_fact, nbval=0, nbret=n3)
            if (n2 .eq. 0) then
                call utmess('F', 'EXTRACTION_7', sk=field_type)
            end if
            nb_poin = -n2
            if ((n2 .ne. 0) .and. (n3 .eq. 0)) then
                nb_spoi = nb_spoi_maxi
            else
                nb_spoi = -n3
            end if
        else
            nb_poin = nb_poin_maxi
            nb_spoi = nb_spoi_maxi
        end if
    else
        nb_poin = 1
        nb_spoi = 1
    end if
!
! - Protection
!
    if (nb_poin .gt. nb_poin_maxi) nb_poin = nb_poin_maxi
    if (nb_spoi .gt. nb_spoi_maxi) nb_spoi = nb_spoi_maxi
!
! - Create lists
!
    call wkvect(list_poin, 'V V I', nb_poin, vi=v_list_poin)
    if (nb_spoi .ne. 0) then
        call wkvect(list_spoi, 'V V I', nb_spoi, vi=v_list_spoi)
    end if
!
! - Set lists
!
    if (type_extr_elem .eq. 'VALE' .and. field_disc .eq. 'ELGA') then
        call getvis(keyw_fact, 'POINT', iocc=i_keyw_fact, nbval=nb_poin, vect=v_list_poin)
        if (nb_spoi .ne. 0) then
            call getvis(keyw_fact, 'SOUS_POINT', iocc=i_keyw_fact, nbval=nb_spoi, &
                        vect=v_list_spoi, nbret=n3)
            if (n3 .eq. 0) then
                do i_spoi = 1, nb_spoi
                    v_list_spoi(i_spoi) = i_spoi
                end do
            end if
        end if
    else
        do i_poin = 1, nb_poin
            v_list_poin(i_poin) = i_poin
        end do
        do i_spoi = 1, nb_spoi
            v_list_spoi(i_spoi) = i_spoi
        end do
    end if
!
end subroutine
