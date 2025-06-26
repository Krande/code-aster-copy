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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine cgExportTableG(cgField, cgTheta, cgTable, cgStat)
!
    use calcG_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/titre.h"
#include "jeveux.h"
!
    type(CalcG_field), intent(in) :: cgField
    type(CalcG_theta), intent(in) :: cgTheta
    type(CalcG_table), intent(in) :: cgTable
    type(CalcG_stat), intent(inout) :: cgStat
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!     Export table container
!
! --------------------------------------------------------------------------------------------------
!
!
    integer(kind=8), parameter :: nb_objet = 2
    character(len=16) :: obje_name(nb_objet)
    character(len=24) :: obje_sdname(nb_objet)
!
    integer(kind=8), parameter :: nbpara = 3
    character(len=19), parameter :: nompar(nbpara) = (/ &
                                    'NOM_OBJET ', 'TYPE_OBJET', 'NOM_SD    '/)
    character(len=19), parameter :: typpar(nbpara) = (/ &
                                    'K16', 'K16', 'K24'/)
!
    integer(kind=8), parameter :: l_nb_obje = 2
    character(len=16), parameter :: l_obje_name(l_nb_obje) = (/ &
                                    'TABLE_G         ', 'CHAM_THETA      '/)
    character(len=16), parameter :: l_obje_type(l_nb_obje) = (/ &
                                    'TABLE_SDASTER   ', 'CHAM_NO_SDASTER '/)
    character(len=80), parameter :: def_title = &
                            'CALCUL DES FACTEURS D''INTENSITE DES CONTRAINTES PAR LA METHODE CALC_G'
    integer(kind=8), parameter :: l_def_title = 70
!
    integer(kind=8) :: i_l_obj, i_obj
    character(len=24) :: vk(3)
    character(len=16) :: obje_type
    aster_logical, parameter :: debug = ASTER_FALSE
    real(kind=8) :: start, finish
!
    call cpu_time(start)
!
    call jemarq()
!
! --- Create table_container to store (calc_g and cham_theta)
!
    obje_name(1) = "TABLE_G"
    obje_sdname(1) = cgTable%table_g
!
    obje_name(2) = "CHAM_THETA"
    obje_sdname(2) = cgTheta%theta_factors
!
! - Create new table_container
!
    call detrsd('TABLE_CONTAINER', cgField%table_out)
    call tbcrsd(cgField%table_out, 'G')
    call tbajpa(cgField%table_out, nbpara, nompar, typpar)
!
! - Loop on objects to add new one
!
    if (debug) then
        print *, "Create table_container in CALC_G"
        print *, "Number of object: ", nb_objet
    end if
!
    do i_obj = 1, nb_objet
!
! ----- Find the type of object
!
        obje_type = ' '
        do i_l_obj = 1, l_nb_obje
            if (l_obje_name(i_l_obj) .eq. obje_name(i_obj)) then
                obje_type = l_obje_type(i_l_obj)
                exit
            end if
        end do
        ASSERT(obje_type .ne. ' ')
!
! ----- Add object (new line)
!
        if (debug) then
            print *, "OBJET NAME", i_obj, " : ", obje_name(i_obj)
            print *, "OBJET TYPE", i_obj, " : ", obje_type
            print *, "SD NAME   ", i_obj, " : ", obje_sdname(i_obj)
        end if
!
        vk(1) = obje_name(i_obj)
        vk(2) = obje_type
        vk(3) = obje_sdname(i_obj)
!
        call tbajli(cgField%table_out, nbpara, nompar, &
                    [0], [0.d0], [dcmplx(0., 0.)], vk, 0)
    end do
!
    call jedema()
!
    call titre(defTitle=def_title, lDefTitle=l_def_title)
!
    call cpu_time(finish)
    cgStat%cgExpTabl = cgStat%cgExpTabl+finish-start
!
end subroutine
