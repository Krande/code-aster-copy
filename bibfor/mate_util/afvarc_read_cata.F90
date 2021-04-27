! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine afvarc_read_cata(varc_cata)
!
use Material_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
type(Mat_DS_VarcListCata), intent(out) :: varc_cata
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! Read list of variables from AFFE_MATERIAU catalog
!
! --------------------------------------------------------------------------------------------------
!
! Out varc_cata        : datastructure for catalog of external state variables
!
!   /!\ The caller is responsible to deallocate the 'varc_cata'
!       content with 'afvarc_free'.
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nmCmpMax = 9
    integer, parameter :: nbExteStatVari = 14
    integer :: iExteStatVari, nb_cmp, i_cmp
    character(len=8) :: varc_name, phys_para
    character(len=16) :: field_type
    character(len=8), parameter :: listExteStatVari(nbExteStatVari) = (/"TEMP    ", "GEOM    ",&
                                                                        "CORR    ", "IRRA    ",&
                                                                        "HYDR    ", "SECH    ",&
                                                                        "EPSA    ", "M_ACIER ",&
                                                                        "M_ZIRC  ", "NEUT1   ",&
                                                                        "NEUT2   ", "NEUT3   ",&
                                                                        "PTOT    ", "DIVU    "/)
    character(len=8), parameter :: listPhysQuantity(nbExteStatVari) = (/"TEMP_R  ", "GEOM_R  ",&
                                                                        "CORR_R  ", "IRRA_R  ",&
                                                                        "HYDR_R  ", "TEMP_R  ",&
                                                                        "EPSI_R  ", "VARI_R  ",&
                                                                        "VARI_R  ", "NEUT_R  ",&
                                                                        "NEUT_R  ", "NEUT_R  ",&
                                                                        "DEPL_R  ", "EPSI_R  "/)
    integer, parameter :: listNbCmp(nbExteStatVari) = (/ 7, 3,&
                                                         1, 1,&
                                                         1, 1,&
                                                         6, 9,&
                                                         5, 1,&
                                                         1, 1,&
                                                         1, 1/)
!
! --------------------------------------------------------------------------------------------------
!

! - Allocate
    varc_cata%nb_varc = nbExteStatVari
    allocate(varc_cata%list_cata_varc(nbExteStatVari))

! - Definition
    do iExteStatVari = 1, nbExteStatVari
! ----- Main properties
        varc_name = listExteStatVari(iExteStatVari)
        phys_para = listPhysQuantity(iExteStatVari)
        nb_cmp    = listNbCmp(iExteStatVari)
        ASSERT(nb_cmp .le. nmCmpMax)
        varc_cata%list_cata_varc(iExteStatVari)%name           = varc_name
        varc_cata%list_cata_varc(iExteStatVari)%type_phys_para = phys_para
        varc_cata%list_cata_varc(iExteStatVari)%nb_cmp         = nb_cmp

! ----- List of components
        allocate(varc_cata%list_cata_varc(iExteStatVari)%list_cmp(nb_cmp))
        if (varc_name .eq. 'TEMP') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "TEMP"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%phys_para_cmp = "TEMP_MIL"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%phys_para_cmp = "TEMP_INF"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(4)%phys_para_cmp = "TEMP_SUP"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(5)%phys_para_cmp = "DTX"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(6)%phys_para_cmp = "DTY"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(7)%phys_para_cmp = "DTZ"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "TEMP"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%varc_cmp = "TEMP_MIL"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%varc_cmp = "TEMP_INF"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(4)%varc_cmp = "TEMP_SUP"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(5)%varc_cmp = "DTX"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(6)%varc_cmp = "DTY"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(7)%varc_cmp = "DTZ"
            ASSERT(nb_cmp .eq. 7)
            field_type = varc_name

        elseif (varc_name .eq. 'NEUT1') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "X1"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "NEUT1"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'NEUT'

        elseif (varc_name .eq. 'NEUT2') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "X1"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "NEUT2"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'NEUT'

        elseif (varc_name .eq. 'NEUT3') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "X1"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "NEUT3"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'NEUT'

        elseif (varc_name .eq. 'GEOM') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "X"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%phys_para_cmp = "Y"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%phys_para_cmp = "Z"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "X"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%varc_cmp = "Y"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%varc_cmp = "Z"
            ASSERT(nb_cmp .eq. 3)
            field_type = 'GEOM'

        elseif (varc_name .eq. 'CORR') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "CORR"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "CORR"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'CORR'

        elseif (varc_name .eq. 'IRRA') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "IRRA"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "IRRA"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'IRRA'

        elseif (varc_name .eq. 'DIVU') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "DIVU"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "DIVU"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'DIVU'

        elseif (varc_name .eq. 'HYDR') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "HYDR"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "HYDR"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'HYDR_ELNO'

        elseif (varc_name .eq. 'SECH') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "TEMP"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "SECH"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'TEMP'

        elseif (varc_name .eq. 'PTOT') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "PTOT"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "PTOT"
            ASSERT(nb_cmp .eq. 1)
            field_type = 'PTOT'

        elseif (varc_name .eq. 'EPSA') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "EPXX"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%phys_para_cmp = "EPYY"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%phys_para_cmp = "EPZZ"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(4)%phys_para_cmp = "EPXY"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(5)%phys_para_cmp = "EPXZ"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(6)%phys_para_cmp = "EPYZ"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "EPSAXX"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%varc_cmp = "EPSAYY"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%varc_cmp = "EPSAZZ"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(4)%varc_cmp = "EPSAXY"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(5)%varc_cmp = "EPSAXZ"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(6)%varc_cmp = "EPSAYZ"
            ASSERT(nb_cmp .eq. 6)
            field_type = 'EPSA'

        elseif (varc_name .eq. 'M_ACIER') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "V1"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%phys_para_cmp = "V2"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%phys_para_cmp = "V3"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(4)%phys_para_cmp = "V4"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(5)%phys_para_cmp = "V5"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(6)%phys_para_cmp = "V6"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(7)%phys_para_cmp = "V7"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(8)%phys_para_cmp = "V8"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(9)%phys_para_cmp = "V9"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "PFERRITE"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%varc_cmp = "PPERLITE"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%varc_cmp = "PBAINITE"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(4)%varc_cmp = "PMARTENS"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(5)%varc_cmp = "PAUSTENI"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(6)%varc_cmp = "PCOLDSUM"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(7)%varc_cmp = "TAUSTE"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(8)%varc_cmp = "TRANSF"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(9)%varc_cmp = "TACIER"
            ASSERT(nb_cmp .eq. 9)
            field_type = 'META_ELNO'

        elseif (varc_name .eq. 'M_ZIRC') then
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%phys_para_cmp = "V1"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%phys_para_cmp = "V2"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%phys_para_cmp = "V3"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(4)%phys_para_cmp = "V4"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(5)%phys_para_cmp = "V5"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(1)%varc_cmp = "ALPHPUR"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(2)%varc_cmp = "ALPHBETA"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(3)%varc_cmp = "BETA"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(4)%varc_cmp = "TZIRC"
            varc_cata%list_cata_varc(iExteStatVari)%list_cmp(5)%varc_cmp = "TEMPS"
            ASSERT(nb_cmp .eq. 5)
            field_type = 'META_ELNO'
        else
            ASSERT(ASTER_FALSE)
        endif

        varc_cata%list_cata_varc(iExteStatVari)%field_type_def = field_type

    end do
!
    if (.false.) then
        do iExteStatVari = 1, nbExteStatVari
            write(6,*) 'Variable de commande :', iExteStatVari
            write(6,*) '> Nom      :', varc_cata%list_cata_varc(iExteStatVari)%name
            write(6,*) '> GRANDEUR :', varc_cata%list_cata_varc(iExteStatVari)%type_phys_para
            write(6,*) '> Default  :', varc_cata%list_cata_varc(iExteStatVari)%field_type_def
            write(6,*) '> NB_CMP   :', varc_cata%list_cata_varc(iExteStatVari)%nb_cmp
            do i_cmp = 1, varc_cata%list_cata_varc(iExteStatVari)%nb_cmp
                write(6,*) '> Nombre de composantes :', i_cmp
                write(6,*) '>> CMP_GD   :',&
                  varc_cata%list_cata_varc(iExteStatVari)%list_cmp(i_cmp)%phys_para_cmp
                write(6,*) '>> CMP_VARC :',&
                  varc_cata%list_cata_varc(iExteStatVari)%list_cmp(i_cmp)%varc_cmp
            end do
        end do
         write(6,*) 'Nombre total de variables de commande :', varc_cata%nb_varc
    endif
!
end subroutine
