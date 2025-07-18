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
subroutine romMultiParaRead(ds_multipara)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/getexm.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jelira.h"
#include "asterfort/utmess.h"
#include "asterfort/romMultiCoefRead.h"
#include "asterfort/romVariParaRead.h"
#include "asterfort/romFieldGetInfo.h"
#include "asterfort/as_allocate.h"
!
    type(ROM_DS_MultiPara), intent(inout) :: ds_multipara
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Read data for multiparametric problems
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_multipara     : datastructure for multiparametric problems
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nb_matr, nb_vect
    integer(kind=8) :: i_matr, i_vect, i_vari_para, nbret, nb_vari_coef, nb_vari_para
    character(len=1)  :: matr_type, vect_type, ktyp, matr_elem_type, vect_elem_type
    character(len=16) :: keywfact, type_vari_coef
    character(len=8) :: matr_asse, vect_asse, gran_name, model
    character(len=19) :: nume_equa
    character(len=24) :: field_name, vect_assez
    type(ROM_DS_MultiCoef) :: ds_multicoef
    type(ROM_DS_VariPara) :: ds_varipara
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM2_20')
    end if
!
! - List of matrix
!
    if (niv .ge. 2) then
        call utmess('I', 'ROM2_25')
    end if
    nb_matr = 0
    nb_vect = 0
    keywfact = 'MATR_ASSE'
    matr_type = 'R'
    call getfac(keywfact, nb_matr)
    ASSERT(nb_matr .gt. 0)
    AS_ALLOCATE(vk8=ds_multipara%matr_name, size=nb_matr)
    AS_ALLOCATE(vk8=ds_multipara%matr_type, size=nb_matr)
    allocate (ds_multipara%matr_coef(nb_matr))
    do i_matr = 1, nb_matr
        matr_asse = ' '
        matr_elem_type = ' '
        if (getexm(keywfact, 'MATRICE') .eq. 1) then
            call getvid(keywfact, 'MATRICE', iocc=i_matr, scal=matr_asse, nbret=nbret)
            ASSERT(nbret .gt. 0)
            call jelira(matr_asse//'           .VALM', 'TYPE', cval=ktyp)
            if (ktyp .eq. 'C') then
                matr_elem_type = 'C'
                matr_type = 'C'
            elseif (ktyp .eq. 'R') then
                matr_elem_type = 'R'
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        call romMultiCoefRead(ds_multicoef, keywfact, i_matr)
        ds_multipara%matr_name(i_matr) = matr_asse
        ds_multipara%matr_type(i_matr) = matr_elem_type
        ds_multipara%matr_coef(i_matr) = ds_multicoef
    end do
    ds_multipara%nb_matr = nb_matr
    AS_ALLOCATE(vk24=ds_multipara%matr_mode_curr, size=nb_matr)
    AS_ALLOCATE(vk24=ds_multipara%prod_matr_mode, size=nb_matr)
    AS_ALLOCATE(vk24=ds_multipara%matr_redu, size=nb_matr)
!
! - List of vectors
!
    if (niv .ge. 2) then
        call utmess('I', 'ROM2_23')
    end if
    keywfact = 'VECT_ASSE'
    vect_type = 'R'
    call getfac(keywfact, nb_vect)
    ASSERT(nb_vect .gt. 0)
    AS_ALLOCATE(vk8=ds_multipara%vect_name, size=nb_vect)
    AS_ALLOCATE(vk8=ds_multipara%vect_type, size=nb_vect)
    allocate (ds_multipara%vect_coef(nb_vect))
    do i_vect = 1, nb_vect
        vect_asse = ' '
        vect_elem_type = ' '
        if (getexm(keywfact, 'VECTEUR') .eq. 1) then
            call getvid(keywfact, 'VECTEUR', iocc=i_vect, scal=vect_asse, nbret=nbret)
            ASSERT(nbret .gt. 0)
            call dismoi('NOM_GD', vect_asse, 'CHAM_NO', repk=gran_name)
            if (gran_name .eq. 'DEPL_C') then
                vect_elem_type = 'C'
                vect_type = 'C'
            elseif (gran_name .eq. 'DEPL_R') then
                vect_elem_type = 'R'
            else
                ASSERT(.false.)
            end if
        end if
        call romMultiCoefRead(ds_multicoef, keywfact, i_vect)
        ds_multipara%vect_name(i_vect) = vect_asse
        ds_multipara%vect_type(i_vect) = vect_type
        ds_multipara%vect_coef(i_vect) = ds_multicoef
    end do
    ds_multipara%nb_vect = nb_vect
    AS_ALLOCATE(vk24=ds_multipara%vect_redu, size=nb_vect)
!
! - Get informations from field
!
    vect_assez = vect_asse
    call dismoi('NUME_EQUA', vect_assez, 'CHAM_NO', repk=nume_equa)
    call dismoi('NOM_MODELE', nume_equa, 'NUME_EQUA', repk=model)
    field_name = 'DEPL'
    call romFieldGetInfo(model, field_name, vect_assez, ds_multipara%field, l_chck_=ASTER_FALSE)
!
! - Set system type
!
    ds_multipara%syst_type = 'R'
    if (vect_type .eq. 'C' .or. matr_type .eq. 'C') then
        ds_multipara%syst_type = 'C'
    end if
!
! - Read data for variations of multiparametric problems
!
    call getvis(' ', 'NB_VARI_COEF', scal=nb_vari_coef, nbret=nbret)
    ASSERT(nbret .eq. 1)
    ds_multipara%nb_vari_coef = nb_vari_coef
    call getvtx(' ', 'TYPE_VARI_COEF', scal=type_vari_coef, nbret=nbret)
    ds_multipara%type_vari_coef = type_vari_coef
    keywfact = 'VARI_PARA'
    call getfac(keywfact, nb_vari_para)
    if (niv .ge. 2) then
        if (nb_vari_para .eq. 0) then
            call utmess('I', 'ROM2_38')
        else
            call utmess('I', 'ROM2_24')
        end if
    end if
    ds_multipara%nb_vari_para = nb_vari_para
    allocate (ds_multipara%vari_para(nb_vari_para))
    do i_vari_para = 1, nb_vari_para
        call romVariParaRead(ds_varipara, keywfact, i_vari_para)
        ds_multipara%vari_para(i_vari_para) = ds_varipara
    end do
!
end subroutine
