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

subroutine ccchuc_norm(norm, model, name_gd, field_in, type_field_in, &
                       field_out)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/calc_coor_elga.h"
#include "asterfort/calc_norm_coef.h"
#include "asterfort/calc_norm_elem.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/chpchd.h"
#include "asterfort/chsut1.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/megeom.h"
#include "asterfort/nopar2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!
    character(len=16), intent(in) :: norm
    character(len=8), intent(in) :: model
    character(len=8), intent(in) :: name_gd
    character(len=19), intent(in) :: field_in
    character(len=4), intent(in) :: type_field_in
    character(len=19), intent(in) :: field_out
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP - CHAM_UTIL - NORME
!
! Compute NORME
!
! --------------------------------------------------------------------------------------------------
!
! In  norm          : type of norm
! In  model         : name of model
! In  name_gd       : name of <GRANDEUR> of input field
! In  type_field_in : type of input field
! In  field_in      : name of <CHAM_ELEM> input field FROM which extract values
! In  field_out     : name of <CHAM_ELEM> output field IN which compute values
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_cmp_max
    parameter(nb_cmp_max=30)
    integer(kind=8) :: iexist
    character(len=19) :: modelLigrel, celmod
    character(len=19) :: field_in_s, field_neut_s, field_neut, field_neut_mod
    integer(kind=8) ::  jchsc
    integer(kind=8) :: nb_elem, nb_cmp, nb_cmp_act
    character(len=24) :: list_cmp, list_cmp_neut, valk(3)
    integer(kind=8) :: j_liscmp_in, j_liscmp_ne
    character(len=24) :: chcoef, chgaus, chgeom, chcalc
    integer(kind=8) :: nb_coef_user
    real(kind=8) :: coef_user(1)
    character(len=4) :: ki
    integer(kind=8) :: icmp, nncp, iret, ibid
    character(len=16) :: option
    character(len=8) :: nopar
    integer(kind=8), pointer :: cesd(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    field_in_s = '&&CCCHUC_NORM.FIELS'
    field_neut_s = '&&CCCHUC_NORM.NEUTS'
    field_neut = '&&CCCHUC_NORM.NEUTR'
    field_neut_mod = '&&CCCHUC_NORM.NEUTM'
    chcoef = '&&CCCHUC_NORM.CHCOEF'
    chgaus = '&&CCCHUC_NORM.CHGAUS'
    chcalc = '&&CCCHUC_NORM.CHCALC'
    list_cmp_neut = '&&CCCHUC_NORM.CMPN'
    nb_coef_user = 0
    coef_user(1) = 0.d0
!
! - Compute <CARTE> with informations on Gauss points
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call exisd('CHAMP', chgaus, iexist)
    if (iexist .eq. 0) then
        call megeom(model, chgeom)
        call calc_coor_elga(model, modelLigrel, chgeom, chgaus)
    end if
!
! - Create <CHAM_ELEM_S> from input field
!
    call celces(field_in, 'V', field_in_s)
    call jeveuo(field_in_s//'.CESD', 'L', vi=cesd)
    call jeveuo(field_in_s//'.CESC', 'L', jchsc)
    nb_elem = cesd(1)
    nb_cmp = cesd(2)
    list_cmp = field_in_s//'.CESC'
!
! - <NEUT_R> components
!
    call wkvect(list_cmp_neut, 'V V K8', nb_cmp, j_liscmp_ne)
    do icmp = 1, nb_cmp
        call codent(icmp, 'G', ki)
        zk8(j_liscmp_ne-1+icmp) = 'X'//ki(1:len(ki))
    end do
!
! - Construction of <CARTE> of <NEUT_R> by selection of components
!
    call calc_norm_coef(model, name_gd, nb_cmp_max, nb_cmp, norm, &
                        'NORM', list_cmp, nb_coef_user, coef_user, chcoef, &
                        chcalc, nb_cmp_act)
!
! - Transform input field in NEUT_R
!
    call jeveuo(list_cmp, 'L', j_liscmp_in)
    call chsut1(field_in_s, 'NEUT_R', nb_cmp, zk8(j_liscmp_in), zk8(j_liscmp_ne), &
                'V', field_neut_s)
!
! - Convert CHAMELEM_S field to CHAMELEM field
!
    if (type_field_in .eq. 'ELNO') then
        option = 'TOU_INI_ELNO'
    else if (type_field_in .eq. 'ELGA') then
        option = 'TOU_INI_ELGA'
    else
        ASSERT(.false.)
    end if
    call nopar2(option, 'NEUT_R', 'OUT', nopar)
    call cescel(field_neut_s, modelLigrel, option, nopar, 'OUI', &
                nncp, 'V', field_neut, 'F', iret)
    ASSERT(iret .eq. 0)
!
! - Change type of field
!
    if (norm .eq. 'L2') then
        option = 'NORME_L2'
    else if (norm .eq. 'FROBENIUS') then
        option = 'NORME_FROB'
    else
        ASSERT(.false.)
    end if
    if (type_field_in .eq. 'ELGA') then
        field_neut_mod = field_neut
    else
        nopar = 'PCHAMPG'
        celmod = '&&PENORM.CELMOD'
        call alchml(modelLigrel, option, nopar, 'V', celmod, &
                    ibid, ' ')
        if (ibid .ne. 0) then
            valk(1) = modelLigrel
            valk(2) = nopar
            valk(3) = option
            call utmess('F', 'UTILITAI3_23', nk=3, valk=valk)
        end if
        call chpchd(field_neut, 'ELGA', celmod, 'OUI', 'V', &
                    field_neut_mod, model)
        call detrsd('CHAMP', celmod)
    end if
!
! - Compute Norm (integration on finite element)
!
    call calc_norm_elem(norm, modelLigrel, chcoef, chgaus, chcalc, &
                        field_neut_mod, field_out)
!
    call jedetr(chcoef)
    call jedetr(list_cmp_neut)
    call detrsd('CHAM_ELEM_S', field_in_s)
    call detrsd('CHAM_ELEM_S', field_neut_s)
    call detrsd('CHAM_ELEM', field_neut)
    call detrsd('CHAM_ELEM', field_neut_mod)
!
end subroutine
