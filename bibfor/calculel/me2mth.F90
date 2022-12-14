! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine me2mth(model_, nb_load, list_name_, cara_elem_,&
                  time_ , temp_  , vect_elem_)
!
implicit none
!
#include "asterfort/assert.h"
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/megeom.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/load_neut_excl.h"
!
!
    character(len=*), intent(in) :: model_
    character(len=*), intent(in) :: temp_
    character(len=*), intent(in) :: cara_elem_
    character(len=*), intent(in) :: vect_elem_
    character(len=*), intent(in) :: time_
    character(len=*), intent(in) :: list_name_(*)
    integer, intent(in) :: nb_load
!
! --------------------------------------------------------------------------------------------------
!
!         CALCUL DE TOUS LES SECONDS MEMBRES ELEMENTAIRES PROVENANT
!         DES CHARGES_THERMIQUES
!
!     ENTREES:
!
!     LES NOMS QUI SUIVENT SONT LES PREFIXES UTILISATEUR K8:
!        MODELZ : NOM DU MODELE
!        NCHAR  : NOMBRE DE CHARGES
!        LCHAR  : LISTE DES CHARGES
!        CARAZ  : CARAC_ELEM
!        TIMEZ  : CHAMPS DE TEMPSR
!        CHTNZ  : CHAM_NO DE TEMPERATURE A L'INSTANT TN
!        VECELZ : NOM DU VEC_ELE (N RESUELEM) PRODUIT
!                 SI VECEL EXISTE DEJA, ON LE DETRUIT.
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: model, cara_elem
    character(len=8) :: load_name
    character(len=19) :: vect_elem
    character(len=24) :: time
    aster_logical :: lfonc
    character(len=8) :: lpain(5), lpaout(1), k8bid
    character(len=16) :: option
    character(len=24) :: lchin(5), lchout(1), ligrmo, ligrch
    character(len=24) :: chgeom
    integer :: i_load, iret, ilires
    character(len=8), pointer :: v_nomo(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    model     = model_
    cara_elem = cara_elem_
    vect_elem = vect_elem_
    time      = time_
    call megeom(model, chgeom)
!
    call jeexin(vect_elem//'.RERR', iret)
    if (iret .gt. 0) then
        call jedetr(vect_elem//'.RERR')
        call jedetr(vect_elem//'.RELR')
    endif
    call memare('G', vect_elem, model, ' ', cara_elem,&
                'CHAR_THER')
!
    lpaout(1) = 'PVECTTR'
    lchout(1) = vect_elem(1:8)//'.VE000'
    ilires = 0
!
!     BOUCLE SUR LES CHARGES POUR CALCULER :
!         ( CHAR_THER_TEXT_F , ISO_FACE ) SUR LE MODELE
!         ( CHAR_THER_FLUN_F , ISO_FACE ) SUR LE MODELE
!         ( CHAR_THER_SOUR_F , ISO      ) SUR LE MODELE
!         ( CHAR_THER_SOUR_F , SOURCE_NO) SUR LE LIGREL(CHARGE)
!         ( THER_DDLI_F      , CAL_TI   ) SUR LE LIGREL(CHARGE)
!
!
    if (nb_load .ne. 0) then
!
! ----- Exclusion of some loads
!
        call load_neut_excl('CALC_VECT_ELEM',list_nbload_ = nb_load, list_name_ = list_name_)
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PTEMPSR'
        lchin(2) = time
        if (model .ne. '        ') then
            ligrmo = model//'.MODELE'
        else
            load_name = list_name_(1)
            call jeveuo(load_name//'.CHTH      .NOMO', 'L', vk8=v_nomo)
            ligrmo = v_nomo(1)//'.MODELE'
        endif
        do i_load = 1, nb_load
            load_name = list_name_(i_load)
            call dismoi('TYPE_CHARGE', load_name, 'CHARGE', repk=k8bid)
            if (k8bid(5:7) .eq. '_FO') then
                lfonc = .true.
            else
                lfonc = .false.
            endif
!
            ligrch = load_name//'.CHTH.LIGREL'

!
!  =====================================================================
!           --  ( CHAR_THER_TEXT_F , ISO_FACE ) SUR LE MODELE
            call exisd('CHAMP_GD', ligrch(1:13)//'.COEFH', iret)
            if (iret .ne. 0) then
                if (lfonc) then
                    option = 'CHAR_THER_TEXT_F'
                    lpain(3) = 'PT_EXTF'
                    lpain(4) = 'PCOEFHF'
                else
                    option = 'CHAR_THER_TEXT_R'
                    lpain(3) = 'PT_EXTR'
                    lpain(4) = 'PCOEFHR'
                endif
                lchin(3) = ligrch(1:13)//'.T_EXT     '
                lchin(4) = ligrch(1:13)//'.COEFH     '
                lpain(5) = 'PTEMPER'
                lchin(5) = temp_
                ilires = ilires + 1
                call codent(ilires, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrmo, 5, lchin,&
                            lpain, 1, lchout, lpaout, 'G',&
                            'OUI')
                call reajre(vect_elem, lchout(1), 'G')
            endif
!  =====================================================================
!           --  ( CHAR_THER_FLUN_F , ISO_FACE ) SUR LE MODELE
            call exisd('CHAMP_GD', ligrch(1:13)//'.FLURE', iret)
            if (iret .ne. 0) then
                if (lfonc) then
                    option = 'CHAR_THER_FLUN_F'
                    lpain(3) = 'PFLUXNF'
                else
                    option = 'CHAR_THER_FLUN_R'
                    lpain(3) = 'PFLUXNR'
                endif
                lchin(3) = ligrch(1:13)//'.FLURE     '
                ilires = ilires + 1
                call codent(ilires, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrmo, 3, lchin,&
                            lpain, 1, lchout, lpaout, 'G',&
                            'OUI')
                call reajre(vect_elem, lchout(1), 'G')
            endif
!  =====================================================================
!           --  ( CHAR_THER_FLUX_  , ISO_FACE ) SUR LE MODELE
            call exisd('CHAMP_GD', ligrch(1:13)//'.FLUR2', iret)
            if (iret .ne. 0) then
                if (lfonc) then
                    option = 'CHAR_THER_FLUX_F'
                    lpain(3) = 'PFLUXVF'
                else
                    option = 'CHAR_THER_FLUX_R'
                    lpain(3) = 'PFLUXVR'
                endif
                lchin(3) = ligrch(1:13)//'.FLUR2     '
                ilires = ilires + 1
                call codent(ilires, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrmo, 3, lchin,&
                            lpain, 1, lchout, lpaout, 'G',&
                            'OUI')
                call reajre(vect_elem, lchout(1), 'G')
            endif
!  =====================================================================
!           --   ( CHAR_THER_SOUR_F , ISO    )  SUR LE MODELE
            call exisd('CHAMP_GD', ligrch(1:13)//'.SOURE', iret)
            if (iret .ne. 0) then
                if (lfonc) then
                    option = 'CHAR_THER_SOUR_F'
                    lpain(3) = 'PSOURCF'
                else
                    option = 'CHAR_THER_SOUR_R'
                    lpain(3) = 'PSOURCR'
                endif
                lchin(3) = ligrch(1:13)//'.SOURE     '
                ilires = ilires + 1
                call codent(ilires, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrmo, 3, lchin,&
                            lpain, 1, lchout, lpaout, 'G',&
                            'OUI')
                call reajre(vect_elem, lchout(1), 'G')
            endif
!  =====================================================================
!           --   ( CHAR_THER_SOUR_R , CALC    )  SUR LE MODELE
            call exisd('CHAMP_GD', ligrch(1:13)//'.SOURC', iret)
            if (iret .ne. 0) then
                if (lfonc) then
                    ASSERT(.false.)
                else
                    option = 'CHAR_THER_SOUR_R'
                    lpain(3) = 'PSOURCR'
                endif
                lchin(3) = ligrch(1:13)//'.SOURC     '
                ilires = ilires + 1
                call codent(ilires, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrmo, 3, lchin,&
                            lpain, 1, lchout, lpaout, 'G',&
                            'OUI')
                call reajre(vect_elem, lchout(1), 'G')
            endif
!  =====================================================================
!           --   ( CHAR_THER_GRAI_  ,  ISO_VOLU  ) SUR LE   MODELE
            call exisd('CHAMP_GD', ligrch(1:13)//'.GRAIN', iret)
            if (iret .ne. 0) then
                if (lfonc) then
                    option = 'CHAR_THER_GRAI_F'
                    lpain(3) = 'PGRAINF'
                else
                    option = 'CHAR_THER_GRAI_R'
                    lpain(3) = 'PGRAINR'
                endif
                lchin(3) = ligrch(1:13)//'.GRAIN     '
                ilires = ilires + 1
                call codent(ilires, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrmo, 3, lchin,&
                            lpain, 1, lchout, lpaout, 'G',&
                            'OUI')
                call reajre(vect_elem, lchout(1), 'G')
            endif
!  =====================================================================
!           --   ( THER_DDLI_F    , CAL_TI   )  SUR LE LIGREL(CHARGE)
            call exisd('CHAMP_GD', ligrch(1:13)//'.CIMPO', iret)
            if (iret .ne. 0) then
                if (lfonc) then
                    option = 'THER_DDLI_F'
                    lpain(3) = 'PDDLIMF'
                else
                    option = 'THER_DDLI_R'
                    lpain(3) = 'PDDLIMR'
                endif
                lchin(3) = ligrch(1:13)//'.CIMPO     '
                ilires = ilires + 1
                call codent(ilires, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrch, 3, lchin,&
                            lpain, 1, lchout, lpaout, 'G',&
                            'OUI')
                call reajre(vect_elem, lchout(1), 'G')
            endif
        end do
    endif
!
    call jedema()
end subroutine
