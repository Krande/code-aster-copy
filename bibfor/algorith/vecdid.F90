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

subroutine vecdid(model, list_load, nume_dof, veelem, veasse)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assvec.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/exisd.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/nmchex.h"
#include "asterfort/vemare.h"
#include "asterfort/reajre.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=24), intent(in) :: model
    character(len=19), intent(in) :: list_load
    character(len=24), intent(in) :: nume_dof
    character(len=19), intent(in) :: veelem(*)
    character(len=19), intent(in) :: veasse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Computation
!
! Compute elementary vector and assemble for DIDI loads
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  list_load        : name of datastructure for list of loads
! In  nume_dof         : name of numbering (NUME_DDL)
! In  veelem           : hat variable for elementary vectors
! In  veasse           : hat variable for vectors
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=1, nbin=3)
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    integer(kind=8) :: iret
    integer(kind=8) :: nchar, nbres, icha
    real(kind=8) :: alpha
    character(len=8) :: nomcha
    character(len=19) :: vect_elem, vect_asse, disp_didi
    character(len=16) :: option
    character(len=1) :: base
    character(len=24) :: masque
    character(len=24) :: ligrch, chalph
    aster_logical :: l_didi
    integer(kind=8), pointer :: infc(:) => null()
    character(len=24), pointer :: lcha(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! --- LISTE DES CHARGES
!
    call jeexin(list_load(1:19)//'.LCHA', iret)
    if (iret .eq. 0) goto 99
!
    call jelira(list_load(1:19)//'.LCHA', 'LONMAX', nchar)
    call jeveuo(list_load(1:19)//'.LCHA', 'L', vk24=lcha)
    call jeveuo(list_load(1:19)//'.INFC', 'L', vi=infc)
!
! --- VERIF SI CHARGE DE TYPE DIRICHLET DIFFERENTIEL
!
    l_didi = .false.
    do icha = 1, nchar
        if (infc(icha+1) .gt. 0 .or. infc(1+3*nchar+2+icha) .eq. 0) then
            l_didi = .true.
        end if
    end do
    if (.not. l_didi) goto 99
!
! --- INITIALISATIONS
!
    call nmchex(veelem, 'VEELEM', 'CNDIDI', vect_elem)
    call nmchex(veasse, 'VEASSE', 'CNDIDI', vect_asse)
    disp_didi = '&&CNCHAR.DIDI'

    base = 'V'
    option = 'MECA_BU_R'
!
! --- INITIALISATION DES CHAMPS POUR CALCUL
!
    call inical(nbin, lpain, lchin, nbout, lpaout, &
                lchout)
!
! --- CONSTRUCTION DU VECTEUR BDIDI.UREF
!
! REM : LE TERME BT.LAMBDA EST EGALEMENT CALCULE. IL EST NUL CAR A CE
!       STADE, LES LAMBDAS SONT NULS.
!
!
! --- ALLOCATION DE LA CARTE DU CONDITIONNEMENT DES LAGRANGES
! REM : A CE STADE, ON FIXE LE COND A 1
!
    alpha = 1.d0
    chalph = '&&VEBUME.CH_NEUT_R'
    call mecact('V', chalph, 'MODELE', model, 'NEUT_R  ', &
                ncmp=1, nomcmp='X1', sr=alpha)
!
! --- PREPARATION DES VECT_ELEM
!
    call jeexin(vect_elem(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        call vemare('V', vect_elem, model(1:8))
    end if
    call jedetr(vect_elem(1:19)//'.RELR')
    call reajre(vect_elem, ' ', 'V')
    masque = vect_elem(1:19)//'.VEXX'
!
! --- BOUCLE SUR LES CHARGES DE TYPE DIRICHLET DIFFERENTIEL
!
    nbres = 0
    do icha = 1, nchar
!
! --- VERIF SI CHARGE DE TYPE DIRICHLET DIFFERENTIEL
!
        if (infc(icha+1) .le. 0 .or. infc(1+3*nchar+2+icha) .eq. 0) then
            cycle
        end if
        nomcha = lcha(icha) (1:8)
        call jeexin(nomcha(1:8)//'.CHME.LIGRE.LIEL', iret)
        if (iret .le. 0) cycle
        call exisd('CHAMP_GD', nomcha(1:8)//'.CHME.CMULT', iret)
        if (iret .le. 0) cycle
!
        call codent(nbres+1, 'D0', masque(12:14))
!
        ligrch = nomcha//'.CHME.LIGRE'
        lpain(1) = 'PDDLMUR'
        lchin(1) = nomcha//'.CHME.CMULT'
        lpain(2) = 'PDDLIMR'
        lchin(2) = disp_didi(1:19)
        lpain(3) = 'PALPHAR'
        lchin(3) = chalph(1:19)
        lpaout(1) = 'PVECTUR'
        lchout(1) = masque(1:19)
        call calcul('S', option, ligrch, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, &
                    'OUI')
        nbres = nbres+1
        call reajre(vect_elem, lchout(1), 'V')
    end do
!
! - Assembly
!
    call assvec('V', vect_asse, 1, vect_elem, [1.d0], nume_dof)

99  continue
!
    call jedema()
end subroutine
