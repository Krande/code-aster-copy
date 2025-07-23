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
subroutine vethbu(model, matasz, loadNameJv, loadInfoJv, &
                  chtni, vebtem)
    implicit none
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/conlag.h"
#include "asterfort/corich.h"
#include "asterfort/gcncon.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/vemare.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: loadNameJv, loadInfoJv
    character(len=24) :: chtni, vebtem
    character(len=*) :: matasz
! ----------------------------------------------------------------------
! CALCUL DES VECTEURS ELEMENTAIRES B * TEMPERATURE
!
! IN  MODELE  : NOM DU MODELE
! IN  MATASZ  : SD MATRICE ASSEMBLEE
! IN  CHARGE  : LISTE DES CHARGES
! IN  INFCHA  : INFORMATIONS SUR LES CHARGEMENTS
! IN  CARELE  : CHAMP DE CARA_ELEM
! IN  MATE    : MATERIAU CODE
! IN  CHTNI   : IEME ITEREE DU CHAMP DE TEMPERATURE
! OUT VEBTEM  : VECTEURS ELEMENTAIRES PROVENANT DE B * TMOINS
!
!
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=1, nbin=3)
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    real(kind=8) :: alpha
    character(len=8) :: nomcha, newnom
    character(len=16) :: option
    character(len=19) :: vecel, matass
    character(len=24) :: ligrch, chalph
    integer(kind=8) :: iret, nchar, icha, jchar, jinf, jdir, ndir
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- ACCES AUX CHARGES
!
    call jeexin(loadNameJv, iret)
    if (iret .ne. 0) then
        call jelira(loadNameJv, 'LONMAX', nchar)
        call jeveuo(loadNameJv, 'L', jchar)
    else
        nchar = 0
    end if
!
! --- ALLOCATION DU VECT_ELEM RESULTAT :
!
    call jeexin(vebtem, iret)
    vecel = '&&VEBUEE           '
    if (iret .eq. 0) then
        vebtem = vecel//'.RELR'
        call vemare('V', vecel, model(1:8))
        call wkvect(vebtem, 'V V K24', nchar, jdir)
    else
        call jeveuo(vebtem, 'E', jdir)
    end if
    call jeecra(vebtem, 'LONUTI', 0)
!
! --- ALLOCATION DE LA CARTE DU CONDITIONNEMENT DES LAGRANGES
!
    matass = matasz
    call conlag(matass, alpha)
    chalph = '&&VETHBU.CH_NEUT_R'
    call mecact('V', chalph, 'MODELE', model, 'NEUT_R  ', &
                ncmp=1, nomcmp='X1', sr=alpha)
!
! --- CALCUL DE L'OPTION B.T
!
    option = 'THER_BU_R'
    if (nchar .gt. 0) then
        call jeveuo(loadInfoJv, 'L', jinf)
        do icha = 1, nchar
            ndir = 0
            if (zi(jinf+icha) .gt. 0) then
                nomcha = zk24(jchar+icha-1) (1:8)
                ligrch = nomcha//'.CHTH.LIGRE'
                lpain(1) = 'PDDLMUR'
                lchin(1) = nomcha//'.CHTH.CMULT'
                lpain(2) = 'PDDLIMR'
                lchin(2) = chtni(1:19)
                lpain(3) = 'PALPHAR'
                lchin(3) = chalph(1:19)
                lpaout(1) = 'PVECTTR'
                lchout(1) = '&&VETHBU.???????'
                call gcncon('.', newnom)
                lchout(1) (10:16) = newnom(2:8)
                call corich('E', lchout(1), ichin_=-1)
                call calcul('S', option, ligrch, nbin, lchin, &
                            lpain, nbout, lchout, lpaout, 'V', &
                            'OUI')
                ndir = ndir+1
                zk24(jdir+ndir-1) = lchout(1)
            end if
            call jeecra(vebtem, 'LONUTI', ndir)
!
        end do
!-
    end if
! FIN ------------------------------------------------------------------
    call jedema()
end subroutine
