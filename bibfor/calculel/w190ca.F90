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

subroutine w190ca(modele, carele, chmar1, chefge, chamfer, chefge0, chmar2)
    implicit none
#include "asterfort/calcul.h"
#include "asterfort/exlim3.h"
#include "asterfort/mecara.h"
    character(len=8) :: modele, carele
    character(len=19) :: chmar1, chmar2, chefge, chamfer, chefge0
!
! ----------------------------------------------------------------------
!     CALCUL DE L'OPTION VERIFICATION DE FERRAILLAGE
!
! IN  MODELE  : NOM DU MODELE
! IN  CARELE  : CARACTERISTIQUES COQUES
! IN  CHMAR1  : CHAMP DE MAR1_R
! IN  CHEFGE  : CHAMP DE EFGE_ELNO
! IN  CHEFGE0 : CHAMP DE EFGE_ELNO DE REFERENCE
! IN  CHAMFER : CHAMP DE FERRAILLAGE
! OUT CHMAR2  : RESULTAT DU CALCUL DE LA MARGE MECANIQUE
!
    character(len=8) :: lpain(6), lpaout(1)
    character(len=16) :: option
    character(len=19) :: chcara(18)
    character(len=19) :: lchin(6), lchout(2), ligrel
!
    call exlim3('AFFE', 'G', modele, ligrel)
    option = 'MARG_ELEM'
!
    call mecara(carele, chcara)
!
    lpain(1) = 'PCACOQU'
    lchin(1) = chcara(7)
    lpain(2) = 'PVFER1'
    lchin(2) = chmar1
    lpain(3) = 'PEFFORR'
    lchin(3) = chefge
    lpain(4) = 'PCAGEPO'
    lchin(4) = chcara(5)
    lpain(5) = 'PEFFOR0'
    lchin(5) = chefge0
    lpain(6) = 'PVFER0'
    lchin(6) = chamfer
!
!
    lpaout(1) = 'PVFER2'
    lchout(1) = chmar2
!
    call calcul('S', option, ligrel, 6, lchin, &
                lpain, 1, lchout, lpaout, 'G', &
                'OUI')
!
end subroutine
