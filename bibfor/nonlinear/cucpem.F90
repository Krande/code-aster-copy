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
subroutine cucpem(deficu, resocu, nbliai)
!
!
    implicit none
#include "jeveux.h"
!
#include "asterfort/compute_ineq_conditions_matrix.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=24) :: deficu, resocu
    integer(kind=8) :: nbliai
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - PENALISATION)
!
! CALCUL DE LA MATRICE DE CONTACT PENALISEE ELEMENTAIRE [E_N*AT]
!
! ----------------------------------------------------------------------
!
!
! IN  DEFICU : SD DE DEFINITION (ISSUE D'AFFE_CHAR_MECA)
! IN  RESOCU : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT
!
!
!
!
    character(len=24) :: apcoef, appoin
    integer(kind=8) :: japcoe, japptr
    character(len=24) :: coefpe
    integer(kind=8) :: jcoef_pena
    character(len=24) :: jeux
    integer(kind=8) :: jjeux
    character(len=24) :: enat
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    appoin = deficu(1:16)//'.POINOE'
    apcoef = resocu(1:14)//'.APCOEF'
    jeux = resocu(1:14)//'.APJEU'
    coefpe = deficu(1:16)//'.COEFPE'
    enat = resocu(1:14)//'.ENAT'
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apcoef, 'L', japcoe)
    call jeveuo(jeux, 'L', jjeux)
    call jeveuo(coefpe, 'L', jcoef_pena)
!
! --- CALCUL DE LA MATRICE DE CONTACT PENALISEE
!
    call compute_ineq_conditions_matrix(enat, nbliai, japptr, japcoe, jjeux, &
                                        jcoef_pena-1, 1, 1)
!
    call jedema()
!
end subroutine
