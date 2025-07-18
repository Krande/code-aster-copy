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
subroutine cfgcac(resoco, tole, neq, nbliai, nbliac)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/calatm.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/r8inir.h"
    character(len=24) :: resoco
    integer(kind=8) :: neq, nbliai, nbliac
    real(kind=8) :: tole
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - GCP)
!
! ACTIVATION DES LIAISONS ET CALCUL DE LA FORCE DE CONTACT
!
! ----------------------------------------------------------------------
!
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT
! IN  NEQ    : NOMBRE D'EQUATIONS
! IN  TOLE   : TOLERANCE POUR DETECTER PRESSION NULLE
! OUT NBLIAC : NOMBRE DE LIAISONS ACTIVES
!
!
!
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iliai, jdecal, nbddl
    character(len=24) :: apcoef, apddl, appoin
    integer(kind=8) :: japcoe, japddl, japptr
    character(len=19) :: mu, atmu, liac
    integer(kind=8) :: jmu, jatmu, jliac
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    appoin = resoco(1:14)//'.APPOIN'
    apcoef = resoco(1:14)//'.APCOEF'
    apddl = resoco(1:14)//'.APDDL'
    mu = resoco(1:14)//'.MU'
    atmu = resoco(1:14)//'.ATMU'
    liac = resoco(1:14)//'.LIAC'
!
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apcoef, 'L', japcoe)
    call jeveuo(apddl, 'L', japddl)
    call jeveuo(mu, 'L', jmu)
    call jeveuo(atmu, 'E', jatmu)
    call jeveuo(liac, 'E', jliac)
!
! --- CALCUL DE ATMU SUR TOUTES LES LIAISONS
!
    call r8inir(neq, 0.d0, zr(jatmu), 1)
    do iliai = 1, nbliai
        jdecal = zi(japptr+iliai-1)
        nbddl = zi(japptr+iliai)-zi(japptr+iliai-1)
        call calatm(neq, nbddl, zr(jmu+iliai-1), zr(japcoe+jdecal), zi(japddl+jdecal), &
                    zr(jatmu))
    end do
!
! --- COMPTE DES LIAISONS ACTIVES ET ACTIVATION
!
    nbliac = 0
    do iliai = 1, nbliai
        if (zr(jmu+iliai-1) .gt. tole) then
            nbliac = nbliac+1
            zi(jliac+nbliac-1) = iliai
        end if
    end do
!
    call jedema()
!
end subroutine
