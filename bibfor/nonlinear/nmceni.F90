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
subroutine nmceni(numedd, depdel, deppr1, deppr2, rho, &
                  eta, f)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=24) :: numedd
    character(len=19) ::  depdel, deppr1, deppr2
    real(kind=8) :: eta, rho, f
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE - SELECTION PARAMETRE)
!
! CALCUL DU PARAMETRE DE SELECTION DE TYPE NORM_INCR_DEPL
!
! ----------------------------------------------------------------------
!
!
! IN  NUMEDD : NUME_DDL
! IN  DEPDEL : INCREMENT DE DEPLACEMENT DEPUIS DEBUT PAS DE TEMPS
! IN  DEPPR1 : INCREMENT DE DEPLACEMENT K-1.F_DONNE
! IN  DEPPR2 : INCREMENT DE DEPLACEMENT K-1.F_PILO
! IN  RHO    : PARAMETRE DE RECHERCHE LINEAIRE
! IN  ETA    : PARAMETRE DE PILOTAGE
! OUT F      : VALEUR DU CRITERE
!
!
!
!
    character(len=19) :: profch
    integer(kind=8) :: neq, i
    integer(kind=8), pointer :: deeq(:) => null()
    real(kind=8), pointer :: depde(:) => null()
    real(kind=8), pointer :: du0(:) => null()
    real(kind=8), pointer :: du1(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    f = 0.d0
!
! --- INFORMATIONS SUR NUMEROTATION
!
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
    call dismoi('NUME_EQUA', depdel, 'CHAM_NO', repk=profch)
    call jeveuo(profch(1:19)//'.DEEQ', 'L', vi=deeq)
!
! --- ACCES AUX VECTEURS SOLUTIONS
!
    call jeveuo(depdel(1:19)//'.VALE', 'L', vr=depde)
    call jeveuo(deppr1(1:19)//'.VALE', 'L', vr=du0)
    call jeveuo(deppr2(1:19)//'.VALE', 'L', vr=du1)

! --- CALCUL DE LA NORME
    do i = 1, neq
        if (deeq(2*i+2) .gt. 0) then
            f = f+(depde(1+i)+rho*du0(1+i)+eta*du1(1+i))**2
        end if
    end do

    call jedema()
end subroutine
