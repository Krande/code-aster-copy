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
! aslint: disable=W0413
!
subroutine nmceai(numedd, depdel, deppr1, deppr2, depold, &
                  sdpilo, rho, eta, f, &
                  indic)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: indic
    character(len=24) :: numedd
    character(len=19) :: sdpilo, depdel, depold, deppr1, deppr2
    real(kind=8) :: eta, rho, f
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE - SELECTION PARAMETRE)
!
! CALCUL DU PARAMETRE DE SELECTION DE TYPE ANGL_INCR_DEPL
!
! ----------------------------------------------------------------------
!
!
! IN  NUMEDD : NUME_DDL
! IN  SDPILO : SD PILOTAGE
! IN  DEPDEL : INCREMENT DE DEPLACEMENT DEPUIS DEBUT PAS DE TEMPS
! IN  DEPOLD : INCREMENT DE DEPLACEMENT PAS DE TEMPS PRECEDENT
! IN  DEPPR1 : INCREMENT DE DEPLACEMENT K-1.F_DONNE
! IN  DEPPR2 : INCREMENT DE DEPLACEMENT K-1.F_PILO
! IN  RHO    : PARAMETRE DE RECHERCHE LINEAIRE
! IN  ETA    : PARAMETRE DE PILOTAGE
! OUT F      : VALEUR DU CRITERE
! OUT INDIC  : 0 CRITERE NON UTILISABLE
!              1 CRITERE UTILISABLE
!
    real(kind=8) :: sca, nodup, coef
    real(kind=8) :: norm_depold
    integer(kind=8) :: jdu1
    integer(kind=8) :: neq, i
    character(len=19) :: selpil
    real(kind=8), pointer :: depde(:) => null()
    real(kind=8), pointer :: depol(:) => null()
    real(kind=8), pointer :: du0(:) => null()
    real(kind=8), pointer :: plsl(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!

!
! --- ACCES VECTEUR DE SELCTION CMP DX/DY/DZ
!
    selpil = sdpilo(1:14)//'.PLSL'
    call jeveuo(selpil(1:19)//'.VALE', 'L', vr=plsl)

!
!
! --- INITIALISATIONS--------------------------------
!
    sca = 0.d0
    nodup = 0.d0
    norm_depold = 0.d0
    f = 0.d0
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
!
! --- ACCES AUX VECTEURS SOLUTIONS
!
    call jeveuo(depdel(1:19)//'.VALE', 'L', vr=depde)
    call jeveuo(deppr1(1:19)//'.VALE', 'L', vr=du0)
    call jeveuo(deppr2(1:19)//'.VALE', 'L', jdu1)
    call jeveuo(depold(1:19)//'.VALE', 'L', vr=depol)
!
! --- CALCUL DE L'ANGLE

    do i = 1, neq
        coef = plsl(i)
        sca = sca+(depol(i)*(depde(i)+rho*du0(i)+eta*zr(jdu1+i-1)))*coef
        nodup = nodup+(depde(i)+rho*du0(i)+eta*zr(jdu1+i-1))**2
        norm_depold = norm_depold+depol(i)**2
    end do
!
    if (nodup .eq. 0.d0 .or. norm_depold .eq. 0.d0) then
        indic = 0
        f = 0.d0
    else
        indic = 1
        f = sca/sqrt(nodup*norm_depold)
    end if
    f = -f
!
    call jedema()
end subroutine
