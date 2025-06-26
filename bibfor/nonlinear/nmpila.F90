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
subroutine nmpila(numedd, sdpilo, dtau, depdel, &
                  ddepl0, ddepl1, nbeffe, eta, pilcvg)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/zerop2.h"
!
    integer(kind=8) :: pilcvg, nbeffe
    character(len=19) :: sdpilo
    character(len=24) :: numedd
    character(len=19) :: ddepl0, ddepl1, depdel
    real(kind=8) :: dtau, eta(2)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE - CALCUL DE ETA)
!
! RESOLUTION DE L'EQUATION DE PILOTAGE PAR LONGUEUR D'ARC
!
! ----------------------------------------------------------------------
!
!
! IN  NUMEDD : NUME_DDL
! IN  SDPILO : SD PILOTAGE
! IN  DEPDEL : INCREMENT DE DEPLACEMENT DEPUIS DEBUT PAS DE TEMPS
! IN  DDEPL0 : INCREMENT DE DEPLACEMENT K-1.F_DONNE
! IN  DDEPL1 : INCREMENT DE DEPLACEMENT K-1.F_PILO
! IN  DTAU   : SECOND MEMBRE DE L'EQUATION DE PILOTAGE
! OUT NBEFFE : NOMBRE DE SOLUTIONS EFFECTIVES
! OUT ETA    : ETA_PILOTAGE
! OUT PILCVG : CODE DE CONVERGENCE POUR LE PILOTAGE
!                -1 : PAS DE CALCUL DU PILOTAGE
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : PAS DE SOLUTION
!                 2 : BORNE ATTEINTE -> FIN DU CALCUL
!
!
!
!
    integer(kind=8) :: i, nrac
    real(kind=8) :: r0, r1, r2, dtau2, rac(2)
    integer(kind=8) :: jdepde
    integer(kind=8) :: neq
    integer(kind=8) :: ifm, niv
    character(len=19) :: chapil
    real(kind=8), pointer :: coef(:) => null()
    real(kind=8), pointer :: dep0(:) => null()
    real(kind=8), pointer :: dep1(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('PILOTAGE', ifm, niv)
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        write (ifm, *) '<PILOTAGE> ...... PILOTAGE PAR LONGUEUR D''ARC'
    end if
!
! --- INITIALISATIONS
!
    pilcvg = -1
    dtau2 = dtau**2
    r0 = -dtau2
    r1 = 0.d0
    r2 = 0.d0
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
!
! --- ACCES OBJETS JEVEUX
!
    call jeveuo(ddepl0(1:19)//'.VALE', 'L', vr=dep0)
    call jeveuo(ddepl1(1:19)//'.VALE', 'L', vr=dep1)
    call jeveuo(depdel(1:19)//'.VALE', 'L', jdepde)
    chapil = sdpilo(1:14)//'.PLCR'
    call jeveuo(chapil(1:19)//'.VALE', 'L', vr=coef)

    do i = 1, neq
        r0 = r0+coef(i)*(zr(jdepde+i-1)+dep0(i))**2
        r1 = r1+coef(i)*(zr(jdepde+i-1)+dep0(i))*dep1(i)
        r2 = r2+coef(i)*dep1(i)*dep1(i)
    end do

!
    r1 = 2.d0*r1
    if (r2 .eq. 0) then
        ASSERT(.false.)
    end if
    if (niv .ge. 2) then
        write (ifm, *) '<PILOTAGE> ....EQUATION X2+BX+C: ', r1/r2, r0/r2
    end if
!
! --- RESOLUTION DE L'EQUATION DE DEGRE DEUX
!
    call zerop2(r1/r2, r0/r2, rac, nrac)
!
    if (nrac .eq. 0) then
        pilcvg = 1
    else if (nrac .eq. 1) then
        pilcvg = 0
        nbeffe = 1
        eta(1) = rac(1)
    else
        pilcvg = 0
        nbeffe = 2
        eta(1) = rac(1)
        eta(2) = rac(2)
    end if
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        write (ifm, *) '<PILOTAGE> ...... SOLUTIONS: ', nrac, rac
    end if
!
    call jedema()
!
end subroutine
