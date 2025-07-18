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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmreso(fonact, cndonn, cnpilo, cncine, solveu, &
                  maprec, matass, depso1, depso2, rescvg)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/isfonc.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmdebg.h"
#include "asterfort/resoud.h"
#include "asterfort/vtzero.h"
!
    integer(kind=8) :: fonact(*)
    character(len=19) :: maprec, matass
    character(len=19) :: solveu, cndonn, cnpilo
    character(len=19) :: cncine, depso1, depso2
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - CALCUL)
!
! RESOLUTION AVEC PILOTAGE  K.U = F0 + ETA.F1
! SUR DDLS PHYSIQUES
!
! ----------------------------------------------------------------------
!
!
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IN  CNDONN : SECOND MEMBRE DONNE
! IN  CNPILO : SECOND MEMBRE PILOTE
! IN  CNCINE : CHAM_NO DE CHARGE CINEMATIQUE
! IN  SOLVEU : SOLVEUR
! IN  MAPREC : MATRICE DE PRECONDITIONNEMENT (GCPC)
! IN  MATASS : MATRICE ASSEMBLEE
! OUT DEPSO1 : SOLUTION DU DU SYSTEME K.U = F EN L'ABSENCE DE PILOTAGE
! OUT DEPSO2 : SOLUTION DU DU SYSTEME K.U = F AVEC PILOTAGE
! OUT RESCVG : CODE RETOUR DE LA RESOLUTION
!
!
!
!
    aster_logical :: lpilo
    integer(kind=8) :: rescvg, ifm, niv
    complex(kind=8) :: c16bid
    character(len=19) :: crgc
    c16bid = dcmplx(0.d0, 0.d0)
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
!
! --- INITIALISATIONS
!
    crgc = '&&RESGRA_GCPC'
    call vtzero(depso1)
    call vtzero(depso2)
!
! --- FONCTIONNALITES ACTIVEES
!
    lpilo = isfonc(fonact, 'PILOTAGE')
!
! --- SECOND MEMBRE FIXE
!
    if (niv .ge. 2) then
        call nmdebg('VECT', cndonn, 6)
    end if
!
! --- SECOND MEMBRE PILOTE
!
    if (lpilo) then
        if (niv .ge. 2) then
            call nmdebg('VECT', cnpilo, 6)
        end if
    end if
!
    if (niv .ge. 2) then
        call nmdebg('MATA', matass, 6)
    end if
!
! --- INVERSION DE LA PARTIE FIXE
!
    call resoud(matass, maprec, solveu, cncine, 0, &
                cndonn, depso1, 'V', [0.d0], [c16bid], &
                crgc, .true._1, -9999, rescvg)
!
! --- ERREUR SANS POSSIBILITE DE CONTINUER
!
    if (rescvg .eq. 1) goto 999
!
! --- INVERSION DE LA PARTIE PILOTEE
!
    if (lpilo) then
        call resoud(matass, maprec, solveu, cncine, 0, &
                    cnpilo, depso2, 'V', [0.d0], [c16bid], &
                    crgc, .true._1, -9999, rescvg)
        if (rescvg .eq. 1) goto 999
    end if
!
! - Print solution
!
    if (niv .ge. 2) then
        call nmdebg('VECT', depso1, 6)
        if (lpilo) then
            call nmdebg('VECT', depso2, 6)
        end if
    end if
!
999 continue
!
    call jedetr(crgc//'.CRTI')
    call jedetr(crgc//'.CRTR')
    call jedetr(crgc//'.CRDE')
!
    call jedema()
end subroutine
