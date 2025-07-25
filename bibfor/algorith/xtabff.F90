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
subroutine xtabff(nbfond, nfon, ndim, fiss, operation)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ltcrsd.h"
#include "asterfort/ltnotb.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
    integer(kind=8) :: ndim, nbfond, nfon
    character(len=8) :: fiss
    character(len=16) :: operation
!
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (PREPARATION)
!
!     CONSTRUCTION DES 2 TABLES SUR LES FONDS DE FISSURE :
!             - TABLE DES COORDONNEES DES FONDS DE FISSURE
!             - TABLE DU NOMBRE DE FONDS DE FISSURE
!
! ----------------------------------------------------------------------
!
!
! I/O FISS   : NOM DE LA FISSURE
! IN  NBFOND : NOMBRE DE FONDS DE FISSURES DETECTES
! IN  NFON   : NOMBRE DE POINTS DE FOND DE FISSURE
! IN  NDIM   : DIMENSION DE L'ESPACE
!
!
!
!
    integer(kind=8) :: ifm, niv, i
    integer(kind=8) :: npara, nfonl, nfondl, vali(2)
    real(kind=8) :: vale(4), r8bid
    character(len=1) :: typar2(3), typar3(6)
    character(len=8) :: k8bid
    character(len=12) :: nopar2(3), nopar3(6)
    character(len=19) :: tabcoo, tabnb
    complex(kind=8) :: c16b
    integer(kind=8), pointer :: fondmult(:) => null()
    real(kind=8), pointer :: fondfiss(:) => null()
    data nopar2/'NUME_FOND', 'COOR_X', 'COOR_Y'/
    data typar2/'I', 'R', 'R'/
    data nopar3/'NUME_FOND', 'NUM_PT', 'ABSC_CURV',&
     &                                       'COOR_X', 'COOR_Y', 'COOR_Z'/
    data typar3/'I', 'I', 'R', 'R', 'R', 'R'/
!
! ----------------------------------------------------------------------
!
    call jemarq()
    c16b = (0.d0, 0.d0)
    r8bid = 0.d0
!
    call infdbg('XFEM', ifm, niv)
!
!     S'IL N'Y A PAS DE FOND DE FISSURE ON SORT
    if (nbfond .eq. 0 .and. operation .ne. 'PROPA_COHESIF') goto 999
!
    call jeveuo(fiss//'.FONDMULT', 'L', vi=fondmult)
    call jeveuo(fiss//'.FONDFISS', 'L', vr=fondfiss)
!
    if (operation .eq. 'PROPA_COHESIF') then
        call jelira(fiss//'.FONDMULT', 'LONUTI', nbfond, k8bid)
        nbfond = nbfond/2
        call jelira(fiss//'.FONDFISS', 'LONUTI', nfon, k8bid)
        nfon = nfon/6
    end if
    call ltcrsd(fiss, 'G')
!
!     ------------------------------------------------------------------
!     CONSTRUCTION DE LA TABLE DES COORDONNEES DES FONDS DE FISSURE
!     ------------------------------------------------------------------
!
    call ltnotb(fiss, 'FOND_FISS', tabcoo)
    call tbcrsd(tabcoo, 'G')
!
    if (ndim .eq. 2) then
        npara = 3
        call tbajpa(tabcoo, npara, nopar2, typar2)
        do i = 1, nbfond
            vali(1) = i
            vale(1) = fondfiss(4*(i-1)+1)
            vale(2) = fondfiss(4*(i-1)+2)
            call tbajli(tabcoo, npara, nopar2, vali, vale, &
                        [c16b], k8bid, 0)
        end do
    else if (ndim .eq. 3) then
        npara = 6
        call tbajpa(tabcoo, npara, nopar3, typar3)
        nfonl = 1
        nfondl = 0
        do i = 1, nfon
            if (fondmult(2*nfondl+1) .eq. i) then
                nfondl = nfondl+1
                nfonl = 1
            else
                nfonl = nfonl+1
            end if
            vali(1) = nfondl
            vali(2) = nfonl
            vale(1) = fondfiss(4*(i-1)+4)
            vale(2) = fondfiss(4*(i-1)+1)
            vale(3) = fondfiss(4*(i-1)+2)
            vale(4) = fondfiss(4*(i-1)+3)
            call tbajli(tabcoo, npara, nopar3, vali, vale, &
                        [c16b], k8bid, 0)
        end do
    else
        ASSERT(.false.)
    end if
!
!     ------------------------------------------------------------------
!     CONSTRUCTION DE LA TABLE DU NOMBRE DE FONDS DE FISSURE
!     ------------------------------------------------------------------
!
    call ltnotb(fiss, 'NB_FOND_FISS', tabnb)
    call tbcrsd(tabnb, 'G')
    call tbajpa(tabnb, 1, 'NOMBRE', 'I')
    call tbajli(tabnb, 1, 'NOMBRE', [nbfond], [r8bid], &
                [c16b], k8bid, 0)
!
999 continue
    call jedema()
end subroutine
