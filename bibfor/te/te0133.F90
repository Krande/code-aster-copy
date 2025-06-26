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
subroutine te0133(option, nomte)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/getElemOrientation.h"
!
    character(len=16) :: option, nomte
!----------------------------------------------------------------------
!
!   - FONCTION REALISEE: CALCUL DU REPERE LOCAL DONNE PAR L'UTILISATEUR
!   - TYPES D'ELEMENT 3D ET C_PLAN, D_PLAN, AXIS (2D)
!
!
!----------------------------------------------------------------------
!
    integer(kind=8) :: jgeom, jrepl1, jrepl2, jrepl3
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano
    integer(kind=8) :: i
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: ux(3), uy(3), uz(3)
    real(kind=8) :: ang(3)
    aster_logical :: interf
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    ASSERT(nnos .le. 8)
!
    interf = lteatt('INTERFACE', 'OUI')
    if (interf) ndim = ndim+1
!
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PREPLO1', 'E', jrepl1)
    call jevech('PREPLO2', 'E', jrepl2)
    if (ndim .eq. 3) call jevech('PREPLO3', 'E', jrepl3)
!
!     CALCUL DES ANGLES NAUTIQUES
    call getElemOrientation(ndim, nnos, jgeom, ang)
!     CALCUL DE LA MATRICE DE PASSAGE GLOBAL->LOCAL
    call matrot(ang, pgl)
!
!     (UX,UY,UZ) : VECTEUR LOCAUX UTILISATEUR DANS LE BASE GLOBALE
!     UX EST LA PREMIERE LIGNE DE PGL
!     UY LA DEUXIEME
!     UZ LA TROISIEME
!
    ux(1) = pgl(1, 1)
    ux(2) = pgl(1, 2)
    ux(3) = pgl(1, 3)
!
    uy(1) = pgl(2, 1)
    uy(2) = pgl(2, 2)
    uy(3) = pgl(2, 3)
!
    uz(1) = pgl(3, 1)
    uz(2) = pgl(3, 2)
    uz(3) = pgl(3, 3)
!
    do i = 1, 2
        zr(jrepl1-1+i) = ux(i)
        zr(jrepl2-1+i) = uy(i)
    end do
!
    if (ndim .eq. 3) then
        zr(jrepl1-1+3) = ux(3)
        zr(jrepl2-1+3) = uy(3)
        do i = 1, 3
            zr(jrepl3-1+i) = uz(i)
        end do
    end if
!
end subroutine
