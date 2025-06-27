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

subroutine d1mamc(fami, mater, instan, poum, kpg, &
                  ksp, angl, nbsig, d1)
!.======================================================================
    implicit none
!
!      D1MAMC :   CALCUL DE L'INVERSE DE LA MATRICE DE HOOKE
!                 POUR LES ELEMENTS ISOPARAMETRIQUES POUR DES
!                 MATERIAUX ISOTROPE, ORTHOTROPE ET ISOTROPE TRANSVERSE
!
!   ARGUMENT        E/S  TYPE         ROLE
!    MATER          IN     I        MATERIAU
!    INSTAN         IN     R        INSTANT DE CALCUL (0 PAR DEFAUT)
!    ANGL(3)        IN     R        ANGLES NAUTIQUE DEFINISSANT LE REPERE
!                                   D'ORTHOTROPIE
!    NBSIG          IN     I        NOMBRE DE CONTRAINTES ASSOCIE A
!                                   L'ELEMENT
!    D1(NBSIG,1)    OUT    R        MATRICE DE HOOKE
!
!
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterfort/d1ma3d.h"
#include "asterfort/d1macp.h"
#include "asterfort/d1madp.h"
#include "asterfort/lteatt.h"
#include "asterfort/utmess.h"
    character(len=*) :: fami, poum
    integer(kind=8) :: kpg, ksp
    integer(kind=8) :: mater, nbsig
    real(kind=8) :: angl(3), d1(nbsig, 1), instan
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
!       ------------------------
! ----  CAS MASSIF 3D ET FOURIER
!       ------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    if (lteatt('DIM_TOPO_MAILLE', '3') .or. lteatt('FOURIER', 'OUI')) then
!
        call d1ma3d(fami, mater, instan, poum, kpg, &
                    ksp, angl, d1)
!
!       ----------------------------------------
! ----  CAS DEFORMATIONS PLANES ET AXISYMETRIQUE
!       ----------------------------------------
    elseif (lteatt('D_PLAN', 'OUI') .or. lteatt('AXIS', 'OUI')) &
        then
!
        call d1madp(fami, mater, instan, poum, kpg, &
                    ksp, angl(1), d1)
!
!       ----------------------
! ----  CAS CONTRAINTES PLANES
!       ----------------------
    else if (lteatt('C_PLAN', 'OUI')) then
!
        call d1macp(fami, mater, instan, poum, kpg, &
                    ksp, angl(1), d1)
!
    else
        call utmess('F', 'ELEMENTS_11')
    end if
!.============================ FIN DE LA ROUTINE ======================
end subroutine
