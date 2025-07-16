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
subroutine te0245(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/getDensity.h"
#include "asterfort/teattr.h"
!
    character(len=*) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!           CALCUL DES TERMES PROPRES A UNE STRUCTURE  (ELEMENT DE BARRE)
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       OPTION  : 'MASS_INER      : CALCUL DES CARACTERISTIQUES DE STRUCTURES
!       NOMTE   :
!        'MECA_BARRE'       : BARRE
!        'MECA_2D_BARRE'    : BARRE 2D
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: lcastr, lmater, lsect, igeom
    integer(kind=8) :: dimModel, ier
    real(kind=8) :: rho, a, xl
    character(len=8) :: attrib
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PMATERC', 'L', lmater)

    call teattr('S', 'DIM_COOR_MODELI', attrib, ier)
    read (attrib, '(I8)') dimModel
!
    call getDensity(zi(lmater), rho)
!
!   recuperation des caracteristiques generales des sections
    call jevech('PCAGNBA', 'L', lsect)
    a = zr(lsect)
!
!   Longueur de l'élément
    if (dimModel .eq. 3) then
        xl = lonele(igeom=igeom)
    else if (dimModel .eq. 2) then
        xl = lonele(dime=2, igeom=igeom)
    else
        ASSERT(ASTER_FALSE)
    end if
!
!   calcul des caracteristiques elementaires
    if (option .eq. 'MASS_INER') then
        call jevech('PMASSINE', 'E', lcastr)
!       masse et cdg de l'element
        if (nomte .eq. 'MECA_BARRE') then
            zr(lcastr) = rho*a*xl
            zr(lcastr+1) = (zr(igeom+4)+zr(igeom+1))/2.d0
            zr(lcastr+2) = (zr(igeom+5)+zr(igeom+2))/2.d0
            zr(lcastr+3) = (zr(igeom+6)+zr(igeom+3))/2.d0
        else if (nomte .eq. 'MECA_2D_BARRE') then
            zr(lcastr) = rho*a*xl
            zr(lcastr+1) = (zr(igeom+3)+zr(igeom+1))/2.d0
            zr(lcastr+2) = (zr(igeom+4)+zr(igeom+2))/2.d0
        end if
!       inertie de l'element
        zr(lcastr+4) = 0.d0
        zr(lcastr+5) = 0.d0
        zr(lcastr+6) = 0.d0
        zr(lcastr+7) = 0.d0
        zr(lcastr+8) = 0.d0
        zr(lcastr+9) = 0.d0
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
