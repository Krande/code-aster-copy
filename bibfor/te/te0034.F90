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
subroutine te0034(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_ME_FR1D3D  '
!                                   'CHAR_ME_FF1D3D  '
!                        ELEMENT  : 'MEBODKT'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    real(kind=8) :: dx, dy, dz, long, fx, fy, fz, mx, my, mz
!
    real(kind=8) :: valpar(4)
    character(len=8) :: nompar(4)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icod1, icod2, icod3, icod4, icod5, icod6, iforc
    integer(kind=8) :: igeom, ino, itpsr, ivectu
!-----------------------------------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PVECTUR', 'E', ivectu)
!
!     -- CALCUL DE LA LONGUEUR DU SEGMENT :
    dx = zr(igeom-1+4)-zr(igeom-1+1)
    dy = zr(igeom-1+5)-zr(igeom-1+2)
    dz = zr(igeom-1+6)-zr(igeom-1+3)
    long = sqrt(dx**2+dy**2+dz**2)
!
!     -- CALCUL DE LA FORCE MOYENNE :
    if (option(11:16) .eq. 'FR1D3D') then
        call jevech('PFR1D3D', 'L', iforc)
        fx = zr(iforc-1+1)
        fy = zr(iforc-1+2)
        fz = zr(iforc-1+3)
        mx = zr(iforc-1+4)
        my = zr(iforc-1+5)
        mz = zr(iforc-1+6)
    else if (option(11:16) .eq. 'FF1D3D') then
        call jevech('PFF1D3D', 'L', iforc)
        call jevech('PINSTR', 'L', itpsr)
!
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
!
        valpar(1) = (zr(igeom-1+1)+zr(igeom-1+4))/2.0d0
        valpar(2) = (zr(igeom-1+2)+zr(igeom-1+5))/2.0d0
        valpar(3) = (zr(igeom-1+3)+zr(igeom-1+6))/2.0d0
        valpar(4) = zr(itpsr)
!
        call fointe('FM', zk8(iforc-1+1), 4, nompar, valpar, &
                    fx, icod1)
        call fointe('FM', zk8(iforc-1+2), 4, nompar, valpar, &
                    fy, icod2)
        call fointe('FM', zk8(iforc-1+3), 4, nompar, valpar, &
                    fz, icod3)
        call fointe('FM', zk8(iforc-1+4), 4, nompar, valpar, &
                    mx, icod4)
        call fointe('FM', zk8(iforc-1+5), 4, nompar, valpar, &
                    my, icod5)
        call fointe('FM', zk8(iforc-1+6), 4, nompar, valpar, &
                    mz, icod6)
    else
        call utmess('F', 'ELEMENTS2_77', sk=option)
    end if
!
!     -- AFFECTATION DU RESULTAT:
!
    do 1, ino = 1, 2
        zr(ivectu-1+(ino-1)*6+1) = fx*long/2.0d0
        zr(ivectu-1+(ino-1)*6+2) = fy*long/2.0d0
        zr(ivectu-1+(ino-1)*6+3) = fz*long/2.0d0
        zr(ivectu-1+(ino-1)*6+4) = mx*long/2.0d0
        zr(ivectu-1+(ino-1)*6+5) = my*long/2.0d0
        zr(ivectu-1+(ino-1)*6+6) = mz*long/2.0d0
1   end do
!
end subroutine
