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

subroutine nirela(irela, jp, gm, gp, am, ap, bp, boa, aa, bb, daa, dbb, dboa, d2boa, iret)
!
!
    implicit none
#include "asterfort/assert.h"
    integer(kind=8), intent(in) :: irela
    real(kind=8), intent(in) :: jp, gm, gp
    real(kind=8), intent(out) :: am, ap, bp, boa, aa, bb, daa, dbb, dboa, d2boa
    integer(kind=8), intent(out) :: iret
!-----------------------------------------------------------------------
!          CALCUL DES OPTIONS DE MECANIQUE NON LINEAIRE
!           GRANDES DEFORMATIONS QUASI-INCOMPRESSIBLES
!
!  INITIALISATION DES FONCTIONS POUR LA RELATION : B(J) = B(A(G))
!-----------------------------------------------------------------------
! IN  IRELA  IDENTIFIANT DE LA RELATION : 1) J = 1 + G
!                                         2) ln(J) = G
!                                         3) J = exp(G)
!                                         4) J^2 = 1 + G
!                                         5) ln(J) = ln(1 + G)
!
! IN  JP     CHANGEMENT DE VOLUME EN T+
! IN  GM     GONFLEMENT EN T-
! IN  GP     GONFLEMENT EN T+
! OUT AM     A(GM)
! OUT AP     A(GP)
! OUT BP     B(JP)
! OUT BOA    B O A EN T+
! OUT AA(G)  DA/DG / A(G)
! OUT BB(J)  J * DB/DJ
! OUT DAA    DAA/DG
! OUT DBB    DBB/DJ
! OUT DBOA   D(B O A)/DG
! OUT D2BOA  D2(B O A)/DG2
! OUT IRET   O (OK) / -1 (ERROR)
!
    iret = 0
    if (jp .le. 0.d0) then
        iret = -1
        go to 999
    elseif (irela == 2 .or. irela == 3) then
        if (gm .ge. 200 .or. gp .ge. 200) then
            iret = -1
            go to 999
        end if
    end if
!
    if (irela .eq. 1) then
!-----------------------------------------------------------------------
!    APPLICATION: A(G) = 1+G   ET   B(J) = J
!-----------------------------------------------------------------------
        am = 1.d0+gm
        ap = 1.d0+gp
        bp = jp
        boa = ap
        aa = 1.d0/ap
        bb = jp
        daa = -1.d0/(1.d0+gp)**2
        dbb = 1.d0
        dboa = 1.d0
        d2boa = 0.d0
    else if (irela .eq. 2) then
!-----------------------------------------------------------------------
!    APPLICATION: A(G) = exp(G)   ET   B(J) = ln(J)
!-----------------------------------------------------------------------
        am = exp(gm)
        ap = exp(gp)
        bp = log(jp)
        boa = gp
        aa = 1.d0
        bb = 1.d0
        daa = 0.d0
        dbb = 0.d0
        dboa = 1.d0
        d2boa = 0.d0
    else if (irela .eq. 3) then
!-----------------------------------------------------------------------
!    APPLICATION: A(G) = exp(G)   ET   B(J) = J
!-----------------------------------------------------------------------
        am = exp(gm)
        ap = exp(gp)
        bp = jp
        boa = ap
        aa = 1.d0
        bb = jp
        daa = 1.d0
        dbb = jp
        dboa = ap
        d2boa = ap
    else if (irela .eq. 4) then
!-----------------------------------------------------------------------
!    APPLICATION: A(G) = SQRT(1+G)   ET   B(J) = J^2
!-----------------------------------------------------------------------
        am = sqrt(1.d0+gm)
        ap = sqrt(1.d0+gp)
        bp = jp**2
        boa = 1.d0+gp
        aa = 0.5d0/ap
        bb = 2.d0*jp**2
        daa = -0.25d0/ap**(3.d0/2.d0)
        dbb = 4.d0*jp
        dboa = 1.d0
        d2boa = 0.d0
    else if (irela .eq. 5) then
!-----------------------------------------------------------------------
!    APPLICATION: A(G) = 1+G   ET   B(J) = ln(J)
!-----------------------------------------------------------------------
        am = 1.d0+gm
        ap = 1.d0+gp
        bp = log(jp)
        boa = log(1.d0+gp)
        aa = 1.d0/ap
        bb = 1.d0
        daa = -1.d0/(ap)**2
        dbb = 0.d0
        dboa = 1.d0/(ap)
        d2boa = -1.d0/(ap)**2
    else
        ASSERT(ASTER_FALSE)
    end if
!
999 continue
!
    if (iret == -1) then
        am = 0.d0
        ap = 0.d0
        bp = 0.d0
        boa = 0.d0
        aa = 0.d0
        bb = 0.d0
        daa = 0.d0
        dbb = 0.d0
        dboa = 0.d0
        d2boa = 0.d0
    end if
end subroutine
