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

subroutine btldth(fami, xi3, nb1, kpg, btild, &
                  wgt, indic, young, nu, alpha, &
                  temper, forthi)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8nnem.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nb1, kpg
    real(kind=8) :: wgt, young, nu, alpha, xi3
    real(kind=8) :: btild(5, 42), forthi(1), vecthr(2)
    integer(kind=8) :: jcou, imoy
    character(len=4) :: fami
    real(kind=8) :: p1xi3, p2xi3, p3xi3
!
!
!     CALCUL DE TEMPERATURE AUX PTS D'INTEGRATION
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, indic, iret1, iret2, iret3, iret4, k
!
    real(kind=8) :: temper, tinf, tmoy, tref, tsup
!-----------------------------------------------------------------------
    p1xi3 = 1-xi3*xi3
    p2xi3 = -xi3*(1-xi3)/2.d0
    p3xi3 = xi3*(1+xi3)/2.d0
    call jevech('PNBSP_I', 'L', jcou)
    imoy = (3*zi(jcou)+1)/2
    call rcvarc(' ', 'TEMP', 'REF', fami, 1, &
                1, tref, iret1)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                1, tinf, iret2)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                imoy, tmoy, iret3)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                3*zi(jcou), tsup, iret4)
    if ((iret2+iret3+iret4) .eq. 0) then
        if ((iret1 .eq. 1) .or. (indic .eq. 0)) then
            call utmess('F', 'COMPOR5_43')
        else
            temper = (tmoy*p1xi3+tinf*p2xi3+tsup*p3xi3)-tref
        end if
    else
        temper = r8nnem()
    end if
!
!
    if (indic .eq. 1) then
!
        vecthr(1) = young*alpha*temper/(1.d0-nu)
        vecthr(2) = vecthr(1)
!
!     CONSTRUCTION DES EFFORTS DUS AUX DILATATIONS THERMIQUES
!
        do i = 1, 5*nb1+2
            forthi(i) = 0.d0
            do k = 1, 2
                forthi(i) = forthi(i)+btild(k, i)*vecthr(k)*wgt
            end do
        end do
    end if
!
!
end subroutine
