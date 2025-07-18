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

subroutine ef0231(nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/ppgan2.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
    character(len=16) :: nomte
! ......................................................................
!     CALCUL DE EFGE_ELNO
!     ------------------------------------------------------------------
!
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
    character(len=8) :: elrefe
    character(len=16) :: nomres(3)
    integer(kind=8) :: icodre(3)
    real(kind=8) :: e, nu, tpg, tgmoy, tgsup, tginf, tref
    real(kind=8) :: x3, eps(5), c, h, epsthe, valres(3)
    real(kind=8) :: e11, e22, k11, k22, ep11, ep22
    real(kind=8) :: dfdx(3), effopg(24)
    real(kind=8) :: jac, r, cosa, sina, cour
    integer(kind=8) :: i, k, kp, igeom, imate, icaco, idepl
    integer(kind=8) :: nno, npg, idfdk, ivf, iret, iret2, iret1, iret3, iret4
    integer(kind=8) :: jcoopg, ip, correc, jdfd2
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ieffor, ipoids, jgano, ndim, nnos
!-----------------------------------------------------------------------
    call elref1(elrefe)
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, &
                     jdfd2=jdfd2, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCACOQU', 'L', icaco)
    call jevech('PDEPLAR', 'L', idepl)
    call jevech('PEFFORR', 'E', ieffor)
    call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, &
                1, tref, iret)
!
    h = zr(icaco)
!JMP  CORREC = CORRECTION DE METRIQUE = 0 (NON) OU 1 (OUI)
!JMP  CORREC = ZR(ICACO+2)
    correc = nint(zr(icaco+2))
!
    do i = 1, npg*6
        effopg(i) = 0.d0
    end do
!
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(igeom), dfdx, &
                    cour, jac, cosa, sina)
!
        do i = 1, 5
            eps(i) = 0.d0
        end do
        r = 0.d0
        do i = 1, nno
            eps(1) = eps(1)+dfdx(i)*zr(idepl+3*i-3)
            eps(2) = eps(2)+dfdx(i)*zr(idepl+3*i-2)
            eps(3) = eps(3)+dfdx(i)*zr(idepl+3*i-1)
            eps(4) = eps(4)+zr(ivf+k+i-1)*zr(idepl+3*i-3)
            eps(5) = eps(5)+zr(ivf+k+i-1)*zr(idepl+3*i-1)
            r = r+zr(ivf+k+i-1)*zr(igeom+2*i-2)
        end do
!
        e11 = eps(2)*cosa-eps(1)*sina
        k11 = eps(3)
        if (nomte .eq. 'MECXSE3') then
            e22 = eps(4)/r
            k22 = -eps(5)*sina/r
        else
            e22 = 0.d0
            k22 = 0.d0
        end if
!
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, &
                    1, tginf, iret1)
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, &
                    2, tgmoy, iret2)
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, &
                    3, tgsup, iret3)
        iret4 = iret1+iret2+iret3
        ASSERT(iret4 .eq. 0 .or. iret4 .eq. 3)
!
!
!---- UTILISATION DE 4 POINTS DE GAUSS DANS L'EPAISSEUR
!---- COMME POUR LA LONGUEUR
!
        do ip = 1, npg
            x3 = zr(jcoopg-1+ip)
            if (iret4 .eq. 0) then
                tpg = tgmoy*(1.d0-(x3)**2)+tgsup*x3*(1.d0+x3)/2.d0- &
                      tginf*x3*(1.d0-x3)/2.d0
            else
                tpg = r8nnem()
            end if
            x3 = x3*h/2.d0
            ep11 = (e11+x3*k11)/(1.d0+(correc*x3*cour))
            nomres(1) = 'E'
            nomres(2) = 'NU'
            nomres(3) = 'ALPHA'
            call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                        ' ', 'ELAS', 1, 'TEMP', [tpg], &
                        2, nomres, valres, icodre, 1)
            call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                        ' ', 'ELAS', 1, 'TEMP', [tpg], &
                        1, nomres(3), valres(3), icodre(3), 0)
            e = valres(1)
            nu = valres(2)
            if (iret4 .eq. 0) then
                if ((icodre(3) .ne. 0) .or. (iret .eq. 1)) then
                    call utmess('F', 'CALCULEL_15')
                else
                    epsthe = (tpg-tref)*valres(3)*e/(1.d0-nu)
                end if
            else
                epsthe = 0.d0
            end if
!
            c = e/(1.d0-nu*nu)
            if (nomte .eq. 'MECXSE3') then
                ep22 = (e22+x3*k22)/(1.d0+(correc*cosa*x3/r))
                effopg(6*(kp-1)+1) = effopg(6*(kp-1)+1)+zr(ipoids-1+ip) &
                                     *(h/2.d0)*(c*(ep11+nu*ep22)-epsthe)
                effopg(6*(kp-1)+2) = effopg(6*(kp-1)+2)+zr(ipoids-1+ip) &
                                     *(h/2.d0)*(c*(nu*ep11+ep22)-epsthe)
                effopg(6*(kp-1)+4) = effopg(6*(kp-1)+4)+zr(ipoids-1+ip) &
                                     *x3*(h/2.d0)*(c*(ep11+nu*ep22)-epsthe)
                effopg(6*(kp-1)+5) = effopg(6*(kp-1)+5)+zr(ipoids-1+ip) &
                                     *x3*(h/2.d0)*(c*(nu*ep11+ep22)-epsthe)
            else
                effopg(6*(kp-1)+1) = effopg(6*(kp-1)+1)+zr(ipoids-1+ip) &
                                     *(h/2.d0)*(c*ep11-epsthe)
                effopg(6*(kp-1)+2) = effopg(6*(kp-1)+2)+zr(ipoids-1+ip) &
                                     *(h/2.d0)*(c*nu*ep11-epsthe)
                effopg(6*(kp-1)+4) = effopg(6*(kp-1)+4)+zr(ipoids-1+ip) &
                                     *x3*(h/2.d0)*(c*ep11-epsthe)
                effopg(6*(kp-1)+5) = effopg(6*(kp-1)+5)+zr(ipoids-1+ip) &
                                     *x3*(h/2.d0)*(c*nu*ep11-epsthe)
            end if
!
        end do
        effopg(6*(kp-1)+3) = 0.d0
        effopg(6*(kp-1)+6) = 0.d0
!
    end do
!
    call ppgan2(jgano, 1, 6, effopg, zr(ieffor))
end subroutine
