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
subroutine usuvus(puusur, vusur, nbinst, temps, isupp, &
                  nbpt, fn, vg, iret)
    implicit none
!     CALCULE LE VOLUME USE
!
! IN  : PUUSUR : PUISSANCE USURE
! OUT : VUSUR  : VOLUME USE
! IN  : NBINST : NOMBRE D'INSTANTS
! IN  : TEMPS  : LES INSTANTS
! VAR : ISUPP  : = 1, CALCULE LE VOLUME USE MOBILE
!                = 2, CALCULE LE VOLUME USE OBSTACLE
!                NE CALCULE PAS LE VOLUME USE OBSTACLE, ISUPP = 0
!-----------------------------------------------------------------------
#include "asterc/getfac.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/iunifi.h"
#include "asterfort/usuban.h"
#include "asterfort/usukwu.h"
#include "asterfort/utmess.h"
    real(kind=8) :: vusur(*), temps(*), para(7), fn(*), vg(*)
    character(len=8) :: k8b
    character(len=24) :: loi, mate
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifires, iret, isupp, n1, n2
    integer(kind=8) :: n3, n4, n5, n6, nbinst, nbpt, nn
!
    real(kind=8) :: puusur, t, v0, w, x1, xa, xb
    real(kind=8) :: xd, xk, xn, xs
!-----------------------------------------------------------------------
    ifires = iunifi('RESULTAT')
!
    call getvtx(' ', 'LOI_USURE', scal=loi, nbret=n1)
    iret = 0
    k8b = ' '
!
! **********************************************************************
!                 M O D E L E     A R C H A R D
! **********************************************************************
!
    if (loi(1:7) .eq. 'ARCHARD') then
        if (isupp .eq. 1) then
            write (ifires, 1000)
            call getvr8('MOBILE', 'COEF_USURE', iocc=1, scal=xk, nbret=n1)
            if (n1 .eq. 0) then
                call getvtx(' ', 'MATER_USURE', scal=mate, nbret=n2)
                call usuban(mate, isupp, para, iret)
                xk = para(1)
            end if
            write (ifires, 2100)
        else if (isupp .eq. 2) then
            call getvr8('OBSTACLE', 'COEF_USURE', iocc=1, scal=xk, nbret=n1)
            if (n1 .eq. 0) then
                call getvtx(' ', 'USURE_OBST', scal=k8b, nbret=n2)
                if (k8b(1:3) .eq. 'OUI') then
                    call getvtx(' ', 'MATER_USURE', scal=mate, nbret=n3)
                    call usuban(mate, isupp, para, iret)
                    xk = para(1)
                else
                    isupp = 0
                    goto 999
                end if
            end if
            write (ifires, 2200)
        end if
        write (ifires, 2010) xk
        do i = 1, nbinst
            vusur(i) = xk*puusur*temps(i)
        end do
!
! **********************************************************************
!                 M O D E L E     K W U _ E P R I
! **********************************************************************
!
    else if (loi(1:8) .eq. 'KWU_EPRI') then
        if (isupp .eq. 1) then
            write (ifires, 1010)
            call getvr8('MOBILE', 'COEF_USURE', iocc=1, scal=para(1), nbret=n1)
            call getvr8('MOBILE', 'COEF_FNOR', iocc=1, scal=para(2), nbret=n2)
            call getvr8('MOBILE', 'COEF_VTAN', iocc=1, scal=para(3), nbret=n3)
            call getvr8('MOBILE', 'COEF_K', iocc=1, scal=para(4), nbret=n4)
            call getvr8('MOBILE', 'COEF_C', iocc=1, scal=para(5), nbret=n5)
            if (n4 .eq. 0) para(4) = 5.d0
            if (n5 .eq. 0) para(5) = 10.d0
            call getvtx(' ', 'MATER_USURE', scal=mate, nbret=n6)
            if (n6 .ne. 0) then
                call usuban(mate, isupp, para, iret)
            end if
            write (ifires, 2100)
        else if (isupp .eq. 2) then
            call getvr8('OBSTACLE', 'COEF_USURE', iocc=1, scal=para(1), nbret=n1)
            call getvr8('OBSTACLE', 'COEF_FNOR', iocc=1, scal=para(2), nbret=n2)
            call getvr8('OBSTACLE', 'COEF_VTAN', iocc=1, scal=para(3), nbret=n3)
            call getvr8('OBSTACLE', 'COEF_K', iocc=1, scal=para(4), nbret=n4)
            call getvr8('OBSTACLE', 'COEF_C', iocc=1, scal=para(5), nbret=n5)
            if (n4 .eq. 0) para(4) = 5.d0
            if (n5 .eq. 0) para(5) = 10.d0
            call getvtx(' ', 'MATER_USURE', scal=mate, nbret=n6)
            if (n6 .ne. 0) then
                call getvtx(' ', 'USURE_OBST', scal=k8b, nbret=n2)
                if (k8b(1:3) .eq. 'OUI') then
                    call usuban(mate, isupp, para, iret)
                else
                    isupp = 0
                    goto 999
                end if
            end if
            nn = n1+n2+n3+n4+n5
            if (nn .eq. 0) then
                isupp = 0
                goto 999
            end if
            write (ifires, 2200)
        end if
        write (ifires, 2010) para(1)
        write (ifires, 2050) para(3)
        write (ifires, 2060) para(2)
        write (ifires, 2070) para(4)
        write (ifires, 2080) para(5)
        call usukwu(nbpt, fn, vg, para, w, &
                    iret)
        if (iret .eq. 10) then
            call utmess('F', 'PREPOST4_85')
        end if
        do i = 1, nbinst
            vusur(i) = para(1)*w*puusur*temps(i)
        end do
!
! **********************************************************************
!                 M O D E L E     E D F _ M Z
! **********************************************************************
!
    else if (loi(1:6) .eq. 'EDF_MZ') then
        if (isupp .eq. 1) then
            write (ifires, 1020)
            call getvr8('MOBILE', 'COEF_S', iocc=1, scal=xs, nbret=n1)
            call getvr8('MOBILE', 'COEF_B', iocc=1, scal=xb, nbret=n2)
            call getvr8('MOBILE', 'COEF_N', iocc=1, scal=xn, nbret=n3)
            call getvr8('MOBILE', 'COEF_USURE', iocc=1, scal=xa, nbret=n4)
            if (n1 .eq. 0) xs = 1.14d-16
            if (n2 .eq. 0) xb = 1.2d0
            if (n3 .eq. 0) xn = 2.44d-08
            if (n4 .eq. 0) xa = 1.d-13
            call getvtx(' ', 'MATER_USURE', scal=mate, nbret=n5)
            if (n5 .ne. 0) then
                call usuban(mate, isupp, para, iret)
                xs = para(1)
                xb = para(2)
                xn = para(3)
                xa = para(4)
            end if
            write (ifires, 2100)
        else if (isupp .eq. 2) then
            call getvr8('OBSTACLE', 'COEF_S', iocc=1, scal=xs, nbret=n1)
            call getvr8('OBSTACLE', 'COEF_B', iocc=1, scal=xb, nbret=n2)
            call getvr8('OBSTACLE', 'COEF_N', iocc=1, scal=xn, nbret=n3)
            call getvr8('OBSTACLE', 'COEF_USURE', iocc=1, scal=xa, nbret=n4)
            if (n1 .eq. 0) xs = 1.14d-16
            if (n2 .eq. 0) xb = 1.2d0
            if (n3 .eq. 0) xn = 2.44d-08
            if (n4 .eq. 0) xa = 1.d-13
            call getvtx(' ', 'MATER_USURE', scal=mate, nbret=n5)
            if (n5 .ne. 0) then
                call getvtx(' ', 'USURE_OBST', scal=k8b, nbret=n6)
                if (k8b(1:3) .eq. 'OUI') then
                    call usuban(mate, isupp, para, iret)
                    xs = para(2)
                    xb = para(3)
                    xn = para(4)
                    xa = para(1)
                else
                    isupp = 0
                    goto 999
                end if
            end if
            call getfac('OBSTACLE', n6)
            nn = n1+n2+n3+n4+n5+n6
            if (nn .eq. 0) then
                isupp = 0
                goto 999
            end if
            write (ifires, 2200)
        end if
        write (ifires, 2010) xa
        write (ifires, 2020) xs
        write (ifires, 2030) xb
        write (ifires, 2040) xn
        v0 = xa*(puusur**xb)
        xd = xs/v0
        if (xd .gt. 1.d0) then
            iret = 10
            call utmess('I', 'PREPOST4_86')
            call utmess('I', 'PREPOST4_87')
            goto 999
        end if
        x1 = (1.d0-xd)/xn
        do i = 1, nbinst
            t = temps(i)
            vusur(i) = v0*(xd*t+x1*(1.d0-exp(-xn*t)))
        end do
!
    end if
!
1000 format(/, '******* MODELE ARCHARD *******')
1010 format(/, '******* MODELE KWU_EPRI *******')
1020 format(/, '******* MODELE EDF_MZ *******')
2100 format(/, '===> COEFFICIENT(S) UTILISE(S) POUR LE MOBILE :')
2200 format(/, '===> COEFFICIENT(S) UTILISE(S) POUR L''OBSTACLE :')
2010 format(1p, 4x, '       COEFFICIENT D''USURE : ', e12.5)
2020 format(1p, 4x, '                     SEUIL : ', e12.5)
2030 format(1p, 4x, '                  EXPOSANT : ', e12.5)
2040 format(1p, 4x, '    TAUX DE RALENTISSEMENT : ', e12.5)
2050 format(1p, 4x, ' COEFFICIENT DE GLISSEMENT : ', e12.5)
2060 format(1p, 4x, '      COEFFICIENT D''IMPACT : ', e12.5)
2070 format(1p, 4x, '               CONSTANTE K : ', e12.5)
2080 format(1p, 4x, '               CONSTANTE C : ', e12.5)
!
999 continue
end subroutine
