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
subroutine usupru(vusurt, vusuro, nbinst, prust)
    implicit none
!     CALCULE LA PROFONDEUR D'USURE
!
! IN  : VUSURT : VOLUME USE TUBE A CHAQUE INSTANT
! IN  : VUSURO : VOLUME USE OBSTACLE A CHAQUE INSTANT
! IN  : NBINST : NOMBRE D'INSTANTS
! OUT : PRUST  : PROFONDEUR D'USURE DU TUBE POUR CHAQUE INSTANT
!-----------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/usubis.h"
#include "asterfort/usufon.h"
#include "asterfort/usunew.h"
#include "asterfort/utmess.h"
    real(kind=8) :: lsup, vusurt(*), vusuro(*), prust(*), para(7)
    character(len=4) :: crit
    character(len=24) :: type, typ1, typ2
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ire1, ire2, iret, n0, n1, n2, n3
    integer(kind=8) :: n4, n5, nbinst, n6
    real(kind=8) :: aimp, angl, cst1, cst2, de, depi, des3
    real(kind=8) :: des5, df, epsi, rapp, rayo
    real(kind=8) :: rayt, resu, un, uns3, uns5, v1, v2
    real(kind=8) :: vulim, x1, x11, x2, xla, zero, para_cont(6)
!-----------------------------------------------------------------------
    type = ' '
    zero = 0.d0
    un = 1.d0
    de = 2.d0
    uns3 = un/3.d0
    des3 = de/3.d0
    uns5 = un/5.d0
    des5 = de/5.d0
    depi = r8depi()
    crit = 'RELA'
    epsi = 1.d-06
    para(1) = zero
    para(2) = zero
    para(3) = zero
    para(4) = zero
    para(5) = zero
    para(6) = zero
    para(7) = zero
    angl = zero
!
    call getvtx(' ', 'CONTACT', scal=type, nbret=n0)
!
    if (n0 .ne. 0) then
        call getvr8(' ', 'RAYON_MOBILE', scal=para_cont(1), nbret=n1)
        call getvr8(' ', 'RAYON_OBST', scal=para_cont(2), nbret=n2)
        call getvr8(' ', 'LARGEUR_OBST', scal=para_cont(3), nbret=n3)
        call getvr8(' ', 'ANGL_INCLI', scal=para_cont(4), nbret=n4)
        call getvr8(' ', 'ANGL_ISTHME', scal=para_cont(5), nbret=n5)
        call getvr8(' ', 'ANGL_IMPACT', scal=para_cont(6), nbret=n6)
    end if
!
!
!     --- TUBE - BARRE ANTI VIBRATOIRE ---
    if (type(1:8) .eq. 'TUBE_BAV') then
        if (n1 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'RAYON_MOBILE'])
        if (n3 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'LARGEUR_OBST'])
        if (n2 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'RAYON_OBST'])
        if (n5 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_ISTHME'])
!
        rayt = para_cont(1)
        lsup = para_cont(3)
!
        if (n6 .ne. 0) then
            aimp = para_cont(6)
            rapp = cos(aimp*r8dgrd())
        else
            rapp = un
        end if
        if (n4 .eq. 0) then
            cst1 = (un/(de*rayt))**uns3
            cst2 = 3.d0/(4.d0*lsup)
            do i = 1, nbinst
                v1 = vusurt(i)*rapp+vusuro(i)*rapp
                v2 = vusurt(i)*rapp/v1
                prust(i) = v2*cst1*((cst2*v1)**des3)
            end do
        else
            angl = para_cont(4)*r8dgrd()
            xla = lsup*angl
            cst1 = (un/(de*rayt))**uns5
            cst2 = 15.d0*angl/8.d0
            para(1) = rayt
            para(2) = lsup
            para(3) = angl
            x1 = xla
            x2 = rayt
            do i = 1, nbinst
                v1 = vusurt(i)*rapp+vusuro(i)*rapp
                v2 = vusurt(i)*rapp/v1
                prust(i) = v2*cst1*((cst2*v1)**des5)
                if (prust(i) .gt. xla) then
                    para(4) = vusurt(i)*rapp
                    para(5) = vusuro(i)*rapp
                    if (prust(i) .ge. x2) then
                        call utmess('A', 'PREPOST4_83')
                        prust(i) = 9999.d0
                        goto 12
                    end if
                    call usunew(type, para, crit, epsi, x1, &
                                x2, resu, iret)
                    if (iret .eq. 0) then
                        prust(i) = resu
                        x1 = resu
                    else
                        prust(i) = 9999.d0
                    end if
                end if
12              continue
            end do
        end if
!
!     --- TUBE - TROU CIRCULAIRE ---
    else if (type(1:12) .eq. 'TUBE_ALESAGE') then
        if (n1 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'RAYON_MOBILE'])
        if (n3 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'LARGEUR_OBST'])
        if (n5 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_ISTHME'])
        if (n6 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_IMPACT'])
!
        rayt = para_cont(1)
        rayo = para_cont(2)
        lsup = para_cont(3)
        angl = para_cont(4)
!
        para(1) = para_cont(1)
        para(2) = para_cont(2)
        para(3) = para_cont(3)
        if (n4 .ne. 0) angl = para_cont(4)*r8dgrd()
        para(4) = angl
!
        if (n2 .eq. 0) then
            do i = 1, nbinst
                prust(i) = vusurt(i)/(depi*lsup*rayt)
            end do
        else
            x1 = zero
            x2 = rayt
            if (n4 .eq. 0) then
                do i = 1, nbinst
                    para(5) = vusurt(i)
                    para(6) = vusuro(i)
                    call usubis(type, para, crit, epsi, x1, &
                                x2, resu, iret)
                    if (iret .eq. 0) then
                        if (resu .ge. x2) then
                            call utmess('A', 'PREPOST4_83')
                            prust(i) = 9999.d0
                            goto 22
                        end if
                        prust(i) = resu
                        x1 = resu
                    else
                        prust(i) = 9999.d0
                    end if
22                  continue
                end do
            else
!            --- CAS 3 OU D < L * THETA ---
                typ1 = 'TUBE_ALESAG_3A'
!            --- CAS 3 OU D > L * THETA ---
                typ2 = 'TUBE_ALESAG_3B'
!
                xla = lsup*angl
                x1 = zero
                x2 = de*rayo
                do i = 1, nbinst
                    para(5) = vusurt(i)
                    para(6) = vusuro(i)
                    call usubis(typ1, para, crit, epsi, x1, &
                                x2, resu, ire1)
                    if (ire1 .eq. 0) then
                        if (resu .gt. xla) then
                            call usubis(typ2, para, crit, epsi, x1, &
                                        x2, resu, ire2)
                            if (ire2 .eq. 0) then
                                prust(i) = resu
                                x1 = resu
                            else
                                prust(i) = 9999.d0
                            end if
                        else
                            prust(i) = resu
                            x1 = resu
                        end if
                    else
                        prust(i) = 9999.d0
                    end if
                end do
            end if
        end if
!
!     --- TUBE - TROU QUADRIFOLIE OU TRIFOLIE ---
    elseif (type(1:11) .eq. 'TUBE_4_ENCO' .or. type(1:11) .eq. &
            'TUBE_3_ENCO') then
!
        if (n1 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'RAYON_MOBILE'])
        if (n2 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'RAYON_OBST'])
        if (n3 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'LARGEUR_OBST'])
        if (n6 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_IMPACT'])
!
        para(1) = para_cont(1)
        para(2) = para_cont(2)
        para(3) = para_cont(3)
        if (n4 .ne. 0) para(4) = para_cont(4)*r8dgrd()
        if (n5 .ne. 0) para(7) = para_cont(5)*r8dgrd()
!
        x1 = zero
        x2 = para(1)
        if (n4 .eq. 0) then
            do i = 1, nbinst
                if (type(1:11) .eq. 'TUBE_4_ENCO') then
                    para(5) = vusurt(i)/de
                    para(6) = vusuro(i)/de
                else
                    para(5) = vusurt(i)
                    para(6) = vusuro(i)
                end if
                call usubis(type, para, crit, epsi, x1, &
                            x2, resu, iret)
                if (iret .eq. 0) then
                    if (resu .ge. x2) then
                        call utmess('A', 'PREPOST4_83')
                        prust(i) = 9999.d0
                        goto 30
                    end if
                    prust(i) = resu
                    x1 = resu
                else
                    prust(i) = 9999.d0
                end if
30              continue
            end do
        else
!           --- CAS 2 OU D < L * THETA ---
            typ1 = 'TUBE_ENCO_2A'
!           --- CAS 2 OU D > L * THETA ---
            typ2 = 'TUBE_ENCO_2B'
!
            xla = para(3)*para(4)
            x1 = zero
            x2 = de*para(2)
            do i = 1, nbinst
                if (type(1:11) .eq. 'TUBE_4_ENCO') then
                    para(5) = vusurt(i)/de
                    para(6) = vusuro(i)/de
                else
                    para(5) = vusurt(i)
                    para(6) = vusuro(i)
                end if
                call usubis(typ1, para, crit, epsi, x1, &
                            x2, resu, ire1)
                if (ire1 .eq. 0) then
                    if (resu .gt. xla) then
                        call usubis(typ2, para, crit, epsi, x1, &
                                    x2, resu, ire2)
                        if (ire2 .eq. 0) then
                            prust(i) = resu
                            x1 = resu
                        else
                            prust(i) = 9999.d0
                        end if
                    else
                        prust(i) = resu
                        x1 = resu
                    end if
                else
                    prust(i) = 9999.d0
                end if
            end do
        end if
!
!     --- TUBE - TUBE ---
    else if (type(1:9) .eq. 'TUBE_TUBE') then
!
        if (n1 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'RAYON_MOBILE'])
        if (n2 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'RAYON_OBST'])
        if (n3 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'LARGEUR_OBST'])
        if (n5 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_ISTHME'])
        if (n6 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_IMPACT'])
!
        rayt = para_cont(1)
        if (n4 .ne. 0) angl = para_cont(4)
!
        cst1 = (un/(de*rayt))**uns5
        cst2 = 15.d0*angl*r8dgrd()/8.d0
        do i = 1, nbinst
            prust(i) = cst1*((cst2*vusurt(i))**des5)
        end do
!
!     --- GRAPPE - ALESAGE ---
    else if (type(1:14) .eq. 'GRAPPE_ALESAGE') then
!
        if (n1 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'RAYON_MOBILE'])
        if (n2 .eq. 0) call utmess('F', 'PREPOST4_1', nk=2, valk=[type, 'RAYON_OBST'])
        if (n3 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'LARGEUR_OBST'])
        if (n4 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_INCLI'])
        if (n5 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_ISTHME'])
        if (n6 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_IMPACT'])
!
        para(1) = para_cont(1)
        para(2) = para_cont(2)
!
        x11 = zero
        x2 = para(2)
        do i = 1, nbinst
            prust(i) = 9999.d0
            para(5) = vusurt(i)
            call usubis(type, para, crit, epsi, x11, &
                        x2, resu, iret)
            if (iret .eq. 0) then
                if (resu .ge. x2) then
                    call utmess('A', 'PREPOST4_83')
                    goto 50
                end if
                prust(i) = resu
                x11 = resu
            end if
50          continue
        end do
!
!     --- GRAPPE - 1 ENCOCHE ---
!     --- GRAPPE - 2 ENCOCHE ---
    elseif (type(1:13) .eq. 'GRAPPE_1_ENCO' .or. type(1:13) .eq. &
            'GRAPPE_2_ENCO') then
!
        if (n1 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'RAYON_MOBILE'])
        if (n2 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'RAYON_OBST'])
        if (n3 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'LARGEUR_OBST'])
        if (n4 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_INCLI'])
        if (n5 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_ISTHME'])
        if (n6 .ne. 0) call utmess('A', 'PREPOST4_2', nk=2, valk=[type, 'ANGL_IMPACT'])
!
        if (type(1:13) .eq. 'GRAPPE_2_ENCO') then
            para(1) = -48.89d+03/11.d0
            para(2) = 106.03d0/11.d0
            para(3) = -0.88d-03/11.d0
        else
            para(1) = -0.5d0*48.89d+03/11.d0
            para(2) = 0.5d0*106.03d0/11.d0
            para(3) = -0.5d0*0.88d-03/11.d0
        end if
        x11 = zero
        x2 = 0.00144d0
        para(5) = zero
        call usufon(type, para, x2, vulim, df)
        do i = 1, nbinst
            prust(i) = 9999.d0
            para(5) = vusurt(i)
            if (vusurt(i) .gt. vulim) goto 62
            call usunew(type, para, crit, epsi, x11, &
                        x2, resu, iret)
            if (iret .eq. 0) then
                if (resu .ge. x2) then
                    call utmess('A', 'PREPOST4_83')
                    goto 62
                end if
                prust(i) = resu
                x11 = resu
            end if
62          continue
        end do
!
    end if
!
end subroutine
