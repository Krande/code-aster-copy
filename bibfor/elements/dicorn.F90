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
subroutine dicorn(irmetg, nbt, neq, iterat, icodma, &
                  ul, dul, utl, sim, varim, &
                  klv, klv2, varip)
! ----------------------------------------------------------------------
    implicit none
#include "asterfort/dicor0.h"
#include "asterfort/dicor2.h"
#include "asterfort/dicor3.h"
#include "asterfort/dicor4.h"
#include "asterfort/dicor5.h"
#include "asterfort/dikfin.h"
#include "asterfort/dikini.h"
#include "asterfort/pmavec.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/vecma.h"
    integer(kind=8) :: irmetg, nbt, neq, iterat, icodma
    real(kind=8) :: ul(neq), dul(neq), utl(neq)
    real(kind=8) :: sim(neq), varim(7)
    real(kind=8) :: klv(nbt), klv2(nbt), varip(7)
!
!     RELATION DE COMPORTEMENT "ASSE_CORN" (CORNIERE).
!
! ----------------------------------------------------------------------
!
! IN  : IRMETG : VAUT 1 SI ON CALCULE L'OPTION "RIGI_MECA_TANG"
!       NBT    : NOMBRE DE VALEURS POUR LA DEMI-MATRICE
!       NEQ    : NOMBRE DE DDL DE L'ELEMENT
!       ITERAT : NUMERO DE L'ITERATION DE NEWTON
!       ICODMA : ADRESSE DU MATERIAU CODE
!       UL     : DEPLACEMENT PRECEDENT REPERE LOCAL (DIM NEQ)
!       DUL    : INCREMENT DE DEPLACEMENT REPERE LOCAL (DIM NEQ)
!       UTL    : DEPLACEMENT COURANT REPERE LOCAL (DIM NEQ)
!       SIM    : EFFORTS GENERALISES A L'INSTANT PRECEDENT (DIM NEQ)
!       VARIM$ : VARIABLES INTERNES A L'INSTANT PRECEDENT (7 VALEURS)
!
! OUT : KLV    :                                (DIM NBT)
!       KLV2   :                                (DIM NBT)
!       VARIP$ : VARIABLES INTERNES REACTUALISEES (7 VALEURS)
!
!***************** DECLARATION DES VARIABLES LOCALES *******************
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, nbpar, nbre1
    real(kind=8) :: a1, a2, c1, c2, dbar1, dbar2, dmsdt
    real(kind=8) :: dmsdt2, dnsdt, dnsdt2, dnsdu, dnsdu2, dry2, dryr
    real(kind=8) :: dryu1, dryu2, du2, dur, dxu1, dxu2, feq1
    real(kind=8) :: feq2, g1, g2, p1, p2, pi, plouf
    real(kind=8) :: rg1, rg2, t2, test, ti, tr2, tt
    real(kind=8) :: ttot, u2, ui, ur2, utot, uu, valpar
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    parameter(nbre1=15)
    real(kind=8) :: nu1, mu1, nu2, mu2, ky, kz, krx, krz, rp0
    real(kind=8) :: si(12), k01(78), k02(78), klc(144), valre1(nbre1)
    integer(kind=8) :: codre1(nbre1), kpg, spt
    character(len=8) :: nompar, nomre1(nbre1), fami, poum
!
!************ FIN DES DECLARATIONS DES VARIABLES LOCALES ***************
!
!****************************** DATA ***********************************
!
    data nomre1/'NU_1', 'MU_1', 'DXU_1', 'DRYU_1', 'C_1',&
     &            'NU_2', 'MU_2', 'DXU_2', 'DRYU_2', 'C_2',&
     &            'KY', 'KZ', 'KRX', 'KRZ', 'R_P0'/
!
! ----------------------------------------------------------------------
! --- DEFINITION DES PARAMETRES
!
    zero = 0.d0
    nbpar = 0
    nompar = ' '
    valpar = 0.d0
    call r8inir(nbre1, zero, valre1, 1)
!
! --- CARACTERISTIQUES DU MATERIAU
!    (LES DEFINITIONS DE DRYU1 ET DRYU2 SYMETRISENT LA MATRICE)
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, icodma, &
                ' ', 'ASSE_CORN', nbpar, nompar, [valpar], &
                nbre1, nomre1, valre1, codre1, 1)
!
    nu1 = valre1(1)
    mu1 = valre1(2)
    dxu1 = valre1(3)
    c1 = valre1(5)
    nu2 = valre1(6)
    mu2 = valre1(7)
    dxu2 = valre1(8)
    c2 = valre1(10)
    ky = valre1(11)
    kz = valre1(12)
    krx = valre1(13)
    krz = valre1(14)
!      DRYU1 = NU1 * DXU1 / MU1
!      DRYU2 = NU2 * DXU2 / MU2
    dryu1 = valre1(4)
    dryu2 = valre1(9)
    rp0 = valre1(15)
!      RP0   = 1.D4
!
! --- CONSTANTES DE LA RELATION DE COMPORTEMENT
!
    a1 = 1.d0
    a2 = 1.d0
    dbar1 = c1**(a1+1.d0)/(1.d0-c1**a1)
    dbar2 = c2**(a2+1.d0)/(1.d0-c2**a2)
!
! --- ECRITURE DANS LE REPERE LOCAL DE K01 ET K02 (MATRICES DE
!     RAIDEUR TANGENTE INITIALES POUR LES DEUX MECANISMES)
!
    call dikini(nbt, nu1, mu1, dxu1, dryu1, &
                nu2, mu2, dxu2, dryu2, ky, &
                kz, krx, krz, k01, k02, &
                rp0)
!
! ======================================================================
!                  DEBUT DU TRAITEMENT DE L'ASSEMBLAGE
! ======================================================================
!
! --- DUR  : INCREMENT DE LONGUEUR DANS L'AXE DE L'ELEMENT
! --- DRYR : INCREMENT DE ROTATION
!
    dur = dul(7)-dul(1)
    dryr = dul(11)-dul(5)
    uu = utl(7)-utl(1)
    tt = utl(11)-utl(5)
    ui = ul(7)-ul(1)
    ti = ul(11)-ul(5)
!      INDECH = 0
!
! -*-*-*-*       TEST POUR SAVOIR SI L'ON DECHARGE OU NON      *-*-*-*-*
!
!      IF ((((UU*DUR).GT.0.D0.AND.(UI*DUR).GE.0.D0).OR.
!     &   ((TT*DRYR).GT.0.D0.AND.(TI*DRYR).GE.0.D0))) INDECH = 1
!
    if (irmetg .ne. 1) then
!
! ======================================================================
!                       FULL_MECA
! ======================================================================
!
        varip(2) = 0.d0
!
! -*-*-*-* TEST POUR DETERMINER LE MECANISME OU L'ON SE TROUVE *-*-*-*-*
!
        if (varim(1) .le. 1.d0 .or. varim(3) .eq. 1.d0) then
!
! ====================================
! ====== ON EST EN MECANISME 1 =======
! ====================================
!
            call vecma(k01, nbt, klc, neq)
            call pmavec('ZERO', neq, klc, dul, si)
            pi = varim(1)
!
! ****** TEST SUR LE NUMERO D'ITERATION
!
            if (iterat .eq. 1) then
!
! ****** CAS DE LA PREMIERE ITERATION
!
                p1 = varim(1)
                g1 = dbar1*p1
                rg1 = 0.5d0*(-g1+sqrt(g1**2+4.d0*g1))
!
! **** TEST SUR LA POSITION PAR RAPPORT A LA SLU1
!
!
! **** ON EST SUR LA SLU1
!
                if (varim(1) .eq. 0.d0) then
                    dnsdu2 = k01(1)
                    dmsdt2 = k01(15)
                else
                    dnsdu2 = rg1*nu1/dxu1/p1
                    if (dur .eq. 0.d0) dnsdu2 = k01(1)
                    dmsdt2 = rg1*mu1/dryu1/p1
                    if (dryr .eq. 0.d0) dmsdt2 = k01(15)
                end if
!
                dnsdt2 = 0.d0
                si(7) = sim(7)+dnsdu2*dur
                si(11) = sim(11)+dmsdt2*dryr
                si(1) = -si(7)
                si(5) = -si(11)
!
                feq1 = sqrt((si(7)/nu1)**2+(si(11)/mu1)**2)
!
! ** TEST DE CHANGEMENT DE MECANISME
!
                if (feq1 .lt. c1) then
!
! ** ON RESTE EN MECANISME 1
!
                    p1 = feq1**2/(1.d0-feq1)/dbar1
                    u2 = p1*dxu1*si(7)/nu1/feq1
                    t2 = p1*dryu1*si(11)/mu1/feq1
                    utot = u2+varim(4)
                    ttot = t2+varim(5)
!
                    if (dur .ne. 0.d0) dnsdu2 = si(7)/utot
                    if (dur .eq. 0.d0) dnsdu2 = k01(1)
                    if (dryr .ne. 0.d0) dmsdt2 = si(11)/ttot
                    if (dryr .eq. 0.d0) dmsdt2 = k01(15)
                    dnsdt2 = 0.d0
                    si(7) = dnsdu2*uu
                    si(11) = dmsdt2*tt
                    varip(1) = p1
                    varip(2) = varim(2)
                    varip(3) = 1.0d0
!
                    call dicor3(k01, dur, dryr, sim, si, &
                                dnsdu, dmsdt, dnsdt)
!
                    do i = 4, 7
                        varip(i) = varim(i)
                    end do
                else
!
! ** ON PASSE EN MECANISME 2
!
                    u2 = ui-varim(4)
                    t2 = ti-varim(5)
                    call dicor4(k02, sim, si, pi, u2, &
                                t2, dxu1, dxu2, dryu1, dryu2, &
                                nu1, nu2, mu1, mu2, feq1, &
                                c1, dbar2, uu, tt, dur, &
                                dryr, p2, utot, ttot, dnsdu, &
                                dmsdt, dnsdt, dnsdu2, dmsdt2, dnsdt2)
                    varip(4) = utot-si(7)/k02(1)
                    varip(5) = ttot-si(11)/k02(15)
                    varip(6) = si(7)
                    varip(7) = si(11)
                    u2 = utot-varim(4)
                    t2 = ttot-varim(5)
                    varip(1) = sqrt((u2/dxu1)**2+(t2/dryu1)**2)
                    varip(2) = p2
                    varip(3) = 2.d0
!
                end if
!
!             ELSE
!
! **** ON EST SOUS LA SLU1
!
!
!             ENDIF
!
            else if (iterat .ge. 2) then
!
! ****** CAS DES ITERATIONS 2 ET SUIVANTES
!
                u2 = uu-varim(4)
                t2 = tt-varim(5)
                varip(1) = sqrt((u2/dxu1)**2+(t2/dryu1)**2)
                p1 = varip(1)
!
                if (p1 .le. 1.d0) then
!
! **** ON RESTE EN MECANISME 1
!
                    g1 = dbar1*p1
                    rg1 = 0.5d0*(-g1+sqrt(g1**2+4.d0*g1))
                    dnsdu2 = rg1*nu1/dxu1/p1
                    if (dur .eq. 0.d0) dnsdu2 = k01(1)
                    dmsdt2 = rg1*mu1/dryu1/p1
                    if (dryr .eq. 0.d0) dmsdt2 = k01(15)
!
                    dnsdt2 = 0.d0
!
!
                    call dicor2(k01, varim(2), p1, dur, dryr, &
                                dxu1, dryu1, rg1, nu1, mu1, &
                                u2, t2, sim, dnsdu2, dmsdt2, &
                                dnsdt2, varip(2), varip(3), si)
!
                    call dicor3(k01, dur, dryr, sim, si, &
                                dnsdu, dmsdt, dnsdt)
                    do i = 4, 7
                        varip(i) = varim(i)
                    end do
!
                else
!
! **** ON PASSE EN MECANISME 2
!
!
                    g1 = dbar1*varim(1)
                    rg1 = 0.5d0*(-g1+sqrt(g1**2+4.d0*g1))
                    u2 = ui-varim(4)
                    t2 = ti-varim(5)
                    call dicor5(k02, sim, p1, pi, u2, &
                                t2, dxu1, dxu2, dryu1, dryu2, &
                                nu1, nu2, mu1, mu2, c1, &
                                dbar2, uu, tt, dur, dryr, &
                                dnsdu, dmsdt, dnsdt, dnsdu2, dmsdt2, &
                                dnsdt2, si, varip(2), varip(3))
                    varip(4) = uu-si(7)/k02(1)
                    varip(5) = tt-si(11)/k02(15)
                    varip(6) = si(7)
                    varip(7) = si(11)
!
                end if
!
            end if
!
!
! -*-*-*-*-*-*-*-*-*-*-*-* FIN DU MECANISME 1 *-*-*-*-*-*-*-*-*-*-*-*-*
!
        else
!
! ====================================
! ====== ON EST EN MECANISME 2 =======
! ====================================
!
            p2 = varim(2)
            varip(1) = varim(1)
            g2 = dbar2*p2
            rg2 = 0.5d0*(-g2+sqrt(g2**2+4.d0*g2))
! ****** TEST SUR LA POSITION PAR RAPPORT A LA SLU2
!
            if (varim(3) .eq. 2.d0) then
!
! ****** ON EST SUR LA SLU2
!
                dnsdu2 = rg2*nu2/dxu2/p2
                dmsdt2 = rg2*mu2/dryu2/p2
                feq2 = sqrt(((sim(7)+dnsdu2*dur)/nu2)**2+((sim(11)+dmsdt2*dryr)/mu2)**2)
                if (feq2 .lt. rg2) then
                    call dicor0(k02, varim(2), varip(2), varip(3), dnsdu, &
                                dmsdt, dnsdt)
                    call dicor0(k02, varim(2), varip(2), varip(3), dnsdu2, &
                                dmsdt2, dnsdt2)
                    do i = 4, 7
                        varip(i) = varim(i)
                    end do
                else
                    if (iterat .eq. 1) then
                        if (feq2 .ge. c2) then
                            call utmess('I', 'ELEMENTS_26')
                        end if
                        si(7) = sim(7)+dnsdu2*dur
                        si(11) = sim(11)+dmsdt2*dryr
                        varip(2) = feq2**2/(1.d0-feq2)/dbar2
                        ur2 = (varip(2)*si(7)/feq2-p2*sim(7)/rg2)/nu2
                        tr2 = (varip(2)*si(11)/feq2-p2*sim(11)/rg2)/mu2
                        u2 = ur2*dxu2
                        t2 = tr2*dryu2
                        utot = u2+ui
                        ttot = t2+ti
!
                        if (dur .ne. 0.d0) dnsdu2 = si(7)/utot
                        if (dur .eq. 0.d0) dnsdu2 = k02(1)
                        if (dryr .ne. 0.d0) dmsdt2 = si(11)/ttot
                        if (dryr .eq. 0.d0) dmsdt2 = k02(15)
                        dnsdt2 = 0.d0
                        varip(4) = utot-si(7)/k02(1)
                        varip(5) = ttot-si(11)/k02(15)
                        varip(6) = si(7)
                        varip(7) = si(11)
                        si(7) = dnsdu2*uu
                        si(11) = dmsdt2*tt
                        call utmess('I', 'ELEMENTS_27')
                        varip(3) = 2.0d0
!
                        call dicor3(k02, dur, dryr, sim, si, &
                                    dnsdu, dmsdt, dnsdt)
                    else
                        u2 = dur+p2*sim(7)*dxu2/rg2/nu2
                        t2 = dryr+p2*sim(11)*dryu2/rg2/mu2
                        varip(2) = sqrt((u2/dxu2)**2+(t2/dryu2)**2)
                        g2 = dbar2*varip(2)
                        feq2 = 0.5d0*(-g2+sqrt(g2**2+4.d0*g2))
                        dnsdu2 = feq2*nu2/dxu2/varip(2)
                        if (dur .eq. 0.d0) dnsdu2 = k02(1)
                        dmsdt2 = feq2*mu2/dryu2/varip(2)
                        if (dryr .eq. 0.d0) dmsdt2 = k02(15)
                        dnsdt2 = 0.d0
                        si(7) = u2*feq2*nu2/dxu2/varip(2)
                        si(11) = t2*feq2*mu2/dryu2/varip(2)
                        call utmess('I', 'ELEMENTS_27')
                        call dicor3(k02, dur, dryr, sim, si, &
                                    dnsdu, dmsdt, dnsdt)
                        varip(3) = 2.d0
                        varip(4) = uu-si(7)/k02(1)
                        varip(5) = tt-si(11)/k02(15)
                        varip(6) = si(7)
                        varip(7) = si(11)
                    end if
                end if
!
            else if (varim(3) .eq. 0.d0) then
!
! ****** ON EST SOUS LA SLU2
!
                feq2 = sqrt(((sim(7)+k02(1)*dur)/nu2)**2+((sim(11)+k02(15)*dryr)/mu2)**2)
!
! **** TEST POUR SAVOIR SI L'ON RESTE SOUS LA SLU
!
                if (feq2 .le. rg2) then
!
! **** ON RESTE SOUS LA SLU2
!
                    si(7) = sim(7)+k02(1)*dur
                    si(11) = sim(11)+k02(15)*dryr
                    test = varim(6)*si(7)/nu2**2+varim(7)*si(11)/mu2**2
                    if (test .lt. 0.d0) then
                        if (iterat .eq. 1) then
                            feq1 = sqrt((si(7)/nu1)**2+(si(11)/mu1)**2)
                            if (feq1 .ge. c1) then
                                call utmess('I', 'ELEMENTS_28')
                                goto 19
                            end if
                            call utmess('I', 'ELEMENTS_29')
!
! ** ON REPASSE EN MECANISME 1
!
                            p1 = feq1**2/(1.d0-feq1)/dbar1
                            u2 = p1*dxu1*si(7)/nu1/feq1
                            t2 = p1*dryu1*si(11)/mu1/feq1
                            utot = u2+varim(4)
                            ttot = t2+varim(5)
                            du2 = utot-ui
                            dry2 = ttot-ti
                            feq2 = sqrt( &
                                   ((sim(7)+k02(1)*du2)/nu2)**2+((sim(11)+k02(15)*dry2)/mu2 &
                                                                 )**2 &
                                   )
                            if (feq2 .gt. rg2) then
                                call utmess('I', 'ELEMENTS_30')
                            end if
!
                            if (dur .ne. 0.d0) dnsdu2 = si(7)/utot
                            if (dur .eq. 0.d0) dnsdu2 = k01(1)
                            if (dryr .ne. 0.d0) dmsdt2 = si(11)/ttot
                            if (dryr .eq. 0.d0) dmsdt2 = k01(15)
                            dnsdt2 = 0.d0
                            si(7) = dnsdu2*uu
                            si(11) = dmsdt2*tt
                            varip(1) = p1
                            varip(2) = varim(2)
                            varip(3) = 1.0d0
!
                            call dicor3(k01, dur, dryr, sim, si, &
                                        dnsdu, dmsdt, dnsdt)
!
                        else
!
! ****** CAS DES ITERATIONS 2 ET SUIVANTES
!
                            u2 = uu-varim(4)
                            t2 = tt-varim(5)
                            varip(1) = sqrt((u2/dxu1)**2+(t2/dryu1)**2)
                            p1 = varip(1)
                            call utmess('I', 'ELEMENTS_29')
!
                            if (p1 .gt. 1.d0) then
                                call utmess('I', 'ELEMENTS_28')
                                goto 19
                            end if
!
! **** ON EST EN MECANISME 1
!
                            g1 = dbar1*p1
                            rg1 = 0.5d0*(-g1+sqrt(g1**2+4.d0*g1))
                            dnsdu2 = rg1*nu1/dxu1/p1
                            if (dur .eq. 0.d0) dnsdu2 = k01(1)
                            dmsdt2 = rg1*mu1/dryu1/p1
                            if (dryr .eq. 0.d0) dmsdt2 = k01(15)
!
                            dnsdt2 = 0.d0
!
!
                            call dicor2(k01, varim(2), p1, dur, dryr, &
                                        dxu1, dryu1, rg1, nu1, mu1, &
                                        u2, t2, sim, dnsdu2, dmsdt2, &
                                        dnsdt2, varip(2), varip(3), si)
!
                            call dicor3(k01, dur, dryr, sim, si, &
                                        dnsdu, dmsdt, dnsdt)
!
                        end if
                        goto 20
                    end if
19                  continue
                    call dicor0(k02, varim(2), varip(2), varip(3), dnsdu, &
                                dmsdt, dnsdt)
                    call dicor0(k02, varim(2), varip(2), varip(3), dnsdu2, &
                                dmsdt2, dnsdt2)
20                  continue
                    do i = 4, 7
                        varip(i) = varim(i)
                    end do
!
                else
!
! **** ON REVIENT SUR LA SLU2
!
                    if (iterat .eq. 1) then
                        si(7) = sim(7)+k02(1)*dur
                        si(11) = sim(11)+k02(15)*dryr
                        varip(2) = feq2**2/(1.d0-feq2)/dbar2
!
                        ur2 = (varip(2)*si(7)/feq2-p2*varim(6)/rg2)/nu2
                        tr2 = (varip(2)*si(11)/feq2-p2*varim(7)/rg2)/mu2
                        u2 = ur2*dxu2
                        t2 = tr2*dryu2
                        utot = u2+ui+(varim(6)-sim(7))/k02(1)
                        ttot = t2+ti+(varim(7)-sim(11))/k02(15)
!
                        if (dur .ne. 0.d0) dnsdu2 = si(7)/utot
                        if (dur .eq. 0.d0) dnsdu2 = k02(1)
                        if (dryr .ne. 0.d0) dmsdt2 = si(11)/ttot
                        if (dryr .eq. 0.d0) dmsdt2 = k02(15)
                        dnsdt2 = 0.d0
                        varip(4) = utot-si(7)/k02(1)
                        varip(5) = ttot-si(11)/k02(15)
                        varip(6) = si(7)
                        varip(7) = si(11)
                        si(7) = dnsdu2*uu
                        si(11) = dmsdt2*tt
!
                        call utmess('I', 'ELEMENTS_27')
                        varip(3) = 2.0d0
                        call dicor3(k02, dur, dryr, sim, si, &
                                    dnsdu, dmsdt, dnsdt)
                    else
!
                        u2 = dur+p2*varim(6)*dxu2/rg2/nu2-(varim(6)-sim(7))/k02(1)
                        t2 = dryr+p2*varim(7)*dryu2/rg2/mu2-(varim(7)-sim(11))/k02(15)
                        varip(2) = sqrt((u2/dxu2)**2+(t2/dryu2)**2)
                        g2 = dbar2*varip(2)
                        feq2 = 0.5d0*(-g2+sqrt(g2**2+4.d0*g2))
                        dnsdu2 = feq2*nu2/dxu2/varip(2)
                        if (dur .eq. 0.d0) dnsdu2 = k02(1)
                        dmsdt2 = feq2*mu2/dryu2/varip(2)
                        if (dryr .eq. 0.d0) dmsdt2 = k02(15)
                        dnsdt2 = 0.d0
                        si(7) = u2*feq2*nu2/dxu2/varip(2)
                        si(11) = t2*feq2*mu2/dryu2/varip(2)
                        call utmess('I', 'ELEMENTS_27')
                        call dicor3(k02, dur, dryr, sim, si, &
                                    dnsdu, dmsdt, dnsdt)
                        varip(3) = 2.d0
                        varip(4) = uu-si(7)/k02(1)
                        varip(5) = tt-si(11)/k02(15)
                        varip(6) = si(7)
                        varip(7) = si(11)
                    end if
!
                end if
!
            end if
!
! -*-*-*-*-*-*-*-*-*-*-*-* FIN DU MECANISME 2 *-*-*-*-*-*-*-*-*-*-*-*-*
!
        end if
!
    else
!
! ======================================================================
!                             RIGI_MECA_TANG
! ======================================================================
!
        if (varim(1) .le. 1.d0 .or. varim(3) .eq. 1.d0) then
            call dicor0(k01, varim(1), varip(1), plouf, dnsdu, &
                        dmsdt, dnsdt)
            call dicor0(k01, varim(2), varip(2), varip(3), dnsdu2, &
                        dmsdt2, dnsdt2)
            varip(3) = varim(3)
            p1 = varim(1)
            g1 = dbar1*p1
            rg1 = 0.5d0*(-g1+sqrt(g1**2+4.d0*g1))
            if (p1 .ne. 0.d0) dnsdu2 = rg1*nu1/dxu1/p1
            if (p1 .ne. 0.d0) dmsdt2 = rg1*mu1/dryu1/p1
        else
            call dicor0(k02, varim(1), varip(1), plouf, dnsdu, &
                        dmsdt, dnsdt)
            call dicor0(k02, varim(2), varip(2), varip(3), dnsdu2, &
                        dmsdt2, dnsdt2)
            varip(3) = varim(3)
            p2 = varim(2)
            g2 = dbar2*p2
            rg2 = 0.5d0*(-g2+sqrt(g2**2+4.d0*g2))
            if (varim(3) .eq. 2.d0) dnsdu2 = rg2*nu2/dxu2/p2
            if (varim(3) .eq. 2.d0) dmsdt2 = rg2*mu2/dryu2/p2
        end if
        do i = 4, 7
            varip(i) = varim(i)
        end do
!
    end if
!
! ======================================================================
!                         PARAMETRES EN SORTIE
! ======================================================================
!
! --- ECRITURE DE LA MATRICE TANGENTE EN REPERE LOCAL
!
    call dikfin(nbt, dnsdu, dnsdt, dmsdt, dnsdu2, &
                dnsdt2, dmsdt2, ky, kz, krx, &
                krz, klv, klv2)
! ----------------------------------------------------------------------
!
end subroutine
