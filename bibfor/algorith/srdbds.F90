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

subroutine srdbds(nmat, mater, i1, devsig, nvi, vint, para, val, tmp, dbetds, dbetdi)

!

!!!
!!! MODELE LKR : CALCUL DE LA DERIVEE DE BPRIME PAR RAPPORT A I1 ET S
!!!

! ===================================================================================
! IN  : VAL           : ENTIER PRECISANT DILATANCE EN PRE(0) OU POST-PIC(1)
!     : TMP           : TEMPERATURE A L'INSTANT "-" OU "+"
!     : NMAT          : DIMENSION TABLE DES PARAMETRES MATERIAU
!     : MATER(NMAT,2) : PARAMETRES MATERIAU A T+DT
!     : DEVSIG(6)     : DEVIATEUR DES CONTRAINTES
!     : I1            : TRACE DES CONTRAINTES
!     : NVI           : NOMBRE DE VARIABLES INTERNES
!     : VINT(NVI)     : VARIABLES INTERNES
!     : PARA(3)       : VECTEUR CONTENANT AXI, SXI ET MXI
! OUT : DBETDS(6)     : DERIVEE DE BPRIME PAR RAPPORT A DEVSIG (DEVIATEUR)
!     : DBETDI        : DERIVEE DE BPRIME PAR RAPPORT A I1 (TRACE SIGMA)
! ===================================================================================

    implicit none

#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/cos3t.h"
#include "asterfort/srdhds.h"
#include "asterfort/srhtet.h"

    !!!
    !!! Variables globales
    !!!

    integer(kind=8) :: nmat, nvi, iret, val
    real(kind=8) :: mater(nmat, 2), devsig(6), i1, para(3)
    real(kind=8) :: vint(nvi), dbetds(6), dbetdi, tmp

    !!!
    !!! Variables locales
    !!!

    integer(kind=8) :: i, ndt, ndi
    real(kind=8) :: pi, pref, sigc, m3, xi2
    real(kind=8) :: m1, a2, s1, qi, fi, a5, m5, s5, ffp
    real(kind=8) :: xi1, fp, ucar, fact6, fact7
    real(kind=8) :: rho_1, rho_2, alres, rcos3t, rtheta
    real(kind=8) :: c, phi, xip, tiers, a0, m00, s0, m0
    real(kind=8) :: sigmin, sigmax, sii, sigcar, alpha, sigtil, sinpsi
    real(kind=8) :: dsinds(6), dsmids(6), dsmads(6), r0c, dhds(6), dscardi
    real(kind=8) :: dsmidi, dsmadi, dsindi, rho_4, dbdsin, fact3, fact4, fact5
    real(kind=8) :: m10, qi0, xi10, xi20, rq, rm, f3p, cp, phip, sigtilp, dscards(6)
    real(kind=8) :: rx1, rx2, rs, trr, dtmp, spre, spos, temp(6), fact1, fact2
    common/tdim/ndt, ndi

    !!!
    !!! Init
    !!!

    alpha = 0.d0
    pi = 0.d0
    pref = 0.d0
    sigc = 0.d0
    m3 = 0.d0
    xi2 = 0.d0
    m1 = 0.d0
    a2 = 0.d0
    s1 = 0.d0
    qi = 0.d0
    fi = 0.d0
    a5 = 0.d0
    m5 = 0.d0
    s5 = 0.d0
    ffp = 0.d0
    xi1 = 0.d0
    fp = 0.d0
    ucar = 0.d0
    fact6 = 0.d0
    fact7 = 0.d0
    rho_1 = 0.d0
    rho_2 = 0.d0
    alres = 0.d0
    rcos3t = 0.d0
    rtheta = 0.d0
    c = 0.d0
    phi = 0.d0
    xip = 0.d0
    tiers = 0.d0
    a0 = 0.d0
    m00 = 0.d0
    s0 = 0.d0
    m0 = 0.d0
    sigmin = 0.d0
    sigmax = 0.d0
    sii = 0.d0
    sigcar = 0.d0
    sigtil = 0.d0
    sinpsi = 0.d0
    dsinds = 0.d0
    dsmids = 0.d0
    dsmads = 0.d0
    r0c = 0.d0
    dhds = 0.d0
    dscardi = 0.d0
    dsmidi = 0.d0
    dsmadi = 0.d0
    dsindi = 0.d0
    rho_4 = 0.d0
    dbdsin = 0.d0
    fact3 = 0.d0
    fact4 = 0.d0
    fact5 = 0.d0
    m10 = 0.d0
    qi0 = 0.d0
    xi10 = 0.d0
    xi20 = 0.d0
    rq = 0.d0
    rm = 0.d0
    f3p = 0.d0
    cp = 0.d0
    phip = 0.d0
    sigtilp = 0.d0
    dscards = 0.d0
    rx1 = 0.d0
    rx2 = 0.d0
    rs = 0.d0
    trr = 0.d0
    dtmp = 0.d0
    spre = 0.d0
    spos = 0.d0
    temp = 0.d0
    fact1 = 0.d0
    fact2 = 0.d0

    !!!
    !!! Calcul de sii
    !!!

    sii = norm2(devsig(1:ndt))

    !!!
    !!! Recup. de para. du modele
    !!!

    !!! a T0
    pi = r8pi()
    pref = mater(1, 2)
    sigc = mater(3, 2)
    a0 = 5.d-1
    a2 = mater(8, 2)
    m00 = mater(9, 2)
    m10 = mater(10, 2)
    qi0 = mater(11, 2)
    xi10 = mater(12, 2)
    xi20 = mater(13, 2)
    fp = mater(15, 2)
    rho_1 = mater(18, 2)
    rho_2 = mater(19, 2)
    rho_4 = mater(20, 2)
    rq = mater(21, 2)
    rm = mater(22, 2)
    rs = mater(23, 2)
    rx1 = mater(24, 2)
    rx2 = mater(25, 2)

    !!! Temperatures
    trr = mater(8, 1)

    !!! a T
    if ((tmp .ge. trr) .and. (trr .gt. 0.d0)) then
        qi = qi0*(1.d0-rq*log(tmp/trr))
        dtmp = tmp-trr
    else
        qi = qi0
        dtmp = 0.d0
    end if

    m0 = m00*exp(-rm*(dtmp**2.d0))
    m1 = m10*exp(-rm*(dtmp**2.d0))
    s0 = (m0*1.d-1/(1.d0-1.d-1**2.d0))**2.d0
    s1 = 1.d0*exp(-rs*(dtmp**2.d0))
    xi1 = xi10*exp(rx1*tmp)
    xi2 = xi20*exp(rx2*tmp)
    fi = qi/sigc
    m3 = m1*fi/(fi**2.d0-s1)

    !!! Calculs de s5, m5 et a5
    if (fp .le. 0.d0) then
        s5 = s0
        m5 = m0
        a5 = a0
    else
        ffp = s0*(1.d0-fp)/m0+fp*s1/m1
        fact1 = ffp*m1*fi**(1.d0/a2)
        fact2 = fi**2.d0-s1+m1*ffp
        !!!
        s5 = fact1/fact2
        m5 = s5/ffp
        a5 = a2
    end if

    !!!
    !!! Recup. de h(theta)
    !!!

    rcos3t = cos3t(devsig, pref, 1.d-8)
    call srhtet(nmat, mater, rcos3t, r0c, rtheta)

    !!!
    !!! Calcul de alpha_res
    !!!

    alres = 1.d0+m3

    !!!
    !!! Calcul de sii
    !!!

    sii = norm2(devsig(1:ndt))

    !!!
    !!! Calcul de sigma_min et sigma_max
    !!!

    tiers = 1.d0/3.d0
    sigmin = i1*tiers-sii*rtheta/sqrt(6.d0)/r0c
    sigmax = i1*tiers+sqrt(2.d0/3.d0)*sii*rtheta/r0c

    if (sigmax .le. r8prem()) sigmax = r8prem()
    if (sigmin .le. r8prem()) sigmin = r8prem()

    !!!
    !!! Calcul de sigma_car
    !!!

    fact4 = sigmin*m5/sigc+s5

    if (fact4 .ge. 0.d0) then
        sigcar = sigmin+sigc*fact4**a5
        if (sigcar .le. r8prem()) sigcar = r8prem()
    else
        sigcar = r8prem()
    end if

    !!!
    !!! Calcul de sigma_tilde et alpha
    !!!

    xip = vint(1)

    if ((xi1 .le. xip) .and. (xip .lt. xi2)) then

        fact3 = 1.d0+para(1)*para(3)*para(2)**(para(1)-1.d0)
        c = sigc*(para(2)**para(1))/2.d0/sqrt(fact3)
        phi = 2.d0*atan(sqrt(fact3))-pi/2.d0
        sigtil = c/tan(phi)

        if (sigmin+sigtil .ge. r8prem()) then
            alpha = (sigmax+sigtil)/(sigmin+sigtil)
        else
            alpha = (sigmax+sigtil)/r8prem()
        end if

    else if (xip .ge. xi2) then

        sigtil = 0.d0

        if (sigmin+sigtil .ge. r8prem()) then
            alpha = (sigmax+sigtil)/(sigmin+sigtil)
        else
            alpha = (sigmax+sigtil)/r8prem()
        end if

    end if

    !!!
    !!! Calcl de sigma_tilde au pic de resistance
    !!!

    f3p = 1.d0+m1/2.d0/s1**(1.d0/2.d0)
    cp = sigc*(s1**(1.d0/2.d0))/2.d0/sqrt(f3p)
    phip = 2.d0*atan(sqrt(f3p))-pi/2.d0
    sigtilp = cp/tan(phip)

    if (sigtilp .le. r8prem()) then
        sigtilp = r8prem()
    end if

    !!!
    !!! Calcul de sin(psi)
    !!!

    if (val .eq. 0) then

        sinpsi = rho_1*(sigmax-sigcar)/(rho_2*sigmax+sigcar)

    else

        if (sigtil .gt. 0.d0) then

            if (sigmax-sigcar .gt. 0.d0) then
                spre = rho_1*(sigmax-sigcar)/(rho_2*sigmax+sigcar)
                spos = rho_1*(alpha-alres)/(rho_4*alpha+alres)
            else
                spre = 0.d0
                spos = rho_1*(alpha-alres)/(rho_4*alpha+alres)
            end if

            sinpsi = spre+(1.d0-sigtil/sigtilp)*spos

        else

            if (sigmax-sigcar .gt. 0.d0) then
                spre = rho_1*(sigmax-sigcar)/(rho_2*sigmax+sigcar)
                spos = rho_1*(alpha-alres)/(rho_4*alpha+alres)
            else
                spre = 0.d0
                spos = rho_1*(alpha-alres)/(rho_4*alpha+alres)
            end if

            sinpsi = spre+spos

        end if
    end if

    if (sinpsi .le. -1.d0) then
        sinpsi = -1.d0
    else if (sinpsi .ge. 1.d0) then
        sinpsi = 1.d0
    end if

    !!!
    !!! Calcul de d(beta')/d(sin(psi))
    !!!

    dbdsin = -6.d0*sqrt(6.d0)/(3.d0-sinpsi)**2.d0

    !!!
    !!! Termes communs
    !!!

    fact5 = rho_1*(1.d0+rho_2)/(rho_2*sigmax+sigcar)**2.d0
    ucar = sigmin*m5/sigc+s5

    if (ucar .gt. 0.d0) then
        fact6 = 1.d0+m5*a5*(ucar**(a5-1.d0))
    else
        fact6 = 1.d0
    end if

    fact7 = rho_1*alres*(1.d0+rho_4)/(rho_4*alpha+alres)**2.d0/(sigmin+sigtil)**2.d0

    !!!
    !!! Calcul de d(sigma_max)/d(s), d(sigma_min)/d(s) et d(sigma_car)/d(s)
    !!!

    call srdhds(nmat, mater, devsig, dhds, iret)

    do i = 1, ndt
        dsmids(i) = -(sii*dhds(i)+rtheta*devsig(i)/sii)/r0c/sqrt(6.d0)
        dsmads(i) = sqrt(2.d0/3.d0)*(sii*dhds(i)+rtheta*devsig(i)/sii)/r0c
        dscards(i) = dsmids(i)*fact6
    end do

    !!!
    !!! Calcul de d(sin(psi))/d(s)
    !!!

    !!! Assemblage d(sin)/d(s) pre-pic et visco.

    if (val .eq. 0) then

        do i = 1, ndt
            dsinds(i) = fact5*(sigcar*dsmads(i)-sigmax*dscards(i))
        end do

    else if (val .eq. 1) then

        if (sigmax-sigcar .gt. 0.d0) then

            do i = 1, ndt
                temp(i) = (sigmin+sigtil)*dsmads(i)-(sigmax+sigtil)*dsmids(i)
                dsinds(i) = fact5*(sigcar*dsmads(i)-sigmax*dscards(i))+ &
                            (1.d0-sigtil/sigtilp)*fact7*temp(i)
            end do

        else

            do i = 1, ndt
                temp(i) = (sigmin+sigtil)*dsmads(i)-(sigmax+sigtil)*dsmids(i)
                dsinds(i) = (1.d0-sigtil/sigtilp)*fact7*temp(i)
            end do

        end if
    end if

    !!!
    !!! Calcul de d(sin(psi))/d(i1)
    !!!

    !!!
    !!! Calcul de d(sigmax)/d(i1) et d(sigmin)/d(i1)
    !!!

    dsmidi = tiers
    dsmadi = tiers
    dscardi = dsmidi*fact6

    !!!
    !!! Assemblage de d(sin)/d(i1)
    !!!

    if (val .eq. 0) then

        dsindi = fact5*(sigcar*dsmadi-sigmax*dscardi)

    else if (val .eq. 1) then

        if (sigmax-sigcar .gt. 0.d0) then

            dsindi = fact5*(sigcar*dsmadi-sigmax*dscardi)+ &
                     (1.d0-sigtil/sigtilp)*fact7*(sigmin-sigmax)/3.d0

        else

            dsindi = (1.d0-sigtil/sigtilp)*fact7*(sigmin-sigmax)/3.d0

        end if
    end if

    !!!
    !!! Calcul de d(beta')/d(s)
    !!!

    dbetds(1:ndt) = dbdsin*dsinds(1:ndt)

    !!!
    !!! Calcul de d(beta')/d(i1)
    !!!

    dbetdi = dbdsin*dsindi

end subroutine
