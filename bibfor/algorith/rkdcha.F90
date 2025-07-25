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
subroutine rkdcha(nvi, vini, coeft, nmat, sigi, &
                  dvin)
    implicit none
!     MODELE  ELASTO-VISCOPLASTIQUE DE CHABOCHE A 26 PARAMETRES
!     A ECROUISSAGE ISOTROPE ET CINEMATIQUE
!     INTEGRATION DE  LA LOI (VISCOCHAB) PAR UNE METHODE DE RUNGE KUTTA
!
!     CETTE ROUTINE FOURNIT LA DERIVEE DE L ENSEMBLE DES VARIABLES
!     INTERNES DU MODELE
!     ----------------------------------------------------------------
    integer(kind=8) :: itens, nvi, nmat
    real(kind=8) :: coeft(nmat), vini(*), dvin(*)
    real(kind=8) :: smx(6), petin(6), petin2(6), sigi(6)
    real(kind=8) :: evi(6), a1v(6), a2v(6), csi(6), devi(6), da1v(6), da2v(6)
    real(kind=8) :: dcsi(6)
    real(kind=8) :: k0, k, n, mr, mu, m1, m2, gamma1, g20, g10, evcum, eta
    real(kind=8) :: drayvi, dqcum
    real(kind=8) :: devcum, detat, grjx2, grjx1, grjeps, grj2v, granqr, granq
    real(kind=8) :: gr, alp
    real(kind=8) :: d2, gamma2, trsig, trest, rayvi, qr0, qm, qcum, q0, gx1, gx2
    real(kind=8) :: xx, xna1v, xna2v, ze, xxn, ai, ak, tempo
    real(kind=8) :: b, c1, c2, critme, critv, d1
    data ze/0.0d0/
!
! --    COEFFICIENTS MATERIAU INELASTIQUE
!
    k0 = coeft(1)
    ak = coeft(2)
    k = coeft(4)
    n = coeft(5)
    alp = coeft(6)
    b = coeft(7)
    mr = coeft(8)
    gr = coeft(9)
    mu = coeft(10)
    qm = coeft(11)
    q0 = coeft(12)
    qr0 = coeft(13)
    eta = coeft(14)
    c1 = coeft(15)
    m1 = coeft(16)
    d1 = coeft(17)
    gx1 = coeft(18)
    g10 = coeft(19)
    c2 = coeft(20)
    m2 = coeft(21)
    d2 = coeft(22)
    gx2 = coeft(23)
    g20 = coeft(24)
    ai = coeft(25)
!
! --  VARIABLES INTERNES
!
    do itens = 1, 6
        evi(itens) = vini(itens)
        a1v(itens) = vini(itens+6)
        a2v(itens) = vini(itens+12)
        csi(itens) = vini(itens+18)
    end do
    rayvi = vini(25)
    qcum = vini(26)
    evcum = vini(27)
!
!       ----------------------------------------------------------------
    trsig = (sigi(1)+sigi(2)+sigi(3))/3.0d0
    grj2v = 0.0d0
    do itens = 1, 6
        smx(itens) = sigi(itens)-(c1*a1v(itens)+c2*a2v(itens))/1.5d0
        if (itens .le. 3) smx(itens) = smx(itens)-trsig
        grj2v = grj2v+smx(itens)**2
    end do
    grj2v = sqrt(1.5d0*grj2v)
    critv = grj2v-rayvi-k
    if (critv .le. 0.0d0) then
        drayvi = 0.0d0
        dqcum = 0.0d0
        devcum = 0.0d0
        do itens = 1, 6
            devi(itens) = 0.0d0
            da1v(itens) = 0.0d0
            da2v(itens) = 0.0d0
            dcsi(itens) = 0.0d0
        end do
    else
        tempo = critv/(k0+ak*rayvi)
        devcum = tempo**n
        if (alp .gt. 1.0d-30) devcum = devcum*exp(alp*tempo**(n+1.0d0))
        gamma1 = ai+(1.0d0-ai)*exp(-b*evcum)
        gamma2 = g20*gamma1
        gamma1 = g10*gamma1
        xna1v = 0.0d0
        xna2v = 0.0d0
        do itens = 1, 6
            petin(itens) = smx(itens)/grj2v
            devi(itens) = 1.5d0*petin(itens)*devcum
            petin(itens) = sqrt(1.5d0)*petin(itens)
            xna1v = xna1v+a1v(itens)*petin(itens)
            xna2v = xna2v+a2v(itens)*petin(itens)
        end do
!
! --    ECROUISSAGE CINEMATIQUE
!
        do itens = 1, 6
            da1v(itens) = d1*a1v(itens)+(1.0d0-d1)*xna1v*petin(itens)
            da1v(itens) = devi(itens)-gamma1*da1v(itens)*devcum
            da2v(itens) = d2*a2v(itens)+(1.0d0-d1)*xna2v*petin(itens)
            da2v(itens) = devi(itens)-gamma2*da2v(itens)*devcum
        end do
        grjx1 = 0.0d0
        grjx2 = 0.0d0
        do itens = 1, 6
            grjx1 = grjx1+a1v(itens)**2
            grjx2 = grjx2+a2v(itens)**2
        end do
        grjx1 = c1*sqrt(grjx1/1.5d0)
        if (grjx1 .gt. 1.0d-30) then
            trest = (grjx1**m1)/grjx1
            do itens = 1, 6
                da1v(itens) = da1v(itens)-gx1*trest*a1v(itens)
            end do
        end if
        grjx2 = c2*sqrt(grjx2/1.5d0)
        if (grjx2 .gt. 1.0d-30) then
            trest = (grjx2**m2)/grjx2
            do itens = 1, 6
                da2v(itens) = da2v(itens)-gx2*trest*a2v(itens)
            end do
        end if
!
! --    ECROUISSAGE ISOTROPE
!
        granq = q0+(qm-q0)*(1.0d0-exp(-2.0d0*mu*qcum))
        granqr = (qm-granq)/qm
        granqr = granq-qr0*(1.0d0-granqr*granqr)
        xx = sign(1.0d0, granqr-rayvi)
        drayvi = b*(granq-rayvi)*devcum
        drayvi = drayvi+gr*xx*(abs(granqr-rayvi))**mr
        grjeps = 0.0d0
        do itens = 1, 6
            grjeps = grjeps+(evi(itens)-csi(itens))**2
        end do
        grjeps = sqrt(grjeps*1.5d0)
        critme = grjeps/1.5d0-qcum
        if (critme .le. 0.0d0) then
            dqcum = 0.0d0
            do itens = 1, 6
                dcsi(itens) = 0.0d0
            end do
        else
            xxn = 0.0d0
            tempo = sqrt(1.5d0)/grjeps
            do itens = 1, 6
                petin2(itens) = tempo*(evi(itens)-csi(itens))
                xxn = xxn+petin(itens)*petin2(itens)
            end do
            if (xxn .le. 0.0d0) then
                dqcum = 0.0d0
                do itens = 1, 6
                    dcsi(itens) = 0.0d0
                end do
            else
                dqcum = eta*xxn*devcum
                tempo = sqrt(1.5d0)*(1.0d0-eta)*xxn*devcum
                do itens = 1, 6
                    dcsi(itens) = tempo*petin2(itens)
                end do
            end if
        end if
    end if
    detat = ze
!
! --    DERIVEES DES VARIABLES INTERNES
!
    do itens = 1, 6
        dvin(itens) = devi(itens)
        dvin(itens+6) = da1v(itens)
        dvin(itens+12) = da2v(itens)
        dvin(itens+18) = dcsi(itens)
    end do
    dvin(25) = drayvi
    dvin(26) = dqcum
    dvin(27) = devcum
    dvin(nvi) = detat
!
end subroutine
