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
subroutine te0311(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/assert.h"
#include "asterfort/coor_cyl.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvad2.h"
#include "asterfort/utmess.h"
#include "asterfort/xdeffk.h"
!
    character(len=16) :: option, nomte
!.......................................................................
!
!  CALCUL DU TAUX DE RESTITUTION D'ENERGIE ELEMENTAIRE
!  ET DES FACTEURS D'INTENSITE DES CONTRAINTES K1 K2 K3
!  BORDS ELEMENTS ISOPARAMETRIQUES 3D
!
!  OPTION : 'CALC_K_G'    (CHARGES REELLES)
!           'CALC_K_G_F'  (CHARGES FONCTIONS)
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!
! VECTEURS DIMENSIONNES POUR  NNO = 9 , NPG = 9
!.......................................................................
!
    integer(kind=8) :: ndim, nno, npg1, compt, iforf
    integer(kind=8) :: ipoids, ivf, idfde, i, j, k, kp, iforc, ind
    integer(kind=8) :: idepl, ipres, ithet, igthet, igeom, ipref, itemps, icode
    integer(kind=8) :: nnos, jgano, ino, jlsn, jlst
!
    real(kind=8) :: a1(3), a2(3), a3(3), i1(3), i2(3), epsi, dfdx(9), dfdy(9)
    real(kind=8) :: coor(18), depl(3), valpar(4)
    real(kind=8) :: a1norm, a3norm, i2norm, divt, tcla, thetx, thety, thetz
    real(kind=8) :: dth1d1, dth2d2, poids, th1, th2, lsng, lstg
    real(kind=8) :: forc, dford1(3), dford2(3), dfor(3), coorg(3)
!                                NNO      3*NNO
    real(kind=8) :: presg, forcg(3), presn(9), forcn(27)
!
    real(kind=8) :: depi
!
    real(kind=8) :: tcla1, tcla2, tcla3
    real(kind=8) :: p(3, 3), invp(3, 3), ffp(27), fkpo(3, 3)
    real(kind=8) :: phig, rg, valres(3), devres(3)
    real(kind=8) :: e, nu, mu
    real(kind=8) :: ka, coeff, coeff3
    real(kind=8) :: u1g(3), u2g(3), u3g(3)
    real(kind=8) :: g, k1, k2, k3, prsc
    integer(kind=8) :: ibalo, imate
    character(len=16) :: nomres(3)
    integer(kind=8) :: icodre(3)
!
    character(len=8) :: nompar(4)
    character(len=4) :: fami
!
    aster_logical :: fonc, l_not_zero
!.......................................................................
!
!
    depi = r8depi()
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    call jevech('PTHETAR', 'L', ithet)
    tcla = 0.d0
    tcla1 = 0.d0
    tcla2 = 0.d0
    tcla3 = 0.d0
    coeff = 0.d0
    coeff3 = 0.d0
    call jevech('PGTHETA', 'E', igthet)
!
!
! - PAS DE CALCUL DE G POUR LES ELEMENTS OU LA VALEUR DE THETA EST NULLE
!
    compt = 0
    epsi = 1.d-10
    do i = 1, nno
        thetx = zr(ithet+3*(i-1)+1-1)
        thety = zr(ithet+3*(i-1)+2-1)
        thetz = zr(ithet+3*(i-1)+3-1)
        if ((abs(thetx) .lt. epsi) .and. (abs(thety) .lt. epsi) .and. (abs(thetz) .lt. epsi)) then
            compt = compt+1
        end if
    end do
    if (compt .eq. nno) goto 999
!
! RECUPERATION CHARGE, MATER...
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
!
    call jevech('PBASLOR', 'L', ibalo)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
!
    call jevech('PMATERC', 'L', imate)
!
    if (option .eq. 'CALC_K_G_XFEM_F') then
        fonc = .true.
        call jevech('PFF2D3D', 'L', iforf)
        call jevech('PPRESSF', 'L', ipref)
        call jevech('PINSTR', 'L', itemps)
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
        valpar(4) = zr(itemps)
    else
        fonc = .false.
        call jevech('PFR2D3D', 'L', iforc)
        call jevech('PPRESSR', 'L', ipres)
    end if
!
!
!
! - SI CHARGE FONCTION RECUPERATION DES VALEURS AUX PG ET NOEUDS
!
    if (fonc) then
        do i = 1, nno
            do j = 1, 3
                valpar(j) = zr(igeom+3*(i-1)+j-1)
            end do
            call fointe('FM', zk8(ipref), 4, nompar, valpar, &
                        presn(i), icode)
            do j = 1, 3
                call fointe('FM', zk8(iforf+j-1), 4, nompar, valpar, &
                            forcn(3*(i-1)+j), icode)
            end do
        end do
    end if
!
! CALCUL DU REPERE LOCAL (A1, A2, A3)
!
    do j = 1, 3
        a1(j) = zr(igeom+3*(2-1)+j-1)-zr(igeom+3*(1-1)+j-1)
        a2(j) = zr(igeom+3*(3-1)+j-1)-zr(igeom+3*(1-1)+j-1)
    end do
!
    a3(1) = a1(2)*a2(3)-a1(3)*a2(2)
    a3(2) = a1(3)*a2(1)-a1(1)*a2(3)
    a3(3) = a1(1)*a2(2)-a1(2)*a2(1)
!
! CALCUL DU REPERE LOCAL ORTHONORME (I1, I2, A3)
!
    i2(1) = a3(2)*a1(3)-a3(3)*a1(2)
    i2(2) = a3(3)*a1(1)-a3(1)*a1(3)
    i2(3) = a3(1)*a1(2)-a3(2)*a1(1)
!
    a1norm = sqrt(a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3))
    i2norm = sqrt(i2(1)*i2(1)+i2(2)*i2(2)+i2(3)*i2(3))
    a3norm = sqrt(a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3))
    do i = 1, 3
        i1(i) = a1(i)/a1norm
        i2(i) = i2(i)/i2norm
        a3(i) = a3(i)/a3norm
    end do
!
    do i = 1, nno
        coor(2*i-1) = 0.d0
        coor(2*i) = 0.d0
        do j = 1, 3
            coor(2*i-1) = coor(2*i-1)+(zr(igeom+3*(i-1)+j-1)-zr(igeom+j-1))*i1(j)
            coor(2*i) = coor(2*i)+(zr(igeom+3*(i-1)+j-1)-zr(igeom+j-1))*i2(j)
        end do
    end do
!
! --- BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg1
!
        k = (kp-1)*nno
!
        do j = 1, 3
            depl(j) = 0.d0
            dford1(j) = 0.d0
            dford2(j) = 0.d0
            dfor(j) = 0.d0
            coorg(j) = 0.d0
        end do
!
!
        th1 = 0.d0
        th2 = 0.d0
        dth1d1 = 0.d0
        dth2d2 = 0.d0
        lsng = 0.d0
        lstg = 0.d0
!
        do i = 1, nno
!
            lsng = lsng+zr(jlsn-1+i)*zr(ivf+k+i-1)
            lstg = lstg+zr(jlst-1+i)*zr(ivf+k+i-1)
!
            do j = 1, 3
                coorg(j) = coorg(j)+zr(ivf+k+i-1)*zr(igeom+3*(i-1)+j-1)
            end do
        end do
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    CALCUL DES COOR. CYL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        p(:, :) = 0.d0
        invp(:, :) = 0.d0
        do ino = 1, nno
            ffp(ino) = zr(ivf-1+nno*(kp-1)+ino)
        end do
        call coor_cyl(3, nno, zr(ibalo), zr(igeom), ffp, &
                      p, invp, rg, phig, l_not_zero)
!
!       ON A PAS PU CALCULER LES DERIVEES DES FONCTIONS SINGULIERES
!       CAR ON SE TROUVE SUR LE FOND DE FISSURE
        ASSERT(l_not_zero)
!
!
        if ((abs(lsng) .lt. 1.0d-8) .and. (lstg .lt. 0.0d0)) then
!
!         ON DETERMINE SI ON EST SUR LA LEVRE X2 > 0 OU
!         SUR LA LEVRE X2 < 0
!
!         PRODUIT SCALAIRE E2 (AXE X2) * A3 (NORMALE DE L'ELEMENT)
            prsc = 0.0d0
            do i = 1, 3
                prsc = prsc+p(i, 2)*a3(i)
            end do
            if (prsc .gt. 0.0d0) then
!            ON EST SUR LA LEVRE X2 < 0
                phig = -1.0d0*abs(phig)
            else
!            ON EST SUR LA LEVRE X2 > 0
                phig = abs(phig)
            end if
        end if
!
!       RECUPERATION DES DONNEES MATERIAU
!
        nomres(1) = 'E'
        nomres(2) = 'NU'
        nomres(3) = 'ALPHA'
!
        call rcvad2(fami, kp, 1, '+', zi(imate), &
                    'ELAS', 3, nomres, valres, devres, &
                    icodre)
        if ((icodre(1) .ne. 0) .or. (icodre(2) .ne. 0)) then
            call utmess('F', 'RUPTURE1_25')
        end if
        if (icodre(3) .ne. 0) then
            valres(3) = 0.d0
            devres(3) = 0.d0
        end if
!
        e = valres(1)
        nu = valres(2)
!
        mu = e/(2.d0*(1.d0+nu))
!
!       COEFFICIENTS DE CALCUL
        ka = 3.d0-4.d0*nu
!
        coeff = e/(1.d0-nu*nu)
        coeff3 = 2.d0*mu
!
!-----------------------------------------------------------------------
!       DEFINITION DU CHAMP SINGULIER AUXILIAIRE U1
!-----------------------------------------------------------------------
!
! --------- champs singuliers
        call xdeffk(ka, mu, rg, phig, 3, &
                    fkpo)
!
        u1g(:) = 0.
        u2g(:) = 0.
        u3g(:) = 0.
        do i = 1, 3
            do ind = 1, 3
                u1g(i) = u1g(i)+p(i, ind)*fkpo(1, ind)
                u2g(i) = u2g(i)+p(i, ind)*fkpo(2, ind)
                u3g(i) = u3g(i)+p(i, ind)*fkpo(3, ind)
            end do
        end do
!
!
        call dfdm2d(nno, kp, ipoids, idfde, coor, &
                    poids, dfdx, dfdy)
!
        if (fonc) then
            do j = 1, 3
                valpar(j) = coorg(j)
            end do
            call fointe('FM', zk8(ipref), 4, nompar, valpar, &
                        presg, icode)
            do j = 1, 3
                call fointe('FM', zk8(iforf+j-1), 4, nompar, valpar, &
                            forcg(j), icode)
            end do
!
            do i = 1, nno
                do j = 1, 3
                    dford1(j) = dford1(j)+(forcn(3*(i-1)+j)-presn(i)*a3(j))*dfdx(i)
                    dford2(j) = dford2(j)+(forcn(3*(i-1)+j)-presn(i)*a3(j))*dfdy(i)
                end do
            end do
        else
            presg = 0.d0
            forcg(1) = 0.d0
            forcg(2) = 0.d0
            forcg(3) = 0.d0
            do i = 1, nno
                presg = presg+zr(ipres+i-1)*zr(ivf+k+i-1)
                do j = 1, 3
                    forcg(j) = forcg(j)+zr(iforc+3*(i-1)+j-1)*zr(ivf+k+i-1)
                end do
            end do
        end if
!
        do i = 1, nno
            do j = 1, 3
                depl(j) = depl(j)+zr(ivf+k+i-1)*zr(idepl+3*(i-1)+j-1)
                th1 = th1+zr(ivf+k+i-1)*zr(ithet+3*(i-1)+j-1)*i1(j)
                th2 = th2+zr(ivf+k+i-1)*zr(ithet+3*(i-1)+j-1)*i2(j)
                dth1d1 = dth1d1+zr(ithet+3*(i-1)+j-1)*i1(j)*dfdx(i)
                dth2d2 = dth2d2+zr(ithet+3*(i-1)+j-1)*i2(j)*dfdy(i)
            end do
        end do
!
        do j = 1, 3
            dfor(j) = dfor(j)+dford1(j)*th1+dford2(j)*th2
        end do
!
        divt = dth1d1+dth2d2
!
        do j = 1, 3
            forc = forcg(j)-presg*a3(j)
            tcla = tcla+poids*(forc*divt+dfor(j))*depl(j)
            tcla1 = tcla1+0.5d0*poids*(forc*divt+dfor(j))*u1g(j)
            tcla2 = tcla2+0.5d0*poids*(forc*divt+dfor(j))*u2g(j)
            tcla3 = tcla3+0.5d0*poids*(forc*divt+dfor(j))*u3g(j)
!
        end do
!
    end do
999 continue
!
    g = tcla
    k1 = tcla1*coeff
    k2 = tcla2*coeff
    k3 = tcla3*coeff3
!
    zr(igthet) = g
    zr(igthet+1) = tcla1*sqrt(coeff)
    zr(igthet+2) = tcla2*sqrt(coeff)
    zr(igthet+3) = tcla3*sqrt(coeff3)
    zr(igthet+4) = k1
    zr(igthet+5) = k2
    zr(igthet+6) = k3
!
end subroutine
