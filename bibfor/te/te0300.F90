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
subroutine te0300(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterc/r8prem.h"
#include "asterfort/coor_cyl.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/provec.h"
#include "asterfort/rcvad2.h"
#include "asterfort/utmess.h"
#include "asterfort/xdeffk.h"
!
    character(len=16) :: option, nomte
!.......................................................................
!
!      CALCUL DES COEFFICIENTS DE CONTRAINTES K1 ET K2
!      BORDS ELEMENTS ISOPARAMETRIQUES 2D AVEC CHARGEMENT DE BORD
!      PRESSION-CISAILLEMENT ET FORCE REPARTIE
!
!      OPTION : 'CALC_K_G_XFEM'  (CHARGES REELLES)
!               'CALC_K_G_XFEM_F' (CHARGES FONCTIONS)
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!
! VECTEURS DIMENSIONNES POUR  NNO = 3 , NPG = 4
!.......................................................................
!
    integer(kind=8) :: nno, nnos, jgano, ndim, npg, kp, compt, i, j, k
    integer(kind=8) :: idepl, ific, ifond, iforc, imate, ipres, ithet
    integer(kind=8) :: ipoids, ivf, idfdk, igeom, itemps
    integer(kind=8) :: iforf, ipref, icode, ino, ind
!
    real(kind=8) :: depi, eps, valres(3), devres(3), valpar(3)
    real(kind=8) :: tcla, tcla1, tcla2, u1s(2), u2s(2), ux, uy
    real(kind=8) :: vf, dfde, dxde, dyde, dsde, poids, dthxde, dthyde, thx, thy
    real(kind=8) :: g, k1, k2, fx, fy, pres, cisa, divthe
    real(kind=8) :: xg, yg, rpol, phi
    real(kind=8) :: cpk, dpk, ck, coefk, dcoefk, ccoefk, cform
    real(kind=8) :: the, dfxde, dfyde, presno, cisano, fxno, fyno
!                                            2*NNO     2*NNO
    real(kind=8) :: presg(2), forcg(2), presn(6), forcn(6)
    real(kind=8) :: basloc(9*6), p(3, 3), invp(3, 3), e1(3), e2(3), e3(3)
    real(kind=8) :: fkpo(3, 3), ffp(9), mu, pt_ree(2), pt_loc(2)
    real(kind=8) :: xno1, xno2, yno1, yno2, d1, d2
!
    integer(kind=8) :: icodre(3)
    character(len=4) :: fami
    character(len=8) :: nompar(3), elrefe
    character(len=16) :: nomres(3)
!
    aster_logical :: fonc, axi, l_not_zero
!.......................................................................
!
    call elref1(elrefe)
    eps = r8prem()
    depi = r8depi()
    axi = .false.
    if (lteatt('AXIS', 'OUI')) axi = .true.
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
    call jevech('PTHETAR', 'L', ithet)
    tcla = 0.d0
    tcla1 = 0.d0
    tcla2 = 0.d0
    call jevech('PGTHETA', 'E', ific)
!
! PAS DE CALCUL DE G POUR LES ELEMENTS OU LA VALEUR DE THETA EST NULLE
!
    compt = 0
    do i = 1, nno
        thx = zr(ithet+2*(i-1))
        thy = zr(ithet+2*(i-1)+1)
        if ((abs(thx) .lt. eps) .and. (abs(thy) .lt. eps)) then
            compt = compt+1
        end if
    end do
    if (compt .eq. nno) goto 110
!
! RECUPERATION CHARGE, MATER...
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PDEPLAR', 'L', idepl)
    call jevech('PFISSR', 'L', ifond)
    if ((option .eq. 'CALC_K_G_XFEM_F') .or. (option .eq. 'G_MODA_F')) then
        fonc = .true.
        call jevech('PFF1D2D', 'L', iforf)
        call jevech('PPRESSF', 'L', ipref)
        call jevech('PINSTR', 'L', itemps)
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'INST'
        valpar(3) = zr(itemps)
    else
        fonc = .false.
        call jevech('PFR1D2D', 'L', iforc)
        call jevech('PPRESSR', 'L', ipres)
    end if
!
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
!
!
! - SI CHARGE FONCTION RECUPERATION DES VALEURS AUX PG ET NOEUDS
!
    if (fonc) then
        do i = 1, nno
            do j = 1, 2
                valpar(j) = zr(igeom+2*(i-1)+j-1)
            end do
            do j = 1, 2
                call fointe('FM', zk8(ipref+j-1), 3, nompar, valpar, &
                            presn(2*(i-1)+j), icode)
                call fointe('FM', zk8(iforf+j-1), 3, nompar, valpar, &
                            forcn(2*(i-1)+j), icode)
            end do
        end do
    end if
!
! --- BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg
        k = (kp-1)*nno
        xg = 0.d0
        yg = 0.d0
        dxde = 0.d0
        dyde = 0.d0
        ux = 0.d0
        uy = 0.d0
        thx = 0.d0
        dfxde = 0.d0
        dfyde = 0.d0
        dthxde = 0.d0
        dthyde = 0.d0
        divthe = 0.d0
!
!
        do i = 1, nno
            vf = zr(ivf+k+i-1)
            dfde = zr(idfdk+k+i-1)
            xg = xg+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
            yg = yg+zr(igeom+2*(i-1)+1)*zr(ivf+k+i-1)
            dxde = dxde+dfde*zr(igeom+2*(i-1))
            dyde = dyde+dfde*zr(igeom+2*(i-1)+1)
            ux = ux+vf*zr(idepl+2*(i-1))
            uy = uy+vf*zr(idepl+2*(i-1)+1)
            thx = thx+vf*zr(ithet+2*(i-1))
            thy = thy+vf*zr(ithet+2*(i-1)+1)
            dthxde = dthxde+dfde*zr(ithet+2*(i-1))
            dthyde = dthyde+dfde*zr(ithet+2*(i-1)+1)
        end do
!
        if (fonc) then
            valpar(1) = xg
            valpar(2) = yg
            do j = 1, 2
                call fointe('FM', zk8(ipref+j-1), 3, nompar, valpar, &
                            presg(j), icode)
                call fointe('FM', zk8(iforf+j-1), 3, nompar, valpar, &
                            forcg(j), icode)
            end do
        else
            presg(1) = 0.d0
            presg(2) = 0.d0
            forcg(1) = 0.d0
            forcg(2) = 0.d0
            do i = 1, nno
                do j = 1, 2
                    presg(j) = presg(j)+zr(ipres+2*(i-1)+j-1)*zr(ivf+k+i-1)
                    forcg(j) = forcg(j)+zr(iforc+2*(i-1)+j-1)*zr(ivf+k+i-1)
                end do
            end do
        end if
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
        dpk = 3.d0-4.d0*valres(2)
        cpk = (3.d0-valres(2))/(1.d0+valres(2))
        cform = (1.d0+valres(2))/(sqrt(depi)*valres(1))
        dcoefk = valres(1)/(1.d0-valres(2)*valres(2))
        ccoefk = valres(1)
        mu = valres(1)/(2.d0*(1.d0+valres(2)))
!
        if (axi) then
            ck = dpk
            coefk = dcoefk
        else
            ck = cpk
            coefk = ccoefk
        end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    CALCUL DES COOR. CYL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        p(:, :) = 0.d0
        invp(:, :) = 0.d0
        do ino = 1, nno
            ffp(ino) = zr(ivf-1+nno*(kp-1)+ino)
            basloc((6*(ino-1)+1):(6*(ino-1)+6)) = zr((ifond-1+1):(ifond-1+6))
        end do
        call coor_cyl(2, nno, basloc, zr(igeom), ffp, &
                      p, invp, rpol, phi, l_not_zero)
! BRICOLAGE POUR CALCULER LE SIGNE DE K2 QUAND NDIM=2
        e1(:) = 0.d0
        e1(1:2) = p(1:2, 1)
        e2(:) = 0.d0
        e2(1:2) = p(1:2, 2)
        call provec(e1, e2, e3)
        p(3, 3) = e3(3)
        invp(3, 3) = e3(3)
! BRICOLAGE POUR DETERMINER LE SIGNE DE LA PHI
        pt_ree = [xg, yg]-zr((ifond-1+1):(ifond-1+2))
!
        pt_loc(:) = 0.d0
        do i = 1, 2
            do ind = 1, 2
                pt_loc(i) = pt_loc(i)+invp(i, ind)*pt_ree(ind)
            end do
        end do
!
        if ((abs(pt_loc(2)) .lt. 1.0d-8) .and. (pt_loc(1) .lt. 0.0d0)) then
!
! ON DETERMINE SI ON EST SUR LA LEVRE X2 > 0 OU
! SUR LA LEVRE X2 < 0
!
            xno1 = zr(igeom)
            yno1 = zr(igeom+1)
            xno2 = zr(igeom+2)
            yno2 = zr(igeom+3)
            d1 = ( &
                 ( &
                 xno1-zr(ifond-1+1))*(xno1-zr(ifond-1+1)))+((yno1-zr(ifond-1+2))*(yno1-zr(ifond-1&
                 &+2) &
                 ) &
                 )
            d2 = ( &
                 ( &
                 xno2-zr(ifond-1+1))*(xno2-zr(ifond-1+1)))+((yno2-zr(ifond-1+2))*(yno2-zr(ifond-1&
                 &+2) &
                 ) &
                 )
            if (d2 .gt. d1) then
                phi = -1.0d0*phi
            else
                phi = abs(phi)
            end if
        end if
!
        if (axi .and. (xg .lt. r8prem())) then
            call utmess('F', 'RUPTURE0_56')
        end if
!
! --------- champs singuliers
        call xdeffk(ck, mu, rpol, phi, 2, &
                    fkpo(1:2, 1:2))
!
        u1s(:) = 0.d0
        u2s(:) = 0.d0
        do i = 1, 2
            do ind = 1, 2
                u1s(i) = u1s(i)+p(i, ind)*fkpo(1, ind)
                u2s(i) = u2s(i)+p(i, ind)*fkpo(2, ind)
            end do
        end do
!        print*,' *** KOR ***'
!        print*,'  - rg, phig',rpol, phi
!        print*,'  - ori',zr((ifond-1+1):(ifond-1+2))
!        print*,'  - e1=',e1
!        print*,'  - e2=',e2
!        print*,'  - u1s',u1s
!        print*,'  - u2s',u2s
!        print*,' ***********'
!
        dsde = sqrt(dxde**2+dyde**2)
!
        pres = presg(1)
        cisa = presg(2)
        fx = forcg(1)-(dyde*pres-dxde*cisa)/dsde
        fy = forcg(2)+(dxde*pres+dyde*cisa)/dsde
!
        if (fonc) then
            do i = 1, nno
                dfde = zr(idfdk+k+i-1)
                presno = presn(2*(i-1)+1)
                cisano = presn(2*(i-1)+2)
                fxno = forcn(2*(i-1)+1)-(dyde*presno-dxde*cisano)/dsde
                fyno = forcn(2*(i-1)+2)+(dxde*presno+dyde*cisano)/dsde
                dfxde = dfxde+dfde*fxno
                dfyde = dfyde+dfde*fyno
            end do
        end if
!
        poids = zr(ipoids+kp-1)
        if (axi) poids = poids*xg
        the = (thx*dxde+thy*dyde)/dsde
        divthe = (dthxde*dxde+dthyde*dyde)/dsde
        if (axi) divthe = divthe+(thx*dsde/xg)
!
        tcla1 = tcla1+poids*((divthe*fx+dfxde*the)*u1s(1)+(divthe*fy+dfyde*the)*u1s(2))
        tcla2 = tcla2+poids*((divthe*fx+dfxde*the)*u2s(1)+(divthe*fy+dfyde*the)*u2s(2))
        tcla = tcla+poids*((divthe*fx+dfxde*the)*ux+(divthe*fy+dfyde*the)*uy)
!
    end do
!
    g = tcla
    k1 = tcla1*coefk/2.d0
    k2 = tcla2*coefk/2.d0
!    if (e3(3) .lt. 0) k2=-k2
!
    zr(ific) = g
    zr(ific+1) = k1/sqrt(coefk)
    zr(ific+2) = k2/sqrt(coefk)
    zr(ific+3) = k1
    zr(ific+4) = k2
!
110 continue
end subroutine
