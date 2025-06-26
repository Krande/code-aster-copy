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

subroutine te0498(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/pronor.h"
#include "asterfort/rcvalb.h"
#include "asterfort/trigom.h"
#include "asterfort/utmess.h"
!
!
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
! Option: ONDE_PLAN
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom, i, j
    integer(kind=8) :: ndim, nno, ipg, npg1, ino, jno
    integer(kind=8) :: idec, jdec, kdec, ldec, ires, imate
    integer(kind=8) :: ii, mater, jinst, indic1, indic2
    integer(kind=8) :: ionde, iondc, ier, nnos, jgano
    real(kind=8) :: jac, nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9)
    real(kind=8) :: valres(5), e, nu, lambda, mu, cp, cs, rho, typer
    real(kind=8) :: taux, tauy, tauz, dirx, diry, dirz
    real(kind=8) :: norm, tanx, tany, norx, nory, norz
    real(kind=8) :: taondx, taondy, taondz, sign_tan
    real(kind=8) :: nux, nuy, nuz, scal, coedir, coef_amor, coef_dse
    real(kind=8) :: nortan, cele, trace, cele2
    real(kind=8) :: param0, param, h, h2, instd, ris, rip, l0, usl0
    real(kind=8) :: sigma(3, 3), epsi(3, 3), grad(3, 3), valfon
    real(kind=8) :: xgg(9), ygg(9), zgg(9), vondn(3), vondt(3), uondn(3), uondt(3)
    real(kind=8) :: a2, b2, sina, cosa, cosg, sing, sinb2, cosb2, rc1c2, ra12, ra13, kr, nr
    real(kind=8) :: xsv, zsv, ysv, dist1, dist2, instd1, instd2, x0, y0, z0, x1, y1, z1
    real(kind=8) :: valfon1, valfon2, param1, param2
    integer(kind=8) :: icodre(5), ndim2
    character(len=2) :: type
    character(len=8) :: fami, poum
    character(len=8) :: nompar(3)
    character(len=8) :: lpar2(2)
    real(kind=8) :: vpar2(2)
    real(kind=8) :: xyzgau(3)
    integer(kind=8) :: idecpg, idecno
    character(len=16), parameter :: nomres(5) = (/'E        ', 'NU       ', &
                                                  'RHO      ', &
                                                  'COEF_AMOR', 'LONG_CARA'/)
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'ONDE_PLAN')
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PONDPLA', 'L', ionde)
    call jevech('PONDPLR', 'L', iondc)
    call jevech('PINSTR', 'L', jinst)
    call jevech('PVECTUR', 'E', ires)
!
    if (zk8(ionde) (1:7) .eq. '&FOZERO') goto 99
!
!     --- INITIALISATION
!
    sigma = 0.d0
!
    mater = zi(imate)
    fami = 'RIGI'
    poum = '+'
    ASSERT(ndim .ne. 3)
    ndim2 = ndim+1
!
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'Z'
!
!     --- CARACTERISTIQUES DE L'ONDE PLANE
!
    dirx = zr(iondc)
    diry = zr(iondc+1)
    dirz = zr(iondc+2)
    typer = zr(iondc+3)
    x0 = zr(iondc+4)
    y0 = zr(iondc+5)
    z0 = zr(iondc+6)
    x1 = zr(iondc+7)
    y1 = zr(iondc+8)
    z1 = zr(iondc+9)
!
    if (typer .eq. 0.d0) type = 'P'
    if (typer .eq. 1.d0) type = 'SV'
    if (typer .eq. 2.d0) type = 'SH'
!
!     --- CALCUL DU VECTEUR DIRECTEUR UNITAIRE DE L'ONDE PLANE
!
    norm = sqrt(dirx**2.d0+diry**2.d0+dirz**2.d0)
    dirx = dirx/norm
    diry = diry/norm
    dirz = dirz/norm

    if (x0 .ne. r8vide()) then
        h = x0*dirx+y0*diry+z0*dirz
        if (x1 .ne. r8vide()) then
            h2 = x1*dirx+y1*diry+z1*dirz
        else
            h2 = r8vide()
        end if
    else
        h = r8vide()
    end if
!
    if (abs(dirz-1.d0) .le. r8prem()) then

        cosa = dirz
        sina = 0.d0
        cosg = 1.d0
        sing = 0.d0

    else

        cosa = dirz

        if (abs(diry) .le. r8prem()) then

            sina = dirx
            cosg = 1.d0
            sing = 0.d0

        else if (abs(dirx) .le. r8prem()) then

            sina = diry
            cosg = 0.d0
            sing = 1.d0

        else

            if (dirx .gt. 0.d0) then

                sina = sin(acos(cosa))

            else

                sina = -sin(acos(cosa))

            end if

            cosg = dirx/sina
            sing = diry/sina

        end if

    end if

! Cette condition est nécéssaire pour DEFI_SOL_EQUI
! qui opère avec l'axe vertical identifié comme l'axe Y

    if (dirz .gt. 0) then

        coef_dse = 1.d0

    else

        coef_dse = 0.d0

    end if

    !     CALCUL DU REPERE ASSOCIE A L'ONDE

    tanx = sing
    tany = -cosg

    if (cosg .gt. 0) then

        sign_tan = tany/abs(tany)

    else

        sign_tan = -tanx/abs(tanx)

    end if
!
    nortan = sqrt(tanx**2.d0+tany**2.d0)

    tanx = tanx/nortan
    tany = tany/nortan
!
    norx = -tany*dirz
    nory = tanx*dirz
    norz = -tanx*diry+tany*dirx
!
!
!     --- CALCUL DES PRODUITS VECTORIELS OMI X OMJ ---
!
    do ino = 1, nno
        i = igeom+3*(ino-1)-1
        do jno = 1, nno
            j = igeom+3*(jno-1)-1
            sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
            sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
            sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
        end do
    end do
    do ipg = 1, npg1
        xgg(ipg) = 0.d0
        ygg(ipg) = 0.d0
        zgg(ipg) = 0.d0
    end do
!
!    write(6,*) 'npg=',npg1,'nno=',nno
    do ipg = 1, npg1
        ldec = (ipg-1)*nno
        do i = 1, nno
            ii = 3*i-2
            xgg(ipg) = xgg(ipg)+zr(igeom+ii-1)*zr(ivf+ldec+i-1)
            ygg(ipg) = ygg(ipg)+zr(igeom+ii)*zr(ivf+ldec+i-1)
            zgg(ipg) = zgg(ipg)+zr(igeom+ii+1)*zr(ivf+ldec+i-1)
        end do
!       write(6,*) 'kp=',ipg,'xgg=',xgg(ipg),'ygg=',ygg(ipg),'zgg=',zgg(ipg)
    end do
!
!     --- BOUCLE SUR LES POINTS DE GAUSS ---
!
    do ipg = 1, npg1
!
! -     Get material properties
!
        idecpg = nno*(ipg-1)-1
! ----- Coordinates for current Gauss point
        xyzgau(:) = 0.d0
        do i = 1, nno
            idecno = ndim2*(i-1)-1
            do j = 1, ndim2
                xyzgau(j) = xyzgau(j)+zr(ivf+i+idecpg)*zr(igeom+j+idecno)
            end do
        end do
!
        call rcvalb(fami, ipg, 1, poum, mater, &
                    ' ', 'ELAS', 3, nompar, xyzgau, &
                    4, nomres, valres, icodre, 1)
!       appel LONG_CARA en iarret = 0
        call rcvalb(fami, ipg, 1, poum, mater, &
                    ' ', 'ELAS', 3, nompar, xyzgau, &
                    1, nomres(5), valres(5), icodre(5), 0)
!
        e = valres(1)
        nu = valres(2)
        rho = valres(3)
        coef_amor = valres(4)
        !
        usl0 = 0.d0
        if (icodre(5) .eq. 0) then
            l0 = valres(5)
            usl0 = 1.d0/l0
        end if
        !
        lambda = e*nu/(1.d0+nu)/(1.d0-2.d0*nu)
        mu = e/2.d0/(1.d0+nu)
        !
        cp = sqrt((lambda+2.d0*mu)/rho)
        cs = sqrt(mu/rho)
        rip = (lambda+2.d0*mu)*usl0
        ris = mu*usl0
!
        if (type .eq. 'P') then
            cele = cp
            cele2 = cs
        else
            cele = cs
            cele2 = cp
        end if

! Calcul du rapport des célérités des ondes de cisaillement (c2) et de compression (c1).
        rc1c2 = sqrt(2.d0+2.d0*nu/(1.d0-2.d0*nu))

! Calcul de l'angle réflexion de l'onde SV réfléchie.

        a2 = asin(sina)

        if (h .ne. r8vide()) then
            if (h2 .ne. r8vide()) then

                if (type .eq. 'P') then
                    b2 = asin(sina/rc1c2)
                    cosb2 = cos(b2)
                    sinb2 = sin(b2)

! Coefficients de réflexions (Calcul intégré au script).
                    kr = 1.d0/rc1c2
                    nr = (kr)**2*sin(2.d0*b2)*sin(2.d0*a2)+cos(2.d0*b2)**2
                    ra12 = (kr**2*sin(2.d0*a2)*sin(2.d0*b2)-cos(2.d0*b2)**2)/nr
                    ra13 = 2.d0*kr*sin(2.d0*a2)*cos(2.d0*b2)/nr
                else if (type .eq. 'SV') then

                    if ((abs(sina*rc1c2)-1.d0) .gt. r8prem()) then
                        call utmess('F', 'ALGORITH2_82')
                    end if

                    b2 = trigom('ASIN', sina*rc1c2*coef_dse)
                    cosb2 = cos(b2)
                    sinb2 = sin(b2)

! Coefficients de réflexions (Calcul intégré au script).
                    kr = 1.d0/rc1c2
                    nr = (kr)**2*sin(2.d0*b2)*sin(2.d0*a2)+cos(2.d0*a2)**2
                    ra12 = (kr**2*sin(2.d0*a2)*sin(2.d0*b2)-cos(2.d0*a2)**2)/nr
                    ra13 = -2.d0*kr*sin(2.d0*a2)*cos(2.d0*a2)/nr
                else if (type .eq. 'SH') then
                    b2 = a2
                    cosb2 = cos(b2)
                    sinb2 = sin(b2)

! Coefficients de réflexions (Calcul intégré au script).
                    ra12 = 1.0d0
                    ra13 = 0.d0
                end if

! Calcul des bons paramètres dist à insérer dans le calcul.

                if (abs(cosa) .gt. 0.d0) then
                    dist1 = x0*sina*cosg+y0*sina*sing+(2.d0*z1-z0)*(-cosa)
                    if (type .eq. 'P') then
                        zsv = z1+cosb2*(z1-z0)/rc1c2/cosa
                        xsv = x0+((sina/cosa)*(z1-z0)-sinb2*(z1-z0)/rc1c2/cosa)*cosg
                        ysv = y0+((sina/cosa)*(z1-z0)-sinb2*(z1-z0)/rc1c2/cosa)*sing
                    else
                        zsv = z1+cosb2*(z1-z0)*rc1c2/cosa
                        xsv = x0+((sina/cosa)*(z1-z0)-sinb2*(z1-z0)*rc1c2/cosa)*cosg
                        ysv = y0+((sina/cosa)*(z1-z0)-sinb2*(z1-z0)*rc1c2/cosa)*sing
                    end if
                    dist2 = xsv*sinb2*cosg+ysv*sinb2*sing+zsv*(-cosb2)
                end if
            end if
        end if
!
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
!        --- CALCUL DU CHARGEMENT PAR ONDE PLANE
!KH          ON SUPPOSE QU'ON RECUPERE UNE VITESSE
        param0 = dirx*xgg(ipg)+diry*ygg(ipg)+dirz*zgg(ipg)
        valfon1 = 0.d0
        valfon2 = 0.d0
        if (h .ne. r8vide()) then
            param = param0-h
            instd = zr(jinst)-param/cele
            if (instd .lt. 0.d0) then
                valfon = 0.d0
            else
                call fointe('F ', zk8(ionde), 1, 'INST', [instd], valfon, ier)
            end if
            if (h2 .ne. r8vide()) then
                if (abs(cosa) .gt. 0.d0) then
                    param1 = dirx*xgg(ipg)+diry*ygg(ipg)-dirz*zgg(ipg)
                    param = param1-dist1
                    instd1 = zr(jinst)-param/cele
                    if (instd1 .lt. 0.d0) then
                        valfon1 = 0.d0
                    else
                        call fointe('F ', zk8(ionde), 1, 'INST', [instd1], valfon1, ier)
                        valfon1 = valfon1*ra12
                    end if
                    param2 = sinb2*cosg*xgg(ipg)+sinb2*sing*ygg(ipg)-cosb2*zgg(ipg)
                    param = param2-dist2
                    instd2 = zr(jinst)-param/cele2
                    if (instd2 .lt. 0.d0) then
                        valfon2 = 0.d0
                    else
                        call fointe('F ', zk8(ionde), 1, 'INST', [instd2], valfon2, ier)
                        valfon2 = valfon2*ra13
                    end if
                else
                    param1 = 2.0d0*(h2-h)-param
                    instd1 = zr(jinst)-param1/cele
                    if (instd1 .lt. 0.d0) then
                        valfon1 = 0.d0
                    else
                        call fointe('F ', zk8(ionde), 1, 'INST', [instd1], valfon1, ier)
                    end if
                end if
            end if
        else
            lpar2(1) = 'X'
            lpar2(2) = 'INST'
            vpar2(1) = 1.0*param0
            vpar2(2) = zr(jinst)
            call fointe('F ', zk8(ionde), 2, lpar2, vpar2, valfon, ier)
            if (type .ne. 'P') then
                valfon = -valfon
            end if
        end if

        valfon = -valfon/cele
        valfon1 = -valfon1/cele
        valfon2 = -valfon2/cele2
!
!        CALCUL DES CONTRAINTES ASSOCIEES A L'ONDE PLANE
!        CALCUL DU GRADIENT DU DEPLACEMENT
        if (type .eq. 'P') then
!
            if (abs(cosa) .gt. 0.d0) then
                grad(1, 1) = dirx*valfon*dirx
                grad(1, 2) = diry*valfon*dirx
                grad(1, 3) = dirz*valfon*dirx

                grad(2, 1) = dirx*valfon*diry
                grad(2, 2) = diry*valfon*diry
                grad(2, 3) = dirz*valfon*diry

                grad(3, 1) = dirx*valfon*dirz
                grad(3, 2) = diry*valfon*dirz
                grad(3, 3) = dirz*valfon*dirz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        grad(1, 1) = grad(1, 1)+dirx*valfon1*dirx
                        grad(1, 1) = grad(1, 1)+sinb2*cosg*valfon2*cosb2*cosg
                        grad(1, 2) = grad(1, 2)+diry*valfon1*dirx
                        grad(1, 2) = grad(1, 2)+sinb2*sing*valfon2*cosb2*cosg
                        grad(1, 3) = grad(1, 3)-dirz*valfon1*dirx
                        grad(1, 3) = grad(1, 3)-cosb2*valfon2*cosb2*cosg

                        grad(2, 1) = grad(2, 1)+dirx*valfon1*diry
                        grad(2, 1) = grad(2, 1)+sinb2*cosg*valfon2*cosb2*sing
                        grad(2, 2) = grad(2, 2)+diry*valfon1*diry
                        grad(2, 2) = grad(2, 2)+sinb2*sing*valfon2*cosb2*sing
                        grad(2, 3) = grad(2, 3)-dirz*valfon1*diry
                        grad(2, 3) = grad(2, 3)-cosb2*valfon2*cosb2*sing

                        grad(3, 1) = grad(3, 1)-dirx*valfon1*dirz
                        grad(3, 1) = grad(3, 1)+sinb2*cosg*valfon2*sinb2
                        grad(3, 2) = grad(3, 2)-diry*valfon1*dirz
                        grad(3, 2) = grad(3, 2)+sinb2*sing*valfon2*sinb2
                        grad(3, 3) = grad(3, 3)+dirz*valfon1*dirz
                        grad(3, 3) = grad(3, 3)-cosb2*valfon2*sinb2
                    end if
                end if

            else
                grad(1, 1) = dirx*(valfon-valfon1)*dirx
                grad(1, 2) = diry*(valfon-valfon1)*dirx
                grad(1, 3) = dirz*(valfon-valfon1)*dirx
!
                grad(2, 1) = dirx*(valfon-valfon1)*diry
                grad(2, 2) = diry*(valfon-valfon1)*diry
                grad(2, 3) = dirz*(valfon-valfon1)*diry
!
                grad(3, 1) = dirx*(valfon-valfon1)*dirz
                grad(3, 2) = diry*(valfon-valfon1)*dirz
                grad(3, 3) = dirz*(valfon-valfon1)*dirz
            end if
!
        else if (type .eq. 'SV') then
!
            if (abs(cosa) .gt. 0.d0) then
                grad(1, 1) = dirx*valfon*norx
                grad(1, 2) = diry*valfon*norx
                grad(1, 3) = dirz*valfon*norx

                grad(2, 1) = dirx*valfon*nory
                grad(2, 2) = diry*valfon*nory
                grad(2, 3) = dirz*valfon*nory

                grad(3, 1) = dirx*valfon*norz
                grad(3, 2) = diry*valfon*norz
                grad(3, 3) = dirz*valfon*norz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        grad(1, 1) = grad(1, 1)-dirx*valfon1*norx
                        grad(1, 1) = grad(1, 1)+sinb2*cosg*valfon2*sinb2*cosg*sign_tan
                        grad(1, 2) = grad(1, 2)-diry*valfon1*norx
                        grad(1, 2) = grad(1, 2)+sinb2*sing*valfon2*sinb2*cosg*sign_tan
                        grad(1, 3) = grad(1, 3)+dirz*valfon1*norx
                        grad(1, 3) = grad(1, 3)-cosb2*valfon2*sinb2*cosg*sign_tan

                        grad(2, 1) = grad(2, 1)-dirx*valfon1*nory
                        grad(2, 1) = grad(2, 1)+sinb2*cosg*valfon2*sinb2*sing*sign_tan
                        grad(2, 2) = grad(2, 2)-diry*valfon1*nory
                        grad(2, 2) = grad(2, 2)+sinb2*sing*valfon2*sinb2*sing*sign_tan
                        grad(2, 3) = grad(2, 3)+dirz*valfon1*nory
                        grad(2, 3) = grad(2, 3)-cosb2*valfon2*sinb2*sing*sign_tan

                        grad(3, 1) = grad(3, 1)+dirx*valfon1*norz
                        grad(3, 1) = grad(3, 1)-sinb2*cosg*valfon2*cosb2*sign_tan
                        grad(3, 2) = grad(3, 2)+diry*valfon1*norz
                        grad(3, 2) = grad(3, 2)-sinb2*sing*valfon2*cosb2*sign_tan
                        grad(3, 3) = grad(3, 3)-dirz*valfon1*norz
                        grad(3, 3) = grad(3, 3)+cosb2*valfon2*cosb2*sign_tan
                    end if
                end if

            else
                grad(1, 1) = dirx*(valfon-valfon1)*norx
                grad(1, 2) = diry*(valfon-valfon1)*norx
                grad(1, 3) = dirz*(valfon-valfon1)*norx
!
                grad(2, 1) = dirx*(valfon-valfon1)*nory
                grad(2, 2) = diry*(valfon-valfon1)*nory
                grad(2, 3) = dirz*(valfon-valfon1)*nory
!
                grad(3, 1) = dirx*(valfon-valfon1)*norz
                grad(3, 2) = diry*(valfon-valfon1)*norz
                grad(3, 3) = dirz*(valfon-valfon1)*norz
            end if
!
        else if (type .eq. 'SH') then
!
            if (abs(cosa) .gt. 0.d0) then
                grad(1, 1) = dirx*valfon*tanx
                grad(1, 1) = grad(1, 1)+dirx*valfon1*tanx
                grad(1, 2) = diry*valfon*tanx
                grad(1, 2) = grad(1, 2)+diry*valfon1*tanx
                grad(1, 3) = dirz*valfon*tanx
                grad(1, 3) = grad(1, 3)-dirz*valfon1*tanx
                grad(2, 1) = dirx*valfon*tany
                grad(2, 1) = grad(2, 1)+dirx*valfon1*tany
                grad(2, 2) = diry*valfon*tany
                grad(2, 2) = grad(2, 2)+diry*valfon1*tany
                grad(2, 3) = dirz*valfon*tany
                grad(2, 3) = grad(2, 3)-dirz*valfon1*tany
            else
                grad(1, 1) = dirx*(valfon-valfon1)*tanx
                grad(1, 2) = diry*(valfon-valfon1)*tanx
                grad(1, 3) = dirz*(valfon-valfon1)*tanx
                grad(2, 1) = dirx*(valfon-valfon1)*tany
                grad(2, 2) = diry*(valfon-valfon1)*tany
                grad(2, 3) = dirz*(valfon-valfon1)*tany
            end if
!
            grad(3, 1) = 0.d0
            grad(3, 2) = 0.d0
            grad(3, 3) = 0.d0
!
        end if
!
!        CALCUL DES DEFORMATIONS
        do indic1 = 1, 3
            do indic2 = 1, 3
                epsi(indic1, indic2) = .5d0*(grad(indic1, indic2)+grad(indic2, indic1))
            end do
        end do
!
!        CALCUL DES CONTRAINTES
        trace = 0.d0
        do indic1 = 1, 3
            trace = trace+epsi(indic1, indic1)
        end do
!
        do indic1 = 1, 3
            do indic2 = 1, 3
                if (indic1 .eq. indic2) then
                    sigma(indic1, indic2) = lambda*trace+2.d0*mu*epsi(indic1, indic2)
                else
                    sigma(indic1, indic2) = 2.d0*mu*epsi(indic1, indic2)
                end if
            end do
        end do
!
        nx = 0.0d0
        ny = 0.0d0
        nz = 0.0d0
!
!        --- CALCUL DE LA NORMALE AU POINT DE GAUSS IPG ---
!
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
                nx = nx+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                ny = ny+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                nz = nz+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
            end do
        end do
!
!        --- LE JACOBIEN EST EGAL A LA NORME DE LA NORMALE ---
!
        jac = sqrt(nx*nx+ny*ny+nz*nz)
!
!        --- CALCUL DE LA NORMALE UNITAIRE ---
!
        nux = nx/jac
        nuy = ny/jac
        nuz = nz/jac
!
!        --- TEST DU SENS DE LA NORMALE PAR RAPPORT A LA DIRECTION
!            DE L'ONDE
!
        scal = nux*dirx+nuy*diry+nuz*dirz
        if (scal .gt. 0.d0) then
            coedir = 1.d0
        else
            coedir = -1.d0
        end if
        coedir = -1.d0
        if (h .ne. r8vide()) then
            coedir = -1.d0
        else
            coedir = 0.d0
        end if
!
!        --- CALCUL DE V.N ---
!
        vondt(1) = 0.d0
        vondt(2) = 0.d0
        vondt(3) = 0.d0
!
        if (type .eq. 'P') then
            if (abs(cosa) .gt. 0.d0) then
                vondt(1) = -cele*valfon*dirx
                vondt(2) = -cele*valfon*diry
                vondt(3) = -cele*valfon*dirz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        vondt(1) = vondt(1)-cele*valfon1*dirx
                        vondt(1) = vondt(1)-cele2*valfon2*cosb2*cosg
                        vondt(2) = vondt(2)-cele*valfon1*diry
                        vondt(2) = vondt(2)-cele2*valfon2*cosb2*sing
                        vondt(3) = vondt(3)+cele*valfon1*dirz
                        vondt(3) = vondt(3)-cele2*valfon2*sinb2
                    end if
                end if

            else
                vondt(1) = -cele*(valfon+valfon1)*dirx
                vondt(2) = -cele*(valfon+valfon1)*diry
                vondt(3) = -cele*(valfon+valfon1)*dirz
            end if
        else if (type .eq. 'SV') then
            if (abs(cosa) .gt. 0.d0) then
                vondt(1) = -cele*valfon*norx
                vondt(2) = -cele*valfon*nory
                vondt(3) = -cele*valfon*norz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        vondt(1) = vondt(1)+cele*valfon1*norx
                        vondt(1) = vondt(1)-cele2*valfon2*sinb2*cosg*sign_tan
                        vondt(2) = vondt(2)+cele*valfon1*nory
                        vondt(2) = vondt(2)-cele2*valfon2*sinb2*sing*sign_tan
                        vondt(3) = vondt(3)-cele*valfon1*norz
                        vondt(3) = vondt(3)+cele2*valfon2*cosb2*sign_tan
                    end if
                end if
            else
                vondt(1) = -cele*(valfon+valfon1)*norx
                vondt(2) = -cele*(valfon+valfon1)*nory
                vondt(3) = -cele*(valfon+valfon1)*norz
            end if
        else if (type .eq. 'SH') then
            if (abs(cosa) .gt. 0.d0) then
                vondt(1) = -cele*valfon*tanx
                vondt(1) = vondt(1)-cele*valfon1*tanx
                vondt(2) = -cele*valfon*tany
                vondt(2) = vondt(2)-cele*valfon1*tany
            else
                vondt(1) = -cele*(valfon+valfon2)*tanx
                vondt(2) = -cele*(valfon+valfon2)*tany
            end if
            vondt(3) = 0.d0
        end if
!
!        --- CALCUL DE LA VITESSE NORMALE ET DE LA VITESSE TANGENTIELLE
        call pronor(nux, nuy, nuz, vondt, vondn)
!
!        --- CALCUL DU VECTEUR CONTRAINTE
!
        taux = -rho*(cp*vondn(1)+cs*vondt(1))*coef_amor
        tauy = -rho*(cp*vondn(2)+cs*vondt(2))*coef_amor
        tauz = -rho*(cp*vondn(3)+cs*vondt(3))*coef_amor
!
        if (zk8(ionde+1) (1:7) .eq. '&FOZERO') goto 98

        uondt(1) = 0.d0
        uondt(2) = 0.d0
        uondt(3) = 0.d0
!
        valfon1 = 0.d0
        valfon2 = 0.d0
        if (h .ne. r8vide()) then
            if (instd .lt. 0.d0) then
                valfon = 0.d0
            else
                call fointe('F ', zk8(ionde+1), 1, 'INST', [instd], valfon, ier)
            end if
            if (h2 .ne. r8vide()) then
                if (abs(cosa) .gt. 0.d0) then
                    if (instd1 .lt. 0.d0) then
                        valfon1 = 0.d0
                    else
                        call fointe('F ', zk8(ionde+1), 1, 'INST', [instd1], valfon1, ier)
                        valfon1 = valfon1*ra12
                    end if
                    if (instd2 .lt. 0.d0) then
                        valfon2 = 0.d0
                    else
                        call fointe('F ', zk8(ionde+1), 1, 'INST', [instd2], valfon2, ier)
                        valfon2 = valfon2*ra13
                    end if
                else
                    if (instd1 .lt. 0.d0) then
                        valfon1 = 0.d0
                    else
                        call fointe('F ', zk8(ionde+1), 1, 'INST', [instd1], valfon1, ier)
                    end if
                end if
            end if
        else
            lpar2(1) = 'X'
            lpar2(2) = 'INST'
            vpar2(1) = 1.0*param0
            vpar2(2) = zr(jinst)
            call fointe('F ', zk8(ionde+1), 2, lpar2, vpar2, valfon, ier)
            if (type .ne. 'P') then
                valfon = -valfon
            end if
        end if
        if (type .eq. 'P') then
            if (abs(cosa) .gt. 0.d0) then
                uondt(1) = valfon*dirx
                uondt(2) = valfon*diry
                uondt(3) = valfon*dirz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        uondt(1) = uondt(1)+valfon1*dirx
                        uondt(1) = uondt(1)+valfon2*cosb2*cosg
                        uondt(2) = uondt(2)+valfon1*diry
                        uondt(2) = uondt(2)+valfon2*cosb2*sing
                        uondt(3) = uondt(3)-valfon1*dirz
                        uondt(3) = uondt(3)+valfon2*sinb2
                    end if
                end if
            else
                uondt(1) = (valfon+valfon1)*dirx
                uondt(2) = (valfon+valfon1)*diry
                uondt(3) = (valfon+valfon1)*dirz
            end if
        else if (type .eq. 'SV') then
            if (abs(cosa) .gt. 0.d0) then
                uondt(1) = valfon*norx
                uondt(2) = valfon*nory
                uondt(3) = valfon*norz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        uondt(1) = uondt(1)-valfon1*norx
                        uondt(1) = uondt(1)+valfon2*sinb2*cosg*sign_tan
                        uondt(2) = uondt(2)-valfon1*nory
                        uondt(2) = uondt(2)+valfon2*sinb2*sing*sign_tan
                        uondt(3) = uondt(3)+valfon1*norz
                        uondt(3) = uondt(3)-valfon2*cosb2*sign_tan
                    end if
                end if
            else
                uondt(1) = (valfon+valfon1)*norx
                uondt(2) = (valfon+valfon1)*nory
                uondt(3) = (valfon+valfon1)*norz
            end if
        else if (type .eq. 'SH') then
            if (abs(cosa) .gt. 0.d0) then
                uondt(1) = valfon*tanx
                uondt(1) = uondt(1)+valfon1*tanx
                uondt(2) = valfon*tany
                uondt(2) = uondt(2)+valfon1*tany
            else
                uondt(1) = (valfon+valfon1)*tanx
                uondt(2) = (valfon+valfon1)*tany
            end if
            uondt(3) = 0.d0
        end if
!        --- CALCUL DES DEPLACEMENTS NORMAL ET TANGENTIEL
        call pronor(nux, nuy, nuz, uondt, uondn)
!        --- CALCUL DU VECTEUR CONTRAINTE
        taux = taux-(rip*uondn(1)+ris*uondt(1))
        tauy = tauy-(rip*uondn(2)+ris*uondt(2))
        tauz = tauz-(rip*uondn(3)+ris*uondt(3))
98      continue
!
!        --- CALCUL DU VECTEUR CONTRAINTE DU A UNE ONDE PLANE
!
        taondx = sigma(1, 1)*nux
        taondx = taondx+sigma(1, 2)*nuy
        taondx = taondx+sigma(1, 3)*nuz
!
        taondy = sigma(2, 1)*nux
        taondy = taondy+sigma(2, 2)*nuy
        taondy = taondy+sigma(2, 3)*nuz
!
        taondz = sigma(3, 1)*nux
        taondz = taondz+sigma(3, 2)*nuy
        taondz = taondz+sigma(3, 3)*nuz
!
!        --- CALCUL DU VECTEUR ELEMENTAIRE
!
        do i = 1, nno
            ii = 3*i-2
            zr(ires+ii-1) = zr(ires+ii-1)+(taux+coedir*taondx)*&
                            &zr(ivf+ldec+i-1)*jac*zr(ipoids+ipg-1)
            zr(ires+ii+1-1) = zr(ires+ii+1-1)+(tauy+coedir*taondy)*&
                              &zr(ivf+ldec+i-1)*jac*zr(ipoids+ipg-1)
            zr(ires+ii+2-1) = zr(ires+ii+2-1)+(tauz+coedir*taondz)*&
                              &zr(ivf+ldec+i-1)*jac*zr(ipoids+ipg-1)
        end do
!
    end do
!
99  continue
!
end subroutine
