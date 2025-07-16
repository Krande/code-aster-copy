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

subroutine te0499(option, nomte)
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
#include "asterfort/rcvalb.h"
#include "asterfort/trigom.h"
#include "asterfort/vff2dn.h"
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
! Elements: 2D
! Option: ONDE_PLAN
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: fami, poum
    integer(kind=8) :: icodre(5), kpg
    character(len=1) :: type
    real(kind=8) :: poids, nx, ny, valres(5), e, nu, lambda, mu, cp, cs
    real(kind=8) :: rho, taux, tauy, nux, nuy, scal, sign_tan, nortan
    real(kind=8) :: sigma(2, 2), epsi(2, 2), grad(2, 2)
    real(kind=8) :: xgg(4), ygg(4), vondn(2), vondt(2), uondn(2), uondt(2)
    real(kind=8) :: taondx, taondy, norx, nory, norz, dirx, diry, dirz, cele, cele2
    real(kind=8) :: trace, norm, jac, coef_amor
    real(kind=8) :: param0, param, h, h2, instd, ris, rip, l0, usl0
    real(kind=8) :: a2, b2, sina, cosa, sinb2, cosb2, rc1c2, ra12, ra13, sing, cosg
    real(kind=8) :: coedir, typer, valfon, tanx, tany
    real(kind=8) :: xsv, zsv, dist1, dist2, instd1, instd2, x0, z0, x1, z1
    real(kind=8) :: kr, nr
    real(kind=8) :: valfon1, valfon2, param1, param2
    integer(kind=8) :: nno, npg, ipoids, ivf, idfde, igeom
    integer(kind=8) :: ivectu, k, i, mater
    integer(kind=8) :: ier, ii, imate, indic1, indic2, iondc, ionde
    integer(kind=8) :: j, jgano, jinst, ndim, nnos, ndim2
    character(len=8) :: nompar(2), lpar2(2)
    real(kind=8) :: vpar2(2)
    real(kind=8) :: xygau(2)
    integer(kind=8) :: idecpg, idecno
    character(len=16), parameter :: nomres(5) = (/'E        ', 'NU       ', &
                                                  'RHO      ', &
                                                  'COEF_AMOR', 'LONG_CARA'/)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'ONDE_PLAN')
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PONDPLA', 'L', ionde)
    call jevech('PONDPLR', 'L', iondc)
    call jevech('PINSTR', 'L', jinst)
    call jevech('PVECTUR', 'E', ivectu)
!
    if (zk8(ionde) (1:7) .eq. '&FOZERO') goto 99
!
!     --- INITIALISATION DE SIGMA
!
    do i = 1, 2
        do j = 1, 2
            sigma(i, j) = 0.d0
        end do
    end do
!
    mater = zi(imate)
    fami = 'RIGI'
    poum = '+'
    ASSERT(ndim .ne. 2)
    ndim2 = ndim+1
!
    nompar(1) = 'X'
    nompar(2) = 'Y'
!
!     --- CARACTERISTIQUES DE L'ONDE PLANE
!
    dirx = zr(iondc)
    dirz = zr(iondc+1)
    typer = zr(iondc+3)
    x0 = zr(iondc+4)
    z0 = zr(iondc+5)
    x1 = zr(iondc+7)
    z1 = zr(iondc+8)
!
    if (typer .eq. 0.d0) type = 'P'
    if (typer .eq. 1.d0) type = 'S'
!
!     --- CALCUL DU VECTEUR DIRECTEUR UNITAIRE DE L'ONDE PLANE
!
    norm = sqrt(dirx**2.d0+dirz**2.d0)
    dirx = dirx/norm
    diry = 0.d0
    dirz = dirz/norm

    if (x0 .ne. r8vide()) then
        h = x0*dirx+z0*dirz
        if (x1 .ne. r8vide()) then
            h2 = x1*dirx+z1*dirz
        else
            h2 = r8vide()
        end if
    else
        h = r8vide()
    end if
!
!     CALCUL DU REPERE ASSOCIE A L'ONDE

    sing = 0.d0
    cosg = 1.d0

    tanx = -sing
    tany = cosg

    sign_tan = tany/abs(tany)

!
    nortan = sqrt(tanx**2.d0+tany**2.d0)

    tanx = tanx/nortan
    tany = tany/nortan

    norx = -tany*dirz
    nory = tanx*dirz
    norz = -tanx*diry+tany*dirx
!
!
    do kpg = 1, npg
        xgg(kpg) = 0.d0
        ygg(kpg) = 0.d0
    end do
!
!    write(6,*) 'npg=',npg,'nno=',nno
    do kpg = 1, npg
        k = (kpg-1)*nno
        do i = 1, nno
            ii = 2*i-1
            xgg(kpg) = xgg(kpg)+zr(igeom+ii-1)*zr(ivf+k+i-1)
            ygg(kpg) = ygg(kpg)+zr(igeom+ii)*zr(ivf+k+i-1)
        end do
!      write(6,*) 'kp=',kp,'xgg=',xgg(kp),'ygg=',ygg(kp)
    end do
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do kpg = 1, npg
!
        idecpg = nno*(kpg-1)-1
! ----- Coordinates for current Gauss point
        xygau(:) = 0.d0
        do i = 1, nno
            idecno = ndim2*(i-1)-1
            do j = 1, ndim2
                xygau(j) = xygau(j)+zr(ivf+i+idecpg)*zr(igeom+j+idecno)
            end do
        end do
!
        call rcvalb(fami, kpg, 1, poum, mater, &
                    ' ', 'ELAS', 2, nompar, xygau, &
                    4, nomres, valres, icodre, 1)
!       appel LONG_CARA en iarret = 0
        call rcvalb(fami, kpg, 1, poum, mater, &
                    ' ', 'ELAS', 2, nompar, xygau, &
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
        sina = dirx
        cosa = dirz
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
                else

                    if ((abs(sina*rc1c2)-1.d0) .gt. r8prem()) then
                        call utmess('F', 'ALGORITH2_82')
                    end if

                    !          b2=asin(sina*rc1c2)
                    b2 = trigom('ASIN', sina*rc1c2)
                    cosb2 = cos(b2)
                    sinb2 = sin(b2)

                    ! Coefficients de réflexions (Calcul intégré au script).
                    kr = 1.d0/rc1c2
                    nr = (kr)**2*sin(2.d0*b2)*sin(2.d0*a2)+cos(2.d0*a2)**2
                    ra12 = (kr**2*sin(2.d0*a2)*sin(2.d0*b2)-cos(2.d0*a2)**2)/nr
                    ra13 = -2.d0*kr*sin(2.d0*a2)*cos(2.d0*a2)/nr
                end if

! Calcul des bons paramètres dist à insérer dans le calcul.

                dist1 = x0*sina+(2.d0*z1-z0)*(-cosa)
                if (type .eq. 'P') then
                    zsv = z1+cosb2*(z1-z0)/rc1c2/cosa
                    xsv = x0+(sina/cosa)*(z1-z0)-sinb2*(z1-z0)/rc1c2/cosa
                else
                    zsv = z1+cosb2*(z1-z0)*rc1c2/cosa
                    xsv = x0+(sina/cosa)*(z1-z0)-sinb2*(z1-z0)*rc1c2/cosa
                end if
                dist2 = xsv*sinb2+zsv*(-cosb2)
            end if
        end if
!
        k = (kpg-1)*nno
!
!        CALCUL DU CHARGEMENT PAR ONDE PLANE
        param0 = dirx*xgg(kpg)+dirz*ygg(kpg)

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
                    param1 = sina*xgg(kpg)-cosa*ygg(kpg)
                    param = param1-dist1
                    instd1 = zr(jinst)-param/cele
                    if (instd1 .lt. 0.d0) then
                        valfon1 = 0.d0
                    else
                        call fointe('F ', zk8(ionde), 1, 'INST', [instd1], valfon1, ier)
                        valfon1 = valfon1*ra12
                    end if
                    param2 = sinb2*xgg(kpg)-cosb2*ygg(kpg)
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

!
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
                grad(1, 2) = dirz*valfon*dirx

                grad(2, 1) = dirx*valfon*dirz
                grad(2, 2) = dirz*valfon*dirz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        grad(1, 1) = grad(1, 1)+dirx*valfon1*dirx
                        grad(1, 1) = grad(1, 1)+sinb2*valfon2*cosb2
                        grad(1, 2) = grad(1, 2)-dirz*valfon1*dirx
                        grad(1, 2) = grad(1, 2)-cosb2*valfon2*cosb2

                        grad(2, 1) = grad(2, 1)-dirx*valfon1*dirz
                        grad(2, 1) = grad(2, 1)+sinb2*valfon2*sinb2
                        grad(2, 2) = grad(2, 2)+dirz*valfon1*dirz
                        grad(2, 2) = grad(2, 2)-cosb2*valfon2*sinb2
                    end if
                end if
            else
                grad(1, 1) = dirx*(valfon-valfon1)*dirx
                grad(1, 2) = dirz*(valfon-valfon1)*dirx
                grad(2, 1) = dirx*(valfon-valfon1)*dirz
                grad(2, 2) = dirz*(valfon-valfon1)*dirz
            end if
        else if (type .eq. 'S') then
!
            if (abs(cosa) .gt. 0.d0) then
                grad(1, 1) = dirx*valfon*norx
                grad(1, 2) = dirz*valfon*norx

                grad(2, 1) = dirx*valfon*norz
                grad(2, 2) = dirz*valfon*norz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        grad(1, 1) = grad(1, 1)-dirx*valfon1*norx
                        grad(1, 1) = grad(1, 1)+sinb2*cosg*valfon2*sinb2*cosg*sign_tan
                        grad(1, 2) = grad(1, 2)+dirz*valfon1*norx
                        grad(1, 2) = grad(1, 2)-cosb2*valfon2*sinb2*cosg*sign_tan

                        grad(2, 1) = grad(2, 1)+dirx*valfon1*norz
                        grad(2, 1) = grad(2, 1)-sinb2*cosg*valfon2*cosb2*sign_tan
                        grad(2, 2) = grad(2, 2)-dirz*valfon1*norz
                        grad(2, 2) = grad(2, 2)+cosb2*valfon2*cosb2*sign_tan
                    end if
                end if
            else
                grad(1, 1) = dirx*(valfon-valfon1)*norx
                grad(1, 2) = dirz*(valfon-valfon1)*norx
                grad(2, 1) = dirx*(valfon-valfon1)*norz
                grad(2, 2) = dirz*(valfon-valfon1)*norz
            end if
!
        end if

!
!        CALCUL DES DEFORMATIONS
        do indic1 = 1, 2
            do indic2 = 1, 2
                epsi(indic1, indic2) = .5d0*(grad(indic1, indic2)+grad(indic2, indic1))
            end do
        end do
!
!        CALCUL DES CONTRAINTES
        trace = 0.d0
        do indic1 = 1, 2
            trace = trace+epsi(indic1, indic1)
        end do
        do indic1 = 1, 2
            do indic2 = 1, 2
                if (indic1 .eq. indic2) then
                    sigma(indic1, indic2) = lambda*trace+2.d0*mu*epsi(indic1, indic2)
                else
                    sigma(indic1, indic2) = 2.d0*mu*epsi(indic1, indic2)
                end if
            end do
        end do
!
        call vff2dn(ndim, nno, kpg, ipoids, idfde, &
                    zr(igeom), nx, ny, poids)
!
        jac = sqrt(nx*nx+ny*ny)
!
!        --- CALCUL DE LA NORMALE UNITAIRE ---
!
        nux = nx/jac
        nuy = ny/jac
!
!        --- TEST DU SENS DE LA NORMALE PAR RAPPORT A LA DIRECTION
!            DE L'ONDE
!
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
!
        if (type .eq. 'P') then
            if (abs(cosa) .gt. 0.d0) then
                vondt(1) = -cele*valfon*dirx
                vondt(2) = -cele*valfon*dirz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        vondt(1) = vondt(1)-cele*valfon1*dirx
                        vondt(1) = vondt(1)-cele2*valfon2*cosb2
                        vondt(2) = vondt(2)+cele*valfon1*dirz
                        vondt(2) = vondt(2)-cele2*valfon2*sinb2
                    end if
                end if

            else
                vondt(1) = -cele*(valfon+valfon1)*dirx
                vondt(2) = -cele*(valfon+valfon1)*dirz
            end if
        else if (type .eq. 'S') then
            if (abs(cosa) .gt. 0.d0) then
                vondt(1) = -cele*valfon*norx
                vondt(2) = -cele*valfon*norz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        vondt(1) = vondt(1)+cele*valfon1*norx
                        vondt(1) = vondt(1)-cele2*valfon2*sinb2*cosg*sign_tan
                        vondt(2) = vondt(2)-cele*valfon1*norz
                        vondt(2) = vondt(2)+cele2*valfon2*cosb2*sign_tan
                    end if
                end if
            else
                vondt(1) = -cele*(valfon+valfon1)*norx
                vondt(2) = -cele*(valfon+valfon1)*norz
            end if
        end if
!
        scal = nux*vondt(1)+nuy*vondt(2)
!
!        --- CALCUL DE LA VITESSE NORMALE ET DE LA VITESSE TANGENCIELLE
!
        vondn(1) = nux*scal
        vondn(2) = nuy*scal
!
        vondt(1) = vondt(1)-vondn(1)
        vondt(2) = vondt(2)-vondn(2)
!
!        --- CALCUL DU VECTEUR CONTRAINTE
!
        taux = -rho*(cp*vondn(1)+cs*vondt(1))*coef_amor
        tauy = -rho*(cp*vondn(2)+cs*vondt(2))*coef_amor
!
        if (zk8(ionde+1) (1:7) .eq. '&FOZERO') goto 98
        uondt(1) = 0.d0
        uondt(2) = 0.d0
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
                uondt(2) = valfon*dirz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        uondt(1) = uondt(1)+valfon1*dirx
                        uondt(1) = uondt(1)+valfon2*cosb2
                        uondt(2) = uondt(2)-valfon1*dirz
                        uondt(2) = uondt(2)+valfon2*sinb2
                    end if
                end if
            else
                uondt(1) = (valfon+valfon1)*dirx
                uondt(2) = (valfon+valfon1)*dirz
            end if
        else if (type .eq. 'S') then
            if (abs(cosa) .gt. 0.d0) then
                uondt(1) = valfon*norx
                uondt(2) = valfon*norz

                if (h .ne. r8vide()) then
                    if (h2 .ne. r8vide()) then
                        uondt(1) = uondt(1)-valfon1*norx
                        uondt(1) = uondt(1)+valfon2*sinb2*cosg*sign_tan
                        uondt(2) = uondt(2)+valfon1*norz
                        uondt(2) = uondt(2)-valfon2*cosb2*sign_tan
                    end if
                end if
            else
                uondt(1) = (valfon+valfon1)*norx
                uondt(2) = (valfon+valfon1)*norz
            end if
        end if
        scal = nux*uondt(1)+nuy*uondt(2)
        uondn(1) = nux*scal
        uondn(2) = nuy*scal
        uondt(1) = uondt(1)-uondn(1)
        uondt(2) = uondt(2)-uondn(2)
!
        taux = taux-(rip*uondn(1)+ris*uondt(1))
        tauy = tauy-(rip*uondn(2)+ris*uondt(2))
98      continue
!
!        --- CALCUL DU VECTEUR CONTRAINTE DU A UNE ONDE PLANE
!
        taondx = sigma(1, 1)*nux
        taondx = taondx+sigma(1, 2)*nuy
!
        taondy = sigma(2, 1)*nux
        taondy = taondy+sigma(2, 2)*nuy

!
!        --- CALCUL DU VECTEUR ELEMENTAIRE
!
        do i = 1, nno
            ii = 2*i-1
            zr(ivectu+ii-1) = zr(ivectu+ii-1)+(taux+coedir*taondx)*&
                                  &zr(ivf+k+i-1)*poids
            zr(ivectu+ii+1-1) = zr(ivectu+ii+1-1)+(tauy+coedir*taondy)*&
                                      &zr(ivf+k+i-1)*poids
        end do
!
    end do
!
99  continue
!
end subroutine
