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
subroutine pieigv(neps, tau, imate, vim, epsm, &
                  epspc, epsdc, typmod, etamin, etamax, &
                  copilo)
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/critev.h"
#include "asterfort/diagp3.h"
#include "asterfort/infniv.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/zerod2.h"
#include "blas/ddot.h"
    character(len=8) :: typmod(2)
    integer(kind=8) :: neps, imate
    real(kind=8) :: epsm(neps), epspc(neps), epsdc(neps), etamin, etamax, tau
    real(kind=8) :: vim(2)
    real(kind=8) :: copilo(2, 2)
! ----------------------------------------------------------------------
!     PILOTAGE LOI DE COMPORTEMENT ENDO_ISOT_BETON EN NON LOCAL
!     GRAD_VARI
!
! IN  NEPS    : DIMENSION DES DEFORMATIONS GENERALISEES
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  TAU     : 2ND MEMBRE DE L'EQUATION F(ETA)=TAU
! IN  IMATE   : NATURE DU MATERIAU
! IN  VIM     : VARIABLES INTERNES EN T-
! IN  EPSM    : CHAMP DE DEFORMATION EN T-
! IN  EPSPC   : INCREMENT DE DEFORMATION DU AUX CHARGES FIXES
! IN  EPSDC   : INCREMENT DE DEFORMATION DU AUX CHARGES PILOTEES
! IN  ETAMIN  : DONNEE UTILISATEUR DU MINIMUM DE ETA
! IN  ETAMAX  : DONNEE UTILISATEUR DU MAXIMUM DE ETA
! OUT COPILO  : COEFFICIENTS DE PILOTAGE :
!                F := COPILO(1,1)+COPILO(2,1)*ETA = TAU
!                F := COPILO(1,2)+COPILO(2,2)*ETA = TAU
! ----------------------------------------------------------------------
!
    aster_logical :: cplan
    integer(kind=8) :: ndim, ndimsi, k, iter, nitmax
    integer(kind=8) :: ifm, niv
    real(kind=8) :: trepsd, coplan, sigeld(6)
    real(kind=8) :: phim, phip, phid
    real(kind=8) :: tr(6), vecp(3, 3), rac2
    real(kind=8) :: fpd, dm, d, eta, epm(3)
    real(kind=8) :: e, nu, lambda, deuxmu, gamma, seuil, trepsm
    real(kind=8) :: k0, k1, sicr, r
    real(kind=8) :: epsp(7), epsd(7), x(4), y(4), z(4)
    real(kind=8) :: epstol
    real(kind=8) :: crit1, crit2
    real(kind=8) :: critp1, critp2
    real(kind=8) :: epsvp
    integer(kind=8) :: icodre(3), kpg, spt
    character(len=16) :: nomres(3)
    character(len=8) :: fami, poum
    real(kind=8) :: valres(3)
!
    real(kind=8) :: epsmax, etasup, etainf, epsnor
    real(kind=8) :: treinf, tresup
    real(kind=8) :: linter, epsto2
    real(kind=8) :: xs, ys, zs
    real(kind=8) :: x1, y1, z1
    real(kind=8) :: x2, y2, z2
!
    real(kind=8) :: kron(6)
    blas_int :: b_incx, b_incy, b_n
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
!----- GET INFO=1,2
    call infniv(ifm, niv)
!
    nitmax = 100
    epstol = 1.d-6
    epsvp = 1.d-6/abs(etamax-etamin)
    epsto2 = 1.d-2
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
!
!
! -- OPTION ET MODELISATION
    cplan = (typmod(1) .eq. 'C_PLAN  ')
    ndim = (neps-2)/3
    ndimsi = 2*ndim
    rac2 = sqrt(2.d0)
!
! -- CAS DE L'ENDOMMAGEMENT SATURE, on ne pilote pas
    if ((nint(vim(2)) .eq. 2)) then
        if (niv .eq. 2) then
            call utmess('I', 'PILOTAGE_2')
        end if
        goto 666
    end if
!
!
! -- LECTURE DES CARACTERISTIQUES THERMOELASTIQUES
    nomres(1) = 'E'
    nomres(2) = 'NU'
    call rcvalb(fami, kpg, spt, poum, imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres, valres, icodre, 1)
    e = valres(1)
    nu = valres(2)
    lambda = e*nu/(1.d0+nu)/(1.d0-2.d0*nu)
    deuxmu = e/(1.d0+nu)
!
!    LECTURE DES CARACTERISTIQUES DE REGULARISATION
    nomres(1) = 'PENA_LAGR'
    call rcvala(imate, ' ', 'NON_LOCAL', 0, ' ', &
                [0.d0], 1, nomres, valres, icodre, &
                0)
    if (icodre(1) .ne. 0) then
        valres(1) = 1.e3
    end if
    r = valres(1)
!
!    LECTURE DES CARACTERISTIQUES D'ENDOMMAGEMENT
    nomres(1) = 'D_SIGM_EPSI'
    nomres(2) = 'SYT'
    nomres(3) = 'SYC'
    call rcvalb(fami, kpg, spt, poum, imate, &
                ' ', 'BETON_ECRO_LINE', 0, ' ', [0.d0], &
                3, nomres, valres, icodre, 0)
    gamma = -e/valres(1)
    k0 = valres(2)**2*(1.d0+gamma)/(2.d0*e)*(1.d0+nu-2.d0*nu**2)/(1.d0+nu)
    if (nu .eq. 0) then
        if (icodre(3) .eq. 0) then
            call utmess('F', 'ALGORITH4_52')
        else
            seuil = k0
        end if
    else
        sicr = sqrt((1.d0+nu-2.d0*nu**2)/(2.d0*nu**2))*valres(2)
        if (icodre(3) .eq. 1) then
            seuil = k0
        else
            if (valres(3) .lt. sicr) then
                call utmess('F', 'ALGORITH4_53')
            else
                k1 = valres(3)*(1.d0+gamma)*nu**2/(1.d0+nu)/(1.d0-2.d0*nu)-k0*e/(1.d0-2.d0*nu)/v&
                     &alres(3)
                trepsm = 0.d0
                do k = 1, ndim
                    trepsm = trepsm+epsm(k)
                end do
                if (trepsm .gt. 0.d0) then
                    trepsm = 0.d0
                end if
                seuil = k0-k1*trepsm
            end if
        end if
    end if
!
!
!
!
!
!
!    ETAT MECANIQUE EN T-
!
!
    dm = vim(1)
    d = dm+tau
    fpd = (1+gamma)/(1+gamma*d)**2
!
! -- CAS DE L'ENDOMMAGEMENT QUI SATURERA, ON NE PILOTE PAS
    if (d .gt. 1.d0) then
        if (niv .eq. 2) then
            call utmess('I', 'PILOTAGE_2')
        end if
        goto 666
    end if
!
!
!
!
! -- CALCUL DES DEFORMATIONS EN PRESENCE DE CONTRAINTES PLANES
!
    if (cplan) then
        coplan = -nu/(1.d0-nu)
        epsm(3) = coplan*(epsm(1)+epsm(2))
        epspc(3) = coplan*(epspc(1)+epspc(2))
        epsdc(3) = coplan*(epsdc(1)+epsdc(2))
    end if
!
!
    do k = 1, 3
        epsp(k) = epsm(k)+epspc(k)
        epsd(k) = epsdc(k)
    end do
!
!
    do k = 4, ndimsi
        epsp(k) = epsm(k)+epspc(k)
        epsd(k) = epsdc(k)
    end do
!
    if (ndimsi .lt. 6) then
        do k = ndimsi+1, 6
            epsp(k) = 0.d0
            epsd(k) = 0.d0
        end do
    end if
!
    phim = epsm(ndimsi+2)+r*epsm(ndimsi+1)
    phip = epspc(ndimsi+2)+r*epspc(ndimsi+1)
    phid = epsdc(ndimsi+2)+r*epsdc(ndimsi+1)
!
    epsp(7) = phim+phip
    epsd(7) = phid
!
!
!
!
!
!
! Calcul du nombre de solutions sur un intervalle raisonnable
!
! - ON COMMENCE DONC PAR REGARDER LES VP DE EPSD
    tr(1) = epsd(1)
    tr(2) = epsd(4)/rac2
    tr(3) = epsd(5)/rac2
    tr(4) = epsd(2)
    tr(5) = epsd(6)/rac2
    tr(6) = epsd(3)
!
!
!
! -- DIAGONALISATION AVEC TRI EN VAL RELATIVE CROISSANT
    call diagp3(tr, vecp, epm)
!
!
! On prend la valeur absolue max des valeurs propres de EPSD
    epsmax = max(abs(epm(1)), abs(epm(3)))
!
! Si les valeurs propres sont trop petites, on ne pilote pas ce point
    if (epsmax .lt. epsvp) goto 666
!
!
!
! on "normalise" les deformations pilotees
!
    trepsd = epsd(1)+epsd(2)+epsd(3)
    do k = 1, ndimsi
        sigeld(k) = lambda*trepsd*kron(k)+deuxmu*epsd(k)
    end do
!
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    epsnor = 1.d0/sqrt(0.5d0*ddot(b_n, epsd, b_incx, sigeld, b_incy))
!
    do k = 1, 7
        epsd(k) = epsd(k)*epsnor
    end do
!
!
!
!
! CALIBRAGE DE L'INTERVALLE DE RECHERCHE POUR EVITER LES DIVERGENCES
! DE L'ALGORITHME DE RECHERCHE
!
! On repercute la normalisation sur les bornes de ETA pour
! definir l'intervalle de recherche
!
    etasup = etamax/epsnor
    etainf = -etamax/epsnor
!
!
! Test sur la valeur de la trace de la deformee pour eta=etainf
! pour s'assurer qu'elle ne diverge pas. On fixe une borne tr(eps)<1
    treinf = epsp(1)+epsp(2)+epsp(3)+etainf*(epsd(1)+epsd(2)+epsd(3))
    if (abs(treinf) .gt. 1.d0) then
!        WRITE(6,*) 'Modification de etainf  :',ETAINF
        etainf = (treinf/abs(treinf)-(epsp(1)+epsp(2)+epsp(3)))/(epsd(1)+epsd(2)+epsd(3))
!        WRITE(6,*) 'devient  :',ETAINF
    end if
!
    eta = etainf
    call critev(epsp, epsd, eta, lambda, deuxmu, &
                fpd, seuil, r*d, crit1, critp1)
!
!
! Test sur la valeur de la trace de la deformee pour eta=etasup
! pour s'assurer qu'elle ne diverge pas. On fixe une borne tr(eps)<1
!
    tresup = epsp(1)+epsp(2)+epsp(3)+etasup*(epsd(1)+epsd(2)+epsd(3))
    if (abs(tresup) .gt. 1.d0) then
!        WRITE(6,*) 'Modification de etasup  :',ETASUP
        etasup = (tresup/abs(tresup)-(epsp(1)+epsp(2)+epsp(3)))/(epsd(1)+epsd(2)+epsd(3))
!        WRITE(6,*) 'devient  :',ETASUP
    end if
!
    eta = etasup
    call critev(epsp, epsd, eta, lambda, deuxmu, &
                fpd, seuil, r*d, crit2, critp2)
!
!
!
! Longueur de l'intervalle
    linter = abs(etasup-etainf)
!
!
!
!###############################################################
! RECHERCHE DES SOLUTIONS SUR L'INTERVALLE [-ETASUP,ETASUP]
!###############################################################
!
!
! CAS A 0 SOLUTION
!
! on reste en dessous du seuil sur l'intervalle
    if ((crit1 .lt. 0.d0) .and. (crit2 .lt. 0.d0)) then
!        WRITE(6,*) 'cas 1'
        goto 666
    end if
!
! on reste au dessus du seuil sur l'intervalle,
!        on utilise la convexite pour le voir
!
    if (((crit1 .gt. 0.d0) .and. (critp1 .gt. (-crit1/linter))) .or. &
        ((crit2 .gt. 0.d0) .and. (critp2 .lt. (crit2/linter)))) then
!        WRITE(6,*) 'cas 2'
        goto 666
    end if
!
!
! CAS A 1 SOLUTION
!
    if ((crit1 .lt. 0.d0) .and. (crit2 .gt. 0.d0)) then
!        WRITE(6,*) 'cas 3'
        x(1) = etainf
        y(1) = crit1
        z(1) = critp1
        x(2) = etasup
        y(2) = crit2
        z(2) = critp2
!
        x(3) = x(1)
        y(3) = y(1)
        z(3) = z(1)
!
        do iter = 1, nitmax
!
            if (abs(y(3)) .le. epstol*seuil) goto 201
            if (abs(z(1)-z(2)) .lt. epsto2*abs(z(2))) then
                x(3) = (-y(3)+z(3)*x(3))/z(3)
                goto 555
            end if
            call zerod2(x, y, z)
555         continue
!
            call critev(epsp, epsd, x(3), lambda, deuxmu, &
                        fpd, seuil, r*d, y(3), z(3))
!
        end do
        call utmess('F', 'PILOTAGE_87')
201     continue
!
        copilo(2, 1) = z(3)/epsnor
        copilo(1, 1) = tau-x(3)*copilo(2, 1)*epsnor
        copilo(1, 2) = r8vide()
        copilo(2, 2) = r8vide()
!
        goto 999
!
    end if
!
    if ((crit1 .gt. 0.d0) .and. (crit2 .lt. 0.d0)) then
!        WRITE(6,*) 'cas 4'
        x(2) = etainf
        y(2) = crit1
        z(2) = critp1
        x(1) = etasup
        y(1) = crit2
        z(1) = critp2
!
        x(3) = x(1)
        y(3) = y(1)
        z(3) = z(1)
        do iter = 1, nitmax
            if (abs(y(3)) .le. epstol*seuil) goto 203
!
            if (abs(z(1)-z(2)) .lt. epsto2*abs(z(2))) then
                x(3) = (-y(3)+z(3)*x(3))/z(3)
                goto 556
            end if
            call zerod2(x, y, z)
556         continue
!
            call critev(epsp, epsd, x(3), lambda, deuxmu, &
                        fpd, seuil, r*d, y(3), z(3))
        end do
        call utmess('F', 'PILOTAGE_87')
203     continue
!
        copilo(2, 1) = z(3)/epsnor
        copilo(1, 1) = tau-x(3)*copilo(2, 1)*epsnor
        copilo(1, 2) = r8vide()
        copilo(2, 2) = r8vide()
!
        goto 999
!
    end if
!
!
!
! CAS A 2 OU 0 SOLUTIONS
!
    if (((crit1 .gt. 0.d0) .and. (critp1 .lt. (-crit1/linter))) .and. &
        ((crit2 .gt. 0.d0) .and. (critp2 .gt. (crit2/linter)))) then
!
!        WRITE(6,*) 'cas 5'
!
!
! il faut chercher s'il y a une valeur dans l'intervalle qui donne une
! valeur du critere negative
! s'il y en a une, il y a 2 solutions, sinon 0 solution
!
! On utilise les tangentes pour aller vers le "minimum"
! on s'arrete quand le critere est negatif, on se fiche
! de trouver exactement le minimum
!
        x1 = etainf
        y1 = crit1
        z1 = critp1
        x2 = etasup
        y2 = crit2
        z2 = critp2
!
        ys = y1
!
        iter = 0
!
750     continue
!
        if (iter .lt. nitmax) then
            xs = (y2-y1+z1*x1-z2*x2)/(z1-z2)
            call critev(epsp, epsd, xs, lambda, deuxmu, &
                        fpd, seuil, r*d, ys, zs)
            if (ys .lt. 0.d0) goto 751
!
            if (zs .gt. 0.d0) then
                x2 = xs
                y2 = ys
                z2 = zs
                linter = x2-x1
                if ((z1 .gt. (-y1/linter)) .or. (z2 .lt. (y2/linter))) then
                    goto 666
                end if
                goto 750
            else
                x1 = xs
                y1 = ys
                z1 = zs
                linter = x2-x1
                if ((z1 .gt. (-y1/linter)) .or. (z2 .lt. (y2/linter))) then
                    goto 666
                end if
                goto 750
            end if
        else
            goto 666
        end if
!
751     continue
!
!
! il y a une solution sur [ETAINF,XS] et une sur [XS,ETASUP]
!
!
! Calcul de la solution sur [XS,ETASUP]
        x(1) = xs
        y(1) = ys
        z(1) = zs
        x(2) = etasup
        y(2) = crit2
        z(2) = critp2
!
        x(3) = x(1)
        y(3) = y(1)
        z(3) = z(1)
!
        do iter = 1, nitmax
            if (abs(y(3)) .le. epstol*seuil) goto 205
!
            if (abs(z(1)-z(2)) .lt. epsto2*abs(z(2))) then
                x(3) = (-y(3)+z(3)*x(3))/z(3)
                goto 557
            end if
            call zerod2(x, y, z)
557         continue
!
            call critev(epsp, epsd, x(3), lambda, deuxmu, &
                        fpd, seuil, r*d, y(3), z(3))
!
        end do
        call utmess('F', 'PILOTAGE_87')
205     continue
!
        copilo(2, 1) = z(3)/epsnor
        copilo(1, 1) = tau-x(3)*copilo(2, 1)*epsnor
!
!
! Calcul de la solution sur [-ETASUP,XS]
        x(1) = xs
        y(1) = ys
        z(1) = zs
        x(2) = etainf
        y(2) = crit1
        z(2) = critp1
!
        x(3) = x(1)
        y(3) = y(1)
        z(3) = z(1)
!
        do iter = 1, nitmax
            if (abs(y(3)) .le. epstol*seuil) goto 207
            if (abs(z(1)-z(2)) .lt. epsto2*abs(z(2))) then
                x(3) = (-y(3)+z(3)*x(3))/z(3)
                goto 558
            end if
            call zerod2(x, y, z)
558         continue
!
            call critev(epsp, epsd, x(3), lambda, deuxmu, &
                        fpd, seuil, r*d, y(3), z(3))
!
        end do
        call utmess('F', 'PILOTAGE_87')
207     continue
!
        copilo(2, 2) = z(3)/epsnor
        copilo(1, 2) = tau-x(3)*copilo(2, 2)*epsnor
!
        goto 999
!
!
    end if
!
666 continue
    copilo(1, 1) = 0.d0
    copilo(2, 1) = 0.d0
    copilo(1, 2) = r8vide()
    copilo(2, 2) = r8vide()
!
!
999 continue
!
end subroutine
