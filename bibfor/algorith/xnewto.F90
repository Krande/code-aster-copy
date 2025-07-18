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

subroutine xnewto(elrefp, name, n, ndime, ptxx, &
                  ndim, tabco, tabls, ipp, ip, &
                  itemax, epsmax, ksi, exit, dekker)
    implicit none
!
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/utmess.h"
#include "asterfort/xdelt0.h"
#include "asterfort/xdelt2.h"
#include "asterfort/xdelt3.h"
#include "asterfort/xintva.h"
    integer(kind=8) :: ndime, ndim, ipp, ip, n(3)
    real(kind=8) :: ptxx(*), tabco(*), tabls(*)
    integer(kind=8) :: itemax
    real(kind=8) :: epsmax, ksi(ndime)
    character(len=6) :: name
    character(len=8) :: elrefp
    real(kind=8), intent(in), optional :: dekker(4*ndime)
    integer(kind=8), intent(inout), optional :: exit(2)
!
!             ALGORITHME DE NEWTON POUR CALCULER LES COORDONNEES
!                DE REFERENCE D'UN POINT PT MILIEU D'UNE ARETE
!
!     ENTREE
!       NUM     : NUMERO DE LA FONCTION A RESOUDRE (DANS XDELT1)
!       NDIM    : DIMENSION TOPOLOGIQUE DU MAILLAGE
!       COORSG  : COORDONNEES DES 3 NOEUDS DE L'ARETE
!       S       : ABSCISSE CURVILIGNE DU POINT SUR L'ARETE
!       ITEMAX  : NOMBRE MAXI D'ITERATIONS DE NEWTON
!       EPSMAX  : RESIDU POUR CONVERGENCE DE NEWTON
!       N       : LES INDICES DES NOEUX D'UNE FACE DANS L'ELEMENT PARENT
!       PMILIE  : LES COORDONNES DES POINTS MILIEUX
!       EXIT    : RETOUR EVENTUEL A UNE DECOUPE PRIMAIRE DIFFERENTE EN CAS
!                 D'ECHEC DU NEWTON
!
!     SORTIE
!       KSI     : COORDONNEES DE REFERENCE DU POINT
!     --------------------------------------------------------------
!
    real(kind=8) :: eps
    real(kind=8) :: test, epsrel, epsabs, refe, itermin
    integer(kind=8) :: iter, i, arete
    real(kind=8) :: zero
    parameter(zero=0.d0)
    real(kind=8) :: dist, dmin, intinf, intsup
    real(kind=8) :: ksi2(ndime), delta(ndime), ksim(ndime)
    data itermin/3/
!
! ------------------------------------------------------------------
!
!
! --- POINT DE DEPART
!
!  ATTENTION: ON SUPPOSE QUE LA FONCTION APPELANTE A DEJA
!  INITIALISE LE NEWTON EN AMONT
    do i = 1, ndime
        ksi2(i) = ksi(i)
    end do
!
    delta(:) = zero
    iter = 0
    epsabs = epsmax/100.d0
    epsrel = epsmax
    dmin = r8gaem()
    dist = 0.d0
    arete = 1
!
    if (present(dekker)) then
        call xintva(name, dekker, ptxx, ndime, intinf, &
                    intsup)
    end if
!
! --- DEBUT DE LA BOUCLE
!
20  continue
!-------------------------
!
!     FAIRE TANT QUE
!
!
! --- CALCUL DE L'INCREMENT
!
    if (name .eq. 'XMILFI') then
        call xdelt2(elrefp, n, ndime, ksi2, ptxx, &
                    ndim, tabco, tabls, ipp, ip, &
                    delta)
    elseif (name .eq. 'XINTAR') then
        call xdelt3(ndim, ksi2, tabls, delta(1))
    else if (name .eq. 'XINTER') then
        call xdelt0(elrefp, ndime, tabls, ptxx, ksi2(1), &
                    delta(1), arete)
    else if (name .eq. 'XMIFIS') then
        call xdelt0(elrefp, ndime, tabls, ptxx, ksi2(1), &
                    delta(1))
    else if (name .eq. 'XCENFI') then
        call xdelt0(elrefp, ndime, tabls, ptxx, ksi2(1), &
                    delta(1))
    end if
!
! --- ACTUALISATION
!
    do i = 1, ndime
        ksi2(i) = ksi2(i)-delta(i)
    end do
!
    iter = iter+1
!
!   ON VERIFIE POUR XMIFIS QUE LE NEWTON RESTE DANS LA FACE TRIA
!   DE RECHERCHE, SINON ON ACTIVE LA METHODE DE DEKKER
    if (name .eq. 'XMIFIS' .or. name .eq. 'XCENFI' .or. name .eq. 'XINTER') then
        if (present(dekker)) then
            if (iter .eq. itemax) then
!   DANS CERTAINS CAS, IL EST IMPOSSIBLE DE TROUVER UN POINT MILIEU SUR LA
!   FISSURE DANS LA FACE DU SOUS ELEMENT, ON EFFECTUE ALORS UNE APPROXIAMTION
!   LINEAIRE DE LA FISSURE SUR CETTE FACE
                if (name .eq. 'XINTER') then
                    ksi2(1) = ksi(1)
                else
                    ksi2(1) = 0.d0
                    if (exit(1) .le. 1) then
                        exit(1) = 1
                    end if
                end if
                goto 30
            end if
            if (ksi2(1) .gt. intsup) then
                ksi2(1) = ksi2(1)+delta(1)
                ksi2(1) = (ksi2(1)+intsup)/2.d0
            else if (ksi2(1) .lt. intinf) then
                ksi2(1) = ksi2(1)+delta(1)
                ksi2(1) = (ksi2(1)+intinf)/2.d0
            end if
        end if
    end if
!
    do i = 1, ndime
        dist = dist+delta(i)*delta(i)
    end do
    dist = sqrt(dist)
!
    if (dist .le. dmin) then
        do i = 1, ndime
            ksim(i) = ksi2(i)
        end do
    end if
!
! --- CALCUL DE LA REFERENCE POUR TEST DEPLACEMENTS
!
    refe = zero
    do i = 1, ndime
        refe = refe+ksi2(i)*ksi2(i)
    end do
    if (refe .le. epsrel) then
        refe = 1.d0
        eps = epsabs
    else
        eps = epsrel
    end if
!
! --- CALCUL POUR LE TEST DE CONVERGENCE
!
    test = zero
    do i = 1, ndime
        test = test+delta(i)*delta(i)
    end do
    test = sqrt(test/refe)
!
! --- EVALUATION DE LA CONVERGENCE
!
    if (iter .lt. itermin) goto 20
!
    if ((test .gt. eps) .and. (iter .lt. itemax)) then
        goto 20
    else if ((iter .ge. itemax) .and. (test .gt. eps)) then
        call utmess('F', 'XFEM_67')
        do i = 1, ndime
            ksi2(i) = ksim(i)
        end do
    end if
!
! --- FIN DE LA BOUCLE
!
    do i = 1, ndime
        ksi2(i) = ksi2(i)-delta(i)
    end do
!
30  continue
!   GESTION DU CAS NDIME<NDIM
    do i = 1, ndime
        ksi(i) = ksi2(i)
    end do
!    write(6,*)'CONVERGENCE DE ',name,' EN ',iter,' ITERATIONS'
end subroutine
