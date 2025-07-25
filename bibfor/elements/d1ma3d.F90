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

subroutine d1ma3d(fami, mater, instan, poum, kpg, &
                  ksp, angl, d1)
!.======================================================================
    implicit none
!
!     D1MA3D  --   CALCUL DE L'INVERSE DE LA MATRICE DE HOOKE
!                  POUR LES ELEMENTS MASSIFS EN 3D OU EN SERIE DE
!                  FOURIER POUR DES MATERIAUX ISOTROPE, ORTHOTROPE
!                  ET ISOTROPE TRANSVERSE
!
!   ARGUMENT        E/S  TYPE         ROLE
!    FAMI           IN     K*       FAMILLE DU POINT DE GAUSS
!    MATER          IN     I        MATERIAU
!    INSTAN         IN     R        INSTANT DE CALCUL (0 PAR DEFAUT)
!    POUM           IN     K1       T ou T+DT
!    KPG            IN     I        POINT DE GAUSS
!    KSP            IN     I        SOUS-POINT DE GAUSS
!    ANGL(3)        IN     R        ANGLES NAUTIQUES DEFINISSANT LE REPERE
!                                   D'ORTHOTROPIE
!    D1(6,6)        OUT    R        INVERSE DE LA MATRICE DE HOOKE
!
!
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/d1pa3d.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
    character(len=*) :: poum, fami
    integer(kind=8) :: kpg, ksp
    real(kind=8) :: angl(3), d1(6, 6), instan
! -----  VARIABLES LOCALES
!-----------------------------------------------------------------------
    integer(kind=8) :: i, irep, j, mater, nbres, nbv
    real(kind=8) :: coef1, coef2, coef3, deux, e, e1, e2
    real(kind=8) :: e3, un, zero
!-----------------------------------------------------------------------
    parameter(nbres=9)
!
    integer(kind=8) :: icodre(nbres)
    character(len=8) :: nompar(2)
    character(len=16) :: nomres(nbres)
    character(len=32) :: phenom
!
    real(kind=8) :: valres(nbres), valpar(1)
    real(kind=8) :: passag(6, 6), d1orth(6, 6), work(6, 6)
    real(kind=8) :: nu, nu12, nu21, nu13, nu23, nu31, nu32
    integer(kind=8) :: nbpar
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! ---- INITIALISATIONS
!      ---------------
    zero = 0.0d0
    un = 1.0d0
    deux = 2.0d0
!
    if (instan .eq. r8vide()) then
        nbpar = 0
    else
        nbpar = 1
        nompar(1) = 'INST'
        valpar(1) = instan
    end if
!
    do i = 1, 6
        do j = 1, 6
            d1(i, j) = zero
            d1orth(i, j) = zero
            work(i, j) = zero
        end do
    end do
!
! ---- RECUPERATION DU TYPE DU MATERIAU DANS PHENOM
!      --------------------------------------------
    call rccoma(mater, 'ELAS', 1, phenom, icodre(1))
!
!      ------------
! ---- CAS ISOTROPE
!      ------------
    if (phenom .eq. 'ELAS') then
!
        nomres(1) = 'E'
        nomres(2) = 'NU'
        nbv = 2
!
! ----   INTERPOLATION DES COEFFICIENTS EN FONCTION DE LA TEMPERATURE
! ----   ET DU TEMPS
!        -----------
        call rcvalb(fami, kpg, ksp, poum, mater, &
                    ' ', phenom, nbpar, nompar, [valpar], &
                    nbv, nomres, valres, icodre, 1)
!
        e = valres(1)
        nu = valres(2)
!
        coef1 = un/e
        coef2 = -nu/e
        coef3 = deux*(un+nu)/e
!
        d1(1, 1) = coef1
        d1(1, 2) = coef2
        d1(1, 3) = coef2
!
        d1(2, 1) = coef2
        d1(2, 2) = coef1
        d1(2, 3) = coef2
!
        d1(3, 1) = coef2
        d1(3, 2) = coef2
        d1(3, 3) = coef1
!
        d1(4, 4) = coef3
        d1(5, 5) = coef3
        d1(6, 6) = coef3
!
!      --------------
! ---- CAS ORTHOTROPE
!      --------------
    else if (phenom .eq. 'ELAS_ORTH') then
!
        nomres(1) = 'E_L'
        nomres(2) = 'E_T'
        nomres(3) = 'E_N'
        nomres(4) = 'NU_LT'
        nomres(5) = 'NU_LN'
        nomres(6) = 'NU_TN'
        nomres(7) = 'G_LT'
        nomres(8) = 'G_LN'
        nomres(9) = 'G_TN'
        nbv = 9
!
! ----   INTERPOLATION DES COEFFICIENTS EN FONCTION DE LA TEMPERATURE
! ----   ET DU TEMPS
!        -----------
        call rcvalb(fami, kpg, ksp, poum, mater, &
                    ' ', phenom, nbpar, nompar, [valpar], &
                    nbv, nomres, valres, icodre, 1)
!
        e1 = valres(1)
        e2 = valres(2)
        e3 = valres(3)
        nu12 = valres(4)
        nu13 = valres(5)
        nu23 = valres(6)
        nu21 = e2*nu12/e1
        nu31 = e3*nu13/e1
        nu32 = e3*nu23/e2
!
        d1orth(1, 1) = un/e1
        d1orth(1, 2) = -nu21/e2
        d1orth(1, 3) = -nu31/e3
        d1orth(2, 2) = un/e2
        d1orth(2, 3) = -nu32/e3
        d1orth(3, 3) = un/e3
        d1orth(2, 1) = d1orth(1, 2)
        d1orth(3, 1) = d1orth(1, 3)
        d1orth(3, 2) = d1orth(2, 3)
!
        d1orth(4, 4) = un/valres(7)
        d1orth(5, 5) = un/valres(8)
        d1orth(6, 6) = un/valres(9)
!
! ----   CALCUL DE LA MATRICE DE PASSAGE DU REPERE D'ORTHOTROPIE AU
! ----   REPERE GLOBAL POUR L'INVERSE DE LA MATRICE DE HOOKE
!        ---------------------------------------------------
        call d1pa3d(angl, irep, passag)
!
! ----   'INVERSE' DU TENSEUR D'ELASTICITE DANS LE REPERE GLOBAL :
! ----    D1_GLOB = PASSAG_T * D1_ORTH * PASSAG
! ----    (ON NE FAIT REELLEMENT LE PRODUIT QUE SI LA MATRICE
! ----     DE PASSAGE N'EST PAS L'IDENTITE)
!        ----------------------------------
        ASSERT((irep .eq. 1) .or. (irep .eq. 0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 6, 6, d1orth, passag, &
                        work, d1)
        else if (irep .eq. 0) then
            do i = 1, 6
                do j = 1, 6
                    d1(i, j) = d1orth(i, j)
                end do
            end do
        end if
!
!      -----------------------
! ---- CAS ISOTROPE-TRANSVERSE
!      -----------------------
    else if (phenom .eq. 'ELAS_ISTR') then
!
        nomres(1) = 'E_L'
        nomres(2) = 'E_N'
        nomres(3) = 'NU_LT'
        nomres(4) = 'NU_LN'
        nomres(5) = 'G_LN'
        nbv = 5
!
! ----   INTERPOLATION DES COEFFICIENTS EN FONCTION DE LA TEMPERATURE
! ----   ET DU TEMPS
!        -----------
        call rcvalb(fami, kpg, ksp, poum, mater, &
                    ' ', phenom, nbpar, nompar, [valpar], &
                    nbv, nomres, valres, icodre, 1)
!
        e1 = valres(1)
        e3 = valres(2)
        nu12 = valres(3)
        nu13 = valres(4)
        nu31 = e3*nu13/e1
!
        d1orth(1, 1) = un/e1
        d1orth(1, 2) = -nu12/e1
        d1orth(1, 3) = -nu31/e3
        d1orth(2, 1) = d1orth(1, 2)
        d1orth(2, 2) = un/e1
        d1orth(2, 3) = -nu31/e3
        d1orth(3, 1) = d1orth(1, 3)
        d1orth(3, 2) = d1orth(2, 3)
        d1orth(3, 3) = un/e3
        d1orth(4, 4) = deux*(un+nu12)/e1
        d1orth(5, 5) = un/valres(5)
        d1orth(6, 6) = d1orth(5, 5)
!
! ----   CALCUL DE LA MATRICE DE PASSAGE DU REPERE D'ORTHOTROPIE AU
! ----   REPERE GLOBAL POUR L'INVERSE DE LA MATRICE DE HOOKE
!        ---------------------------------------------------
        call d1pa3d(angl, irep, passag)
!
! ----   'INVERSE' DU TENSEUR D'ELASTICITE DANS LE REPERE GLOBAL :
! ----    D_GLOB = PASSAG_T * D_ORTH * PASSAG
! ----    (ON NE FAIT REELLEMENT LE PRODUIT QUE SI LA MATRICE
! ----     DE PASSAGE N'EST PAS L'IDENTITE)
!        ----------------------------------
        ASSERT((irep .eq. 1) .or. (irep .eq. 0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 6, 6, d1orth, passag, &
                        work, d1)
        else if (irep .eq. 0) then
            do i = 1, 6
                do j = 1, 6
                    d1(i, j) = d1orth(i, j)
                end do
            end do
        end if
!
    else
        call utmess('F', 'ELEMENTS_15', sk=phenom)
    end if
!.============================ FIN DE LA ROUTINE ======================
end subroutine
