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
subroutine d1madp(fami, mater, instan, poum, kpg, &
                  ksp, angl, d1)
!.======================================================================
    implicit none
!
!      D1MADP --   CALCUL DE L'INVERSE DE LA MATRICE DE HOOKE
!                  POUR LES ELEMENTS MASSIFS 2D
!                  EN DEFORMATIONS PLANES OU AXISYMETRIQUES
!                  POUR DES MATERIAUX ISOTROPE, ORTHOTROPE
!                  ET ISOTROPE TRANSVERSE
!
!   ARGUMENT        E/S  TYPE         ROLE
!    MATER          IN     I        MATERIAU
!    INSTAN         IN     R        INSTANT DE CALCUL (0 PAR DEFAUT)
!    ANGL           IN     R        ANGLE DEFINISSANT LE REPERE
!                                   D'ORTHOTROPIE
!    D1(4,4)        OUT    R        INVERSE MATRICE DE HOOKE
!
!
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterfort/assert.h"
#include "asterfort/d1pa2d.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
    character(len=*) :: fami, poum
    integer(kind=8) :: kpg, ksp
    real(kind=8) :: angl, d1(4, *), instan
! -----  VARIABLES LOCALES
!-----------------------------------------------------------------------
    integer(kind=8) :: i, irep, j, mater, nbres, nbv
    real(kind=8) :: deux, e, e1, e2, e3, un, zero
!
!-----------------------------------------------------------------------
    parameter(nbres=7)
!
    integer(kind=8) :: icodre(nbres)
    character(len=8) :: nompar
    character(len=16) :: nomres(nbres)
    character(len=32) :: phenom
!
    real(kind=8) :: valres(nbres), valpar
    real(kind=8) :: passag(4, 4), d1orth(4, 4), work(4, 4)
    real(kind=8) :: nu, nu12, nu13, nu21, nu23, nu31, nu32
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! ---- INITIALISATIONS
!      ---------------
    zero = 0.0d0
    un = 1.0d0
    deux = 2.0d0
!
    nompar = 'INST'
    valpar = instan
!
    do i = 1, 4
        do j = 1, 4
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
! ----   INTERPOLATION DES COEFFICIENTS EN FONCTION DU TEMPS
!        -----------
        call rcvalb(fami, kpg, ksp, poum, mater, &
                    ' ', phenom, 1, nompar, [valpar], &
                    nbv, nomres, valres, icodre, 1)
!
        e = valres(1)
        nu = valres(2)
!
!          D1(1,1) =  (UN - NU*NU)/E
!          D1(1,2) = -(NU + NU*NU)/E
!          D1(1,3) =  D1(1,2)
        d1(1, 1) = un/e
        d1(1, 2) = -nu/e
        d1(1, 3) = -nu/e
!
        d1(2, 1) = d1(1, 2)
        d1(2, 2) = d1(1, 1)
        d1(2, 3) = d1(1, 2)
!
        d1(3, 1) = d1(1, 2)
        d1(3, 2) = d1(1, 2)
        d1(3, 3) = d1(1, 1)
!
        d1(4, 4) = deux*(1+nu)/e
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
        nbv = 7
!
! ----   INTERPOLATION DES COEFFICIENTS EN FONCTION DU TEMPS
!        -----------
        call rcvalb(fami, kpg, ksp, poum, mater, &
                    ' ', phenom, 1, nompar, [valpar], &
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
!
! ----   CALCUL DE LA MATRICE DE PASSAGE DU REPERE D'ORTHOTROPIE AU
! ----   REPERE GLOBAL POUR L'INVERSE DE LA MATRICE DE HOOKE
!        ---------------------------------------------------
        call d1pa2d(angl, irep, passag)
!
! ----   'INVERSE' DU TENSEUR D'ELASTICITE DANS LE REPERE GLOBAL :
! ----    D1_GLOB = PASSAG_T * D1_ORTH * PASSAG
! ----    (ON NE FAIT REELLEMENT LE PRODUIT QUE SI LA MATRICE
! ----     DE PASSAGE N'EST PAS L'IDENTITE)
!        ----------------------------------
        ASSERT((irep .eq. 1) .or. (irep .eq. 0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 4, 4, d1orth, passag, &
                        work, d1)
        else if (irep .eq. 0) then
            do i = 1, 4
                do j = 1, 4
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
        nbv = 4
!
! ----   INTERPOLATION DES COEFFICIENTS EN FONCTION DU TEMPS
!        -----------
        call rcvalb(fami, kpg, ksp, poum, mater, &
                    ' ', phenom, 1, nompar, [valpar], &
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
        d1orth(2, 2) = d1orth(1, 1)
        d1orth(2, 3) = d1orth(1, 3)
        d1orth(3, 1) = d1orth(1, 3)
        d1orth(3, 2) = d1orth(2, 3)
        d1orth(3, 3) = un/e3
        d1orth(4, 4) = deux*(un+nu12)/e1
!
! ----   CALCUL DE LA MATRICE DE PASSAGE DU REPERE D'ORTHOTROPIE AU
! ----   REPERE GLOBAL POUR L'INVERSE DE LA MATRICE DE HOOKE
!        ---------------------------------------------------
        call d1pa2d(angl, irep, passag)
!
! ----   'INVERSE' DU TENSEUR D'ELASTICITE DANS LE REPERE GLOBAL :
! ----    D1_GLOB = PASSAG_T * D1_ORTH * PASSAG
! ----    (ON NE FAIT REELLEMENT LE PRODUIT QUE SI LA MATRICE
! ----     DE PASSAGE N'EST PAS L'IDENTITE)
!        ----------------------------------
        ASSERT((irep .eq. 1) .or. (irep .eq. 0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 4, 4, d1orth, passag, &
                        work, d1)
        else if (irep .eq. 0) then
            do i = 1, 4
                do j = 1, 4
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
