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

subroutine lceib1(fami, kpg, ksp, imate, &
                  ndim, epsm, sref, sechm, hydrm, &
                  t, lambda, deuxmu, epsthe, kdess, &
                  bendo, gamma, seuil)
!
    implicit none
!
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/Behaviour_type.h"
    character(len=*) :: fami
    integer(kind=8) :: imate, ndim, t(3, 3), kpg, ksp
    real(kind=8) :: epsm(6), lambda, deuxmu, epsthe(2), kdess, bendo
    real(kind=8) :: gamma, seuil
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT ENDO_ISOT_BETON - INITIALISATION
!
! IN  COMPOR     : NOM DE LA LOI DE COMPORTEMENT
! IN  IMATE      : CODE MATERIAU
! IN  EPSM       : DEFORMATION AU TEMPS MOINS
! IN  TM         : TEMPERATURE A T-
! IN  TREF       : TEMPERATURE DE REFERENCE
! IN  SREF       : SECHAGE DE REFEERNCE
! IN  SECHM      : SECHAGE AU TEMPS -
! IN  HYDRM      : HYDRATATION AU TEMPS-
! OUT T          : TENSEUR DE PLACEMENT (PASSAGE VECT -> MATRICE)
! OUT LAMBDA
! OUT DEUXMU
! OUT ALPHA
! OUT KDESS
! OUT BENDO
! OUT GAMMA
! OUT SEUIL
! ----------------------------------------------------------------------
!
    integer(kind=8) :: icodre(3)
    character(len=16) :: nomres(3)
    integer(kind=8) :: i, k, ndimsi
    real(kind=8) :: valres(3), e, nu
    real(kind=8) :: sref, sechm, hydrm
    real(kind=8) :: k0, k1, sicr, trepsm, eps(6), kron(6)
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
!
    ndimsi = 2*ndim
    t(1, 1) = 1
    t(1, 2) = 4
    t(1, 3) = 5
    t(2, 1) = 4
    t(2, 2) = 2
    t(2, 3) = 6
    t(3, 1) = 5
    t(3, 2) = 6
    t(3, 3) = 3
!
!    LECTURE DES CARACTERISTIQUES DU MATERIAU
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres, valres, icodre, 1)
    call verift(fami, kpg, ksp, '-', imate, &
                epsth_=epsthe(1))
    call verift(fami, kpg, ksp, '+', imate, &
                epsth_=epsthe(2))
!
    e = valres(1)
    nu = valres(2)
!
    lambda = e*nu/(1.d0+nu)/(1.d0-2.d0*nu)
    deuxmu = e/(1.d0+nu)
!
!    LECTURE DES CARACTERISTIQUES DE RETRAIT ENDOGENE ET DESSICCATION
    nomres(1) = 'B_ENDOGE'
    nomres(2) = 'K_DESSIC'
    call rcvalb(fami, 1, 1, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres, valres, icodre, 0)
    if (icodre(1) .ne. 0) valres(1) = 0.d0
    if (icodre(2) .ne. 0) valres(2) = 0.d0
    bendo = valres(1)
    kdess = valres(2)
!
!    LECTURE DES CARACTERISTIQUES D'ENDOMMAGEMENT
    nomres(1) = 'D_SIGM_EPSI'
    nomres(2) = 'SYT'
    nomres(3) = 'SYC'
    call rcvalb(fami, 1, 1, '+', imate, &
                ' ', 'BETON_ECRO_LINE', 0, ' ', [0.d0], &
                3, nomres, valres, icodre, 0)
    if ((icodre(1) .ne. 0) .or. (icodre(2) .ne. 0)) then
        call utmess('F', 'ALGORITH4_51')
    end if
    gamma = -e/valres(1)
    k0 = valres(2)**2*(1.d0+gamma)/(2.d0*e)*(1.d0+nu-2.d0*nu**2)/( &
         1.d0+nu)
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
                k1 = valres(3)*(1.d0+gamma)*nu**2/(1.d0+nu)/(1.d0- &
                                                             2.d0*nu)-k0*e/(1.d0-2.d0*nu)/valres(3)
!      PASSAGE AUX DEFORMATIONS ELASTIQUES
                call r8inir(6, 0.d0, eps, 1)
                do k = 1, ndimsi
                    eps(k) = epsm(k)-(epsthe(1)-kdess*(sref-sechm)-bendo*hydrm) &
                            &*kron(k)
                end do
                trepsm = 0.d0
                do i = 1, ndim
                    trepsm = trepsm+eps(i)
                end do
                if (trepsm .gt. 0.d0) then
                    trepsm = 0.d0
                end if
                seuil = k0-k1*trepsm
            end if
        end if
    end if
!
end subroutine
