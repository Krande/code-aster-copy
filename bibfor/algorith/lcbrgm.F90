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
subroutine lcbrgm(ndim, typmod, imate, epsm, deps, &
                  vim, option, sig, vip, dsidpt, &
                  codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvala.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
! LOI DE COMPORTEMENT ELASTIQUE ENDO HETEROGENE
! (AVEC REGULARIS. DES CONTRAINTES)
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : NATURE DU MATERIAU
! IN  EPSM    : DEFORMATION EN T-
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  VIM     : VARIABLES INTERNES EN T-
! IN  OPTION  : OPTION DEMANDEE
!                 RIGI_MECA_TANG ->     DSIDEP
!                 FULL_MECA      -> SIG DSIDEP VIP
!                 RAPH_MECA      -> SIG        VIP
! OUT SIG     : CONTRAINTE
! OUT VIP     : VARIABLES INTERNES
!                 1   -> VALEUR DE L'ENDOMMAGEMENT
!                 2   -> ELASTIQUE (0) OU POINTE (1) RUPT AMORCAGE (2)
!                     -> RUPT PROPAGATION (3)
!                 3   -> CONTRAINTE RUPT AMORCAGE
!                 4   -> CONTRAINTE RUPT PROPAGATION
!                 5   -> NUMERO ELEMENT POINTE 1
!                 6   -> NUMERO ELEMENT POINTE 2 (SI RUPT AMORCAGE)
!                 7   -> IT DE NEWTON DE RUPTURE
!                 8   -> IT DE NEWTON COURANTE
!                 9   -> COORX POINTE DE FISSURE (APRES RUPT PROPA)
!                 10  -> COORY POINTE DE FISSURE (APRES RUPT PROPA)
! OUT DSIDPT  : MATRICE TANGENTE
! ----------------------------------------------------------------------
!
    character(len=8) :: typmod(*)
    character(len=16) :: option
    integer(kind=8) :: ndim, imate, codret
    real(kind=8) :: epsm(12), deps(12), vim(*)
    real(kind=8) :: sig(6), vip(*), dsidpt(6, 6, 2)
! ----------------------------------------------------------------------
!
!
!
!
    aster_logical :: cplan, resi, rigi
    integer(kind=8) :: ndimsi, k, l, etat
!
    real(kind=8) :: eps(6), epsr(6), treps, sigel(6)
    real(kind=8) :: kron(6)
    real(kind=8) :: fd, d, dm, e, nu, lambda, deuxmu
!
    integer(kind=8) :: icodre(2)
    character(len=16) :: nomres(2)
    real(kind=8) :: valres(2)
!
    real(kind=8) :: dmax
    blas_int :: b_incx, b_incy, b_n
    parameter(dmax=0.999999d0)
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!
!
! ======================================================================
!                            INITIALISATION
! ======================================================================
!
! -- OPTION ET MODELISATION
!
    codret = 0
    resi = option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
!
    cplan = (typmod(1) .eq. 'C_PLAN  ')
    ndimsi = 2*ndim
!
! -- LECTURE DES CARACTERISTIQUES ELASTIQUES
!
    nomres(1) = 'E'
    nomres(2) = 'NU'
    call rcvala(imate, ' ', 'ELAS', 0, ' ', &
                [0.d0], 2, nomres, valres, icodre, &
                1)
!
    e = valres(1)
    nu = valres(2)
    lambda = e*nu/(1.d0+nu)/(1.d0-2.d0*nu)
    deuxmu = e/(1.d0+nu)
!
! -- DEFORMATIONS
!
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm, b_incx, eps, b_incy)
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm(7), b_incx, epsr, b_incy)
    if (resi) then
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps, b_incx, eps, &
                   b_incy)
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps(7), b_incx, epsr, &
                   b_incy)
    end if
    do k = 1, ndimsi
        dsidpt(k, k, 1) = 0.d0
    end do
!
!
!
! ======================================================================
!                         CONTRAINTES ELASTIQUES
! ======================================================================
!
! -- CALCUL DES CONTRAINTES ELASTIQUES
!
    treps = eps(1)+eps(2)+eps(3)
    do k = 1, ndimsi
        sigel(k) = lambda*treps*kron(k)+deuxmu*eps(k)
    end do
!
! ======================================================================
!                 INTEGRATION DE LA LOI DE COMPORTEMENT
! ======================================================================
!
    if (resi) then
        dm = vim(1)
        etat = nint(vip(2))
        vip(3) = vim(3)
        vip(4) = vim(4)
!
        if (etat .eq. 3) then
            d = dmax
            etat = 3
        else if (etat .eq. 2) then
            d = dmax
            etat = 2
        else if (etat .eq. 1) then
            d = dm
            etat = 1
        else if (etat .eq. 0) then
            d = dm
            etat = 0
        end if
!
        do k = 1, ndimsi
            sig(k) = (1-d)*sigel(k)
        end do
!
        vip(1) = d
        vip(2) = etat
!
    else
        d = vim(1)
        etat = nint(vim(2))
    end if
!
!
! ======================================================================
!                            MATRICE TANGENTE
! ======================================================================
!
    if (rigi) then
        call r8inir(72, 0.d0, dsidpt, 1)
        fd = 1.d0-d
        do k = 1, 3
            do l = 1, 3
                dsidpt(k, l, 1) = fd*lambda
            end do
        end do
        do k = 1, ndimsi
            dsidpt(k, k, 1) = dsidpt(k, k, 1)+fd*deuxmu
        end do
        if (cplan) then
            do k = 1, ndimsi
                if (k .ne. 3) then
                    do l = 1, ndimsi
                        if (l .ne. 3) then
                            dsidpt(k, l, 1) = dsidpt(k, l, 1)-1.d0/dsidpt(3, 3, 1)*dsidpt(k, 3, 1&
                                              &)*dsidpt(3, l, 1)
                        end if
                    end do
                end if
            end do
        end if
    end if
!
end subroutine
