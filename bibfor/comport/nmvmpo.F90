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
subroutine nmvmpo(fami, npg, nno, option, nc, &
                  xl, wgauss, icodma, sect, u, &
                  du, contm, contp, fl, klv)
!
!
    implicit none
    character(len=*) :: fami, option
    integer(kind=8) :: npg, nno, nc, icodma
    real(kind=8) :: xl, sect(*), u(nno*nc), du(nno*nc), fl(nno*nc), klv(*)
    real(kind=8) :: contm(npg*nc), contp(npg*nc), wgauss(npg)
!
#include "asterf_types.h"
#include "asterfort/jsd1ff.h"
#include "asterfort/mavec.h"
#include "asterfort/moytem.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utbtab.h"
#include "asterfort/verifm.h"
#include "blas/dscal.h"
!
! --------------------------------------------------------------------------------------------------
!
!                       COMPORTEMENT ELAS POUR LES MECA_POU_D_TG
!
! --------------------------------------------------------------------------------------------------
!
!   IN :
!       fami        : famille de point de gauss
!       npg         : nombre de point de gauss
!       nno         : nombre de noeuds
!       option      : RAPH_MECA  FULL_MECA RIGI_MECA_TANG
!       nc          : nombre de composantes
!       xl          : longueur de l element
!       wgauss      : poids des points de Gauss
!       icodma      : adresse du materiau code
!       sect        : caracteristiques de la section
!       u           : vecteur deplacement a l'instant precedent
!       du          : vecteur accroissement de deplacement
!       contm       : contraintes a l'instant precedent
!   OUT :
!       contp       : contraintes à l'instant actuel
!       fl          : force nodale = bt*contp
!       klv         : matrice de rigidité tangente
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: codres(2), itemp, iret
    character(len=2) :: nomres(2)
    aster_logical :: vecteu, matric
    integer(kind=8) :: dimklv, kp, kk, i, j, k
    real(kind=8) :: eps(nc), deps(nc), fg(nno*nc), sigp(nc), sigm(nc)
    real(kind=8) :: e, nu, g, phiy, phiz, xls2, epsthf, epsthd
    real(kind=8) :: aa, xiy, xiz, alfay, alfaz, xjx, xjg
    real(kind=8) :: valres(2)
    real(kind=8) :: temp
    real(kind=8) :: hoel(nc, nc), hota(nc, nc), d1b(nc, nno*nc)
    real(kind=8) :: work(nc, nno*nc), rg0(nno*nc, nno*nc)
!   pour la thermique
    real(kind=8) :: temm, em, num, f, df
    blas_int :: b_incx, b_n
!
! --------------------------------------------------------------------------------------------------
!
    xls2 = xl/2.d0
    dimklv = 2*nc*(2*nc+1)/2
!
!   booléens pratiques
    matric = option .eq. 'FULL_MECA' .or. option .eq. 'RIGI_MECA_TANG'
    vecteu = option .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
!
    fl(1:nno*nc) = 0.0d0
!
    hoel(:, :) = 0.0d0
    hota(:, :) = 0.0d0
    d1b(:, :) = 0.0d0
    work(:, :) = 0.0d0
    rg0(:, :) = 0.0d0
    fg(:) = 0.0d0
!
!   Température
    call verifm(fami, npg, 1, '-', icodma, &
                epsthf, iret)
    call verifm(fami, npg, 1, 'T', icodma, &
                epsthd, iret)
    itemp = 0
    if (iret .eq. 0) itemp = 1
    nomres(1) = 'E'
    nomres(2) = 'NU'
!   Thermique à T+
    call moytem(fami, npg, 1, '+', temp, &
                iret)
    call rcvalb(fami, 1, 1, '+', icodma, &
                ' ', 'ELAS', 1, 'TEMP', [temp], &
                2, nomres, valres, codres, 1)
    e = valres(1)
    nu = valres(2)
    g = e/(2.d0*(1.d0+nu))
!   thermique à T-
    call moytem(fami, npg, 1, '-', temm, &
                iret)
    call rcvalb(fami, 1, 1, '-', icodma, &
                ' ', 'ELAS', 1, 'TEMP', [temm], &
                2, nomres, valres, codres, 1)
    em = valres(1)
    num = valres(2)
!
!   caractéristiques de la section
    aa = sect(1)
    xiy = sect(2)
    xiz = sect(3)
    alfay = sect(4)
    alfaz = sect(5)
    xjx = sect(8)
    xjg = sect(9)
!
!   matériau integré sur la section
    hoel(1, 1) = e*aa
    hoel(2, 2) = g*aa/alfay
    hoel(3, 3) = g*aa/alfaz
    phiy = e*xiz*12.d0*alfay/(xl*xl*g*aa)
    phiz = e*xiy*12.d0*alfaz/(xl*xl*g*aa)
    hoel(4, 4) = g*xjx
    hoel(5, 5) = e*xiy
    hoel(6, 6) = e*xiz
    hoel(7, 7) = e*xjg
!
!   boucle sur les points de gauss
    do kp = 1, npg
!       calcul de d1b ( epsi = d1b * u )
        call jsd1ff(kp, xl, phiy, phiz, d1b)
!       calcul de eps, deps et sigm (effort au pt de gauss)
!       et de dsigm = incrément d'effort élastique
        eps(:) = 0.d0
        deps(:) = 0.0d0
        sigm(:) = 0.d0
        do i = 1, nc
            do j = 1, nno*nc
                eps(i) = eps(i)+d1b(i, j)*u(j)
                deps(i) = deps(i)+d1b(i, j)*du(j)
            end do
            sigm(i) = contm(nc*(kp-1)+i)*e/em
        end do
        if ((epsthd .ne. 0.d0) .and. (itemp .ne. 0)) then
            f = epsthf
            df = epsthd
            eps(1) = eps(1)-f
            deps(1) = deps(1)-df
        end if
!       Élasticité
        do i = 1, nc
            hota(i, i) = hoel(i, i)
            sigp(i) = sigm(i)+hoel(i, i)*deps(i)
        end do
!       calcul de bt*h*b
        if (matric) then
            b_n = to_blas_int(nc*nc)
            b_incx = to_blas_int(1)
            call dscal(b_n, xls2, hota, b_incx)
            b_n = to_blas_int(nc*nc)
            b_incx = to_blas_int(1)
            call dscal(b_n, wgauss(kp), hota, b_incx)
            call utbtab('CUMU', nc, nno*nc, hota, d1b, &
                        work, rg0)
        end if
!       Les contraintes "+" et fl
        if (vecteu) then
            do i = 1, nc
                contp(nc*(kp-1)+i) = sigp(i)
            end do
            do k = 1, nno*nc
                do kk = 1, nc
                    fl(k) = fl(k)+xls2*sigp(kk)*d1b(kk, k)*wgauss(kp)
                end do
            end do
        end if
    end do
!
    if (matric) then
        call mavec(rg0, nno*nc, klv, dimklv)
    end if
!
end subroutine
