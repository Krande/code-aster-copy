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
subroutine te0327(option, nomte)
    implicit none
!....................................................................
!   CALCUL DES TERMES ELEMENTAIRES D'AMORTISSEMENT AJOUTE
!     OPTION : AMOR_AJOU
!....................................................................
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/subacv.h"
#include "asterfort/sumetr.h"
!
    character(len=16) :: nomte, option
    real(kind=8) :: sx(9, 9), sy(9, 9), sz(9, 9), jac(9)
    real(kind=8) :: nx(9), ny(9), nz(9), norm(3, 9), acc(3, 9)
    real(kind=8) :: flufn(9), acloc(3, 8), cova(3, 3), cnva(3, 3)
    real(kind=8) :: a(2, 2), metr(2, 2), e1(3, 9), e2(3, 9), jc
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom, imattt
    integer(kind=8) :: ndim, nno, ipg, npg1
    integer(kind=8) :: idec, jdec, kdec, ldec, nnos, jgano
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacce, idim, idir, ii, ij, ino
    integer(kind=8) :: j, jdir, jj, jno, k
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    call jevech('PACCELR', 'L', iacce)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATTTR', 'E', imattt)
!
    do i = 1, nno
        acloc(1, i) = 0.0d0
        acloc(2, i) = 0.0d0
        acloc(3, i) = 0.0d0
    end do
    k = 0
    do i = 1, nno
        do idim = 1, 3
            k = k+1
            acloc(idim, i) = zr(iacce+k-1)
        end do
    end do
!
    do ipg = 1, npg1
        acc(1, ipg) = 0.0d0
        acc(2, ipg) = 0.0d0
        acc(3, ipg) = 0.0d0
    end do
!
    do ipg = 1, npg1
        ldec = (ipg-1)*nno
!
        do i = 1, nno
            acc(1, ipg) = acc(1, ipg)+acloc(1, i)*zr(ivf+ldec+i-1)
            acc(2, ipg) = acc(2, ipg)+acloc(2, i)*zr(ivf+ldec+i-1)
            acc(3, ipg) = acc(3, ipg)+acloc(3, i)*zr(ivf+ldec+i-1)
        end do
    end do
!
! --- CALCUL DES PRODUITS VECTORIELS OMI X OMJ
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
!
! --- BOUCLE SUR LES POINTS DE GAUSS
!
    do ipg = 1, npg1
!
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
        nx(ipg) = 0.0d0
        ny(ipg) = 0.0d0
        nz(ipg) = 0.0d0
!
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
!
                nx(ipg) = nx(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                ny(ipg) = ny(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                nz(ipg) = nz(ipg)+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
!
            end do
        end do
!
! ------ CALCUL DU JACOBIEN AU POINT DE GAUSS IPG
!
        jac(ipg) = sqrt(nx(ipg)*nx(ipg)+ny(ipg)*ny(ipg)+nz(ipg)*nz(ipg))
!
! ------ CALCUL DE LA NORMALE UNITAIRE
!
        norm(1, ipg) = nx(ipg)/jac(ipg)
        norm(2, ipg) = ny(ipg)/jac(ipg)
        norm(3, ipg) = nz(ipg)/jac(ipg)
    end do
!
! --- CALCUL DES VECTEURS E1, E2 TANGENTS A L'ELEMENT NON NORMALISES
!
    do ipg = 1, npg1
        kdec = (ipg-1)*nno*ndim
        e1(1, ipg) = 0.0d0
        e1(2, ipg) = 0.0d0
        e1(3, ipg) = 0.0d0
!
        e2(1, ipg) = 0.0d0
        e2(2, ipg) = 0.0d0
        e2(3, ipg) = 0.0d0
!
        do j = 1, nno
            idec = (j-1)*ndim
!
            e1(1, ipg) = e1(1, ipg)+zr(igeom+3*(j-1)-1+1)*zr(idfdx+ &
                                                             kdec+idec)
            e1(2, ipg) = e1(2, ipg)+zr(igeom+3*(j-1)-1+2)*zr(idfdx+ &
                                                             kdec+idec)
            e1(3, ipg) = e1(3, ipg)+zr(igeom+3*(j-1)-1+3)*zr(idfdx+ &
                                                             kdec+idec)
!
            e2(1, ipg) = e2(1, ipg)+zr(igeom+3*(j-1)-1+1)*zr(idfdy+ &
                                                             kdec+idec)
            e2(2, ipg) = e2(2, ipg)+zr(igeom+3*(j-1)-1+2)*zr(idfdy+ &
                                                             kdec+idec)
            e2(3, ipg) = e2(3, ipg)+zr(igeom+3*(j-1)-1+3)*zr(idfdy+ &
                                                             kdec+idec)
        end do
!
! ------ CONSTITUTION DE LA BASE COVARIANTE
!
        do i = 1, 3
            cova(i, 1) = e1(i, ipg)
            cova(i, 2) = e2(i, ipg)
        end do
!
! ------ ON CALCULE LE TENSEUR METRIQUE
!
        call sumetr(cova, metr, jc)
!
! ------ CALCUL DE LA BASE CONTRAVARIANTE
!
        call subacv(cova, metr, jc, cnva, a)
!
! ------ CALCUL DU FLUX FLUIDE NORMAL AUX POINTS DE GAUSS
!
        flufn(ipg) = acc(1, ipg)*norm(1, ipg)+acc(2, ipg)*norm(2, ipg)+acc(3, ipg)*norm(3, ipg)
!
! CALCUL DE LA MATRICE DES PRODUITS SCALAIRES DES GRADIENTS SURFACIQUES
! DES FONCTIONS DE FORMES
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, i
                jdec = (j-1)*ndim
                ij = (i-1)*i/2+j
                do ii = 1, 2
                    idir = ii-1
                    do jj = 1, 2
                        jdir = jj-1
!
                        zr(imattt+ij-1) = zr(imattt+ij-1)+jac(ipg)*zr(ipoids+ipg-1)*flufn(ipg)&
                                         & *zr(idfdx+idir+kdec+idec)*zr(idfdx+jdir+kdec+jdec) &
                                         &*a(ii, jj)
                    end do
                end do
            end do
        end do
    end do
!
end subroutine
