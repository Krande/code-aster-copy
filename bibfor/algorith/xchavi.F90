! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine xchavi(actpoi, jbasc, jffis, jfon, jvit, &
                  jbeta, ndim, nfonn, sifval)
!
! person_in_charge: patrick.massin at edf.fr
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "blas/ddot.h"
!
! Propagation XFEM avec éléments cohésifs
! Calculer l'avancée qui produit le nouveau front
!  à partir de l'ancien
!
! In actpoi => Position du premier point du morceau
!              de front à traiter
! In jbasc => base covariante au front de fissure (ancien)
! In jfiss => coordonnées points ancien front
! In jfon => coordonnées points nouveau front
! Out jvit => avancée qui produit le nouveau front
! In ndim => dimension
! In nfonn => nombre de points nouveau front
! In sifval => nombre de points ancien front
!
    integer :: actpoi
    real(kind=8) :: b(3), beta1, ci(3), cosb
    real(kind=8) :: dir(3)
    integer :: i
    real(kind=8) :: poitot, poiav, poids
    integer :: ipt, j, jbasc, jbeta, jffis, jfon, jvit
    real(kind=8) :: loncar, m(3), mi(3), mtast, pi
    integer :: ndim, nfonn
    real(kind=8) :: n(3), t(3)
    integer :: sifval, nbptfo, ibid
    real(kind=8) :: sinb, tast(3), vecv(3)
    real(kind=8) :: vitn, vnor, vpnt, linf, lprop, lcalc
    aster_logical :: linter
    blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------
!
!
!       LONGUEUR CARACTERISTIQUE MAILLAGE
!       SERT POUR LA LONGUEUR D INFLUENCE
    loncar = (zr(jffis-1+4*(actpoi+sifval-1)+4)-zr(jffis-1+4*actpoi+4))/sifval
    loncar = abs(loncar)
    pi = r8pi()
    call getvis(' ', 'NB_POINT_FOND', scal=nbptfo, nbret=ibid)
    linf = (zr(jffis-1+4*(actpoi+sifval-1)+4)-zr(jffis-1+4*actpoi+4))/nbptfo
    linf = abs(linf)
!
    do i = 1, sifval
!
!       RECUP POINT ET BASE ANCIEN FRONT
        do j = 1, ndim
            m(j) = zr(jffis-1+4*(i-1)+j)
            n(j) = zr(jbasc-1+2*ndim*(i-1)+j)
            t(j) = zr(jbasc-1+2*ndim*(i-1)+ndim+j)
        end do
!
!       ON REORTHOGONALISE (cf XPRVIT)
        call normev(n, mtast)
        ASSERT(mtast .gt. 0.d0)
        call provec(n, t, b)
        call normev(b, mtast)
        ASSERT(mtast .gt. 0.d0)
        call provec(b, n, tast)
        call normev(tast, mtast)
        ASSERT(mtast .gt. 0.d0)
        do j = 1, ndim
            t(j) = tast(j)
        end do
!
!       LONGUEUR INFLUENCE EN DUR
        lcalc = (1.d-5)*loncar
        lprop = (1.d-3)*loncar
!
!       INITIALISATIONS
        poitot = 0.d0
        vpnt = 0.d0
        beta1 = 0.d0
        linter = .false.
!
!       BOUCLE SUR LES POINTS DU NOUVEAU FOND
!       QUI NE SONT PAS ORDONNES
        do ipt = 1, nfonn
!
!           DISTANCE PT NOUVEAU FRONT/PT ANCIEN FRONT
            do j = 1, ndim
                ci(j) = zr(jfon-1+11*(ipt-1)+j)
                mi(j) = ci(j)-m(j)
            end do
!
!           DISTANCE PT PLAN (SIGNEE)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            vitn = ddot(b_n, mi, b_incx, b, b_incy)
!
!           CALCUL DU VECTEUR VITESSE DS LE PLAN VECV
!           ET SA DIRECTION NORMALISEE DIR
            do j = 1, ndim
                vecv(j) = mi(j)-vitn*b(j)
                dir(j) = vecv(j)
            end do
            call normev(dir, vnor)
!
!           SI DISTANCE < DISTANCE INFLUENCE
!           ON RENTRE CETTE VITESSE DANS LA MOYENNE
!           POIDS GAUSSIENS TRONQUES A [-2*Linf,+2*Linf]
            poids = 1.d0/(sqrt(2*pi)*linf)*(exp(-vitn*vitn/(2*linf*linf))-0.05)
            if (poids .gt. 0.d0) then
                linter = .true.
                poiav = poitot
                poitot = poitot+poids
                b_n = to_blas_int(3)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                if ((ddot(b_n, vecv, b_incx, t, b_incy)) .le. lcalc) then
                    vpnt = vpnt*poiav/poitot
                else
                    vpnt = vpnt*poiav/poitot+vnor*poids/poitot
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    sinb = ddot(b_n, dir, b_incx, n, b_incy)
                    b_n = to_blas_int(3)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    cosb = ddot(b_n, dir, b_incx, t, b_incy)
                    beta1 = beta1*poiav/poitot+(atan2(sinb, cosb))*poids/poitot
                end if
            end if
!
        end do
!
!       STOCKAGE VITESSE ET ANGLE
!       PROPAGATION MINIMALE LPROP
!
        if (.not. linter) then
            zr(jvit-1+i) = -2.d0
        else if (vpnt .gt. lprop) then
            zr(jvit-1+i) = vpnt
        else
            zr(jvit-1+i) = 0.d0
        end if
        zr(jbeta-1+i) = beta1
!
    end do
!
! --- SI PAS DE POINT TROUVE, ON PREND LA VITESSE LA PLUS PROCHE
!
    do i = 1, sifval
        if (zr(jvit-1+i) .lt. -1.d0) then
!            Il faudrait ajouter un message pour l'utilisateur ici
!            pour dire qu'on prend la vitesse la plus proche
            if (i .lt. sifval/2) then
                do j = i+1, sifval/2
                    if (zr(jvit-1+j) .ge. 0.d0) zr(jvit-1+i) = zr(jvit-1+j)
                    if (zr(jvit-1+j) .ge. 0.d0) zr(jbeta-1+i) = zr(jbeta-1+j)
                    if (zr(jvit-1+j) .ge. 0.d0) goto 500
                end do
500             continue
            else if (i .gt. sifval/2) then
                do j = i-1, sifval/2, -1
                    if (zr(jvit-1+j) .ge. 0.d0) zr(jvit-1+i) = zr(jvit-1+j)
                    if (zr(jvit-1+j) .ge. 0.d0) zr(jbeta-1+i) = zr(jbeta-1+j)
                    if (zr(jvit-1+j) .ge. 0.d0) goto 505
                end do
505             continue
            end if
        end if
    end do
end subroutine
