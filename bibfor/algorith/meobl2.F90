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
subroutine meobl2(eps, b, d, deltab, deltad, &
                  mult, lambda, mu, ecrob, ecrod, &
                  alpha, k1, k2, bdim, dsidep)
!
!
    implicit none
!
#include "asterfort/ceobfb.h"
#include "asterfort/ceobfd.h"
#include "asterfort/dfbdb.h"
#include "asterfort/dfbde.h"
#include "asterfort/dfddd.h"
#include "asterfort/dfdde.h"
#include "asterfort/dfmdf.h"
#include "asterfort/mgauss.h"
#include "asterfort/r8inir.h"
    real(kind=8) :: eps(6), b(6), d, dsidep(6, 6)
    real(kind=8) :: deltab(6), deltad, mult
    real(kind=8) :: lambda, mu, alpha, k1, k2, ecrob, ecrod
    integer(kind=8) :: bdim
!
!--CALCUL DE LA MATRICE TANGENTE POUR LA LOI ENDO_ORTHO_BETON
!-------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, iret
    real(kind=8) :: rac2, nofbm, un, det, deux, fb(6)
    real(kind=8) :: treps, fd, dfmf(3, 3), tdfbdb(6, 6), tdfbde(6, 6)
    real(kind=8) :: tdfdde(6), tdfddd, interd(3), intert(6), interg(3)
    real(kind=8) :: psi(3, 6), ksi(3, 3), iksi(3, 3), matb(3, 6), matd(6)
    real(kind=8) :: fbs(3), deltas(3)
    real(kind=8) :: fbsm(6), sdfbdb(3, 3), sdfbde(3, 6)
    real(kind=8) :: coupl, dcrit(6)
!
    deux = 2.d0
    rac2 = sqrt(deux)
    un = 1.d0
!
!-------------------------------------------------------
!-------------------------------------------------------
!----CALCUL DE FB: FORCE THERMO ASSOCIEE A
!-------------------ENDOMMAGEMENT ANISOTROPE DE TRACTION
!
    call ceobfb(b, eps, lambda, mu, ecrob, &
                bdim, fb, nofbm, fbsm)
!
    fbs(1) = fb(1)
    fbs(2) = fb(2)
    fbs(3) = fb(4)
!
    deltas(1) = deltab(1)
    deltas(2) = deltab(2)
    deltas(3) = deltab(4)
!
!----CALCUL DE FD: PARTIE POSITIVE DE LA FORCE THERMO ASSOCIEE A
!-------------------ENDOMMAGEMENT ISOTROPE DE COMPRESSION
!
    call ceobfd(d, eps, lambda, mu, ecrod, &
                fd)
!
!---CALCUL DE DERIVEES UTILES----------------------------------
!
    call dfmdf(3, fbs, dfmf)
!
!----CALCUL DE LA DERIVEE DU SEUIL---------------------
!
    treps = eps(1)+eps(2)+eps(3)
    if (treps .gt. 0.d0) then
        treps = 0.d0
    end if
    dcrit(1) = -k1*(-treps/k2/(un+(-treps/k2)**deux)&
     &           +atan2(-treps/k2, un))
    dcrit(2) = -k1*(-treps/k2/(un+(-treps/k2)**deux)&
     &           +atan2(-treps/k2, un))
    dcrit(3) = -k1*(-treps/k2/(un+(-treps/k2)**deux)&
     &           +atan2(-treps/k2, un))
    dcrit(4) = 0.d0
    dcrit(5) = 0.d0
    dcrit(6) = 0.d0
!
    call dfbdb(3, b, eps, deux*mu, lambda, &
               ecrob, tdfbdb)
    call dfbde(3, b, eps, deux*mu, lambda, &
               tdfbde)
!
    sdfbdb(1, 1) = tdfbdb(1, 1)
    sdfbdb(1, 2) = tdfbdb(1, 2)
    sdfbdb(1, 3) = tdfbdb(1, 4)
    sdfbdb(2, 1) = tdfbdb(2, 1)
    sdfbdb(2, 2) = tdfbdb(2, 2)
    sdfbdb(2, 3) = tdfbdb(2, 4)
    sdfbdb(3, 1) = tdfbdb(4, 1)
    sdfbdb(3, 2) = tdfbdb(4, 2)
    sdfbdb(3, 3) = tdfbdb(4, 4)
!
    do i = 1, 6
        sdfbde(1, i) = tdfbde(1, i)
        sdfbde(2, i) = tdfbde(2, i)
        sdfbde(3, i) = tdfbde(4, i)
    end do
!
    fbsm(3) = rac2*fbsm(3)
    deltas(3) = deltas(3)*rac2
!
    call dfdde(eps, d, 3, lambda, mu, &
               tdfdde)
    call dfddd(eps, d, 3, lambda, mu, &
               ecrod, tdfddd)
!
    nofbm = fbsm(1)**2+fbsm(2)**2+fbsm(3)**2
!
    coupl = sqrt(alpha*nofbm+(un-alpha)*fd**deux)
    call r8inir(36, 0.d0, dsidep, 1)
!
    if ((fd .ne. 0.d0) .and. (nofbm .ne. 0.d0)) then
!
!---CALCUL DE DBDE ET DDDE-------------------------------------
!
!---CALCUL DE KSI ET PSI
!
        call r8inir(3, 0.d0, interd, 1)
        call r8inir(3, 0.d0, interg, 1)
        call r8inir(6, 0.d0, intert, 1)
        call r8inir(18, 0.d0, psi, 1)
        call r8inir(9, 0.d0, ksi, 1)
!
        do i = 1, 6
            intert(i) = (un-alpha)*fd*tdfdde(i)-coupl*dcrit(i)
            do j = 1, 3
                do k = 1, 3
                    intert(i) = intert(i)+alpha*fbsm(k)*dfmf(k, j)* &
                                sdfbde(j, i)
                end do
            end do
        end do
!
        do i = 1, 3
            interg(i) = deltas(i)/fd-alpha*fbsm(i)/(un-alpha)/fd/tdfddd
            do j = 1, 3
                do k = 1, 3
                    ksi(i, j) = ksi(i, j)+alpha*deltad*dfmf(i, k)*sdfbdb(k, &
                                                                         j)
                    interd(i) = interd(i)+alpha*fbsm(k)*dfmf(k, j)* &
                                sdfbdb(j, i)
                end do
            end do
            do j = 1, 6
                do k = 1, 3
                    psi(i, j) = psi(i, j)-alpha*deltad*dfmf(i, k)*sdfbde(k, &
                                                                         j)
                end do
            end do
        end do
!
        do i = 1, 3
            ksi(i, i) = ksi(i, i)-(un-alpha)*fd
        end do
!
        do i = 1, 3
            do j = 1, 3
                ksi(i, j) = ksi(i, j)+interg(i)*interd(j)
            end do
            do j = 1, 6
                psi(i, j) = psi(i, j)-interg(i)*intert(j)+(un-alpha)* &
                            deltas(i)*tdfdde(j)
            end do
        end do
!
        call r8inir(9, 0.d0, iksi, 1)
        do i = 1, 3
            iksi(i, i) = 1.d0
        end do
!
        call mgauss('NFVP', ksi, iksi, 3, 3, &
                    3, det, iret)
!
!-- ! ksi n est plus disponible
!
        call r8inir(18, 0.d0, matb, 1)
        call r8inir(6, 0.d0, matd, 1)
!
        do i = 1, 6
            matd(i) = -intert(i)/(un-alpha)/fd/tdfddd
            do j = 1, 3
                do k = 1, 3
                    matb(j, i) = matb(j, i)+iksi(j, k)*psi(k, i)
                    matd(i) = matd(i)-interd(j)*iksi(j, k)*psi(k, i) &
                              /(un-alpha)/fd/tdfddd
                end do
!            WRITE(6,*) 'MB(',J,',',I,')=',MATB(J,I),';'
            end do
        end do
!
        do i = 1, 6
            do j = 1, 6
                dsidep(i, j) = -tdfdde(i)*matd(j)
!         WRITE(6,*) 'DID(',I,',',J,')=', DSIDEP(I,J),';'
                do k = 1, 3
                    dsidep(i, j) = dsidep(i, j)-sdfbde(k, i)*matb(k, j)
                end do
            end do
        end do
!
    else if ((fd .eq. 0.d0) .and. (nofbm .ne. 0.d0)) then
!
        call r8inir(9, 0.d0, ksi, 1)
        call r8inir(18, 0.d0, psi, 1)
!
        do i = 1, 3
            do j = 1, 3
                ksi(i, j) = -fbsm(i)*fbsm(j)/nofbm
                do k = 1, 3
                    ksi(i, j) = ksi(i, j)-alpha*mult*dfmf(i, k)*sdfbdb(k, j)
                end do
            end do
            do j = 1, 6
                psi(i, j) = psi(i, j)-fbsm(i)*alpha*mult/coupl*dcrit(j)
                do k = 1, 3
                    psi(i, j) = psi(i, j)+alpha*mult*dfmf(i, k)*sdfbde(k, j)
                end do
            end do
        end do
!
        do i = 1, 3
            ksi(i, i) = ksi(i, i)+1
        end do
!
        call r8inir(9, 0.d0, iksi, 1)
        do i = 1, 3
            iksi(i, i) = 1.d0
        end do
!
        call mgauss('NFVP', ksi, iksi, 3, 3, &
                    3, det, iret)
!
        call r8inir(18, 0.d0, matb, 1)
!
        do i = 1, 3
            do j = 1, 6
                do k = 1, 3
                    matb(i, j) = matb(i, j)+iksi(i, k)*psi(k, j)
                end do
            end do
        end do
!
        do i = 1, 6
            do j = 1, 6
                do k = 1, 3
                    dsidep(i, j) = dsidep(i, j)-sdfbde(k, i)*matb(k, j)
                end do
            end do
        end do
!
    else if ((fd .ne. 0.d0) .and. (nofbm .eq. 0.d0)) then
!
        do i = 1, 6
            do j = 1, 6
                dsidep(i, j) = -tdfdde(i)*(-tdfdde(j)+coupl/(un-alpha) &
                                           *dcrit(j)/fd)/tdfddd
            end do
        end do
!
    end if
!
end subroutine
