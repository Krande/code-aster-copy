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
subroutine meobl3(eps, b, d, deltab, deltad, &
                  mult, lambda, mu, ecrob, ecrod, &
                  alpha, k1, k2, bdim, dsidep)
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
    real(kind=8) :: rac2, nofbm, un, det, deux
    real(kind=8) :: fb(6), fbm(6)
    real(kind=8) :: treps, fd
    real(kind=8) :: dfbmdf(6, 6), tdfbdb(6, 6), tdfbde(6, 6)
    real(kind=8) :: tdfdde(6), tdfddd
    real(kind=8) :: interd(6), intert(6), interg(6)
    real(kind=8) :: psi(6, 6), ksi(6, 6), iksi(6, 6)
    real(kind=8) :: matb(6, 6), matd(6)
    real(kind=8) :: coupl, dcrit(6)
!
    un = 1.d0
    deux = 2.d0
    rac2 = sqrt(deux)
!
!-------------------------------------------------------
!-------------------------------------------------------
!----CALCUL DE FB: FORCE THERMO ASSOCIEE A
!-------------------ENDOMMAGEMENT ANISOTROPE DE TRACTION
!
    call ceobfb(b, eps, lambda, mu, ecrob, &
                bdim, fb, nofbm, fbm)
!
!----CALCUL DE FD: PARTIE POSITIVE DE LA FORCE THERMO ASSOCIEE A
!-------------------ENDOMMAGEMENT ISOTROPE DE COMPRESSION
!
    call ceobfd(d, eps, lambda, mu, ecrod, &
                fd)
!
!---CALCUL DE DERIVEES UTILES----------------------------------
!
    call dfmdf(6, fb, dfbmdf)
    call dfbdb(3, b, eps, deux*mu, lambda, &
               ecrob, tdfbdb)
    call dfbde(3, b, eps, deux*mu, lambda, &
               tdfbde)
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
    do i = 4, 6
        fbm(i) = rac2*fbm(i)
        deltab(i) = deltab(i)*rac2
    end do
!
    call dfdde(eps, d, 3, lambda, mu, &
               tdfdde)
    call dfddd(eps, d, 3, lambda, mu, &
               ecrod, tdfddd)
!
    nofbm = fbm(1)**2+fbm(2)**2+fbm(3)**2+fbm(4)**2&
     &        +fbm(5)**2+fbm(6)**2
!
    coupl = sqrt(alpha*nofbm+(un-alpha)*fd**deux)
!
    call r8inir(36, 0.d0, dsidep, 1)
!
    if ((fd .ne. 0.d0) .and. (nofbm .ne. 0.d0)) then
!
!---CALCUL DE DBDE ET DDDE-------------------------------------
!
!---CALCUL DE KSI ET PSI
!
        call r8inir(6, 0.d0, interd, 1)
        call r8inir(6, 0.d0, interg, 1)
        call r8inir(6, 0.d0, intert, 1)
        call r8inir(36, 0.d0, psi, 1)
        call r8inir(36, 0.d0, ksi, 1)
!
        do i = 1, 6
            interg(i) = deltab(i)/fd-alpha*fbm(i)/(un-alpha)/fd/tdfddd
            intert(i) = (un-alpha)*fd*tdfdde(i)-coupl*dcrit(i)
            do j = 1, 6
                do k = 1, 6
                    ksi(i, j) = ksi(i, j)+alpha*deltad*dfbmdf(i, k)*tdfbdb( &
                                k, j)
                    interd(i) = interd(i)+alpha*fbm(k)*dfbmdf(k, j)* &
                                tdfbdb(j, i)
                    psi(i, j) = psi(i, j)-alpha*deltad*dfbmdf(i, k)*tdfbde( &
                                k, j)
                    intert(i) = intert(i)+alpha*fbm(k)*dfbmdf(k, j)* &
                                tdfbde(j, i)
                end do
            end do
        end do
!
        do i = 1, 6
            ksi(i, i) = ksi(i, i)-(un-alpha)*fd
        end do
!
        do i = 1, 6
            do j = 1, 6
                ksi(i, j) = ksi(i, j)+interg(i)*interd(j)
                psi(i, j) = psi(i, j)-interg(i)*intert(j)+(un-alpha)* &
                            deltab(i)*tdfdde(j)
            end do
        end do
!
        call r8inir(36, 0.d0, iksi, 1)
        do i = 1, 6
            iksi(i, i) = 1.d0
        end do
!
        call mgauss('NFVP', ksi, iksi, 6, 6, &
                    6, det, iret)
!
!-- ! ksi n est plus disponible
!
        call r8inir(36, 0.d0, matb, 1)
        call r8inir(6, 0.d0, matd, 1)
!
        do i = 1, 6
            matd(i) = -intert(i)/(un-alpha)/fd/tdfddd
            do j = 1, 6
                do k = 1, 6
                    matb(i, j) = matb(i, j)+iksi(i, k)*psi(k, j)
                    matd(i) = matd(i)-interd(j)*iksi(j, k)*psi(k, i) &
                              /(un-alpha)/fd/tdfddd
                end do
            end do
        end do
!
        do i = 1, 6
            do j = 1, 6
                dsidep(i, j) = -tdfdde(i)*matd(j)
                do k = 1, 6
                    dsidep(i, j) = dsidep(i, j)-tdfbde(k, i)*matb(k, j)
                end do
!             WRITE(6,*) 'tang(',I,',',J,')=',DSIDEP(I,J)
            end do
        end do
!
    else if ((fd .eq. 0.d0) .and. (nofbm .ne. 0.d0)) then
!
! 567     CONTINUE
!
        call r8inir(36, 0.d0, ksi, 1)
        call r8inir(36, 0.d0, psi, 1)
!
        do i = 1, 6
            do j = 1, 6
                ksi(i, j) = -fbm(i)*fbm(j)/nofbm
                psi(i, j) = psi(i, j)-fbm(i)*alpha*mult/coupl*dcrit(j)
                do k = 1, 6
                    ksi(i, j) = ksi(i, j)-alpha*mult*dfbmdf(i, k)*tdfbdb(k, &
                                                                         j)
                    psi(i, j) = psi(i, j)+alpha*mult*dfbmdf(i, k)*tdfbde(k, &
                                                                         j)
                end do
            end do
        end do
!
        do i = 1, 6
            ksi(i, i) = ksi(i, i)+1
        end do
!
        call r8inir(36, 0.d0, iksi, 1)
        do i = 1, 6
            iksi(i, i) = 1.d0
        end do
!
        call mgauss('NFVP', ksi, iksi, 6, 6, &
                    6, det, iret)
!
        call r8inir(36, 0.d0, matb, 1)
!
        do i = 1, 6
            do j = 1, 6
                do k = 1, 6
                    matb(i, j) = matb(i, j)+iksi(i, k)*psi(k, j)
                end do
            end do
        end do
!
        do i = 1, 6
            do j = 1, 6
                do k = 1, 6
                    dsidep(i, j) = dsidep(i, j)-tdfbde(k, i)*matb(k, j)
                end do
!             WRITE(6,*) 'tang(',I,',',J,')=',DSIDEP(I,J)
            end do
        end do
!
    else if ((fd .ne. 0.d0) .and. (nofbm .eq. 0.d0)) then
!
! 568     CONTINUE
!
        do i = 1, 6
            do j = 1, 6
                dsidep(i, j) = -tdfdde(i)*(-tdfdde(j)+coupl/(un-alpha) &
                                           *dcrit(j)/fd)/tdfddd
!             WRITE(6,*) 'tang(',I,',',J,')=',DSIDEP(I,J)
            end do
        end do
!
    end if
!
end subroutine
