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
subroutine lceob2(intmax, tole, eps, bm, dm, &
                  lambda, mu, alpha, ecrob, ecrod, &
                  seuil, bdim, b, d, mult, &
                  elas, dbloq, iret)
!
!
    implicit none
#include "asterf_types.h"
#include "asterfort/ceobfb.h"
#include "asterfort/ceobfd.h"
#include "asterfort/dfbdb.h"
#include "asterfort/dfmdf.h"
#include "asterfort/mgauss.h"
#include "asterfort/r8inir.h"
    real(kind=8) :: eps(6)
    real(kind=8) :: bm(6), dm, b(6), d, mult
    real(kind=8) :: lambda, mu, alpha, seuil, ecrob, ecrod
    real(kind=8) :: tole
!
    integer(kind=8) :: intmax, iret, bdim
!
    aster_logical :: elas, dbloq
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT DU MODELE D'ENDOMMAGEMENT ANISOTROPE
!     ROUTINE DE RESOLUTION DU SYSTEME NON LINEAIRE
!     ALGORITHME DE NEWTON
!
!
!
!  IN INTMAX  : NBRE D'ITERATION MAX POUR LE NEWTON LOCAL
!  IN TOLE    : RESIDU TOLERE POUR LE NEWTON LOCAL
!  IN  BDIM   : DIMENSION DE L'ESPACE
!  IN  CRIT   : CRITERES DE CONVERGENCE LOCAUX
!  IN  EPSM   : DEFORMATION EN T- REPERE GLOBAL
!  IN  DEPS   : INCREMENT DE DEFORMATION
!  IN  BM DM  : VARIABLES INTERNES EN T-
!  IN LAMBDA  : /
!  IN MU      : / COEFFICIENTS DE LAME
!  IN  ALPHA  : /
!  IN  ECROB  : /
!  IN  ECROD  : / PARAMETRES DU MODELE
!  IN  SEUIL  : SEUIL DU CRITERE D'ENDOMMAGEMENT
!  IN  BDIM   : DIMENSION DE L ESPACE
!
! OUT  B D    : VARIABLES INTERNES EN T+
! OUT MULT    : MULTIPLICATEUR PLASTIQUE DU PRINCIPE DE NORMALITE
! OUT ELAS    : ELASTIQUE OU DISSIPATION?
! OUT DBLOQ   : BLOQUAGE DE L'ENDOMMAGEMENT DE COMPRESSION
! OUT IRET    : CODE RETOUR
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, compte, iret1
!
    real(kind=8) :: bs(3), bms(3), dbs(3)
    real(kind=8) :: fb(6), fd, dd, ddg, resb(3)
    real(kind=8) :: rac2, rtemp2, rtemp3, delta1(3), delta2
    real(kind=8) :: tolc, det, tata, normrb, rtemp, crit
    real(kind=8) :: mte1(3, 3), mte2(6, 6), mte2s(3, 3)
    real(kind=8) :: fbs(3), fbsm(6)
    real(kind=8) :: ksi(3, 3), iksi(3, 3), toti(3, 3), ide(3, 3)
    real(kind=8) :: teme(3, 3), coupl
    real(kind=8) :: resd, dfddd, psi
    real(kind=8) :: inter1, inter2, inter3, inter4
    real(kind=8) :: deux, un
!
    deux = 2.d0
    rac2 = sqrt(deux)
    tolc = seuil*tole
    un = 1.d0
    compte = 0
    mult = 0.d0
!
    do i = 1, 6
        b(i) = bm(i)
    end do
    d = dm
!
!-------------------------------------------------------
!-------------------------------------------------------
!----CALCUL DE FB: FORCE THERMO ASSOCIEE A
!-------------------ENDOMMAGEMENT ANISOTROPE DE TRACTION
!
    call ceobfb(b, eps, lambda, mu, ecrob, &
                bdim, fb, rtemp, fbsm)
!
    fbs(1) = fb(1)
    fbs(2) = fb(2)
    fbs(3) = fb(4)
!
    bs(1) = b(1)
    bs(2) = b(2)
    bs(3) = b(4)
!
    bms(1) = bm(1)
    bms(2) = bm(2)
    bms(3) = bm(4)
!
!----CALCUL DE FD: PARTIE POSITIVE DE LA FORCE THERMO ASSOCIEE A
!-------------------ENDOMMAGEMENT ISOTROPE DE COMPRESSION
!
    if (dbloq) then
        fd = 0.d0
    else
        call ceobfd(d, eps, lambda, mu, ecrod, &
                    fd)
    end if
!
!----CALCUL DU CRITERE-------------------------------------
    coupl = sqrt(alpha*rtemp+(un-alpha)*fd**2)
    crit = coupl-seuil
!
!
    elas = .false.
!
    if (crit .le. tolc) then
        elas = .true.
        goto 999
!
    else
!
        do i = 1, 3
            resb(i) = -bs(i)+bms(i)+alpha*mult*fbsm(i)
        end do
        resd = -d+dm+(un-alpha)*mult*fd
        resb(3) = rac2*resb(3)
!
        tata = 0.d0
        do i = 1, 3
            tata = tata+resb(i)*resb(i)
        end do
!
        normrb = sqrt(tata)
!
        ddg = 0.d0
!
!--------------------------------------------------------
!--BOUCLE DU NEWTON SUR LES VARIABLES INTERNES-----------
!--------------------------------------------------------
!
!
38      continue
        if (((crit .gt. tolc) .or. (normrb .gt. tole) .or. (abs(resd) .gt. tole))) then
            if ((compte .lt. intmax) .and. (coupl .ne. 0.d0)) then
! Rajout du test sur COUPL (fiche 15020) : lorsque c'est le cas,
! la derivee du residu est une matrice singuliere et le systeme ne
! peut etre resolu. On sort pour enclencher la decoupe du pas de temps
!
                call dfmdf(3, fbs, mte1)
!
                call dfbdb(3, b, eps, deux*mu, lambda, &
                           ecrob, mte2)
!
                mte2s(1, 1) = mte2(1, 1)
                mte2s(1, 2) = mte2(1, 2)
                mte2s(1, 3) = mte2(1, 4)
                mte2s(2, 1) = mte2(2, 1)
                mte2s(2, 2) = mte2(2, 2)
                mte2s(2, 3) = mte2(2, 4)
                mte2s(3, 1) = mte2(4, 1)
                mte2s(3, 2) = mte2(4, 2)
                mte2s(3, 3) = mte2(4, 4)
!
                dfddd = 0.d0
!
                if ((.not. dbloq) .and. (fd .ne. 0.d0)) then
                    dfddd = -(fd+deux*ecrod)/(un-d)
                end if
!
                call r8inir(9, 0.d0, ksi, 1)
                do i = 1, 3
                    do j = 1, 3
                        do k = 1, 3
                            ksi(i, j) = ksi(i, j)-mult*alpha*mte1(i, k)* &
                                        mte2s(k, j)
                        end do
                    end do
                end do
!
                do i = 1, 3
                    ksi(i, i) = ksi(i, i)+un
                end do
!
                do i = 1, 3
                    do j = 1, 3
                        toti(i, j) = ksi(i, j)
                    end do
                end do
!
                do i = 1, 3
                    do j = 1, 3
                        if (i .eq. j) then
                            ide(i, j) = 1.d0
                        else
                            ide(i, j) = 0.d0
                        end if
                    end do
                end do
                call r8inir(9, 0.d0, teme, 1)
                do i = 1, 3
                    do j = 1, 3
                        teme(i, j) = ide(i, j)
                    end do
                end do
                call mgauss('NFVP', toti, teme, 3, 3, &
                            3, det, iret1)
                call r8inir(9, 0.d0, iksi, 1)
                do i = 1, 3
                    do j = 1, 3
                        iksi(i, j) = teme(i, j)
                    end do
                end do
!
                psi = un-mult*(un-alpha)*dfddd
!
                call r8inir(3, 0.d0, delta1, 1)
                do i = 1, 3
                    do j = 1, 3
                        if (j .eq. 3) then
                            rtemp2 = rac2
                        else
                            rtemp2 = 1.d0
                        end if
                        delta1(i) = delta1(i)+alpha/coupl*rtemp2*fbsm(j) &
                                    *mte2s(j, i)
                    end do
                end do
!
                delta2 = (un-alpha)/coupl*fd*dfddd
!
                inter1 = 0.d0
                inter3 = 0.d0
                do i = 1, 3
                    do j = 1, 3
                        if (j .eq. 3) then
                            rtemp2 = rac2
                        else
                            rtemp2 = 1.d0
                        end if
                        inter1 = inter1+delta1(i)*iksi(i, j)*resb(j)
                        inter3 = inter3+alpha*rtemp2*delta1(i)*iksi(i, &
                                                                    j)*fbsm(j)
                    end do
                end do
!
                inter2 = delta2/psi*resd
                inter4 = delta2/psi*(un-alpha)*fd
!
                ddg = -(crit+inter1+inter2)/(inter3+inter4)
!
                dd = resd/psi+ddg*(un-alpha)*fd/psi
                call r8inir(3, 0.d0, dbs, 1)
                do i = 1, 3
                    do j = 1, 3
                        if (i .eq. 3) then
                            rtemp2 = 1/rac2
                        else
                            rtemp2 = 1.d0
                        end if
                        if (j .eq. 3) then
                            rtemp3 = rac2
                        else
                            rtemp3 = 1.d0
                        end if
                        dbs(i) = dbs(i)+rtemp2*iksi(i, j)*(resb(j)+ddg* &
                                                           alpha*fbsm(j)*rtemp3)
                    end do
                end do
!
                do i = 1, 3
                    bs(i) = bs(i)+dbs(i)
                end do
                d = d+dd
                compte = compte+1
                mult = mult+ddg
!
!----CALCUL DE FB DANS NEWTON---------------------------
!
                call r8inir(6, 0.d0, b, 1)
!
                b(1) = bs(1)
                b(2) = bs(2)
                b(4) = bs(3)
!
                call ceobfb(b, eps, lambda, mu, ecrob, &
                            bdim, fb, rtemp, fbsm)
!
                fbs(1) = fb(1)
                fbs(2) = fb(2)
                fbs(3) = fb(4)
!
!----CALCUL DE FD: PARTIE POSITIVE DE LA FORCE THERMO ASSOCIEE A
!-------------------ENDOMMAGEMENT ISOTROPE DE COMPRESSION
!
                if (dbloq) then
                    fd = 0.d0
                else
                    call ceobfd(d, eps, lambda, mu, ecrod, &
                                fd)
                end if
!
!----CALCUL DU CRITERE-------------------------------------
                coupl = sqrt(alpha*rtemp+(un-alpha)*fd**deux)
                crit = coupl-seuil
!
                do i = 1, 3
                    resb(i) = -bs(i)+bms(i)+alpha*mult*fbsm(i)
                end do
                resd = -d+dm+(un-alpha)*mult*fd
                resb(3) = rac2*resb(3)
!
                tata = 0.d0
                do i = 1, 3
                    tata = tata+resb(i)*resb(i)
                end do
!
                normrb = sqrt(tata)
!
                goto 38
            else
                iret = 1
                goto 999
            end if
        end if
!
    end if
999 continue
!
!
end subroutine
