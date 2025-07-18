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
subroutine lceob3(intmax, tole, eps, bm, dm, &
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
    real(kind=8) :: fb(6), db(6), fd, dd, fbm(6), resb(6)
    real(kind=8) :: rac2, un, deux
    real(kind=8) :: rtemp2, rtemp3, delta1(6), delta2
    real(kind=8) :: ddg, tolc, det, tata, normrb, rtemp, crit
    real(kind=8) :: mte1(6, 6), mte2(6, 6)
    real(kind=8) :: ksi(6, 6), iksi(6, 6), toti(6, 6), ide(6, 6)
    real(kind=8) :: teme(6, 6), coupl
    real(kind=8) :: resd, dfddd, psi
    real(kind=8) :: inter1, inter2, inter3, inter4
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
!
    d = dm
!
!-------------------------------------------------------
!-------------------------------------------------------
!----CALCUL DE FB: FORCE THERMO ASSOCIEE A
!-------------------ENDOMMAGEMENT ANISOTROPE DE TRACTION
!
    call ceobfb(b, eps, lambda, mu, ecrob, &
                bdim, fb, rtemp, fbm)
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
!
    coupl = sqrt(alpha*rtemp+(1-alpha)*fd**2)
    crit = coupl-seuil
!
    elas = .false.
!
    if (crit .le. tolc) then
        elas = .true.
        goto 999
!
    else
!
        do i = 1, 6
            resb(i) = -b(i)+bm(i)+alpha*mult*fbm(i)
        end do
        resd = -d+dm+(1-alpha)*mult*fd
        do i = 4, 6
            resb(i) = rac2*resb(i)
        end do
!
        tata = 0.d0
        do i = 1, 6
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
38      continue
        if (((crit .gt. tolc) .or. (normrb .gt. tole) .or. (abs(resd) .gt. tole))) then
            if ((compte .lt. intmax) .and. (coupl .ne. 0.d0)) then
! Rajout du test sur COUPL (fiche 15020) : lorsque c'est le cas,
! la derivee du residu est une matrice singuliere et le systeme ne
! peut etre resolu. On sort pour enclencher la decoupe du pas de temps
!
                call dfmdf(6, fb, mte1)
                call dfbdb(3, b, eps, deux*mu, lambda, &
                           ecrob, mte2)
!
                dfddd = 0.d0
!
                if ((.not. dbloq) .and. (fd .ne. 0.d0)) then
                    dfddd = -(fd+deux*ecrod)/(un-d)
                end if
!
                call r8inir(36, 0.d0, ksi, 1)
                do i = 1, 6
                    do j = 1, 6
                        do k = 1, 6
                            ksi(i, j) = ksi(i, j)-mult*alpha*mte1(i, k)* &
                                        mte2(k, j)
                        end do
                    end do
                end do
!
                do i = 1, 6
                    ksi(i, i) = ksi(i, i)+1
                end do
!
                do i = 1, 6
                    do j = 1, 6
                        toti(i, j) = ksi(i, j)
                    end do
                end do
                do i = 1, 6
                    do j = 1, 6
                        if (i .eq. j) then
                            ide(i, j) = 1.d0
                        else
                            ide(i, j) = 0.d0
                        end if
                    end do
                end do
                call r8inir(36, 0.d0, teme, 1)
                do i = 1, 6
                    do j = 1, 6
                        teme(i, j) = ide(i, j)
                    end do
                end do
                call mgauss('NFVP', toti, teme, 6, 6, &
                            6, det, iret1)
                call r8inir(36, 0.d0, iksi, 1)
                do i = 1, 6
                    do j = 1, 6
                        iksi(i, j) = teme(i, j)
                    end do
                end do
!
                psi = 1-mult*(1-alpha)*dfddd
!
                call r8inir(6, 0.d0, delta1, 1)
                do i = 1, 6
                    do j = 1, 6
                        if (j .ge. 4) then
                            rtemp2 = rac2
                        else
                            rtemp2 = 1.d0
                        end if
                        delta1(i) = delta1(i)+alpha/coupl*rtemp2*fbm(j)* &
                                    mte2(j, i)
                    end do
                end do
!
                delta2 = (1-alpha)/coupl*fd*dfddd
!
                inter1 = 0.d0
                inter3 = 0.d0
                do i = 1, 6
                    do j = 1, 6
                        if (j .ge. 4) then
                            rtemp2 = rac2
                        else
                            rtemp2 = 1.d0
                        end if
                        inter1 = inter1+delta1(i)*iksi(i, j)*resb(j)
                        inter3 = inter3+alpha*rtemp2*delta1(i)*iksi(i, &
                                                                    j)*fbm(j)
                    end do
                end do
!
                inter2 = delta2/psi*resd
                inter4 = delta2/psi*(1-alpha)*fd
!
                ddg = -(crit+inter1+inter2)/(inter3+inter4)
!
                call r8inir(6, 0.d0, db, 1)
!
                dd = resd/psi+ddg*(1-alpha)*fd/psi
                do i = 1, 6
                    do j = 1, 6
                        if (i .ge. 4) then
                            rtemp2 = 1/rac2
                        else
                            rtemp2 = 1.d0
                        end if
                        if (j .ge. 4) then
                            rtemp3 = rac2
                        else
                            rtemp3 = 1.d0
                        end if
                        db(i) = db(i)+rtemp2*iksi(i, j)*(resb(j)+ddg* &
                                                         alpha*fbm(j)*rtemp3)
                    end do
                end do
!
                do i = 1, 6
                    b(i) = b(i)+db(i)
                end do
                d = d+dd
!
                compte = compte+1
                mult = mult+ddg
!
!----CALCUL DE FB DANS NEWTON---------------------------
!
                call ceobfb(b, eps, lambda, mu, ecrob, &
                            bdim, fb, rtemp, fbm)
!
!----CALCUL DE FD DANS NEWTON----------------------------
                if (dbloq) then
                    fd = 0.d0
                else
                    call ceobfd(d, eps, lambda, mu, ecrod, &
                                fd)
                end if
!
!----CALCUL DU CRITERE-------------------------------------
                coupl = sqrt(alpha*rtemp+(1-alpha)*fd**2)
                crit = coupl-seuil
!
                do i = 1, 6
                    resb(i) = -b(i)+bm(i)+alpha*mult*fbm(i)
                end do
                resd = -d+dm+(1-alpha)*mult*fd
                do i = 4, 6
                    resb(i) = rac2*resb(i)
                end do
!
                tata = 0.d0
                do i = 1, 6
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
end subroutine
