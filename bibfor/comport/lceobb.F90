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
subroutine lceobb(intmax, toler, epsm, deps, bm, &
                  dm, lambda, mu, alpha, ecrob, &
                  ecrod, rk, rk1, rk2, b, &
                  d, mult, elas, dbloq, iret)
!
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/diago3.h"
#include "asterfort/lceob1.h"
#include "asterfort/lceob2.h"
#include "asterfort/lceob3.h"
#include "asterfort/r8inir.h"
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: bm(6), dm, b(6), d, mult
    real(kind=8) :: lambda, mu, alpha, rk, rk1, rk2, ecrob, ecrod
    real(kind=8) :: toler
    integer(kind=8) :: intmax, iret
    aster_logical :: elas, dbloq
!
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT DU MODELE D'ENDOMMAGEMENT ANISOTROPE
!     ROUTINE DE DECOUPAGE DE L'INCREMENT DE CHARGE LORSQUE
!     L ENDOMMAGEMENT APPROCHE DE 1
!
!  IN INTMAX  : NBRE D'ITERATION MAX POUR LE NEWTON LOCAL
!  IN TOLER   : RESIDU TOLERE POUR LE NEWTON LOCAL
!  IN  NDIM    : DIMENSION DE L'ESPACE
!  IN  TYPMOD  : TYPE DE MODELISATION
!  IN  IMATE   : NATURE DU MATERIAU
!  IN  CRIT   : CRITERES DE CONVERGENCE LOCAUX
!  IN  EPSM    : DEFORMATION EN T- REPERE GLOBAL
!  IN  DEPS    : INCREMENT DE DEFORMATION
!  IN  BM DM     : VARIABLES INTERNES EN T-
!  IN LAMBDA     : /
!  IN MU        : / COEFFICIENTS DE LAME
!  IN  ALPHA    : /
!  IN  ECROB    : /
!  IN  ECROD    : /
!  IN  RK    : /
!  IN  RK1    : /
!  IN  RK2    : / PARAMETRES DU MODELE
!
! OUT  B D     : VARIABLES INTERNES EN T+
! OUT MULT     : MULTIPLICATEUR PLASTIQUE DU PRINCIPE DE NORMALITE
! OUT ELAS     : ELASTIQUE OU DISSIPATION?
! OUT DBLOQ  : BLOQUAGE DE L'ENDOMMAGEMENT DE COMPRESSION
! OUT IRET     : CODE RETOUR
! ----------------------------------------------------------------------
!
    aster_logical :: reinit, tot1, tot2, tot3
    integer(kind=8) :: i, j, k, l
    integer(kind=8) :: bdim, compte, t(3, 3)
!
    real(kind=8) :: tolb, un, deux
    real(kind=8) :: valbm(3), vecbm(3, 3), valbr(3), vecbr(3, 3)
    real(kind=8) :: valb(3), vecb(3, 3), seuil
    real(kind=8) :: bmr(6), br(6), epsr(6)
    real(kind=8) :: interm(3, 3), binter(6), epi(6)
    real(kind=8) :: epsi(6), epst(6), epsf(6), trepsm
!
    un = 1.d0
    deux = 2.d0
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
    tolb = 1.d-2
    compte = 0
!-------------------------------------------------
! -- DEFORMATIONS
!-------------------------------------------------
    call r8inir(6, 0.d0, epsi, 1)
    call r8inir(6, 0.d0, epsf, 1)
    call r8inir(6, 0.d0, epst, 1)
!
    do k = 1, 6
        epsf(k) = epsm(k)+deps(k)
        epsi(k) = epsm(k)
        epst(k) = (epsf(k)+epsi(k))/deux
    end do
!
    reinit = .false.
!
999 continue
    if (( &
        ( &
        (epsi(1) .ne. epsf(1)) .or. (epsi(2) .ne. epsf(2)) .or. (epsi(3) .ne. epsf(3)) .or. &
        (epsi(4) .ne. epsf(4)) .or. (epsi(5) .ne. epsf(5)) .or. (epsi(6) .ne. epsf(6)) &
        ) &
        .or. reinit &
        ) &
        .and. (compte .le. 100)) then
!
        reinit = .false.
        compte = compte+1
!
        if (compte .eq. 100) then
            iret = 0
            goto 9999
        end if
!
        call diago3(bm, vecbm, valbm)
        bdim = 3
        do i = 1, 3
            if (valbm(i)-tolb .le. 0.d0) then
                bdim = bdim-1
            end if
        end do
!
        trepsm = epst(1)+epst(2)+epst(3)
        if (trepsm .gt. 0.d0) then
            trepsm = 0.d0
        end if
!
        seuil = rk-rk1*trepsm*(atan2(-trepsm/rk2, un))
!
!----CAS OU LES 3 VALEURS PROPRES SONT NON NULLES---------------------
        if (bdim .eq. 3) then
!
            call lceob3(intmax, toler, epst, bm, dm, &
                        lambda, mu, alpha, ecrob, ecrod, &
                        seuil, bdim, b, d, mult, &
                        elas, dbloq, iret)
!
            call diago3(b, vecb, valb)
            reinit = .false.
            if (compte .lt. 100) then
                do i = 1, 3
                    if ((valb(i) .lt. 0) .or. (d .gt. 1.d0)) then
                        reinit = .true.
                    else
                        if (valb(i)-tolb .le. 0.d0) then
                            valb(i) = tolb-r8prem()
                        end if
                        if (un-d-tolb .le. 0.d0) then
                            d = un-tolb+r8prem()
                            dbloq = .true.
                        end if
                    end if
                end do
!
                if (reinit) then
                    do i = 1, 6
                        epst(i) = (epst(i)+epsi(i))/deux
                    end do
                    goto 999
                else
                    call r8inir(6, 0.d0, b, 1)
                    call r8inir(6, 0.d0, bm, 1)
                    do i = 1, 3
                        do j = i, 3
                            do k = 1, 3
                                b(t(i, j)) = b(t(i, j))+vecb(i, k)*valb(k)* &
                                             vecb(j, k)
                                bm(t(i, j)) = bm(t(i, j))+vecb(i, k)*valb( &
                                              k)*vecb(j, k)
                            end do
                        end do
                    end do
                    dm = d
                    do i = 1, 6
                        epsi(i) = epst(i)
                        epst(i) = epsf(i)
                    end do
                    goto 999
                end if
            else
                do i = 1, 3
                    if ((valb(i) .lt. 0) .and. (abs(valb(i))-tolb .le. 0.d0)) then
                        valb(i) = tolb-r8prem()
                    end if
                end do
                if (abs(un-d)-tolb .le. 0.d0) then
                    d = un-tolb+r8prem()
                    dbloq = .true.
                end if
                call r8inir(6, 0.d0, b, 1)
                call r8inir(6, 0.d0, bm, 1)
                do i = 1, 3
                    do j = i, 3
                        do k = 1, 3
                            b(t(i, j)) = b(t(i, j))+vecb(i, k)*valb(k)* &
                                         vecb(j, k)
                            bm(t(i, j)) = bm(t(i, j))+vecb(i, k)*valb(k)* &
                                          vecb(j, k)
                        end do
                    end do
                end do
                dm = d
                do i = 1, 6
                    epsi(i) = epst(i)
                    epst(i) = epsf(i)
                end do
!
            end if
!
!----CAS OU 1 VALEUR PROPRE EST NULLE---------------------------------
!
        else if (bdim .eq. 2) then
!
            call r8inir(9, 0.d0, interm, 1)
            call r8inir(6, 0.d0, epi, 1)
            do i = 1, 3
                do l = 1, 3
                    do k = 1, 3
                        interm(i, l) = interm(i, l)+vecbm(k, i)*epst(t(k, l) &
                                                                     )
                    end do
                    do j = i, 3
                        epi(t(i, j)) = epi(t(i, j))+interm(i, l)*vecbm(l, j)
                    end do
                end do
            end do
            tot1 = .false.
            tot2 = .false.
            tot3 = .false.
            call r8inir(6, 0.d0, bmr, 1)
            if (valbm(1)-tolb .le. 0.d0) then
!
                bmr(1) = valbm(2)
                bmr(2) = valbm(3)
                epsr(1) = epi(2)
                epsr(2) = epi(3)
                epsr(3) = epi(1)
                epsr(4) = epi(6)
                epsr(5) = epi(4)
                epsr(6) = epi(5)
                tot1 = .true.
            else if (valbm(2)-tolb .le. 0.d0) then
!
                bmr(1) = valbm(3)
                bmr(2) = valbm(1)
                epsr(1) = epi(3)
                epsr(2) = epi(1)
                epsr(3) = epi(2)
                epsr(4) = epi(5)
                epsr(5) = epi(6)
                epsr(6) = epi(4)
                tot2 = .true.
!
            else if (valbm(3)-tolb .le. 0.d0) then
!
                bmr(1) = valbm(1)
                bmr(2) = valbm(2)
                epsr(1) = epi(1)
                epsr(2) = epi(2)
                epsr(3) = epi(3)
                epsr(4) = epi(4)
                epsr(5) = epi(5)
                epsr(6) = epi(6)
                tot3 = .true.
!
            end if
!
            call lceob2(intmax, toler, epsr, bmr, dm, &
                        lambda, mu, alpha, ecrob, ecrod, &
                        seuil, bdim, br, d, mult, &
                        elas, dbloq, iret)
!
            call diago3(br, vecbr, valbr)
!
            if (compte .lt. 100) then
!
                do i = 1, 2
                    if (valbr(i) .lt. 0) then
                        reinit = .true.
                    end if
                    if (valbr(i)-tolb .le. 0.d0) then
                        valbr(i) = tolb-r8prem()
                    end if
                end do
                if (d .gt. 1.d0) then
                    reinit = .true.
                end if
                if (un-d-tolb .le. 0.d0) then
                    d = un-tolb+r8prem()
                    dbloq = .true.
                end if
            else
!
                reinit = .false.
                do i = 1, 2
                    if (valbr(i)-tolb .le. 0.d0) then
                        valbr(i) = tolb-r8prem()
                    end if
                end do
                if (d-(un-tolb) .ge. 0.d0) then
                    d = un-tolb+r8prem()
                    dbloq = .true.
                end if
!
            end if
!
            if (reinit) then
                do i = 1, 6
                    epst(i) = (epst(i)+epsi(i))/2
                end do
                goto 999
            else
!
                call r8inir(6, 0.d0, br, 1)
                do i = 1, 3
                    do j = i, 3
                        do k = 1, 3
                            br(t(i, j)) = br(t(i, j))+vecbr(i, k)*valbr(k)* &
                                          vecbr(j, k)
                        end do
                    end do
                end do
!
                if (tot1) then
                    binter(1) = tolb-r8prem()
                    binter(2) = br(1)
                    binter(3) = br(2)
                    binter(4) = 0.d0
                    binter(5) = 0.d0
                    binter(6) = br(4)
                else if (tot2) then
                    binter(1) = br(2)
                    binter(2) = tolb-r8prem()
                    binter(3) = br(1)
                    binter(4) = 0.d0
                    binter(5) = br(4)
                    binter(6) = 0.d0
                else if (tot3) then
                    binter(1) = br(1)
                    binter(2) = br(2)
                    binter(3) = tolb-r8prem()
                    binter(4) = br(4)
                    binter(5) = 0.d0
                    binter(6) = 0.d0
                end if
!
                call r8inir(9, 0.d0, interm, 1)
                call r8inir(6, 0.d0, b, 1)
                call r8inir(6, 0.d0, bm, 1)
                do i = 1, 3
                    do l = 1, 3
                        do k = 1, 3
                            interm(i, l) = interm(i, l)+vecbm(i, k)*binter( &
                                           t(k, l))
                        end do
                        do j = i, 3
                            b(t(i, j)) = b(t(i, j))+interm(i, l)*vecbm(j, l)
                            bm(t(i, j)) = bm(t(i, j))+interm(i, l)*vecbm(j, &
                                                                         l)
                        end do
                    end do
                end do
                dm = d
!
                do i = 1, 6
                    epsi(i) = epst(i)
                    epst(i) = epsf(i)
                end do
                goto 999
            end if
!
!---- CAS OU 2 VALEURS PROPRES SONT NULLES-----------------------------
!
        else if (bdim .eq. 1) then
!
            call r8inir(9, 0.d0, interm, 1)
            call r8inir(6, 0.d0, epi, 1)
            do i = 1, 3
                do l = 1, 3
                    do k = 1, 3
                        interm(i, l) = interm(i, l)+vecbm(k, i)*epst(t(k, l) &
                                                                     )
                    end do
                    do j = i, 3
                        epi(t(i, j)) = epi(t(i, j))+interm(i, l)*vecbm(l, j)
                    end do
                end do
            end do
!
            tot1 = .false.
            tot2 = .false.
            tot3 = .false.
            call r8inir(6, 0.d0, bmr, 1)
            if (valbm(1)-tolb .gt. 0.d0) then
                bmr(1) = valbm(1)
                epsr(1) = epi(1)
                epsr(2) = epi(2)
                epsr(3) = epi(3)
                epsr(4) = epi(4)
                epsr(5) = epi(5)
                epsr(6) = epi(6)
                tot1 = .true.
!
            else if (valbm(2)-tolb .gt. 0.d0) then
                bmr(1) = valbm(2)
                epsr(1) = epi(2)
                epsr(2) = epi(3)
                epsr(3) = epi(1)
                epsr(4) = epi(6)
                epsr(5) = epi(4)
                epsr(6) = epi(5)
                tot2 = .true.
!
            else if (valbm(3)-tolb .gt. 0.d0) then
                bmr(1) = valbm(3)
                epsr(1) = epi(3)
                epsr(2) = epi(1)
                epsr(3) = epi(2)
                epsr(4) = epi(5)
                epsr(5) = epi(6)
                epsr(6) = epi(4)
                tot3 = .true.
            end if
!
            call lceob1(intmax, toler, epsr, bmr, dm, &
                        lambda, mu, alpha, ecrob, ecrod, &
                        seuil, bdim, br, d, mult, &
                        elas, dbloq, iret)
!
            if (compte .lt. 100) then
!
                if (br(1) .lt. 0) then
                    reinit = .true.
                end if
                if (br(1)-tolb .le. 0.d0) then
                    br(1) = tolb-r8prem()
                end if
                if (d .gt. 1.d0) then
                    reinit = .true.
                end if
                if (un-d-tolb .le. 0.d0) then
                    d = un-tolb+r8prem()
                    dbloq = .true.
                end if
!
            else
!
                reinit = .false.
                if (br(1)-tolb .le. 0.d0) then
                    br(1) = tolb-r8prem()
                end if
                if (d-(un-tolb) .ge. 0.d0) then
                    d = un-tolb+r8prem()
                    dbloq = .true.
                end if
!
            end if
!
            if (reinit) then
                do i = 1, 6
                    epst(i) = (epst(i)+epsi(i))/2
                end do
                goto 999
!
            else
                valb(1) = tolb-r8prem()
                valb(2) = tolb-r8prem()
                valb(3) = tolb-r8prem()
                if (tot1) valb(1) = br(1)
                if (tot2) valb(2) = br(1)
                if (tot3) valb(3) = br(1)
                call r8inir(6, 0.d0, b, 1)
                call r8inir(6, 0.d0, bm, 1)
                do i = 1, 3
                    do j = i, 3
                        do k = 1, 3
                            b(t(i, j)) = b(t(i, j))+vecbm(i, k)*valb(k)* &
                                         vecbm(j, k)
                            bm(t(i, j)) = bm(t(i, j))+vecbm(i, k)*valbm(k)* &
                                          vecbm(k, j)
                        end do
                    end do
                end do
                dm = d
                do i = 1, 6
                    epsi(i) = epst(i)
                    epst(i) = epsf(i)
                end do
                goto 999
            end if
        end if
!
    end if
9999 continue
!
end subroutine
