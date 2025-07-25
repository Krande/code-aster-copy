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
subroutine preml0(n1, n2, diag, col, delg, &
                  prno, deeq, nec, p, q, &
                  lbd1, lbd2, rl, rl1, rl2, &
                  nrl, lt, lmat)
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: n1, diag(0:*), col(*)
    integer(kind=8) :: delg(*), prno(*), deeq(*), nec, lbd1(n1), lbd2(n1)
    integer(kind=8) :: rl(4, *), rl1(*), rl2(*)
    integer(kind=8) :: p(*), q(*)
!     VARIABLES LOCALES
    integer(kind=8) :: nrl, lt, n2, ino, num, nobl, i, j, lmat, i2, iddl, ier, ifm, niv
    integer(kind=8) :: idiai, idiai1, ii, li, iconne, nfois, vali(3)
    aster_logical :: nivdbg
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    nivdbg = .false.
    nfois = 0
    iconne = 0
    call infniv(ifm, niv)
!
!---------------------------------------------INITIALISATIONS
    diag(0) = 0
    do i = 1, n1
        p(i) = 0
        lbd1(i) = 0
        lbd2(i) = 0
        q(i) = 0
        rl1(i) = 0
        rl2(i) = 0
    end do
!---------------------------------------------CALCUL DE ADJNC1
    lmat = diag(n1)
    n2 = 0
    nrl = 0
    do iddl = 1, n1
        if (delg(iddl) .eq. 0) then
            n2 = n2+1
            p(iddl) = n2
            q(n2) = iddl
        else
            ino = deeq(2*iddl-1)
            if (ino .ne. 0) then
!     IDDL EST UN LAGRANGE DE BLOCAGE
                num = -deeq(2*iddl)
                if (num .eq. 0) then
                    vali(1) = iddl
                    vali(2) = ino
                    vali(3) = num
                    call utmess('F', 'ALGELINE5_31', ni=3, vali=vali)
                end if
                nobl = prno((nec+2)*(ino-1)+1)
!     RECHERCHE DE NOBL : NUMERO DU DDL BLOQUE
!     DO WHILE (DEEQ(2*NOBL).NE.NUM)
20              continue
                if (deeq(2*nobl) .ne. num) then
                    nobl = nobl+1
                    ASSERT(nobl .le. n1)
                    goto 20
!     FIN DO WHILE
                end if
                if (delg(iddl) .eq. -1) then
                    if (lbd1(nobl) .ne. 0) nfois = nfois+1
                    lbd1(nobl) = iddl
                else if (delg(iddl) .eq. -2) then
                    if (lbd2(nobl) .ne. 0) nfois = nfois+1
                    lbd2(nobl) = iddl
                else
                    vali(1) = delg(iddl)
                    call utmess('F', 'ALGELINE5_32', si=vali(1))
                end if
                if (nfois .gt. 0) then
                    vali(1) = nobl
                    call utmess('F', 'ALGELINE5_33', si=vali(1))
                end if
            else
!     IDDL EST UN LAGRANGE DE RELATION LINEAIRE
!     POUR CHQE REL. LIN. I,ON A
!     RL(2,I) = LAMBDA2 ET RL(1,I) = LAMBDA1
                if (delg(iddl) .eq. -2) then
                    nrl = nrl+1
                    rl(2, nrl) = iddl
!     RL(1,NRL) SERA DEFINI DANS PREMLC, COMME LE NO DE COLONNE
!     DU 1ER TERME DE LA LIGNE RL(2,NRL).
                end if
            end if
        end if
    end do
!     CALCUL DE LA TAILLE DE LA LISTE
    lt = 0
    do i = 1, nrl
        i2 = rl(2, i)
        lt = lt+(diag(i2)-diag(i2-1))
    end do
!     ON MAJORE LT POUR LES PETITS CAS-TESTS
    if (lt .le. 10) then
        lt = lt**2
    else
        lt = lt*10
    end if
!
!     VERIFICATION DES CONNEXIONS DES LAGRANGES
    if (nivdbg) then
        do i = 1, n1
            li = lbd1(i)
            if (li .ne. 0) then
                idiai1 = diag(li-1)+1
                idiai = diag(li)
                if (idiai1 .lt. idiai) then
!
                    write (ifm, *) 'LE DDL BLOQUE: ', i, ' A POUR LAMBDA1: ', lbd1(i)
                    write (ifm, *) 'LE DDL BLOQUE: ', i, ' A POUR LAMBDA2: ', lbd2(i)
                    write (ifm, *) 'LE LAMBDA1 ', lbd1(i),&
     &               ' A POUR VOISIN INATTENDUS '
                    do j = idiai1, idiai-1
                        write (ifm, *) 'LE DDL ', col(j)
                        iconne = iconne+1
                    end do
                end if
                do ii = li+1, n1
                    idiai1 = diag(ii-1)+1
                    idiai = diag(ii)
                    do j = idiai1, idiai
                        if (col(j) .eq. li) then
                            if (ii .ne. i .and. ii .ne. lbd2(i)) then
                                write (ifm, *) 'LE DDL BLOQUE: ', i, &
                                    ' A POUR LAMBDA1: ', lbd1(i)
                                write (ifm, *) 'LE DDL BLOQUE: ', i, &
                                    ' A POUR LAMBDA2: ', lbd2(i)
                                write (ifm, *) 'LE LAMBDA1 ', lbd1(i), &
                                    ' A POUR VOISIN INATTENDU', ii
                                iconne = iconne+1
                            end if
                        end if
                    end do
!
                end do
            end if
        end do
        if (iconne .gt. 0) then
            call utmess('A', 'ALGELINE5_53')
            write (ifm, *) 2*iconne, ' TERMES SUPPLEMENTAIRES DANS'&
     &'    LA MATRICE INITIALE'
        end if
    end if
!
    if (niv .eq. 2) then
        ier = 0
        do i = 1, n1
            if (lbd1(i) .ne. 0) then
!            WRITE (IFM,*) 'LE DDL BLOQUE: ',I,' A POUR LAMBDA1: ',
!     &        LBD1(I)
!            WRITE (IFM,*) 'LE DDL BLOQUE: ',I,' A POUR LAMBDA2: ',
!     &        LBD2(I)
!            IF (LBD2(I).EQ.0) IER = 1
            else if (lbd2(i) .ne. 0) then
                ier = 1
            end if
            if (ier .eq. 1) then
                vali(1) = i
                vali(2) = lbd1(i)
                vali(3) = lbd1(i)
                call utmess('F', 'ALGELINE5_34', ni=3, vali=vali)
            end if
!
        end do
    end if
end subroutine
