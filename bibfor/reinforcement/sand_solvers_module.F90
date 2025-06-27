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

module sand_solvers_module

    implicit none
    public :: solver_sandcas1, solver_sandcas2, solver_optimum

contains

!---------------------------------------------------------------------
!! SUBROUTINE 1
!---------------------------------------------------------------------
    subroutine solver_sandcas1(N_INF, N_SUP, N_TOT, ITER, MONO, p, &
                               nSX_INF, nSY_INF, nSX_SUP, nSY_SUP, nS_TOT, &
                               tINF, tSUP, ncMAX_INF, ncMIN_INF, ncMAX_SUP, ncMIN_SUP, &
                               RESIDU, AngleINF, AngleSUP)

#include "asterfort/wkvect.h"
#include "asterfort/juveca.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jedetr.h"

!Variables principales
        integer(kind=8) :: N_INF
        integer(kind=8) :: N_SUP
        integer(kind=8) :: N_TOT
        integer(kind=8) :: ITER
        logical :: MONO
        character(20) :: p(15)
        real(kind=8), pointer :: nSX_INF(:)
        real(kind=8), pointer :: nSY_INF(:)
        real(kind=8), pointer :: nSX_SUP(:)
        real(kind=8), pointer :: nSY_SUP(:)
        real(kind=8), pointer :: nS_TOT(:)
        real(kind=8), pointer :: tINF(:)
        real(kind=8), pointer :: tSUP(:)
        real(kind=8), pointer :: ncMAX_INF(:)
        real(kind=8), pointer :: ncMIN_INF(:)
        real(kind=8), pointer :: ncMAX_SUP(:)
        real(kind=8), pointer :: ncMIN_SUP(:)
        real(kind=8), pointer :: RESIDU(:)
        real(kind=8), pointer :: AngleINF(:)
        real(kind=8), pointer :: AngleSUP(:)

!Variables intermediaires de calcul
        integer(kind=8) :: N1, N2, indx1, indx2, indx3, COUNT_SOL, i, j, k, iBIS, N_Theta_ADD
        integer(kind=8) :: COUNT_Q, INDICE, Q
        logical :: COND_nSX_INF, COND_nSY_INF, COND_nSX_SUP, COND_nSY_SUP, COND_insert
        real(kind=8) :: Theta_INTER_nSX_INF, Theta_INTER_nSX_SUP
        real(kind=8) :: Theta_INTER_nSY_INF, Theta_INTER_nSY_SUP
        real(kind=8) :: angle1, angle2, angle3, suppl, suppl1, suppl2, suppl3, suppl4
        real(kind=8) :: Theta_INTER(4, 4), a12, a13
        real(kind=8), pointer :: Theta_ADD(:) => null()
        real(kind=8) :: debug1, debug2

!Consideration de l'axe d'iteration pour recherche des racines eventuelles
        !ITER=1 ==> INF'
        !ITER=2 ==> SUP'
        if (ITER .eq. 1) then
            N1 = N_SUP
            N2 = N_INF
        elseif (ITER .eq. 2) then
            N1 = N_INF
            N2 = N_SUP
        end if

!Consideration supplémentaire pour AngleINF et AngleSUP
        if (MONO .eqv. (.true.)) then
            !SAND_CAS = 2
            !ITER = 2
            !AngleINF(1:N_SUP) et 'N_INF'=N1=1
            !
            !SAND_CAS = 3
            !ITER = 1
            !AngleSUP(1:N_INF) et 'N_SUP'=N1=1
            N1 = 1
        end if

        do i = 1, N1
            j = 1
            do while (j .lt. N2)
                COND_nSX_SUP = .false.
                COND_nSX_INF = .false.
                COND_nSY_SUP = .false.
                COND_nSY_INF = .false.
                COND_insert = .false.

                if (ITER .eq. 1) then
                    indx1 = i+N1*(j-1)
                    indx2 = i+N1*j
                elseif (ITER .eq. 2) then
                    indx1 = j+N2*(i-1)
                    indx2 = (j+1)+N2*(i-1)
                end if

                if ((nSX_SUP(indx1) .ne. (-1.d0)) .and. (nSX_INF(indx1) .ne. (-1.d0)) &
                     & .and. (nSY_SUP(indx1) .ne. (-1.d0)) .and. (nSY_INF(indx1) .ne. (-1.d0)) &
                     & .and. (nSX_SUP(indx2) .ne. (-1.d0)) .and. (nSX_INF(indx2) .ne. (-1.d0)) &
                     & .and. (nSY_SUP(indx2) .ne. (-1.d0)) .and. (nSY_INF(indx2) .ne. (-1.d0))) then

                    ! Exclude cases where we have zeros

                    debug1 = nSX_SUP(indx1)
                    debug2 = nSX_SUP(indx2)
                    if (((debug1*debug2) .lt. 0) &
                       & .and. (abs(debug1) .le. epsilon(debug1)) &
                       & .and. (abs(debug2) .le. epsilon(debug2))) then
                        COND_nSX_SUP = .true.
                    end if

                    debug1 = nSX_INF(indx1)
                    debug2 = nSX_INF(indx2)
                    if (((debug1*debug2) .lt. 0) &
                       & .and. (abs(debug1) .le. epsilon(debug1)) &
                       & .and. (abs(debug2) .le. epsilon(debug2))) then
                        COND_nSX_INF = .true.
                    end if

                    debug1 = nSY_SUP(indx1)
                    debug2 = nSY_SUP(indx2)
                    if (((debug1*debug2) .lt. 0) &
                       & .and. (abs(debug1) .le. epsilon(debug1)) &
                       & .and. (abs(debug2) .le. epsilon(debug2))) then
                        COND_nSY_SUP = .true.
                    end if

                    debug1 = nSY_INF(indx1)
                    debug2 = nSY_INF(indx2)
                    if (((debug1*debug2) .lt. 0) &
                       & .and. (abs(debug1) .le. epsilon(debug1)) &
                       & .and. (abs(debug2) .le. epsilon(debug2))) then
                        COND_nSY_INF = .true.
                    end if

                    if ((COND_nSX_SUP .eqv. (.true.)) .or. (COND_nSX_INF .eqv. (.true.)) &
                      & .or. (COND_nSY_SUP .eqv. (.true.)) .or. (COND_nSY_INF .eqv. (.true.))) then
                        COND_insert = .true.
                    end if

                    if (COND_insert .eqv. (.true.)) then

                        COUNT_SOL = 0
                        if (COND_nSX_SUP .eqv. (.true.)) then
                            COUNT_SOL = COUNT_SOL+1
                        end if
                        if (COND_nSX_INF .eqv. (.true.)) then
                            COUNT_SOL = COUNT_SOL+1
                        end if
                        if (COND_nSY_SUP .eqv. (.true.)) then
                            COUNT_SOL = COUNT_SOL+1
                        end if
                        if (COND_nSY_INF .eqv. (.true.)) then
                            COUNT_SOL = COUNT_SOL+1
                        end if

                        N2 = N2+COUNT_SOL
                        N_TOT = N1*N2

                        do k = 1, 12
                            call juveca(p(k), N_TOT)
                        end do

                        call jeveuo(p(1), 'E', vr=nSX_SUP)
                        call jeveuo(p(2), 'E', vr=nSY_SUP)
                        call jeveuo(p(3), 'E', vr=nSX_INF)
                        call jeveuo(p(4), 'E', vr=nSY_INF)
                        call jeveuo(p(5), 'E', vr=nS_TOT)
                        call jeveuo(p(6), 'E', vr=tSUP)
                        call jeveuo(p(7), 'E', vr=tINF)
                        call jeveuo(p(8), 'E', vr=ncMAX_SUP)
                        call jeveuo(p(9), 'E', vr=ncMIN_SUP)
                        call jeveuo(p(10), 'E', vr=ncMAX_INF)
                        call jeveuo(p(11), 'E', vr=ncMIN_INF)
                        call jeveuo(p(12), 'E', vr=RESIDU)

                        if (ITER .eq. 1) then
                            call juveca(p(14), N2)
                            call jeveuo(p(14), 'E', vr=AngleINF)
                            do k = 0, (N2-j-(COUNT_SOL+1))
                                AngleINF(N2-k) = AngleINF(N2-k-COUNT_SOL)
                            end do

                        elseif (ITER .eq. 2) then
                            call juveca(p(13), N2)
                            call jeveuo(p(13), 'E', vr=AngleSUP)
                            do k = 0, (N2-j-(COUNT_SOL+1))
                                AngleSUP(N2-k) = AngleSUP(N2-k-COUNT_SOL)
                            end do

                        end if

                        do iBIS = 1, N1
                        do k = 0, (N2-j-(COUNT_SOL+1))
                            if (ITER .eq. 1) then
                                indx1 = iBIS+N1*(N2-k-1)
                                indx2 = iBIS+N1*(N2-k-COUNT_SOL-1)
                            elseif (ITER .eq. 2) then
                                indx1 = (N2-k)+N2*(iBIS-1)
                                indx2 = (N2-k-COUNT_SOL)+N2*(iBIS-1)
                            end if
                            nSX_SUP(indx1) = nSX_SUP(indx2)
                            nSX_INF(indx1) = nSX_INF(indx2)
                            nSY_SUP(indx1) = nSY_SUP(indx2)
                            nSY_INF(indx1) = nSY_INF(indx2)
                            nS_TOT(indx1) = nS_TOT(indx2)
                            tSUP(indx1) = tSUP(indx2)
                            tINF(indx1) = tINF(indx2)
                            ncMAX_SUP(indx1) = ncMAX_SUP(indx2)
                            ncMIN_SUP(indx1) = ncMIN_SUP(indx2)
                            ncMAX_INF(indx1) = ncMAX_INF(indx2)
                            ncMIN_INF(indx1) = ncMIN_INF(indx2)
                            RESIDU(indx1) = RESIDU(indx2)
                        end do
                        end do

                        !NOW WE MUST CONSTRUCT THE NEW INSTANCES

                        Theta_INTER_nSX_SUP = -1.d0
                        Theta_INTER_nSX_INF = -1.d0
                        Theta_INTER_nSY_SUP = -1.d0
                        Theta_INTER_nSY_INF = -1.d0

                        if (ITER .eq. 1) then
                            indx1 = i+N1*(j-1)
                            indx2 = i+N1*((j+COUNT_SOL+1)-1)
                        elseif (ITER .eq. 2) then
                            indx1 = j+N2*(i-1)
                            indx2 = (j+COUNT_SOL+1)+N2*(i-1)
                        end if

                        if (ITER .eq. 1) then
                            angle1 = AngleINF(j)
                            angle2 = AngleINF(j+COUNT_SOL+1)
                        elseif (ITER .eq. 2) then
                            angle1 = AngleSUP(j)
                            angle2 = AngleSUP(j+COUNT_SOL+1)
                        end if

                        if (COND_nSX_SUP .eqv. (.true.)) then
                            suppl = (angle2-angle1)/(1+Abs(nSX_SUP(indx2))/Abs(nSX_SUP(indx1)))
                            Theta_INTER_nSX_SUP = angle1+suppl
                        end if
                        if (COND_nSX_INF .eqv. (.true.)) then
                            suppl = (angle2-angle1)/(1+Abs(nSX_INF(indx2))/Abs(nSX_INF(indx1)))
                            Theta_INTER_nSX_INF = angle1+suppl
                        end if
                        if (COND_nSY_SUP .eqv. (.true.)) then
                            suppl = (angle2-angle1)/(1+Abs(nSY_SUP(indx2))/Abs(nSY_SUP(indx1)))
                            Theta_INTER_nSY_SUP = angle1+suppl
                        end if
                        if (COND_nSY_INF .eqv. (.true.)) then
                            suppl = (angle2-angle1)/(1+Abs(nSY_INF(indx2))/Abs(nSY_INF(indx1)))
                            Theta_INTER_nSY_INF = angle1+suppl
                        end if

                        Theta_INTER(1, 1) = Theta_INTER_nSX_SUP
                        Theta_INTER(2, 1) = Theta_INTER_nSX_INF
                        Theta_INTER(3, 1) = Theta_INTER_nSY_SUP
                        Theta_INTER(4, 1) = Theta_INTER_nSY_INF
                        Theta_INTER(1, 2) = 1.d0
                        Theta_INTER(2, 2) = 2.d0
                        Theta_INTER(3, 2) = 3.d0
                        Theta_INTER(4, 2) = 4.d0
                        Theta_INTER(1, 3) = -1.d0
                        Theta_INTER(2, 3) = -1.d0
                        Theta_INTER(3, 3) = -1.d0
                        Theta_INTER(4, 3) = -1.d0

                        N_Theta_ADD = 2*COUNT_SOL
                        call wkvect(p(15), ' V V R ', N_Theta_ADD, vr=Theta_ADD)
                        do k = 1, N_Theta_ADD
                            Theta_ADD(k) = -1.d0
                        end do

                        do k = 1, COUNT_SOL
                            COUNT_Q = 0
                            INDICE = 0
                            do Q = 1, 4
                                if (((COUNT_Q .eq. 0) .or. (Theta_INTER(Q, 1) .le. Theta_ADD(k))) &
                                   & .and. (Theta_INTER(Q, 1) .ne. (-1.d0)) &
                                   & .and. (Theta_INTER(Q, 3) .eq. (-1.d0))) then
                                    COUNT_Q = COUNT_Q+1
                                    Theta_ADD(k) = Theta_INTER(Q, 1)
                                    Theta_ADD(k+COUNT_SOL) = Theta_INTER(Q, 2)
                                    INDICE = Q
                                end if
                            end do
                            if (INDICE .gt. 0) then
                                Theta_INTER(INDICE, 3) = 1.d0
                            end if
                        end do

                        do k = 1, COUNT_SOL

                            if (ITER .eq. 1) then
                                AngleINF(j+k) = Theta_ADD(k)
                                angle3 = AngleINF(j+k)
                                indx3 = i+N1*((j+k)-1)
                            elseif (ITER .eq. 2) then
                                AngleSUP(j+k) = Theta_ADD(k)
                                angle3 = AngleSUP(j+k)
                                indx3 = (j+k)+N2*(i-1)
                            end if

                            suppl1 = nSX_SUP(indx2)-nSX_SUP(indx1)
                            suppl1 = suppl1/(angle2-angle1)
                            suppl1 = suppl1*(angle3-angle1)

                            suppl2 = nSX_INF(indx2)-nSX_INF(indx1)
                            suppl2 = suppl2/(angle2-angle1)
                            suppl2 = suppl2*(angle3-angle1)

                            suppl3 = nSY_SUP(indx2)-nSY_SUP(indx1)
                            suppl3 = suppl3/(angle2-angle1)
                            suppl3 = suppl3*(angle3-angle1)

                            suppl4 = nSY_INF(indx2)-nSY_INF(indx1)
                            suppl4 = suppl4/(angle2-angle1)
                            suppl4 = suppl4*(angle3-angle1)

                            if (Theta_ADD(k+COUNT_SOL) .eq. (1.d0)) then
                                nSX_SUP(indx3) = 0
                                nSX_INF(indx3) = nSX_INF(indx1)+suppl2
                                nSY_SUP(indx3) = nSY_SUP(indx1)+suppl3
                                nSY_INF(indx3) = nSY_INF(indx1)+suppl4
                            elseif (Theta_ADD(k+COUNT_SOL) .eq. (2.d0)) then
                                nSX_SUP(indx3) = nSX_SUP(indx1)+suppl1
                                nSX_INF(indx3) = 0
                                nSY_SUP(indx3) = nSY_SUP(indx1)+suppl3
                                nSY_INF(indx3) = nSY_INF(indx1)+suppl4
                            elseIf (Theta_ADD(k+COUNT_SOL) .eq. (3.d0)) then
                                nSX_SUP(indx3) = nSX_SUP(indx1)+suppl1
                                nSX_INF(indx3) = nSX_INF(indx1)+suppl2
                                nSY_SUP(indx3) = 0
                                nSY_INF(indx3) = nSY_INF(indx1)+suppl4
                            elseif (Theta_ADD(k+COUNT_SOL) .eq. (4.d0)) then
                                nSX_SUP(indx3) = nSX_SUP(indx1)+suppl1
                                nSX_INF(indx3) = nSX_INF(indx1)+suppl2
                                nSY_SUP(indx3) = nSY_SUP(indx1)+suppl3
                                nSY_INF(indx3) = 0
                            end if
                            a12 = (angle2-angle1)
                            a13 = (angle3-angle1)
                            nS_TOT(indx3) = Abs(nSX_SUP(indx3))+Abs(nSX_INF(indx3))
                            nS_TOT(indx3) = nS_TOT(indx3)+Abs(nSY_SUP(indx3))+Abs(nSY_INF(indx3))
                            suppl = (tSUP(indx2)-tSUP(indx1))
                            tSUP(indx3) = tSUP(indx1)+(suppl/a12)*a13
                            suppl = (tINF(indx2)-tINF(indx1))
                            tINF(indx3) = tINF(indx1)+(suppl/a12)*a13
                            suppl = (ncMAX_SUP(indx2)-ncMAX_SUP(indx1))
                            ncMAX_SUP(indx3) = ncMAX_SUP(indx1)+(suppl/a12)*a13
                            suppl = (ncMIN_SUP(indx2)-ncMIN_SUP(indx1))
                            ncMIN_SUP(indx3) = ncMIN_SUP(indx1)+(suppl/a12)*a13
                            suppl = (ncMAX_INF(indx2)-ncMAX_INF(indx1))
                            ncMAX_INF(indx3) = ncMAX_INF(indx1)+(suppl/a12)*a13
                            suppl = (ncMIN_INF(indx2)-ncMIN_INF(indx1))
                            ncMIN_INF(indx3) = ncMIN_INF(indx1)+(suppl/a12)*a13
                            suppl = (RESIDU(indx2)-RESIDU(indx1))
                            RESIDU(indx3) = RESIDU(indx1)+(suppl/a12)*a13
                        end do

                        do iBIS = 1, N1
                            if (iBIS .ne. i) then
                                do k = 1, COUNT_SOL

                                    if (ITER .eq. 1) then
                                        indx1 = iBIS+N1*(j-1)
                                        indx2 = iBIS+N1*((j+COUNT_SOL+1)-1)
                                        indx3 = iBIS+N1*((j+k)-1)
                                        angle3 = AngleINF(j+k)
                                    elseif (ITER .eq. 2) then
                                        indx1 = j+N2*(iBIS-1)
                                        indx2 = (j+COUNT_SOL+1)+N2*(iBIS-1)
                                        indx3 = (j+k)+N2*(iBIS-1)
                                        angle3 = AngleSUP(j+k)
                                    end if

                                    a12 = (angle2-angle1)
                                    a13 = (angle3-angle1)
                                    suppl1 = ((nSX_SUP(indx2)-nSX_SUP(indx1))/a12)*a13
                                    suppl2 = ((nSX_INF(indx2)-nSX_INF(indx1))/a12)*a13
                                    suppl3 = ((nSY_SUP(indx2)-nSY_SUP(indx1))/a12)*a13
                                    suppl4 = ((nSY_INF(indx2)-nSY_INF(indx1))/a12)*a13

                                    if ((nSX_SUP(indx1) .ne. (-1.d0)) &
                                         &.and. (nSX_SUP(indx2) .ne. (-1.d0))) then
                                        nSX_SUP(indx3) = nSX_SUP(indx1)+suppl1
                                    else
                                        nSX_SUP(indx3) = -1.d0
                                    end if
                                    if ((nSX_INF(indx1) .ne. (-1.d0)) &
                                         &.and. (nSX_INF(indx2) .ne. (-1.d0))) then
                                        nSX_INF(indx3) = nSX_INF(indx1)+suppl2
                                    else
                                        nSX_INF(indx3) = -1.d0
                                    end if
                                    if ((nSY_SUP(indx1) .ne. (-1.d0)) &
                                         &.and. (nSY_SUP(indx2) .ne. (-1.d0))) then
                                        nSY_SUP(indx3) = nSY_SUP(indx1)+suppl3
                                    else
                                        nSY_SUP(indx3) = -1.d0
                                    end if
                                    if ((nSY_INF(indx1) .ne. (-1.d0)) &
                                         &.and. (nSY_INF(indx2) .ne. (-1.d0))) then
                                        nSY_INF(indx3) = nSY_INF(indx1)+suppl4
                                    else
                                        nSY_INF(indx3) = -1.d0
                                    end if

                                    if ((nSX_SUP(indx3) .ne. (-1.d0)) &
                                         &.and. (nSX_INF(indx3) .ne. (-1.d0)) &
                                         &.and. (nSY_SUP(indx3) .ne. (-1.d0)) &
                                         &.and. (nSY_INF(indx3) .ne. (-1.d0))) then
                                        nS_TOT(indx3) = Abs(nSX_SUP(indx3))
                                        nS_TOT(indx3) = nS_TOT(indx3)+Abs(nSX_INF(indx3))
                                        nS_TOT(indx3) = nS_TOT(indx3)+Abs(nSY_SUP(indx3))
                                        nS_TOT(indx3) = nS_TOT(indx3)+Abs(nSY_INF(indx3))
                                        suppl = (tSUP(indx2)-tSUP(indx1))
                                        tSUP(indx3) = tSUP(indx1)+(suppl/a12)*a13
                                        suppl = (tINF(indx2)-tINF(indx1))
                                        tINF(indx3) = tINF(indx1)+(suppl/a12)*a13
                                        suppl = (ncMAX_SUP(indx2)-ncMAX_SUP(indx1))
                                        ncMAX_SUP(indx3) = ncMAX_SUP(indx1)+(suppl/a12)*a13
                                        suppl = (ncMIN_SUP(indx2)-ncMIN_SUP(indx1))
                                        ncMIN_SUP(indx3) = ncMIN_SUP(indx1)+(suppl/a12)*a13
                                        suppl = (ncMAX_INF(indx2)-ncMAX_INF(indx1))
                                        ncMAX_INF(indx3) = ncMAX_INF(indx1)+(suppl/a12)*a13
                                        suppl = (ncMIN_INF(indx2)-ncMIN_INF(indx1))
                                        ncMIN_INF(indx3) = ncMIN_INF(indx1)+(suppl/a12)*a13
                                    else
                                        nS_TOT(indx3) = -1.d0
                                        tSUP(indx3) = -1.d0
                                        tINF(indx3) = -1.d0
                                        ncMAX_SUP(indx3) = -1.d0
                                        ncMIN_SUP(indx3) = -1.d0
                                        ncMAX_INF(indx3) = -1.d0
                                        ncMIN_INF(indx3) = -1.d0
                                    end if
                                end do
                            end if
                        end do

                        j = j+COUNT_SOL+1
                        call jedetr(p(15))

                    else

                        j = j+1

                    end if
                    !Distinction pour COND_insert

                else

                    j = j+1

                end if
                !Distinction si "-1" existe ou pas

            end do
            !Pour iteration sur N2

        end do
!Pour iteration sur N1

    end subroutine solver_sandcas1

!---------------------------------------------------------------------
!! SUBROUTINE 2
!---------------------------------------------------------------------
    subroutine solver_optimum(N_INF, N_SUP, MONO, &
                              nSX_INF, nSY_INF, nSX_SUP, nSY_SUP, nS_TOT, &
                              tINF, tSUP, ncMAX_INF, ncMIN_INF, ncMAX_SUP, ncMIN_SUP, &
                              AngleINF, AngleSUP, &
                              fyd, ferrcomp, ferrsyme, slsyme, &
                              vect, ierr)

#include "asterfort/wkvect.h"
#include "asterfort/juveca.h"
#include "asterfort/jeveuo.h"

!Variables principales
        integer(kind=8) :: N_INF
        integer(kind=8) :: N_SUP
        logical :: MONO
        real(kind=8), pointer :: nSX_INF(:)
        real(kind=8), pointer :: nSY_INF(:)
        real(kind=8), pointer :: nSX_SUP(:)
        real(kind=8), pointer :: nSY_SUP(:)
        real(kind=8), pointer :: nS_TOT(:)
        real(kind=8), pointer :: tINF(:)
        real(kind=8), pointer :: tSUP(:)
        real(kind=8), pointer :: ncMAX_INF(:)
        real(kind=8), pointer :: ncMIN_INF(:)
        real(kind=8), pointer :: ncMAX_SUP(:)
        real(kind=8), pointer :: ncMIN_SUP(:)
        real(kind=8), pointer :: AngleINF(:)
        real(kind=8), pointer :: AngleSUP(:)
        real(kind=8) :: fyd
        integer(kind=8) :: ferrcomp
        integer(kind=8) :: ferrsyme
        real(kind=8) :: slsyme
        real(kind=8) :: vect(20)
        integer(kind=8) :: ierr

!VECT =
!1:dnsxi,2:dnsxs,3:dnsyi,4:dnsys,&
!5:etsxi,6:etsxs,7:etsyi,8:etsys,&
!9:snsxi,10:snsxs,11:snsyi,12:snsys,&
!13:ncmaxi,14:ncmini,15:ncmaxs,16:ncmins,&
!17:t_inf,18:t_sup,19:theta_inf,20:theta_sup

!Variables intermediaires de calcul
        real(kind=8) :: dnsxi, dnsxs, dnsyi, dnsys
        real(kind=8) :: etsxi, etsxs, etsyi, etsys
        real(kind=8) :: snsxi, snsxs, snsyi, snsys
        real(kind=8) :: ncmaxi, ncmini, ncmaxs, ncmins
        real(kind=8) :: t_inf, t_sup, theta_inf, theta_sup
        integer(kind=8) :: COUNT_SOL, indx, i, j, INDICEi, INDICEj, indx_INDICE
        logical :: cond_found
        real(kind=8) :: CalcX, CalcY, nS_TOT_opt, a12, a13

        COUNT_SOL = 0
        do j = 1, N_INF
        do i = 1, N_SUP

            indx = i+N_SUP*(j-1)
            if (nS_TOT(indx) .ne. (-1.d0)) then
                if (ferrcomp .eq. 0) then
                    if ((nSX_SUP(indx) .lt. 0) .or. (nSY_SUP(indx) .lt. 0) &
                         & .or. (nSX_INF(indx) .lt. 0) .or. (nSY_INF(indx) .lt. 0)) then
                        cond_found = .false.
                        ierr = 2
                    else
                        cond_found = .true.
                    end if
                else
                    cond_found = .true.
                end if
                if ((cond_found) .eqv. (.true.)) then
                    CalcX = abs(abs(nSX_SUP(indx))-abs(nSX_INF(indx)))
                    CalcY = abs(abs(nSY_SUP(indx))-abs(nSY_INF(indx)))
                    if (ferrsyme .eq. 1) then
                        if ((CalcX .le. slsyme) .and. (CalcY .le. slsyme)) then
                            cond_found = .true.
                        else
                            cond_found = .false.
                            ierr = 3
                        end if
                    else
                        cond_found = .true.
                    end if
                end if
            else
                cond_found = .false.
            end if

            if (cond_found .eqv. (.true.)) then
                COUNT_SOL = COUNT_SOL+1
                if (COUNT_SOL .eq. 1) then
                    INDICEi = i
                    INDICEj = j
                    indx_INDICE = INDICEi+N_SUP*(INDICEj-1)
                    nS_TOT_opt = nS_TOT(indx)
                else
                    if (nS_TOT(indx) .lt. nS_TOT_opt) then
                        INDICEi = i
                        INDICEj = j
                        indx_INDICE = INDICEi+N_SUP*(INDICEj-1)
                        nS_TOT_opt = nS_TOT(indx)
                    elseif (nS_TOT(indx) .eq. nS_TOT_opt) then
                        a12 = (tSUP(indx)+tINF(indx))
                        a13 = (tSUP(indx_INDICE)+tINF(indx_INDICE))
                        if (a12 .lt. a13) then
                            INDICEi = i
                            INDICEj = j
                            indx_INDICE = INDICEi+N_SUP*(INDICEj-1)
                        end if
                    end if
                end if
            end if
        end do
        end do

        if (COUNT_SOL .gt. 0) then
            dnsxs = Abs(nSX_SUP(indx_INDICE))/fyd
            dnsys = Abs(nSY_SUP(indx_INDICE))/fyd
            dnsxi = Abs(nSX_INF(indx_INDICE))/fyd
            dnsyi = Abs(nSY_INF(indx_INDICE))/fyd

            if (nSX_SUP(indx_INDICE) .ge. 0) then
                etsxs = 0.d0
            else
                etsxs = 1.d0
            end if

            if (nSY_SUP(indx_INDICE) .ge. 0) then
                etsys = 0.d0
            else
                etsys = 1.d0
            end if

            if (nSX_INF(indx_INDICE) .ge. 0) then
                etsxi = 0.d0
            else
                etsxi = 1.d0
            end if

            if (nSY_INF(indx_INDICE) .ge. 0) then
                etsyi = 0.d0
            else
                etsyi = 1.d0
            end if

            t_sup = tSUP(indx_INDICE)
            t_inf = tINF(indx_INDICE)

            theta_sup = AngleSUP(INDICEi)
            if (MONO .eqv. (.true.)) then
                theta_inf = AngleINF(INDICEi)
            else
                theta_inf = AngleINF(INDICEj)
            end if

            if (abs(t_sup) .ge. epsilon(t_sup)) then
                ncmaxs = ncMAX_SUP(indx_INDICE)/t_sup
                ncmins = ncMIN_SUP(indx_INDICE)/t_sup
            else
                ncmaxs = 0.d0
                ncmins = 0.d0
            end if

            if (abs(t_inf) .ge. epsilon(t_inf)) then
                ncmaxi = ncMAX_INF(indx_INDICE)/t_inf
                ncmini = ncMIN_INF(indx_INDICE)/t_inf
            else
                ncmaxi = 0.d0
                ncmini = 0.d0
            end if

            if (abs(dnsxs) .ge. epsilon(dnsxs)) then
                snsxs = fyd
            else
                snsxs = 0.d0
            end if
            if (abs(dnsys) .ge. epsilon(dnsys)) then
                snsys = fyd
            else
                snsys = 0.d0
            end if

            if (abs(dnsxi) .ge. epsilon(dnsxi)) then
                snsxi = fyd
            else
                snsxi = 0.d0
            end if
            if (abs(dnsyi) .ge. epsilon(dnsyi)) then
                snsyi = fyd
            else
                snsyi = 0.d0
            end if

        else

            dnsxs = -1.d0
            dnsys = -1.d0
            dnsxi = -1.d0
            dnsyi = -1.d0

            etsxs = 2
            etsys = 2
            etsxi = 2
            etsyi = 2

            t_sup = -1.d0
            t_inf = -1.d0
            theta_sup = -1.d0
            theta_inf = -1.d0

            ncmaxs = -1.d0
            ncmins = -1.d0
            ncmaxi = -1.d0
            ncmini = -1.d0

            snsxs = -1.d0
            snsys = -1.d0
            snsxi = -1.d0
            snsyi = -1.d0

            if (ierr .eq. 0) then
                ierr = 1
            end if

        end if

        vect(1) = dnsxi
        vect(2) = dnsxs
        vect(3) = dnsyi
        vect(4) = dnsys
        vect(5) = etsxi
        vect(6) = etsxs
        vect(7) = etsyi
        vect(8) = etsys
        vect(9) = snsxi
        vect(10) = snsxs
        vect(11) = snsyi
        vect(12) = snsys
        vect(13) = ncmaxi
        vect(14) = ncmini
        vect(15) = ncmaxs
        vect(16) = ncmins
        vect(17) = t_inf
        vect(18) = t_sup
        vect(19) = theta_inf
        vect(20) = theta_sup

    end subroutine solver_optimum

!---------------------------------------------------------------------
!! SUBROUTINE 3
!---------------------------------------------------------------------
    subroutine solver_sandcas2(ITER, yINF, ySUP, ht, effrts, fcd, fcd1, cond109, &
                               AngleSUP, tINF, &
                               tSUP, AngleINF, Residu, &
                               nSX_INF, nSY_INF, nSX_SUP, nSY_SUP, nS_TOT, &
                               ncMAX_INF, ncMIN_INF, ncMAX_SUP, ncMIN_SUP)

#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
#include "asterc/r8rddg.h"
#include "asterfort/mgauss.h"

!Variables principales
        integer(kind=8) :: ITER
        real(kind=8) :: yINF
        real(kind=8) :: ySUP
        real(kind=8) :: ht
        real(kind=8) :: effrts(6)
        real(kind=8) :: fcd
        real(kind=8) :: fcd1
        integer(kind=8) :: cond109
        real(kind=8) :: AngleSUP
        real(kind=8) :: tINF
        real(kind=8) :: tSUP
        real(kind=8) :: AngleINF
        real(kind=8) :: Residu
        real(kind=8) :: nSX_INF
        real(kind=8) :: nSY_INF
        real(kind=8) :: nSX_SUP
        real(kind=8) :: NSY_SUP
        real(kind=8) :: nS_TOT
        real(kind=8) :: ncMAX_INF
        real(kind=8) :: ncMIN_INF
        real(kind=8) :: ncMAX_SUP
        real(kind=8) :: ncMIN_SUP

!Variables intermediaires de calcul
        real(kind=8) :: Nxx, Nyy, Nxy, Mxx, Myy, Mxy
        real(kind=8) :: fc, a00, b00, c00, DELTA, xA, xB
        real(kind=8) :: Calc, Calc1, Calc2
        logical :: COND_xA, COND_xB
        real(kind=8) :: ncX_INF, ncY_INF, ncXY_INF, alpha, fcd2
        real(kind=8) ::nC_SUP, mC_SUP, theta_sup, pi
        real(kind=8) :: Ds(4, 4), SOL(4), det, denum
        integer(kind=8) :: iret

        Nxx = effrts(1)
        Nyy = effrts(2)
        Nxy = effrts(3)
        Mxx = effrts(4)
        Myy = effrts(5)
        Mxy = effrts(6)
        pi = r8pi()

        fc = fcd1

        Calc1 = abs(AngleSUP)
        Calc2 = abs(abs(AngleSUP)-90)
        theta_sup = AngleSUP*r8dgrd()
        if ((Calc1 .ge. epsilon(Calc1)) .and. (Calc2 .ge. epsilon(Calc2))) then
            a00 = -fc*sin(theta_sup)*cos(theta_sup)
            b00 = fc*sin(theta_sup)*cos(theta_sup)*(2*ht-tINF)
            c00 = (ht-tINF)*Nxy+2*Mxy

            if (abs(a00) .ge. epsilon(a00)) then
                DELTA = b00*b00-4*a00*c00
                if (DELTA .ge. 0) then
                    xA = (-b00+sqrt(DELTA))/(2*a00)
                    xB = (-b00-sqrt(DELTA))/(2*a00)
                    if ((xA .ge. 0) .and. (xA .le. (0.5*ht))) then
                        COND_xA = .true.
                    else
                        COND_xA = .false.
                    end if
                    if ((xB .ge. 0) .and. (xB .le. (0.5*ht))) then
                        COND_xB = .true.
                    else
                        COND_xB = .false.
                    end if
                else
                    COND_xA = .false.
                    COND_xB = .false.
                end if

                if ((COND_xA .eqv. (.true.)) .and. (COND_xB .eqv. (.true.))) then
                    tSUP = min(xA, xB)
                elseif (COND_xA .eqv. (.true.)) then
                    tSUP = xA
                elseif (COND_xB .eqv. (.true.)) then
                    tSUP = xB
                else
                    tSUP = -1.d0
                end if

            else

                if (abs(b00) .ge. epsilon(b00)) then
                    tSUP = -c00/b00
                    if ((tSUP .lt. 0) .or. (tSUP .gt. (0.5*ht))) then
                        tSUP = -1.d0
                    end if
                else
                    tSUP = -1.d0
                end if

            end if

        else
            Calc = Mxy+Nxy*0.5*(ht-tINF)
            if (abs(Calc) .lt. epsilon(Calc)) then
                tSUP = 0.d0
            else
                tSUP = -1.d0
            end if

        end if

        if (tSUP .ne. (-1.d0)) then

            ncXY_INF = Nxy+fc*sin(theta_sup)*cos(theta_sup)*tSUP
            nC_SUP = -tSUP*fc
            mC_SUP = 0.5*(ht-tSUP)*nC_SUP

            ncMAX_SUP = -nC_SUP
            ncMIN_SUP = 0

            Ds(1, 1) = 1
            Ds(1, 2) = 0
            Ds(1, 3) = 1
            Ds(1, 4) = 0
            Ds(2, 1) = 0
            Ds(2, 2) = 1
            Ds(2, 3) = 0
            Ds(2, 4) = 1
            !ITER=1 ==> REINF in SUP (SAND_CAS=2)
            !ITER=2 ==> REINF in INF (SAND_CAS=3)
            if (ITER .eq. 1) then
                Ds(3, 1) = ySUP
                Ds(3, 3) = -0.5*(ht-tINF)
                Ds(4, 2) = ySUP
                Ds(4, 4) = -0.5*(ht-tINF)
            elseif (ITER .eq. 2) then
                Ds(3, 1) = -ySUP
                Ds(3, 3) = 0.5*(ht-tINF)
                Ds(4, 2) = -ySUP
                Ds(4, 4) = 0.5*(ht-tINF)
            end if
            Ds(3, 2) = 0
            Ds(3, 4) = 0
            Ds(4, 1) = 0
            Ds(4, 3) = 0

            SOL(1) = Nxx-nC_SUP*((sin(theta_sup))**2)
            SOL(2) = Nyy-nC_SUP*((cos(theta_sup))**2)
            if (ITER .eq. 1) then
                SOL(3) = Mxx-mC_SUP*((sin(theta_sup))**2)
                SOL(4) = Myy-mC_SUP*((cos(theta_sup))**2)
            elseif (ITER .eq. 2) then
                SOL(3) = Mxx+mC_SUP*((sin(theta_sup))**2)
                SOL(4) = Myy+mC_SUP*((cos(theta_sup))**2)
            end if

            call mgauss('NFSP', Ds, SOL, 4, 4, 1, det, iret)

            nSX_SUP = SOL(1)
            nSY_SUP = SOL(2)
            nSX_INF = 0
            nSY_INF = 0
            ncX_INF = SOL(3)
            ncY_INF = SOL(4)

            ncMAX_INF = -(0.5*(ncX_INF+ncY_INF)-sqrt(((0.5*ncX_INF-0.5*ncY_INF)**2)+(ncXY_INF**2)))
            ncMIN_INF = -(0.5*(ncX_INF+ncY_INF)+sqrt(((0.5*ncX_INF-0.5*ncY_INF)**2)+(ncXY_INF**2)))

            if (cond109 .eq. 1) then
                alpha = 0.5*(ncX_INF+ncY_INF)
                alpha = alpha+sqrt(((0.5*ncX_INF-0.5*ncY_INF)**2)+(ncXY_INF**2))
                denum = 0.5*(ncX_INF+ncY_INF)-sqrt(((0.5*ncX_INF-0.5*ncY_INF)**2)+(ncXY_INF**2))
                alpha = alpha/denum
                alpha = abs(alpha)
                fcd2 = 0.85*fcd*(1+3.8*alpha)/((1+alpha)**2)
            else
                fcd2 = fcd
            end if

            fc = fcd2
            Residu = 0.5*(ncX_INF+ncY_INF)
            Residu = Residu-sqrt(((0.5*ncX_INF-0.5*ncY_INF)**2)+(ncXY_INF**2))
            Residu = Residu+tINF*fc

            if (abs(ncXY_INF) .ge. epsilon(ncXY_INF)) then
                AngleINF = atan((-tINF*fc-ncX_INF)/(ncXY_INF))
            elseif (abs(ncX_INF) .ge. abs(ncY_INF)) then
                AngleINF = pi/2.0
            else
                AngleINF = 0.d0
            end if

            if (AngleINF .ge. 0) then
                AngleINF = 90-AngleINF*r8rddg()
            else
                AngleINF = -90-AngleINF*r8rddg()
            end if

            nS_TOT = abs(nSX_SUP)+abs(nSX_INF)+abs(nSY_SUP)+abs(nSY_INF)

        else

            nSX_SUP = -1.d0
            nSX_INF = 0.d0
            nSY_SUP = -1.d0
            nSY_INF = 0.d0
            nS_TOT = -1.d0
            Residu = -1.d0
            AngleINF = 0.d0
            ncMAX_SUP = -1.d0
            ncMIN_SUP = -1.d0
            ncMAX_INF = -1.d0
            ncMIN_INF = -1.d0

        end if

    end subroutine solver_sandcas2

end module sand_solvers_module
