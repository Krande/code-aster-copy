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

subroutine mcordo(dpstrs, pstrs, pstra, dirprj, &
                  edge, apex, codret)
!***********************************************************************
!
! OBJECT: RE-ORDER MATRICEX AND VECTORS ACOORDING TO:
!
!   . PSTRA1 > PSTRA2 > PSTRA3   => EDGE=0 APEX=0
!   . PSTRA1!= PSTRA2 = PSTRA3   => EDGE=1 APEX=0
!   . PSTRA1 = PSTRA2 = PSTRA3   => EDGE=0 APEX=1
! ----------------------------------------------------------------------
!
!     LOI DE COMPORTEMENT DE MOHR-COULOMB
!
! IN  DPSTRS  : D Y_PRIN / D X_PRIN (3,3)-MATRIX
! IN  PSTRS   : VALEURS PRINCIPALES DE SIGMA (3)
! IN  PSTRA   : VALEURS PRINCIPALES DE EPSILON (3)
! IN  DIRPRJ  : DIRECTIONS PRINCIPALES
!
! OUT EDGE    : MATRICE TANGENTE COHERENTE REACTUALISEE
! OUT APEX    : MATRICE TANGENTE COHERENTE REACTUALISEE
! OUT CODRET  : CODE RETOUR
!               = | 0: OK
!                 | 1: NOOK
!
! ----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
! ======================================================================
!
    real(kind=8) :: dpstrs(3, 3)
    real(kind=8) :: pstrs(3)
    real(kind=8) :: pstra(3)
    real(kind=8) :: dirprj(3, 3)
    real(kind=8) :: edge
    real(kind=8) :: apex
    integer(kind=8) :: codret
!
! Declaration of integer type variables
    integer(kind=8) :: i, j, mdim, ndim
!
! Declaration of integer type variables
!     aster_logical :: epflag
    aster_logical :: lorder
!
    parameter(mdim=3, ndim=6)
!
! Declaration of vector and matrix type variables
    integer(kind=8) :: iorder(mdim)
    real(kind=8) :: pst1, pst2, pst3, refe, tbidon(ndim, mdim), vbidon(mdim), r0, r1, r2
    real(kind=8) :: r3, r4, small, tol, dmax1, pstmin, pstmax
!
!
    data r0, r1, r2, r3, r4, small, tol/&
     &    0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 1.d-06, 1.d-10/
!
! Declaration of Common space variables
!     common / debug / epflag
!
    lorder = .false.
    pst1 = pstra(1)
    pst2 = pstra(2)
    pst3 = pstra(3)
    refe = dmax1(abs(pst1), abs(pst2), abs(pst3))*small
! Re-ordering principal components if two repeated eigenvalues
! such that x1!=x2=x3
    if (abs(pst1-pst3) .lt. refe) then
        iorder(1) = 2
        iorder(2) = 1
        iorder(3) = 3
        lorder = .true.
        edge = r1
        apex = r0
    end if
    if (abs(pst1-pst2) .lt. refe .and. .not. lorder) then
        iorder(1) = 3
        iorder(2) = 1
        iorder(3) = 2
        lorder = .true.
        edge = r1
        apex = r0
    else if (abs(pst2-pst3) .lt. refe .and. .not. lorder) then
        lorder = .false.
        edge = r1
        apex = r0
    elseif (abs(pst2-pst3) .lt. refe .and. abs(pst1-pst2) .lt. refe .and. &
            lorder) then
        lorder = .false.
        edge = r0
        apex = r1
    else if (.not. lorder) then
        iorder(1) = 1
        iorder(3) = 1
        pstmax = pst1
        pstmin = pst1
        do i = 2, 3
            if (pstra(i) .ge. pstmax) then
                iorder(1) = i
                pstmax = pstra(i)
            end if
            if (pstra(i) .lt. pstmin) then
                iorder(3) = i
                pstmin = pstra(i)
            end if
        end do
        if (iorder(1) .ne. 1 .and. iorder(3) .ne. 1) iorder(2) = 1
        if (iorder(1) .ne. 2 .and. iorder(3) .ne. 2) iorder(2) = 2
        if (iorder(1) .ne. 3 .and. iorder(3) .ne. 3) iorder(2) = 3
        if (iorder(1) .ne. 1 .or. iorder(2) .ne. 2 .or. iorder(3) .ne. 3) lorder = .true.
        edge = r0
        apex = r0
    else
        codret = 1
        goto 999
    end if
!
!
! Re-order PSTRA and DPRSTS
    if (lorder) then
!
        do i = 1, mdim
            vbidon(i) = pstra(iorder(i))
            do j = 1, mdim
                tbidon(i, j) = dpstrs(iorder(i), iorder(j))
            end do
        end do
!
        do i = 1, mdim
            pstra(i) = vbidon(i)
            do j = 1, mdim
                dpstrs(i, j) = tbidon(i, j)
            end do
        end do
!
! Re-order PSTRS and DIRPRJ
        do i = 1, mdim
            vbidon(i) = pstrs(iorder(i))
            do j = 1, mdim
                tbidon(j, i) = dirprj(j, iorder(i))
            end do
        end do
!
        do i = 1, mdim
            pstrs(i) = vbidon(i)
            do j = 1, mdim
                dirprj(j, i) = tbidon(j, i)
            end do
        end do
!
    end if
999 continue
end subroutine
