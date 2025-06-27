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
subroutine btdmsn(ind, nb1, intsn, npgsr, xr, &
                  btdm, btdf, btds, btild)
!
    implicit none
!
    integer(kind=8) :: nb1, intsn, npgsr
    real(kind=8) :: xr(*), btdm1, btds1
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42), btdf(3, 42), btild(5, 42)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, ind, j, k, l
!-----------------------------------------------------------------------
    l = 702
!
!     CALCUL DE BTILDMN, BTILDSN AUX PTS DE GAUSS NORMAL
!            M=MEMBRANE, S=CISAILLEMENT, N=NORMAL
!
!        BTILDMN = SOMME MKBARRE BTILDMR   OU  K=1,NPGSR  R=REDUIT
!        BTILDSN = SOMME MKBARRE BTILDSR   OU  K=1,NPGSR  R=REDUIT
!        (AUX PTS DE GAUSS NORMAL)
!
!     MKBARRE = FONCTIONS DE FORME ASSOCIEES AUX PTS DE GAUSS REDUITS
!
    i1 = l+4*(intsn-1)
!
    if (ind .eq. 0) then
!
!     INTEGRATION UNIFORME (ICI REDUITE)
!     BTILD =  BTDF + BTDM  : BTDF, BTDM , BTDS
!              BTDS           OBTENUES PAR INTEGRATION REDUITE
!
        do i = 1, 3
            do j = 1, 5*nb1+2
                btild(i, j) = btdf(i, j)+btdm(intsn, i, j)
            end do
        end do
!
        do i = 1, 2
            do j = 1, 5*nb1+2
                btild(3+i, j) = btds(intsn, i, j)
            end do
        end do
!
    else if (ind .eq. 1) then
!
!     INTEGRATION SELECTIVE
!     BTILD =  BTDF + BTDM  : BTDF, BTDM , BTDS
!              BTDS           BTDM, BTDS INTEGRATION REDUITE
!
!
        do j = 1, 5*nb1+2
            do i = 1, 2
                btdm1 = 0.d0
                btds1 = 0.d0
                do k = 1, npgsr
                    btdm1 = btdm1+xr(i1+k)*btdm(k, i, j)
                    btds1 = btds1+xr(i1+k)*btds(k, i, j)
                end do
!
!                               BTILDMN + BTILDFN
!     CONSTRUCTION DE BTILD =
!                               BTILDSN
!
!
                btild(i, j) = btdm1+btdf(i, j)
                btild(i+3, j) = btds1
            end do
            btdm1 = 0.d0
            do k = 1, npgsr
                btdm1 = btdm1+xr(i1+k)*btdm(k, 3, j)
            end do
!
!                               BTILDMN + BTILDFN
!     CONSTRUCTION DE BTILD =
!                               BTILDSN
!
!
            btild(3, j) = btdm1+btdf(3, j)
        end do
!
    end if
!
end subroutine
