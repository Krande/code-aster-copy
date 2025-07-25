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

subroutine glrc_change_rep_mat(vmp, vfp, dspdep, dsidep)
!
    implicit none
    real(kind=8) :: dspdep(6, *), dsidep(6, *), vmp(2, 2), vfp(2, 2)
    integer(kind=8) :: t(2, 2), i, j, k, l, n, m
! ----------------------------------------------------------------------
!
!   ROUTINE DE CHANGEMENT DE REPERE DE LA MATRICE TANGENTE POUR GLRC_LC
!
! IN:
!       VMP(2,2)    : MATRICE DES VECTEURS PROPRES DE MEMBRANE
!       VFP(2,2)    : MATRICE DES VECTEURS PROPRES DE FLEXION
!       DSPDEP(4,4) : MATRICE TANGENTE DANS LE REPERE PROPRE
!
! OUT:
!       DSIDEP(6,6) : MATRICE TANGENTE DANS LE REPERE GLOBAL
! ----------------------------------------------------------------------
    real(kind=8) :: rtemp, rtpf, rtmf, rtfm
! ----------------------------------------------------------------------
!
! MATRICE DE PASSAGE TENSEUR D'ORDRE 4 >> TENSEUR D'ORDRE 2
    t(1, 1) = 1
    t(1, 2) = 3
    t(2, 1) = 3
    t(2, 2) = 2
!
    do i = 1, 2
        do j = i, 2
            do k = 1, 2
                do l = 1, 2
                    if (t(i, j) .ge. t(k, l)) then
                        rtemp = 0.d0
                        rtpf = 0.d0
                        rtmf = 0.d0
                        rtfm = 0.d0
                        do m = 1, 2
                            do n = 1, 2
                                rtemp = rtemp+dspdep(n, m)*vmp(k, m)*vmp(i, n)*vmp(j, n)*vm&
                                        &p(l, m)
!
                                rtpf = rtpf+dspdep(n+3, m+3)*vfp(k, m)*vfp(i, n)*vfp(j, n)*&
                                       & vfp(l, m)
!
                                rtmf = rtmf+dspdep(n, m+3)*vfp(k, m)*vmp(i, n)*vmp(j, n)*vf&
                                       &p(l, m)
!
                                rtfm = rtfm+dspdep(n+3, m)*vmp(k, m)*vfp(i, n)*vfp(j, n)*vm&
                                       &p(l, m)
                            end do
                        end do
!
                        rtemp = rtemp+dspdep(3, 3)*vmp(i, 1)*vmp(j, 2)*vmp(k, 1)*vmp(l, 2)
                        rtemp = rtemp+dspdep(3, 3)*vmp(i, 2)*vmp(j, 1)*vmp(k, 2)*vmp(l, 1)
!
                        rtpf = rtpf+dspdep(6, 6)*vfp(i, 1)*vfp(j, 2)*vfp(k, 1)*vfp(l, 2)
                        rtpf = rtpf+dspdep(6, 6)*vfp(i, 2)*vfp(j, 1)*vfp(k, 2)*vfp(l, 1)
!
                        dsidep(t(i, j), t(k, l)) = dsidep(t(i, j), t(k, l))+rtemp
                        dsidep(t(i, j)+3, t(k, l)+3) = dsidep(t(i, j)+3, t(k, l)+3)+rtpf
                        dsidep(t(i, j)+3, t(k, l)) = dsidep(t(i, j)+3, t(k, l))+rtfm
                        if (t(i, j) .ne. t(k, l)) then
                            dsidep(t(k, l)+3, t(i, j)) = dsidep(t(k, l)+3, t(i, j))+rtmf
                        end if
                    end if
                end do
            end do
        end do
    end do
!
    dsidep(3, 3) = dsidep(3, 3)*0.5d0
    dsidep(6, 6) = dsidep(6, 6)*0.5d0
    dsidep(6, 3) = dsidep(6, 3)*0.5d0
!
    do i = 1, 6
        do j = i+1, 6
            dsidep(i, j) = dsidep(j, i)
        end do
    end do
!
end subroutine
