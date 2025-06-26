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
subroutine rc32s0b(seis, sig, trescamax)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rctres.h"
    real(kind=8) :: sig(6), trescamax
    real(kind=8) :: seis(6*12)
!
!
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3200
!     CALCUL AVEC TOUTES LES POSSIBILITES DE SIGNE
!     POUR LES CMP DE SEISME
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: i0, i1, e0(2), i2, i3, i4, i5, i6, i7
    integer(kind=8) :: i8, i9, i10, i11, i12, j, jinfois, typseis
    real(kind=8) :: st(6), tresca
!
! DEB ------------------------------------------------------------------
!
    trescamax = 0.d0
    do i0 = 1, 2
        i1 = 2*(i0-2)+1
        e0(i0) = i1
    end do
!
    call jeveuo('&&RC3200.SEIS_INFOI', 'L', jinfois)
    typseis = zi(jinfois+3)
!
    if (typseis .eq. 1) then
!
        do i1 = 1, 2
            do i2 = 1, 2
                do i3 = 1, 2
                    do i4 = 1, 2
                        do i5 = 1, 2
                            do i6 = 1, 2
                                do i7 = 1, 2
                                    do i8 = 1, 2
                                        do i9 = 1, 2
                                            do i10 = 1, 2
                                                do i11 = 1, 2
                                                    do i12 = 1, 2
                                                        do j = 1, 6
!
                                                            st(j) = sig(j)+e0(i1)*seis(j)+e0(i2)&
                                                                    &*seis(6+j)+e0(i3)*seis(2*6&
                                                                    &+j)+e0(i4)*seis(3*6+j)+e0&
                                                                    &(i5)*seis(4*6+j)+e0(i6)*sei&
                                                                    &s(5*6+j)+e0(i7)*seis(6*6+j&
                                                                    &)+e0(i8)*seis(7*6+j)+e0(i&
                                                                    &9)*seis(8*6+j)+e0(i10)*seis&
                                                                    &(9*6+j)+e0(i11)*seis(10*6+&
                                                                    &j)+e0(i12)*seis(11*6+j)
!
                                                        end do
                                                        call rctres(st, tresca)
                                                        trescamax = max(tresca, trescamax)
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
!
    else
!
        do i1 = 1, 2
            do i2 = 1, 2
                do i3 = 1, 2
                    do i4 = 1, 2
                        do i5 = 1, 2
                            do i6 = 1, 2
                                do j = 1, 6
!
                                    st(j) = sig(j)+e0(i1)*seis(j)+e0(i2)*seis(6+j)+e0(i3)*seis(2&
                                            &*6+j)+e0(i4)*seis(3*6+j)+e0(i5)*seis(4*6+j)+e0(i6&
                                            &)*seis(5*6+j)
!
                                end do
                                call rctres(st, tresca)
                                trescamax = max(tresca, trescamax)
!
                            end do
                        end do
                    end do
                end do
            end do
        end do
!
    end if
!
end subroutine
