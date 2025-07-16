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
subroutine metaGetParaElas(poum, fami, kpg, ksp, j_mater, &
                           e_, deuxmu_, mu_, troisk_, &
                           em_, deuxmum_, mum_, troiskm_)
!
    implicit none
!
#include "asterfort/rcvalb.h"
#include "asterfort/Metallurgy_type.h"
!
    character(len=1), intent(in) :: poum
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp
    integer(kind=8), intent(in) :: j_mater
    real(kind=8), optional, intent(out) :: e_, deuxmu_, mu_, troisk_
    real(kind=8), optional, intent(out) :: em_, deuxmum_, mum_, troiskm_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get elastic parameters
!
! --------------------------------------------------------------------------------------------------
!
! In  poum         : '-' or '+' for parameters evaluation (previous or current)
! In  fami         : Gauss family for integration point rule
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  j_mater      : coded material address
! Out e            : young modulus for "poum"
! Out deuxmu       : Lame ratio for "poum"
! Out mu           : Lame ratio for "poum"
! Out troisk       : compressibility modulus for "poum"
! Out em           : young modulus for "-"
! Out deuxmum      : Lame ratio for "-"
! Out mum          : Lame ratio for "-"
! Out troiskm      : compressibility modulus for "-"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_resu_max = 2
    real(kind=8) :: resu_vale(nb_resu_max)
    integer(kind=8) :: codret(nb_resu_max)
    character(len=16) :: resu_name(nb_resu_max)
    integer(kind=8) :: nb_resu
    real(kind=8) :: em, num, deuxmum, mum, troiskm
    real(kind=8) :: e, nu, deuxmu, mu, troisk
!
! --------------------------------------------------------------------------------------------------
!
    nb_resu = 2
    resu_name(1) = 'E'
    resu_name(2) = 'NU'
!
    call rcvalb(fami, kpg, ksp, '-', j_mater, &
                ' ', 'ELAS_META', 0, ' ', [0.d0], &
                nb_resu, resu_name, resu_vale, codret, 2)
    em = resu_vale(1)
    num = resu_vale(2)
    deuxmum = em/(1.d0+num)
    mum = deuxmum/2.d0
    troiskm = em/(1.d0-2.d0*num)
    call rcvalb(fami, kpg, ksp, poum, j_mater, &
                ' ', 'ELAS_META', 0, ' ', [0.d0], &
                nb_resu, resu_name, resu_vale, codret, 2)
    e = resu_vale(1)
    nu = resu_vale(2)
    deuxmu = e/(1.d0+nu)
    mu = deuxmu/2.d0
    troisk = e/(1.d0-2.d0*nu)
!
    if (present(e_)) then
        e_ = e
    end if
    if (present(deuxmu_)) then
        deuxmu_ = deuxmu
    end if
    if (present(mu_)) then
        mu_ = mu
    end if
    if (present(troisk_)) then
        troisk_ = troisk
    end if
!
    if (present(em_)) then
        em_ = em
    end if
    if (present(deuxmum_)) then
        deuxmum_ = deuxmum
    end if
    if (present(mum_)) then
        mum_ = mum
    end if
    if (present(troiskm_)) then
        troiskm_ = troiskm
    end if
!
end subroutine
