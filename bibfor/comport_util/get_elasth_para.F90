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

subroutine get_elasth_para(fami, j_mater, poum, ipg, ispg, &
                           elas_type, elas_keyword, materi_, temp_vale_, &
                           alpha, alpha_l, alpha_t, alpha_n, &
                           z_h_r_, deps_ch_tref_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
!
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: j_mater
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: ipg
    integer(kind=8), intent(in) :: ispg
    integer(kind=8), intent(in) :: elas_type
    character(len=16), intent(in) :: elas_keyword
    character(len=8), optional, intent(in) :: materi_
    real(kind=8), optional, intent(in) :: temp_vale_
    real(kind=8), optional, intent(out) :: alpha(2)
    real(kind=8), optional, intent(out) :: alpha_l
    real(kind=8), optional, intent(out) :: alpha_t
    real(kind=8), optional, intent(out) :: alpha_n
    real(kind=8), optional, intent(out) :: z_h_r_
    real(kind=8), optional, intent(out) :: deps_ch_tref_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility
!
! Get elastic parameters for thermic dilatation
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  j_mater      : coded material address
! In  poum         : '-' or '+' for parameters evaluation (previous or current temperature)
! In  ipg          : current point gauss
! In  ispg         : current "sous-point" gauss
! In  elas_type    : Type of elasticity
!                       1 - Isotropic
!                       2 - Orthotropic
!                       3 - Transverse isotropic
! In  elas_keyword : keyword factor linked to type of elasticity parameters
! In  materi       : name of material if multi-material Gauss point (PMF)
! In  temp_vale    : specifi temperature (example: mean temperature for structural elements)
! Out alpha        : thermic dilatation ratio (isotropic)
!                     if   META -> alpha(1) for hot phasis and alpha(2) for cold phasis
!                     else alpha(1) only
! Out alpha_l      : thermic dilatation ratio - Direction L (Orthotropic/Transverse isotropic)
! Out alpha_t      : thermic dilatation ratio - Direction T (Orthotropic)
! Out alpha_n      : thermic dilatation ratio - Direction N (Orthotropic/Transverse isotropic)
! Out z_h_r        : characterizes the reference metallurgical phase :
!                    z_h_r = 1 --> reference phasis = hot phase
!                    z_h_r = 0 --> reference phasis = cold phase
! Out deps_ch_tref : Compactness difference between the hot phase and the cold phase
!                    at the reference temperature
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbresm
    parameter(nbresm=4)
    integer(kind=8) :: icodre(nbresm)
    character(len=16) :: nomres(nbresm)
    real(kind=8) :: valres(nbresm)
!
    character(len=8) :: para_name, materi
    character(len=24) :: valk(3)
    real(kind=8) :: para_vale
    integer(kind=8) :: nbres, nb_para, i
    real(kind=8) :: alpha_c, alpha_f, alpha_a
    integer(kind=8) :: iadzi, iazk24
    real(kind=8) :: z_h_r
    real(kind=8) :: deps_ch_tref
!
! --------------------------------------------------------------------------------------------------
!
    nb_para = 0
    para_name = ' '
    para_vale = 0.d0
    materi = ' '
    if (present(materi_)) then
        materi = materi_
    end if
    if (present(temp_vale_)) then
        nb_para = 1
        para_vale = temp_vale_
        para_name = 'TEMP'
    end if
!
! - Get parameters
!
    if (elas_type .eq. 1) then
        if (elas_keyword .eq. 'ELAS_HYPER') then
            call utmess('F', 'COMPOR5_6')
        elseif (elas_keyword .eq. 'ELAS_META') then
            nbres = 4
            nomres(1) = 'C_ALPHA'
            nomres(2) = 'F_ALPHA'
            nomres(3) = 'PHASE_REFE'
            nomres(4) = 'EPSF_EPSC_TREF'
            call rcvalb(fami, ipg, ispg, poum, j_mater, &
                        materi, elas_keyword, nb_para, para_name, [para_vale], &
                        nbres, nomres, valres, icodre, 1)
            alpha_c = valres(1)
            alpha_f = valres(2)
            z_h_r = valres(3)
            deps_ch_tref = valres(4)
            if (present(alpha)) then
                alpha(1) = alpha_c
                alpha(2) = alpha_f
            end if
            if (present(z_h_r_)) then
                z_h_r_ = z_h_r
            end if
            if (present(deps_ch_tref_)) then
                deps_ch_tref_ = deps_ch_tref
            end if
        else
            nbres = 1
            nomres(1) = 'ALPHA'
            call rcvalb(fami, ipg, ispg, poum, j_mater, &
                        materi, elas_keyword, nb_para, para_name, [para_vale], &
                        nbres, nomres, valres, icodre, 1)
            alpha_a = valres(1)
            alpha(1) = alpha_a
            alpha(2) = 0.d0
        end if
    elseif (elas_type .eq. 2) then
        nbres = 3
        nomres(1) = 'ALPHA_L'
        nomres(2) = 'ALPHA_T'
        nomres(3) = 'ALPHA_N'
        call rcvalb(fami, ipg, ispg, poum, j_mater, &
                    materi, elas_keyword, nb_para, para_name, [para_vale], &
                    nbres, nomres, valres, icodre, 1)
        alpha_l = valres(1)
        alpha_t = valres(2)
        alpha_n = valres(3)
    elseif (elas_type .eq. 3) then
        nbres = 2
        nomres(1) = 'ALPHA_L'
        nomres(2) = 'ALPHA_N'
        call rcvalb(fami, ipg, ispg, poum, j_mater, &
                    materi, elas_keyword, nb_para, para_name, [para_vale], &
                    nbres, nomres, valres, icodre, 1)
        alpha_l = valres(1)
        alpha_n = valres(2)
    else
        ASSERT(.false.)
    end if
!
! - Test
!
    do i = 1, nbres
        if (icodre(i) .ne. 0) then
            call tecael(iadzi, iazk24)
            valk(1) = zk24(iazk24-1+3)
            valk(2) = 'TEMP'
            valk(3) = nomres(i)
            call utmess('F', 'COMPOR5_32', nk=3, valk=valk)
        end if
    end do
!
end subroutine
