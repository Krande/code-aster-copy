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

subroutine verift(fami, kpg, ksp, poum, j_mater, &
                  materi_, iret_, epsth_, epsth_anis_, epsth_meta_, &
                  temp_prev_, temp_curr_, temp_refe_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/get_elasth_para.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/metaGetType.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecael.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: j_mater
    character(len=8), optional, intent(in) :: materi_
    integer(kind=8), optional, intent(out) :: iret_
    real(kind=8), optional, intent(out) :: epsth_
    real(kind=8), optional, intent(out) :: epsth_anis_(3)
    real(kind=8), optional, intent(out) :: epsth_meta_
    real(kind=8), optional, intent(out) :: temp_prev_
    real(kind=8), optional, intent(out) :: temp_curr_
    real(kind=8), optional, intent(out) :: temp_refe_
!
! --------------------------------------------------------------------------------------------------
!
! Compute thermic dilatation
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  j_mater      : coded material address
! In  poum         : parameters evaluation
!                     '-' for previous temperature
!                     '+' for current temperature
!                     'T' for current and previous temperature => epsth is increment
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  materi       : name of material if multi-material Gauss point (PMF)
! In  iret         : 0 if temperature defined
!                    1 if not
! Out epsth_       : thermic dilatation
! Out epsth_anis_  : non-isotropic thermic dilatation
! Out epsth_meta_  : multiphasic thermic dilatation
! Out temp_prev    : previous temperature
! Out temp_curr    : current temperature
! Out temp_refe    : reference temperature
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elem_name, materi
    integer(kind=8) :: iret_temp_prev, iret_temp_curr, iret_temp, iret_temp_refe
    real(kind=8) :: temp_prev, temp_refe, temp_curr
    real(kind=8) :: epsth, epsth_anis(3), epsth_meta
    real(kind=8) :: alpha_p(2)
    real(kind=8) :: alpha_l_p, alpha_t_p, alpha_n_p
    real(kind=8) :: alpha_c(2)
    real(kind=8) :: alpha_l_c, alpha_t_c, alpha_n_c
    real(kind=8) :: zcold_p, zhot_p, zcold_c, zhot_c, epsth_meta_h, epsth_meta_c
    real(kind=8) :: z_h_r, deps_ch_tref
    integer(kind=8) :: iadzi, iazk24
    integer(kind=8) :: elas_id, iret_cmp, icompo, meta_type, nb_phasis
    character(len=16) :: elas_keyword, rela_comp
!
! --------------------------------------------------------------------------------------------------
!
    materi = ' '
    if (present(materi_)) then
        materi = materi_
    end if
!
    iret_temp = 0
    iret_temp_prev = 0
    iret_temp_curr = 0
    iret_temp_refe = 0
    epsth = 0.d0
    epsth_anis(1:3) = 0.d0
    epsth_meta = 0.d0
    temp_refe = 0.d0
    temp_curr = 0.d0
    temp_prev = 0.d0
!
! - No temperature -> thermic strain is zero
!
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                ksp, temp_curr, iret_temp)
    if (iret_temp .ne. 0) then
        goto 999
    end if
!
! - Get reference temperature
!
    call rcvarc(' ', 'TEMP', 'REF', fami, kpg, &
                ksp, temp_refe, iret_temp_refe)
    if (iret_temp_refe .eq. 1) then
        call tecael(iadzi, iazk24)
        elem_name = zk24(iazk24-1+3) (1:8)
        call utmess('F', 'COMPOR5_8', sk=elem_name)
    end if
!
! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
!
    call get_elas_id(j_mater, elas_id, elas_keyword)
!
! - Get temperatures
!
    if (poum .eq. 'T' .or. poum .eq. '-') then
        call rcvarc(' ', 'TEMP', '-', fami, kpg, &
                    ksp, temp_prev, iret_temp_prev)
    end if
    if (poum .eq. 'T' .or. poum .eq. '+') then
        call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                    ksp, temp_curr, iret_temp_curr)
    end if
!
! - Get elastic parameters for thermic dilatation
!
    if (poum .eq. 'T' .or. poum .eq. '-') then
        if (iret_temp_prev .eq. 0) then
            call get_elasth_para(fami, j_mater, '-', kpg, ksp, &
                                 elas_id, elas_keyword, materi_=materi, &
                                 alpha=alpha_p, &
                                 alpha_l=alpha_l_p, &
                                 alpha_t=alpha_t_p, &
                                 alpha_n=alpha_n_p)
        end if
    end if
    if (poum .eq. 'T' .or. poum .eq. '+') then
        if (iret_temp_curr .eq. 0) then
            call get_elasth_para(fami, j_mater, '+', kpg, ksp, &
                                 elas_id, elas_keyword, materi_=materi, &
                                 alpha=alpha_c, &
                                 alpha_l=alpha_l_c, &
                                 alpha_t=alpha_t_c, &
                                 alpha_n=alpha_n_c)
        end if
    end if
!
! - Check non-isotropic material
!
    if (elas_id .ne. 1) then
        if (.not. present(epsth_anis_)) then
            call tecael(iadzi, iazk24)
            elem_name = zk24(iazk24-1+3) (1:8)
            call utmess('F', 'COMPOR5_9', sk=elem_name)
        end if
    end if
!
! - Check metallurgical material
!
    if (elas_keyword .eq. 'ELAS_META') then
        if (.not. present(epsth_meta_)) then
            call tecael(iadzi, iazk24)
            elem_name = zk24(iazk24-1+3) (1:8)
            call utmess('F', 'COMPOR5_10', sk=elem_name)
        end if
!
        call metaGetType(meta_type, nb_phasis)
!
        if (poum .eq. 'T' .or. poum .eq. '-') then
            if (iret_temp_prev .eq. 0) then
                call metaGetPhase(fami, '-', kpg, ksp, meta_type, nb_phasis, &
                                  zcold_=zcold_p, &
                                  zhot_=zhot_p)
            end if
        end if
!
        if (poum .eq. 'T' .or. poum .eq. '+') then
            if (iret_temp_prev .eq. 0) then
                call metaGetPhase(fami, '+', kpg, ksp, meta_type, nb_phasis, &
                                  zcold_=zcold_c, &
                                  zhot_=zhot_c)
            end if
        end if
!   - Check behavior name (used only in metallurgical case)
        call tecach('NNO', 'PCOMPOR', 'L', iret_cmp, iad=icompo)
        if (iret_cmp .eq. 0) then
            rela_comp = zk16(icompo)
        else
            rela_comp = 'Unknown'
        end if
    end if
!
! - Compute thermic strain
!
    if (poum .eq. 'T') then
        if (iret_temp_prev+iret_temp_curr .eq. 0) then
            if (elas_id .eq. 1) then
                if (elas_keyword .eq. 'ELAS_META') then
                    epsth_meta_h = zhot_c*alpha_c(1)*(temp_curr-temp_refe)- &
                                   zhot_p*alpha_p(1)*(temp_prev-temp_refe)
                    epsth_meta_c = zcold_c*alpha_c(2)*(temp_curr-temp_refe)- &
                                   zcold_p*alpha_p(2)*(temp_prev-temp_refe)
                    if (rela_comp .ne. 'META_LEMA_ANI') then
                        call get_elasth_para(fami, j_mater, '+', kpg, ksp, &
                                             elas_id, elas_keyword, materi_=materi, &
                                             z_h_r_=z_h_r, deps_ch_tref_=deps_ch_tref)
                        epsth_meta_h = epsth_meta_h+(1-z_h_r)*deps_ch_tref*(zhot_p-zhot_c)
                        epsth_meta_c = epsth_meta_c+z_h_r*deps_ch_tref*(zcold_c-zcold_p)
                    end if
                    epsth_meta = epsth_meta_h+epsth_meta_c
                else
                    epsth = alpha_c(1)*(temp_curr-temp_refe)-alpha_p(1)*(temp_prev-temp_refe)
                end if
            elseif (elas_id .eq. 2) then
                epsth_anis(1) = alpha_l_c*(temp_curr-temp_refe)-alpha_l_p*(temp_prev-temp_refe)
                epsth_anis(2) = alpha_t_c*(temp_curr-temp_refe)-alpha_t_p*(temp_prev-temp_refe)
                epsth_anis(3) = alpha_n_c*(temp_curr-temp_refe)-alpha_n_p*(temp_prev-temp_refe)
            elseif (elas_id .eq. 3) then
                epsth_anis(1) = alpha_l_c*(temp_curr-temp_refe)-alpha_l_p*(temp_prev-temp_refe)
                epsth_anis(2) = alpha_n_c*(temp_curr-temp_refe)-alpha_n_p*(temp_prev-temp_refe)
            else
                ASSERT(.false.)
            end if
        end if
    else if (poum .eq. '-') then
        if (iret_temp_prev .eq. 0) then
            if (elas_id .eq. 1) then
                if (elas_keyword .eq. 'ELAS_META') then
                    epsth_meta_h = zhot_p*alpha_p(1)*(temp_prev-temp_refe)
                    epsth_meta_c = zcold_p*alpha_p(2)*(temp_prev-temp_refe)
                    if (rela_comp .ne. 'META_LEMA_ANI') then
                        call get_elasth_para(fami, j_mater, '+', kpg, ksp, &
                                             elas_id, elas_keyword, materi_=materi, &
                                             z_h_r_=z_h_r, deps_ch_tref_=deps_ch_tref)
                        epsth_meta_h = epsth_meta_h-zhot_p*(1-z_h_r)*deps_ch_tref
                        epsth_meta_c = epsth_meta_c+zcold_p*z_h_r*deps_ch_tref
                    end if
                    epsth_meta = epsth_meta_h+epsth_meta_c
                else
                    epsth = alpha_p(1)*(temp_prev-temp_refe)
                end if
            elseif (elas_id .eq. 2) then
                epsth_anis(1) = alpha_l_p*(temp_prev-temp_refe)
                epsth_anis(2) = alpha_t_p*(temp_prev-temp_refe)
                epsth_anis(3) = alpha_n_p*(temp_prev-temp_refe)
            elseif (elas_id .eq. 3) then
                epsth_anis(1) = alpha_l_p*(temp_prev-temp_refe)
                epsth_anis(2) = alpha_n_p*(temp_prev-temp_refe)
            else
                ASSERT(.false.)
            end if
        end if
    else if (poum .eq. '+') then
        if (iret_temp_curr .eq. 0) then
            if (elas_id .eq. 1) then
                if (elas_keyword .eq. 'ELAS_META') then
                    epsth_meta_h = zhot_c*alpha_c(1)*(temp_curr-temp_refe)
                    epsth_meta_c = zcold_c*alpha_c(2)*(temp_curr-temp_refe)
                    if (rela_comp .ne. 'META_LEMA_ANI') then
                        call get_elasth_para(fami, j_mater, '+', kpg, ksp, &
                                             elas_id, elas_keyword, materi_=materi, &
                                             z_h_r_=z_h_r, deps_ch_tref_=deps_ch_tref)
                        epsth_meta_h = epsth_meta_h-zhot_c*(1-z_h_r)*deps_ch_tref
                        epsth_meta_c = epsth_meta_c+zcold_c*z_h_r*deps_ch_tref
                    end if
                    epsth_meta = epsth_meta_h+epsth_meta_c
                else
                    epsth = alpha_c(1)*(temp_curr-temp_refe)
                end if
            elseif (elas_id .eq. 2) then
                epsth_anis(1) = alpha_l_c*(temp_curr-temp_refe)
                epsth_anis(2) = alpha_t_c*(temp_curr-temp_refe)
                epsth_anis(3) = alpha_n_c*(temp_curr-temp_refe)
            elseif (elas_id .eq. 3) then
                epsth_anis(1) = alpha_l_c*(temp_curr-temp_refe)
                epsth_anis(2) = alpha_n_c*(temp_curr-temp_refe)
            else
                ASSERT(.false.)
            end if
        end if
    else
        ASSERT(.false.)
    end if
!
999 continue
!
! - Output temperature
!
    if (present(temp_refe_)) then
        temp_refe_ = temp_refe
    end if
    if (present(temp_prev_)) then
        temp_prev_ = temp_prev
    end if
    if (present(temp_curr_)) then
        temp_curr_ = temp_curr
    end if
!
! - Output strains
!
    if (present(epsth_meta_)) then
        epsth_meta_ = epsth_meta
    end if
    if (present(epsth_anis_)) then
        epsth_anis_(1:3) = epsth_anis(1:3)
    end if
    if (present(epsth_)) then
        epsth_ = epsth
    end if
!
! - Output error
!
    if (present(iret_)) then
        iret_ = 0
        if ((iret_temp_prev+iret_temp_curr) .ne. 0) then
            iret_ = 1
        end if
        if (iret_temp .ne. 0) then
            iret_ = 1
        end if
    end if
!
end subroutine
