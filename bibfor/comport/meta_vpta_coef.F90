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
subroutine meta_vpta_coef(metaRela, metaGlob, &
                          lgpg, fami, kpg, j_mater, &
                          l_temp, temp, meta_type, nb_phasis, phas_prev, &
                          phas_curr, zcold_curr, young, deuxmu, coef, &
                          trans)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/metaGetMechanism.h"
#include "asterfort/metaGetParaPlasTransf.h"
#include "asterfort/metaGetParaVisc.h"
#include "asterfort/metaGetParaMixture.h"
#include "asterfort/metaGetParaHardTrac.h"
#include "asterfort/metaGetParaHardLine.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: metaRela, metaGlob
    integer(kind=8), intent(in) :: lgpg
    character(len=4), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: j_mater
    aster_logical, intent(in) :: l_temp
    real(kind=8), intent(in) :: temp
    integer(kind=8), intent(in) :: meta_type
    integer(kind=8), intent(in) :: nb_phasis
    real(kind=8), intent(in) :: phas_prev(*)
    real(kind=8), intent(in) :: phas_curr(*)
    real(kind=8), intent(in) :: zcold_curr
    real(kind=8), intent(in) :: young
    real(kind=8), intent(in) :: deuxmu
    real(kind=8), intent(out) :: coef
    real(kind=8), intent(out) :: trans
!
! --------------------------------------------------------------------------------------------------
!
! Metallurgy - Comportment
!
! Compute coefficient for second member
! Effect of command variables phasis variation on transformation plasticity
!
! --------------------------------------------------------------------------------------------------
!
! In  metaRela      : behaviour for each phase
! In  metaGlob      : global behaviour
! In  lgpg          : length of integration point
! In  fami          : integration point type
! In  kpg           : integration point number
! In  j_mater       : coded material address
! In  l_temp        : .true. if temperature command variable is affected
! In  temp          : temperature
! In  meta_type     : type of metallurgy
! In  nb_phasis     : total number of phasis (cold and hot)
! In  phas_prev     : previous phasis
! In  phas_curr     : current phasis
! In  zcold_curr    : sum of cold phasis
! In  young         : Young modulus
! In  deuxmu        : e/(1.d0+nu)
! Out coef          : coefficient for command variable second member
! Out trans         : transformation for command variable second member
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: j_vari
    integer(kind=8) :: i_phasis, i_phasis_c, ksp, nb_phasis_c
    real(kind=8) :: epsp(5), h0(5)
    real(kind=8) :: kpt(4), fpt(4)
    real(kind=8) :: eta(5), n(5), unsurn(5), c(5), m(5)
    real(kind=8) :: rprim, deltaz(8), fmel, coef_hard
    aster_logical :: l_visc, l_elas, l_plas_tran, l_hard_line
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_phasis .le. 5)
    nb_phasis_c = nb_phasis-1
    trans = 0.d0
    coef = 1.d0
    ksp = 1
    fpt(:) = 0.d0
    kpt(:) = 0.d0
    eta(:) = 0.d0
    n(:) = 0.d0
    unsurn(:) = 0.d0
    c(:) = 0.d0
    m(:) = 0.d0
    epsp(:) = 0.d0
    h0(:) = 0.d0
!
! - Cumulated plastic strain
!
    call jevech('PVARIPR', 'L', j_vari)
    do i_phasis = 1, nb_phasis
        epsp(i_phasis) = zr(j_vari+lgpg*(kpg-1)-1+i_phasis)
        deltaz(i_phasis) = (phas_curr(i_phasis)-phas_prev(i_phasis))
    end do
!
! - Is elastic ?
!
    if (meta_type .eq. META_STEEL) then
        l_elas = zr(j_vari+lgpg*(kpg-1)+(nb_phasis+1)) .lt. 0.5d0
    elseif (meta_type .eq. META_ZIRC) then
        l_elas = zr(j_vari+lgpg*(kpg-1)+nb_phasis) .lt. 0.5d0
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Mechanisms of comportment law
!
    call metaGetMechanism(metaRela, metaGlob, &
                          l_visc=l_visc, &
                          l_hard_line=l_hard_line, l_plas_tran=l_plas_tran)
!
! - Transformation plasticity parameters
!
    if (l_plas_tran) then
        call metaGetParaPlasTransf('+', fami, kpg, ksp, j_mater, &
                                   meta_type, nb_phasis, deltaz, zcold_curr, &
                                   kpt, fpt)
    end if
!
! - Visco-plasticity parameters
!
    if (l_visc) then
        call metaGetParaVisc('+', fami, kpg, ksp, j_mater, &
                             meta_type, nb_phasis, eta, n, unsurn, &
                             c, m)
    end if
!
! - Compute Sum(iphase) [kpt * fpt] on cold phasis
!
    trans = 0.d0
    do i_phasis_c = 1, nb_phasis_c
        if (deltaz(i_phasis_c) .gt. 0) then
            trans = trans+kpt(i_phasis_c)*fpt(i_phasis_c)*deltaz(i_phasis_c)
        end if
    end do
!
! - Compute coefficient
!
    if (l_elas .or. l_visc) then
        coef = 1.d0
    else
!
! ----- Mixing law: yield
!
        call metaGetParaMixture('+', fami, kpg, ksp, j_mater, &
                                l_visc, meta_type, nb_phasis, zcold_curr, fmel=fmel)
!
! ----- Get point on hardening curve
!
        if (l_hard_line) then
            coef_hard = 1.d0
            call metaGetParaHardLine('+', fami, kpg, ksp, j_mater, &
                                     meta_type, nb_phasis, &
                                     young, coef_hard, h0)
        else
            call metaGetParaHardTrac(j_mater, meta_type, nb_phasis, &
                                     l_temp, temp, &
                                     epsp, h0)
        end if
!
! ----- Compute coefficient
!
        rprim = 0.d0
        if (zcold_curr .gt. 0.d0) then
            do i_phasis_c = 1, nb_phasis_c
                rprim = rprim+phas_curr(i_phasis_c)*h0(i_phasis_c)
            end do
            rprim = rprim/zcold_curr
        else
            rprim = 0.d0
        end if
        rprim = (1.d0-fmel)*h0(nb_phasis)+fmel*rprim
        coef = 1.d0-(1.5d0*deuxmu)/(1.5d0*deuxmu+rprim)
    end if
!
end subroutine
