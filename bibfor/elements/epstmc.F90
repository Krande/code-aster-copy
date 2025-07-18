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

subroutine epstmc(fami, ndim, instan, poum, kpg, &
                  ksp, angl_naut, j_mater, option, &
                  epsi_varc)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/calc_epth_elga.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
!
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: ndim
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: j_mater
    real(kind=8), intent(in) :: angl_naut(3)
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: instan
    real(kind=8), intent(out) :: epsi_varc(6)
!
! --------------------------------------------------------------------------------------------------
!
! Compute variable commands strains (themrics, drying, etc.)
!
! For isoparametric elements
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  ndim         : dimension of space
! In  poum         : parameters evaluation
!                     '-' for previous temperature
!                     '+' for current temperature
!                     'T' for current and previous temperature
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  j_mater      : coded material address
! In  angl_naut    : nautical angles (for non-isotropic materials)
! In  instan       : current time
! In  option       : name of option to compute
! Out epsi_varc    : command variables strains
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbres = 3
    integer(kind=8) :: icodre(nbres)
    character(len=16) :: nomres(nbres)
    real(kind=8) :: valres(nbres)
!
    integer(kind=8) :: nbv, elas_id, nbpar
    real(kind=8) :: biot, e
    character(len=8) :: nompar
    character(len=32) :: phenom
    real(kind=8) :: valpar, bendog, kdessi
    real(kind=8) :: troisk, nu
    real(kind=8) :: hydr, sech, sref, ptot
    integer(kind=8) :: k, iret
    character(len=16) :: elas_keyword
    character(len=6), parameter :: epsa(6) = (/'EPSAXX', 'EPSAYY', 'EPSAZZ', &
                                               'EPSAXY', 'EPSAXZ', 'EPSAYZ'/)
!
! --------------------------------------------------------------------------------------------------
!
    if (instan .eq. r8vide()) then
        nbpar = 0
    else
        nbpar = 1
        nompar = 'INST'
        valpar = instan
    end if
    epsi_varc(1:6) = 0.d0
    biot = 0.d0
    bendog = 0.d0
    kdessi = 0.d0
!
! - Get command variables
!
    call rcvarc(' ', 'HYDR', poum, fami, kpg, &
                ksp, hydr, iret)
    if (iret .eq. 1) hydr = 0.d0
    call rcvarc(' ', 'SECH', poum, fami, kpg, &
                ksp, sech, iret)
    if (iret .eq. 1) sech = 0.d0
    call rcvarc(' ', 'PTOT', poum, fami, kpg, &
                ksp, ptot, iret)
    if (iret .eq. 1) ptot = 0.d0
    call rcvarc(' ', 'SECH', 'REF', fami, 1, &
                1, sref, iret)
    if (iret .eq. 1) sref = 0.d0
!
    if (option(11:14) .eq. 'HYDR') then
!
! ----- Hydric strains
!
        if (hydr .ne. 0.d0) then
            call get_elas_id(j_mater, elas_id, elas_keyword)
            nomres(1) = 'B_ENDOGE'
            nbv = 1
            call rcvalb(fami, kpg, ksp, poum, j_mater, &
                        ' ', elas_keyword, nbpar, nompar, [valpar], &
                        nbv, nomres, valres, icodre, 0)
            if (icodre(1) .eq. 0) then
                bendog = valres(1)
                epsi_varc(1) = -bendog*hydr
                epsi_varc(2) = -bendog*hydr
                epsi_varc(3) = -bendog*hydr
            else
                call utmess('I', 'COMPOR5_12')
            end if
        end if
    else if (option(11:14) .eq. 'PTOT') then
!
! ----- Fluid pressure strain
!
        if (ptot .ne. 0.d0) then
            phenom = 'THM_DIFFU'
            nomres(1) = 'BIOT_COEF'
            nbv = 1
            call rcvalb(fami, kpg, ksp, poum, j_mater, &
                        ' ', phenom, nbpar, nompar, [valpar], &
                        nbv, nomres, valres, icodre, 0)
            if (icodre(1) .eq. 0) then
                biot = valres(1)
            else
                biot = 0.d0
                call utmess('I', 'COMPOR5_13')
            end if
            call get_elas_id(j_mater, elas_id, elas_keyword)
            call get_elas_para(fami, j_mater, poum, kpg, ksp, &
                               elas_id, elas_keyword, &
                               time=instan, &
                               e_=e, nu_=nu)
            troisk = e/(1.d0-2.d0*nu)
            epsi_varc(1) = biot/troisk*ptot
            epsi_varc(2) = epsi_varc(1)
            epsi_varc(3) = epsi_varc(1)
        end if
    else if (option(11:14) .eq. 'SECH') then
!
! ----- Drying strains
!
        call get_elas_id(j_mater, elas_id, elas_keyword)
        nomres(1) = 'K_DESSIC'
        nbv = 1
        call rcvalb(fami, kpg, ksp, poum, j_mater, &
                    ' ', elas_keyword, nbpar, nompar, [valpar], &
                    nbv, nomres, valres, icodre, 0)
        if (icodre(1) .eq. 0) then
            kdessi = valres(1)
        else
            kdessi = 0.d0
        end if
        epsi_varc(1) = -kdessi*(sref-sech)
        epsi_varc(2) = -kdessi*(sref-sech)
        epsi_varc(3) = -kdessi*(sref-sech)
    else if (option(11:14) .eq. 'EPSA') then
!
! ----- User strains for command variables
!
        do k = 1, 6
            call rcvarc(' ', epsa(k), poum, fami, kpg, &
                        ksp, epsi_varc(k), iret)
            if (iret .eq. 1) epsi_varc(k) = 0.d0
        end do
    else
!
! ----- Thermic strains
!
        call calc_epth_elga(fami, ndim, poum, kpg, ksp, &
                            j_mater, angl_naut, epsi_varc)
    end if
!
end subroutine
