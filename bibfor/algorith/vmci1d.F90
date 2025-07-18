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

subroutine vmci1d(fami, kpg, ksp, imate, em, &
                  ep, sigm, deps, vim, option, &
                  materi, sigp, vip, dsde)
!
! person_in_charge: jean-luc.flejou at edf.fr
    implicit none
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: kpg, ksp, imate
    real(kind=8) :: ep, em, sigm, deps, sigp, dsde
    real(kind=8) :: vim(*), vip(*)
    character(len=16) :: option
    character(len=*) :: fami, materi
!
! --------------------------------------------------------------------------------------------------
!
!           PLASTICITE VON MISES CINEMATIQUE LINEAIRE EN 1D
!              FORTEMENT INSPIRE DE NM1DCI
!  IN
!        FAMI   : FAMILLE DU POINT DE GAUSS
!        KPG    : NUMERO DU POINT DE GAUSS
!        KSP    : NUMERO DU SOUS-POINT DE GAUSS / FIBRE
!        IMATE  : POINTEUR MATERIAU CODE
!        EM     : MODULE D YOUNG MOINS
!        EP     : MODULE D YOUNG PLUS
!        SIGM   : CONTRAINTE AU TEMPS MOINS
!        DEPS   : DEFORMATION TOTALE PLUS - DEFORMATION MOINS
!                       - INCREMENT DEFORMATION THERMIQUE
!        VIM    : VARIABLE INTERNES MOINS
!        OPTION : OPTION DE CALCUL
!  OUT
!        SIGP   : CONTRAINTES PLUS
!        VIP    : VARIABLE INTERNES PLUS
!        DSDE   : DSIG/DEPS
! --------------------------------------------------------------------------------------------------
!     Variables internes
!       icels : critère sigma
!       icelu : critère epsi
!       iepsq : déformation équivalente
!       iplas : indicateur plastique
!       idiss : dissipation plastique
!       iwthe : dissipation thermodynamique
!       i..m  : ecrouissage cinematique
! --------------------------------------------------------------------------------------------------
!   index des variables internes
!           'CRITSIG', 'CRITEPS', 'EPSPEQ', 'INDIPLAS', 'DISSIP', 'DISSTHER',
!           'XCINXX',  'XCINYY',  'XCINZZ', 'XCINXY', 'XCINXZ', 'XCINYZ',
    integer(kind=8), parameter :: icels = 1, icelu = 2, iepsq = 3, iplas = 4, idiss = 5, iwthe = 6
    integer(kind=8), parameter :: ixxm = 7
    integer(kind=8), parameter :: nbvari = 12
! --------------------------------------------------------------------------------------------------
    real(kind=8)        :: sigy, sieleq, sige, dp, etm, etp, xp, xm, hm, hp, sgels, epelu
    character(len=16)   :: valkm(3)
    integer(kind=8)             :: icodre(4)
    real(kind=8)        :: valres(4)
    character(len=16)   :: nomecl(4)
!
    data nomecl/'D_SIGM_EPSI', 'SY', 'SIGM_LIM', 'EPSI_LIM'/
! --------------------------------------------------------------------------------------------------
!   instant -
    call rcvalb(fami, kpg, ksp, '-', imate, materi, 'ECRO_LINE', 0, ' ', [0.d0], &
                1, nomecl, valres, icodre, 1)
    etm = valres(1)
    hm = em*etm/(em-etm)
!   instant +
    call rcvalb(fami, kpg, ksp, '+', imate, materi, 'ECRO_LINE', 0, ' ', [0.d0], &
                4, nomecl, valres, icodre, 1)
!   vérification que SIGM_LIM, EPSI_LIM sont présents
    if (icodre(3)+icodre(4) .ne. 0) then
        valkm(1) = 'VMIS_CINE_GC'
        valkm(2) = nomecl(3)
        valkm(3) = nomecl(4)
        call utmess('F', 'COMPOR1_76', nk=3, valk=valkm)
    end if
    etp = valres(1)
    sigy = valres(2)
    sgels = valres(3)
    epelu = valres(4)
!
    hp = ep*etp/(ep-etp)
    xm = vim(ixxm)
!
    sige = ep*(sigm/em+deps)-hp*xm/hm
    sieleq = abs(sige)
!
! --------------------------------------------------------------------------------------------------
!   calcul : EPSP, P , SIG
    if ((option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RAPH_MECA')) then
        vip(1:nbvari) = vim(1:nbvari)
        if (sieleq .le. sigy) then
            vip(iplas) = 0.d0
            dsde = ep
            dp = 0.d0
            xp = hp*xm/hm
            sigp = ep*(sigm/em+deps)
            vip(ixxm) = xp
            vip(icelu) = (sigm/em+deps)/epelu
            vip(iepsq) = (sigm/em+deps)
        else
            vip(iplas) = 1.d0
            dp = (sieleq-sigy)/(ep+hp)
            if (option .eq. 'FULL_MECA_ELAS') then
                dsde = ep
            else
                dsde = etp
            end if
            xp = hp*xm/hm+hp*dp*sige/sieleq
            sigp = xp+sigy*sige/sieleq
            vip(ixxm) = xp
            vip(icelu) = ((sigp-sigy)/etp+sigy/ep)/epelu
            vip(iepsq) = ((sigp-sigy)/etp+sigy/ep)
        end if
        vip(icels) = sigp/sgels
!       dissipation thermodynamique
        vip(iwthe) = vim(iwthe)+sigy*dp
!       dissipation irréversible
        vip(idiss) = vim(idiss)+(dsde*deps-(sigp-sigm))*deps/2.0d0
    end if
    if (option(1:10) .eq. 'RIGI_MECA_') then
        if ((vim(iplas) .lt. 0.5d0) .or. (option .eq. 'RIGI_MECA_ELAS')) then
            dsde = ep
        else
            dsde = etp
        end if
    end if
end subroutine
