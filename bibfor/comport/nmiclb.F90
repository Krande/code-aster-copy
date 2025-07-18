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
subroutine nmiclb(fami, kpg, ksp, option, rela_comp, &
                  imate, xlong0, aire, tmoins, tplus, &
                  dlong0, effnom, vim, effnop, vip, &
                  klv, fono, epsm, carcri, codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcimpl.h"
#include "asterfort/nm1dci.h"
#include "asterfort/nm1dco.h"
#include "asterfort/nm1dis.h"
#include "asterfort/relax_acier_cable.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imate, neq, nbt, kpg, ksp, codret
    parameter(neq=6, nbt=21)
!
    real(kind=8) :: xlong0, aire, tmoins, tplus, dlong0, carcri(CARCRI_SIZE), epsm
    real(kind=8) :: effnom, vim(*), effnop, vip(*), fono(neq), klv(nbt)
!
    character(len=16) :: rela_comp, option
    character(len=*) :: fami
!
! --------------------------------------------------------------------------------------------------
!
!
!    TRAITEMENT DE LA RELATION DE COMPORTEMENT -ELASTOPLASTICITE-
!    ECROUISSAGE ISOTROPE ET CINEMATIQUE- LINEAIRE - VON MISES-
!    POUR UN MODELE BARRE ELEMENT MECA_BARRE
!
!
! --------------------------------------------------------------------------------------------------
!
! IN  : IMATE : POINTEUR MATERIAU CODE
!       XLONG0 : LONGUEUR DE L'ELEMENT DE BARRE AU REPOS
!       aire   : SECTION DE LA BARRE
!       TMOINS : INSTANT PRECEDENT
!       TPLUS  : INSTANT COURANT
!       DLONG0 : INCREMENT D'ALLONGEMENT DE L'ELEMENT
!       EFFNOM : EFFORT NORMAL PRECEDENT
!       TREF   : TEMPERATURE DE REFERENCE
!       TEMPM  : TEMPERATURE IMPOSEE A L'INSTANT PRECEDENT
!       TEMPP  : TEMPERATURE IMPOSEE A L'INSTANT COURANT
!       OPTION : OPTION DEMANDEE (R_M_T,FULL OU RAPH_MECA)
! OUT : EFFNOP : CONTRAINTE A L'INSTANT ACTUEL
!       VIP    : VARIABLE INTERNE A L'INSTANT ACTUEL
!       FONO   : FORCES NODALES COURANTES
!       KLV    : MATRICE TANGENTE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)       :: codres(1)
    real(kind=8)  :: sigm, deps, depsth, depsm, em, ep
    real(kind=8)  :: sigp, xrig, val(1), dsde
    aster_logical :: isot, cine, elas, corr, implex, isotli, relax
!
! --------------------------------------------------------------------------------------------------
!
    elas = .false.
    isot = .false.
    cine = .false.
    corr = .false.
    implex = option .eq. 'RIGI_MECA_IMPLEX' .or. option .eq. 'RAPH_MECA_IMPLEX'
    isotli = .false.
    relax = .false.
    if (rela_comp .eq. 'ELAS') then
        elas = .true.
    else if ((rela_comp .eq. 'VMIS_ISOT_LINE') .or. &
             (rela_comp .eq. 'VMIS_ISOT_TRAC')) then
        isot = .true.
        if (rela_comp .eq. 'VMIS_ISOT_LINE') then
            isotli = .true.
        end if
    else if (rela_comp .eq. 'VMIS_CINE_LINE') then
        cine = .true.
    else if (rela_comp .eq. 'CORR_ACIER') then
        corr = .true.
    else if (rela_comp .eq. 'RELAX_ACIER') then
        relax = .true.
    end if

    if (implex) then
        if ((.not. elas) .and. (.not. isotli)) then
            call utmess('F', 'POUTRE0_49', sk=rela_comp)
        end if
    end if
!
    klv = 0.d0
    fono = 0.d0
!
!   Récupération des caractéristiques
    deps = dlong0/xlong0
    sigm = effnom/aire
!
    if (isot .and. (.not. implex)) then
!       Caractéristiques élastiques a t-
        call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        em = val(1)
!       Caractéristiques élastiques a t+
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        ep = val(1)
!
        call verift(fami, kpg, ksp, 'T', imate, epsth_=depsth)
        depsm = deps-depsth
        call nm1dis(fami, kpg, ksp, imate, em, ep, sigm, depsm, vim, option, &
                    rela_comp, ' ', sigp, vip, dsde)
    else if (cine) then
!       Caractéristiques élastiques a t-
        call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        em = val(1)
!       Caractéristiques élastiques a t+
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        ep = val(1)
!
        call verift(fami, kpg, ksp, 'T', imate, epsth_=depsth)
        depsm = deps-depsth
        call nm1dci(fami, kpg, ksp, imate, em, ep, sigm, depsm, vim, option, &
                    ' ', sigp, vip, dsde)
    else if (relax) then
        call verift(fami, kpg, ksp, 'T', imate, epsth_=depsth)
        depsm = deps-depsth
        call relax_acier_cable(fami, kpg, ksp, imate, sigm, epsm, depsm, vim, option, &
                               ' ', sigp, vip, dsde)
    else if (elas) then
!       Caractéristiques élastiques a t-
        call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        em = val(1)
!       Caractéristiques élastiques a t+
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        ep = val(1)
!
        dsde = ep
        vip(1) = 0.d0
        call verift(fami, kpg, ksp, 'T', imate, epsth_=depsth)
        sigp = ep*(sigm/em+deps-depsth)
    else if (corr) then
!       Caractéristiques élastiques a t-
        call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        em = val(1)
!       Caractéristiques élastiques a t+
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        ep = val(1)
!
        call nm1dco(fami, kpg, ksp, option, imate, ' ', ep, sigm, epsm, deps, &
                    vim, sigp, vip, dsde, carcri, codret)
    else if (implex) then
!       Caractéristiques élastiques a t-
        call rcvalb(fami, kpg, ksp, '-', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        em = val(1)
!       Caractéristiques élastiques a t+
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', &
                    0, ' ', [0.d0], 1, 'E', val, codres, 1)
        ep = val(1)
!
        call lcimpl(fami, kpg, ksp, imate, em, ep, sigm, tmoins, tplus, deps, &
                    vim, option, sigp, vip, dsde)
    else
        ASSERT(ASTER_FALSE)
    end if
!
!   Calcul du coefficient non nul de la matrice tangente
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
        xrig = dsde*aire/xlong0
        klv(1) = xrig
        klv(7) = -xrig
        klv(10) = xrig
    end if
!
!   Calcul des forces nodales
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        effnop = sigp*aire
        fono(1) = -effnop
        fono(4) = effnop
    end if
!
    if (implex) then
        effnop = sigp*aire
        fono(1) = -effnop
        fono(4) = effnop
    end if
!
end subroutine
